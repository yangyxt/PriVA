"""R5/R6 atomic plans, parser fallback and three original product boundaries."""
import ast
import copy
import csv
import importlib.util
from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import acmg_scoring
import priva_rust_gene_weights as bridge
import priva_rust_rank_setup as setup
import priva_rust_pedigree_context as pedigree
import priva_rust_control_multiply as multiply
import priva_rust_ranking_compatibility as compatibility
import priva_rust_terminal_order as terminal


def entry():
    data = pd.DataFrame({'SYMBOL': ['A', 'B', 'C', 'A'], 'sort_index': [.95, .475, -0., .12345678901234568],
        'zygosity_inheritance_mechanism_compatibility': [.9, 1., .9, .9], 'payload': ['original'] * 4})
    data.index = pd.Index(['fourth', 'second', 'third', 'first'], name='source rows')
    data.columns.name = 'source columns'
    data.attrs = {'scope': ['actual post-R4 frame']}
    data.flags.allows_duplicate_labels = False
    return data


def mask(data, nullable=False):
    values = (pd.arrays.BooleanArray(np.array([True, True, False, True]), np.array([False, True, False, True]))
              if nullable else np.array([True, False, True, False]))
    result = pd.Series(values, index=data.index, name='ACMG_quant_score')
    result.attrs = {'scope': ['actual R4 mask']}
    return result


def original(frame, dispensable, relevant, penalty):
    """Execute only the preserved original R5/R6 statements on current buffers."""
    tree = ast.parse(Path(acmg_scoring.__file__).read_text())
    ranker = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'sort_and_rank_variants')
    fallback = next(node for node in ranker.body if isinstance(node, ast.If) and ast.unparse(node.test) == 'native_weights is None')
    wrapper = ast.parse('def original_weights(df, dispensable_gene_list, relevant_gene_list, compatibility_penalty_rows):\n    pass\n').body[0]
    wrapper.body = copy.deepcopy(fallback.body) + ast.parse('return locals()').body
    namespace = dict(acmg_scoring.__dict__)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[wrapper], type_ignores=[])), '<original R5/R6 test excerpt>', 'exec'), namespace)
    checkpoints = {}
    lines = {wrapper.body[1].lineno: 'dispensable', wrapper.body[2].lineno: 'relevant', wrapper.body[-1].lineno: 'compatibility'}
    def trace(frame, event, arg):
        if frame.f_code is namespace['original_weights'].__code__ and event == 'line' and frame.f_lineno in lines:
            checkpoints[lines[frame.f_lineno]] = frame.f_locals['df']['sort_index'].to_numpy(copy=True)
        return trace
    previous = sys.gettrace()
    try:
        sys.settrace(trace)
        result = namespace['original_weights'](frame, dispensable, relevant, penalty)
    finally:
        sys.settrace(previous)
    return result, checkpoints


def response(command, frame, penalty, dispensable, relevant, change=None):
    """Encode actual independent Python products in the existing wire format."""
    data, selected_mask = pickle.loads(pickle.dumps((frame, penalty), protocol=5))
    state, checkpoints = original(data, dispensable, relevant, selected_mask)
    names = list(frame.columns); planned = list(names); steps = []
    gene_lists = {}
    for role, path, target, factor, line in [('dispensable', dispensable, 'dispensable_gene_index', .9, 604),
                                           ('relevant', relevant, 'relevant_gene_index', 1.2, 611)]:
        if path is not None:
            gene_lists[role] = pd.read_table(path, header=None, names=['SYMBOL'], comment='#')['SYMBOL'].tolist()
            position = planned.index(target) if target in planned else len(planned)
            if position == len(planned): planned.append(target)
            chosen = frame['SYMBOL'].isin(state[role + '_genes']).to_numpy()
            steps.append((role, target, line, position, chosen, np.where(chosen, factor, 1.)))
    chosen = selected_mask.to_numpy(dtype=bool, na_value=False)
    steps.append(('compatibility', 'zygosity_inheritance_mechanism_compatibility', 613,
                  names.index('zygosity_inheritance_mechanism_compatibility'), chosen,
                  frame['zygosity_inheritance_mechanism_compatibility'].to_numpy()))
    rows = []
    for step, (role, target, line, position, chosen, multipliers) in enumerate(steps):
        for row, (selected, factor, value) in enumerate(zip(chosen, multipliers.view(np.uint64), checkpoints[role].view(np.uint64))):
            rows.append(dict(step=str(step), source_position=str(row), index_kind='opaque', index_value=str(row),
                             selected=str(int(selected)), multiplier_bits=f'{factor:016x}', score_bits=f'{value:016x}'))
    meta = [('axis', 'caller-owned'), ('score_column_position', str(names.index('sort_index'))), ('step_count', str(len(steps)))]
    for step, (_, target, line, position, _, _) in enumerate(steps):
        meta.extend([(f'step.{step}.name', target), (f'step.{step}.line', str(line)),
                     (f'step.{step}.multiplier_column_position', str(position))])
    for role, path in [('dispensable', dispensable), ('relevant', relevant)]:
        meta.append((role + '.enabled', str(path is not None).lower()))
        if path is not None:
            values = gene_lists[role]
            meta.append((role + '.count', str(len(values))))
            meta.extend((f'{role}.gene.{i}', value) for i, value in enumerate(values))
    if change: change(rows, meta)
    fields = ['step', 'source_position', 'index_kind', 'index_value', 'selected', 'multiplier_bits', 'score_bits']
    with Path(command[command.index('--output') + 1]).open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(fields)
        writer.writerows([row.get(field) for field in fields] + row.get('_extra', []) for row in rows)
    with Path(command[command.index('--metadata-output') + 1]).open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['name', 'value']); writer.writerows(meta)
    return SimpleNamespace(returncode=0)


def publish(data, penalty, plan, dispensable, relevant):
    if dispensable is not None:
        data['dispensable_gene_index'] = 1.0
        data.loc[bridge.gene_membership(data, plan['dispensable_selected']), 'dispensable_gene_index'] = .9
        data.loc[:, 'sort_index'] = plan['dispensable_product']
    if relevant is not None:
        data['relevant_gene_index'] = 1.0
        data.loc[bridge.gene_membership(data, plan['relevant_selected']), 'relevant_gene_index'] = 1.2
        data.loc[:, 'sort_index'] = plan['relevant_product']
    data.loc[penalty, 'sort_index'] = plan['compatibility_product']
    return data


def rank_frame():
    result = pd.DataFrame({'BIOTYPE': ['protein_coding'] * 4, 'ACMG_quant_score': [.95] * 4,
        'var_plausible_patho_mechs': ['dominant_LOF', '', 'dominant_DN', ''] ,
        'variant_effect': ['predicted_LOF_high_confidence'] * 4, 'variant_lof_score': [1] * 4,
        'chrom': ['1'] * 4, 'pos': [1, 2, 3, 4], 'ref': ['A'] * 4, 'alt': ['G'] * 4,
        'SYMBOL': ['A', 'B', 'C', 'A'], 'variant_id': ['a', 'b', 'c', 'd']})
    return result


class GeneWeightsBridgeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='r5-r6-test-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.dispensable = self.root / 'dispensable.txt'; self.dispensable.write_text('# gene\nA\nB\nA\n')
        self.relevant = self.root / 'relevant.txt'; self.relevant.write_text('B\nC\n')

    def planned(self, data, penalty, dispensable=None, relevant=None, change=None):
        with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs:
                response(command, data, penalty, dispensable, relevant, change)):
            return bridge.try_gene_weights(data, dispensable, relevant, penalty, 'native')

    def test_default_and_terminal_only_never_enter_gene_bridge(self):
        with patch.object(bridge, 'try_gene_weights') as native, patch.object(terminal, 'try_terminal_order', return_value=None):
            wanted = acmg_scoring.sort_and_rank_variants(rank_frame(), None, None, dispensable_gene_list=None)
            got = acmg_scoring.sort_and_rank_variants(rank_frame(), None, None, dispensable_gene_list=None, terminal_order_binary='tail')
        native.assert_not_called(); pd.testing.assert_frame_equal(got, wanted, check_exact=True)

    def test_plan_atomic_bool_nullable_and_three_products_exact(self):
        for nullable in (False, True):
            data = entry(); penalty = mask(data, nullable)
            expected, _ = original(*pickle.loads(pickle.dumps((data, self.dispensable, self.relevant, penalty), protocol=5)))
            before = pickle.dumps((data, penalty), protocol=5)
            plan = self.planned(data, penalty, self.dispensable, self.relevant)
            self.assertIsNotNone(plan); self.assertEqual(pickle.dumps((data, penalty), protocol=5), before)
            actual = publish(data, penalty, plan, self.dispensable, self.relevant)
            pd.testing.assert_frame_equal(actual, expected['df'], check_exact=True)
            self.assertEqual(actual.attrs, expected['df'].attrs)
            np.testing.assert_array_equal(actual['sort_index'].to_numpy().view(np.uint64), expected['df']['sort_index'].to_numpy().view(np.uint64))
            self.assertEqual(plan['dispensable_genes'], {'A', 'B'})
            self.assertEqual(plan['relevant_genes'], {'B', 'C'})

    def test_none_disabled_empty_file_enabled_and_existing_target_positions(self):
        empty = self.root / 'empty.txt'; empty.write_text('')
        for dispensable, relevant in ((None, None), (empty, None), (None, self.relevant), (self.dispensable, self.relevant)):
            data = entry(); data.insert(1, 'relevant_gene_index', ['old'] * 4); penalty = mask(data)
            expected, _ = original(*pickle.loads(pickle.dumps((data, dispensable, relevant, penalty), protocol=5)))
            plan = self.planned(data, penalty, dispensable, relevant); self.assertIsNotNone(plan)
            pd.testing.assert_frame_equal(publish(data, penalty, plan, dispensable, relevant), expected['df'], check_exact=True)
            self.assertEqual('dispensable_genes' in plan, dispensable is not None)

    def test_empty_object_mask_and_no_symbol_when_gene_lists_disabled(self):
        data = entry().iloc[:0].drop(columns='SYMBOL'); penalty = pd.Series([], dtype=object, index=data.index)
        plan = self.planned(data, penalty); self.assertIsNotNone(plan)
        self.assertEqual(len(plan['compatibility_product']), 0)
        data = entry().drop(columns='SYMBOL'); penalty = mask(data)
        self.assertIsNotNone(self.planned(data, penalty))

    def test_raw_nullable_mask_and_payload_placeholders_are_transport_exact(self):
        class Payload:
            def __str__(self): raise AssertionError('payload str callback')
            def __repr__(self): raise AssertionError('payload repr callback')
        data = entry(); data['payload'] = [Payload() for _ in range(4)]; penalty = mask(data, True)
        # Build independent output before putting opaque objects in the frame.
        clean = entry()
        def native(command, **kwargs):
            with Path(command[command.index('--penalty-input') + 1]).open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([(row['raw_data'], row['is_null']) for row in rows], [('1', '0'), ('1', '1'), ('0', '0'), ('1', '1')])
            with Path(command[command.index('--input') + 1]).open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([(row['payload_kind'], row['payload']) for row in rows], [('null', '')] * 4)
            self.assertEqual((kwargs['timeout'], kwargs['check']), (300, False))
            return response(command, clean, mask(clean, True), self.dispensable, self.relevant)
        with patch.object(bridge.subprocess, 'run', side_effect=native):
            self.assertIsNotNone(bridge.try_gene_weights(data, self.dispensable, self.relevant, penalty, 'native'))

    def test_duplicate_gapped_axes_and_series_metadata_survive_publication(self):
        for index in (pd.Index([7, 7, -2, 7], name='duplicate'), pd.Index([8, -2, 90, 4], name='gapped')):
            data = entry(); data.flags.allows_duplicate_labels = True; data.index = index
            data['SYMBOL'].attrs = {'source': ['cached metadata']}
            penalty = mask(data)
            expected, _ = original(*pickle.loads(pickle.dumps((data, self.dispensable, self.relevant, penalty), protocol=5)))
            plan = self.planned(data, penalty, self.dispensable, self.relevant); self.assertIsNotNone(plan)
            pd.testing.assert_frame_equal(publish(data, penalty, plan, self.dispensable, self.relevant), expected['df'], check_exact=True)

    def test_unsupported_axis_dtype_versions_and_custom_metadata_defer(self):
        cases = []
        for target in ('SYMBOL', 'sort_index', 'zygosity_inheritance_mechanism_compatibility'):
            data = entry(); data[target] = data[target].astype('category'); cases.append((data, mask(data)))
            data = entry().drop(columns=target); cases.append((data, mask(data)))
        data = entry(); data.index = pd.MultiIndex.from_tuples([('a', i) for i in range(4)]); cases.append((data, mask(data)))
        data = entry(); penalty = mask(data); penalty.index = data.index.copy(); cases.append((data, penalty))
        data = entry(); cases.append((data, mask(data).astype(object)))
        data = entry(); data.attrs = {'custom': object()}; cases.append((data, mask(data)))
        data = entry(); data['SYMBOL'].attrs = {'custom': object()}; cases.append((data, mask(data)))
        data = entry(); data.flags.allows_duplicate_labels = True
        data = pd.concat([data, data[['sort_index']]], axis=1); cases.append((data, mask(data)))
        for data, penalty in cases:
            with self.subTest(dtypes=list(map(str, data.dtypes))), patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_gene_weights(data, self.dispensable, self.relevant, penalty, 'native')); native.assert_not_called()
        for module in (pd, np):
            data = entry(); penalty = mask(data)
            with patch.object(module, '__version__', 'unqualified'), patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_gene_weights(data, None, None, penalty, 'native')); native.assert_not_called()

    def test_unconsumed_symbol_is_not_read_when_both_lists_disabled(self):
        class Symbol:
            def __str__(self): raise AssertionError('unconsumed SYMBOL')
        data = entry(); data['SYMBOL'] = [Symbol()] * 4; penalty = mask(data)
        clean = entry()
        with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs: response(command, clean, mask(clean), None, None)):
            self.assertIsNotNone(bridge.try_gene_weights(data, None, None, penalty, 'native'))

    def test_custom_paths_links_compression_and_missing_files_keep_original_parser_route(self):
        link = self.root / 'link.txt'; link.symlink_to(self.dispensable.name)
        compressed = self.root / 'compressed.txt.gz'; compressed.write_text('A\n')
        for value in ('', False, 7, pd.NA, object(), self.root / 'missing.txt', link, compressed):
            data = entry(); penalty = mask(data)
            with self.subTest(value=type(value).__name__), patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_gene_weights(data, value, None, penalty, 'native')); native.assert_not_called()

    def test_malformed_step_rows_mask_multipliers_and_metadata_are_atomic(self):
        changes = [lambda rows, meta: rows.reverse(), lambda rows, meta: rows.pop(),
            lambda rows, meta: rows[0].update(step='01'), lambda rows, meta: rows[0].update(source_position='1'),
            lambda rows, meta: rows[0].update(index_value='wrong'), lambda rows, meta: rows[0].update(selected='2'),
            lambda rows, meta: rows[0].update(multiplier_bits='0000000000000000'),
            lambda rows, meta: rows[-1].update(selected='1'), lambda rows, meta: rows[-1].update(multiplier_bits='0000000000000000'),
            lambda rows, meta: rows[0].update(score_bits='nan'), lambda rows, meta: rows[0].update(_extra=['unexpected']),
            lambda rows, meta: meta.reverse(), lambda rows, meta: meta.append(meta[-1]),
            lambda rows, meta: meta.__setitem__(0, ('axis', 'changed')),
            lambda rows, meta: meta.__setitem__(2, ('step_count', '1')),
            lambda rows, meta: meta.__setitem__(5, ('step.0.multiplier_column_position', '0')),
            lambda rows, meta: meta.__setitem__(-1, ('relevant.gene.1', 'NA'))]
        for change in changes:
            data = entry(); penalty = mask(data); before = pickle.dumps((data, penalty), protocol=5)
            with self.subTest(change=change):
                self.assertIsNone(self.planned(data, penalty, self.dispensable, self.relevant, change))
                self.assertEqual(pickle.dumps((data, penalty), protocol=5), before)

    def test_gene_file_metadata_race_rejects_the_whole_plan(self):
        data = entry(); penalty = mask(data); before = pickle.dumps((data, penalty), protocol=5)
        def native(command, **kwargs):
            result = response(command, data, penalty, self.dispensable, self.relevant)
            self.relevant.write_text('C\nB\n'); return result
        with patch.object(bridge.subprocess, 'run', side_effect=native):
            self.assertIsNone(bridge.try_gene_weights(data, self.dispensable, self.relevant, penalty, 'native'))
        self.assertEqual(pickle.dumps((data, penalty), protocol=5), before)

    def test_timeout_nonzero_and_resource_failures_preserve_atomicity(self):
        data = entry(); penalty = mask(data); before = pickle.dumps((data, penalty), protocol=5)
        for error in (FileNotFoundError('missing'), subprocess.TimeoutExpired('native', 300), RuntimeError('transport')):
            with patch.object(bridge.subprocess, 'run', side_effect=error):
                self.assertIsNone(bridge.try_gene_weights(data, None, None, penalty, 'native'))
        with patch.object(bridge.subprocess, 'run', return_value=SimpleNamespace(returncode=1)):
            self.assertIsNone(bridge.try_gene_weights(data, None, None, penalty, 'native'))
        for error in (MemoryError, RecursionError):
            with patch.object(bridge.subprocess, 'run', side_effect=error), self.assertRaises(error):
                bridge.try_gene_weights(data, None, None, penalty, 'native')
        self.assertEqual(pickle.dumps((data, penalty), protocol=5), before)

    def test_full_ranker_fallback_executes_once_on_actual_post_R4_buffers(self):
        events, current = [], []
        original_control, original_masks, original_read = acmg_scoring.compute_control_common, acmg_scoring._variant_mechanism_masks, pd.read_table
        def control(*args):
            events.append('C1'); result = original_control(*args); current.append(result); return result
        def masks(data): events.append('R4'); return original_masks(data)
        def genes(data, first, second, penalty, binary):
            events.append('R5/R6-attempt'); self.assertIs(data, current[0]); self.assertIn(bridge._COMPATIBILITY, data); return None
        def read(path, **kwargs): events.append('parse-dispensable' if path == self.dispensable else 'parse-relevant'); return original_read(path, **kwargs)
        expected = acmg_scoring.sort_and_rank_variants(rank_frame(), None, None, dispensable_gene_list=self.dispensable, relevant_gene_list=self.relevant)
        with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(pedigree, 'try_pedigree_context', return_value=None), patch.object(
                multiply, 'try_control_multiply', return_value=None), patch.object(compatibility, 'try_ranking_compatibility', return_value=None), patch.object(
                bridge, 'try_gene_weights', side_effect=genes), patch.object(terminal, 'try_terminal_order', return_value=None), patch.object(
                acmg_scoring, 'compute_control_common', side_effect=control), patch.object(acmg_scoring, '_variant_mechanism_masks', side_effect=masks), patch.object(pd, 'read_table', side_effect=read):
            actual = acmg_scoring.sort_and_rank_variants(rank_frame(), None, None, dispensable_gene_list=self.dispensable, relevant_gene_list=self.relevant, ranking_binary='native')
        pd.testing.assert_frame_equal(actual, expected, check_exact=True)
        self.assertEqual(events, ['C1', 'R4', 'R5/R6-attempt', 'parse-dispensable', 'parse-relevant'])

    def test_full_ranker_native_publication_keeps_three_separate_assignments(self):
        expected = acmg_scoring.sort_and_rank_variants(rank_frame(), None, None,
            dispensable_gene_list=self.dispensable, relevant_gene_list=self.relevant)
        events, current = [], []
        original_control, actual_try = acmg_scoring.compute_control_common, bridge.try_gene_weights
        loc_type = type(rank_frame().loc); original_assign = loc_type.__setitem__
        def control(*args):
            result = original_control(*args); current.append(result); return result
        def genes(data, first, second, penalty, binary):
            self.assertIs(data, current[0]); events.append('plan')
            with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs:
                    response(command, data, penalty, first, second)):
                result = actual_try(data, first, second, penalty, binary)
            self.assertIsNotNone(result); return result
        def assign(indexer, key, value):
            if events and indexer.obj is current[0] and isinstance(key, tuple) and key[1] == 'sort_index':
                events.append('sort_index-assignment')
            return original_assign(indexer, key, value)
        with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(pedigree, 'try_pedigree_context', return_value=None), patch.object(
                multiply, 'try_control_multiply', return_value=None), patch.object(compatibility, 'try_ranking_compatibility', return_value=None), patch.object(
                bridge, 'try_gene_weights', side_effect=genes), patch.object(terminal, 'try_terminal_order', return_value=None), patch.object(
                acmg_scoring, 'compute_control_common', side_effect=control), patch.object(loc_type, '__setitem__', assign):
            actual = acmg_scoring.sort_and_rank_variants(rank_frame(), None, None,
                dispensable_gene_list=self.dispensable, relevant_gene_list=self.relevant, ranking_binary='native')
        pd.testing.assert_frame_equal(actual, expected, check_exact=True)
        self.assertEqual(events, ['plan', 'sort_index-assignment', 'sort_index-assignment', 'sort_index-assignment'])

    def test_second_read_failure_keeps_first_gene_product_and_warning_order(self):
        outcomes = []
        for native in (False, True):
            events, current = [], []; original_control = acmg_scoring.compute_control_common
            def control(*args):
                result = original_control(*args); current.append(result); events.append('C1'); return result
            def genes(*args): events.append('R5/R6-attempt'); return None
            with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(pedigree, 'try_pedigree_context', return_value=None), patch.object(
                    multiply, 'try_control_multiply', return_value=None), patch.object(compatibility, 'try_ranking_compatibility', return_value=None), patch.object(
                    bridge, 'try_gene_weights', side_effect=genes), patch.object(acmg_scoring, 'compute_control_common', side_effect=control), patch.object(
                    acmg_scoring.logger, 'warning', side_effect=lambda *args: events.append('warning')):
                with self.assertRaises(FileNotFoundError) as failure:
                    acmg_scoring.sort_and_rank_variants(rank_frame(), None, None, dispensable_gene_list=self.dispensable,
                        relevant_gene_list=self.root / 'missing.txt', **({'ranking_binary': 'native'} if native else {}))
            self.assertEqual(events, ['warning', 'C1'] + (['R5/R6-attempt'] if native else []))
            self.assertIn('dispensable_gene_index', current[0]); self.assertNotIn('relevant_gene_index', current[0])
            outcomes.append((failure.exception.args, current[0]))
        self.assertEqual(outcomes[0][0], outcomes[1][0]); pd.testing.assert_frame_equal(outcomes[0][1], outcomes[1][1], check_exact=True)

    def test_parser_projection_requires_one_complete_enabled_replacement(self):
        spec = importlib.util.spec_from_file_location('r5_parser_guard', ROOT / 'rust/tests/terminal_order_bridge_compare.py')
        harness = importlib.util.module_from_spec(spec); spec.loader.exec_module(harness)
        project = harness._gene_parser_projection
        for roles in ([], ['dispensable'], ['relevant'], ['dispensable', 'relevant']):
            self.assertEqual(project([], roles, [True], roles), roles)
            self.assertIs(project(roles, [], [False], roles), roles)
        for values in ((['dispensable'], ['dispensable'], [True], ['dispensable']),
                       ([], [], [True], ['dispensable']), ([], ['relevant'], [True], ['dispensable']),
                       ([], [], [False, False], []), ([], [], [], ['other'])):
            with self.assertRaises(AssertionError): project(*values)


if __name__ == '__main__':
    unittest.main()
