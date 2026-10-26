"""R4 atomic plans, original current-buffer fallback and publication order."""
import ast
import copy
import csv
import importlib.util
from pathlib import Path
import pickle
import subprocess
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import acmg_scoring
import priva_rust_ranking_compatibility as bridge
import priva_rust_control_multiply as multiply
import priva_rust_pedigree_context as pedigree
import priva_rust_rank_setup as setup
import priva_rust_terminal_order as terminal


def entry():
    data = pd.DataFrame({'var_plausible_patho_mechs': ['', 'dominant_LOF', 'x_linked_recessive_LOF', 'dominant_DN'],
        'variant_effect': ['predicted_LOF_high_confidence'] * 4, 'variant_lof_score': [1] * 4,
        'ACMG_quant_score': [.95] * 4, 'sort_index': [.475] * 4, 'chrom': ['X'] * 4})
    data.index = pd.Index(['fourth', 'second', 'third', 'first'], name='source rows')
    data.columns.name = 'source columns'
    data.attrs = {'scope': ['current post-R3 frame']}
    data.flags.allows_duplicate_labels = False
    return data


def mask(data, nullable=False):
    values = (pd.arrays.BooleanArray(np.array([True, True, False, True]), np.array([False, True, False, True]))
              if nullable else np.array([True, False, True, True]))
    result = pd.Series(values, index=data.index, name='PROBAND')
    result.attrs = {'scope': ['actual R2 mask']}
    result.flags.allows_duplicate_labels = False
    return result


def source_plan(data, het, male=True, female=False):
    """Execute the preserved original R4 statements on private test buffers."""
    tree = ast.parse(Path(acmg_scoring.__file__).read_text())
    ranker = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'sort_and_rank_variants')
    fallback = next(node for node in ranker.body if isinstance(node, ast.If)
                    and ast.unparse(node.test) == 'native_compatibility is None')
    wrapper = ast.parse('def original_r4(df, proband_is_male, proband_is_female, proband_het):\n    pass\n').body[0]
    wrapper.body = copy.deepcopy(fallback.body) + ast.parse('return locals()').body
    namespace = dict(acmg_scoring.__dict__)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[wrapper], type_ignores=[])), '<original R4 test excerpt>', 'exec'), namespace)
    result = namespace['original_r4'](data, male, female, het)
    return {name: result[name] for name in ('mechanism_masks', 'chrom', 'x_locus', 'has_single_alt_lof_compatible',
        'exact_gof', 'asserted_lof_effect', 'dn_without_hi_or_ambiguous', 'heterozygous_lof_incompatible',
        'high_acmg_lof', 'compatibility_penalty_rows', 'pext_modulation_rows')}


def response(command, plan, change=None):
    """Encode independent original-Python outputs as the fixed native contract."""
    with Path(command[command.index('--input') + 1]).open(newline='') as stream:
        source = csv.reader(stream); header = next(source); input_rows = list(source)
    names = header[3::2]
    masks = plan['mechanism_masks']
    roles = ('chrom', 'x_locus', 'has_single_alt_lof_compatible', 'exact_gof', 'asserted_lof_effect',
             'dn_without_hi_or_ambiguous', 'heterozygous_lof_incompatible', 'high_acmg_lof',
             'compatibility_penalty_rows', 'compatibility_penalty_rows_null', 'pext_modulation_rows')
    fields = ['source_position', 'index_kind', 'index_value', 'chrom_hex', 'compatibility_bits',
              *masks, *[name for name in roles if name != 'chrom']]
    penalty = plan['compatibility_penalty_rows']
    nullable = str(penalty.dtype) == 'boolean'
    raw = penalty.array._data if nullable else penalty.to_numpy(copy=False)
    null = penalty.array._mask if nullable else np.zeros(len(penalty), dtype=bool)
    rows = []
    for position, record in enumerate(input_rows):
        row = dict(source_position=str(position), index_kind=record[0], index_value=record[1],
            chrom_hex=plan['chrom'].iloc[position].encode().hex(),
            compatibility_bits='3feccccccccccccd' if raw[position] and not null[position] else '3ff0000000000000')
        row.update((name, str(int(value.iloc[position]))) for name, value in masks.items())
        for name in roles[1:]:
            value = null[position] if name == 'compatibility_penalty_rows_null' else (
                    raw[position] if name == 'compatibility_penalty_rows' else plan[name].iloc[position])
            row[name] = str(int(value))
        rows.append(row)
    target = 'zygosity_inheritance_mechanism_compatibility'
    message = ('Zygosity/inheritance/mechanism compatibility signals for ranking: '
        f"hub_profile_rows={masks['modern_profile'].sum()}, asserted_lof_effect={pd.Series(plan['asserted_lof_effect'], index=penalty.index).sum()}, "
        f"dominant_DN_without_HI_or_ambiguous={plan['dn_without_hi_or_ambiguous'].sum()}, "
        f"heterozygous_lof_incompatible={plan['heterozygous_lof_incompatible'].sum()}, compatibility_penalty_rows={penalty.sum()}")
    meta = [('target_source_column_position', str(names.index(target) if target in names else len(names))),
            ('penalty_dtype', str(penalty.dtype)), ('axis', 'caller-owned'), ('log_message', message)]
    for name in (*masks, *roles):
        value = penalty if name == 'compatibility_penalty_rows_null' else masks[name] if name in masks else plan[name]
        if name in ('chrom', 'x_locus'):
            origin = 'chrom' if 'chrom' in names else None
        elif name in ('is_exact_gof', 'is_exact_dn', 'is_exact_nonlof', 'exact_gof'):
            origin = None
        elif name in ('is_predicted_lof', 'is_uncertain', 'asserted_lof_effect', 'heterozygous_lof_incompatible'):
            origin = 'variant_effect'
        elif name in ('high_acmg_lof', 'pext_modulation_rows', 'compatibility_penalty_rows', 'compatibility_penalty_rows_null'):
            origin = 'ACMG_quant_score'
        else:
            origin = 'var_plausible_patho_mechs'
        meta.extend(((name + '.series_name_kind', 'none' if value.name is None else 'str'),
                     (name + '.series_name', '' if value.name is None else value.name),
                     (name + '.metadata_origin', 'fresh' if origin is None else 'column:' + str(names.index(origin)))))
    if change:
        change(rows, meta)
    with Path(command[command.index('--output') + 1]).open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(fields)
        writer.writerows([row.get(field) for field in fields] + row.get('_extra', []) for row in rows)
    with Path(command[command.index('--metadata-output') + 1]).open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['name', 'value']); writer.writerows(meta)
    return SimpleNamespace(returncode=0)


def rank_frame():
    data = entry()
    for name, values in dict(BIOTYPE=['protein_coding'] * 4, pos=[1, 2, 3, 4], ref=['A'] * 4,
                              alt=['G'] * 4, SYMBOL=['G'] * 4, variant_id=['a', 'b', 'c', 'd']).items():
        data[name] = values
    return data


def rank(data, **kwargs):
    return acmg_scoring.sort_and_rank_variants(data, None, None, dispensable_gene_list=None, **kwargs)


class CompatibilityBridgeTests(unittest.TestCase):
    def test_default_and_r8_only_never_enter_r4_adapter(self):
        with patch.object(bridge, 'try_ranking_compatibility') as native, patch.object(
                terminal, 'try_terminal_order', return_value=None):
            wanted = rank(rank_frame())
            pd.testing.assert_frame_equal(rank(rank_frame(), terminal_order_binary='tail'), wanted, check_exact=True)
        native.assert_not_called()

    def test_plan_is_atomic_and_preserves_names_metadata_aliases_and_raw_bits(self):
        for nullable in (False, True):
            data = entry(); het = mask(data, nullable); expected = source_plan(data, het)
            before = pickle.dumps((data, het), protocol=5)
            with self.subTest(nullable=nullable), patch.object(bridge.subprocess, 'run',
                    side_effect=lambda command, **kwargs: response(command, expected)) as native:
                actual = bridge.try_ranking_compatibility(data, True, False, het, 'native')
            self.assertIsNotNone(actual)
            self.assertEqual(list(actual), list(expected))
            for name in expected:
                if name == 'mechanism_masks':
                    self.assertEqual(list(actual[name]), list(expected[name]))
                    for key in expected[name]:
                        pd.testing.assert_series_equal(actual[name][key], expected[name][key], check_exact=True)
                        self.assertEqual(actual[name][key].attrs, expected[name][key].attrs)
                else:
                    pd.testing.assert_series_equal(actual[name], expected[name], check_exact=True)
                    self.assertEqual(actual[name].attrs, expected[name].attrs)
                    self.assertEqual(actual[name].flags.allows_duplicate_labels, expected[name].flags.allows_duplicate_labels)
            self.assertIs(actual['exact_gof'], actual['mechanism_masks']['is_exact_gof'])
            self.assertIs(actual['asserted_lof_effect'], actual['mechanism_masks']['is_predicted_lof'])
            self.assertIs(actual['pext_modulation_rows'], actual['high_acmg_lof'])
            self.assertIs(actual['chrom'].index, data.index)
            if nullable:
                for name in ('_data', '_mask'):
                    np.testing.assert_array_equal(getattr(actual['compatibility_penalty_rows'].array, name),
                                                  getattr(expected['compatibility_penalty_rows'].array, name))
            self.assertEqual(pickle.dumps((data, het), protocol=5), before)
            self.assertEqual((native.call_args.kwargs['timeout'], native.call_args.kwargs['check']), (300, False))

    def test_nullable_mask_transport_retains_hidden_data(self):
        data = entry(); het = mask(data, True); expected = source_plan(data, het)
        def native(command, **kwargs):
            with Path(command[command.index('--heterozygosity-input') + 1]).open(newline='') as stream:
                records = list(csv.DictReader(stream))
            self.assertEqual([(row['raw_data'], row['is_null']) for row in records],
                             [('1', '0'), ('1', '1'), ('0', '0'), ('1', '1')])
            self.assertEqual(command[command.index('--heterozygosity-dtype') + 1], 'boolean')
            return response(command, expected)
        with patch.object(bridge.subprocess, 'run', side_effect=native):
            self.assertIsNotNone(bridge.try_ranking_compatibility(data, True, False, het, 'native'))

    def test_empty_object_context_and_missing_chrom_preserve_fresh_metadata(self):
        data = entry().iloc[:0].drop(columns='chrom')
        het = pd.Series([], dtype=object, index=data.index)
        expected = source_plan(data, het)
        with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs: response(command, expected)):
            result = bridge.try_ranking_compatibility(data, False, False, het, 'native')
        self.assertIsNotNone(result)
        self.assertEqual(str(result['compatibility_penalty_rows'].dtype), 'bool')
        self.assertEqual(result['chrom'].attrs, {}); self.assertTrue(result['chrom'].flags.allows_duplicate_labels)
        self.assertIsNone(result['chrom'].name)

    def test_payload_objects_are_not_transported_or_called(self):
        class Payload:
            def __str__(self): raise AssertionError('payload string callback')
            def __repr__(self): raise AssertionError('payload repr callback')
            def __deepcopy__(self, memo): raise AssertionError('payload deepcopy callback')
        data = entry(); data['payload'] = [Payload() for _ in range(len(data))]
        het = mask(data); expected = source_plan(data, het)
        def native(command, **kwargs):
            with Path(command[command.index('--input') + 1]).open(newline='') as stream:
                records = list(csv.DictReader(stream))
            self.assertEqual([(row['payload_kind'], row['payload']) for row in records], [('null', '')] * 4)
            return response(command, expected)
        with patch.object(bridge.subprocess, 'run', side_effect=native):
            self.assertIsNotNone(bridge.try_ranking_compatibility(data, True, False, het, 'native'))

    def test_unsupported_versions_axes_labels_and_consumed_storage_defer_before_native(self):
        cases = []
        for name in ('var_plausible_patho_mechs', 'variant_effect', 'ACMG_quant_score'):
            data = entry().drop(columns=name); cases.append((data, mask(data)))
        for name in ('chrom', 'ACMG_quant_score', 'variant_lof_score'):
            data = entry(); data[name] = data[name].astype('category'); cases.append((data, mask(data)))
        data = entry(); data.flags.allows_duplicate_labels = True
        data = pd.concat([data, data[['variant_effect']]], axis=1); cases.append((data, mask(data)))
        data = entry(); het = mask(data); het.index = data.index.copy(); cases.append((data, het))
        data = entry(); het = mask(data).astype(object); cases.append((data, het))
        for data, het in cases:
            with self.subTest(columns=list(data.columns)), patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_ranking_compatibility(data, True, False, het, 'native')); native.assert_not_called()
        data = entry(); het = mask(data)
        for module in (pd, np):
            with patch.object(module, '__version__', 'unsupported'), patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_ranking_compatibility(data, True, False, het, 'native')); native.assert_not_called()

    def test_numpy_sex_scalar_and_custom_consumed_cells_defer_without_callbacks(self):
        data = entry(); het = mask(data)
        with patch.object(bridge.subprocess, 'run') as native:
            self.assertIsNone(bridge.try_ranking_compatibility(data, np.bool_(True), False, het, 'native')); native.assert_not_called()
        class Consumed:
            def __str__(self): raise AssertionError('custom consumed callback')
        data['variant_effect'] = [Consumed()] * 4
        with patch.object(bridge.subprocess, 'run') as native:
            self.assertIsNone(bridge.try_ranking_compatibility(data, True, False, het, 'native')); native.assert_not_called()

    def test_malformed_output_is_rejected_without_partial_publication(self):
        changes = [lambda rows, meta: rows.reverse(), lambda rows, meta: rows.pop(),
            lambda rows, meta: rows[0].update(source_position='01'),
            lambda rows, meta: rows[0].update(index_value='wrong'),
            lambda rows, meta: rows[0].update(chrom_hex='ff'),
            lambda rows, meta: rows[0].update(compatibility_bits='0000000000000000'),
            lambda rows, meta: rows[0].update(high_acmg_lof='2'),
            lambda rows, meta: rows[0].update(_extra=['unexpected']),
            lambda rows, meta: rows[0].update(exact_gof='1'),
            lambda rows, meta: rows[0].update(compatibility_penalty_rows_null='1')]
        for change in changes:
            data = entry(); het = mask(data); expected = source_plan(data, het)
            before = pickle.dumps((data, het), protocol=5)
            with self.subTest(change=change), patch.object(bridge.subprocess, 'run',
                    side_effect=lambda command, **kwargs: response(command, expected, change)):
                self.assertIsNone(bridge.try_ranking_compatibility(data, True, False, het, 'native'))
            self.assertEqual(pickle.dumps((data, het), protocol=5), before)
            self.assertNotIn('zygosity_inheritance_mechanism_compatibility', data)

    def test_metadata_target_names_origins_and_INFO_are_all_required(self):
        for key in ('target_source_column_position', 'penalty_dtype', 'axis', 'log_message',
                    'is_predicted_lof.series_name', 'chrom.metadata_origin', 'high_acmg_lof.series_name_kind'):
            data = entry(); het = mask(data); expected = source_plan(data, het)
            def change(rows, meta):
                position = next(i for i, row in enumerate(meta) if row[0] == key)
                meta[position] = (key, 'changed')
            with self.subTest(key=key), patch.object(bridge.subprocess, 'run',
                    side_effect=lambda command, **kwargs: response(command, expected, change)):
                self.assertIsNone(bridge.try_ranking_compatibility(data, True, False, het, 'native'))
        for change in (lambda rows, meta: meta.append(meta[-1]), lambda rows, meta: meta.reverse()):
            with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs: response(command, expected, change)):
                self.assertIsNone(bridge.try_ranking_compatibility(data, True, False, het, 'native'))

    def test_native_errors_and_timeouts_are_atomic_but_resource_errors_propagate(self):
        data = entry(); het = mask(data); before = pickle.dumps((data, het), protocol=5)
        for error in (FileNotFoundError('missing'), subprocess.TimeoutExpired('native', 300), RuntimeError('native failure')):
            with self.subTest(error=error), patch.object(bridge.subprocess, 'run', side_effect=error):
                self.assertIsNone(bridge.try_ranking_compatibility(data, True, False, het, 'native'))
        with patch.object(bridge.subprocess, 'run', return_value=SimpleNamespace(returncode=1)):
            self.assertIsNone(bridge.try_ranking_compatibility(data, True, False, het, 'native'))
        for error in (MemoryError, RecursionError):
            with patch.object(bridge.subprocess, 'run', side_effect=error), self.assertRaises(error):
                bridge.try_ranking_compatibility(data, True, False, het, 'native')
        self.assertEqual(pickle.dumps((data, het), protocol=5), before)

    def test_current_frame_fallback_runs_original_helper_once_after_C1_R3(self):
        data = rank_frame(); expected = rank(data); events, current = [], []
        original_control, original_masks = acmg_scoring.compute_control_common, acmg_scoring._variant_mechanism_masks
        def control(*args):
            events.append('C1'); result = original_control(*args); current.append(result); return result
        def product(actual, binary):
            events.append('R3'); self.assertIs(actual, current[0]); return None
        def compatibility(actual, male, female, het, binary):
            events.append('R4-attempt'); self.assertIs(actual, current[0]); return None
        def mechanisms(actual):
            events.append('R4-original'); self.assertIs(actual, current[0]); return original_masks(actual)
        with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(pedigree, 'try_pedigree_context', return_value=None), patch.object(
                multiply, 'try_control_multiply', side_effect=product), patch.object(bridge, 'try_ranking_compatibility', side_effect=compatibility), patch.object(
                acmg_scoring, 'compute_control_common', side_effect=control), patch.object(acmg_scoring, '_variant_mechanism_masks', side_effect=mechanisms), patch.object(
                terminal, 'try_terminal_order', return_value=None):
            actual = rank(data, ranking_binary='native')
        pd.testing.assert_frame_equal(actual, expected, check_exact=True)
        self.assertEqual(events, ['C1', 'R3', 'R4-attempt', 'R4-original'])

    def test_accepted_plan_keeps_original_INFO_before_both_assignments(self):
        data = rank_frame(); expected = rank(data); events, current = [], []
        original_control, original_setitem = acmg_scoring.compute_control_common, pd.DataFrame.__setitem__
        loc_type = type(data.loc)
        original_loc_setitem = loc_type.__setitem__
        def control(*args):
            events.append('C1'); result = original_control(*args); current.append(result); return result
        def compatibility(actual, male, female, het, binary):
            self.assertIs(actual, current[0]); self.assertNotIn(bridge._TARGET, actual)
            events.append('R4-plan'); return source_plan(actual, het, male, female)
        def assign(actual, key, value):
            if key == bridge._TARGET:
                events.append(('assignment', value))
            return original_setitem(actual, key, value)
        def assign_loc(indexer, key, value):
            if isinstance(key, tuple) and key[1] == bridge._TARGET:
                events.append(('loc assignment', value))
            return original_loc_setitem(indexer, key, value)
        with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(pedigree, 'try_pedigree_context', return_value=None), patch.object(
                multiply, 'try_control_multiply', return_value=None), patch.object(bridge, 'try_ranking_compatibility', side_effect=compatibility), patch.object(
                acmg_scoring, 'compute_control_common', side_effect=control), patch.object(acmg_scoring.logger, 'info', side_effect=lambda *args: events.append('INFO')), patch.object(
                pd.DataFrame, '__setitem__', assign), patch.object(loc_type, '__setitem__', assign_loc), patch.object(
                terminal, 'try_terminal_order', return_value=None):
            actual = rank(data, ranking_binary='native')
        pd.testing.assert_frame_equal(actual, expected, check_exact=True)
        self.assertEqual(events, ['C1', 'R4-plan', 'INFO', ('assignment', 1.0), ('loc assignment', 0.9)])
        self.assertEqual(actual[bridge._TARGET].tolist(), expected[bridge._TARGET].tolist())

    def test_original_helper_failure_and_warning_order_survive_rejected_plan(self):
        failures = []
        for native in (False, True):
            data = rank_frame().drop(columns='variant_effect'); events = []
            original = acmg_scoring.compute_control_common
            def control(*args): events.append('C1'); return original(*args)
            with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(pedigree, 'try_pedigree_context', return_value=None), patch.object(
                    multiply, 'try_control_multiply', return_value=None), patch.object(bridge, 'try_ranking_compatibility', return_value=None), patch.object(
                    acmg_scoring, 'compute_control_common', side_effect=control), patch.object(acmg_scoring.logger, 'warning', side_effect=lambda *args: events.append('warning')), patch.object(
                    acmg_scoring.logger, 'info', side_effect=lambda *args: events.append('INFO')):
                with self.assertRaises(KeyError) as error:
                    rank(data, **({'ranking_binary': 'native'} if native else {}))
            failures.append((type(error.exception), error.exception.args))
            self.assertEqual(events, ['warning', 'C1'])
        self.assertEqual(failures[0], failures[1])

    def test_INFO_failure_precedes_compatibility_publication(self):
        data = rank_frame(); current = []
        original = acmg_scoring.compute_control_common
        def control(*args):
            result = original(*args); current.append(result); return result
        def compatibility(actual, male, female, het, binary):
            return source_plan(actual, het, male, female)
        with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(pedigree, 'try_pedigree_context', return_value=None), patch.object(
                multiply, 'try_control_multiply', return_value=None), patch.object(bridge, 'try_ranking_compatibility', side_effect=compatibility), patch.object(
                acmg_scoring, 'compute_control_common', side_effect=control), patch.object(acmg_scoring.logger, 'info', side_effect=RuntimeError('INFO failure')):
            with self.assertRaisesRegex(RuntimeError, 'INFO failure'):
                rank(data, ranking_binary='native')
        self.assertEqual(len(current), 1); self.assertNotIn(bridge._TARGET, current[0])

    def test_call_projection_requires_one_native_replacement_at_exact_source_position(self):
        spec = importlib.util.spec_from_file_location('R4_source_call_guard', ROOT / 'rust/tests/terminal_order_bridge_compare.py')
        harness = importlib.util.module_from_spec(spec); spec.loader.exec_module(harness)
        project = harness._source_call_projection
        self.assertEqual(project(['compute_control_common', 'coerce_pext_value'],
            [dict(position=1, name='_variant_mechanism_masks')], [True]),
            ['compute_control_common', '_variant_mechanism_masks', 'coerce_pext_value'])
        raw = ['compute_control_common', '_variant_mechanism_masks']
        self.assertIs(project(raw, [], [False]), raw)
        invalid = [([], [], [True]), (raw, [dict(position=1, name='_variant_mechanism_masks')], [True]),
            (['compute_control_common'], [dict(position=0, name='_variant_mechanism_masks')], [True]),
            (['compute_control_common'], [dict(position=1, name='other')], [True]),
            (['compute_control_common', 'other'], [dict(position=1, name='_variant_mechanism_masks')], [True]),
            (raw, [], [False, False])]
        for arguments in invalid:
            with self.subTest(arguments=arguments), self.assertRaises(AssertionError): project(*arguments)


if __name__ == '__main__':
    unittest.main()
