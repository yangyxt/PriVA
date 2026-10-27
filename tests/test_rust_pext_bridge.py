"""R7 atomic planning, actual map boxing and ordered current-buffer publication."""
import ast
import copy
import csv
from contextlib import ExitStack, contextmanager
from pathlib import Path
import pickle
import struct
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
import priva_rust_pext as bridge
import priva_rust_rank_setup as setup
import priva_rust_pedigree_context as pedigree
import priva_rust_control_multiply as multiply
import priva_rust_ranking_compatibility as compatibility
import priva_rust_gene_weights as weights
import priva_rust_terminal_order as terminal


def entry():
    result = pd.DataFrame({'PEXT_MEAN': ['.02', '.1', None, '-0&0'],
        'sort_index': [.95, .475, -0., .12345678901234568], 'payload': ['untouched'] * 4})
    result.index = pd.Index(['fourth', 'second', 'third', 'first'], name='source rows')
    result.columns.name = 'source columns'
    result.attrs = {'scope': ['actual post-gene-weight frame']}
    result.flags.allows_duplicate_labels = False
    return result


def mask(frame, nullable=False):
    array = (pd.arrays.BooleanArray(np.array([True, True, False, True]), np.array([False, True, False, True]))
             if nullable else np.array([True, False, True, False]))
    result = pd.Series(array, index=frame.index, name='ACMG_quant_score')
    result.attrs = {'scope': ['R4 high-ACMG PEXT mask']}
    return result


def original(frame, penalty, column='PEXT_MEAN', cutoff=.1, floor=.8, shape=.5):
    """Execute the preserved source arithmetic and both assignment statements."""
    tree = ast.parse(Path(acmg_scoring.__file__).read_text())
    definitions = [node for node in ast.walk(tree) if isinstance(node, ast.If)
                   and ast.unparse(node.test) == 'native_pext is None']
    assert len(definitions) == 1
    wrapper = ast.parse('def original_pext(df, pext_modulation_rows, pext_col_to_use, '
        'pext_low_expression_cutoff, pext_penalty_floor, pext_penalty_shape):\n    pass\n').body[0]
    wrapper.body = copy.deepcopy(definitions[0].body) + ast.parse('return locals()').body
    namespace = dict(acmg_scoring.__dict__)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[wrapper], type_ignores=[])),
        '<original R7 test excerpt>', 'exec'), namespace)
    return namespace['original_pext'](frame, penalty, column, cutoff, floor, shape)


def response(command, frame, penalty, column='PEXT_MEAN', tissues='', cutoff=.1, floor=.8, shape=.5, change=None):
    """Independent original Python outputs encoded in the qualified native wire format."""
    data, chosen = pickle.loads(pickle.dumps((frame, penalty), protocol=5))
    state = original(data, chosen, column, cutoff, floor, shape)
    numeric, scaled = state['pext_numeric'], state['scaled_pext']
    power = np.power(scaled, shape)
    product = (1.0 - floor) * power
    summed = floor + product
    outputs = [numeric, scaled, power, product, summed, data['pext_sort_index'], data['sort_index']]
    bits = [np.asarray(value, dtype=np.float64).view(np.uint64) for value in outputs]
    selected = penalty.to_numpy(dtype=bool, na_value=False)
    rows = [dict(zip(bridge._FIELDS, [str(i), 'opaque', str(i), str(int(selected[i]))]
                + [f'{value[i]:016x}' for value in bits])) for i in range(len(frame))]
    listed = [value.strip() for value in tissues.split(',') if value.strip()] if tissues else None
    first = f"PEXT_{listed[0].replace('-', '_')}" if listed else None
    names = list(frame.columns)
    meta = [('axis', 'caller-owned'), ('applied', 'true'), ('column', column),
        ('used_tissue', str(first is not None and first in names).lower()),
        ('tissue_list_present', str(listed is not None).lower()), ('tissue_count', str(len(listed or []))),
        ('first_tissue_column_present', str(first is not None).lower()), ('first_tissue_column', first or ''),
        ('source_column_position', str(names.index(column))),
        ('target_column_position', str(names.index('pext_sort_index') if 'pext_sort_index' in names else len(names))),
        ('score_column_position', str(names.index('sort_index')))]
    meta += [(f'tissue.{i}', value) for i, value in enumerate(listed or [])]
    if change:
        change(rows, meta)
    with Path(command[command.index('--output') + 1]).open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(bridge._FIELDS)
        writer.writerows([row.get(name) for name in bridge._FIELDS] + row.get('_extra', []) for row in rows)
    with Path(command[command.index('--metadata-output') + 1]).open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['name', 'value']); writer.writerows(meta)
    return SimpleNamespace(returncode=0)


def publish(frame, penalty, plan):
    frame['pext_sort_index'] = plan['pext_sort_index']
    frame.loc[penalty, 'sort_index'] = bridge.pext_score_product(frame, penalty, plan)
    return frame


def rank_frame():
    return pd.DataFrame({'BIOTYPE': ['protein_coding'] * 4, 'ACMG_quant_score': [.95, .2, .95, .8],
        'var_plausible_patho_mechs': ['dominant_LOF'] * 4,
        'variant_effect': ['predicted_LOF_high_confidence'] * 4, 'variant_lof_score': [1] * 4,
        'chrom': ['1'] * 4, 'pos': [1, 2, 3, 4], 'ref': ['A'] * 4, 'alt': ['G'] * 4,
        'SYMBOL': ['A', 'B', 'C', 'A'], 'variant_id': ['a', 'b', 'c', 'd'], 'PEXT_MEAN': [.02] * 4})


@contextmanager
def preceding_fallbacks():
    with ExitStack() as stack:
        for module, name in ((setup, 'try_rank_setup'), (pedigree, 'try_pedigree_context'),
                (multiply, 'try_control_multiply'), (compatibility, 'try_ranking_compatibility'),
                (weights, 'try_gene_weights'), (terminal, 'try_terminal_order')):
            stack.enter_context(patch.object(module, name, return_value=None))
        yield


class PextBridgeTests(unittest.TestCase):
    def planned(self, frame, penalty, column='PEXT_MEAN', tissues='', cutoff=.1, floor=.8, shape=.5, change=None):
        with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs:
                response(command, frame, penalty, column, tissues, cutoff, floor, shape, change)):
            return bridge.try_pext(frame, column, tissues, cutoff, floor, shape, penalty, 'native')

    def test_default_and_terminal_only_do_not_enter_pext_bridge(self):
        with patch.object(bridge, 'try_pext') as native, patch.object(terminal, 'try_terminal_order', return_value=None):
            wanted = acmg_scoring.sort_and_rank_variants(rank_frame(), None, None, dispensable_gene_list=None)
            got = acmg_scoring.sort_and_rank_variants(rank_frame(), None, None,
                dispensable_gene_list=None, terminal_order_binary='tail')
        native.assert_not_called(); pd.testing.assert_frame_equal(got, wanted, check_exact=True)

    def test_atomic_plan_nullable_bits_series_metadata_and_two_assignments(self):
        for nullable in (False, True):
            data = entry(); data['PEXT_MEAN'].attrs = {'source': ['cached column metadata']}
            penalty = mask(data, nullable)
            before = pickle.dumps((data, penalty), protocol=5)
            expected_data = entry(); expected_data['PEXT_MEAN'].attrs = {'source': ['cached column metadata']}
            wanted = original(expected_data, mask(expected_data, nullable))
            plan = self.planned(data, penalty); self.assertIsNotNone(plan)
            self.assertEqual(pickle.dumps((data, penalty), protocol=5), before)
            self.assertIs(plan['pext_numeric'].name, data['PEXT_MEAN'].name)
            for name in ('pext_numeric', 'scaled_pext'):
                pd.testing.assert_series_equal(plan[name], wanted[name], check_exact=True)
                self.assertEqual(plan[name].attrs, wanted[name].attrs)
            actual = publish(data, penalty, plan)
            pd.testing.assert_frame_equal(actual, wanted['df'], check_exact=True)
            np.testing.assert_array_equal(actual['sort_index'].to_numpy().view(np.uint64),
                                          wanted['df']['sort_index'].to_numpy().view(np.uint64))

    def test_first_tissue_only_then_mean_and_existing_target_position(self):
        for tissues, column in [('Brain-Cortex,Liver', 'PEXT_Brain_Cortex'), ('Missing,Liver', 'PEXT_MEAN'),
                                (' , , ', 'PEXT_MEAN'), ('\u001cLiver\u001f,Brain-Cortex', 'PEXT_Liver')]:
            data = entry(); data['PEXT_Brain_Cortex'] = [.01] * 4; data['PEXT_Liver'] = [.03] * 4
            data.insert(1, 'pext_sort_index', ['previous'] * 4); penalty = mask(data)
            wanted = original(*pickle.loads(pickle.dumps((data, penalty), protocol=5)), column)
            plan = self.planned(data, penalty, column, tissues); self.assertIsNotNone(plan)
            pd.testing.assert_frame_equal(publish(data, penalty, plan), wanted['df'], check_exact=True)

    def test_transport_uses_actual_public_series_map_boxing(self):
        for dtype, values, kinds in [
                ('object', [np.bool_(False), False, np.float32(.02), np.int64(-1)], ['numpy_bool', 'bool', 'float32', 'int']),
                ('bool', [True, False, True, False], ['bool'] * 4),
                ('boolean', [True, False, pd.NA, True], ['bool', 'bool', 'pd_na', 'bool']),
                ('float32', [.02, np.nan, -0., .1], ['float'] * 4),
                ('Float64', [.02, pd.NA, -0., .1], ['float'] * 4),
                ('string', ['.02', pd.NA, '-0', '.1'], ['str', 'pd_na', 'str', 'str'])]:
            data = entry(); data['PEXT_MEAN'] = pd.Series(values, index=data.index, dtype=dtype); penalty = mask(data)
            def native(command, **kwargs):
                with Path(command[command.index('--coercion-input') + 1]).open(newline='') as stream:
                    rows = list(csv.DictReader(stream))
                self.assertEqual([row['value_kind'] for row in rows], kinds)
                self.assertEqual([row['source_position'] for row in rows], list(map(str, range(4))))
                return response(command, data, penalty)
            with self.subTest(dtype=dtype), patch.object(bridge.subprocess, 'run', side_effect=native):
                self.assertIsNotNone(bridge.try_pext(data, 'PEXT_MEAN', '', .1, .8, .5, penalty, 'native'))

    def test_raw_nullable_mask_and_unconsumed_payload_are_not_coerced(self):
        class Payload:
            def __str__(self): raise AssertionError('unconsumed payload str')
            def __repr__(self): raise AssertionError('unconsumed payload repr')
        data = entry(); penalty = mask(data, True); clean = entry(); clean_mask = mask(clean, True)
        data['payload'] = [Payload()] * 4
        def native(command, **kwargs):
            with Path(command[command.index('--mask-input') + 1]).open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([(row['raw_data'], row['is_null']) for row in rows], [('1', '0'), ('1', '1'), ('0', '0'), ('1', '1')])
            with Path(command[command.index('--input') + 1]).open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([(row['payload_kind'], row['payload']) for row in rows], [('null', '')] * 4)
            return response(command, clean, clean_mask)
        with patch.object(bridge.subprocess, 'run', side_effect=native):
            self.assertIsNotNone(bridge.try_pext(data, 'PEXT_MEAN', '', .1, .8, .5, penalty, 'native'))

    def test_duplicate_and_gapped_axes_keep_original_loc_metadata(self):
        for index in (pd.Index([7, 7, -2, 7], name='duplicate'), pd.RangeIndex(9, 1, -2, name='reverse')):
            data = entry(); data.flags.allows_duplicate_labels = True; data.index = index; penalty = mask(data)
            wanted = original(*pickle.loads(pickle.dumps((data, penalty), protocol=5)))
            plan = self.planned(data, penalty); self.assertIsNotNone(plan)
            pd.testing.assert_frame_equal(publish(data, penalty, plan), wanted['df'], check_exact=True)

    def test_unsupported_options_axes_dtypes_and_metadata_keep_python(self):
        cases = []
        for name in ('sort_index', 'PEXT_MEAN'):
            data = entry(); data[name] = data[name].astype('category'); cases.append((data, mask(data)))
            data = entry().drop(columns=name); cases.append((data, mask(data)))
        data = entry(); data.index = pd.MultiIndex.from_tuples([('a', i) for i in range(4)]); cases.append((data, mask(data)))
        data = entry(); penalty = mask(data); penalty.index = data.index.copy(); cases.append((data, penalty))
        data = entry(); cases.append((data, mask(data).astype(object)))
        data = entry(); data.attrs = {'custom': object()}; cases.append((data, mask(data)))
        data = entry(); data['PEXT_MEAN'].attrs = {'custom': object()}; cases.append((data, mask(data)))
        data = entry().iloc[:0]; cases.append((data, mask(entry()).iloc[:0]))
        data = entry(); data.flags.allows_duplicate_labels = True
        data = pd.concat([data, data[['PEXT_MEAN']]], axis=1); cases.append((data, mask(data)))
        for data, penalty in cases:
            with patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_pext(data, 'PEXT_MEAN', '', .1, .8, .5, penalty, 'native')); native.assert_not_called()
        for value in (None, 1, np.float64(.1), np.nan, np.inf, -1.):
            data = entry(); penalty = mask(data)
            with patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_pext(data, 'PEXT_MEAN', '', value, .8, .5, penalty, 'native')); native.assert_not_called()
        for module in (pd, np):
            data = entry(); penalty = mask(data)
            with patch.object(module, '__version__', 'unqualified'), patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_pext(data, 'PEXT_MEAN', '', .1, .8, .5, penalty, 'native')); native.assert_not_called()

    def test_custom_objects_and_dtype_getters_are_rejected_without_callbacks(self):
        callbacks = []
        class Scalar:
            def __str__(self): callbacks.append('str'); raise AssertionError('custom scalar')
        class Dtype(pd.StringDtype):
            def __str__(self): callbacks.append('dtype str'); return super().__str__()
        class Array(pd.arrays.StringArray):
            @property
            def dtype(self): callbacks.append('dtype getter'); return Dtype()
        for values in ([Scalar()] * 4, pd.Series(Array(np.array(['.1'] * 4, dtype=object)), index=entry().index)):
            data = entry(); data['PEXT_MEAN'] = values; penalty = mask(data); callbacks.clear()
            with patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_pext(data, 'PEXT_MEAN', '', .1, .8, .5, penalty, 'native'))
            native.assert_not_called(); self.assertEqual(callbacks, [])

    def test_malformed_rows_masks_bits_and_metadata_reject_atomically(self):
        changes = [lambda rows, meta: rows.reverse(), lambda rows, meta: rows.pop(),
            lambda rows, meta: rows[0].update(source_position='1'), lambda rows, meta: rows[0].update(index_value='wrong'),
            lambda rows, meta: rows[0].update(selected='2'), lambda rows, meta: rows[0].update(pext_numeric_bits='nan'),
            lambda rows, meta: rows[1].update(score_bits='0000000000000000'),
            lambda rows, meta: rows[0].update(_extra=['unexpected']), lambda rows, meta: meta.reverse(),
            lambda rows, meta: meta.append(meta[-1]), lambda rows, meta: meta.__setitem__(1, ('applied', 'false')),
            lambda rows, meta: meta.__setitem__(9, ('target_column_position', '0'))]
        for change in changes:
            data = entry(); penalty = mask(data); before = pickle.dumps((data, penalty), protocol=5)
            self.assertIsNone(self.planned(data, penalty, change=change))
            self.assertEqual(pickle.dumps((data, penalty), protocol=5), before)

    def test_transport_errors_nonzero_and_resource_failures_preserve_atomicity(self):
        data = entry(); penalty = mask(data); before = pickle.dumps((data, penalty), protocol=5)
        for error in (FileNotFoundError('missing'), subprocess.TimeoutExpired('native', 300), RuntimeError('transport')):
            with patch.object(bridge.subprocess, 'run', side_effect=error):
                self.assertIsNone(bridge.try_pext(data, 'PEXT_MEAN', '', .1, .8, .5, penalty, 'native'))
        with patch.object(bridge.subprocess, 'run', return_value=SimpleNamespace(returncode=1)):
            self.assertIsNone(bridge.try_pext(data, 'PEXT_MEAN', '', .1, .8, .5, penalty, 'native'))
        for error in (MemoryError, RecursionError):
            with patch.object(bridge.subprocess, 'run', side_effect=error), self.assertRaises(error):
                bridge.try_pext(data, 'PEXT_MEAN', '', .1, .8, .5, penalty, 'native')
        self.assertEqual(pickle.dumps((data, penalty), protocol=5), before)

    def test_no_selected_column_skips_options_validation_and_adapter(self):
        data = rank_frame().drop(columns='PEXT_MEAN')
        with preceding_fallbacks(), patch.object(bridge, 'try_pext') as native:
            acmg_scoring.sort_and_rank_variants(data, None, None, dispensable_gene_list=None,
                pext_low_expression_cutoff=None, pext_penalty_floor='invalid', pext_penalty_shape=-1., ranking_binary='native')
        native.assert_not_called()

    def test_original_selection_log_precedes_conditional_validation_failure(self):
        events = []
        with preceding_fallbacks(), patch.object(bridge, 'try_pext') as native, patch.object(acmg_scoring.logger, 'info',
                side_effect=lambda message, *args: events.append(message)):
            with self.assertRaisesRegex(ValueError, '^pext_low_expression_cutoff must be > 0$'):
                acmg_scoring.sort_and_rank_variants(rank_frame(), None, None,
                    dispensable_gene_list=None, pext_low_expression_cutoff=0., ranking_binary='native')
        native.assert_not_called(); self.assertEqual(events[-1], 'Using PEXT_MEAN column for expression-aware sorting')

    def test_native_publication_uses_high_acmg_mask_and_two_assignments_once(self):
        wanted = acmg_scoring.sort_and_rank_variants(rank_frame(), None, None, dispensable_gene_list=None)
        actual_try, original_control = bridge.try_pext, acmg_scoring.compute_control_common
        events, current = [], []
        def control(*args):
            events.append('C1'); result = original_control(*args); current.append(result); return result
        def native(data, column, tissues, cutoff, floor, shape, penalty, binary):
            self.assertIs(data, current[0]); np.testing.assert_array_equal(penalty.to_numpy(), [True, False, True, False])
            self.assertEqual(data['zygosity_inheritance_mechanism_compatibility'].tolist(), [1.] * 4)
            events.append('R7-plan')
            with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs:
                    response(command, data, penalty, column, tissues, cutoff, floor, shape)):
                plan = actual_try(data, column, tissues, cutoff, floor, shape, penalty, binary)
            self.assertIsNotNone(plan); return plan
        cls = type(rank_frame().loc); original_loc = cls.__setitem__; original_set = pd.DataFrame.__setitem__
        def assign(indexer, key, value):
            if events and events[-1] == 'target' and isinstance(key, tuple) and key[1] == 'sort_index': events.append('score')
            return original_loc(indexer, key, value)
        def target(frame, key, value):
            if events and events[-1] == 'R7-plan' and key == 'pext_sort_index': events.append('target')
            return original_set(frame, key, value)
        with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(pedigree, 'try_pedigree_context', return_value=None), patch.object(
                multiply, 'try_control_multiply', return_value=None), patch.object(compatibility, 'try_ranking_compatibility', return_value=None), patch.object(
                weights, 'try_gene_weights', return_value=None), patch.object(terminal, 'try_terminal_order', return_value=None), patch.object(
                bridge, 'try_pext', side_effect=native), patch.object(acmg_scoring, 'compute_control_common', side_effect=control), patch.object(
                cls, '__setitem__', assign), patch.object(pd.DataFrame, '__setitem__', target):
            actual = acmg_scoring.sort_and_rank_variants(rank_frame(), None, None, dispensable_gene_list=None, ranking_binary='native')
        self.assertEqual(events, ['C1', 'R7-plan', 'target', 'score'])
        pd.testing.assert_frame_equal(actual, wanted, check_exact=True)

    def test_second_assignment_failure_keeps_first_target_and_current_score(self):
        outcomes = []
        for enabled in (False, True):
            current = []; actual_try = bridge.try_pext; original_control = acmg_scoring.compute_control_common
            cls = type(rank_frame().loc); original_assign = cls.__setitem__
            def control(*args):
                result = original_control(*args); current.append(result); return result
            def native(data, column, tissues, cutoff, floor, shape, penalty, binary):
                with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs:
                        response(command, data, penalty, column, tissues, cutoff, floor, shape)):
                    return actual_try(data, column, tissues, cutoff, floor, shape, penalty, binary)
            def fail(indexer, key, value):
                if current and indexer.obj is current[0] and 'pext_sort_index' in indexer.obj and isinstance(key, tuple) and key[1] == 'sort_index':
                    raise RuntimeError('original assignment failure')
                return original_assign(indexer, key, value)
            with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(pedigree, 'try_pedigree_context', return_value=None), patch.object(
                    multiply, 'try_control_multiply', return_value=None), patch.object(compatibility, 'try_ranking_compatibility', return_value=None), patch.object(
                    weights, 'try_gene_weights', return_value=None), patch.object(terminal, 'try_terminal_order', return_value=None), patch.object(
                    bridge, 'try_pext', side_effect=native), patch.object(acmg_scoring, 'compute_control_common', side_effect=control), patch.object(cls, '__setitem__', fail):
                with self.assertRaisesRegex(RuntimeError, '^original assignment failure$'):
                    acmg_scoring.sort_and_rank_variants(rank_frame(), None, None, dispensable_gene_list=None,
                        **({'ranking_binary': 'native'} if enabled else {}))
            self.assertIn('pext_sort_index', current[0]); outcomes.append(current[0])
        pd.testing.assert_frame_equal(outcomes[0], outcomes[1], check_exact=True)


if __name__ == '__main__':
    unittest.main()
