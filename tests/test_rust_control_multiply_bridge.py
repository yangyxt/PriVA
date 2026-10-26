"""R3 atomic product, current-buffer fallback and original control ordering."""
import csv
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
import priva_rust_control_multiply as bridge
import priva_rust_pedigree_context as pedigree
import priva_rust_rank_setup as setup
import priva_rust_terminal_order as terminal


class OpaquePayload:
    def __str__(self):
        raise AssertionError('unconsumed payload callback')
    def __repr__(self):
        raise AssertionError('unconsumed payload callback')
    def __deepcopy__(self, _memo):
        raise AssertionError('unconsumed payload callback')


def entry():
    result = pd.DataFrame({'payload': [OpaquePayload(), OpaquePayload()],
                           'sort_index': [1., -0.], 'control_common_index': [.5, .5]})
    result.index = pd.Index(['second', 'first'], name='opaque rows')
    result.columns.name = 'original columns'
    result.attrs = {'nested': {'scope': ['current control return']}}
    result.flags.allows_duplicate_labels = False
    return result


def frame():
    return pd.DataFrame(dict(ACMG_quant_score=[.96, .7], BIOTYPE=['protein_coding'] * 2,
        chrom=['1'] * 2, pos=[1, 2], ref=['A'] * 2, alt=['G'] * 2, variant_id=['a', 'b'],
        SYMBOL=['A'] * 2, var_plausible_patho_mechs=[''] * 2, variant_effect=['uncertain'] * 2,
        control_AN=[100., 100.], control_AC=[30., 0.], control_AF=[.3, 0.], control_nhomalt=[0., 0.]))


def rank(data, **kwargs):
    return acmg_scoring.sort_and_rank_variants(data, None, None, dispensable_gene_list=None, **kwargs)


def response(command, *, bits=('3fe0000000000000', '8000000000000000'), target=1, change=None):
    fields = ['source_position', 'index_kind', 'index_value', 'float_bits']
    rows = [dict(source_position=str(i), index_kind='opaque', index_value=str(i), float_bits=value)
            for i, value in enumerate(bits)]
    meta = [('target_source_column_position', str(target)), ('dtype', 'float64'), ('axis', 'caller-owned')]
    if change:
        change(rows, meta)
    with Path(command[command.index('--output') + 1]).open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(fields)
        writer.writerows([row.get(field) for field in fields] + row.get('_extra', []) for row in rows)
    with Path(command[command.index('--metadata-output') + 1]).open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['name', 'value']); writer.writerows(meta)
    return SimpleNamespace(returncode=0)


class ControlMultiplyBridgeTests(unittest.TestCase):
    def test_default_and_r8_only_never_enter_r3_adapter(self):
        with patch.object(bridge, 'try_control_multiply') as product, patch.object(
                terminal, 'try_terminal_order', return_value=None):
            wanted = rank(frame())
            pd.testing.assert_frame_equal(rank(frame(), terminal_order_binary='tail'), wanted, check_exact=True)
        product.assert_not_called()

    def test_fallback_uses_current_control_return_once_with_same_frame(self):
        data = frame(); before = pickle.dumps(data, protocol=5)
        expected = rank(pickle.loads(before))
        events, current = [], []
        original_control, original_masks = acmg_scoring.compute_control_common, acmg_scoring._variant_mechanism_masks
        def control(*args):
            events.append('C1'); result = original_control(*args); current.append(result); return result
        def product(actual, binary):
            events.append('R3'); self.assertIs(actual, current[0]); self.assertEqual(binary, 'native')
            self.assertIn('control_common_index', actual)
            return None
        def masks(actual):
            events.append('R4'); self.assertIs(actual, current[0]); return original_masks(actual)
        with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(
                pedigree, 'try_pedigree_context', return_value=None), patch.object(
                bridge, 'try_control_multiply', side_effect=product), patch.object(
                acmg_scoring, 'compute_control_common', side_effect=control), patch.object(
                acmg_scoring, '_variant_mechanism_masks', side_effect=masks), patch.object(
                terminal, 'try_terminal_order', return_value=None):
            actual = rank(data, ranking_binary='native')
        pd.testing.assert_frame_equal(actual, expected, check_exact=True)
        self.assertEqual(events, ['C1', 'R3', 'R4'])
        self.assertEqual(pickle.dumps(data, protocol=5), before)

    def test_product_is_atomic_and_preserves_bits_axes_attrs_flags_and_payload(self):
        data = entry(); before = pickle.dumps(data, protocol=5)
        seen = []
        def native(command, **kwargs):
            self.assertEqual(command[1], 'control-multiply')
            self.assertEqual(command[command.index('--row-index-kind') + 1], 'flat')
            with Path(command[command.index('--input') + 1]).open(newline='') as stream:
                records = list(csv.DictReader(stream))
            self.assertEqual([(row['payload_kind'], row['payload']) for row in records], [('null', '')] * 2)
            seen.append(kwargs)
            return response(command)
        with patch.object(bridge.subprocess, 'run', side_effect=native):
            product = bridge.try_control_multiply(data, 'native')
        self.assertIsNotNone(product)
        self.assertEqual(product.to_numpy(copy=False).view(np.uint64).tolist(), [0x3fe0000000000000, 0x8000000000000000])
        self.assertIs(product.index, data.index); self.assertIsNone(product.name)
        self.assertEqual(product.attrs, data.attrs); self.assertIsNot(product.attrs, data.attrs)
        self.assertEqual(product.flags.allows_duplicate_labels, data.flags.allows_duplicate_labels)
        self.assertEqual(pickle.dumps(data, protocol=5), before)
        self.assertEqual((seen[0]['timeout'], seen[0]['check']), (300, False))

    def test_accepted_product_keeps_post_control_frame_and_original_assignment(self):
        data = frame().drop(columns=['control_AN', 'control_AC', 'control_AF', 'control_nhomalt'])
        wanted = rank(data)
        current, events = [], []
        original_control, original_masks = acmg_scoring.compute_control_common, acmg_scoring._variant_mechanism_masks
        def control(*args):
            events.append('C1'); result = original_control(*args); current.append(result); return result
        def product(actual, _binary):
            events.append('R3'); self.assertIs(actual, current[0])
            return pd.Series([.96, .7], index=actual.index)
        def masks(actual):
            events.append('R4'); self.assertIs(actual, current[0])
            self.assertEqual(actual['sort_index'].tolist(), [.96, .7])
            return original_masks(actual)
        with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(
                pedigree, 'try_pedigree_context', return_value=None), patch.object(
                bridge, 'try_control_multiply', side_effect=product), patch.object(
                acmg_scoring, 'compute_control_common', side_effect=control), patch.object(
                acmg_scoring, '_variant_mechanism_masks', side_effect=masks), patch.object(
                terminal, 'try_terminal_order', return_value=None):
            actual = rank(data, ranking_binary='native')
        pd.testing.assert_frame_equal(actual, wanted, check_exact=True)
        self.assertEqual(events, ['C1', 'R3', 'R4'])

    def test_actual_control_failure_precedes_native_attempt(self):
        data = frame().drop(columns='control_AN')
        failures = []
        for native in (False, True):
            events = []
            original = acmg_scoring.compute_control_common
            def control(*args):
                events.append('C1'); return original(*args)
            with self.subTest(native=native), patch.object(setup, 'try_rank_setup', return_value=None), patch.object(
                    pedigree, 'try_pedigree_context', return_value=None), patch.object(
                    acmg_scoring.logger, 'warning', side_effect=lambda *args: events.append('warning')), patch.object(
                    acmg_scoring, 'compute_control_common', side_effect=control), patch.object(
                    bridge, 'try_control_multiply') as product:
                with self.assertRaises(AttributeError) as error:
                    rank(data, **({'ranking_binary': 'native'} if native else {}))
            failures.append((type(error.exception), error.exception.args))
            self.assertEqual(events, ['warning', 'C1']); product.assert_not_called()
        self.assertEqual(failures[0], failures[1])

    def test_missing_post_control_labels_retain_original_left_lookup_failure(self):
        data = frame(); returned = []
        def control(actual, _proband):
            result = actual.drop(columns='sort_index'); returned.append(result); return result
        with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(
                pedigree, 'try_pedigree_context', return_value=None), patch.object(
                acmg_scoring, 'compute_control_common', side_effect=control) as actual_control, patch.object(
                bridge.subprocess, 'run') as native:
            with self.assertRaises(KeyError) as error:
                rank(data, ranking_binary='native')
        self.assertEqual(error.exception.args, ('sort_index',)); self.assertEqual(actual_control.call_count, 1)
        native.assert_not_called(); self.assertNotIn('sort_index', returned[0])

    def test_empty_product_retains_float64_and_metadata(self):
        data = entry().iloc[:0]
        with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs: response(command, bits=())):
            product = bridge.try_control_multiply(data, 'native')
        self.assertIsNotNone(product); self.assertEqual(len(product), 0); self.assertEqual(str(product.dtype), 'float64')
        pd.testing.assert_index_equal(product.index, data.index, exact=True)
        self.assertEqual(product.attrs, data.attrs); self.assertFalse(product.flags.allows_duplicate_labels)

    def test_unsupported_axes_dtypes_and_labels_fall_back_before_native(self):
        cases = []
        data = entry(); data.index = pd.MultiIndex.from_tuples([('a', 1), ('b', 2)]); cases.append(data)
        for operand in ('sort_index', 'control_common_index'):
            for dtype in ('float32', '>f8', 'Int64', 'Float64', 'object'):
                data = entry(); data[operand] = pd.Series([1, 2], index=data.index, dtype=dtype); cases.append(data)
            data = entry().drop(columns=operand); cases.append(data)
            data = entry(); data.flags.allows_duplicate_labels = True
            cases.append(pd.concat([data, data[[operand]]], axis=1))
        for data in cases:
            with self.subTest(columns=list(data.columns)), patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_control_multiply(data, 'native')); native.assert_not_called()

    def test_two_nan_payload_route_is_an_atomic_fallback(self):
        data = entry()
        data['sort_index'] = np.asarray([0x7ff0000000000123, 0x8000000000000000], dtype=np.uint64).view(np.float64)
        data['control_common_index'] = np.asarray([0xfff8000000000456, 0x3ff0000000000000], dtype=np.uint64).view(np.float64)
        before = pickle.dumps(data, protocol=5)
        with patch.object(bridge.subprocess, 'run') as native:
            self.assertIsNone(bridge.try_control_multiply(data, 'native')); native.assert_not_called()
        self.assertEqual(pickle.dumps(data, protocol=5), before)

    def test_real_large_frame_retains_original_numerical_execution_route(self):
        data = pd.DataFrame({'sort_index': np.zeros(1_000_001), 'control_common_index': np.ones(1_000_001)})
        with patch.object(bridge.subprocess, 'run') as native:
            self.assertIsNone(bridge.try_control_multiply(data, 'native')); native.assert_not_called()

    def test_custom_metadata_and_unqualified_versions_do_not_run_callbacks(self):
        data = entry(); data.attrs['custom'] = OpaquePayload()
        with patch.object(bridge.subprocess, 'run') as native:
            self.assertIsNone(bridge.try_control_multiply(data, 'native')); native.assert_not_called()
        for library in (pd, np):
            with patch.object(library, '__version__', 'unqualified'), patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_control_multiply(entry(), 'native')); native.assert_not_called()

    def test_nonzero_timeout_and_os_failure_leave_current_frame_unchanged(self):
        data = entry(); before = pickle.dumps(data, protocol=5)
        for outcome in (SimpleNamespace(returncode=1), subprocess.TimeoutExpired('native', 300), FileNotFoundError('native')):
            with self.subTest(outcome=type(outcome).__name__), patch.object(bridge.subprocess, 'run') as native:
                if isinstance(outcome, Exception): native.side_effect = outcome
                else: native.return_value = outcome
                self.assertIsNone(bridge.try_control_multiply(data, 'native'))
            self.assertEqual(pickle.dumps(data, protocol=5), before)

    def test_malformed_product_or_metadata_never_publishes(self):
        data = entry(); before = pickle.dumps(data, protocol=5)
        changes = [lambda rows, meta: rows.reverse(), lambda rows, meta: rows.pop(),
                   lambda rows, meta: rows[0].update(index_value='wrong'),
                   lambda rows, meta: rows[0].update(float_bits='not-bits'),
                   lambda rows, meta: rows[0].update(_extra=['ignored value']),
                   lambda rows, meta: meta.__setitem__(0, ('target_source_column_position', '0')),
                   lambda rows, meta: meta.__setitem__(1, ('dtype', 'Float64')),
                   lambda rows, meta: meta.__setitem__(2, ('axis', 'wrong')),
                   lambda rows, meta: meta.append(meta[0])]
        for change in changes:
            with self.subTest(change=changes.index(change)), patch.object(bridge.subprocess, 'run',
                    side_effect=lambda command, **kwargs: response(command, change=change)):
                self.assertIsNone(bridge.try_control_multiply(data, 'native'))
            self.assertEqual(pickle.dumps(data, protocol=5), before)

    def test_product_construction_failure_is_atomic(self):
        data = entry(); before = pickle.dumps(data, protocol=5)
        with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs: response(command)), patch.object(
                bridge.pd, 'Series', side_effect=RuntimeError('private product construction failed')):
            self.assertIsNone(bridge.try_control_multiply(data, 'native'))
        self.assertEqual(pickle.dumps(data, protocol=5), before)

    def test_existing_numerical_environment_is_preserved(self):
        seen = []
        def native(command, **kwargs):
            seen.append(kwargs['env']['POLARS_MAX_THREADS']); return response(command)
        with patch.dict(bridge.os.environ, {'POLARS_MAX_THREADS': '3'}), patch.object(
                bridge.subprocess, 'run', side_effect=native):
            self.assertIsNotNone(bridge.try_control_multiply(entry(), 'native'))
        self.assertEqual(seen, ['3'])

    def test_process_resource_errors_propagate(self):
        for error in (MemoryError('resource'), RecursionError('depth')):
            with self.subTest(error=type(error).__name__), patch.object(bridge.subprocess, 'run', side_effect=error):
                with self.assertRaises(type(error)):
                    bridge.try_control_multiply(entry(), 'native')


if __name__ == '__main__':
    unittest.main()
