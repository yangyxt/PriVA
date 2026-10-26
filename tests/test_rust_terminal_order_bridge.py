"""Check opt-in routing and atomic fallback independently of a native build."""
from pathlib import Path
import copy
import csv
import pickle
import struct
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

import acmg_scoring
import priva_rust_terminal_order as bridge


def frame():
    return pd.DataFrame(dict(ACMG_quant_score=[.96, .94], BIOTYPE=['protein_coding'] * 2,
        chrom=['2', '1'], pos=[20, 10], ref=['A', 'G'], alt=['G', 'A'],
        variant_id=['second', 'first'], SYMBOL=['A', 'B'],
        var_plausible_patho_mechs=['', ''], variant_effect=['uncertain', 'uncertain']))


def rank(data, **kwargs):
    return acmg_scoring.sort_and_rank_variants(data, None, None, dispensable_gene_list=None, **kwargs)


def state(data):
    """Snapshot public axes/attrs/flags and every current raw cell buffer."""
    return (type(data.index), data.index.names, data.index.to_numpy().tobytes(),
            type(data.columns), data.columns.names, data.columns.to_numpy().tobytes(),
            copy.deepcopy(data.attrs), data.flags.allows_duplicate_labels,
            [(str(data.iloc[:, i].dtype), data.iloc[:, i].to_numpy().tobytes())
             for i in range(data.shape[1])])


def native_response(command, data, rows=None):
    """Small valid protocol response, used only for injecting malformed fields."""
    if rows is None:
        rows = [[str(i), 'opaque', str(i), 'float', struct.pack('>d', score).hex(), str(i + 1)]
                for i, score in enumerate(data['sort_index'])]
    with Path(command[command.index('--output') + 1]).open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['source_position', 'index_kind', 'index_value', 'max_kind', 'max_value', 'variant_rank'])
        writer.writerows(rows)
    with Path(command[command.index('--metadata-output') + 1]).open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['record', 'position', 'name', 'dtype'])
        for position, (name, dtype) in enumerate(zip(data.columns, data.dtypes)):
            writer.writerow(['column', position, name, str(dtype)])
        writer.writerow(['column', len(data.columns), 'variant_rank', 'int64'])
        writer.writerow(['temporary', len(data.columns), 'max_variant_score', 'float64'])
        writer.writerow(['axis', '', 'caller-owned', 'opaque'])
    return SimpleNamespace(returncode=0)


class TerminalOrderBridgeTests(unittest.TestCase):
    def test_default_ranker_does_not_enter_native_adapter(self):
        with patch.object(bridge, 'try_terminal_order') as native:
            self.assertEqual(rank(frame())['variant_rank'].tolist(), [1, 2])
        native.assert_not_called()

    def test_fallback_receives_current_weighted_frame(self):
        data = frame()
        data['control_AF'] = [.9, 0.]
        data['control_AN'] = [100, 100]
        data['control_AC'] = [95, 0]
        data['control_nhomalt'] = [0, 0]
        seen = []
        def fallback(current, binary):
            self.assertEqual(binary, 'explicit-path')
            seen.append(current)
            self.assertIsNot(current, data)
            self.assertNotEqual(current['sort_index'].tolist(), data['ACMG_quant_score'].tolist())
            self.assertNotIn('max_variant_score', current)
            return None
        wanted = rank(pickle.loads(pickle.dumps(data)))
        with patch.object(bridge, 'try_terminal_order', fallback):
            actual = rank(data, terminal_order_binary='explicit-path')
        pd.testing.assert_frame_equal(actual, wanted, check_exact=True)
        self.assertEqual(len(seen), 1)
        self.assertNotIn('sort_index', data)
        self.assertNotIn('max_variant_score', data)

    def test_failed_process_keeps_entry_unchanged_and_cleans_temporaries(self):
        current = frame()
        current['sort_index'] = current['ACMG_quant_score']
        wanted = current.copy(deep=True)
        calls = []
        def fail(command, **kwargs):
            calls.append(Path(command[command.index('--input') + 1]).parent)
            raise OSError('unavailable executable')
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(bridge.subprocess, 'run', fail), patch.object(bridge.tempfile, 'tempdir', directory):
                self.assertIsNone(bridge.try_terminal_order(current, 'unavailable'))
            self.assertEqual(list(Path(directory).iterdir()), [])
        pd.testing.assert_frame_equal(current, wanted, check_exact=True)
        self.assertEqual(len(calls), 1)
        self.assertFalse(calls[0].exists())

    def test_malformed_native_result_falls_back_to_original_tail(self):
        class Result:
            returncode = 0
        def malformed(command, **kwargs):
            Path(command[command.index('--output') + 1]).write_text('wrong_header\n')
            Path(command[command.index('--metadata-output') + 1]).write_text('wrong_header\n')
            return Result()
        data = frame()
        wanted = rank(pickle.loads(pickle.dumps(data)))
        with patch.object(bridge.subprocess, 'run', malformed):
            actual = rank(data, terminal_order_binary='malformed')
        pd.testing.assert_frame_equal(actual, wanted, check_exact=True)

    def test_original_failure_before_r8_does_not_call_adapter(self):
        data = frame().drop(columns='BIOTYPE')
        with patch.object(bridge, 'try_terminal_order') as native:
            with self.assertRaisesRegex(KeyError, 'BIOTYPE'):
                rank(data, terminal_order_binary='unused')
        native.assert_not_called()

    def test_original_r8_exception_is_reexecuted_after_native_failure(self):
        calls = []
        original = bridge.try_terminal_order
        def observe(current, binary):
            calls.append(current)
            return original(current, binary)
        data = frame().drop(columns='pos')
        with self.assertRaises(KeyError) as expected:
            rank(pickle.loads(pickle.dumps(data)))
        with patch.object(bridge, 'try_terminal_order', observe):
            with self.assertRaises(type(expected.exception)) as actual:
                rank(data, terminal_order_binary='unused')
        self.assertEqual(actual.exception.args, expected.exception.args)
        self.assertEqual(len(calls), 1)
        self.assertNotIn('max_variant_score', calls[0])

    def test_unsupported_consumed_dtypes_stay_in_python(self):
        for dtype in ['category', 'complex128']:
            with self.subTest(dtype=dtype):
                current = frame()
                current['sort_index'] = current['ACMG_quant_score']
                current['chrom'] = current['chrom'].astype(dtype) if dtype == 'category' else pd.Series([1j, 2j])
                with patch.object(bridge.subprocess, 'run') as native:
                    self.assertIsNone(bridge.try_terminal_order(current, 'unused'))
                native.assert_not_called()

    def test_unqualified_runtime_retains_python(self):
        current = frame()
        current['sort_index'] = current['ACMG_quant_score']
        with patch.object(pd, '__version__', 'other'), patch.object(bridge.subprocess, 'run') as native:
            self.assertIsNone(bridge.try_terminal_order(current, 'unused'))
        native.assert_not_called()

    def test_timeout_and_nonzero_exit_leave_original_buffers_and_axes(self):
        for failure in ['timeout', 'nonzero']:
            with self.subTest(failure=failure):
                current = frame()
                current['sort_index'] = current['ACMG_quant_score']
                current.index = pd.Index(['same', 'same'], name='original-index')
                current.columns.name = 'original-columns'
                current.attrs = {'scope': ['unchanged']}
                before = state(current)
                calls = []
                def fail(command, **kwargs):
                    calls.append(kwargs['timeout'])
                    if failure == 'timeout':
                        raise subprocess.TimeoutExpired(command, kwargs['timeout'])
                    return SimpleNamespace(returncode=1)
                with patch.object(bridge.subprocess, 'run', fail):
                    self.assertIsNone(bridge.try_terminal_order(current, 'failed-native'))
                self.assertEqual(calls, [300])
                self.assertEqual(state(current), before)

    def test_invalid_permutations_and_rank_scalars_leave_current_buffers(self):
        for field, value in [('source_position', '0'), ('source_position', '2'),
                             ('variant_rank', 'NaN'), ('variant_rank', '1.5'),
                             ('variant_rank', '-1'), ('variant_rank', str(2**63))]:
            with self.subTest(field=field, value=value):
                current = frame()
                current['sort_index'] = current['ACMG_quant_score']
                before = state(current)
                def malformed(command, **kwargs):
                    rows = [[str(i), 'opaque', str(i), 'float', struct.pack('>d', score).hex(), str(i + 1)]
                            for i, score in enumerate(current['sort_index'])]
                    rows[1][0 if field == 'source_position' else 5] = value
                    return native_response(command, current, rows)
                with patch.object(bridge.subprocess, 'run', malformed):
                    self.assertIsNone(bridge.try_terminal_order(current, 'malformed-native'))
                self.assertEqual(state(current), before)

    def test_failure_after_private_staging_leaves_live_frame_untouched(self):
        current = frame()
        current['sort_index'] = current['ACMG_quant_score']
        before = state(current)
        stages = []
        def staged_failure(data, *args, **kwargs):
            self.assertIsNot(data, current)
            self.assertIn('max_variant_score', data)
            stages.append(True)
            raise ValueError('injected failure after private staging')
        with patch.object(bridge.subprocess, 'run', lambda command, **kwargs: native_response(command, current)):
            with patch.object(pd.DataFrame, 'drop', staged_failure):
                self.assertIsNone(bridge.try_terminal_order(current, 'malformed-native'))
        self.assertEqual(stages, [True])
        self.assertEqual(state(current), before)

    def test_custom_metadata_is_rejected_before_user_copy_callbacks(self):
        callbacks = []
        class Callback:
            def __deepcopy__(self, memo):
                callbacks.append(True)
                raise ValueError('caller-owned callback')
        for scope in ['attrs', 'index-label', 'index-name', 'column-name']:
            with self.subTest(scope=scope):
                current = frame()
                current['sort_index'] = current['ACMG_quant_score']
                value = Callback()
                if scope == 'attrs':
                    current.attrs['custom'] = value
                elif scope == 'index-label':
                    current.index = pd.Index([value, 'other'])
                elif scope == 'index-name':
                    current.index.name = value
                else:
                    current.columns.name = value
                with patch.object(bridge.subprocess, 'run') as native:
                    self.assertIsNone(bridge.try_terminal_order(current, 'unused'))
                native.assert_not_called()
        self.assertEqual(callbacks, [])

    def test_numpy_metadata_subclass_cannot_run_copy_callbacks(self):
        callbacks = []
        class CustomFloat(np.float64):
            def __deepcopy__(self, memo):
                callbacks.append(True)
                raise ValueError('caller-owned scalar callback')
        current = frame()
        current['sort_index'] = current['ACMG_quant_score']
        current.attrs['custom'] = CustomFloat(1.)
        with patch.object(bridge.subprocess, 'run') as native:
            self.assertIsNone(bridge.try_terminal_order(current, 'unused'))
        native.assert_not_called()
        self.assertEqual(callbacks, [])

    def test_consumed_scalar_subclasses_cannot_run_conversion_callbacks(self):
        callbacks = []
        class CustomString(str):
            def encode(self, *args, **kwargs):
                callbacks.append(True)
                raise ValueError('caller-owned encode callback')
        current = frame()
        current['sort_index'] = current['ACMG_quant_score']
        current['chrom'] = pd.Series([CustomString('2'), '1'], dtype=object)
        with patch.object(bridge.subprocess, 'run') as native:
            self.assertIsNone(bridge.try_terminal_order(current, 'unused'))
        native.assert_not_called()
        self.assertEqual(callbacks, [])


if __name__ == '__main__':
    unittest.main()
