"""R1 dispatch, source-order and atomic-fallback controls without a native build."""
import copy
import csv
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
import priva_rust_rank_setup as bridge
import priva_rust_terminal_order as terminal


def frame():
    return pd.DataFrame(dict(ACMG_quant_score=[.96, .94, .8],
        BIOTYPE=['protein_coding', 'other', 'protein_coding'],
        chrom=['2', '1', '1'], pos=[20, 10, 9], ref=['A', 'G', 'C'], alt=['G', 'A', 'T'],
        variant_id=['second', 'excluded', 'first'], SYMBOL=['A', 'B', 'C'],
        var_plausible_patho_mechs=['', '', ''], variant_effect=['uncertain'] * 3))


def rank(data, **kwargs):
    return acmg_scoring.sort_and_rank_variants(data, None, None,
        dispensable_gene_list=None, **kwargs)


def native_response(command, data, positions=(0, 2)):
    with Path(command[command.index('--output') + 1]).open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['source_position', 'index_kind', 'index_value'])
        writer.writerows([position, 'opaque', position] for position in positions)
    names, dtypes = list(data.columns), [str(dtype) for dtype in data.dtypes]
    columns, drops = bridge._column_plan(names, dtypes)
    with Path(command[command.index('--metadata-output') + 1]).open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            'record', 'position', 'name', 'dtype', 'source_column', 'operation'])
        writer.writeheader()
        writer.writerows(columns + drops + [dict(record='axis', position='',
            name='caller-owned', dtype='opaque', source_column='', operation='')])
    return SimpleNamespace(returncode=0)


class RankSetupBridgeTests(unittest.TestCase):
    def test_default_enters_neither_native_adapter(self):
        with patch.object(bridge, 'try_rank_setup') as setup, patch.object(terminal, 'try_terminal_order') as tail:
            self.assertEqual(rank(frame())['variant_rank'].tolist(), [1, 2])
        setup.assert_not_called(); tail.assert_not_called()

    def test_fallback_sees_unmodified_entry_frame_before_control(self):
        data = frame()
        data['sort_index'] = 'old'
        data['haplo_insuf_index'] = 'old'
        before = pickle.dumps(data, protocol=5)
        seen = []
        def fallback(current, binary):
            self.assertIs(current, data)
            self.assertEqual(binary, 'one-native-binary')
            self.assertEqual(pickle.dumps(current, protocol=5), before)
            self.assertNotIn('control_common_index', current)
            seen.append(True)
            return None
        with patch.object(bridge, 'try_rank_setup', fallback), patch.object(terminal, 'try_terminal_order', return_value=None):
            actual = rank(data, ranking_binary='one-native-binary')
        wanted = rank(pickle.loads(before))
        pd.testing.assert_frame_equal(actual, wanted, check_exact=True)
        self.assertEqual(pickle.dumps(data, protocol=5), before)
        self.assertEqual(seen, [True])

    def test_native_plan_gathers_axes_values_and_original_mutation_order(self):
        data = frame()
        data.index = pd.Index(['same', None, 'same'], dtype=object, name='labels')
        data.columns.name = 'columns'
        data.attrs = {'nested': {'values': [1, 'x']}}
        data['haplo_insuf_index'] = 'old'
        data.insert(1, 'zygosity_inheritance_mechanism_index', 'old')
        data.insert(3, 'sort_index', 'old')
        before = pickle.dumps(data, protocol=5)
        with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs: native_response(command, data)):
            actual = bridge.try_rank_setup(data, 'native')
        expected = data.copy()
        for name in bridge._STALE:
            expected = expected.drop(columns=[name])
        expected['sort_index'] = expected.loc[:, 'ACMG_quant_score']
        expected = expected.loc[expected['BIOTYPE'] == 'protein_coding', :]
        pd.testing.assert_frame_equal(actual, expected, check_exact=True)
        self.assertEqual(actual.attrs, expected.attrs)
        self.assertEqual(pickle.dumps(data, protocol=5), before)

    def test_original_entry_failures_reexecute_original_source(self):
        for missing in ('ACMG_quant_score', 'BIOTYPE'):
            data = frame().drop(columns=[missing])
            with self.subTest(missing=missing):
                with self.assertRaises(KeyError) as original:
                    rank(data)
                with patch.object(bridge.subprocess, 'run') as native:
                    with self.assertRaises(type(original.exception)) as actual:
                        rank(data, ranking_binary='native')
                self.assertEqual(actual.exception.args, original.exception.args)
                native.assert_not_called()

    def test_native_nonzero_timeout_and_bad_rows_preserve_input(self):
        data = frame()
        for outcome in (SimpleNamespace(returncode=1), subprocess.TimeoutExpired('native', 300),
                        FileNotFoundError('native')):
            before = pickle.dumps(data, protocol=5)
            with self.subTest(outcome=type(outcome).__name__), patch.object(bridge.subprocess, 'run') as native:
                if isinstance(outcome, Exception):
                    native.side_effect = outcome
                else:
                    native.return_value = outcome
                self.assertIsNone(bridge.try_rank_setup(data, 'native'))
            self.assertEqual(pickle.dumps(data, protocol=5), before)
        for positions in ((2, 0), (0, 0), (-1,), (len(data),)):
            with self.subTest(positions=positions), patch.object(bridge.subprocess, 'run',
                    side_effect=lambda command, **kwargs: native_response(command, data, positions)):
                self.assertIsNone(bridge.try_rank_setup(data, 'native'))
            self.assertEqual(pickle.dumps(data, protocol=5), before)

    def test_bad_plan_and_failure_after_private_copy_fall_back_atomically(self):
        data = frame()
        before = pickle.dumps(data, protocol=5)
        def bad(command, **kwargs):
            result = native_response(command, data)
            path = Path(command[command.index('--metadata-output') + 1])
            path.write_text(path.read_text().replace('copy_score', 'retain'))
            return result
        with patch.object(bridge.subprocess, 'run', bad):
            self.assertIsNone(bridge.try_rank_setup(data, 'native'))
        self.assertEqual(pickle.dumps(data, protocol=5), before)
        with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs: native_response(command, data)), \
                patch.object(pd.DataFrame, 'drop', side_effect=ValueError('private gather failed')):
            data['haplo_insuf_index'] = 'old'
            before = pickle.dumps(data, protocol=5)
            self.assertIsNone(bridge.try_rank_setup(data, 'native'))
        self.assertEqual(pickle.dumps(data, protocol=5), before)

    def test_unsupported_dtypes_and_custom_metadata_do_not_enter_native(self):
        examples = []
        data = frame(); data['BIOTYPE'] = data['BIOTYPE'].astype('category'); examples.append(data)
        data = frame(); data['payload'] = [[1], [2], [3]]; examples.append(data)
        class Callback:
            def __deepcopy__(self, memo):
                raise AssertionError('custom metadata callback ran')
        data = frame(); data.attrs['custom'] = Callback(); examples.append(data)
        for data in examples:
            with patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_rank_setup(data, 'native'))
            native.assert_not_called()

    def test_nullable_types_and_unrelated_duplicate_columns_are_transported(self):
        data = frame()
        data['ACMG_quant_score'] = pd.Series([.96, pd.NA, -.0], dtype='Float64')
        data['unused'] = pd.Series([True, pd.NA, False], dtype='boolean')
        data = pd.concat([data, data[['unused']]], axis=1)
        with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs: native_response(command, data)):
            actual = bridge.try_rank_setup(data, 'native')
        self.assertEqual(actual.columns.tolist().count('unused'), 2)
        self.assertEqual(str(actual['sort_index'].dtype), 'Float64')
        self.assertEqual(np.signbit(actual['sort_index'].iloc[-1]), True)

    def test_unqualified_versions_and_consumed_duplicate_labels_stay_python(self):
        for name in ('BIOTYPE', 'ACMG_quant_score', 'sort_index'):
            data = frame()
            if name not in data:
                data[name] = 'old'
            data = pd.concat([data, data[[name]]], axis=1)
            with patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_rank_setup(data, 'native'))
            native.assert_not_called()
        with patch.object(pd, '__version__', 'unqualified'), patch.object(bridge.subprocess, 'run') as native:
            self.assertIsNone(bridge.try_rank_setup(frame(), 'native'))
        native.assert_not_called()


if __name__ == '__main__':
    unittest.main()
