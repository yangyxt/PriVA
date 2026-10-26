"""R2 atomic publication, scalar transport and original failure-order controls."""
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
import priva_rust_pedigree_context as bridge
import priva_rust_rank_setup as setup
import priva_rust_terminal_order as terminal


def frame(values=('0/1:10', '0/0:1', '1/1'), dtype='object'):
    n = len(values)
    result = pd.DataFrame(dict(ACMG_quant_score=[.96] * n, BIOTYPE=['protein_coding'] * n,
        chrom=['1'] * n, pos=list(range(n)), ref=['A'] * n, alt=['G'] * n,
        variant_id=[str(i) for i in range(n)], SYMBOL=['A'] * n,
        var_plausible_patho_mechs=[''] * n, variant_effect=['uncertain'] * n))
    result['PROBAND'] = pd.Series(values, dtype=dtype)
    return result


def pedigree(sex='2'):
    return pd.DataFrame({'#FamilyID': ['OTHER', 'F', 'F'], 'IndividualID': ['ABSENT', 'PROBAND', 'SECOND'],
                         'Phenotype': ['2', '2', 2], 'Sex': ['1', sex, '1']})


def response(command, *, values=(0, 1, 0), nulls=None, selected=1, male='0', female='1',
             warning='none', dtype='bool', change=None):
    nulls = [0] * len(values) if nulls is None else nulls
    rows = [dict(source_position=str(i), index_kind='opaque', index_value=str(i),
                 value=str(value), is_null=str(missing)) for i, (value, missing) in enumerate(zip(values, nulls))]
    present = selected is not None
    meta = list(zip(bridge._META_NAMES, ['1' if present else '0', 'PROBAND' if present else '',
        str(selected) if present else '', male, female, dtype, '1' if present else '0', warning, 'caller-owned']))
    if change:
        change(rows, meta)
    with Path(command[command.index('--output') + 1]).open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['source_position', 'index_kind', 'index_value', 'value', 'is_null'])
        writer.writeheader(); writer.writerows(rows)
    with Path(command[command.index('--metadata-output') + 1]).open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['name', 'value']); writer.writerows(meta)
    return SimpleNamespace(returncode=0)


def rank(data, ped=None, family=None, **kwargs):
    return acmg_scoring.sort_and_rank_variants(data, ped, family, dispensable_gene_list=None, **kwargs)


class PedigreeContextBridgeTests(unittest.TestCase):
    def test_default_and_r8_only_never_enter_r2_adapter(self):
        with patch.object(bridge, 'try_pedigree_context') as context, patch.object(
                terminal, 'try_terminal_order', return_value=None):
            expected = rank(frame())
            pd.testing.assert_frame_equal(rank(frame(), terminal_order_binary='tail'), expected, check_exact=True)
        context.assert_not_called()

    def test_fallback_boundary_after_r1_before_original_control(self):
        data, ped = frame(), pedigree()
        data.loc[1, 'BIOTYPE'] = 'other'
        data['sort_index'] = 'old'; data['haplo_insuf_index'] = 'old'
        before = pickle.dumps((data, ped), protocol=5)
        events = []
        original = acmg_scoring.compute_control_common
        def context(current, actual_ped, family, binary):
            events.append('R2')
            self.assertIs(actual_ped, ped); self.assertEqual((family, binary), ('F', 'one-binary'))
            self.assertEqual(len(current), 2)
            self.assertNotIn('haplo_insuf_index', current)
            self.assertNotIn('control_common_index', current)
            pd.testing.assert_series_equal(current['sort_index'], current['ACMG_quant_score'], check_names=False)
            return None
        def control(*args, **kwargs):
            events.append('C1'); return original(*args, **kwargs)
        with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(
                bridge, 'try_pedigree_context', side_effect=context), patch.object(
                acmg_scoring, 'compute_control_common', side_effect=control), patch.object(
                terminal, 'try_terminal_order', return_value=None):
            actual = rank(data, ped, 'F', ranking_binary='one-binary')
        wanted = rank(*pickle.loads(before), 'F')
        pd.testing.assert_frame_equal(actual, wanted, check_exact=True)
        self.assertEqual(events, ['R2', 'C1'])
        self.assertEqual(pickle.dumps((data, ped), protocol=5), before)

    def test_selected_source_row_and_mask_metadata_are_retained(self):
        data, ped = frame(), pedigree()
        data.index = pd.MultiIndex.from_tuples([('a', 1), ('a', 1), ('b', 9)], names=['first', 'second'])
        data.columns.name = 'sample columns'; data.attrs = {'nested': {'scope': ['raw']}}
        ped.index = pd.Index([7, 7, 3], name='PED source'); ped.columns.name = 'PED columns'
        ped.attrs = {'source': ['pedigree']}
        before = pickle.dumps((data, ped), protocol=5)
        with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs: response(command)):
            actual = bridge.try_pedigree_context(data, ped, 'F', 'native')
        self.assertIsNotNone(actual)
        pd.testing.assert_series_equal(actual['proband_row'], ped.iloc[[1]].iloc[0], check_exact=True)
        self.assertIs(actual['proband'], ped.iloc[[1]].iloc[0]['IndividualID'])
        expected = data['PROBAND'].str.count('1') == 1
        pd.testing.assert_series_equal(actual['proband_het'], expected, check_exact=True)
        self.assertEqual(actual['proband_het'].attrs, expected.attrs)
        self.assertEqual(actual['proband_het'].flags.allows_duplicate_labels, expected.flags.allows_duplicate_labels)
        self.assertEqual(pickle.dumps((data, ped), protocol=5), before)

    def test_nullable_raw_mask_bits_and_label_flags_are_retained(self):
        data, ped = frame(('0/1:10', '0/0:1', pd.NA), 'string'), pedigree()
        data.attrs = {'scope': ['caller']}; data.flags.allows_duplicate_labels = False
        expected = data['PROBAND'].str.count('1') == 1
        before = pickle.dumps((data, ped), protocol=5)
        with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs: response(
                command, values=(0, 1, 1), nulls=(0, 0, 1), dtype='boolean')):
            actual = bridge.try_pedigree_context(data, ped, 'F', 'native')
        pd.testing.assert_series_equal(actual['proband_het'], expected, check_exact=True)
        # Inspect the raw storage already protected by the pinned original R2
        # oracle; assert_series_equal alone does not compare masked data bits.
        np.testing.assert_array_equal(actual['proband_het'].array._data, expected.array._data)
        np.testing.assert_array_equal(actual['proband_het'].array._mask, expected.array._mask)
        self.assertFalse(actual['proband_het'].flags.allows_duplicate_labels)
        self.assertEqual(pickle.dumps((data, ped), protocol=5), before)

    def test_disabled_pedigree_empty_mask_has_original_dtype_and_metadata(self):
        data = frame(())
        data.attrs = {'ignored': 'not inherited'}; data.flags.allows_duplicate_labels = False
        with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs: response(
                command, values=(), selected=None, male='0', female='0', warning='no_pedigree', dtype='object')):
            actual = bridge.try_pedigree_context(data, None, None, 'native')
        expected = pd.Series([False] * len(data), index=data.index)
        pd.testing.assert_series_equal(actual['proband_het'], expected, check_exact=True)
        self.assertEqual(actual['proband_het'].attrs, {})
        self.assertTrue(actual['proband_het'].flags.allows_duplicate_labels)

    def test_large_python_family_scalar_is_not_rounded(self):
        data, ped = frame(), pedigree()
        ped['#FamilyID'] = pd.Series([0, 2**53 + 1, float(2**53)], dtype=object)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bridge._input(data, ped, 2**53 + 1, root)
            with (root / 'family.csv').open(newline='') as stream:
                self.assertEqual(list(csv.reader(stream)), [['kind', 'value'], ['int', str(2**53 + 1)]])
            with (root / 'pedigree.csv').open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual((rows[1]['#FamilyID_kind'], rows[1]['#FamilyID']), ('int', str(2**53 + 1)))

    def test_numpy_family_scalars_and_object_family_cells_fall_back(self):
        for value in (np.int64(1), np.float64(1), np.bool_(True), np.str_('F')):
            with self.subTest(scalar=type(value).__name__), patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_pedigree_context(frame(), pedigree(), value, 'native'))
                native.assert_not_called()
        ped = pedigree(); ped['#FamilyID'] = pd.Series(['OTHER', np.int64(2**53 + 1), 0], dtype=object)
        with patch.object(bridge.subprocess, 'run') as native:
            self.assertIsNone(bridge.try_pedigree_context(frame(), ped, float(2**53), 'native'))
        native.assert_not_called()

    def test_unsupported_contexts_and_custom_objects_fall_back_without_callbacks(self):
        class Custom:
            def __str__(self):
                raise AssertionError('custom scalar must remain with original source')
            def __deepcopy__(self, _memo):
                raise AssertionError('custom metadata must remain with original source')
        cases = []
        data = frame(); data.loc[0, 'PROBAND'] = Custom(); cases.append(data)
        data = frame(); data.attrs = {'custom': Custom()}; cases.append(data)
        data = frame(); data['PROBAND'] = pd.Categorical(data['PROBAND']); cases.append(data)
        for data in cases:
            with patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_pedigree_context(data, pedigree(), 'F', 'native'))
                native.assert_not_called()
        for library in (pd, np):
            with patch.object(library, '__version__', 'unqualified'), patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(bridge.try_pedigree_context(frame(), pedigree(), 'F', 'native'))
                native.assert_not_called()

    def test_original_warning_precedes_missing_sample_failure_once(self):
        data, ped = frame().drop(columns='PROBAND'), pedigree('unknown')
        for native in (False, True):
            events = []
            def warning(template, *args):
                events.append(('warning', template, args))
            with self.subTest(native=native), patch.object(acmg_scoring.logger, 'warning', side_effect=warning), patch.object(
                    setup, 'try_rank_setup', return_value=None), patch.object(
                    bridge.subprocess, 'run', return_value=SimpleNamespace(returncode=1)), patch.object(
                    acmg_scoring, 'compute_control_common') as control:
                with self.assertRaises(KeyError) as error:
                    rank(data, ped, 'F', **({'ranking_binary': 'native'} if native else {}))
                events.append(('failure', error.exception.args))
            control.assert_not_called()
            self.assertEqual(len(events), 2)
            self.assertEqual(events[0][0], 'warning')
            self.assertEqual(events[0][2], ('unknown', 'PROBAND'))
            self.assertEqual(events[1], ('failure', ('PROBAND',)))

    def test_accepted_warning_occurs_before_control_and_tail_override_is_retained(self):
        data, ped = frame(), pedigree('unknown')
        events, original = [], acmg_scoring.compute_control_common
        accepted = dict(proband='PROBAND', proband_row=ped.iloc[[1]].iloc[0],
            proband_het=data['PROBAND'].str.count('1') == 1,
            proband_is_male=False, proband_is_female=False, source_warning='unknown_sex')
        def control(*args, **kwargs):
            events.append('C1'); return original(*args, **kwargs)
        with patch.object(setup, 'try_rank_setup', return_value=None), patch.object(
                bridge, 'try_pedigree_context', return_value=accepted) as context, patch.object(
                acmg_scoring.logger, 'warning', side_effect=lambda *args: events.append('warning')), patch.object(
                acmg_scoring, 'compute_control_common', side_effect=control), patch.object(
                terminal, 'try_terminal_order', return_value=None) as tail:
            rank(data, ped, 'F', ranking_binary='integrated', terminal_order_binary='explicit-tail')
        self.assertEqual(context.call_args.args[-1], 'integrated')
        self.assertEqual(tail.call_args.args[-1], 'explicit-tail')
        self.assertEqual(events[:2], ['warning', 'C1'])

    def test_nonzero_timeout_and_os_failure_leave_both_inputs_unchanged(self):
        data, ped = frame(), pedigree()
        before = pickle.dumps((data, ped), protocol=5)
        for outcome in (SimpleNamespace(returncode=1), subprocess.TimeoutExpired('native', 300), FileNotFoundError('native')):
            with self.subTest(outcome=type(outcome).__name__), patch.object(bridge.subprocess, 'run') as native:
                if isinstance(outcome, Exception): native.side_effect = outcome
                else: native.return_value = outcome
                self.assertIsNone(bridge.try_pedigree_context(data, ped, 'F', 'native'))
            self.assertEqual(pickle.dumps((data, ped), protocol=5), before)

    def test_malformed_native_context_never_publishes_partial_result(self):
        data, ped = frame(), pedigree()
        before = pickle.dumps((data, ped), protocol=5)
        changes = [lambda rows, meta: rows.reverse(),
                   lambda rows, meta: rows[0].update(index_value='wrong'),
                   lambda rows, meta: rows[0].update(value='2'),
                   lambda rows, meta: rows[0].update(is_null='1'),
                   lambda rows, meta: meta.append(meta[0]),
                   lambda rows, meta: meta.__setitem__(-1, ('axis', 'wrong')),
                   lambda rows, meta: meta.__setitem__(2, ('selected_pedigree_position', '-1')),
                   lambda rows, meta: meta.__setitem__(1, ('proband', 'OTHER')),
                   lambda rows, meta: meta.__setitem__(5, ('mask_dtype', 'boolean')),
                   lambda rows, meta: meta.__setitem__(7, ('source_warning', 'unknown_sex'))]
        for change in changes:
            with self.subTest(change=changes.index(change)), patch.object(bridge.subprocess, 'run',
                    side_effect=lambda command, **kwargs: response(command, change=change)):
                self.assertIsNone(bridge.try_pedigree_context(data, ped, 'F', 'native'))
            self.assertEqual(pickle.dumps((data, ped), protocol=5), before)

    def test_failure_after_selected_row_reconstruction_is_atomic(self):
        data, ped = frame(('1', None), 'string'), pedigree()
        before = pickle.dumps((data, ped), protocol=5)
        with patch.object(bridge.subprocess, 'run', side_effect=lambda command, **kwargs: response(
                command, values=(1, 1), nulls=(0, 1), dtype='boolean')), patch.object(
                pd.arrays, 'BooleanArray', side_effect=RuntimeError('private mask construction failed')):
            self.assertIsNone(bridge.try_pedigree_context(data, ped, 'F', 'native'))
        self.assertEqual(pickle.dumps((data, ped), protocol=5), before)

    def test_process_resource_errors_propagate(self):
        for error in (MemoryError('resource'), RecursionError('depth')):
            with self.subTest(error=type(error).__name__), patch.object(bridge.subprocess, 'run', side_effect=error):
                with self.assertRaises(type(error)):
                    bridge.try_pedigree_context(frame(), pedigree(), 'F', 'native')


if __name__ == '__main__':
    unittest.main()
