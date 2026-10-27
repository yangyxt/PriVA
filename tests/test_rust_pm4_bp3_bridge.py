"""Portable source, package, raw-buffer and fallback controls for PM4/BP3."""
import hashlib
import importlib.util
import io
import json
import logging
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / 'scripts'
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(ROOT / 'tests/fixtures'))
import priva_rust_pm4_bp3 as adapter
import priva_rust_pm4_bp3_source as guards
import priva_rust_pm4_bp3_packages as package_guards
import pm4_bp3_cases as cases
import pm4_bp3_state as state

class Admission(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        reference = SCRIPT_DIR / 'acmg_pm4_bp3_protein_length.py'
        spec = importlib.util.spec_from_file_location('acmg_pm4_bp3_protein_length_guard_tests', reference)
        cls.reference = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.reference)
        cls.binding, cls.digest = adapter._package_profile()
        cls.temporary = tempfile.TemporaryDirectory(prefix='priva-pm4-bp3-test-')
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.scratch = Path(cls.temporary.name)

    def setUp(self):
        self.reference.logger.handlers = [logging.StreamHandler(io.StringIO())]
        self.reference.logger.propagate = False
        binary = Path(sys.executable).resolve()
        self.prototype = adapter._Bridge(self.reference, binary=binary,
            binary_sha256=adapter.sha(binary), scratch=self.scratch,
            package_binding=self.binding, package_binding_sha256=adapter.sha(self.binding))

    def admit(self, frame, pvs, loc): return self.prototype._admit(frame, whole=True, pvs=pvs, loc=loc)

    def test_actual_original_source_and_full_code_guard(self): guards.reference_guard(self.reference)

    def test_positive_admission_does_not_populate_cold_cache(self):
        df, pvs, loc = cases.base()
        before = state.frame(df)
        result = self.admit(df, pvs, loc)
        self.assertEqual(result['ref'], ['A', 'AAA', 'AAAA', 'AA'])
        self.assertEqual(before, state.frame(df))

    def test_warm_and_stale_cache_values_and_aliases_are_preserved(self):
        for name in ('whole-warm-cache', 'whole-stale-ref-cache', 'whole-stale-id-cache', 'whole-stale-consequence-cache'):
            with self.subTest(name=name):
                fixture = cases.build(name); df = fixture['frame']
                before = state.frame(df); aliases = {key: id(value) for key, value in df._item_cache.items()}
                self.admit(df, fixture['pvs'], fixture['loc'])
                self.assertEqual(before, state.frame(df))
                self.assertEqual(aliases, {key: id(value) for key, value in df._item_cache.items()})

    def test_nullable_raw_data_and_mask_are_not_changed(self):
        fixture = cases.build('whole-nullable-booleans'); df = fixture['frame']
        before = state.frame(df)
        self.admit(df, fixture['pvs'], fixture['loc'])
        self.assertEqual(before, state.frame(df))

    def test_late_unsupported_input_does_not_prime_columns(self):
        df, pvs, loc = cases.base(); df['span_functional_domains'] = np.arange(4)
        before = state.frame(df)
        with self.assertRaises(ValueError): self.admit(df, pvs, loc)
        self.assertEqual(before, state.frame(df))

    def test_custom_scalar_is_rejected_without_length_callback(self):
        fixture = cases.build('whole-ref-callback-raises'); cases.CALLBACKS.clear()
        with self.assertRaises(ValueError): self.admit(fixture['frame'], fixture['pvs'], fixture['loc'])
        self.assertEqual(cases.CALLBACKS, [])

    def test_changed_helper_is_rejected_before_frame_access(self):
        df, pvs, loc = cases.base(); before = state.frame(df)
        with mock.patch.object(self.reference, 'find_overlaps_bedtools_efficient', lambda *args, **kwargs: set()):
            with self.assertRaises(ValueError): self.admit(df, pvs, loc)
        self.assertEqual(before, state.frame(df))

    def test_changed_default_object_is_not_compared(self):
        class EqualProbe:
            def __eq__(self, other): raise AssertionError('custom equality invoked')
        function = self.reference.find_overlaps_bedtools_efficient
        before = function.__defaults__
        try:
            function.__defaults__ = (EqualProbe(),)
            with self.assertRaises(ValueError): guards.reference_guard(self.reference)
        finally: function.__defaults__ = before

    def test_changed_pandas_getter_rejected_before_callback(self):
        df, pvs, loc = cases.base(); seen = []
        def getter(*args, **kwargs): seen.append('getter'); raise AssertionError('unexpected getter')
        with mock.patch.object(pd.DataFrame, '__getitem__', getter):
            with self.assertRaises(ValueError): self.admit(df, pvs, loc)
        self.assertEqual(seen, [])
        self.assertEqual(df._item_cache, {})

    def test_custom_logger_is_rejected_before_cache_access(self):
        df, pvs, loc = cases.base()
        self.reference.logger.handlers.append(cases.RaiseAtLog(1))
        cases.CALLBACKS.clear()
        with self.assertRaises(ValueError): self.admit(df, pvs, loc)
        self.assertEqual(df._item_cache, {})
        self.assertEqual(cases.CALLBACKS, [])

    def test_numeric_id_requires_exact_original(self):
        fixture = cases.build('whole-integer-ids')
        with self.assertRaises(ValueError): self.admit(fixture['frame'], fixture['pvs'], fixture['loc'])

    def test_binary_digest_change_is_rejected_without_cache_mutation(self):
        df, pvs, loc = cases.base(); self.prototype.binary_sha256 = '0' * 64
        with self.assertRaises(ValueError): self.admit(df, pvs, loc)
        self.assertEqual(df._item_cache, {})

    def test_python_default_returns_before_source_or_frame_access(self):
        self.prototype.binary = None
        with mock.patch.object(guards, 'reference_guard', side_effect=AssertionError('source guard should not run')):
            with self.assertRaises(ValueError): self.prototype._admit(object())

    def test_public_default_never_reads_source_or_prepares_a_tail(self):
        sentinel = object()
        with mock.patch.object(guards, 'regular_bytes', side_effect=AssertionError('default source read')):
            prototype = adapter._Bridge(self.reference, binary=None, binary_sha256=None,
                scratch=self.scratch, package_binding=self.binding, package_binding_sha256=None)
            with mock.patch.object(self.reference, 'PM4_BP3_criteria', return_value=sentinel) as original:
                self.assertIs(prototype.PM4_BP3_criteria(object(), object(), 'unchanged', object()), sentinel)
                original.assert_called_once()
        self.assertIsNone(prototype._tail)

    def test_actual_logger_filename_function_and_line_are_original(self):
        records = []
        class Capture(logging.Handler):
            def emit(self, record): records.append(record)
        self.reference.logger.handlers = [Capture()]
        self.reference.logger.setLevel(logging.DEBUG)
        emitters = self.prototype._original_log_sites()
        values = [np.zeros(4, dtype=np.bool_) for _ in range(5)] + [np.zeros(4, dtype=np.int64), np.zeros(4, dtype=np.int64), np.zeros(4, dtype=np.bool_)]
        for function, line in emitters:
            function(*values)
            self.assertEqual(records[-1].name, self.reference.__name__)
            self.assertEqual(records[-1].pathname, self.reference.__file__)
            self.assertEqual(records[-1].funcName, 'PM4_BP3_criteria')
            self.assertEqual(records[-1].lineno, line)
        self.assertEqual(len(records), 9)

    def test_post_intersection_tail_keeps_original_function_name(self):
        tail = self.prototype._original_tail()
        self.assertEqual(tail.__code__.co_name, 'PM4_BP3_criteria')
        self.assertEqual(tail.__code__.co_filename, self.reference.__file__)

    def test_bed_end_keeps_source_string_nullable_storage(self):
        df, _pvs, _loc = cases.base()
        df['ref'] = df['ref'].astype('string[python]')
        expected = df['pos'] - 1 + df.get('ref', 'A').str.len()
        before = state.frame(df); aliases = {key: id(value) for key, value in df._item_cache.items()}
        values = np.array([11, 23, 34, 42], dtype=np.int64)
        actual = self.prototype._bed_end(df, values)
        self.assertEqual(state.array(actual), state.array(expected._mgr.blocks[0].values))
        self.assertEqual(before, state.frame(df))
        self.assertEqual(aliases, {key: id(value) for key, value in df._item_cache.items()})

    def test_bed_end_leaves_source_numpy_storage_and_alias_unchanged(self):
        df, _pvs, _loc = cases.base()
        expected = df['pos'] - 1 + df.get('ref', 'A').str.len()
        values = np.array([11, 23, 34, 42], dtype=np.int64)
        before = state.frame(df)
        actual = self.prototype._bed_end(df, values)
        self.assertIs(actual, values)
        self.assertEqual(state.array(actual), state.array(expected._mgr.blocks[0].values))
        self.assertEqual(before, state.frame(df))

class CanonicalAdmission(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        reference = SCRIPT_DIR / 'acmg_pm4_bp3_protein_length.py'
        spec = importlib.util.spec_from_file_location('acmg_pm4_bp3_protein_length_guard_tests', reference)
        cls.reference = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.reference)
        cls.binding, cls.digest = adapter._package_profile()
        cls.temporary = tempfile.TemporaryDirectory(prefix='priva-pm4-bp3-test-')
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.scratch = Path(cls.temporary.name)

    def setUp(self):
        self.reference.logger.handlers = [logging.StreamHandler(io.StringIO())]
        self.reference.logger.propagate = False
        self.reference.logger.setLevel(logging.DEBUG)
        binary = Path(sys.executable).resolve()
        self.prototype = adapter._Bridge(self.reference, binary=binary,
            binary_sha256=adapter.sha(binary), scratch=self.scratch,
            package_binding=self.binding, package_binding_sha256=self.digest)

    def package(self): package_guards.guard(self.binding, self.digest)

    def test_pristine_named_source_and_closure_profile_passes(self): self.package()

    def test_selected_complete_profile_passes_and_other_binary_profile_fails(self):
        self.package()
        profiles = [(path, digest) for _extension, path, digest in adapter._PACKAGE_PROFILES]
        self.assertEqual(len(profiles), 2)
        other = [(path, digest) for path, digest in profiles if path != self.binding]
        self.assertEqual(len(other), 1)
        with self.assertRaisesRegex(ValueError, 'NumPy extension changed'):
            package_guards.guard(*other[0])

    def test_unrecognized_binary_is_not_admitted_using_matching_versions(self):
        with mock.patch.object(adapter, 'sha', return_value='0' * 64):
            with self.assertRaisesRegex(ValueError, 'unqualified NumPy binary profile'):
                adapter._package_profile()

    def test_changed_profile_is_rejected_using_its_pinned_digest(self):
        proof = json.loads(self.binding.read_bytes())
        proof['source_files']['pandas.core.series'] = '0' * 64
        with tempfile.TemporaryDirectory() as temporary:
            altered = Path(temporary) / 'proof.json'
            altered.write_text(json.dumps(proof))
            with self.assertRaisesRegex(ValueError, 'package proof changed'):
                package_guards.guard(altered, self.digest)

    def test_wrong_canonical_method_is_not_trusted_as_current_identity(self):
        with mock.patch.object(pd.Series, 'fillna', pd.Series.isin):
            with self.assertRaises(ValueError): self.package()

    def test_contains_allowed_types_closure_is_source_derived(self):
        from pandas.core.strings.accessor import StringMethods
        function = StringMethods.contains
        cell = dict(zip(function.__code__.co_freevars, function.__closure__))['allowed_types']
        previous = cell.cell_contents
        try:
            cell.cell_contents = set()
            with self.assertRaises(ValueError): self.package()
        finally: cell.cell_contents = previous

    def test_operator_flag_and_wrapped_function_are_source_bound(self):
        function = pd.Series.__and__
        cells = dict(zip(function.__code__.co_freevars, function.__closure__))
        previous = cells['is_cmp'].cell_contents
        try:
            cells['is_cmp'].cell_contents = True
            with self.assertRaises(ValueError): self.package()
        finally: cells['is_cmp'].cell_contents = previous
        previous = function.__dict__['__wrapped__']
        try:
            function.__dict__['__wrapped__'] = pd.Series.isin
            with self.assertRaises(ValueError): self.package()
        finally: function.__dict__['__wrapped__'] = previous

    def test_custom_default_is_rejected_without_equality_or_repr(self):
        calls = []
        class Probe:
            def __eq__(self, other): calls.append('eq'); raise AssertionError('custom equality')
            def __repr__(self): calls.append('repr'); raise AssertionError('custom repr')
        function = pd.Series.fillna; previous = function.__kwdefaults__
        try:
            function.__kwdefaults__ = dict(previous, downcast=Probe())
            with self.assertRaises(ValueError): self.package()
        finally: function.__kwdefaults__ = previous
        self.assertEqual(calls, [])

    def test_c_constructor_changed_before_fresh_guard_import_is_rejected(self):
        calls = []
        def fake(*args, **kwargs): calls.append('array'); raise AssertionError('constructor callback')
        root = SCRIPT_DIR
        with mock.patch.object(np, 'array', fake):
            spec = importlib.util.spec_from_file_location('pm4_fresh_package_guard', root / 'priva_rust_pm4_bp3_packages.py')
            fresh = importlib.util.module_from_spec(spec); spec.loader.exec_module(fresh)
            with self.assertRaises(ValueError): fresh.guard(self.binding, self.digest)
        self.assertEqual(calls, [])

    def test_changed_str_descriptor_is_not_called_during_admission(self):
        calls = []
        def callback(*args): calls.append('str'); raise AssertionError('accessor callback')
        df, pvs, loc = cases.base(); before = state.frame(df)
        with mock.patch.object(pd.Series, 'str', property(callback)):
            with self.assertRaises(ValueError): self.prototype._admit(df, whole=True, pvs=pvs, loc=loc)
        self.assertEqual(calls, []); self.assertEqual(before, state.frame(df))

    def test_normal_cached_string_accessor_preserves_owner_and_cache(self):
        df, pvs, loc = cases.base(); accessor = df['ref'].str
        before = state.frame(df); aliases = {key: id(value) for key, value in df._item_cache.items()}
        result = self.prototype._admit(df, whole=True, pvs=pvs, loc=loc)
        self.assertEqual(result['ref'], ['A', 'AAA', 'AAAA', 'AA'])
        self.assertEqual(before, state.frame(df)); self.assertIs(df._item_cache['ref'].__dict__['str'], accessor)
        self.assertEqual(aliases, {key: id(value) for key, value in df._item_cache.items()})

    def test_cached_series_method_override_is_not_called(self):
        calls = []
        df, pvs, loc = cases.base(); series = df['DOMAINS']
        def callback(*args, **kwargs): calls.append('fillna'); raise AssertionError('cached callback')
        series.__dict__['fillna'] = callback
        before = state.frame(df)
        with self.assertRaises(ValueError): self.prototype._admit(df, whole=True, pvs=pvs, loc=loc)
        self.assertEqual(calls, []); self.assertEqual(before, state.frame(df))

    def test_logger_class_instance_and_factory_changes_are_early(self):
        calls = []
        def callback(*args, **kwargs): calls.append('logger'); raise AssertionError('log callback')
        df, pvs, loc = cases.base(); before = state.frame(df)
        with mock.patch.object(logging.Logger, 'info', callback):
            with self.assertRaises(ValueError): self.prototype._admit(df, whole=True, pvs=pvs, loc=loc)
        with mock.patch.object(self.reference.logger, 'info', callback):
            with self.assertRaises(ValueError): self.prototype._admit(df, whole=True, pvs=pvs, loc=loc)
        with mock.patch.object(logging, '_logRecordFactory', callback):
            with self.assertRaises(ValueError): self.prototype._admit(df, whole=True, pvs=pvs, loc=loc)
        self.assertEqual(calls, []); self.assertEqual(before, state.frame(df))

class PublicFacade(unittest.TestCase):
    def test_none_branch_calls_original_once_without_setup_or_source_reads(self):
        sentinel = object()
        with mock.patch.object(adapter, '_optional_bridge', side_effect=AssertionError('default setup')):
            with mock.patch.object(guards, 'regular_bytes', side_effect=AssertionError('default source read')):
                with mock.patch.object(adapter._reference, 'PM4_BP3_criteria', return_value=sentinel) as original:
                    self.assertIs(adapter.PM4_BP3_criteria(1, 2, 'regions', 4), sentinel)
                    original.assert_called_once_with(1, 2, 'regions', 4)
                with mock.patch.object(adapter._reference, 'find_overlaps_bedtools_efficient', return_value=sentinel) as original:
                    self.assertIs(adapter.find_overlaps_bedtools_efficient(1, 'regions', method='all'), sentinel)
                    original.assert_called_once_with(1, 'regions', method='all')

    def test_missing_binary_calls_original_once(self):
        with tempfile.TemporaryDirectory() as temporary:
            binary = str(Path(temporary) / 'missing')
            sentinel = object()
            with mock.patch.object(adapter._reference, 'PM4_BP3_criteria', return_value=sentinel) as original:
                self.assertIs(adapter.PM4_BP3_criteria(1, 2, 'regions', 4, ranking_binary=binary), sentinel)
                original.assert_called_once_with(1, 2, 'regions', 4)

    def test_mature_error_is_not_replayed_and_private_scratch_is_cleaned(self):
        with tempfile.TemporaryDirectory() as temporary:
            scratch = Path(temporary) / 'transport'; scratch.mkdir()
            (scratch / 'private.tsv').write_text('retained until return')
            error = RuntimeError('mature failure')
            bridge = mock.Mock()
            bridge.PM4_BP3_criteria.side_effect = error
            with mock.patch.object(adapter, '_optional_bridge', return_value=(bridge, scratch)):
                with mock.patch.object(adapter._reference, 'PM4_BP3_criteria', side_effect=AssertionError('replay')):
                    with self.assertRaises(RuntimeError) as caught:
                        adapter.PM4_BP3_criteria(1, 2, 'regions', 4, ranking_binary='tool')
            self.assertIs(caught.exception, error)
            bridge.PM4_BP3_criteria.assert_called_once_with(1, 2, 'regions', 4)
            self.assertFalse(scratch.exists())

    def test_transport_uses_the_two_actual_nested_commands(self):
        with tempfile.TemporaryDirectory() as temporary:
            binary = Path(sys.executable).resolve()
            binding, digest = adapter._package_profile()
            bridge = adapter._Bridge(adapter._reference, binary=binary,
                binary_sha256=adapter.sha(binary), scratch=temporary,
                package_binding=binding, package_binding_sha256=digest)
            commands = []
            def execute(argv, **kwargs):
                commands.append(argv)
                Path(argv[argv.index('--output') + 1]).write_text('ordinal\tvalue\n0\t0\n')
                return type('Completed', (), dict(returncode=0, stdout=b'', stderr=b''))()
            with mock.patch.object(adapter.subprocess, 'run', side_effect=execute):
                bridge._run('prepare-bed', dict(ordinal=[0]))
                bridge._run('pm4-bp3', dict(ordinal=[0]), overlap=set())
            self.assertEqual([argv[:3] for argv in commands],
                [[str(binary), 'pm4-bp3', 'prepare-bed'], [str(binary), 'pm4-bp3', 'assign']])


if __name__ == '__main__': unittest.main()
