"""Portable admission/atomicity controls for the optional posterior seam.

Transport mocks here exercise guards and host publication kinds; real native
numerical and caller comparisons are separate source-pinned qualifications.
"""
import ast
import builtins
import copy
import csv
import importlib.util
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
from types import CodeType, FunctionType, SimpleNamespace
import unittest
from unittest.mock import patch
import warnings

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import acmg_scoring as scoring
import priva_rust_posterior as bridge


def matrix(dtype='int64', *, production=False, index=None):
    current = pd.DataFrame({name: pd.array([0, 0], dtype=dtype) for name in bridge.CRITERIA}, index=index)
    if production:
        current['variant_id'] = ['a', 'b']
        for name in bridge.AUDITS:
            current[name] = np.array([0, 1], dtype=np.int64)
    current.attrs = {'nested': ['retained', 7]}
    current.index.name = 'caller_rows'
    current.columns.name = 'caller_columns'
    return current


def raw_state(current):
    blocks = []
    for block in current._mgr.blocks:
        values = block.values
        if type(values) is np.ndarray:
            blocks.append((id(values), values.dtype.str, values.shape, values.tobytes()))
        elif type(values) in (pd.arrays.IntegerArray, pd.arrays.FloatingArray, pd.arrays.BooleanArray):
            blocks.append((id(values), values._data.tobytes(), values._mask.tobytes()))
        else:
            blocks.append((id(values),))
    cache = current._item_cache
    cached = tuple((name, id(value), id(value._mgr), id(value.index), value.name,
                    id(value._attrs), copy.deepcopy(value._attrs), value.flags.allows_duplicate_labels,
                    tuple((id(block.values), block.values.tobytes()) for block in value._mgr.blocks))
                   for name, value in cache.items())
    maps = tuple(None if value is None else (id(value), value.tobytes())
                 for value in (current._mgr._blknos, current._mgr._blklocs))
    return (blocks, id(current.index), id(current.columns), current.flags.allows_duplicate_labels,
            id(cache), cached, maps)


def invoke(current, function=None, binary='fixture-native'):
    return bridge.try_posterior_probability(current, binary,
        source_function=scoring.calculate_posterior_probability if function is None else function)


def fake_native(argv, **kwargs):
    with Path(argv[argv.index('--input') + 1]).open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    with Path(argv[argv.index('--output') + 1]).open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['source_position', 'probability_kind', 'probability_value', 'ACMG_class'])
        for row in rows:
            if row['BA1'] == '5':
                writer.writerow([row['source_position'], 'int', '0', 'Benign'])
            else:
                writer.writerow([row['source_position'], 'float64', struct.pack('>d', 0.1).hex(),
                                 'Uncertain Significance'])
    return SimpleNamespace(returncode=0)


class PosteriorBridge(unittest.TestCase):
    def test_cold_warm_and_stale_caches_aliases_and_manager_maps_survive_every_atomic_route(self):
        for cache_kind in ('cold', 'warm', 'stale'):
            for route in ('success', 'native_failure', 'late_rejection'):
                with self.subTest(cache=cache_kind, route=route), pd.option_context('mode.copy_on_write', False):
                    current = matrix('object')
                    alias = current._mgr.blocks[0].values[0]
                    if cache_kind == 'warm':
                        retained = current['PS1']
                    elif cache_kind == 'stale':
                        retained = pd.Series([3, 3], index=current.index, name='PS1')
                        current._item_cache['PS1'] = retained
                    if route == 'late_rejection':
                        # Direct raw mutation avoids accidentally populating or
                        # invalidating the caller's existing item cache.
                        current._mgr.blocks[0].values[27, 1] = 'unsupported'
                    observed = []
                    def transport(argv, **kwargs):
                        with Path(argv[argv.index('--input') + 1]).open(newline='') as stream:
                            observed.extend(csv.DictReader(stream))
                        return SimpleNamespace(returncode=73) if route == 'native_failure' else fake_native(argv, **kwargs)
                    before = raw_state(current)
                    with patch.object(bridge.subprocess, 'run', side_effect=transport) as native:
                        plan = invoke(current)
                    self.assertEqual(raw_state(current), before)
                    self.assertTrue(np.shares_memory(alias, current._mgr.blocks[0].values))
                    if cache_kind != 'cold':
                        self.assertIs(current._item_cache['PS1'], retained)
                    if route == 'success':
                        self.assertIsNotNone(plan)
                    else:
                        self.assertIsNone(plan)
                    self.assertEqual(native.call_count, int(route != 'late_rejection'))
                    if observed:
                        self.assertEqual([row['PS1'] for row in observed], ['0', '0'])

    def test_genuine_source_callable_reaches_transport_and_preserves_kinds_metadata_and_buffers(self):
        current = matrix(production=True, index=pd.Index([8, 3]))
        current.loc[3, 'BA1'] = 5
        before, expected = raw_state(current), current.copy(deep=True)
        with patch.object(bridge.subprocess, 'run', side_effect=fake_native) as run:
            plan = invoke(current)
        self.assertEqual(run.call_count, 1)
        self.assertEqual(run.call_args.args[0][1], 'posterior-probability')
        self.assertIs(type(plan), tuple)
        self.assertIs(type(plan[0][0]), np.float64)
        self.assertIs(type(plan[1][0]), int)
        self.assertEqual(plan, ((np.float64(0.1), 'Uncertain Significance'), (0, 'Benign')))
        self.assertEqual(raw_state(current), before)
        pd.testing.assert_frame_equal(current, expected, check_exact=True)
        self.assertEqual(current.attrs, expected.attrs)

    def test_all_BA_and_mixed_tuple_publications_retain_original_pandas_dtype_inference(self):
        for all_benign in (True, False):
            current = matrix()
            current.loc[0, 'BA1'] = 5
            if all_benign:
                current.loc[1, 'BA1'] = 5
            with patch.object(bridge.subprocess, 'run', side_effect=fake_native):
                plan = invoke(current)
            expected = zip(*current.apply(scoring.calculate_posterior_probability, axis=1))
            actual = zip(*plan)
            left, right = pd.DataFrame(index=current.index), pd.DataFrame(index=current.index)
            left['score'], left['class'] = tuple(expected)
            right['score'], right['class'] = tuple(actual)
            pd.testing.assert_frame_equal(left, right, check_exact=True)
            self.assertEqual(str(right['score'].dtype), 'int64' if all_benign else 'float64')

    def test_standard_zero_strength_dtypes_and_duplicate_row_positions_are_admitted(self):
        for dtype in ('object', 'bool', 'boolean', 'int32', 'uint64', 'float32', 'float64',
                      'Int32', 'Int64', 'UInt64', 'Float32', 'Float64'):
            current = matrix(dtype, index=pd.Index([7, 7]))
            before = raw_state(current)
            with self.subTest(dtype=dtype), patch.object(bridge.subprocess, 'run', side_effect=fake_native) as run:
                self.assertIsNotNone(invoke(current))
                self.assertEqual(run.call_count, 1)
                self.assertEqual(raw_state(current), before)

    def test_exact_float_strength_encodings_include_negative_zero_without_normalizing_source(self):
        for kind in (float, np.float32, np.float64):
            for value in (-0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0):
                self.assertEqual(bridge._strength(kind(value)), int(value))

    def test_arbitrary_duplicate_missing_and_reordered_columns_fall_back_without_native(self):
        original = matrix()
        candidates = [original.assign(extra=1), original.assign(PS_extra=3), original.iloc[:, :-1],
                      original[original.columns[::-1]], pd.concat([original, original[['PS1']]], axis=1)]
        for current in candidates:
            before = raw_state(current)
            with patch.object(bridge.subprocess, 'run') as run:
                self.assertIsNone(invoke(current))
            run.assert_not_called()
            self.assertEqual(raw_state(current), before)

    def test_empty_current_matrix_retains_original_empty_apply_zip_failure(self):
        current = matrix().iloc[:0]
        with patch.object(bridge.subprocess, 'run') as run:
            self.assertIsNone(invoke(current))
        run.assert_not_called()
        with self.assertRaises(ValueError):
            probabilities, classes = zip(*current.apply(scoring.calculate_posterior_probability, axis=1))

    def test_missing_nonintegral_out_of_domain_and_string_cells_are_not_coerced(self):
        for value in (None, pd.NA, float('nan'), float('inf'), -1, 6, 1.5, '1'):
            current = matrix('object')
            current.loc[0, 'PS1'] = value
            before = raw_state(current)
            with self.subTest(value=repr(value)), patch.object(bridge.subprocess, 'run') as run:
                self.assertIsNone(invoke(current))
            run.assert_not_called()
            self.assertEqual(raw_state(current), before)

    def test_signaling_NaNs_are_rejected_without_warning_callback_or_error_state_change(self):
        values = [np.array([0x7ff0000000000001], dtype=np.uint64).view(np.float64)[0],
                  np.array([0x7f800001], dtype=np.uint32).view(np.float32)[0],
                  struct.unpack('>d', bytes.fromhex('7ff0000000000001'))[0]]
        old, old_callback = np.geterr(), np.geterrcall()
        events = []
        try:
            np.seterr(all='call')
            np.seterrcall(lambda *args: events.append(args))
            for value in values:
                current = matrix('object')
                current.loc[0, 'PS1'] = value
                before = raw_state(current)
                with warnings.catch_warnings(record=True) as caught, patch.object(bridge.subprocess, 'run') as run:
                    warnings.simplefilter('always')
                    self.assertIsNone(invoke(current))
                run.assert_not_called()
                self.assertEqual(caught, [])
                self.assertEqual(events, [])
                self.assertEqual(np.geterr(), dict.fromkeys(old, 'call'))
                self.assertEqual(raw_state(current), before)
        finally:
            np.seterr(**old)
            np.seterrcall(old_callback)

    def test_custom_scalar_is_rejected_before_conversion_or_class_getters(self):
        events, active = [], [False]
        class Unknown:
            def __getattribute__(self, name):
                if active[0] and name == '__class__':
                    events.append(name)
                    raise RuntimeError('custom class getter')
                return object.__getattribute__(self, name)
            def __int__(self):
                events.append('int')
                raise RuntimeError('custom integer conversion')
        current = matrix('object')
        current.loc[0, 'PS1'] = Unknown()
        before = raw_state(current)
        active[0] = True
        with patch.object(bridge.subprocess, 'run') as run:
            self.assertIsNone(invoke(current))
        run.assert_not_called()
        self.assertEqual(events, [])
        self.assertEqual(raw_state(current), before)

    def test_custom_extension_storage_is_rejected_before_dtype_getters(self):
        active, events = [False], []
        class Unknown(pd.arrays.IntegerArray):
            @property
            def dtype(self):
                if active[0]:
                    events.append('dtype')
                    raise RuntimeError('custom dtype callback')
                return super().dtype
        current = matrix()
        current['PS1'] = Unknown(np.array([0, 0], dtype=np.int64), np.zeros(2, dtype=bool))
        before = raw_state(current)
        active[0] = True
        with patch.object(bridge.subprocess, 'run') as run:
            self.assertIsNone(invoke(current))
        run.assert_not_called()
        self.assertEqual(events, [])
        self.assertEqual(raw_state(current), before)

    def test_custom_raw_numpy_storage_is_rejected_without_dtype_getters(self):
        events = []
        class Unknown(np.ndarray):
            def __getattribute__(self, name):
                if name == 'dtype':
                    events.append('dtype')
                    raise RuntimeError('custom raw dtype callback')
                return super().__getattribute__(name)
        current = matrix('Int64')
        current._mgr.blocks[0].values._data = np.array([0, 0], dtype=np.int64).view(Unknown)
        with patch.object(bridge.subprocess, 'run') as run:
            self.assertIsNone(invoke(current))
        run.assert_not_called()
        self.assertEqual(events, [])

    def test_nonstring_ID_nonboolean_audit_and_custom_metadata_retain_original_boxing(self):
        current = matrix(production=True)
        bad_id = current.copy(); bad_id['variant_id'] = [1, 2]
        bad_audit = current.copy(); bad_audit[bridge.AUDITS[0]] = ['0', '1']
        bad_attrs = current.copy(); bad_attrs.attrs = {'custom': object()}
        for current in (bad_id, bad_audit, bad_attrs):
            with patch.object(bridge.subprocess, 'run') as run:
                self.assertIsNone(invoke(current))
            run.assert_not_called()

    def test_changed_actual_helper_code_rejects_even_with_original_filename(self):
        original = scoring.calculate_posterior_probability
        source = compile('def calculate_posterior_probability(row, prior_probability=.1, exp_base=2, odds_pvst=350):\n return 1, "Pathogenic"\n',
                         original.__code__.co_filename, 'exec')
        code = next(value for value in source.co_consts if isinstance(value, CodeType))
        namespace = dict(original.__globals__)
        changed = FunctionType(code, namespace, original.__name__, original.__defaults__)
        namespace[original.__name__] = changed
        with patch.object(bridge.subprocess, 'run') as run:
            self.assertIsNone(invoke(matrix(), changed))
        run.assert_not_called()

    def test_changed_defaults_are_rejected_without_altering_the_actual_function(self):
        function = scoring.calculate_posterior_probability
        original = function.__defaults__
        try:
            for defaults in ((0.2, 2, 350), (0.1, 2.0, 350), (0.1, 2, 351)):
                function.__defaults__ = defaults
                with patch.object(bridge.subprocess, 'run') as run:
                    self.assertIsNone(invoke(matrix()))
                run.assert_not_called()
        finally:
            function.__defaults__ = original

    def test_changed_conversion_namespace_or_numpy_power_never_executes_replacement(self):
        events = []
        def unexpected(*args):
            events.append(args)
            raise RuntimeError('changed helper binding')
        with patch.dict(scoring.__dict__, {'sum': unexpected}), patch.object(bridge.subprocess, 'run') as run:
            self.assertIsNone(invoke(matrix()))
        run.assert_not_called()
        with patch.object(np, 'power', unexpected), patch.object(bridge.subprocess, 'run') as run:
            self.assertIsNone(invoke(matrix()))
        run.assert_not_called()
        self.assertEqual(events, [])

    def test_changed_public_apply_to_dict_boxing_row_generator_and_raw_accessor_fail_before_frame_access(self):
        from pandas.core import apply, series
        from pandas.core.dtypes import cast
        from pandas.core.internals import blocks, managers
        events = []
        def unexpected(*args, **kwargs):
            events.append('mutated API called')
            raise RuntimeError('mutated pandas API')
        cases = ((pd.DataFrame, 'apply', unexpected), (pd.Series, 'to_dict', unexpected),
                 (cast, 'maybe_box_native', unexpected), (series, 'maybe_box_native', unexpected),
                 (apply.FrameColumnApply, 'series_generator', property(unexpected)),
                 (managers.BlockManager, 'iget_values', unexpected),
                 (blocks.ExtensionBlock, 'iget', unexpected))
        for owner, name, changed in cases:
            current = matrix()
            before = raw_state(current)
            storage = bridge._storage
            with self.subTest(name=name), patch.object(owner, name, changed), \
                    patch.object(bridge, '_storage', wraps=storage) as accessed, \
                    patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(invoke(current))
            accessed.assert_not_called()
            native.assert_not_called()
            self.assertEqual(events, [])
            self.assertEqual(raw_state(current), before)

    def test_changed_apply_and_to_dict_defaults_or_wrapper_closure_are_rejected(self):
        current = matrix()
        wrapper = pd.Series.to_dict
        original = wrapper.__wrapped__
        allow_cell = dict(zip(wrapper.__code__.co_freevars, wrapper.__closure__))['num_allow_args']
        apply_defaults, dict_defaults, allowed = pd.DataFrame.apply.__defaults__, original.__defaults__, allow_cell.cell_contents
        try:
            for changed in ('apply', 'to_dict', 'closure'):
                pd.DataFrame.apply.__defaults__ = ((1,) + apply_defaults[1:]) if changed == 'apply' else apply_defaults
                original.__defaults__ = (list,) if changed == 'to_dict' else dict_defaults
                allow_cell.cell_contents = 0 if changed == 'closure' else allowed
                with self.subTest(changed=changed), patch.object(bridge, '_storage') as accessed, \
                        patch.object(bridge.subprocess, 'run') as native:
                    self.assertIsNone(invoke(current))
                accessed.assert_not_called()
                native.assert_not_called()
        finally:
            pd.DataFrame.apply.__defaults__, original.__defaults__, allow_cell.cell_contents = apply_defaults, dict_defaults, allowed

    def test_mutated_method_code_is_rejected_even_with_original_object_globals_and_debug_identity(self):
        function = pd.DataFrame.apply
        original = function.__code__
        replacement = compile('def changed(*args, **kwargs):\n return None\n', original.co_filename, 'exec')
        code = next(value for value in replacement.co_consts if type(value) is CodeType)
        function.__code__ = code.replace(co_name=original.co_name, co_qualname=original.co_qualname,
                                        co_firstlineno=original.co_firstlineno)
        try:
            with patch.object(bridge, '_storage') as accessed, patch.object(bridge.subprocess, 'run') as native:
                self.assertIsNone(invoke(matrix()))
            accessed.assert_not_called()
            native.assert_not_called()
        finally:
            function.__code__ = original

    def test_fresh_import_rejects_a_different_numpy_ufunc_binding(self):
        current = matrix()
        with patch.object(np, 'power', np.add):
            spec = importlib.util.spec_from_file_location('fresh_posterior_wrong_ufunc', ROOT / 'scripts/priva_rust_posterior.py')
            fresh = importlib.util.module_from_spec(spec); spec.loader.exec_module(fresh)
            with patch.object(fresh.subprocess, 'run') as native:
                self.assertIsNone(fresh.try_posterior_probability(current, 'unused',
                    source_function=scoring.calculate_posterior_probability))
        native.assert_not_called()

    def test_custom_executable_path_does_not_add_fspath_callbacks(self):
        events = []
        class Unknown:
            def __fspath__(self):
                events.append('fspath')
                raise RuntimeError('custom executable path')
        with patch.object(bridge.subprocess, 'run') as run:
            self.assertIsNone(invoke(matrix(), binary=Unknown()))
        run.assert_not_called()
        self.assertEqual(events, [])

    def test_nonzero_timeout_and_missing_native_output_leave_current_buffers_untouched(self):
        for result in (SimpleNamespace(returncode=73), subprocess.TimeoutExpired('fixture', 1), SimpleNamespace(returncode=0)):
            current = matrix()
            before = raw_state(current)
            kwargs = {'side_effect': result} if isinstance(result, Exception) else {'return_value': result}
            with patch.object(bridge.subprocess, 'run', **kwargs) as run:
                self.assertIsNone(invoke(current))
            self.assertEqual(run.call_count, 1)
            self.assertEqual(raw_state(current), before)

    def test_malformed_native_results_are_rejected_atomically(self):
        ordinary = ['0,float64,3fb999999999999a,Uncertain Significance',
                    '1,float64,3fb999999999999a,Uncertain Significance']
        bad_results = [
            ['0,float64,7ff0000000000000,Benign', ordinary[1]],
            ['0,float64,8000000000000000,Benign', ordinary[1]],
            ['0,float64,3ff0000000000001,Pathogenic', ordinary[1]],
            ['0,float64,notbits,Benign', ordinary[1]],
            ['0,int,1,Benign', ordinary[1]],
            ['0,unknown,0,Benign', ordinary[1]],
            ['0,float64,3fb999999999999a,Unknown', ordinary[1]],
            [ordinary[1], ordinary[0]], [ordinary[0]], [ordinary[0], ordinary[0]],
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'out.csv'
            for rows in bad_results:
                path.write_text('source_position,probability_kind,probability_value,ACMG_class\n' + '\n'.join(rows) + '\n')
                with self.assertRaises(ValueError):
                    bridge._gather(path, 2)
            path.write_text('bad_header\n0\n')
            with self.assertRaises(ValueError):
                bridge._gather(path, 2)


class PosteriorMissingAPIs(unittest.TestCase):
    def test_fresh_import_reaches_source_then_actual_version_or_storage_admission(self):
        current = matrix()
        for version in ('unsupported-fixture', '2.2.3'):
            with self.subTest(version=version), patch.object(pd, '__version__', version), \
                    patch.object(pd.arrays, 'NumpyExtensionArray'):
                del pd.arrays.NumpyExtensionArray
                spec = importlib.util.spec_from_file_location('fresh_posterior_missing_API', ROOT / 'scripts/priva_rust_posterior.py')
                fresh = importlib.util.module_from_spec(spec); spec.loader.exec_module(fresh)
                counts = {'guard': 0, 'admission': 0, 'storage': 0}
                guard, admission, storage = fresh._source_guard, fresh._admit, fresh._storage
                def actual_guard(*args):
                    result = guard(*args); counts['guard'] += 1; return result
                def actual_admission(*args):
                    counts['admission'] += 1; return admission(*args)
                def actual_storage(*args):
                    counts['storage'] += 1; return storage(*args)
                before = raw_state(current)
                with patch.object(fresh, '_source_guard', side_effect=actual_guard), \
                        patch.object(fresh, '_admit', side_effect=actual_admission), \
                        patch.object(fresh, '_storage', side_effect=actual_storage), \
                        patch.object(fresh.subprocess, 'run') as run:
                    self.assertIsNone(fresh.try_posterior_probability(current, 'unused',
                        source_function=scoring.calculate_posterior_probability))
                run.assert_not_called()
                self.assertEqual(counts, {'guard': 1, 'admission': 1, 'storage': int(version == '2.2.3')})
                self.assertEqual(raw_state(current), before)


class PosteriorReductionSource(unittest.TestCase):
    @staticmethod
    def harness():
        spec = importlib.util.spec_from_file_location('posterior_reduction_contract', ROOT / 'rust/tests/posterior_bridge_compare.py')
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        return module

    def test_complete_source_derived_count_ledger_and_ignored_value_multiplicities(self):
        harness = self.harness()
        ledger = harness.reduction_ledger(ROOT / 'scripts/acmg_scoring.py')
        self.assertEqual(ledger['pathogenic_category_vectors'], 4845)
        self.assertEqual(ledger['benign_category_vectors'], 1820)
        self.assertEqual(ledger['joint_category_vectors'], 8817900)
        self.assertEqual(ledger['BA_containing_category_vectors'], 6613425)
        self.assertEqual(ledger['no_BA_numeric_representatives'], 232560)
        self.assertEqual(len(ledger['records']['no_BA_exponents']), 48)
        self.assertEqual(ledger['ignored_source_values'], {'pathogenic': [0, 5], 'benign': [0, 4]})
        self.assertEqual(sum(int(item['multiplicity']) for item in ledger['records']['pathogenic']), 6 ** 16)
        self.assertEqual(sum(int(item['multiplicity']) for item in ledger['records']['benign']), 6 ** 12)

    def test_reduction_rejects_mutated_missing_duplicate_async_or_decorated_scalar_source(self):
        harness = self.harness()
        source = (ROOT / 'scripts/acmg_scoring.py').read_text()
        node = harness.definition(source, 'calculate_posterior_probability')
        segment = ast.get_source_segment(source, node)
        candidates = [segment.replace('n == 3', 'n == 4', 1),
                      segment.replace('1/exp_base**3', '1/exp_base**2'),
                      segment.replace('c.startswith("P")', 'c.startswith("B")', 1),
                      source.replace(segment, '', 1), source + '\n' + segment,
                      source.replace('def calculate_posterior_probability(', 'async def calculate_posterior_probability(', 1),
                      source.replace('def calculate_posterior_probability(', '@changed\ndef calculate_posterior_probability(', 1)]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'oracle.py'
            for index, changed in enumerate(candidates):
                path.write_text(changed)
                with self.subTest(index=index), self.assertRaises(AssertionError):
                    harness.reduction_ledger(path)


if __name__ == '__main__':
    unittest.main()
