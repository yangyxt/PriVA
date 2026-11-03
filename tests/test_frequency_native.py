"""Focused frequency protocol, failure-frontier and actual native controls.

Set PRIVA_FREQUENCY_NATIVE_BINARY and PRIVA_FREQUENCY_NATIVE_SHA256 to run the
actual executable controls. A skipped executable class is not native coverage.
These small synthetic controls do not rerun the saved original case catalogue.
"""
import csv
import io
import os
from pathlib import Path
import struct
import subprocess
import sys
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import numpy as np
import pandas as pd
import priva_rust_frequency as bridge


def response(kind, dtype, cells, labels=(), sequence=0):
    stream = io.StringIO(newline=''); writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
    writer.writerow(('result', sequence, kind, dtype, len(cells), bridge.PROFILE))
    for i, label in enumerate(labels): writer.writerow(('index', 0, i, 'label', label, ''))
    for i, value in enumerate(cells): writer.writerow(('cell', 0, i, *bridge._cell(value), ''))
    return stream.getvalue().encode()


class FrequencyWireTests(unittest.TestCase):
    def test_bits_unicode_and_order_survive_encoding(self):
        for value in (0., -0., float('inf'), float('-inf'), float('nan')):
            self.assertEqual(struct.pack('>d', bridge._uncell(*bridge._cell(value))), struct.pack('>d', value))
        value = pd.Series(['λ\t"x"\ny', 'z'], index=['b', 'a'])
        rows = list(csv.reader(io.StringIO(bridge._request(4, 'lower', (value,), 'BS1_criteria@579').decode()), delimiter='\t'))
        self.assertEqual(rows[0][1:5], ['4', 'lower', '1', 'BS1_criteria@579'])
        self.assertEqual([r[4] for r in rows if r[0] == 'cell'], list(value))

    def test_response_sequence_denominator_and_profile_are_strict(self):
        raw = response('array', 'int64', [0, 0])
        # Avoid numpy.testing's import-time lscpu subprocess in this tiny
        # wire-only fixture; the array type, shape and exact values stay explicit.
        observed = bridge._result(raw, 0, (2,), 'zeros')
        self.assertIsInstance(observed, np.ndarray)
        self.assertEqual(observed.shape, (2,))
        self.assertEqual(observed.tolist(), [0, 0])
        for changed in (raw.replace(b'\t0\tarray', b'\t1\tarray', 1),
                        raw.replace(b'int64\t2', b'int64\t3', 1),
                        raw.replace(bridge.PROFILE.encode(), b'unknown')):
            with self.assertRaises(bridge.TransportError): bridge._result(changed, 0, (2,), 'zeros')

    def test_nonconsecutive_foreign_axis_and_wrong_kind_are_rejected(self):
        x = pd.Series([1.], index=['a'])
        raws = (response('series', 'bool', [True], ['sforeign']),
                response('array', 'bool', [True]),
                response('series', 'bool', [True], ['sa']).replace(b'cell\t0\t0', b'cell\t0\t1'))
        for raw in raws:
            with self.assertRaises(bridge.TransportError): bridge._result(raw, 0, (x,), 'notna')

    def test_duplicate_and_multiindex_labels_are_preserved(self):
        index = pd.MultiIndex.from_tuples([('x', 1), ('x', 1), ('x', 2)], names=['group', 'row'])
        x = pd.Series([1., 2., 3.], index=index, name='frequency')
        labels = bridge._value(x)[2]
        self.assertEqual(labels[0], labels[1])
        result = bridge._result(response('series', 'bool', [True, False, True], labels), 0, (x,), 'notna')
        self.assertIs(result.index, index)
        self.assertEqual(result.name, 'frequency')

    def test_unmasked_nullable_nan_is_not_changed_to_na(self):
        values = pd.arrays.FloatingArray(np.array([float('nan')]), np.array([False]))
        with self.assertRaises(bridge.Unsupported): bridge._request(0, 'isna', (pd.Series(values),), 'BA1_criteria@791')

    def test_unknown_scalar_never_calls_string_or_float(self):
        class Unknown:
            def __str__(self): raise AssertionError('string conversion')
            def __float__(self): raise AssertionError('float conversion')
        with self.assertRaises(bridge.Unsupported): bridge._cell(Unknown())

    def test_corrupt_native_response_never_becomes_mature_fallback(self):
        session = bridge.NativeSession('unused', '0' * 64)
        # Package admission has separate bound-source controls below. This
        # test isolates propagation of a corrupt wire response.
        runtime = bridge._Runtime(session, SimpleNamespace(operator=lambda *_args: None), 'unused', mode='mixed')
        with patch.object(session, 'request', side_effect=bridge.TransportError('corrupt')):
            with patch.object(runtime, '_mature') as fallback:
                with self.assertRaises(bridge.TransportError): runtime.binary(149, 'control_false_neg_rate', 'le', pd.Series([1.]), .1)
                fallback.assert_not_called()

    def test_mature_failure_is_not_replayed(self):
        session = bridge.NativeSession('unused', '0' * 64)
        runtime = bridge._Runtime(session, object(), 'unused', mode='mixed')
        calls = []
        def fail(): calls.append('one'); raise ValueError('retained failure')
        # This forwarding/failure fixture calls only fail() above. Attribute
        # its generated bridge frame to synthetic source, so zero-science
        # observation never confuses it with an original procedure frame.
        with patch.object(runtime, '_filename', return_value='<frequency-wire-control>'), \
                patch.object(session, 'request', side_effect=bridge.Unsupported('outside domain')):
            with self.assertRaisesRegex(ValueError, 'retained failure'):
                runtime._native(149, 'control_false_neg_rate', 'le', (1,), fail, (), {})
        self.assertEqual(calls, ['one'])
        self.assertEqual(runtime._callers[('control_false_neg_rate', 149)].co_filename,
                         '<frequency-wire-control>')

    def test_native_math_refuses_mixed_without_calling_it(self):
        session = bridge.NativeSession('unused', '0' * 64)
        runtime = bridge._Runtime(session, object(), 'unused', mode='native_math')
        with patch.object(session, 'request', side_effect=bridge.Unsupported('outside domain')):
            with self.assertRaises(bridge.Unsupported):
                runtime._native(149, 'control_false_neg_rate', 'le', (1,), lambda: self.fail('must not run'), (), {})

    def test_invalid_package_pin_is_rejected_before_binding_read(self):
        with self.assertRaisesRegex(ValueError, 'SHA'):
            bridge.packages.Guard('must-not-be-opened', 'not-a-digest', object())

    def test_grouped_response_rejects_duplicates_foreign_keys_and_types(self):
        args = (pd.Series([.1, .2]), pd.Series(['a', 'b']))
        invalid = (['a', .1, 'a', .2], ['a', .1], ['a', .1, 'c', .2],
                   ['b', .2, 'a', .1], [True, .1, 'b', .2],
                   ['a', 1, 'b', .2], ['a'])
        for cells in invalid:
            with self.subTest(cells=cells):
                with self.assertRaises(bridge.TransportError):
                    bridge._result(response('pairs', 'object', cells), 0, args, 'group_max_pairs')
        self.assertEqual(bridge._result(response('pairs', 'object', ['a', .1, 'b', .2]),
                                       0, args, 'group_max_pairs'), {'a': .1, 'b': .2})

    def test_corrupt_grouped_result_never_replays_mature_max(self):
        session = bridge.NativeSession('unused', '0' * 64)
        runtime = bridge._Runtime(session, object(), 'unused', mode='mixed')
        args = (pd.Series([.1, .2]), pd.Series(['a', 'b']))
        raw = response('pairs', 'object', ['a', .1, 'a', .2])
        def invalid(*_args, **_kwargs): return bridge._result(raw, 0, args, 'group_max_pairs')
        with patch.object(session, 'request', side_effect=invalid), patch.object(runtime, '_mature') as mature:
            with self.assertRaises(bridge.TransportError):
                runtime._native(65, '_gene_max_patho_hom_freq', 'group_max_pairs', args,
                                lambda: self.fail('aggregate replay'), (), {})
            mature.assert_not_called()


class FrequencyOwnedChildTests(unittest.TestCase):
    """Real tiny owned Python children exercise transport cleanup only.

    No original criterion, SciPy call or Rust build is involved. These are
    lifecycle controls, separate from the pinned actual-native column tests.
    """
    def child(self, behavior):
        body = """import os,signal,struct,sys,time
behavior=sys.argv[1]
if behavior=='hang': signal.signal(signal.SIGTERM,signal.SIG_IGN)
def read(n):
 out=b''
 while len(out)<n:
  part=os.read(0,n-len(out))
  if not part: return out
  out+=part
 return out
header=read(4)
if header:
 read(struct.unpack('>I',header)[0])
 if behavior!='noack':
  value=('result\\t0\\tscalar\\tscalar\\t1\\t'+sys.argv[2]+'\\ncell\\t0\\t0\\tstr\\tclosed\\t\\n').encode()
  os.write(1,struct.pack('>I',len(value))+value)
 if behavior=='stderr':
  value=b'd'*200000
  while value: value=value[os.write(2,value):]
 if behavior in ('hang','noack'):
  while True: time.sleep(1)
"""
        process = subprocess.Popen([sys.executable, '-B', '-c', body, behavior, bridge.PROFILE],
                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, bufsize=0)
        os.set_blocking(process.stdin.fileno(), False)
        session = bridge.NativeSession('unused', '0' * 64, timeout=1.)
        session.process = process
        def emergency_cleanup():
            if process.poll() is None:
                process.kill(); process.wait(timeout=3)
            for stream in (process.stdin, process.stdout, process.stderr):
                if not stream.closed: stream.close()
        self.addCleanup(emergency_cleanup)
        return session, process

    def test_stdin_close_error_still_reaps_and_retains_failure(self):
        session, process = self.child('eof')
        original_close = process.stdin.close
        def broken_close(): original_close(); raise OSError('synthetic close failure')
        with patch.object(process.stdin, 'close', side_effect=broken_close), patch.object(bridge, 'CLOSE_PHASE_SECONDS', .15):
            with self.assertRaises(bridge.TransportError): session.close(failed=True)
        self.assertIsNotNone(process.returncode)
        self.assertTrue(session.metrics['native_close']['reaped'])
        self.assertEqual(session.metrics['native_close']['errors'][0]['phase'], 'stdin-close')

    def test_full_diagnostic_pipe_drains_without_false_wait_timeout(self):
        session, process = self.child('stderr')
        with patch.object(bridge, 'CLOSE_PHASE_SECONDS', 1.):
            with self.assertRaisesRegex(bridge.TransportError, 'extra output'): session.close()
        value = session.metrics['native_close']
        self.assertEqual(process.returncode, 0)
        self.assertEqual(value['termination_attempts'], [])
        self.assertEqual(value['streams']['stderr']['bytes'], 200000)
        self.assertEqual(len(bytes.fromhex(value['streams']['stderr']['retained_hex'])), bridge.CLOSE_DIAGNOSTIC_BYTES)
        self.assertTrue(value['streams']['stderr']['truncated'])

    def test_hung_owned_child_is_boundedly_reaped(self):
        session, process = self.child('hang')
        began = time.monotonic()
        with patch.object(bridge, 'CLOSE_PHASE_SECONDS', .15):
            with self.assertRaises(bridge.TransportError): session.close()
        self.assertLess(time.monotonic() - began, 3.)
        self.assertIsNotNone(process.returncode)
        self.assertEqual(session.metrics['native_close']['termination_attempts'], ['terminate', 'kill'])
        self.assertTrue(session.metrics['native_close']['reaped'])

    def test_quit_timeout_is_primary_after_owned_child_cleanup(self):
        session, process = self.child('noack')
        with patch.object(bridge, 'CLOSE_PHASE_SECONDS', .15):
            with self.assertRaisesRegex(bridge.TransportError, 'native response timeout'): session.close()
        self.assertIsNotNone(process.returncode)
        self.assertTrue(session.metrics['native_close']['reaped'])

    def test_scientific_exception_survives_stdin_cleanup_failure(self):
        session, process = self.child('eof')
        failure = ValueError('primary scientific failure')
        def fail(): raise failure
        owner = bridge.FrequencyBridge.__new__(bridge.FrequencyBridge)
        owner.functions = {'BA1_criteria': fail}; owner.used = False
        owner.session = session; owner.runtime = SimpleNamespace(metrics=session.metrics)
        original_close = process.stdin.close
        def broken_close(): original_close(); raise OSError('synthetic cleanup failure')
        with patch.object(process.stdin, 'close', side_effect=broken_close), patch.object(bridge, 'CLOSE_PHASE_SECONDS', .15):
            with self.assertRaises(ValueError) as result: owner.run('BA1_criteria')
        self.assertIs(result.exception, failure)
        self.assertTrue(session.metrics['native_close']['reaped'])
        self.assertEqual(session.metrics['close_failure']['type'], 'TransportError')

    def test_wait_error_is_retained_after_actual_reap(self):
        session, process = self.child('eof')
        actual_wait = process.wait
        def failed_wait(*args, **kwargs):
            actual_wait(*args, **kwargs)
            raise OSError('synthetic wait observation failure')
        with patch.object(process, 'wait', side_effect=failed_wait), patch.object(bridge, 'CLOSE_PHASE_SECONDS', .15):
            with self.assertRaises(bridge.TransportError): session.close(failed=True)
        self.assertTrue(session.metrics['native_close']['reaped'])
        self.assertEqual(session.metrics['native_close']['errors'][-1]['phase'], 'final-wait')

    def test_primary_quit_error_survives_unexpected_cleanup_fault(self):
        session, process = self.child('noack')
        with patch.object(session, '_finish_child', side_effect=OSError('synthetic cleanup fault')), patch.object(bridge, 'CLOSE_PHASE_SECONDS', .15):
            with self.assertRaisesRegex(bridge.TransportError, 'native response timeout'): session.close()
        self.assertIsNotNone(process.returncode)
        self.assertTrue(session.metrics['native_close_emergency']['reaped'])


@unittest.skipUnless(os.environ.get('PRIVA_FREQUENCY_PACKAGE_BINDING') and os.environ.get('PRIVA_FREQUENCY_PACKAGE_SHA256'),
                     'requires separately reviewed exact package source binding')
class FrequencyPackageDeclarationTests(unittest.TestCase):
    def setUp(self):
        self.guard = bridge.packages.Guard(os.environ['PRIVA_FREQUENCY_PACKAGE_BINDING'],
            os.environ['PRIVA_FREQUENCY_PACKAGE_SHA256'], bridge._reference)

    def test_mutated_live_isnull_slot_cannot_authorize_itself(self):
        value = pd.Series([1., np.nan])
        captured = value.isnull
        self.assertTrue(self.guard.method_is_canonical(captured, 'isnull'))
        with patch.object(pd.Series, 'isnull', lambda self: self):
            self.assertFalse(self.guard.method_is_canonical(value.isnull, 'isnull'))
            # Original bound-method lookup precedes an argument callback that
            # might mutate the class. The captured implementation is still
            # checked against source, not against that new mutable slot.
            self.assertTrue(self.guard.method_is_canonical(captured, 'isnull'))
        self.guard.check(complete=True)

    def test_undeclared_method_never_bypasses_its_live_operation(self):
        value = pd.Series([1.])
        self.assertFalse(self.guard.method_is_canonical(value.tolist, 'tolist'))
        with self.assertRaisesRegex(ValueError, 'source declaration'):
            self.guard.declared_method('pandas.core.series:Series.tolist')

    def test_grouping_and_dictionary_bypasses_have_source_declarations(self):
        grouped = pd.Series([1., 2.]).groupby(pd.Series(['a', 'b']))
        self.assertTrue(self.guard.method_is_canonical(grouped.max, 'max'))
        self.guard.declared_method('pandas.core.series:Series.to_dict')
        with patch.object(pd.Series, 'to_dict', lambda self: {'forged': 0.}):
            with self.assertRaises(ValueError):
                self.guard.declared_method('pandas.core.series:Series.to_dict')

    def test_gzip_named_default_is_source_bound_and_mutations_reject(self):
        import gzip
        original_node = self.guard.context['gzip']['definitions']['open'][0]
        original_tree = bridge.guards.ast.dump(original_node, include_attributes=True)
        self.guard.check(complete=True)
        with patch.object(gzip, '_COMPRESS_LEVEL_BEST', 8):
            with self.assertRaisesRegex(ValueError, 'global differs'):
                self.guard.check()
        with patch.object(gzip, '_COMPRESS_LEVEL_BEST', 9.0):
            with self.assertRaisesRegex(ValueError, 'global differs'):
                self.guard.check()
        defaults = gzip.open.__defaults__
        changed = (defaults[0], 8, *defaults[2:])
        with patch.object(gzip.open, '__defaults__', changed):
            with self.assertRaisesRegex(ValueError, 'source defaults changed'):
                self.guard.check()
        with patch.object(gzip, 'open', lambda *_args, **_kwargs: self.fail('decoder must not run')):
            with self.assertRaisesRegex(ValueError, 'executable function changed'):
                self.guard.check()
        self.guard.check(complete=True)
        self.assertEqual(bridge.guards.ast.dump(original_node, include_attributes=True), original_tree)

    def test_captured_foreign_numeric_or_numpy_callable_is_rejected(self):
        foreign = lambda *args, **kwargs: None
        with self.assertRaises(ValueError): self.guard.numeric(foreign)
        # The live NumPy slot is canonical; it cannot authorize this different
        # callable captured before an argument restored that slot.
        with self.assertRaises(ValueError): self.guard.array_function('where', foreign)

    def test_mutated_series_math_slot_cannot_reach_native_math(self):
        value = pd.Series([1.])
        self.guard.operator('mul', value, 2.)
        with patch.object(pd.Series, '__mul__', lambda self, other: self):
            with self.assertRaises(ValueError): self.guard.operator('mul', value, 2.)
        self.guard.check(complete=True)


@unittest.skipUnless(os.environ.get('PRIVA_FREQUENCY_NATIVE_BINARY') and os.environ.get('PRIVA_FREQUENCY_NATIVE_SHA256'),
                     'requires explicitly pinned actual frequency executable')
class FrequencyActualNativeTests(unittest.TestCase):
    def setUp(self):
        self.session = bridge.NativeSession(os.environ['PRIVA_FREQUENCY_NATIVE_BINARY'], os.environ['PRIVA_FREQUENCY_NATIVE_SHA256'])
    def tearDown(self): self.session.close(failed=sys.exc_info()[0] is not None)
    def request(self, name, *args):
        return self.session.request(name, args, site='control_false_neg_rate@149')

    def test_round_ties_even(self):
        np.testing.assert_array_equal(self.request('round', np.array([.5, 1.5, 2.5, -.5, -1.5])), [0., 2., 2., -0., -2.])

    def test_integer_arithmetic(self):
        x = np.array([1, 2, 3], dtype=np.int64)
        np.testing.assert_array_equal(self.request('mul', x, 2), [2, 4, 6])
        np.testing.assert_array_equal(self.request('sub', x, 1), [0, 1, 2])

    def test_nullable_decision_is_distinct_from_false(self):
        x = pd.Series([.001, .5, np.nan])
        comparison = self.request('le', x, .01)
        nullable = self.request('cast_boolean', comparison)
        result = self.request('mask_missing', nullable, self.request('isna', x))
        self.assertEqual(str(result.dtype), 'boolean')
        self.assertEqual(result.iloc[:2].tolist(), [True, False])
        self.assertIs(result.iloc[2], pd.NA)

    def test_inclusive_threshold_boundary(self):
        x = pd.Series([np.nextafter(.01, 0), .01, np.nextafter(.01, 1)])
        self.assertEqual(self.request('le', x, .01).tolist(), [True, True, False])

    def test_boolean_three_state_logic(self):
        a = pd.Series([True, False, pd.NA], dtype='boolean')
        b = pd.Series([pd.NA, pd.NA, True], dtype='boolean')
        result = self.request('and', a, b)
        self.assertIs(result.iloc[0], pd.NA)
        self.assertFalse(result.iloc[1])
        self.assertIs(result.iloc[2], pd.NA)

    def test_exact_axes_and_position_are_not_confused(self):
        a = pd.Series([1., 2.], index=['b', 'a'])
        b = pd.Series([1., 2.], index=['a', 'b'])
        with self.assertRaises(bridge.Unsupported): self.request('gt', a, b)
        self.assertEqual(self.request('gt', a, np.array([0., 3.])).tolist(), [True, False])

    def test_integer_text_coercion_and_explicit_decimal_frontier(self):
        result = self.request('numeric', pd.Series(['1', '-2', '+3']))
        self.assertEqual(str(result.dtype), 'int64')
        self.assertEqual(result.tolist(), [1, -2, 3])
        with self.assertRaises(bridge.Unsupported): self.request('numeric', pd.Series(['0.0001']))

    def test_invalid_integer_cast_retains_mature_failure_frontier(self):
        with self.assertRaises(bridge.Unsupported): self.request('cast_int64', pd.Series([np.nan]))
        with self.assertRaises(bridge.Unsupported): self.request('cast_int64', np.array([np.inf]))

    def test_condition_tokens_are_exact(self):
        values = pd.Series(['dominant', 'non_dominant', 'recessive;dominant', 'dominant\n', 'dominant\nx'])
        self.assertEqual(self.request('condition_token', values, 'dominant').tolist(), [True, False, True, True, False])

    def test_typed_string_membership_keeps_all_missing_kinds_false(self):
        values = pd.Series(['MT', None, np.nan, pd.NA, 'M', 'X'], dtype=object)
        result = self.request('isin', values, ['M', 'MT'])
        self.assertEqual(result.tolist(), [True, False, False, False, True, False])
        self.assertEqual(str(result.dtype), 'bool')

    def test_group_max_preserves_sorted_keys_and_nan(self):
        values = pd.Series([.1, .3, np.nan, .2])
        keys = pd.Series(['z', 'a', 'x', 'z'])
        result = self.request('group_max_pairs', values, keys)
        self.assertEqual(list(result), ['a', 'x', 'z'])
        self.assertEqual(result['a'], .3); self.assertEqual(result['z'], .2)
        self.assertTrue(np.isnan(result['x']))

    def test_cache_map_retains_missing_default_and_source_order(self):
        result = self.request('cache_map_zero', pd.Series(['b', 'absent', 'a']), ['a', 'b'], [.1, .2])
        self.assertEqual(result.tolist(), [.2, 0., .1])

    def test_exon_lookup_keeps_two_literal_keys(self):
        result = self.request('exon_map', pd.Series(['f1', 'f1', 'f2']), pd.Series(['1', '2', '1']),
                              ['f1', 'f2'], ['1', '1'], [.01, .02])
        self.assertEqual(result.iloc[0], .01); self.assertTrue(np.isnan(result.iloc[1])); self.assertEqual(result.iloc[2], .02)

    def test_selected_assignment_does_not_change_other_rows(self):
        result = self.request('assign_array', np.array([1, 2, 3]), np.array([False, True, False]), 0)
        np.testing.assert_array_equal(result, [1, 0, 3])

    def test_empty_arrays_and_explicit_none_filter(self):
        result = self.request('zeros', 0)
        self.assertEqual(result.dtype, np.dtype('int64')); self.assertEqual(result.shape, (0,))
        self.assertEqual(self.request('filter_none', [None, True, True]), [False, True, True])


if __name__ == '__main__': unittest.main()
