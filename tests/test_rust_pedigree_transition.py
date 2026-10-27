"""Focused v2 key/cache/envelope controls; fake wire implements no science.

Actual native kernel tests are in transition.rs and require the separately
reviewed affected build. These controls do not qualify full-file parity.
"""
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import priva_rust_pedigree_filtration as bridge
from priva_rust_pedigree_transition import TransitionSession, VERSION, _Memo, _Profile


def context():
    return {'patients.count': '1', 'patients.0': 'p', 'controls.count': '0', 'proband': 'p',
        'father': '', 'mother': '', 'sex': '1', 'father.original_kind': 'none',
        'mother.original_kind': 'none', 'sex.original_kind': 'str'}


def initial_state():
    return {'phase': 'Idle', 'index': '0', 'chrom_type': '',
        'proband.kind': 'none', 'father.kind': 'none', 'mother.kind': 'none'}


class FakeWire:
    """Programmed transport only: returns supplied state and a fixed outcome."""
    def __init__(self, binary, *, timeout=300, protocol=None):
        assert protocol == VERSION
        self.path, self.identity, self.sha256 = Path(binary), ('fixture-only',), 'fixture-only'
        self.starts = 1
        self.sequence = 1
        self.operations = {'hello': 1}
        self.closed = False
        self.status = None
        self.calls = []
        self.reply = None
        read, self.write = os.pipe()
        self.process = SimpleNamespace(stdout=os.fdopen(read, 'rb', buffering=0), poll=lambda: self.status)

    def request(self, operation, **fields):
        self.sequence += 1
        self.operations[operation] = self.operations.get(operation, 0) + 1
        self.calls.append((operation, fields))
        if self.reply is not None:
            return self.reply(operation, fields)
        if operation == 'kernel-init':
            return dict(value='ready', **fields, **{'state.' + key: value for key, value in initial_state().items()})
        if operation == 'transition':
            return {**{key: value for key, value in fields.items() if key.startswith('state.')},
                    'outcome.status': 'ok', 'outcome.value': 'true'}
        if operation == 'quit':
            return {'value': 'closed'}
        raise AssertionError('unplanned fake operation')

    def close(self, *, suppress=False, _request=None):
        if self.closed:
            return
        try:
            if self.status is None:
                (_request or self.request)('quit')
        except bridge.NativeProtocolError:
            if not suppress:
                raise
        finally:
            self.closed = True
            self.process.stdout.close()
            os.close(self.write)


def executable(root, mode):
    """Real child/pipe fixture; fixed metadata reply is NOT a scientific kernel."""
    path = root / ('transport-' + mode)
    text = r'''#!PYTHON
import csv,io,os,struct,sys,time
MODE=MODE_LITERAL
def read():
 h=sys.stdin.buffer.read(8)
 if not h:return None
 n,=struct.unpack('>Q',h)
 return dict(list(csv.reader(io.StringIO(sys.stdin.buffer.read(n).decode())))[1:])
def write(rows):
 s=io.StringIO(newline='');w=csv.writer(s,lineterminator='\n');w.writerow(['name','value']);w.writerows(rows)
 b=s.getvalue().encode();sys.stdout.buffer.write(struct.pack('>Q',len(b))+b);sys.stdout.buffer.flush()
while True:
 q=read()
 if q is None:break
 rows=[('status','ok'),('sequence',q['sequence'])]
 if q['op']=='hello':rows.append(('value',q['protocol']))
 elif q['op']=='transition':
  if MODE=='hang':time.sleep(30)
  if MODE=='corrupt':rows.extend([('outcome.value','true'),('outcome.value','duplicate')])
  elif MODE=='state':rows.extend([('state.foreign','bad'),('outcome.status','ok'),('outcome.value','true')])
  elif MODE=='source-error':rows.extend([('outcome.status','source-error'),('outcome.class','TypeError'),('outcome.message','source fixture error')])
  else:rows.extend([('outcome.status','ok'),('outcome.value','true')])
 elif q['op']=='quit':rows.append(('value','closed'))
 else:sys.exit(38)
 write(rows)
 if q['op']=='transition' and MODE=='dead':sys.exit(37)
 if q['op']=='quit':
  if MODE=='extra':sys.stdout.buffer.write(b'extra');sys.stdout.buffer.flush()
  if MODE=='nonzero':sys.exit(37)
  break
'''.replace('PYTHON', sys.executable, 1).replace('MODE_LITERAL', repr(mode), 1)
    path.write_text(text)
    path.chmod(0o700)
    return path


class TransitionCacheControls(unittest.TestCase):
    def session(self):
        with patch.object(bridge, 'NativeSession', FakeWire):
            client = TransitionSession(bridge, sys.executable)
        self.addCleanup(client.close, suppress=True)
        return client

    def test_reuse_is_distinct_from_native_execution_and_returns_fresh_values(self):
        client = self.session()
        first = client.request('genotype', predicate='has', **{'gt.kind': 'none'})
        first.pop('value')
        self.assertEqual(client.request('genotype', predicate='has', **{'gt.kind': 'none'}), {'value': 'true'})
        self.assertEqual(client.native_transitions, {'genotype': 1})
        self.assertEqual(client.cache_hits, {'genotype': 1})
        self.assertEqual(client.completed_operations['genotype'], 2)
        self.assertEqual(client.wire.operations['transition'], 1)

    def test_complete_typed_message_distinguishes_container_allele_and_value(self):
        client = self.session()
        for container in ('tuple', 'list'):
            for kind, value in (('bool', '1'), ('int', '1'), ('str', '1'), ('int', '2')):
                fields = {'gt.kind': container, 'gt.count': '1', 'gt.0.kind': kind, 'gt.0.value': value}
                client.request('genotype', predicate='has', **fields)
                client.request('genotype', predicate='has', **fields)
        self.assertEqual(client.native_transitions, {'genotype': 8})
        self.assertEqual(client.cache_hits, {'genotype': 8})

    def test_context_state_and_source_profile_are_complete_key_components(self):
        client = self.session()
        client.request('filter-init', **context())
        for phase in ('Idle', 'Done'):
            state = initial_state();state['phase'] = phase
            client._state = tuple(sorted(state.items()))
            client.request('genotype', predicate='has', **{'gt.kind': 'none'})
        client._context = tuple(sorted(dict(context(), **{'father.original_kind': 'str'}).items()))
        client.request('genotype', predicate='has', **{'gt.kind': 'none'})
        namespace = {};binding = object()
        client.admit_namespace(namespace, (binding,))
        client.request('genotype', predicate='has', **{'gt.kind': 'none'})
        client.request('genotype', predicate='has', **{'gt.kind': 'none'})
        self.assertEqual(client.native_transitions, {'genotype': 4})
        self.assertEqual(client.cache_hits, {'genotype': 1})
        self.assertIs(client._namespace_profiles[id(namespace)][0], namespace)

    def test_stored_raw_gt_type_is_part_of_state_key(self):
        client = self.session();client.request('filter-init', **context())
        for kind in ('tuple', 'list'):
            state = initial_state()
            state.update({'proband.kind': kind, 'proband.count': '1', 'proband.0.kind': 'int', 'proband.0.value': '1'})
            client._state = tuple(sorted(state.items()))
            client.request('genotype', predicate='has', **{'gt.kind': 'none'})
        self.assertEqual(client.native_transitions, {'genotype': 2})

    def test_changed_bindings_refuse_reuse_before_a_new_wire_request(self):
        client = self.session();namespace = {};first = object()
        client.admit_namespace(namespace, (first,))
        client.request('genotype', predicate='has', **{'gt.kind': 'none'})
        count = client.wire.sequence
        with self.assertRaisesRegex(bridge.NativeProtocolError, 'profile changed'):
            client.admit_namespace(namespace, (object(),))
        self.assertEqual(client.wire.sequence, count)

    def test_source_errors_are_fresh_and_state_is_applied_before_raise(self):
        client = self.session();client.request('filter-init', **context())
        before = initial_state();before['phase'] = 'Proband';client._state = tuple(sorted(before.items()))
        after = dict(before, **{'proband.kind': 'int', 'proband.value': '7'})
        def error_reply(operation, fields):
            return {**{'state.' + key: value for key, value in after.items()}, 'outcome.status': 'source-error',
                'outcome.class': 'TypeError', 'outcome.message': 'same source error'}
        client.wire.reply = error_reply
        errors = []
        for _ in range(2):
            client._state = tuple(sorted(before.items()))
            try:client.request('filter-reply', **{'gt.kind': 'int', 'gt.value': '7'})
            except TypeError as error:errors.append(error)
            self.assertEqual(dict(client._state), after)
        self.assertEqual([error.args for error in errors], [('same source error',)] * 2)
        self.assertIsNot(errors[0], errors[1])
        self.assertEqual(client.native_transitions, {'filter-reply': 1})
        self.assertEqual(client.cache_hits, {'filter-reply': 1})
        self.assertEqual(client.completed_operations['filter-reply'], 2)
        client.wire.reply = None

    def test_record_order_is_checked_on_each_begin_before_cache(self):
        client = self.session();client.request('filter-init', **context())
        def begin(operation, fields):
            state = initial_state();state.update(phase='DebugDrop')
            return {**{'state.' + key: value for key, value in state.items()},
                'outcome.status': 'ok', 'outcome.action': 'debug-no-patient'}
        client.wire.reply = begin
        client.request('filter-begin', ordinal='1')
        client._state = tuple(sorted(initial_state().items()))
        client.request('filter-begin', ordinal='2')
        self.assertEqual(client.cache_hits, {'filter-begin': 1})
        self.assertEqual(client._ordinal, 2)
        before = client.wire.sequence
        with self.assertRaises(bridge.NativeProtocolError):client.request('filter-begin', ordinal='2')
        self.assertEqual(client.wire.sequence, before)

    def test_lru_entry_and_byte_bounds_and_oversize_never_fallback(self):
        memo = _Memo(entries=2, budget=10000)
        memo.put(('a',), ('one',));memo.put(('b',), ('two',))
        memo.get(('a',));memo.put(('c',), ('three',))
        self.assertIsNone(memo.get(('b',)))
        self.assertEqual(memo.evictions, 1)
        self.assertLessEqual(memo.bytes, memo.budget)
        memo = _Memo(entries=2, budget=1);memo.put(('a',), ('one',))
        self.assertEqual(memo.oversize, 1);self.assertFalse(memo.values)
        profile = _Profile(('complete',));self.assertRaises(AttributeError, setattr, profile, 'contents', ())

    def test_oversize_session_entry_runs_native_each_time(self):
        client = self.session();client._memo = _Memo(budget=1)
        for _ in range(2):client.request('genotype', predicate='has', **{'gt.kind': 'none'})
        self.assertEqual(client.native_transitions, {'genotype': 2})
        self.assertFalse(client.cache_hits)
        self.assertEqual(client._memo.oversize, 2)

    def test_invalid_complete_state_is_not_cached_and_poisoned(self):
        client = self.session()
        client.wire.reply = lambda operation, fields: {'state.foreign': 'bad', 'outcome.status': 'ok', 'outcome.value': 'true'}
        with self.assertRaises(bridge.NativeProtocolError):client.request('genotype', predicate='has', **{'gt.kind': 'none'})
        self.assertFalse(client._memo.values)
        with self.assertRaises(bridge.NativeProtocolError):client.request('genotype', predicate='has', **{'gt.kind': 'none'})
        self.assertEqual(client.wire.sequence, 2)

    def test_canonical_zero_missing_and_arbitrary_integer_tags_validate(self):
        client = self.session();client.request('filter-init', **context())
        for kind, value in (('bool', '0'), ('int', '0'), ('int', '-2'), ('int', '9' * 5000), ('str', 'α\nA')):
            state = initial_state()
            state.update({'proband.kind': 'list', 'proband.count': '2', 'proband.0.kind': kind,
                'proband.0.value': value, 'proband.1.kind': 'none'})
            self.assertEqual(dict(client._validate_state(state)), state)
        for value in ('-0', '+1', '01', ''):
            state = initial_state();state.update({'proband.kind': 'int', 'proband.value': value})
            with self.assertRaises(bridge.NativeProtocolError):client._validate_state(state)


class ActualTransportControls(unittest.TestCase):
    def make(self, root, mode, timeout=1):
        client = TransitionSession(bridge, executable(root, mode), timeout=timeout)
        self.addCleanup(client.close, suppress=True)
        return client

    def test_dead_child_is_checked_even_on_a_cached_key(self):
        with tempfile.TemporaryDirectory() as folder:
            client = self.make(Path(folder), 'dead')
            client.request('genotype', predicate='has', **{'gt.kind': 'none'})
            client.wire.process.wait(timeout=2)
            with self.assertRaises(bridge.NativeProtocolError):client.request('genotype', predicate='has', **{'gt.kind': 'none'})
            self.assertFalse(client.cache_hits)

    def test_corrupt_miss_state_and_timeout_fail_closed_without_an_entry(self):
        with tempfile.TemporaryDirectory() as folder:
            for mode in ('corrupt', 'state', 'hang'):
                with self.subTest(mode=mode):
                    client = self.make(Path(folder), mode, timeout=0.1)
                    with self.assertRaises(bridge.NativeProtocolError):client.request('genotype', predicate='has', **{'gt.kind': 'none'})
                    self.assertFalse(client._memo.values)
                    client.close(suppress=True)

    def test_extra_output_and_nonzero_close_remain_failures(self):
        with tempfile.TemporaryDirectory() as folder:
            for mode in ('extra', 'nonzero'):
                with self.subTest(mode=mode):
                    client = self.make(Path(folder), mode)
                    with self.assertRaises(bridge.NativeProtocolError):client.close()

    def test_real_pipe_source_error_is_reused_as_a_fresh_exception(self):
        with tempfile.TemporaryDirectory() as folder:
            client = self.make(Path(folder), 'source-error')
            errors = []
            for _ in range(2):
                try:client.request('genotype', predicate='has', **{'gt.kind': 'none'})
                except TypeError as error:errors.append(error)
            self.assertEqual([error.args for error in errors], [('source fixture error',)] * 2)
            self.assertIsNot(errors[0], errors[1])
            self.assertEqual(client.wire.operations['transition'], 1)
            self.assertEqual(client.cache_hits, {'genotype': 1})
            client.close()


if __name__ == '__main__':
    unittest.main()
