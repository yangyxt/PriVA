"""Focused native-procedure transport/effect tests, not scientific parity.

No original frequency body, SciPy score or native binary is called here.
Fixtures exercise the actual new wire/host code; scripted native actions are
explicitly protocol fixtures and cannot establish full procedure equivalence.
"""
from pathlib import Path
import struct
import sys
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import priva_rust_frequency_procedure as subject


def action(kind, payload=(), *, sequence=0, token=0, line=263, operation=2):
    effect = kind in subject._EFFECTS
    site = f'gnomAD_rare_AF@{line}:{operation}' if effect or kind == 'unsupported' else ''
    rows = [('action', sequence, kind, token if effect else '', site, subject.PROFILE), *payload]
    if kind != 'closed':
        rows += [('meta', 'native-operations', '0', '', '', ''),
                 ('meta', 'mature-operations', '0', '', '', '')]
    return subject.Action(subject._encode(rows), sequence)


def obj(role, value, ordinal=0):
    return ('object', role, ordinal, 'handle', str(value), '')


def meta(name, value):
    return ('meta', name, value, '', '', '')


class HostProtocolTests(unittest.TestCase):
    def host(self):
        registry = subject.Registry()
        guard = subject.VisitedPackages('/explicitly/absent/binding', '0' * 64)
        return registry, subject.HostEffects(registry, guard, __file__)

    def test_registry_does_not_observe_an_argument(self):
        class Opaque:
            def __repr__(self): raise AssertionError('repr')
            def __eq__(self, other): raise AssertionError('equality')
            def __iter__(self): raise AssertionError('iteration')
        registry = subject.Registry()
        value = Opaque()
        handle = registry.add(value)
        self.assertIs(registry.get(handle), value)
        self.assertIs(registry.get(registry.add(value)), value)

    def test_foreign_and_zero_handles_fail(self):
        registry = subject.Registry()
        for value in (0, 1, 12):
            with self.assertRaises(subject.TransportError): registry.get(value)

    def test_literal_float_bits_and_signed_zero_survive(self):
        registry = subject.Registry()
        for value in ('8000000000000000', '7ff8000000000001', '3ee4f8b588e368f1'):
            self.assertEqual(struct.pack('>d', registry.object('float', value)).hex(), value)

    def test_noncanonical_literals_are_transport_errors(self):
        registry = subject.Registry()
        for tag, value in (('float', '7FF8000000000000'), ('int', '+1'), ('int', '-0'),
                           ('int', '01'), ('none', 'None'), ('handle', '00')):
            with self.assertRaises(subject.TransportError): registry.object(tag, value)

    def test_action_sequence_profile_and_width_are_strict(self):
        samples = [
            [('action', '1', 'closed', '', '', subject.PROFILE)],
            [('action', '0', 'closed', '', '', 'foreign')],
            [('action', '0', 'closed', '', '', subject.PROFILE, '')],
        ]
        for rows in samples:
            with self.assertRaises(subject.TransportError):
                raw = subject._encode(rows) if len(rows[0]) == 6 else b'a\tb\tc\td\te\tf\tg\n'
                subject.Action(raw, 0)

    def test_duplicate_or_extra_metadata_is_rejected(self):
        for payload in ([meta('class', 'float'), meta('class', 'float'), obj('object', 1)],
                        [meta('class', 'float'), meta('extra', 'x'), obj('object', 1)]):
            with self.assertRaises(subject.TransportError): action('isinstance', payload, line=262, operation=0)

    def test_unknown_sites_and_method_names_are_rejected(self):
        with self.assertRaises(subject.TransportError):
            action('getitem', [obj('owner', 1), obj('key', 2)], line=270)
        with self.assertRaises(subject.TransportError):
            action('getattr', [obj('owner', 1), meta('name', 'to_pickle')])

    def test_closed_and_terminal_actions_have_no_effect_fields(self):
        self.assertEqual(action('closed').kind, 'closed')
        with self.assertRaises(subject.TransportError): action('closed', [meta('extra', 'x')])
        with self.assertRaises(subject.TransportError):
            action('return', [obj('result', 1), obj('argument', 1)])

    def test_package_guard_is_lazy_for_other_cutoff_types(self):
        registry, host = self.host()
        handle = registry.add(None)
        with mock.patch.object(subject.columns.guards, 'regular_bytes', side_effect=AssertionError('early source read')):
            first = host.execute(action('isinstance', [obj('object', handle), meta('class', 'float')], line=262, operation=0))
            second = host.execute(action('isinstance', [obj('object', handle), meta('class', 'str')], line=266, operation=1))
        self.assertEqual(first[0][1:3], ('truth', '0'))
        self.assertEqual(second[0][1:3], ('truth', '0'))
        self.assertIsNone(host.packages.binding)

    def test_getter_failure_keeps_the_actual_exception_object(self):
        failure = KeyError(('missing', 3))
        events = []
        class Frame:
            def __getitem__(self, key):
                events.append(key)
                raise failure
        registry, host = self.host()
        owner, key = registry.add(Frame()), registry.add('AN')
        result = host.execute(action('getitem', [obj('owner', owner), obj('key', key)]))
        self.assertEqual(events, ['AN'])
        self.assertEqual(result[0][1], 'raised')
        self.assertIs(registry.get(result[0][2]), failure)
        self.assertIsNone(host.packages.binding)

    def test_getattr_keeps_the_captured_bound_method(self):
        class Column:
            def fillna(self, value): return value
        registry, host = self.host()
        column = Column()
        result = host.execute(action('getattr', [obj('owner', registry.add(column)), meta('name', 'fillna')]))
        method = registry.get(result[0][2])
        self.assertIs(method.__self__, column)
        self.assertIs(method.__func__, Column.fillna)

    def test_where_capture_does_not_evaluate_any_operand(self):
        registry, host = self.host()
        operation = subject.columns._reference.np.where
        result = host.execute(action('bind-where'))
        self.assertIs(registry.get(result[0][2]), operation)
        self.assertIsNone(host.packages.binding)

    def test_builtin_float_inspection_does_not_coerce_subclasses(self):
        class Float(float):
            def __float__(self): raise AssertionError('coercion')
        registry, host = self.host()
        result = host.execute(action('inspect-float', [obj('object', registry.add(Float(.1)))]))
        self.assertEqual(result[0][1:3], ('other-scalar', ''))

    def test_one_mature_call_retains_arguments_and_failure(self):
        registry, host = self.host()
        marker, calls, failure = object(), [], RuntimeError('fixture failure')
        def callback(value):
            calls.append(value)
            raise failure
        result = host.execute(action('mature', [meta('invoke', 'bound'), meta('reason', 'explicit fixture callback'),
            obj('callable', registry.add(callback)), obj('argument', registry.add(marker))]))
        self.assertEqual(len(calls), 1)
        self.assertIs(calls[0], marker)
        self.assertIs(registry.get(result[0][2]), failure)

    def test_callback_transport_named_exception_is_still_the_original_failure(self):
        registry, host = self.host()
        failure = subject.TransportError('actual callback exception')
        def callback(value): raise failure
        result = host.execute(action('mature', [meta('invoke', 'bound'), meta('reason', 'explicit fixture callback'),
            obj('callable', registry.add(callback)), obj('argument', registry.add(None))]))
        self.assertEqual(result[0][1], 'raised')
        self.assertIs(registry.get(result[0][2]), failure)

    def test_code_identity_includes_foreign_filename_and_line_table(self):
        def callback(value): return value
        code = callback.__code__
        foreign = code.replace(co_filename='foreign-context.py')
        self.assertEqual(code, foreign)
        self.assertNotEqual(subject._code_fields(code), subject._code_fields(foreign))
        altered = code.replace(co_linetable=b'')
        self.assertNotEqual(subject._code_fields(code), subject._code_fields(altered))

    def test_original_looking_unknown_call_shape_is_rejected(self):
        with self.assertRaises(subject.TransportError):
            action('mature', [meta('invoke', 'bound'), meta('reason', 'fixture'), obj('callable', 1),
                obj('argument', 2, 0), obj('argument', 3, 1)])

    def test_kernel_operands_cannot_substitute_captured_arguments(self):
        with self.assertRaises(subject.TransportError):
            action('prepare-kernel', [meta('invoke', 'binary'), meta('operator', 'le'), meta('kernel', 'le'),
                obj('left', 1), obj('right', 2), obj('operand', 1, 0), obj('operand', 3, 1)])

    def test_foreign_materialization_fails_before_boxing(self):
        registry, host = self.host()
        payload = [meta('invoke', 'binary'), meta('operator', 'le'), meta('kernel', 'le'),
            obj('left', 1), obj('right', 2), obj('operand', 1, 0), obj('operand', 2, 1),
            ('native', 'array', 'bool', '1', '', ''), ('cell', '0', '0', 'bool', '1', '')]
        with mock.patch.object(subject.columns, '_result', side_effect=AssertionError('foreign boxing')):
            with self.assertRaises(subject.TransportError): host.execute(action('materialize', payload))

    def test_inherited_owned_child_cleanup_is_reused(self):
        self.assertIs(subject.ProcedureSession._finish_child, subject.columns.NativeSession._finish_child)
        self.assertIs(subject.ProcedureSession.close, subject.columns.NativeSession.close)
        self.assertIs(subject.ProcedureSession._write, subject.columns.NativeSession._write)
        self.assertIs(subject.ProcedureSession._read, subject.columns.NativeSession._read)

    def test_column_operation_cannot_enter_procedure_session(self):
        session = subject.ProcedureSession('/absent/ELF', '0' * 64)
        with self.assertRaises(subject.TransportError): session.request('fillna', [], site='gnomAD_rare_AF@263')
        self.assertIsNone(session.process)


if __name__ == '__main__':
    unittest.main()
