"""Focused source/transport controls; no biological inputs or model execution.

Actual native parity is a separately pinned saved-oracle procedure. These unit
controls cannot stand in for that execution or activate a default backend.
"""
import ast
import copy
import csv
import dis
import importlib.util
import io
import logging
import operator
from pathlib import Path
import struct
import sys
from types import FunctionType, SimpleNamespace
import unittest
from unittest.mock import patch
import warnings

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import numpy as np
import pandas as pd
import acmg_pp3_bp4_bp7_insilico as reference
import priva_rust_insilico as bridge
import priva_rust_insilico_source as guards


def result_bytes(kind, dtype, cells, labels=(), sequence=0):
    stream = io.StringIO(newline='')
    writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
    writer.writerow(('result', sequence, kind, dtype, len(cells), bridge.PROFILE))
    for i, label in enumerate(labels): writer.writerow(('index', 0, i, 'label', label, ''))
    for i, value in enumerate(cells): writer.writerow(('cell', 0, i, *bridge._cell(value), ''))
    return stream.getvalue().encode()


class InsilicoWireTests(unittest.TestCase):
    def test_float_bits_and_quoted_unicode_are_lossless(self):
        for value in (0.0, -0.0, float('inf'), float('-inf'), float('nan')):
            tag, raw = bridge._cell(value)
            self.assertEqual(struct.pack('>d', bridge._uncell(tag, raw)), struct.pack('>d', value))
        value = pd.Series(['a\tb', 'λ\n"x"'], name='quoted')
        rows = list(csv.reader(io.StringIO(bridge._request(0, 'as_str', [value]).decode()), delimiter='\t'))
        self.assertEqual([row[4] for row in rows if row[0] == 'cell'], value.tolist())

    def test_response_denominator_and_sequence_cannot_be_waived(self):
        raw = result_bytes('array', 'int64', [0, 1])
        np.testing.assert_array_equal(bridge._result(raw, 0, [2], 'zeros'), [0, 1])
        with self.assertRaises(bridge.TransportError): bridge._result(raw, 1, [2], 'zeros')
        bad = raw.replace(b'int64\t2\t', b'int64\t3\t', 1)
        with self.assertRaises(bridge.TransportError): bridge._result(bad, 0, [2], 'zeros')

    def test_nonconsecutive_and_extra_result_rows_are_rejected(self):
        raw = result_bytes('array', 'bool', [False]).replace(b'cell\t0\t0', b'cell\t0\t1')
        with self.assertRaises(bridge.TransportError): bridge._result(raw, 0, [], 'eq')
        with self.assertRaises(bridge.TransportError): bridge._result(result_bytes('scalar', 'scalar', [1, 2]), 0, [], 'sum')

    def test_unknown_objects_are_not_stringified_or_coerced(self):
        class Object:
            def __str__(self): raise AssertionError('must not call')
        with self.assertRaises(bridge.Unsupported): bridge._cell(Object())

    def test_signed_index_encoding_preserves_actual_order(self):
        values = [-20, -1, 0, 2, 10, 100]
        self.assertEqual([bridge._unlabel(key) for key in sorted(map(bridge._label, values))], values)

    def test_aligned_series_and_positional_array_keep_distinct_wire_kinds(self):
        pvs = pd.Series([0, 3], index=['b', 'a'])
        series = bridge._value(pvs)
        array = bridge._value(np.asarray(pvs))
        self.assertEqual(series[:3], ('series', 'int64', ['sb', 'sa']))
        self.assertEqual(array[:3], ('array', 'int64', []))

    def test_corrupt_response_types_are_always_transport_errors(self):
        with self.assertRaises(bridge.TransportError): bridge._uncell('int', 'not-an-int')
        for raw in (b'\xff', b'"incomplete', result_bytes('array', 'unknown_dtype', [1])):
            with self.assertRaises(bridge.TransportError): bridge._result(raw, 0, [1], 'zeros')
        session = bridge.NativeSession('unused', '0' * 64)
        runtime = bridge._Runtime(session, reference)
        def corrupt(*_args): return bridge._result(b'\xff', 0, [True, False], 'and')
        with patch.object(session, 'request', side_effect=corrupt), patch.object(runtime, 'fallback') as fallback:
            with self.assertRaises(bridge.TransportError): runtime.binary(32, 'PP3_BP4_criteria', 'and', True, False)
            fallback.assert_not_called()

    def test_result_kind_and_foreign_axis_are_rejected(self):
        with self.assertRaises(bridge.TransportError):
            bridge._result(result_bytes('array', 'bool', [True]), 0, [pd.Series([1.])], 'notna')
        with self.assertRaises(bridge.TransportError):
            bridge._result(result_bytes('series', 'bool', [True], ['sforeign']), 0, [pd.Series([1.])], 'notna')

    def test_partial_writes_are_completed_without_replaying_a_frame(self):
        session = bridge.NativeSession('unused', '0' * 64)
        session.process = SimpleNamespace(poll=lambda: None, stdin=SimpleNamespace(fileno=lambda: 42))
        pieces = []
        def write(_fd, buffer):
            piece = bytes(buffer[:3]); pieces.append(piece); return len(piece)
        with patch.object(bridge.time, 'monotonic', return_value=0), patch.object(bridge.select, 'select', return_value=([], [True], [])), patch.object(bridge.os, 'write', side_effect=write):
            session._write(b'0123456789', 1)
        self.assertEqual(b''.join(pieces), b'0123456789')
        self.assertEqual(list(map(len, pieces)), [3, 3, 3, 1])

    def test_no_progress_write_and_deadline_fail_closed(self):
        session = bridge.NativeSession('unused', '0' * 64)
        session.process = SimpleNamespace(poll=lambda: None, stdin=SimpleNamespace(fileno=lambda: 42))
        with patch.object(bridge.time, 'monotonic', return_value=0), patch.object(bridge.select, 'select', return_value=([], [True], [])), patch.object(bridge.os, 'write', return_value=0):
            with self.assertRaises(bridge.TransportError): session._write(b'x', 1)
        with patch.object(bridge.time, 'monotonic', return_value=2), patch.object(bridge.os, 'write') as write:
            with self.assertRaises(bridge.TransportError): session._write(b'x', 1)
            write.assert_not_called()


class SourceFrontierTests(unittest.TestCase):
    def test_transform_leaves_all_getters_loggers_and_branch_statements(self):
        original = guards.source_tree(Path(reference.__file__).read_bytes())
        changed = [bridge._NativeCalls().visit(copy.deepcopy(node)) for node in original]
        def frontiers(nodes):
            getters, logs, branches = [], [], []
            for node in ast.walk(ast.Module(body=nodes, type_ignores=[])):
                if isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name) and node.value.id in ('df', 'raw_values'):
                    getters.append((node.lineno, ast.dump(node, include_attributes=False)))
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and isinstance(node.func.value, ast.Name) and node.func.value.id == 'logger':
                    logs.append((node.lineno, node.func.attr))
                if isinstance(node, (ast.If, ast.Return)):
                    branches.append((type(node).__name__, node.lineno))
            return sorted(getters), sorted(logs), sorted(branches)
        self.assertEqual(frontiers(original), frontiers(changed))
        class RuntimeOperations(ast.NodeVisitor):
            def __init__(self): self.operations = []
            def visit_FunctionDef(self, node):
                for statement in node.body: self.visit(statement)
            def visit_BinOp(self, node):
                if type(node.op) in (ast.BitAnd, ast.BitOr): self.operations.append(node)
                self.generic_visit(node)
        before, after = RuntimeOperations(), RuntimeOperations()
        before.visit(ast.Module(body=original, type_ignores=[]))
        after.visit(ast.Module(body=changed, type_ignores=[]))
        self.assertTrue(before.operations)
        self.assertEqual(after.operations, [])
        def definitions(nodes):
            return [(node.name, ast.dump(node.args, include_attributes=False),
                     ast.dump(node.returns, include_attributes=False) if node.returns else None,
                     [ast.dump(value, include_attributes=False) for value in node.decorator_list])
                    for node in ast.walk(ast.Module(body=nodes, type_ignores=[]))
                    if type(node) is ast.FunctionDef]
        self.assertEqual(definitions(original), definitions(changed))

    def test_actual_nested_annotation_keeps_type_and_definition_timing(self):
        source = guards.source_tree(Path(reference.__file__).read_bytes())
        nested = next(node for function in source for node in ast.walk(function)
                      if type(node) is ast.FunctionDef and node.name == '_coerce_score_column')
        session = bridge.NativeSession('unused', '0' * 64)
        runtime = bridge._Runtime(session, reference)
        observations = []
        for lowered in (False, True):
            tree = ast.parse("def annotation_fixture():\n events.append('before')\n events.append('after')\n return _coerce_score_column\n")
            tree.body[0].body.insert(1, copy.deepcopy(nested))
            if lowered: tree = bridge._NativeCalls().visit(tree)
            namespace = dict(vars(reference), events=[])
            namespace[bridge._GLOBAL] = runtime
            exec(compile(ast.fix_missing_locations(tree), '<annotation-definition-fixture>', 'exec'), namespace)
            function = namespace['annotation_fixture']
            instructions = {item.offset: item for item in dis.get_instructions(function)}
            def trace(frame, event, _argument):
                if frame.f_code is function.__code__:
                    frame.f_trace_opcodes = True
                    instruction = instructions.get(frame.f_lasti)
                    if event == 'opcode' and instruction is not None and instruction.opname == 'BINARY_OP' and instruction.argrepr == '|':
                        namespace['events'].append('annotation-union')
                    return trace
            previous = sys.gettrace()
            caller = sys._getframe()
            previous_opcodes = caller.f_trace_opcodes
            try:
                with patch.object(session, 'request') as request, patch.object(runtime, 'fallback') as fallback:
                    # CPython 3.12 requires opcode monitoring primed before
                    # installing the hook; otherwise the first call can miss it.
                    caller.f_trace_opcodes = True
                    sys.settrace(trace)
                    defined = function()
                    request.assert_not_called()
                    fallback.assert_not_called()
            finally:
                sys.settrace(previous)
                caller.f_trace_opcodes = previous_opcodes
            self.assertEqual(defined.__annotations__['missing_default'], float | None)
            self.assertIs(defined.__annotations__['return'], pd.Series)
            observations.append((namespace['events'], defined.__annotations__, defined.__defaults__))
        self.assertEqual(observations[0], observations[1])
        self.assertEqual(observations[1][0], ['before', 'annotation-union', 'after'])

    def test_annotation_failures_stop_definition_before_runtime_table_operations(self):
        failure = RuntimeError('annotation failure')
        for fail_at in ('argument', 'return', None):
            observations = []
            for lowered in (False, True):
                events = []
                class Annotation:
                    def __init__(self, label): self.label = label
                    def __or__(self, _other):
                        events.append(self.label)
                        if self.label == fail_at: raise failure
                        return self.label
                runtime = SimpleNamespace(binary=lambda *args: events.append(('runtime', args[2])) or args[-2])
                namespace = dict(argument=Annotation('argument'), result=Annotation('return'), events=events)
                namespace[bridge._GLOBAL] = runtime
                tree = ast.parse("def definition(left: argument | None, right=None) -> result | None:\n events.append('body')\n return (left | right) & left\n")
                if lowered: tree = bridge._NativeCalls().visit(tree)
                try:
                    exec(compile(ast.fix_missing_locations(tree), '<annotation-failure-fixture>', 'exec'), namespace)
                except RuntimeError as error:
                    self.assertIs(error, failure)
                    self.assertNotIn('definition', namespace)
                else:
                    self.assertIsNone(fail_at)
                    self.assertEqual(namespace['definition'].__annotations__, {'left': 'argument', 'return': 'return'})
                    if lowered:
                        namespace['definition'](True, False)
                        self.assertEqual(events[-3:], ['body', ('runtime', 'or'), ('runtime', 'and')])
                        events = events[:-3]
                observations.append(events)
            self.assertEqual(observations[0], observations[1])
            expected = ['argument'] if fail_at == 'argument' else ['argument', 'return']
            self.assertEqual(observations[0], expected)

    def test_default_route_calls_original_once_without_child(self):
        marker = object()
        with patch.object(reference, 'BP7_criteria', return_value=marker) as original, patch.object(bridge.subprocess, 'Popen') as child:
            self.assertIs(bridge.BP7_criteria(marker), marker)
            original.assert_called_once_with(marker)
            child.assert_not_called()

    def test_source_unknown_whole_profile_falls_back_before_child(self):
        marker = object()
        call = bridge.Bridge(binary='absent', binary_sha256='0' * 64)
        with patch.object(bridge, 'transformed', side_effect=ValueError('unsupported test profile')), patch.object(reference, 'BP7_criteria', return_value=marker) as original, patch.object(bridge.subprocess, 'Popen') as child:
            self.assertIs(call.call('BP7_criteria', marker), marker)
            original.assert_called_once_with(marker)
            child.assert_not_called()
            self.assertEqual(call.last_metrics['mode'], 'original')

    def test_transport_fault_never_calls_mature_fallback(self):
        session = bridge.NativeSession('unused', '0' * 64)
        runtime = bridge._Runtime(session, reference)
        with patch.object(session, 'request', side_effect=bridge.TransportError('dead')), patch.object(runtime, 'fallback') as fallback:
            with self.assertRaises(bridge.TransportError): runtime.binary(32, 'PP3_BP4_criteria', 'and', True, False)
            fallback.assert_not_called()

    def test_explicit_unsupported_operation_executes_once_with_same_objects(self):
        session = bridge.NativeSession('unused', '0' * 64)
        runtime = bridge._Runtime(session, reference)
        left, right = object(), object()
        seen = []
        def original(a, b): seen.append((a, b)); return a
        with patch.object(session, 'request', side_effect=bridge.Unsupported('type')):
            self.assertIs(runtime.calculate(32, 'PP3_BP4_criteria', 'eq', [left, right], original), left)
        self.assertEqual(len(seen), 1)
        self.assertIs(seen[0][0], left)
        self.assertIs(seen[0][1], right)

    def test_mature_failure_is_not_retried(self):
        session = bridge.NativeSession('unused', '0' * 64)
        runtime = bridge._Runtime(session, reference)
        failure = KeyError('original key')
        seen = []
        def original(value): seen.append(value); raise failure
        with patch.object(session, 'request', side_effect=bridge.Unsupported('type')):
            with self.assertRaises(KeyError) as caught: runtime.calculate(133, 'BP7_criteria', 'eq', [True], original)
        self.assertIs(caught.exception, failure)
        self.assertEqual(seen, [True])


class CompatibilityFrontierTests(unittest.TestCase):
    """Real mature ownership effects, without a criterion or ELF invocation."""

    def runtime(self):
        session = bridge.NativeSession('unused', '0' * 64)
        return session, bridge._Runtime(session, reference)

    def test_floating_boolean_fill_preserves_actual_dtype_warning_and_failure(self):
        for dtype in ('float64', 'Float64'):
            for replacement in (False, np.bool_(False)):
                observations = []
                for lowered in (False, True):
                    receiver = pd.Series([1.0, float('nan')], dtype=dtype)
                    session, runtime = self.runtime()
                    with warnings.catch_warnings(record=True) as emitted:
                        warnings.simplefilter('always')
                        with patch.object(session, 'request') as request:
                            try:
                                result = (runtime.method(77, 'PP3_BP4_criteria', receiver, 'fillna', replacement)
                                          if lowered else receiver.fillna(replacement))
                                outcome = ('returned', str(result.dtype), result.tolist())
                            except (ValueError, TypeError) as error:
                                outcome = ('raised', type(error), str(error))
                            request.assert_not_called()
                    observations.append((outcome, [(type(item.message), str(item.message)) for item in emitted]))
                    if lowered:
                        self.assertEqual(runtime.session.metrics['mature_fallbacks'], {'fillna': 1})
                self.assertEqual(observations[0], observations[1])

    def test_numeric_floating_fill_stays_native(self):
        receiver = pd.Series([1.0, float('nan')])
        marker = object()
        session, runtime = self.runtime()
        with patch.object(session, 'request', return_value=marker) as request:
            self.assertIs(runtime.method(62, '_coerce_score_column', receiver, 'fillna', 0), marker)
        self.assertIs(request.call_args.args[1][0], receiver)
        self.assertEqual(request.call_args.args[1][1], 0)
        self.assertFalse(session.metrics['mature_fallbacks'])

    def test_LoF_warning_frontier_matches_actual_CALL_and_co_lines(self):
        instructions = list(dis.get_instructions(reference.PP3_BP4_criteria))
        ordinal = next(i for i, item in enumerate(instructions)
                       if item.opname == 'LOAD_ATTR' and item.argval == 'fillna'
                       and item.positions.lineno == 39)
        call = next(item for item in instructions[ordinal:] if item.opname == 'CALL')
        self.assertEqual(call.positions.lineno, 39)
        self.assertTrue(any(start <= call.offset < end and line == 39
                            for start, end, line in reference.PP3_BP4_criteria.__code__.co_lines()))
        tree = ast.Module(body=copy.deepcopy(guards.source_tree(Path(reference.__file__).read_bytes())), type_ignores=[])
        lowered = bridge._NativeCalls().visit(tree)
        sites = [node.args[0].value for node in ast.walk(lowered)
                 if type(node) is ast.Call and type(node.func) is ast.Attribute
                 and node.func.attr == 'method' and len(node.args) > 3
                 and type(node.args[3]) is ast.Constant and node.args[3].value == 'fillna']
        self.assertIn(39, sites)
        session, runtime = self.runtime()
        with warnings.catch_warnings(record=True) as emitted:
            warnings.simplefilter('always')
            with patch.object(session, 'request', side_effect=bridge.Unsupported('object downcast')):
                runtime.method(39, 'PP3_BP4_criteria', pd.Series([1, 2], dtype=object), 'fillna', '')
        self.assertTrue(emitted)
        self.assertTrue(all(item.filename == reference.__file__ and item.lineno == 39 for item in emitted))

    def test_assignment_uses_actual_NumPy_mask_ownership_at_its_frontier(self):
        def fixture():
            index = pd.Index(list('abcdefg'), name='synthetic')
            return np.zeros(7, dtype=int), pd.Series([True, False, True, False, False, True, False], index=index)
        expected, old_mask = fixture()
        expected[old_mask] = 1
        actual, mask = fixture()
        session, runtime = self.runtime()
        received = []
        def request(operation, args):
            self.assertEqual(operation, 'assign')
            received.append(args)
            array, positional, value = args
            self.assertIs(type(positional), np.ndarray)
            self.assertTrue(np.shares_memory(positional, mask._mgr.blocks[0].values))
            result = array.copy(); result[positional] = value
            return result
        with patch.object(session, 'request', side_effect=request):
            result = runtime.assign(126, 'PP3_BP4_criteria', actual, mask, 1)
        np.testing.assert_array_equal(result, expected)
        self.assertEqual(len(received), 1)
        self.assertEqual([(key, type(value)) for key, value in vars(mask.index).get('_cache', {}).items()],
                         [(key, type(value)) for key, value in vars(old_mask.index).get('_cache', {}).items()])
        self.assertFalse(session.metrics['mature_fallbacks'])

    def test_assignment_unsupported_does_not_repeat_Series_conversion(self):
        mask = pd.Series([True, False], index=pd.Index(['a', 'b'], name='synthetic'))
        array = np.zeros(2, dtype=int)
        session, runtime = self.runtime()
        convert = np.asarray
        with patch.object(bridge.np, 'asarray', wraps=convert) as ownership:
            with patch.object(session, 'request', side_effect=bridge.Unsupported('assignment')):
                result = runtime.assign(126, 'PP3_BP4_criteria', array, mask, 1)
        self.assertIs(result, array)
        np.testing.assert_array_equal(result, [1, 0])
        # Series.__array__ itself calls np.asarray on its raw ndarray. Count
        # the actual Series ownership frontier, not every nested use of the
        # global NumPy symbol; a repeated Series conversion is still refused.
        self.assertEqual([type(call.args[0]) for call in ownership.call_args_list],
                         [pd.Series, np.ndarray])
        self.assertIs(ownership.call_args_list[0].args[0], mask)
        self.assertIs(ownership.call_args_list[1].args[0], mask._mgr.blocks[0].values)
        self.assertEqual(session.metrics['mature_fallbacks'], {'assign': 1})

    def test_different_Series_axes_use_actual_alignment_once_with_same_operands(self):
        def fixture():
            return (pd.Series([True, False, True], index=pd.RangeIndex(3)),
                    pd.Series([False, True, False], index=pd.RangeIndex(2, -1, -1)))
        old_left, old_right = fixture()
        expected = old_left & old_right
        left, right = fixture()
        session, runtime = self.runtime()
        seen = []
        def original(a, b):
            seen.append((a, b))
            return operator.and_(a, b)
        with patch.object(session, 'request') as request:
            result = runtime.calculate(113, 'PP3_BP4_criteria', 'and', [left, right], original)
            request.assert_not_called()
        self.assertEqual(len(seen), 1)
        self.assertIs(seen[0][0], left); self.assertIs(seen[0][1], right)
        pd.testing.assert_series_equal(result, expected)
        for current, old in ((left, old_left), (right, old_right)):
            self.assertEqual(vars(current.index).get('_cache', {}), vars(old.index).get('_cache', {}))
        self.assertEqual(session.metrics['mature_fallbacks'], {'and': 1})

    def test_matching_Series_axes_stay_native_without_alignment_caches(self):
        left = pd.Series([True, False], index=pd.RangeIndex(2))
        right = pd.Series([False, True], index=pd.RangeIndex(2))
        marker = object()
        session, runtime = self.runtime()
        with patch.object(session, 'request', return_value=marker) as request:
            self.assertIs(runtime.binary(113, 'PP3_BP4_criteria', 'or', left, right), marker)
        request.assert_called_once()
        self.assertEqual(vars(left.index).get('_cache', {}), {})
        self.assertEqual(vars(right.index).get('_cache', {}), {})
        self.assertFalse(session.metrics['mature_fallbacks'])


class PackageBindingTests(unittest.TestCase):
    """Canonical declarations, including slots outside the shared PM profile."""

    def test_canonical_current_bindings_are_admitted_without_original_call(self):
        self.assertEqual(guards.reference_guard(reference), Path(reference.__file__).read_bytes())

    def test_foreign_to_numeric_is_not_accepted_against_its_own_source(self):
        with patch.object(pd, 'to_numeric', lambda *args, **kwargs: None):
            with self.assertRaises(ValueError): guards.reference_guard(reference)

    def test_foreign_concat_is_rejected(self):
        with patch.object(pd, 'concat', lambda *args, **kwargs: None):
            with self.assertRaises(ValueError): guards.reference_guard(reference)

    def test_same_code_with_changed_numeric_defaults_is_rejected(self):
        original = pd.to_numeric
        replacement = FunctionType(original.__code__, original.__globals__,
            original.__name__, ('coerce', None), original.__closure__)
        replacement.__kwdefaults__ = original.__kwdefaults__
        with patch.object(pd, 'to_numeric', replacement):
            with self.assertRaises(ValueError): guards.reference_guard(reference)

    def test_foreign_le_operator_is_rejected(self):
        with patch.object(pd.Series, '__le__', lambda *args, **kwargs: None):
            with self.assertRaises(ValueError): guards.reference_guard(reference)

    def test_foreign_strip_is_rejected(self):
        from pandas.core.strings.accessor import StringMethods
        with patch.object(StringMethods, 'strip', lambda *args, **kwargs: None):
            with self.assertRaises(ValueError): guards.reference_guard(reference)

    def test_foreign_logger_warning_is_rejected(self):
        with patch.object(logging.Logger, 'warning', lambda *args, **kwargs: None):
            with self.assertRaises(ValueError): guards.reference_guard(reference)

    def test_actual_custom_handler_and_formatter_are_retained(self):
        seen = []
        class Formatter(logging.Formatter):
            def format(self, record):
                seen.append(('format', record.getMessage()))
                return 'actual:' + record.getMessage()
        class Handler(logging.Handler):
            def emit(self, record): seen.append(('emit', self.format(record)))
        handler = Handler(); handler.setFormatter(Formatter())
        logger = reference.logger
        with patch.object(logger, 'handlers', [handler]), patch.object(logger, 'propagate', False), patch.object(logger, 'disabled', False), patch.object(logger, 'level', logging.WARNING):
            guards.reference_guard(reference)
            logger.warning('frontier')
            guards.reference_guard(reference)
        self.assertEqual(seen, [('format', 'frontier'), ('emit', 'actual:frontier')])

    def test_callback_binding_mutation_forces_only_remaining_operations_mature(self):
        original = pd.to_numeric
        seen = []
        class Handler(logging.Handler):
            def emit(self, record):
                seen.append(record.getMessage())
                pd.to_numeric = lambda *args, **kwargs: None
        session = bridge.NativeSession('unused', '0' * 64)
        runtime = bridge._Runtime(session, reference)
        logger = reference.logger
        try:
            with patch.object(logger, 'handlers', [Handler()]), patch.object(logger, 'propagate', False), patch.object(logger, 'disabled', False), patch.object(logger, 'level', logging.WARNING):
                guards.reference_guard(reference)
                self.assertIsNone(runtime.after_log(logger.warning('mutate')))
            self.assertTrue(runtime.force_mature)
            self.assertTrue(runtime.binding_failures)
            with patch.object(session, 'request') as request:
                self.assertIs(runtime.calculate(32, 'PP3_BP4_criteria', 'eq', [True], lambda value: value), True)
                request.assert_not_called()
            self.assertEqual(seen, ['mutate'])
        finally: pd.to_numeric = original

    def test_callback_failure_prevents_after_log_and_any_remaining_operation(self):
        failure = RuntimeError('actual callback failure')
        class Handler(logging.Handler):
            def emit(self, record): raise failure
        runtime = bridge._Runtime(bridge.NativeSession('unused', '0' * 64), reference)
        logger = reference.logger
        with patch.object(logger, 'handlers', [Handler()]), patch.object(logger, 'propagate', False), patch.object(logger, 'disabled', False), patch.object(logger, 'level', logging.WARNING):
            guards.reference_guard(reference)
            with patch.object(runtime, 'after_log') as after:
                with self.assertRaises(RuntimeError) as caught:
                    runtime.after_log(logger.warning('fail'))
                self.assertIs(caught.exception, failure)
                after.assert_not_called()


if __name__ == '__main__': unittest.main()
