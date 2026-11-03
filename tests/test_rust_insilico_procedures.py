"""New Rust-control host profile: source/effect/transport controls only.

These controls execute no original criterion, nested score-helper body, ELF or
model. Actual Rust procedure parity is a separate saved-oracle qualification.
"""
from collections import Counter
import logging
import operator
from pathlib import Path
import struct
import sys
from types import CodeType, SimpleNamespace
import unittest
from unittest.mock import patch
import warnings

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import numpy as np
import pandas as pd
import acmg_pp3_bp4_bp7_insilico as reference
import priva_rust_insilico as columns
import priva_rust_insilico_procedures as procedures


def host():
    raw = Path(reference.__file__).read_bytes()
    top = compile(raw, str(reference.__file__), 'exec', dont_inherit=True)
    codes = {code.co_name: code for code in top.co_consts if type(code) is CodeType}
    codes['_coerce_score_column'] = next(code for code in codes['PP3_BP4_criteria'].co_consts
        if type(code) is CodeType and code.co_name == '_coerce_score_column')
    session = procedures.ProcedureSession('unused', '0' * 64)
    return procedures._Host(reference, raw, codes, session)


def request(host, action, details, function='BP7_criteria', line=134):
    return host.service([('host', str(host.sequence), action, function, str(line), procedures.PROFILE), *details])


def call_rows(operation, target, args, function='', receiver='', keywords=()):
    return [('call', operation, target, str(function), str(receiver), ''),
        *(('argument', str(i), str(value), '', '', '') for i, value in enumerate(args)),
        *(('keyword', str(i), name, str(value), '', '') for i, (name, value) in enumerate(keywords))]


class ProcedureHostTests(unittest.TestCase):
    def test_frame_has_no_eager_admission_and_actual_getter_error_is_preserved(self):
        events = []
        failure = RuntimeError('actual lazy getter')
        class Frame:
            def __getitem__(self, key):
                events.append(key)
                if key == 'Conservation': raise failure
                return object()
        service = host(); frame = service.put(Frame())
        for key, line in (('Consequence', 133), ('splicing_len_changing', 134)):
            request(service, 'column', [('object', '0', str(frame), '', '', ''),
                ('text', '0', key, '', '', '')], line=line)
        with self.assertRaises(RuntimeError) as actual:
            request(service, 'column', [('object', '0', str(frame), '', '', ''),
                ('text', '0', 'Conservation', '', '', '')], line=135)
        self.assertIs(actual.exception, failure)
        self.assertEqual(events, ['Consequence', 'splicing_len_changing', 'Conservation'])

    def test_duplicate_columns_are_actual_getter_objects_not_whole_frame_rejection(self):
        service = host()
        frame = pd.DataFrame([['synonymous_variant', 'missense_variant']], columns=['Consequence', 'Consequence'])
        handle = service.put(frame)
        status, result = request(service, 'column', [('object', '0', str(handle), '', '', ''),
            ('text', '0', 'Consequence', '', '', '')], line=133)
        self.assertEqual(status, 'object')
        self.assertIs(type(service.get(result[0][2])), pd.DataFrame)
        self.assertEqual(frame._item_cache, {})

    def test_list_pvs1_comparison_uses_one_actual_operator_at_its_frontier(self):
        service = host(); left, right = service.put([0, 3]), service.put(3)
        details = call_rows('lt', 'lt', [left, right])
        status, payload = request(service, 'prepare', details, function='PP3_BP4_criteria', line=103)
        self.assertEqual((status, payload), ('mature', []))
        with self.assertRaises(TypeError) as actual:
            request(service, 'mature', details, function='PP3_BP4_criteria', line=103)
        self.assertEqual(str(actual.exception), "'<' not supported between instances of 'list' and 'int'")
        self.assertIsNone(service.pending)
        with self.assertRaises(procedures.TransportError):
            request(service, 'mature', details, function='PP3_BP4_criteria', line=103)

    def test_callable_confidence_is_mature_map_not_original_function_fallback(self):
        service = host(); series = service.put(pd.Series(['reviewed_by_expert_panel', None]))
        calls = []
        confidence = service.put(lambda value: calls.append(value) or 3)
        method = service.put(service.objects[series].map)
        service.attributes[method] = (series, 'map')
        ignore = service.put('ignore')
        details = call_rows('map', 'bound', [series, confidence], method, series, [('na_action', ignore)])
        status, _ = request(service, 'prepare', details, function='PP3_BP4_criteria', line=91)
        self.assertEqual(status, 'mature')
        status, result = request(service, 'mature', details, function='PP3_BP4_criteria', line=91)
        self.assertEqual(calls, ['reviewed_by_expert_panel'])
        self.assertEqual(service.get(result[0][2]).iloc[0], 3)
        self.assertEqual(service.session.metrics['original_entries'], 0)

    def test_float_boolean_fill_keeps_actual_pandas_warning_and_source_line(self):
        service = host(); value = service.put(pd.Series([float('nan')]))
        method = service.put(service.objects[value].fillna); service.attributes[method] = (value, 'fillna')
        false = service.put(False)
        details = call_rows('fillna', 'bound', [value, false], method, value)
        status, _ = request(service, 'prepare', details, function='PP3_BP4_criteria', line=74)
        self.assertEqual(status, 'mature')
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always')
            _status, result = request(service, 'mature', details, function='PP3_BP4_criteria', line=74)
        self.assertIs(service.get(result[0][2]).iloc[0], np.False_)
        self.assertEqual(len(captured), 1)
        self.assertEqual(captured[0].lineno, 74)
        self.assertEqual(captured[0].filename, reference.__file__)

    def test_handler_callback_executes_once_and_its_exception_is_not_replayed(self):
        service = host(); calls = []; failure = RuntimeError('handler failure')
        class Handler(logging.Handler):
            def emit(self, record): calls.append(record.getMessage()); raise failure
        handler = Handler(); logger = logging.Logger('procedure-source-callback', level=logging.INFO)
        logger.addHandler(handler); logger.propagate = False
        method = service.put(logger.info); message = service.put('source message')
        details = call_rows('log', 'bound', [message], method)
        self.assertEqual(request(service, 'prepare', details, function='PP3_BP4_criteria', line=104)[0], 'mature')
        with self.assertRaises(RuntimeError) as actual:
            request(service, 'mature', details, function='PP3_BP4_criteria', line=104)
        self.assertIs(actual.exception, failure); self.assertEqual(calls, ['source message'])
        self.assertIsNone(service.pending)

    def test_full_code_tree_rejects_filename_and_line_table_changes(self):
        def example(): return 1
        code = example.__code__
        filename = code.replace(co_filename='changed')
        self.assertEqual(filename, code)  # CPython equality is incomplete provenance.
        self.assertNotEqual(procedures._code_tree(filename), procedures._code_tree(code))
        changed_lines = code.replace(co_linetable=code.co_linetable + b'\x00')
        self.assertNotEqual(procedures._code_tree(changed_lines), procedures._code_tree(code))

    def test_literals_use_actual_original_constant_identity(self):
        service = host()
        original = next(value for value in reference.BP7_criteria.__code__.co_consts
                        if value == 'synonymous_variant')
        result = service._constant(('BP7_criteria', 133), 'synonymous_variant')
        self.assertIs(result, original)
        with self.assertRaises(procedures.TransportError):
            service._constant(('BP7_criteria', 133), 'unrepresented policy constant')

    def test_helper_definition_evaluates_annotations_without_calling_original_body(self):
        service = host(); marker = object()
        helper = service._effect(('PP3_BP4_criteria', 44), 'define_score_helper', [marker], [])
        self.assertIs(helper.__closure__[0].cell_contents, marker)
        self.assertEqual(helper.__annotations__['missing_default'], float | None)
        self.assertIs(helper.__annotations__['return'], pd.Series)
        self.assertIsNone(helper.__defaults__[0])
        self.assertEqual(service.session.metrics['original_entries'], 0)

    def test_missingness_methods_prepare_and_commit_without_mature_replay(self):
        for operation, expected in (('isna', [False, True]), ('notna', [True, False])):
            with self.subTest(operation=operation):
                service = host()
                series = pd.Series([1., np.nan], index=['left', 'right'], name='score')
                value = service.put(series)
                status, payload = request(service, 'attribute', [
                    ('object', '0', str(value), '', '', ''),
                    ('text', '0', operation, '', '', '')])
                self.assertEqual(status, 'object')
                method = payload[0][2]
                details = call_rows(operation, 'bound', [value], method, value)
                status, payload = request(service, 'prepare', details)
                self.assertEqual(status, 'prepared')
                self.assertEqual(payload, [('operands', '0',
                    columns._request(service.sequence, operation, [series]).decode(), '', '', '')])
                # Exercise the actual host protocol with a typed kernel-result
                # fixture; no ELF or original criterion is executed here.
                native = procedures._bytes([
                    ('result', str(service.sequence), 'series', 'bool', '2', columns.PROFILE),
                    *(('index', '0', str(i), 'label', columns._label(label), '')
                      for i, label in enumerate(('left', 'right'))),
                    *(('cell', '0', str(i), 'bool', str(int(result)), '')
                      for i, result in enumerate(expected))]).decode()
                with patch.object(service, '_mature') as mature:
                    status, payload = request(service, 'commit', [
                        *details, ('native', '0', native, '', '', '')])
                    mature.assert_not_called()
                result = service.get(payload[0][2])
                self.assertEqual(status, 'object')
                self.assertEqual(result.tolist(), expected)
                self.assertIs(result.index, series.index)
                self.assertEqual(result.name, series.name)
                self.assertEqual(service.session.metrics['native_operations'], Counter({operation: 1}))
                self.assertFalse(service.session.metrics['mature_operations'])

    def test_visited_operand_binding_change_refuses_native_without_evaluating_again(self):
        service = host(); value = service.put(pd.Series([1., 2.]))
        bound = service.put(service.objects[value].notna); service.attributes[bound] = (value, 'notna')
        details = call_rows('notna', 'bound', [value], bound, value)
        with patch.object(pd.Series, 'notna', side_effect=AssertionError('must use already bound method')):
            status, _ = request(service, 'prepare', details)
            self.assertEqual(status, 'mature')
            _status, result = request(service, 'mature', details)
        self.assertEqual(service.get(result[0][2]).tolist(), [True, True])

    def test_nullable_boolean_assignment_keeps_one_actual_ownership_conversion(self):
        for dtype in (bool, 'boolean'):
            with self.subTest(dtype=dtype):
                service = host()
                array = np.zeros(2, dtype=int)
                mask = pd.Series([True, False], dtype=dtype)
                handles = [service.put(value) for value in (array, mask, 1)]
                details = call_rows('assign', 'assignment', handles)
                self.assertEqual(request(service, 'prepare', details, line=137)[0], 'prepared')
                native = procedures._bytes([
                    ('result', str(service.sequence), 'array', 'int64', '2', columns.PROFILE),
                    ('cell', '0', '0', 'int', '1', ''),
                    ('cell', '0', '1', 'int', '0', '')]).decode()
                at = service._at; calls = []
                def forward(site, function, args, kwargs=None):
                    calls.append(function)
                    return at(site, function, args, kwargs)
                with patch.object(service, '_at', side_effect=forward), patch.object(service, '_mature') as mature:
                    status, payload = request(service, 'commit', [*details, ('native', '0', native, '', '', '')], line=137)
                    mature.assert_not_called()
                self.assertEqual(calls, [np.asarray, operator.setitem])
                self.assertEqual(status, 'object')
                self.assertIs(service.get(payload[0][2]), array)
                self.assertEqual(array.tolist(), [1, 0])
                self.assertEqual(service.session.metrics['native_operations'], Counter(assign=1))

    def test_missing_nullable_mask_cannot_publish_a_fabricated_native_assignment(self):
        service = host(); array = np.zeros(2, dtype=int)
        mask = pd.Series([True, pd.NA], dtype='boolean')
        details = call_rows('assign', 'assignment', [service.put(value) for value in (array, mask, 1)])
        self.assertEqual(request(service, 'prepare', details, line=137)[0], 'prepared')
        # The Rust kernel must request the actual mature assignment for a
        # missing mask. A fabricated completed result is transport corruption.
        native = procedures._bytes([
            ('result', str(service.sequence), 'array', 'int64', '2', columns.PROFILE),
            ('cell', '0', '0', 'int', '1', ''),
            ('cell', '0', '1', 'int', '0', '')]).decode()
        with patch.object(service, '_mature') as mature:
            with self.assertRaisesRegex(procedures.TransportError, 'missing Boolean mask'):
                request(service, 'commit', [*details, ('native', '0', native, '', '', '')], line=137)
            mature.assert_not_called()
        self.assertEqual(array.tolist(), [0, 0])
        self.assertFalse(service.session.metrics['native_operations'])
        self.assertIsNone(service.source_error)

    def test_explicit_release_drops_only_requested_source_temporaries_in_order(self):
        service = host(); events = []
        class Temporary:
            def __init__(self, label): self.label = label
            def __del__(self): events.append(self.label)
        left, right = service.put(Temporary('left')), service.put(Temporary('right'))
        status, payload = request(service, 'release', [('object', '0', str(left), '', '', ''),
            ('object', '1', str(right), '', '', '')])
        self.assertEqual((status, payload), ('released', []))
        self.assertEqual(events, ['left', 'right'])
        with self.assertRaises(procedures.TransportError):
            request(service, 'release', [('object', '0', str(left), '', '', '')])


class ProcedureTransportTests(unittest.TestCase):
    def test_quoted_protocol_and_typed_bits_are_lossless(self):
        rows = [('host-result', '0', 'object', procedures.PROFILE, '', ''),
                ('object', '0', '3', '', '', '')]
        self.assertEqual(procedures._rows(procedures._bytes(rows)), [list(row) for row in rows])
        self.assertEqual(procedures._rows(procedures._bytes([('x', 'λ\t\n', '', '', '', '')]))[0][1], 'λ\t\n')
        for value in ('true', '-1', '١', str(1 << 64)):
            with self.assertRaises(procedures.TransportError): procedures._number(value)

    def test_corrupt_unknown_sequence_and_foreign_handles_fail_closed(self):
        service = host()
        for rows in ([('host', '1', 'truth', 'BP7_criteria', '133', procedures.PROFILE)],
                     [('host', '0', 'truth', 'BP7_criteria', '133', procedures.PROFILE),
                      ('object', '0', '999', '', '', '')]):
            with self.assertRaises(procedures.TransportError): service.service(rows)
        for raw in (b'\xff', b'"unterminated', b'a\tb\n'):
            with self.assertRaises(procedures.TransportError): procedures._rows(raw)

    def test_malformed_native_result_is_transport_failure_never_mature_fallback(self):
        service = host(); left, right = service.put(pd.Series([1.])), service.put(2)
        details = call_rows('lt', 'lt', [left, right])
        self.assertEqual(request(service, 'prepare', details)[0], 'prepared')
        with patch.object(service, '_mature') as mature:
            with self.assertRaises(procedures.TransportError):
                request(service, 'commit', [*details, ('native', '0', 'corrupt', '', '', '')])
            mature.assert_not_called()

    def test_legacy_unsupported_commit_is_impossible_transport_status(self):
        service = host(); left, right = service.put(pd.Series([1.])), service.put(2)
        details = call_rows('lt', 'lt', [left, right])
        self.assertEqual(request(service, 'prepare', details)[0], 'prepared')
        raw = columns._request(0, 'lt', [1, 2])
        # This is a valid legacy error frame but not a possible procedure-v2
        # commit: the driver must ask for the single mature operation instead.
        frame = procedures._bytes([('unsupported', '0', 'legacy kernel rejection', '', '', columns.PROFILE)]).decode()
        with patch.object(service, '_mature') as mature:
            with self.assertRaisesRegex(procedures.TransportError, 'impossible unsupported native commit'):
                request(service, 'commit', [*details, ('native', '0', frame, '', '', '')])
            mature.assert_not_called()
        self.assertIsNone(service.source_error)

    def test_scope_has_no_whole_original_call_or_activation(self):
        import ast
        tree = ast.parse(Path(procedures.__file__).read_bytes())
        direct = [node for node in ast.walk(tree) if type(node) is ast.Call
                  and type(node.func) is ast.Attribute and node.func.attr in ('PP3_BP4_criteria', 'BP7_criteria')]
        self.assertEqual(direct, [])
        self.assertFalse(hasattr(procedures, 'install'))


if __name__ == '__main__': unittest.main()
