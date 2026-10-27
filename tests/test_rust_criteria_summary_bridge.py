"""Criteria summary atomic plans and exact historical fallback/publication."""
import ast
import copy
import csv
import gzip
import hashlib
import importlib.util
import logging
import os
from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
from types import SimpleNamespace
from typing import Dict, Tuple
import unittest
from unittest.mock import patch
import warnings

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'rust/tests'))
import acmg_scoring as scoring
import priva_rust_criteria_summary as bridge
from criteria_summary_source import definition, historical_summary_source, SUMMARY_SHA256
import criteria_summary_bridge_compare as summary_harness
import terminal_order_bridge_compare as ranking_harness


def original():
    text = (ROOT / 'scripts/acmg_scoring.py').read_text()
    names = ('fix_clingen_term', 'translate_strength_level', 'create_criteria_summary')
    nodes = [n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef) and n.name in names]
    nodes.append(ast.parse(historical_summary_source(text)).body[0])
    namespace = dict(pd=pd, np=np, Dict=Dict, Tuple=Tuple, os=os, gzip=gzip, pickle=pickle,
                     logger=scoring.logger)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), 'historical-summary', 'exec'), namespace)
    return namespace['summarize_acmg_criteria']


def frame():
    value = pd.DataFrame({'variant_id': ['a', 'b'], 'untouched': ['β\n"tab\t', 'second'],
                          'ACMG_quant_score': [.94, -.0], 'ACMG_class': ['old', 'old']},
                         index=pd.Index(['second', 'first'], name='caller rows'))
    value.columns.name = 'caller columns'
    value.attrs = {'nested': {'keep': ['original']}}
    value.flags.allows_duplicate_labels = False
    return value


def strengths(dtype='int64'):
    values = {name: pd.array([0, 0], dtype=dtype) for name in bridge.CRITERIA}
    values['PM2'] = pd.array([1, 1] if dtype in ('bool', 'boolean') else [2, 1], dtype=dtype)
    values['PP3'] = pd.array([0, 1], dtype=dtype)
    return values


def matrix(data, values):
    result = pd.DataFrame({name: values[name] for name in bridge.CRITERIA}, index=data.index)
    result['variant_id'] = data['variant_id']
    return result


def response(command, *, curated=False, mutate=None, **_kwargs):
    if curated:
        summaries, override = ['PVS1_Very Strong;BA1_Stand Alone', 'PM2_Supporting;PP3'], ['1', '0']
    else:
        summaries, override = ['PM2', 'PM2_Supporting;PP3'], ['0', '0']
    rows = []
    for i in range(2):
        row = dict(source_position=str(i), summary_kind='str', summary=summaries[i], overridden=override[i])
        row.update({name: '0' for name in bridge.CRITERIA})
        if curated and i == 0:
            row.update(PVS1='4', BA1='5')
        else:
            row.update(PM2=str(2 - i), PP3=str(i))
        rows.append(row)
    fields = ['source_position', 'summary_kind', 'summary', 'overridden', *bridge.CRITERIA]
    if mutate:
        mutate(fields, rows)
    with Path(command[command.index('--output') + 1]).open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    return SimpleNamespace(returncode=0)


def same_frame(test, left, right):
    pd.testing.assert_frame_equal(left, right, check_exact=True)
    test.assertEqual(left.attrs, right.attrs)
    test.assertEqual(left.flags.allows_duplicate_labels, right.flags.allows_duplicate_labels)


_PUBLICATION_EVENTS = []
_PUBLICATION_RHS_FAILURE = None


def publication_rhs(name):
    _PUBLICATION_EVENTS.append(name)
    if _PUBLICATION_RHS_FAILURE == name:
        raise RuntimeError('fixture RHS failed: ' + name)
    return pd.Series(['summary'] if name == 'summary' else [3], index=pd.RangeIndex(1))


def publication_single(df, criteria_dict, map_path, enabled=True):
    criteria_matrix = pd.DataFrame({'PVS1': [0], 'PS1': [0]}, index=df.index)
    overridden = pd.Series([True], index=df.index)
    df['ACMG_criteria'] = publication_rhs('summary')
    for col in ('PVS1', 'PS1'):
        criteria_matrix.loc[overridden, col] = publication_rhs(col)
    return df, criteria_matrix


def publication_multiline(df, criteria_dict, map_path, enabled=True):
    criteria_matrix = pd.DataFrame({'PVS1': [0], 'PS1': [0]}, index=df.index)
    overridden = pd.Series([True], index=df.index)
    df['ACMG_criteria'] = publication_rhs('summary')
    for col in ('PVS1', 'PS1'):
        criteria_matrix.loc[overridden, col] = publication_rhs(
            col
        )
    return df, criteria_matrix


class CriteriaSummaryBridge(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)

    def evidence(self, mapping, compressed=False):
        path = self.folder / ('map.pkl.gz' if compressed else 'map.pkl')
        opener = gzip.open if compressed else open
        with opener(path, 'wb') as stream:
            pickle.dump(mapping, stream)
        return str(path)

    def compare(self, values=None, evidence=None, *, enabled=True, data=None, fake=None):
        values = strengths() if values is None else values
        data = frame() if data is None else data
        path = self.evidence({} if evidence is None else evidence)
        left, right = data.copy(deep=True), data.copy(deep=True)
        expected = original()(left, values, path, enabled)
        with patch.object(bridge.subprocess, 'run', side_effect=response if fake is None else fake):
            actual = scoring.summarize_acmg_criteria(right, values, path, enabled, summary_binary='native')
        self.assertIs(actual[0], right)
        same_frame(self, expected[0], actual[0]); same_frame(self, expected[1], actual[1])
        return actual

    def test_exact_historical_bytes_restore_and_keyword_default(self):
        source = (ROOT / 'scripts/acmg_scoring.py').read_text()
        segment = historical_summary_source(source)
        self.assertEqual(hashlib.sha256(segment.encode()).hexdigest(), SUMMARY_SHA256)
        self.assertEqual(historical_summary_source(segment), segment)

    def test_definition_guard_rejects_missing_duplicate_async_decorated_and_nested_only(self):
        literal = 'def summarize_acmg_criteria():\n    pass\n'
        for source in ['', literal * 2, 'async ' + literal, '@wrap\n' + literal,
                       'def outside():\n    ' + literal.replace('\n', '\n    ')]:
            with self.subTest(source=source), self.assertRaises(AssertionError):
                definition(ast.parse(source))

    def test_guard_rejects_skipped_reordered_or_changed_native_publication(self):
        source = (ROOT / 'scripts/acmg_scoring.py').read_text()
        changes = [
            ("df['ACMG_criteria'] = _summary_plan.summary", 'pass'),
            ('for col in criteria_order:\n                    criteria_matrix.loc',
             'for col in reversed(criteria_order):\n                    criteria_matrix.loc'),
            ("_summary_plan.updates[col]", "_summary_plan.updates['PVS1']"),
            ('if _summary_plan is not None:', 'if True:'),
            ('return df, criteria_matrix\n    # Add', 'return criteria_matrix, df\n    # Add'),
            ('summary_binary=None', 'summary_binary=True'),
            ('if summary_binary is not None:', 'if summary_binary:'),
            ('_summary_plan = try_criteria_summary', '_summary_plan = skip_summary'),
        ]
        for old, new in changes:
            self.assertIn(old, source)
            with self.subTest(old=old), self.assertRaises(AssertionError):
                historical_summary_source(source.replace(old, new, 1))

    def test_guard_rejects_historical_change_even_when_dispatch_is_valid(self):
        source = (ROOT / 'scripts/acmg_scoring.py').read_text()
        for old, new in [("criteria_dict.get(name, np.zeros", "criteria_dict.get(name, np.ones"),
                         ("if not criteria_matrix[col].isin([0, 1, 2, 3, 4, 5]).all():", 'if False:'),
                         ('clingen_map_pkl.endswith(".gz")', 'clingen_map_pkl.endswith(".zip")')]:
            with self.subTest(old=old), self.assertRaises(AssertionError):
                historical_summary_source(source.replace(old, new, 1))

    def test_native_plan_is_atomic_and_keeps_nonmonotonic_axes_metadata(self):
        data = frame(); current = matrix(data, strengths())
        before, old_matrix = data.copy(deep=True), current.copy(deep=True)
        with patch.object(bridge.subprocess, 'run', side_effect=response) as run:
            result = bridge.try_criteria_summary(data, current, {}, 'native')
        self.assertIsNotNone(result); self.assertEqual(run.call_count, 1)
        self.assertEqual(result.summary.index.tolist(), ['second', 'first'])
        same_frame(self, data, before); same_frame(self, current, old_matrix)
        self.compare()

    def test_curated_summary_then_original_ordered_integer_column_updates(self):
        self.compare(evidence={'a': 'PVS1_Very, BA1_Stand'},
                     fake=lambda cmd, **kwargs: response(cmd, curated=True))

    def test_real_nonstring_nan_evidence_is_encoded_as_silence(self):
        data = frame(); current = matrix(data, strengths())
        captured = []
        def check(command, **kwargs):
            with Path(command[command.index('--input') + 1]).open(newline='') as stream:
                captured.extend(csv.DictReader(stream))
            return response(command)
        with patch.object(bridge.subprocess, 'run', side_effect=check):
            plan = bridge.try_criteria_summary(data, current, {'a': np.nan, 'b': 17, 'unused': None}, 'native')
        self.assertIsNotNone(plan)
        self.assertEqual([(r['evidence_kind'], r['evidence']) for r in captured], [('none', '')] * 2)
        self.compare(evidence={'a': np.nan, 'b': 17})

    def test_float_nullable_and_object_strength_publication_matches_source(self):
        for dtype in ['float64', 'float32', 'Int64', 'Float64', 'object']:
            with self.subTest(dtype=dtype):
                self.compare(values=strengths(dtype), evidence={'a': 'PVS1_Very, BA1_Stand'},
                             fake=lambda cmd, **kwargs: response(cmd, curated=True))

    def test_boolean_strengths_defer_curated_assignment_warnings_and_failures(self):
        for dtype in ['bool', 'boolean']:
            values = strengths(dtype)
            data = frame(); current = matrix(data, values)
            with patch.object(bridge.subprocess, 'run') as run:
                self.assertIsNone(bridge.try_criteria_summary(data, current, {'a': 'PS1'}, 'native'))
            run.assert_not_called()
            path = self.evidence({'a': 'PS1'})
            outcomes = []
            for function, kw in [(original(), {}), (scoring.summarize_acmg_criteria, {'summary_binary': 'native'})]:
                copy_df = data.copy(deep=True)
                with warnings.catch_warnings(record=True) as seen:
                    warnings.simplefilter('always')
                    try:
                        result = function(copy_df, values, path, **kw)
                        outcome = ('ok', result[1])
                    except Exception as error:
                        outcome = (type(error).__name__, str(error))
                outcomes.append((outcome, copy_df, [(type(w.message).__name__, str(w.message)) for w in seen]))
            self.assertEqual(outcomes[0][0][0], outcomes[1][0][0])
            if outcomes[0][0][0] == 'ok':
                same_frame(self, outcomes[0][0][1], outcomes[1][0][1])
            else:
                self.assertEqual(outcomes[0][0][1], outcomes[1][0][1])
            same_frame(self, outcomes[0][1], outcomes[1][1])
            self.assertEqual(outcomes[0][2], outcomes[1][2])

    def test_duplicate_index_columns_and_nonstring_ids_fall_back_without_launch(self):
        cases = []
        data = frame(); data.flags.allows_duplicate_labels = True; data.index = ['same', 'same']; cases.append(data)
        data = frame(); data.flags.allows_duplicate_labels = True
        data = pd.concat([data, data[['untouched']]], axis=1); cases.append(data)
        data = frame(); data['variant_id'] = [1, 2]; cases.append(data)
        for data in cases:
            current = matrix(data, strengths())
            with self.subTest(columns=list(data.columns)), patch.object(bridge.subprocess, 'run') as run:
                self.assertIsNone(bridge.try_criteria_summary(data, current, {}, 'native'))
                run.assert_not_called()

    def test_custom_map_keys_or_strength_objects_are_not_coerced(self):
        class Untouched:
            def __int__(self):
                raise AssertionError('custom int invoked')
            def __str__(self):
                raise AssertionError('custom string invoked')
        data = frame(); current = matrix(data, strengths())
        cases = [(current, {17: 'PS1'}), (current, {'a': Untouched()})]
        modified = current.astype({'PM2': object}); modified.iloc[0, modified.columns.get_loc('PM2')] = Untouched()
        cases.append((modified, {}))
        for index, (current, mapping) in enumerate(cases):
            with self.subTest(case=index), patch.object(bridge.subprocess, 'run', side_effect=response) as run:
                result = bridge.try_criteria_summary(data, current, mapping, 'native')
            self.assertIsNone(result); run.assert_not_called()

    def test_custom_map_callbacks_warnings_and_error_time_frame_match_original(self):
        path = self.evidence({})
        for fail in (False, True):
            outcomes = []
            for function, kwargs in [(original(), {}), (scoring.summarize_acmg_criteria,
                                                       {'summary_binary': 'native'})]:
                events = []
                class NonString:
                    @property
                    def __class__(self):
                        events.append('__class__')
                        if len(events) == 2:
                            warnings.warn('custom map callback', UserWarning)
                            if fail:
                                raise RuntimeError('custom map callback failed')
                        return object
                data = frame(); returned = None; error = None
                mapping = {'a': NonString()}
                with patch.object(scoring.pickle, 'load', new=lambda stream: mapping), \
                        patch.object(bridge.subprocess, 'run') as run, warnings.catch_warnings(record=True) as seen:
                    warnings.simplefilter('always')
                    try:
                        returned = function(data, strengths(), path, **kwargs)
                    except Exception as failure:
                        error = (type(failure).__name__, str(failure))
                    run.assert_not_called()
                outcomes.append((list(events), error, data,
                    [(type(w.message).__name__, str(w.message)) for w in seen], returned))
            self.assertEqual(outcomes[0][0:2], outcomes[1][0:2])
            self.assertEqual(len(outcomes[0][0]), 2)
            self.assertEqual(outcomes[0][3], outcomes[1][3])
            same_frame(self, outcomes[0][2], outcomes[1][2])
            if not fail:
                same_frame(self, outcomes[0][4][1], outcomes[1][4][1])

    def test_custom_extension_storage_is_rejected_before_dtype_callbacks(self):
        events = []
        armed = False
        class CustomIntegerArray(pd.arrays.IntegerArray):
            @property
            def dtype(self):
                if armed:
                    events.append('dtype')
                    raise RuntimeError('custom dtype getter')
                return super().dtype
        data = frame(); values = strengths()
        values['PVS1'] = CustomIntegerArray(np.array([1, 0], dtype='int64'), np.zeros(2, dtype=bool))
        current = matrix(data, values)
        armed = True
        with patch.object(bridge.subprocess, 'run') as run:
            self.assertIsNone(bridge.try_criteria_summary(data, current, {}, 'native'))
        self.assertEqual(events, []); run.assert_not_called()
        class CustomNdarray(np.ndarray):
            def __getattribute__(self, name):
                if name == 'dtype':
                    events.append('raw dtype')
                    raise RuntimeError('custom raw dtype getter')
                return super().__getattribute__(name)
        array = pd.array([1, 0], dtype='Int64')
        array._data = np.array([1, 0], dtype='int64').view(CustomNdarray)
        with self.assertRaises(ValueError):
            bridge._ordinary_array(array)
        self.assertEqual(events, [])

    def test_extension_dtype_callback_warnings_and_failure_match_original(self):
        for fail in (False, True):
            outcomes = []
            for function, kwargs in [(original(), {}), (scoring.summarize_acmg_criteria,
                                                       {'summary_binary': 'native'})]:
                events = []
                armed = False
                class CustomIntegerArray(pd.arrays.IntegerArray):
                    @property
                    def dtype(self):
                        if armed:
                            events.append(len(events) + 1)
                            if len(events) == 13:
                                warnings.warn('extension dtype callback', UserWarning)
                                if fail:
                                    raise RuntimeError('extension dtype callback failed')
                        return super().dtype
                data = pd.DataFrame({'variant_id': ['a']})
                values = {'PVS1': CustomIntegerArray(np.array([1], dtype='int64'), np.zeros(1, dtype=bool))}
                error = None; returned = None; armed = True
                with patch.object(bridge.subprocess, 'run') as run, warnings.catch_warnings(record=True) as seen:
                    warnings.simplefilter('always')
                    try:
                        returned = function(data, values, 'UNREAD', False, **kwargs)
                    except Exception as failure:
                        error = (type(failure).__name__, str(failure))
                    run.assert_not_called()
                armed = False
                outcomes.append((list(events), error, data,
                    [(type(w.message).__name__, str(w.message)) for w in seen], returned))
            self.assertEqual(outcomes[0][0:2], outcomes[1][0:2])
            self.assertEqual(outcomes[0][0], list(range(1, 14)))
            self.assertEqual(outcomes[0][3], [('UserWarning', 'extension dtype callback')])
            self.assertEqual(outcomes[0][3], outcomes[1][3])
            same_frame(self, outcomes[0][2], outcomes[1][2])
            if not fail:
                same_frame(self, outcomes[0][4][1], outcomes[1][4][1])

    def test_empty_frame_and_unqualified_versions_do_not_launch(self):
        data = frame(); current = matrix(data, strengths())
        with patch.object(bridge.subprocess, 'run') as run:
            self.assertIsNone(bridge.try_criteria_summary(data.iloc[:0], current.iloc[:0], {}, 'native'))
            with patch.object(pd, '__version__', 'unsupported'):
                self.assertIsNone(bridge.try_criteria_summary(data, current, {}, 'native'))
        run.assert_not_called()

    def test_runtime_guard_rejects_changed_scalar_helpers_and_publication_code(self):
        data = frame(); current = matrix(data, strengths())
        for name in ['create_criteria_summary', 'fix_clingen_term', 'translate_strength_level',
                     'summarize_acmg_criteria']:
            with self.subTest(name=name), patch.object(scoring, name, lambda *args, **kwargs: 'changed'):
                with patch.object(bridge.subprocess, 'run') as run:
                    self.assertIsNone(bridge.try_criteria_summary(data, current, {}, 'native'))
                    run.assert_not_called()
        function = scoring.fix_clingen_term
        previous_code = function.__code__
        try:
            function.__code__ = (lambda value: 'changed').__code__
            with patch.object(bridge.subprocess, 'run') as run:
                self.assertIsNone(bridge.try_criteria_summary(data, current, {}, 'native'))
                run.assert_not_called()
        finally:
            function.__code__ = previous_code

    def test_changed_helper_uses_current_python_semantics_once_without_native_publication(self):
        path = self.evidence({'a': 'PS1'})
        for name, changed in [('create_criteria_summary', lambda *args, **kwargs: 'PM2'),
                              ('fix_clingen_term', lambda value: 'PM2'),
                              ('translate_strength_level', lambda *args: 5)]:
            with self.subTest(name=name), patch.object(scoring, name, changed):
                left, right = frame(), frame()
                expected = scoring.summarize_acmg_criteria(left, strengths(), path)
                with patch.object(bridge.subprocess, 'run') as run:
                    actual = scoring.summarize_acmg_criteria(right, strengths(), path, summary_binary='native')
                    run.assert_not_called()
                same_frame(self, expected[0], actual[0]); same_frame(self, expected[1], actual[1])

    def test_native_nonzero_timeout_and_missing_output_leave_inputs_unchanged(self):
        data = frame(); current = matrix(data, strengths()); before = data.copy(deep=True)
        for failure in [SimpleNamespace(returncode=2), subprocess.TimeoutExpired('native', 1),
                        SimpleNamespace(returncode=0), OSError('missing executable')]:
            kwargs = {'side_effect': failure} if isinstance(failure, Exception) else {'return_value': failure}
            with self.subTest(failure=type(failure).__name__), patch.object(bridge.subprocess, 'run', **kwargs):
                self.assertIsNone(bridge.try_criteria_summary(data, current, {}, 'native'))
            same_frame(self, before, data)

    def test_transport_rejects_changed_order_rows_schema_strengths_and_nulls(self):
        mutations = [lambda f, r: r.reverse(), lambda f, r: r.pop(),
            lambda f, r: r.append(dict(r[-1])), lambda f, r: f.reverse(),
            lambda f, r: r[0].update(overridden='true'),
            lambda f, r: r[0].update(overridden='1'),
            lambda f, r: r[0].update(summary_kind='nan'),
            lambda f, r: r[0].update(summary=''),
            lambda f, r: r[0].update(PM2='6'), lambda f, r: r[0].update(PM2='1')]
        data = frame(); current = matrix(data, strengths())
        for mutate in mutations:
            with self.subTest(mutate=mutate), patch.object(bridge.subprocess, 'run',
                    side_effect=lambda cmd, **kwargs: response(cmd, mutate=mutate)):
                self.assertIsNone(bridge.try_criteria_summary(data, current, {}, 'native'))

    def test_default_and_override_false_do_not_call_adapter_or_load_map(self):
        data = frame()
        with patch.object(bridge, 'try_criteria_summary', side_effect=AssertionError('unexpected adapter')):
            scoring.summarize_acmg_criteria(data, strengths(), 'unread-missing-map', False)
        data = frame()
        with patch.object(pickle, 'load', side_effect=AssertionError('unexpected map read')):
            with patch.object(bridge.subprocess, 'run', side_effect=response):
                scoring.summarize_acmg_criteria(data, strengths(), 'unread-missing-map', False,
                                                summary_binary='native')

    def test_matrix_validation_and_missing_map_fail_before_native_attempt(self):
        bad = strengths(); bad['PVS1'] = [6, 0]
        for values, category in [(bad, ValueError), (strengths(), FileNotFoundError)]:
            with self.subTest(category=category), patch.object(bridge, 'try_criteria_summary') as attempt:
                with self.assertRaises(category):
                    scoring.summarize_acmg_criteria(frame(), values, 'missing-map', summary_binary='native')
                attempt.assert_not_called()

    def test_malformed_clingen_fallback_preserves_summary_mutation_at_error(self):
        path = self.evidence({'a': 'XPS1', 'b': 'PVS1_Unknown'})
        results = []
        for function, kwargs in [(original(), {}), (scoring.summarize_acmg_criteria, {'summary_binary': 'native'})]:
            data = frame()
            with patch.object(bridge.subprocess, 'run', return_value=SimpleNamespace(returncode=2)):
                with self.assertRaises(KeyError) as caught:
                    function(data, strengths(), path, **kwargs)
            results.append((data, str(caught.exception)))
        same_frame(self, results[0][0], results[1][0]); self.assertEqual(results[0][1], results[1][1])
        self.assertEqual(results[1][0]['ACMG_criteria'].tolist(), ['XPS1', 'PVS1_Unknown'])



class SummaryProofAndPublication(unittest.TestCase):
    def test_proof_view_restores_exact_summary_and_preserves_every_other_byte(self):
        source = (ROOT / 'scripts/acmg_scoring.py').read_text()
        current = ast.get_source_segment(source, definition(ast.parse(source)))
        expected = source.replace(current, historical_summary_source(source), 1)
        self.assertEqual(ranking_harness.summary_proof_source(source), expected)
        changed = source.replace('prior_probability=0.1', 'prior_probability=0.2', 1)
        self.assertNotEqual(changed, source)
        changed_view = ranking_harness.summary_proof_source(changed)
        self.assertIn('prior_probability=0.2', changed_view)
        self.assertNotEqual(ast.dump(ast.parse(expected)), ast.dump(ast.parse(changed_view)))

    def test_proof_view_rejects_relocated_nested_duplicate_and_missing_hook(self):
        source = (ROOT / 'scripts/acmg_scoring.py').read_text()
        function = definition(ast.parse(source))
        hook = next(n for n in function.body if isinstance(n, ast.If)
                    and ast.unparse(n.test) == 'summary_binary is not None')
        lines = source.splitlines(keepends=True)
        block = lines[hook.lineno - 1:hook.end_lineno]
        del lines[hook.lineno - 1:hook.end_lineno]
        removed = ''.join(lines)
        original = definition(ast.parse(removed))
        publication = next(n for n in original.body if isinstance(n, ast.Assign)
                           and ast.unparse(n.targets[0]) == "df['ACMG_criteria']")
        loading = next(n for n in original.body if isinstance(n, ast.If)
                       and ast.unparse(n.test) == 'apply_clingen_override')
        changes = [('missing', removed)]
        for name, position, indentation in [('after_publication', publication.end_lineno, ''),
                ('before_loading', loading.lineno - 1, ''),
                ('nested', loading.body[0].lineno - 1, '    ')]:
            changed = removed.splitlines(keepends=True)
            changed[position:position] = [indentation + line for line in block]
            changes.append((name, ''.join(changed)))
        duplicate = source.splitlines(keepends=True)
        duplicate[hook.lineno - 1:hook.lineno - 1] = block
        changes.append(('duplicate', ''.join(duplicate)))
        for name, changed in changes:
            ast.parse(changed)
            with self.subTest(name=name), self.assertRaises(AssertionError):
                ranking_harness.summary_proof_source(changed)

    def test_opcode_observer_records_only_completed_stores_and_keeps_partial_errors(self):
        global _PUBLICATION_RHS_FAILURE
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            map_path = folder / 'unused-map.pkl'
            with map_path.open('wb') as stream: pickle.dump({}, stream)
            case = dict(name='publication-control', frame=pd.DataFrame({'variant_id': ['a'],
                        'ACMG_criteria': ['untouched']}), criteria={'PVS1': [0], 'PS1': [0]},
                        map_path=str(map_path), enabled=True, phases=True)
            dummy_bridge = SimpleNamespace(try_criteria_summary=lambda *a, **k: None)
            configurations = [
                ('single', publication_single, None, None, ['summary', 'column:PVS1', 'column:PS1'], [3, 3], 'summary'),
                ('multiline', publication_multiline, None, None, ['summary', 'column:PVS1', 'column:PS1'], [3, 3], 'summary'),
                ('RHS_summary', publication_multiline, 'summary', None, [], [0, 0], 'untouched'),
                ('RHS_column', publication_multiline, 'PS1', None, ['summary', 'column:PVS1'], [3, 0], 'summary'),
                ('STORE_before', publication_multiline, None, 'before_column', ['summary', 'column:PVS1'], [3, 0], 'summary'),
                ('STORE_after', publication_multiline, None, 'after_column', ['summary', 'column:PVS1'], [3, 3], 'summary'),
                ('summary_STORE_after', publication_multiline, None, 'after_summary', [], [0, 0], 'summary'),
            ]
            real_loc = pd.core.indexing._LocIndexer.__setitem__
            real_df = pd.DataFrame.__setitem__
            for tag, function, rhs_failure, store_failure, expected, cells, summary in configurations:
                _PUBLICATION_EVENTS.clear(); _PUBLICATION_RHS_FAILURE = rhs_failure
                def loc(self_, key, value):
                    if store_failure in ('before_column', 'after_column') and key[1] == 'PS1':
                        if store_failure == 'after_column': real_loc(self_, key, value)
                        raise RuntimeError('fixture STORE failed')
                    return real_loc(self_, key, value)
                def setitem(self_, key, value):
                    if store_failure == 'after_summary' and key == 'ACMG_criteria':
                        real_df(self_, key, value)
                        raise RuntimeError('fixture STORE failed')
                    return real_df(self_, key, value)
                artifact = folder / (tag + '.pkl')
                code = function.__code__; previous_trace = sys.gettrace()
                with self.subTest(tag=tag), patch.object(pd.core.indexing._LocIndexer, '__setitem__', loc), \
                        patch.object(pd.DataFrame, '__setitem__', setitem):
                    result = summary_harness.execute(function, logging.getLogger('acmg_scoring'),
                        case, 'default', Path('/unused'), dummy_bridge, artifact)
                self.assertIs(function.__code__, code)
                self.assertIs(sys.gettrace(), previous_trace)
                self.assertEqual([x['phase'] for x in result['outcome']['checkpoints']], expected)
                saved = pickle.loads(artifact.read_bytes())
                self.assertEqual(saved['matrix'].iloc[0].tolist(), cells)
                self.assertEqual(saved['frame']['ACMG_criteria'].tolist(), [summary])
                if rhs_failure or store_failure:
                    self.assertEqual(result['outcome']['error']['type'], 'RuntimeError')
                else:
                    self.assertIsNone(result['outcome']['error'])
                    self.assertEqual(_PUBLICATION_EVENTS, ['summary', 'PVS1', 'PS1'])
            _PUBLICATION_RHS_FAILURE = None


class SummaryMissingAPIs(unittest.TestCase):
    def test_fresh_import_reaches_actual_admission_without_getters_native_or_mutation(self):
        active, callbacks = [False], []
        class Unknown(pd.arrays.IntegerArray):
            @property
            def dtype(self):
                if active[0]:
                    callbacks.append('dtype')
                    raise RuntimeError('unsupported dtype callback')
                return super().dtype
        data = frame()
        ordinary = matrix(data, strengths())
        values = strengths()
        values['PVS1'] = Unknown(np.array([1, 1], dtype=np.int64), np.zeros(2, dtype=bool))
        unknown = matrix(data, values)
        before = summary_harness.fingerprint(data)
        cases = [('ordinary', ordinary), ('unknown', unknown)]
        states = {name: summary_harness.fingerprint(value) for name, value in cases}
        functions = tuple(getattr(scoring, name) for name in bridge._FUNCTION_PINS)
        for version in ('unsupported-fixture', '2.2.3'):
            with self.subTest(version=version), patch.object(pd, '__version__', version), \
                    patch.object(pd.arrays, 'NumpyExtensionArray'):
                del pd.arrays.NumpyExtensionArray
                spec = importlib.util.spec_from_file_location('fresh_summary_missing_API', ROOT / 'scripts/priva_rust_criteria_summary.py')
                fresh = importlib.util.module_from_spec(spec); spec.loader.exec_module(fresh)
                self.assertEqual(fresh._EXTENSION_DTYPES, ())
                self.assertEqual(fresh._ARRAY_TYPES, ())
                for name, current in cases:
                    guard, admission, storage = fresh._source_guard, fresh._admit, fresh._ordinary_matrix_storage
                    counts = {'guard_completed': 0, 'admission': 0, 'storage': 0}
                    failures = []
                    def actual_guard(*args, **kwargs):
                        result = guard(*args, **kwargs); counts['guard_completed'] += 1; return result
                    def actual_admission(*args, **kwargs):
                        counts['admission'] += 1
                        try: return admission(*args, **kwargs)
                        except Exception as exc:
                            failures.append((type(exc).__name__, str(exc))); raise
                    def actual_storage(*args, **kwargs):
                        counts['storage'] += 1; return storage(*args, **kwargs)
                    callbacks.clear(); active[0] = True
                    try:
                        with patch.object(fresh, '_source_guard', side_effect=actual_guard), \
                                patch.object(fresh, '_admit', side_effect=actual_admission), \
                                patch.object(fresh, '_ordinary_matrix_storage', side_effect=actual_storage), \
                                patch.object(fresh.subprocess, 'run') as native:
                            result = fresh.try_criteria_summary(data, current, {}, 'unused', source_functions=functions)
                            native.assert_not_called()
                    finally: active[0] = False
                    self.assertIsNone(result)
                    self.assertEqual(callbacks, [])
                    self.assertEqual(counts, {'guard_completed': 1, 'admission': 1,
                                             'storage': 0 if version == 'unsupported-fixture' else 1})
                    self.assertEqual(len(failures), 1)
                    if version == 'unsupported-fixture':
                        self.assertEqual(failures[0], ('ValueError', 'unqualified pandas/NumPy context'))
                    elif name == 'ordinary':
                        self.assertEqual(failures[0], ('ValueError', 'unsupported criterion representation requires Python'))
                    else:
                        self.assertEqual(failures[0][0], 'AttributeError')
                        self.assertIn('NumpyExtensionArray', failures[0][1])
                    self.assertEqual(summary_harness.fingerprint(data), before)
                    self.assertEqual(summary_harness.fingerprint(current), states[name])

if __name__ == '__main__':
    unittest.main()
