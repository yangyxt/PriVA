#!/usr/bin/env python3
"""Source-pinned whole-summary comparison, including error-time mutations.

prepare --candidate ROOT --reference ROOT --output NEW_DIR
worker --candidate ROOT --reference ROOT --inputs DIR --mode reference|default|native|fallback
       --profile fresh|source-bootstrap --binary PATH --output NEW_DIR
compare --inputs DIR --workers DIR --output NEW.json

Only assigned criteria summaries and ClinGen replacements execute here. No
posterior, criterion assignment, ranking, external scientific cache or model.
All paths are explicit; this reusable harness has no host-specific data gate.
"""
import argparse
import ast
import contextlib
import dis
import gzip
import hashlib
import importlib.util
import json
import logging
import os
from pathlib import Path
import pickle
import sys
import traceback
from typing import Dict, Tuple
import warnings
from unittest.mock import patch

import numpy as np
import pandas as pd

from criteria_summary_source import historical_summary_source
from ranking_control_oracle import encode


CRITERIA = ['PVS1', 'PS1', 'PS2', 'PS3', 'PS4', 'PM1', 'PM2', 'PM3', 'PM4',
            'PM5', 'PM6', 'PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'BA1', 'BS1',
            'BS2', 'BS3', 'BS4', 'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7']
REFERENCE_SHA = '9111a76e3fc94a064e14d97b535bb1dd2132f7c61836f8e8892477d4c5e69646'
BASE_CANDIDATE_SHA = '67e53cc9f8c38ef82ab170517396a633f16017af3501725c0687aa7ddc952ffb'
FUNCTIONS = {'fix_clingen_term', 'translate_strength_level', 'create_criteria_summary',
             'summarize_acmg_criteria'}
MAP_VALUE_EVENTS = []
EXTENSION_VALUE_EVENTS = []
EXTENSION_CALLBACK_ACTIVE = False
EXTENSION_CALLBACK_FAIL = False


class CallbackEvidence:
    """Unsupported evidence whose original isinstance calls are observable."""
    def __init__(self, fail=False):
        self.fail, self.calls = fail, 0

    def __reduce__(self):
        return type(self), (self.fail,)

    @property
    def __class__(self):
        self.calls += 1
        MAP_VALUE_EVENTS.append(self.calls)
        if self.calls == 2:
            warnings.warn('custom map callback', UserWarning)
            if self.fail:
                raise RuntimeError('custom map callback failed')
        return object


class CallbackIntegerArray(pd.arrays.IntegerArray):
    """Unknown storage whose dtype callbacks belong to the original function."""
    @property
    def dtype(self):
        if EXTENSION_CALLBACK_ACTIVE:
            EXTENSION_VALUE_EVENTS.append(len(EXTENSION_VALUE_EVENTS) + 1)
            if len(EXTENSION_VALUE_EVENTS) == 13:
                warnings.warn('extension dtype callback', UserWarning)
                if EXTENSION_CALLBACK_FAIL:
                    raise RuntimeError('extension dtype callback failed')
        return super().dtype


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def write_json(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, sort_keys=True, indent=2); stream.write('\n')


def fingerprint(value):
    # Observation must not invoke this fixture's production-only callback. The
    # original function still sees every real dtype access, including errors.
    global EXTENSION_CALLBACK_ACTIVE
    previous = EXTENSION_CALLBACK_ACTIVE
    EXTENSION_CALLBACK_ACTIVE = False
    try:
        return digest(json.dumps(encode(value), sort_keys=True, separators=(',', ':')).encode())
    finally:
        EXTENSION_CALLBACK_ACTIVE = previous



def publication_stores(function):
    """Locate actual publication STORE_SUBSCR instructions by exact target span.

    A multiline RHS can emit the assignment's first line twice before its
    STORE. A line event alone cannot identify completed frame publication.
    """
    source = Path(function.__code__.co_filename).read_text()
    definitions = [n for n in ast.parse(source).body
                   if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
                   and n.name == function.__name__]
    assert len(definitions) == 1 and type(definitions[0]) is ast.FunctionDef
    assert not definitions[0].decorator_list
    instructions = list(dis.get_instructions(function))
    stores = {}
    for node in ast.walk(definitions[0]):
        if not isinstance(node, ast.Assign) or len(node.targets) != 1: continue
        target = node.targets[0]
        label = ast.unparse(target)
        if label == "df['ACMG_criteria']": kind = 'summary'
        elif label == 'criteria_matrix.loc[overridden, col]': kind = 'column'
        else: continue
        span = (target.lineno, target.end_lineno, target.col_offset, target.end_col_offset)
        matches = [i for i in instructions
                   if i.opname == 'STORE_SUBSCR' and tuple(i.positions) == span]
        assert len(matches) == 1 and matches[0].offset not in stores
        stores[matches[0].offset] = kind
    assert set(stores.values()) == {'summary', 'column'}
    return stores


def source_proof(reference, candidate):
    old = (reference / 'scripts/acmg_scoring.py').read_text()
    new = (candidate / 'scripts/acmg_scoring.py').read_text()
    assert digest(old.encode()) == REFERENCE_SHA
    old_segment, restored = historical_summary_source(old), historical_summary_source(new)
    assert old_segment == restored
    function = next(n for n in ast.parse(new).body if isinstance(n, ast.FunctionDef)
                    and n.name == 'summarize_acmg_criteria')
    candidate_segment = ast.get_source_segment(new, function)
    assert new.count(candidate_segment) == 1
    assert digest(new.replace(candidate_segment, restored, 1).encode()) == BASE_CANDIDATE_SHA
    generators = []
    for path in [reference / 'seed-source/differential_scoring.py', candidate / 'rust/tests/differential_scoring.py']:
        text = path.read_text()
        definitions = [n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef) and n.name == 'rows_for_fixture']
        assert len(definitions) == 1
        generators.append(ast.get_source_segment(text, definitions[0]))
    assert generators[0] == generators[1]
    for name in FUNCTIONS - {'summarize_acmg_criteria'} | {'calculate_posterior_probability', 'compute_control_common'}:
        definitions = []
        for text in [old, new]:
            nodes = [n for n in ast.parse(text).body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name == name]
            assert len(nodes) == 1 and isinstance(nodes[0], ast.FunctionDef) and not nodes[0].decorator_list
            definitions.append(ast.get_source_segment(text, nodes[0]))
        assert definitions[0] == definitions[1], name
    return dict(original_whole_sha256=REFERENCE_SHA, restored_candidate_whole_sha256=BASE_CANDIDATE_SHA,
        candidate_whole_sha256=digest(new.encode()), exact_historical_summary_sha256=digest(restored.encode()),
        original_posterior_and_C1_C5_unchanged=True,
        original_differential_fixture_generator_sha256=digest(generators[0].encode()))


def reference_function(reference):
    path = reference / 'scripts/acmg_scoring.py'; text = path.read_text()
    assert digest(text.encode()) == REFERENCE_SHA
    namespace = dict(pd=pd, np=np, Dict=Dict, Tuple=Tuple, os=os, gzip=gzip, pickle=pickle,
                     __name__='acmg_scoring', __package__='', logger=logging.getLogger('acmg_scoring'))
    nodes = [n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef) and n.name in FUNCTIONS]
    assert len(nodes) == 4
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), namespace)
    return namespace['summarize_acmg_criteria'], namespace['logger']


def cases(candidate):
    def seed(name, count=3, **kw):
        data = pd.DataFrame({'variant_id': ['a', 'b', 'c'][:count],
            'untouched': ['β\ttab\nquote"', 'second', 'third'][:count],
            'ACMG_quant_score': np.array([.94, -0., np.nan])[:count],
            'ACMG_class': ['retain-class'] * count, 'ACMG_criteria': ['retain-until-summary'] * count})
        data.attrs = {'nested': {'caller': ['keep']}}; data.columns.name = 'original columns'
        values = {name: np.zeros(count, dtype=np.int64) for name in CRITERIA}
        if count:
            values['PVS1'][0] = 4; values['PM2'][:] = [1, 0, 2][:count]
        return dict(name=name, frame=data, criteria=values, mapping={}, enabled=True,
                    map_kind='plain', phases=True, **kw)
    result = [seed('empty_map'), seed('override_off_missing_map')]
    result[-1].update(enabled=False, map_kind='missing')
    case = seed('production_integer_arrays_and_default_unassigned_columns')
    for name in ('PS4', 'PP4'): del case['criteria'][name]
    case['mapping'] = {'a': 'PVS1_Very', 'b': np.nan}; result.append(case)
    for name, mapping in [
        ('curated_and_nonstring', {'a': 'PS1,PM2_Supporting', 'b': np.nan, 'c': 17, 'unused': None}),
        ('bare_and_first_suffix', {'a': 'PS1_Supporting,PS1', 'b': 'PS1_Supporting,PS1_Strong'}),
        ('split_and_unicode', {'a': 'PVS1_Very_extra, BA1_Standalone', 'b': '\x1c\x1f', 'c': ' \tPM2_Moderate, \x1cPP3\x1f '}),
        ('comma_only_and_unknown', {'a': ', ,', 'b': 'UNKNOWN', 'c': 'PS1;PM2_Supporting'}),
        ('no_matching_map_entries', {'not-in-frame': 'PS1'}),
        ('nonstring_keys', {1: 'PS1', None: 'PM2', 'a': 'PS1'}),
        ('nan_and_nonstring_values', {'a': np.nan, 'b': pd.NA, 'c': [1, 2]}),
        ('nan_text_is_explicit_code', {'a': 'NaN', 'b': 'None', 'c': ''}),
    ]:
        case = seed(name); case['mapping'] = mapping; result.append(case)
    for evidence in ['PS1_Unknown', 'XPS1', 'PS10', 'PS1_', 'PS1_ Supporting', 'PVS1_Unknown', 'BP7_Unknown']:
        case = seed('malformed_' + evidence); case['mapping'] = {'a': evidence}; result.append(case)
    case = seed('column_phase_precedes_row_phase'); case['mapping'] = {'a': 'XPS1', 'b': 'PVS1_Unknown'}; result.append(case)
    for kind in ['missing', 'gzip', 'bad_pickle']:
        case = seed('map_' + kind); case.update(map_kind=kind, mapping={'a': 'PS1'}); result.append(case)
    for fail in (False, True):
        case = seed('custom_map_callback_' + str(fail))
        case['mapping'] = {'a': CallbackEvidence(fail)}; result.append(case)
    for dtype in ['bool', 'boolean', 'int32', 'uint64', 'float32', 'float64', 'Int64', 'Float64', 'object']:
        for curated in [False, True]:
            case = seed('strength_' + dtype + ('_curated' if curated else '_unmapped'))
            case['criteria'] = {name: pd.array(np.minimum(value, 1) if dtype in ('bool', 'boolean') else value, dtype=dtype)
                                for name, value in case['criteria'].items()}
            if curated: case['mapping'] = {'a': 'PS1, PM2_Supporting'}
            result.append(case)
    for value in [6, -1, 1.5, np.nan, pd.NA, '1', complex(1, 0)]:
        case = seed('invalid_strength_' + repr(value)); case['criteria']['PVS1'] = [value, 0, 0]; result.append(case)
    case = seed('invalid_validation_before_missing_map'); case.update(map_kind='missing'); case['criteria']['PVS1'] = [6, 0, 0]; result.append(case)
    case = seed('empty', count=0); result.append(case)
    case = seed('nonmonotonic_unique_index'); case['frame'].index = pd.Index(['c', 'a', 'b'], name='rows'); result.append(case)
    case = seed('duplicate_index'); case['frame'].index = pd.Index(['same', 'other', 'same'], name='rows'); case['mapping'] = {'a': 'PS1'}; result.append(case)
    case = seed('multiindex'); case['frame'].index = pd.MultiIndex.from_tuples([('b', 2), ('a', 1), ('b', 2)], names=['x', 'y']); result.append(case)
    case = seed('flags_false'); case['frame'].flags.allows_duplicate_labels = False; result.append(case)
    case = seed('integer_variant_ids'); case['frame']['variant_id'] = [2**53 + 1, 1, 2]; case['mapping'] = {2**53 + 1: 'PS1', 1: 'PM2'}; result.append(case)
    case = seed('missing_variant_id'); case['frame'] = case['frame'].drop(columns=['variant_id']); result.append(case)
    case = seed('duplicate_variant_columns'); case['frame'] = pd.concat([case['frame'], case['frame'][['variant_id']]], axis=1); result.append(case)
    case = seed('misaligned_series_strengths'); case['criteria']['PVS1'] = pd.Series([4, 1, 0], index=[2, 0, 1]); case['frame'].index = [2, 0, 1]; result.append(case)
    case = seed('nullable_missing_strength'); case['criteria']['PVS1'] = pd.array([4, pd.NA, 0], dtype='Int64'); result.append(case)
    case = seed('signed_zero_strengths'); case['criteria']['PVS1'] = np.array([-0., 0., 4.]); result.append(case)
    spec = importlib.util.spec_from_file_location('source_matched_scoring_seeds', candidate / 'rust/tests/differential_scoring.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    bulk, evidence = module.rows_for_fixture()
    for enabled in [False, True]:
        result.append(dict(name='historical_bulk_override_' + str(enabled), frame=bulk,
            criteria={name: bulk[name].to_numpy() for name in CRITERIA}, mapping=evidence,
            enabled=enabled, map_kind='plain', phases=False))
    # Append these cases so every earlier prepared case keeps its ordinal.
    for fail in (False, True):
        case = seed('custom_extension_callback_' + str(fail), count=1)
        case['enabled'] = False
        case['criteria']['PVS1'] = CallbackIntegerArray(np.array([1], dtype='int64'), np.zeros(1, dtype=bool))
        case['extension_callback_fail'] = fail
        result.append(case)
    assert len({case['name'] for case in result}) == len(result)
    return result


def prepare(args):
    proof = source_proof(args.reference, args.candidate)
    args.output.mkdir(parents=True, exist_ok=False)
    inventory = []
    for i, case in enumerate(cases(args.candidate)):
        folder = args.output / f'{i:03d}'; folder.mkdir()
        path = folder / ('map.pkl.gz' if case['map_kind'] == 'gzip' else 'map.pkl')
        if case['map_kind'] == 'bad_pickle': path.write_bytes(b'not a pickle')
        elif case['map_kind'] != 'missing':
            with (gzip.open if case['map_kind'] == 'gzip' else open)(path, 'wb') as stream:
                pickle.dump(case['mapping'], stream, protocol=5)
        case['map_path'] = str(path)
        raw = pickle.dumps(case, protocol=5)
        (folder / 'case.pkl').write_bytes(raw)
        inventory.append(dict(index=i, name=case['name'], rows=len(case['frame']),
            phase_trace=case['phases'], sha256=digest(raw),
            map_sha256=digest(path.read_bytes()) if path.exists() else None))
    write_json(args.output / 'inputs.json', dict(source_proof=proof, cases=inventory))
    print(json.dumps(dict(cases=len(inventory), rows=sum(c['rows'] for c in inventory), source_proof=proof)))


def execute(function, logger, case, mode, binary, bridge, artifact):
    global EXTENSION_CALLBACK_ACTIVE, EXTENSION_CALLBACK_FAIL
    EXTENSION_CALLBACK_ACTIVE = False
    value = pickle.loads(pickle.dumps(case, protocol=5))
    data, values = value['frame'], value['criteria']
    before_values = fingerprint(values)
    path = Path(value['map_path']); map_before = path.read_bytes() if path.exists() else None
    MAP_VALUE_EVENTS.clear()
    EXTENSION_VALUE_EVENTS.clear()
    EXTENSION_CALLBACK_FAIL = value.get('extension_callback_fail', False)
    logs, completed, pending, final_matrix = [], [], [None], [None]
    attempts, native_successes = [0], [0]
    real = bridge.try_criteria_summary
    def dispatch(*args, **kwargs):
        attempts[0] += 1
        plan = None if mode == 'fallback' else real(*args, **kwargs)
        native_successes[0] += plan is not None
        return plan
    class Log(logging.Handler):
        def emit(self, record):
            logs.append([record.levelname, record.name, record.getMessage()])
    handler = Log(); previous_level = logger.level; logger.addHandler(handler); logger.setLevel(logging.INFO)
    stores = publication_stores(function)
    def state(locals_):
        return dict(df=fingerprint(locals_['df']),
            matrix=fingerprint(locals_['criteria_matrix']) if 'criteria_matrix' in locals_ else None)
    def trace(pyframe, event, arg):
        if pyframe.f_code is not function.__code__: return None
        if event == 'call':
            pyframe.f_trace_opcodes = True
            pyframe.f_trace_lines = False
        elif event == 'opcode':
            # This event precedes the current instruction. A pending STORE has
            # therefore completed; exceptions clear it before any publication.
            if pending[0] is not None:
                if value['phases']: completed.append(dict(phase=pending[0], **state(pyframe.f_locals)))
                pending[0] = None
            kind = stores.get(pyframe.f_lasti)
            if kind is not None:
                pending[0] = 'summary' if kind == 'summary' else 'column:' + pyframe.f_locals['col']
        elif event == 'exception': pending[0] = None
        elif event == 'return':
            if 'criteria_matrix' in pyframe.f_locals: final_matrix[0] = pyframe.f_locals['criteria_matrix']
        return trace
    returned, error, raw_traceback = None, None, None
    try:
        with warnings.catch_warnings(record=True) as seen, patch.object(bridge, 'try_criteria_summary', dispatch):
            warnings.simplefilter('always')
            previous_trace = sys.gettrace()
            observer_frame = sys._getframe()
            previous_opcodes = observer_frame.f_trace_opcodes
            # Python 3.12 enables opcode events only if an existing frame has
            # opted in before settrace. Restore this observer flag afterwards.
            observer_frame.f_trace_opcodes = True
            sys.settrace(trace)
            try:
                kwargs = {} if mode in ('reference', 'default') else dict(summary_binary=binary)
                EXTENSION_CALLBACK_ACTIVE = True
                returned = function(data, values, str(path), value['enabled'], **kwargs)
            except Exception as exc:
                error = dict(type=type(exc).__name__, message=str(exc), args=encode(exc.args))
                raw_traceback = traceback.format_exc()
            finally:
                EXTENSION_CALLBACK_ACTIVE = False
                observer_frame.f_trace_opcodes = previous_opcodes
                sys.settrace(previous_trace)
            warning_events = [[type(w.message).__name__, str(w.message)] for w in seen]
    finally:
        logger.removeHandler(handler); logger.setLevel(previous_level)
    assert fingerprint(values) == before_values
    assert (path.read_bytes() if path.exists() else None) == map_before
    outcome = dict(error=error, warnings=warning_events, logs=logs, checkpoints=completed,
        df=fingerprint(data), matrix=fingerprint(final_matrix[0]),
        custom_map_callbacks=list(MAP_VALUE_EVENTS),
        custom_extension_callbacks=list(EXTENSION_VALUE_EVENTS),
        returned_same_df=returned is not None and returned[0] is data,
        returned_frames=None if returned is None else [fingerprint(x) for x in returned],
        criteria_inputs_unchanged=True, map_bytes_unchanged=True)
    with artifact.open('xb') as stream:
        pickle.dump(dict(outcome=outcome, frame=data, matrix=final_matrix[0], returned=returned,
                         raw_traceback=raw_traceback, checkpoints=completed), stream, protocol=5)
    return dict(name=value['name'], rows=len(data), outcome=outcome,
                native_attempts=attempts[0], native_successes=native_successes[0])


def worker(args):
    proof = source_proof(args.reference, args.candidate)
    sys.path.insert(0, str(args.candidate / 'scripts'))
    import acmg_scoring
    import priva_rust_criteria_summary as bridge
    original, original_logger = reference_function(args.reference)
    assert original.__module__ == acmg_scoring.__name__ == 'acmg_scoring'
    assert original.__globals__['__package__'] == acmg_scoring.__package__ == ''
    assert original_logger.name == acmg_scoring.logger.name == 'acmg_scoring'
    index = json.loads((args.inputs / 'inputs.json').read_bytes())
    args.output.mkdir(parents=True, exist_ok=False)
    reports = []
    for case in index['cases']:
        path = args.inputs / f"{case['index']:03d}" / 'case.pkl'; raw = path.read_bytes()
        assert digest(raw) == case['sha256']
        value = pickle.loads(raw)
        if args.profile == 'source-bootstrap':
            execute(original, original_logger, value, 'reference', args.binary, bridge,
                    args.output / f"{case['index']:03d}.bootstrap.pkl")
        function, logger = (original, original_logger) if args.mode == 'reference' else (acmg_scoring.summarize_acmg_criteria, acmg_scoring.logger)
        reports.append(execute(function, logger, value, args.mode, args.binary, bridge,
                               args.output / f"{case['index']:03d}.result.pkl"))
    write_json(args.output / 'result.json', dict(mode=args.mode, profile=args.profile,
        source_proof=proof, cases=reports, pandas=pd.__version__, numpy=np.__version__,
        execution_context=dict(module='acmg_scoring', package='', logger='acmg_scoring',
            historical_functions_compiled_in_original_order=True, ordered_logger_names_protected=True),
        native={k: os.environ.get(k) for k in ['SLURM_JOB_ID', 'SLURM_STEP_ID', 'SLURMD_NODENAME']}))
    print(json.dumps(dict(mode=args.mode, profile=args.profile, cases=len(reports),
        native_attempts=sum(r['native_attempts'] for r in reports),
        native_successes=sum(r['native_successes'] for r in reports))))


def compare(args):
    mismatches, comparisons, rows, native_successes = [], 0, 0, 0
    for profile in ['fresh', 'source-bootstrap']:
        baseline = json.loads((args.workers / (profile + '-reference') / 'result.json').read_bytes())
        for mode in ['default', 'native', 'fallback']:
            actual = json.loads((args.workers / (profile + '-' + mode) / 'result.json').read_bytes())
            assert len(baseline['cases']) == len(actual['cases'])
            for old, new in zip(baseline['cases'], actual['cases']):
                assert old['name'] == new['name']
                comparisons += 1; rows += old['rows']; native_successes += new['native_successes']
                if old['outcome'] != new['outcome']:
                    mismatches.append(dict(profile=profile, mode=mode, case=old['name'],
                        reference=old['outcome'], actual=new['outcome']))
                if mode in ('default', 'fallback'): assert new['native_successes'] == 0
                if mode == 'default': assert new['native_attempts'] == 0
    result = dict(passed=not mismatches, exact_whole_function_comparisons=comparisons,
        compared_rows=rows, native_successful_plans=native_successes,
        mismatch_count=len(mismatches), mismatches=mismatches,
        profiles=['fresh', 'source-bootstrap'], modes=['default', 'native', 'fallback'],
        posterior_or_criteria_assignment_or_ranking_executed=False)
    write_json(args.output, result)
    print(json.dumps({k: v for k, v in result.items() if k != 'mismatches'}))
    assert not mismatches, 'all mismatches retained in comparison result'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='procedure', required=True)
    for name in ['prepare', 'worker']:
        p = sub.add_parser(name)
        for option in ['candidate', 'reference', 'output']:
            p.add_argument('--' + option, type=Path, required=True)
        if name == 'worker':
            p.add_argument('--inputs', type=Path, required=True)
            p.add_argument('--binary', type=Path, required=True)
            p.add_argument('--mode', choices=['reference', 'default', 'native', 'fallback'], required=True)
            p.add_argument('--profile', choices=['fresh', 'source-bootstrap'], required=True)
    p = sub.add_parser('compare')
    for option in ['inputs', 'workers', 'output']: p.add_argument('--' + option, type=Path, required=True)
    args = parser.parse_args(); globals()[args.procedure](args)


if __name__ == '__main__':
    main()
