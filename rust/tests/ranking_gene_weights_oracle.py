#!/usr/bin/env python3
"""Capture unchanged R5 then R6, original ranker L599-616, before Rust work.

capture --full-oracle ROOT --compatibility-oracle ROOT --output NEW_DATA_ROOT
inspect --oracle ROOT
Observe real full-ranker L599/L621 buffers and verbatim source excerpts. Keep
actual pandas parsing, symbol/set membership, three separate multiplications,
partial mutation, warnings and first failure. No earlier-stage replay in the
stage-only current-buffer entry point; no scientific or production changes.
"""
import argparse
import ast
import copy
import csv
import gzip
import hashlib
import importlib.util
import json
import logging
import os
from pathlib import Path
import pickle
import struct
import sys
import warnings

import numpy as np
import pandas as pd

D = Path('/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade')
FULL_PINS = '9e63fb73e49c9f8ce700b0abd3e588702d2489c5430059e76c5d5b723eabb46a'
FULL_CAPTURE = 'db919b02f8738da60f85d16e0a7ed2bb53f8eaaf6047ed8f4db912fe4938508c'
R4_MANIFEST = '0b8fe35cad05d9c960423430a0fa149844c9e62d66c4e69d14fd6eb8c2ce653e'
FUNCTION_SHA = '8aa29d9a4f83ff5629bbdf5e6279a43a374d71ca87fe04fe9db88253cc57d2ca'
INPUTS = ['df', 'dispensable_gene_list', 'relevant_gene_list', 'compatibility_penalty_rows']
LOCALS = ['dispensable_genes', 'relevant_genes']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True, allow_nan=False)
        stream.write('\n')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def excerpt(path, scoring):
    text = path.read_text()
    source = next(node for node in ast.parse(text).body
                  if isinstance(node, ast.FunctionDef) and node.name == 'sort_and_rank_variants')
    assert hashlib.sha256(ast.get_source_segment(text, source).encode()).hexdigest() == FUNCTION_SHA
    body = [copy.deepcopy(node) for node in source.body if 599 <= node.lineno <= 616]
    assert [node.lineno for node in body] == [599, 606, 613] and body[-1].end_lineno == 616
    result = ast.Return(value=ast.Call(func=ast.Name(id='locals', ctx=ast.Load()), args=[], keywords=[]),
                        lineno=617, col_offset=4, end_lineno=617, end_col_offset=19)
    function = ast.FunctionDef(name='original_R5_R6_excerpt', args=ast.arguments(posonlyargs=[],
        args=[ast.arg(arg=name) for name in INPUTS], kwonlyargs=[], kw_defaults=[], defaults=[]),
        body=body + [result], decorator_list=[], lineno=599, col_offset=0, end_lineno=617, end_col_offset=19)
    namespace = dict(scoring.__dict__)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[])), str(path), 'exec'), namespace)
    return namespace['original_R5_R6_excerpt'], dict(lines=[599, 616], full_function_sha256=FUNCTION_SHA,
        statement_ast_sha256=hashlib.sha256(ast.dump(ast.Module(body=body, type_ignores=[])).encode()).hexdigest())


def selected(values):
    return {name: values[name] for name in INPUTS + LOCALS if name in values}


class Logs(logging.Handler):
    def __init__(self, helper):
        super().__init__()
        self.helper, self.records = helper, []

    def emit(self, record):
        self.records.append(dict(level=record.levelname, message=record.getMessage(),
                                 template=record.msg, args=self.helper.encode(record.args)))


def run_current(inputs, original, scoring, helper, entry_observation=None):
    """Execute only original L599-616 on these exact current buffers."""
    entry = helper.encode(inputs) if entry_observation is None else entry_observation
    before = pickle.dumps(inputs, protocol=5)
    events = []

    def trace(frame, event, arg):
        if frame.f_code is original.__code__ and (
                event in ['exception', 'return'] or event == 'line' and frame.f_lineno in [606, 613]):
            events.append(dict(event=event, line=frame.f_lineno,
                               snapshot=pickle.dumps(selected(frame.f_locals), protocol=5)))
        return trace

    observed = Logs(helper)
    scoring.logger.addHandler(observed)
    previous = sys.gettrace()
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            try:
                sys.settrace(trace)
                result = original(*inputs)
            except Exception as error:
                exception, result = dict(type=type(error).__name__, message=str(error)), None
            else:
                exception = None
            finally:
                sys.settrace(previous)
    finally:
        scoring.logger.removeHandler(observed)
    after = pickle.dumps(inputs, protocol=5)
    output = selected(result) if result is not None else None
    info = dict(exception=exception, input_before_sha256=hashlib.sha256(before).hexdigest(),
        input_after_sha256=hashlib.sha256(after).hexdigest(), input_pickle_unchanged=before == after,
        input_before=entry, input_after=helper.encode(inputs), output=helper.encode(output),
        source_logs=observed.records,
        warnings=[dict(category=item.category.__name__, message=str(item.message), source_line=item.lineno) for item in caught],
        source_events=[dict(event=item['event'], line=item['line'],
                            locals=helper.encode(pickle.loads(item['snapshot']))) for item in events],
        first_failure_line=next((item['line'] for item in events if item['event'] == 'exception'), None),
        returned_same_frame=output is not None and output['df'] is inputs[0])
    return info, before, after, pickle.dumps(output, protocol=5) if output is not None else None


def cell(value):
    if value is None:
        return ['null', '']
    if value is pd.NA:
        return ['pd_na', '']
    if isinstance(value, (bool, np.bool_)):
        return ['bool', str(int(value))]
    if isinstance(value, (int, np.integer)):
        return ['uint' if value > 2**63 - 1 else 'int', str(value)]
    if isinstance(value, (float, np.floating)):
        return ['float', struct.pack('>d', value).hex()]
    if isinstance(value, str):
        return ['str', value.encode().hex()]
    return ['unsupported', type(value).__module__ + '.' + type(value).__qualname__]


def export(root, name, inputs, original, scoring, helper, origin, expected=None, expected_error=None):
    folder = root / name
    folder.mkdir()
    prepared = pickle.dumps(inputs, protocol=5)
    (folder / 'prepared-input.pkl').write_bytes(prepared)
    result, before, after, output = run_current(pickle.loads(prepared), original, scoring, helper)
    repeated, rb, ra, _ = run_current(pickle.loads(prepared), original, scoring, helper)
    assert result == repeated and (before, after) == (rb, ra), name
    (folder / 'input.pkl').write_bytes(before)
    (folder / 'input-after.pkl').write_bytes(after)
    write(folder / 'result.json', result)
    if output is not None:
        (folder / 'expected.pkl').write_bytes(output)
    frame, dispensable, relevant, mask = pickle.loads(before)
    with (folder / 'input.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value'] + [item for label in frame.columns for item in [str(label) + '_kind', str(label)]])
        for row in range(len(frame)):
            writer.writerow(['opaque', helper.canonical(helper.encode(frame.index[row])).hex()] + [
                item for col in range(frame.shape[1]) for item in cell(frame.iloc[row, col])])
    with (folder / 'penalty.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value', 'value_kind', 'value', 'raw_data', 'is_null'])
        if isinstance(mask, pd.Series):
            for row, value in enumerate(mask):
                array = mask.array
                raw = bool(array._data[row]) if hasattr(array, '_mask') else None
                missing = bool(array._mask[row]) if hasattr(array, '_mask') else False
                writer.writerow(['opaque', helper.canonical(helper.encode(mask.index[row])).hex(), *cell(value),
                                 '' if raw is None else int(raw), int(missing)])
    report = dict(name=name, origin=origin, rows=len(frame), columns=list(frame.columns), dtypes=[str(x) for x in frame.dtypes],
        index_kind='multi' if isinstance(frame.index, pd.MultiIndex) else 'flat', penalty_type=type(mask).__name__,
        penalty_dtype=str(mask.dtype) if isinstance(mask, pd.Series) else None,
        penalty_index_equals_frame=isinstance(mask, pd.Series) and mask.index.equals(frame.index),
        gene_arguments=[helper.encode(dispensable), helper.encode(relevant)],
        original_exception=result['exception'], first_failure_line=result['first_failure_line'],
        input_before_sha256=result['input_before_sha256'], input_after_sha256=result['input_after_sha256'],
        source_input_pickle_unchanged=result['input_pickle_unchanged'], warnings=result['warnings'], source_logs=result['source_logs'])
    if expected is not None:
        assert result['exception'] is None
        actual = pickle.loads(output)
        assert {key: helper.encode(value) for key, value in actual.items() if key in ['df'] + LOCALS} == expected, name
        report['full_source_R7_entry_exact'] = True
    if expected_error is not None:
        assert result['exception'] == expected_error, (name, result['exception'], expected_error)
        report['full_source_R5_R6_failure_exact'] = True
    return report


def base(n=4):
    return pd.DataFrame({'SYMBOL': ['A', 'B', 'C', 'A'][:n], 'sort_index': [0.95, 0.475, -0., 0.12345678901234568][:n],
                         'zygosity_inheritance_mechanism_compatibility': [.9, 1., .9, .9][:n], 'untouched': ['u'] * n})


def extras(root):
    files = root / 'gene-lists'
    files.mkdir()
    contents = {'simple': 'A\nB\n', 'other': 'B\nC\n', 'empty': '', 'comments': '#SYMBOL\n# ignored\n\n',
        'spaces': ' A\nA \n \t\n', 'case': 'a\nA\nA\nSYMBOL\n', 'inline_comment': 'A#tail\nB\n',
        'quoted': '"A#B"\n"C"\n', 'tabs': 'A\tB\nC\tD\n', 'bad_rows': 'A\nB\tC\n',
        'numeric': '1\n2\n', 'leading_zero': '01\n002\n', 'float': '1e3\n-inf\ninf\nNaN\n',
        'na': 'NA\nN/A\nNaN\nnull\nNone\n<NA>\n', 'mixed_na': 'A\nNA\n',
        'unicode': 'β\nＧ\nA\n', 'crlf': 'A\r\nB\r\n', 'bom': '\ufeffA\nB\n',
        'quoted_newline': '"A\nB"\nC\n', 'whitespace_only': ' \n\t\n\n'}
    paths = {}
    for name, value in contents.items():
        path = files / (name + '.txt')
        path.write_bytes(value.encode())
        paths[name] = str(path)
    bad = files / 'invalid-utf8.txt'
    bad.write_bytes(b'A\n\xff\n')
    paths['invalid_utf8'] = str(bad)
    compressed = files / 'simple.txt.gz'
    compressed.write_bytes(gzip.compress(b'A\nB\n', mtime=0))
    paths['gzip'] = str(compressed)
    missing = str(files / 'missing.txt')
    def mask(frame):
        return pd.Series([i % 2 == 0 for i in range(len(frame))], index=frame.index, name='ACMG_quant_score')
    for name, first, second in [('none', None, None), ('both', paths['simple'], paths['other']),
        ('dispensable_only', paths['simple'], None), ('relevant_only', None, paths['other']),
        ('same_file', paths['simple'], paths['simple']), ('dispensable_missing', missing, paths['simple']),
        ('relevant_missing_after_dispensable', paths['simple'], missing), ('empty_path', '', None),
        ('false_path', False, None), ('integer_path', 7, None), ('pdNA_path', pd.NA, None),
        ('path_object', Path(paths['simple']), None)]:
        frame = base()
        yield name, (frame, first, second, mask(frame))
    for name, path in paths.items():
        frame = base()
        yield 'parse_' + name, (frame, path, paths['other'], mask(frame))
    for dtype, values in [('object', ['A', None, np.nan, pd.NA]), ('string', ['A', None, 'NA', 'B']),
        ('object_mixed', ['A', 1, True, 1.]), ('int64', [1, 2, 3, 1]), ('float64', [1., np.nan, 2., 0.]),
        ('Int64', [1, None, 2, 0]), ('boolean', [True, None, False, True]), ('category', ['A', None, 'B', 'C'])]:
        frame = base()
        frame['SYMBOL'] = pd.Series(values, dtype='object' if dtype == 'object_mixed' else dtype)
        for list_name in ['simple', 'numeric', 'na', 'mixed_na']:
            yield 'symbols_' + dtype + '_list_' + list_name, (frame.copy(), paths[list_name], paths['other'], mask(frame))
    raw_nan = [struct.unpack('>d', bytes.fromhex(value))[0] for value in
               ['7ff8000000000001', 'fff8000000001234', '7ff0000000000001', 'fff0000000000001']]
    scores = [('float64', [np.inf, -np.inf, np.nan, -0.]), ('float64_payloads', raw_nan),
        ('float32', [.95, .475, -0., np.nan]), ('int64', [1, 0, -1, 2]), ('uint64', [2**64-1, 0, 1, 2]),
        ('bool', [True, False, True, False]), ('Float64', [.95, None, -0., np.inf]),
        ('Int64', [1, None, 0, -1]), ('object', [.95, None, 'a', True]), ('string', ['.95', None, '0', '1'])]
    for dtype, values in scores:
        frame = base()
        frame['sort_index'] = pd.Series(values, dtype='float64' if dtype == 'float64_payloads' else dtype)
        yield 'scores_' + dtype, (frame, paths['simple'], paths['other'], mask(frame))
    for name, values in [('nan', raw_nan), ('infinity', [np.inf, -np.inf, -0., 1.]),
        ('object', ['x', .9, None, True]), ('nullable', [.9, None, .9, 1.])]:
        frame = base()
        frame['zygosity_inheritance_mechanism_compatibility'] = pd.Series(values, dtype='Float64' if name == 'nullable' else None)
        yield 'compatibility_' + name, (frame, None, None, mask(frame))
    frame = base()
    frame['sort_index'] = raw_nan
    frame['zygosity_inheritance_mechanism_compatibility'] = raw_nan[::-1]
    yield 'two_nan_operands', (frame, None, None, mask(frame))
    for name in ['SYMBOL', 'sort_index', 'zygosity_inheritance_mechanism_compatibility']:
        for active in [True, False]:
            frame = base().drop(columns=[name])
            yield 'missing_' + name + '_' + str(active), (frame, paths['simple'] if active else None, None, mask(frame))
        frame = base()
        frame = pd.concat([frame, frame[[name]]], axis=1)
        yield 'duplicate_' + name, (frame, paths['simple'], paths['other'], mask(frame))
    for target in ['dispensable_gene_index', 'relevant_gene_index']:
        for dtype, values in [('float64', [7.] * 4), ('int64', [7] * 4), ('string', ['old'] * 4), ('Float64', [7., None, 7., 7.])]:
            frame = base()
            frame.insert(1, target, pd.Series(values, dtype=dtype))
            yield 'existing_' + target + '_' + dtype, (frame, paths['simple'], paths['other'], mask(frame))
            yield 'disabled_' + target + '_' + dtype, (frame.copy(), None, None, mask(frame))
        frame = base()
        frame[target] = 7.
        frame = pd.concat([frame, frame[[target]]], axis=1)
        yield 'duplicate_' + target, (frame, paths['simple'], paths['other'], mask(frame))
    for name, index in [('duplicate', pd.Index([7, 7, -2, 7], name='duplicate')),
        ('gapped', pd.Index([8, -2, 90, 4], name='gapped')),
        ('multi', pd.MultiIndex.from_tuples([('a', 1), ('a', 1), ('b', 2), ('a', 3)], names=['left', 'right']))]:
        frame = base()
        frame.index = index
        frame.attrs = {'scope': 'R5-R6 fixture'}
        frame.columns.name = 'column-axis'
        selected_mask = mask(frame)
        selected_mask.attrs = {'scope': 'R4 mask'}
        yield 'axis_' + name, (frame, paths['simple'], paths['other'], selected_mask)
    for dtype, values in [('bool', [True, False, True, False]), ('boolean', [True, None, False, None]),
        ('object', [True, False, None, pd.NA]), ('int64', [1, 0, 1, 0]), ('string', ['1', None, '0', 'True'])]:
        frame = base()
        yield 'penalty_dtype_' + dtype, (frame, paths['simple'], paths['other'], pd.Series(values, dtype=dtype))
    for bits in [[False] * 4, [True] * 4]:
        frame = base()
        selected_mask = pd.Series(pd.arrays.BooleanArray(np.array(bits), np.array([False, True, False, True])))
        yield 'penalty_masked_raw_' + str(int(bits[0])), (frame, paths['simple'], paths['other'], selected_mask)
    for index, name in [(pd.Index([3, 2, 1, 0]), 'reordered'), (pd.Index([0, 1, 2, 7]), 'missing_extra'),
        (pd.Index([0, 0, 1, 2]), 'different_duplicates')]:
        frame = base()
        yield 'unaligned_penalty_' + name, (frame, paths['simple'], paths['other'], pd.Series([True] * 4, index=index))
    for size in [0, 1]:
        for enabled in [True, False]:
            frame = base(size)
            yield 'size_' + str(size) + '_' + str(enabled), (frame, paths['simple'] if enabled else None,
                                                            paths['other'] if enabled else None, mask(frame))
    frame = base()
    frame.flags.allows_duplicate_labels = False
    selected_mask = mask(frame)
    selected_mask.flags.allows_duplicate_labels = False
    yield 'attrs_flags', (frame, paths['simple'], paths['other'], selected_mask)


def capture(args):
    native = {key: os.environ.get(key) for key in ['SLURM_JOB_ID', 'SLURM_STEP_ID', 'SLURMD_NODENAME']}
    assert native['SLURM_JOB_ID'] and native['SLURM_STEP_ID']
    for path in [args.full_oracle, args.compatibility_oracle, args.output]:
        assert path.resolve().is_relative_to(D)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'capture-harness.py').write_bytes(Path(__file__).read_bytes())
    write(args.output / 'attempt.json', dict(native_environment=native, harness_sha256=sha(Path(__file__))))
    assert sha(args.full_oracle / 'source-pins.json') == FULL_PINS and sha(args.full_oracle / 'capture.json') == FULL_CAPTURE
    assert sha(args.compatibility_oracle / 'manifest.json') == R4_MANIFEST
    for name, wanted in json.loads((args.compatibility_oracle / 'manifest.json').read_text()).items():
        assert sha(args.compatibility_oracle / name) == wanted
    helper = load('frozen_full_R5_R6_oracle', args.full_oracle / 'source/ranking_control_oracle.py')
    helper.verify_sources(args.full_oracle, FULL_PINS)
    scoring, _, _, _ = helper.original_modules(args.full_oracle / 'source/scripts')
    original, source = excerpt(args.full_oracle / 'source/scripts/acmg_scoring.py', scoring)
    write(args.output / 'runtime-identities.json', helper.runtime_identity())
    captured = json.loads((args.full_oracle / 'capture.json').read_text())
    reports, preceding = [], []
    for case in captured['cases']:
        if case['kind'] != 'ranking':
            continue
        folder = args.full_oracle / case['folder']
        r4_boundary = args.compatibility_oracle / ('full_' + case['name']) / 'actual-full-L599-boundary.pkl'
        if not r4_boundary.exists():
            preceding.append(dict(name=case['name'], first_failure=case['first_uncaught_failure']))
            continue
        data, pedigree = pickle.loads((folder / 'prepared-input.pkl').read_bytes())
        helper.encode(data)
        helper.encode(pedigree)
        bound = json.loads((folder / 'traced/arguments.json').read_text())['arguments']
        kwargs = {key['value']: helper.decode_atom(value) for key, value in bound['items'] if key['value'] not in {'df', 'ped_df'}}
        entries, boundaries = [], []
        def trace(frame, event, arg):
            if frame.f_code is scoring.sort_and_rank_variants.__code__ and event == 'line':
                if frame.f_lineno == 599:
                    entries.append(pickle.dumps(tuple(frame.f_locals[name] for name in INPUTS), protocol=5))
                if frame.f_lineno == 621:
                    boundaries.append(pickle.dumps(selected(frame.f_locals), protocol=5))
            return trace
        previous = sys.gettrace()
        with warnings.catch_warnings(record=True):
            warnings.simplefilter('always')
            try:
                sys.settrace(trace)
                full_output = scoring.sort_and_rank_variants(data, pedigree, **kwargs)
            except Exception as error:
                assert case['status'] == 'error'
                assert dict(type=type(error).__name__, message=str(error)) == dict(
                    type=case['exception']['type'].rsplit('.', 1)[-1], message=case['exception']['message'])
            else:
                assert case['status'] == 'value'
                assert helper.encode(full_output) == json.loads((folder / 'traced/outcome.json').read_text())['output'], case['name']
            finally:
                sys.settrace(previous)
        assert len(entries) == 1, case['name']
        entry = pickle.loads(entries[0])
        retained = pickle.loads(r4_boundary.read_bytes())
        assert helper.encode(entry[0]) == helper.encode(retained['df']), case['name']
        assert helper.encode(entry[3]) == helper.encode(retained['compatibility_penalty_rows']), case['name']
        expected, expected_error = None, None
        if boundaries:
            assert len(boundaries) == 1
            boundary = pickle.loads(boundaries[0])
            expected = {key: helper.encode(value) for key, value in boundary.items() if key in ['df'] + LOCALS}
            events = json.loads((folder / 'traced/trace.json').read_text())['events']
            following = next(event for event in events if event['function'] == 'sort_and_rank_variants' and event['line'] > 616)
            closed_state = json.loads((args.full_oracle / 'states' / (following['state'] + '.json')).read_text())
            assert expected == {key: closed_state[key] for key in expected}, case['name']
        else:
            assert case['first_uncaught_failure']['phase'] in ['R5', 'R6']
            expected_error = dict(type=case['exception']['type'].rsplit('.', 1)[-1], message=case['exception']['message'])
        reports.append(export(args.output, 'full_' + case['name'], entry, original, scoring, helper,
            dict(kind='actual-full-source-R5-entry', full_case=case['folder'], scientific_arguments_sha256=sha(folder / 'traced/arguments.json'),
                 retained_R4_boundary_sha256=sha(r4_boundary)), expected, expected_error))
        if boundaries:
            (args.output / ('full_' + case['name']) / 'actual-full-L621-boundary.pkl').write_bytes(boundaries[0])
    for name, inputs in extras(args.output):
        reports.append(export(args.output, 'extra_' + name, inputs, original, scoring, helper,
                              dict(kind='unchanged-original-R5-R6-statements')))
    write(args.output / 'oracle.json', dict(source_commit=helper.COMMIT, source=source, harness_sha256=sha(Path(__file__)),
        native_environment=native, full_oracle=str(args.full_oracle), full_source_pins_sha256=FULL_PINS, full_capture_sha256=FULL_CAPTURE,
        compatibility_oracle=str(args.compatibility_oracle), compatibility_oracle_manifest_sha256=R4_MANIFEST,
        cases=reports, preceding_full_source_failures=preceding,
        scope='R5 then R6 actual pandas parse/set membership and separate multiplication order on current buffers; original scientific source unchanged',
        scientific_acceptance=False, production_backend_selected=False))
    write(args.output / 'manifest.json', {str(path.relative_to(args.output)): sha(path) for path in sorted(args.output.rglob('*')) if path.is_file()})
    print(json.dumps(dict(cases=len(reports), outputs=sum(case['original_exception'] is None for case in reports),
        source_failures=sum(case['original_exception'] is not None for case in reports), preceding_failures=len(preceding),
        manifest_sha256=sha(args.output / 'manifest.json'), native_environment=native)), flush=True)


def inspect(args):
    assert args.oracle.resolve().is_relative_to(D)
    data = json.loads((args.oracle / 'oracle.json').read_text())
    for name, wanted in json.loads((args.oracle / 'manifest.json').read_text()).items():
        assert sha(args.oracle / name) == wanted
    print(json.dumps(dict(cases=len(data['cases']), full_source_cases=sum(case['origin']['kind'] == 'actual-full-source-R5-entry' for case in data['cases']),
        outputs=sum(case['original_exception'] is None for case in data['cases']),
        source_failures=[dict(name=case['name'], exception=case['original_exception'], line=case['first_failure_line'])
                         for case in data['cases'] if case['original_exception'] is not None],
        manifest_sha256=sha(args.oracle / 'manifest.json'))), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    command = commands.add_parser('capture')
    for name in ['full-oracle', 'compatibility-oracle', 'output']:
        command.add_argument('--' + name, type=Path, required=True)
    commands.add_parser('inspect').add_argument('--oracle', type=Path, required=True)
    args = parser.parse_args()
    (capture if args.command == 'capture' else inspect)(args)


if __name__ == '__main__':
    main()
