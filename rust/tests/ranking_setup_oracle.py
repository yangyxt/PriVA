#!/usr/bin/env python3
"""Capture unchanged original R1 copy/drop/assign/filter behavior before Rust.

capture --full-oracle CLOSED_ROOT --output NEW_DATA_ROOT
Four original source statements only, stopping before PED/C1/SciPy. Existing
full-ranker R2-entry checkpoints independently bind the boundary. No science edit.
"""
import argparse
import ast
import copy
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import struct

import numpy as np
import pandas as pd

D = Path('/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade')
FULL_PINS = '9e63fb73e49c9f8ce700b0abd3e588702d2489c5430059e76c5d5b723eabb46a'
FULL_CAPTURE = 'db919b02f8738da60f85d16e0a7ed2bb53f8eaaf6047ed8f4db912fe4938508c'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True, allow_nan=False); stream.write('\n')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def cell(value):
    if value is None: return ['null', '']
    if value is pd.NA: return ['pd_na', '']
    if isinstance(value, (bool, np.bool_)): return ['bool', str(int(value))]
    if isinstance(value, (int, np.integer)): return ['uint' if value > 2**63-1 else 'int', str(value)]
    if isinstance(value, (float, np.floating)): return ['float', struct.pack('>d', value).hex()]
    if isinstance(value, str): return ['str', value.encode().hex()]
    raise TypeError('outside primitive transport: ' + str(type(value)))


def original_r1(path):
    text = path.read_text()
    function = next(node for node in ast.parse(text).body if isinstance(node, ast.FunctionDef) and node.name == 'sort_and_rank_variants')
    body = [copy.deepcopy(node) for node in function.body if 492 <= node.lineno <= 498]
    assert len(body) == 4 and [node.lineno for node in body] == [492, 493, 497, 498]
    # Original instructions retain their filenames/line numbers. The wrapper's
    # only added operation returns the captured R1 boundary, before original R2.
    boundary = ast.Return(value=ast.Name(id='df', ctx=ast.Load()), lineno=499, col_offset=4, end_lineno=499, end_col_offset=13)
    wrapper = ast.FunctionDef(name='original_R1_excerpt', args=ast.arguments(posonlyargs=[], args=[ast.arg(arg='df')],
        kwonlyargs=[], kw_defaults=[], defaults=[]), body=body+[boundary], decorator_list=[], lineno=492, col_offset=0,
        end_lineno=499, end_col_offset=13)
    scope = {}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[wrapper], type_ignores=[])), str(path), 'exec'), scope)
    return scope['original_R1_excerpt'], dict(lines=[492, 498], statement_ast_sha256=hashlib.sha256(
        ast.dump(ast.Module(body=body, type_ignores=[])).encode()).hexdigest(),
        full_function_sha256=hashlib.sha256(ast.get_source_segment(text, function).encode()).hexdigest())


def extras(helper):
    values = ['protein_coding', 'protein_coding ', ' protein_coding', 'Protein_coding', 'PROTEIN_CODING',
              'protein_coding\n', '', None, pd.NA, np.nan, 0, 1, True, False, np.inf]
    frame = helper.base_frame(len(values)); frame['BIOTYPE'] = pd.Series(values, dtype=object)
    yield 'exact_token_and_primitive_object', frame
    for dtype, values in [('string', ['protein_coding', None, '', 'protein_coding\n']),
        ('Int64', [1, None, 0, -1]), ('UInt64', [2**64-1, None, 0, 2**63]),
        ('Float64', [1., None, np.inf, -0.]), ('boolean', [True, None, False, True]),
        ('int64', [1, 0, -1, 2]), ('uint64', [2**64-1, 0, 1, 2**63]),
        ('float64', [np.nan, np.inf, -0., 1.]), ('bool', [True, False, True, False])]:
        frame = helper.base_frame(4); frame['BIOTYPE'] = pd.Series(values, dtype=dtype)
        yield 'biotype_dtype_' + dtype, frame
    for dtype, values in [('object', ['x', None, np.nan, True]), ('string', ['x', None, '', 'y']),
        ('Int64', [2**63-1, None, -1, 0]), ('UInt64', [2**64-1, None, 0, 2**63]),
        ('Float64', [np.inf, None, -0., 1.]), ('float64', [np.inf, np.nan, -0., np.nextafter(1., 2.)])]:
        frame = helper.base_frame(4, BIOTYPE=['protein_coding', 'other', 'protein_coding', 'protein_coding'])
        frame['ACMG_quant_score'] = pd.Series(values, dtype=dtype)
        yield 'raw_score_copy_' + dtype, frame
    for keep in [[True]*4, [False]*4, [True, False, True, False]]:
        frame = helper.base_frame(4, BIOTYPE=['protein_coding' if x else 'other' for x in keep])
        for name, index in [('range', pd.RangeIndex(10, 18, 2, name='range')),
            ('duplicate', pd.Index([7, 7, -2, 100], name='gapped')),
            ('object', pd.Index(['x', None, 'x', ''], dtype=object, name='opaque')),
            ('multi', pd.MultiIndex.from_tuples([('x', 1), ('x', 1), ('z', -2), ('a', 9)], names=['one', 'two']))]:
            example = frame.copy(deep=True); example.index = index
            example.columns.name = 'column-axis'; example.attrs['retained'] = 'original attrs'
            yield 'index_' + name + '_' + ''.join(str(int(x)) for x in keep), example
    for n in [0, 1, 4]:
        frame = helper.base_frame(n)
        frame.insert(0, 'sort_index', ['stale']*n)
        frame.insert(1, 'haplo_insuf_index', ['legacy']*n)
        frame.insert(3, 'zygosity_inheritance_mechanism_index', ['legacy']*n)
        frame['variant_rank'] = 77
        yield 'column_positions_' + str(n), frame
    for missing in [['ACMG_quant_score'], ['BIOTYPE'], ['ACMG_quant_score','BIOTYPE']]:
        for n in [0, 1]:
            yield 'missing_' + '_'.join(missing) + '_' + str(n), helper.base_frame(n).drop(columns=missing)
    for name in ['ACMG_quant_score', 'BIOTYPE', 'sort_index', 'unused']:
        frame = helper.base_frame(2)
        if name not in frame: frame[name] = ['stale']*2
        frame = pd.concat([frame, frame[[name]]], axis=1)
        yield 'duplicate_' + name, frame


def export(output, name, frame, function, helper, origin, expected=None, expected_exception=None):
    folder = output / name; folder.mkdir()
    prepared = pickle.dumps(frame, protocol=5)
    (folder / 'input.pkl').write_bytes(prepared)
    write(folder / 'input.json', helper.encode(frame))
    original = pickle.loads(prepared)
    before = pickle.dumps(original, protocol=5)
    try:
        result = function(original)
    except Exception as error:
        result = None; exception = dict(type=type(error).__name__, message=str(error))
    else:
        exception = None
    after = pickle.dumps(original, protocol=5)
    (folder / 'input-after.pkl').write_bytes(after)
    write(folder / 'input-after.json', helper.encode(original))
    report = dict(name=name, rows=len(frame), columns=list(frame.columns), dtypes=[str(x) for x in frame.dtypes], origin=origin,
        original_exception=exception, input_semantically_unchanged=helper.encode(frame)==helper.encode(original),
        input_pickle_unchanged=before==after, input_before_sha256=hashlib.sha256(before).hexdigest(),
        input_after_sha256=hashlib.sha256(after).hexdigest())
    assert report['input_semantically_unchanged'], name
    with (folder / 'input.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind','index_value']+[x for col in frame.columns for x in [col+'_kind',col]])
        for row in range(len(frame)):
            writer.writerow(['opaque',helper.canonical(helper.encode(frame.index[row])).hex()]+[
                item for col in range(frame.shape[1]) for item in cell(frame.iloc[row,col])])
    if result is not None:
        (folder / 'expected.pkl').write_bytes(pickle.dumps(result, protocol=5))
        write(folder / 'expected.json', helper.encode(result))
        positions = {token: i for i,token in enumerate(frame['__oracle_row_position'])}
        report.update(permutation=[positions[token] for token in result['__oracle_row_position']],
            output_columns=list(result.columns), output_dtypes=[str(x) for x in result.dtypes])
        if expected is not None:
            assert helper.encode(result) == expected, name
            report['actual_full_ranker_checkpoint_exact'] = True
    if expected_exception is not None:
        assert exception == expected_exception, (name, exception, expected_exception)
        report['actual_full_ranker_R1_failure_exact'] = True
    return report


def capture(args):
    native = {name: os.environ.get(name) for name in ['SLURM_JOB_ID','SLURM_STEP_ID','SLURMD_NODENAME']}
    assert native['SLURM_JOB_ID'] and native['SLURM_STEP_ID']
    assert args.full_oracle.resolve().is_relative_to(D) and args.output.resolve().is_relative_to(D)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'capture-harness.py').write_bytes(Path(__file__).read_bytes())
    write(args.output / 'attempt.json', dict(native_environment=native, source_sha256=sha(Path(__file__))))
    assert sha(args.full_oracle / 'source-pins.json') == FULL_PINS and sha(args.full_oracle / 'capture.json') == FULL_CAPTURE
    helper = load('frozen_full_ranking_oracle', args.full_oracle / 'source/ranking_control_oracle.py')
    helper.verify_sources(args.full_oracle, FULL_PINS)
    function, source = original_r1(args.full_oracle / 'source/scripts/acmg_scoring.py')
    assert source['full_function_sha256'] == '8aa29d9a4f83ff5629bbdf5e6279a43a374d71ca87fe04fe9db88253cc57d2ca'
    write(args.output / 'runtime-identities.json', helper.runtime_identity())
    reports = []
    for case in json.loads((args.full_oracle / 'capture.json').read_text())['cases']:
        if case['kind'] != 'ranking': continue
        folder = args.full_oracle / case['folder']
        frame, _ = pickle.loads((folder / 'prepared-input.pkl').read_bytes())
        trace = json.loads((folder / 'traced/trace.json').read_text())
        events = [event for event in trace['events'] if event['function']=='sort_and_rank_variants' and event['phase']=='R2']
        expected, error = None, None
        if events:
            state = json.loads((args.full_oracle / 'states' / (events[0]['state']+'.json')).read_text())
            expected = state['df']
        elif case['first_uncaught_failure'] and case['first_uncaught_failure']['phase']=='R1':
            error = dict(type=case['exception']['type'].rsplit('.',1)[-1], message=case['exception']['message'])
        else:
            raise AssertionError('No actual R1 boundary or R1 failure: '+case['name'])
        reports.append(export(args.output, 'full_'+case['name'], frame, function, helper,
            dict(kind='actual-closed-full-ranker-R2-entry-or-R1-failure', full_case=case['folder'], full_capture_sha256=FULL_CAPTURE),
            expected=expected, expected_exception=error))
    for name, frame in extras(helper):
        reports.append(export(args.output, 'extra_'+name, frame, function, helper, dict(kind='verbatim-R1-statements-only')))
    write(args.output / 'oracle.json', dict(native_environment=native, created_utc=helper.utc(), source_commit=helper.COMMIT,
        full_oracle=str(args.full_oracle), full_source_pins_sha256=FULL_PINS, full_capture_sha256=FULL_CAPTURE,
        harness_sha256=sha(Path(__file__)), source=source, cases=reports, scientific_acceptance=False, production_backend_selected=False,
        scope='R1 statements only; no PED/C1/SciPy/scientific numeric operation executed'))
    write(args.output / 'manifest.json', {str(path.relative_to(args.output)): sha(path) for path in sorted(args.output.rglob('*')) if path.is_file()})
    print(json.dumps(dict(cases=len(reports), values=sum(x['original_exception'] is None for x in reports),
        source_failures=sum(x['original_exception'] is not None for x in reports), oracle_sha256=sha(args.output/'oracle.json'),
        manifest_sha256=sha(args.output/'manifest.json'), native_environment=native)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest='command', required=True).add_parser('capture')
    command.add_argument('--full-oracle', type=Path, required=True); command.add_argument('--output', type=Path, required=True)
    capture(parser.parse_args())


if __name__ == '__main__':
    main()
