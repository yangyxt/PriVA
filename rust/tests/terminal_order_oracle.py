#!/usr/bin/env python3
"""Capture actual R8-entry/final full-ranker evidence before the opt-in Rust port.

capture --full-oracle CLOSED_ROOT --output NEW_DATA_ROOT
The closed full-ranking oracle supplies original source checkpoints. Extra cases
execute its unchanged complete ranker; a verbatim R8 AST excerpt checks the
boundary and captures R8-only failure/mutation behavior. No scientific edits.
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
import sys

import numpy as np
import pandas as pd

Y = Path('/stornext/snfs170/next-gen/scratch/yangyxt')
D = Path('/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade')
FULL_PINS = '9e63fb73e49c9f8ce700b0abd3e588702d2489c5430059e76c5d5b723eabb46a'
FULL_CAPTURE = 'db919b02f8738da60f85d16e0a7ed2bb53f8eaaf6047ed8f4db912fe4938508c'
FULL_HARNESS = '7173862933e0b92201707231275ff9e041fa8d49ccf60e33c47c841ef1f5b839'
KEYS = ['chrom', 'pos', 'ref', 'alt', 'sort_index', 'ACMG_quant_score', 'variant_id']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True, allow_nan=False)
        stream.write('\n')


def cell(value):
    """Explicit primitives only. Index labels instead travel as opaque JSON tokens."""
    if value is None or value is pd.NA:
        return ['null', '']
    if isinstance(value, (bool, np.bool_)):
        return ['bool', str(int(value))]
    if isinstance(value, (int, np.integer)):
        return ['uint' if value > 2**63 - 1 else 'int', str(value)]
    if isinstance(value, (float, np.floating)):
        return ['float', struct.pack('>d', value).hex()]
    if isinstance(value, str):
        return ['str', value.encode().hex()]
    raise TypeError('unsupported primitive transport: ' + str(type(value)))


def atom(value):
    kind = value['kind']
    if kind in ['none', 'pd.NA']:
        return None if kind == 'none' else pd.NA
    if kind == 'str':
        return value['value']
    if kind == 'int':
        return int(value['decimal'])
    if kind == 'bool':
        return value['value']
    if kind == 'float':
        return struct.unpack('>d', bytes.fromhex(value['ieee_f64_be']))[0]
    if kind == 'numpy_scalar':
        return np.frombuffer(bytes.fromhex(value['storage_hex']), dtype=value['dtype'])[0]
    raise TypeError('unsupported source scalar: ' + kind)


def encoded_frame(value):
    """Restore only declared primitive columns; retain the original encoding too."""
    columns = [atom(x) for x in value['columns']['values']]
    frame = pd.concat([pd.Series([atom(x) for x in item['values']], dtype=item['dtype']['text'])
                       for item in value['data']], axis=1)
    frame.columns = columns
    index = value['index']
    frame.index = pd.Index([atom(x) for x in index['values']], dtype=index['dtype']['text'],
                          name=atom(index['names']['values'][0]))
    return frame


def r8_function(source):
    text = source.read_text()
    node = next(x for x in ast.parse(text).body if isinstance(x, ast.FunctionDef)
                and x.name == 'sort_and_rank_variants')
    body = [copy.deepcopy(x) for x in node.body if x.lineno >= 668]
    assert body[0].lineno == 668 and body[-1].lineno == 697
    function = ast.FunctionDef(name='original_R8_excerpt',
        args=ast.arguments(posonlyargs=[], args=[ast.arg(arg='df')], kwonlyargs=[], kw_defaults=[], defaults=[]),
        body=body, decorator_list=[], lineno=668, col_offset=0, end_lineno=697, end_col_offset=20)
    module = ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[]))
    namespace = {'pd': pd}
    exec(compile(module, str(source), 'exec'), namespace)
    return namespace['original_R8_excerpt'], dict(lines=[668, 697],
        statement_ast_sha256=hashlib.sha256(ast.dump(ast.Module(body=body, type_ignores=[])).encode()).hexdigest())


def extras(helper):
    def base(n, **changes):
        return helper.base_frame(n, **changes)
    payload_nan = struct.unpack('>d', bytes.fromhex('7ff8000000000123'))[0]
    yield 'signed_zero_equal_group', base(6, pos=[1, 1, 2, 2, 3, 3],
        ACMG_quant_score=[-0., 0., 0., -0., payload_nan, np.nan])
    yield 'nonfinite_group_reductions', base(8, pos=[1, 1, 2, 2, 3, 3, 4, 4],
        ACMG_quant_score=[np.nan, np.inf, -np.inf, np.nan, payload_nan, np.nan, -np.inf, np.inf])
    yield 'null_all_coordinate_keys', base(5, chrom=[None, '1', '1', '1', '1'],
        pos=[1., np.nan, 1., 1., 1.], ref=['A', 'A', None, 'A', 'A'], alt=['G', 'G', 'G', None, 'G'])
    yield 'all_six_key_ties', base(24, pos=[1]*24, variant_id=['b', 'a', None, 'b', 'c', np.nan]*4)
    yield 'ids_object_missing', base(6, variant_id=['same', None, np.nan, pd.NA, 'same', ''])
    for name, values, dtype in [
        ('float', [1., np.nan, -0., 0., np.inf, -np.inf], 'float64'),
        ('int', [1, 7, 1, -1, 0, -1], 'int64'),
        ('uint', [2**64-1, 2**63, 1, 2**64-1, 0, 1], 'uint64'),
        ('bool', [True, False, True, False, True, False], 'bool'),
        ('nullable_string', ['x', None, 'x', '', None, 'y'], 'string')]:
        frame = base(6); frame['variant_id'] = pd.Series(values, dtype=dtype)
        yield 'ids_' + name, frame
    frame = base(6, pos=[0., -0., np.nan, 1.e-5, 1.e16, np.inf]).drop(columns='variant_id')
    yield 'absent_id_float_coordinate', frame
    frame = base(3, chrom=['a:b', 'a', 'a'], pos=['c', 'b:c', 'z']).drop(columns='variant_id')
    yield 'absent_id_colon_collision', frame
    for name, values, dtype in [('string', ['z', None, 'a', 'z'], 'string'),
                                ('Int64', [1, None, -1, 1], 'Int64'),
                                ('UInt64', [2**64-1, None, 0, 2**64-1], 'UInt64')]:
        frame = base(4); frame['pos'] = pd.Series(values, dtype=dtype)
        yield 'coordinate_' + name, frame
    for name, index in [('object', pd.Index(['x', None, 'x', ''], dtype=object, name='opaque')),
                        ('multi', pd.MultiIndex.from_tuples([('b', 2), ('a', 9), ('b', 2), ('z', -3)], names=['a', 'b'])),
                        ('range', pd.RangeIndex(9, 1, -2, name='reverse')),
                        ('float', pd.Index([np.nan, -0., 0., np.nan], name='float-index'))]:
        frame = base(4, ACMG_quant_score=[.7, .9, .9, .8]); frame.index = index
        yield 'index_' + name, frame
    frame = base(3)
    frame.insert(0, 'variant_rank', ['stale']*3)
    frame.insert(2, 'max_variant_score', ['stale']*3)
    frame.insert(1, 'sort_index', ['stale']*3)
    frame.columns.name = 'kept-column-axis'; frame.attrs['context'] = 'opaque attrs'
    yield 'existing_column_positions', frame
    yield 'empty', base(0)


def export_case(output, name, frame, expected, helper, excerpt, origin):
    folder = output / name
    folder.mkdir()
    prepared = pickle.dumps(frame, protocol=5)
    (folder / 'entry.pkl').write_bytes(prepared)
    write(folder / 'entry.json', helper.encode(frame))
    working = pickle.loads(prepared)
    maxima = []
    def trace(pyframe, event, arg):
        if pyframe.f_code is excerpt.__code__ and event == 'line' and pyframe.f_lineno == 672:
            maxima.append(helper.encode(pyframe.f_locals['max_scores']))
        return trace
    previous = sys.gettrace()
    try:
        sys.settrace(trace)
        result = excerpt(working)
    except Exception as error:
        exception = dict(type=type(error).__name__, message=str(error))
        result = None
    else:
        exception = None
    finally:
        sys.settrace(previous)
    write(folder / 'entry-after.json', helper.encode(working))
    write(folder / 'maxima.json', maxima)
    report = dict(name=name, origin=origin, columns=list(frame.columns), dtypes=[str(x) for x in frame.dtypes],
                  rows=len(frame), original_exception=exception, entry_sha256=hashlib.sha256(prepared).hexdigest(),
                  full_ranker_agrees=None, index_opaque=True)
    with (folder / 'input.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value'] + [x for name in frame.columns for x in [name+'_kind', name]])
        for row in range(len(frame)):
            token = helper.canonical(helper.encode(frame.index[row])).hex()
            writer.writerow(['opaque', token] + [v for col in range(frame.shape[1]) for v in cell(frame.iloc[row, col])])
    if result is not None:
        (folder / 'expected.pkl').write_bytes(pickle.dumps(result, protocol=5))
        write(folder / 'expected.json', helper.encode(result))
        # Row identity is the original fixture token, independent of duplicate index labels.
        positions = {token: i for i, token in enumerate(frame['__oracle_row_position'])}
        order = [positions[token] for token in result['__oracle_row_position']]
        report.update(permutation=order, ranks=[int(x) for x in result['variant_rank']],
                      output_columns=list(result.columns), output_dtypes=[str(x) for x in result.dtypes])
        if expected is not None:
            pd.testing.assert_frame_equal(result, expected, check_exact=True)
            assert helper.encode(result) == helper.encode(expected), name
            report['full_ranker_agrees'] = True
    return report


def capture(args):
    assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_STEP_ID')
    assert args.full_oracle.resolve().is_relative_to(D) and args.output.resolve().is_relative_to(D)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'capture-harness.py').write_bytes(Path(__file__).read_bytes())
    write(args.output / 'attempt.json', dict(native_environment={name: os.environ.get(name) for name in
        ['SLURM_JOB_ID', 'SLURM_STEP_ID', 'SLURMD_NODENAME']}, source_sha256=sha(Path(__file__))))
    assert sha(args.full_oracle / 'source-pins.json') == FULL_PINS
    assert sha(args.full_oracle / 'capture.json') == FULL_CAPTURE
    harness = args.full_oracle / 'source/ranking_control_oracle.py'
    assert sha(harness) == FULL_HARNESS
    spec = importlib.util.spec_from_file_location('frozen_full_oracle', harness)
    helper = importlib.util.module_from_spec(spec); spec.loader.exec_module(helper)
    helper.verify_sources(args.full_oracle, FULL_PINS)
    scoring, _, _, _ = helper.original_modules(args.full_oracle / 'source/scripts')
    excerpt, details = r8_function(args.full_oracle / 'source/scripts/acmg_scoring.py')
    write(args.output / 'runtime-identities.json', helper.runtime_identity())
    reports = []
    full = json.loads((args.full_oracle / 'capture.json').read_text())
    for item in full['cases']:
        if item['kind'] != 'ranking' or item['status'] != 'value':
            continue
        source = args.full_oracle / item['folder']
        trace = json.loads((source / 'traced/trace.json').read_text())
        # The original observer deduplicates identical states. The before-672
        # checkpoint contains the untouched R8-entry df plus computed maxima;
        # groupby/transform have not yet written max_variant_score into df.
        events = [event for event in trace['events'] if event['function'] == 'sort_and_rank_variants'
                  and event['event'] == 'line' and event['line'] == 672]
        assert len(events) == 1, item['name']
        state = json.loads((args.full_oracle / 'states' / (events[0]['state'] + '.json')).read_text())
        frame = encoded_frame(state['df'])
        expected = pickle.loads((source / 'traced/output.pkl').read_bytes())
        # Recover axes/attrs from the actual output without changing source column/row data.
        frame.columns.name = expected.columns.name
        frame.attrs = copy.deepcopy(expected.attrs)
        reports.append(export_case(args.output, 'full_' + item['name'], frame, expected, helper, excerpt,
            dict(kind='closed-full-ranker-checkpoint', folder=item['folder'], entry_state=events[0]['state'],
                 full_capture_sha256=FULL_CAPTURE)))
    for name, frame in extras(helper):
        entries = []
        def observe(pyframe, event, arg):
            if pyframe.f_code is scoring.sort_and_rank_variants.__code__ and event == 'line' and pyframe.f_lineno == 668:
                entries.append(pickle.dumps(pyframe.f_locals['df'], protocol=5))
            return observe
        before = helper.encode(frame)
        previous = sys.gettrace()
        try:
            sys.settrace(observe)
            expected = scoring.sort_and_rank_variants(frame, None, None, dispensable_gene_list=None)
        finally:
            sys.settrace(previous)
        assert before == helper.encode(frame), name
        assert len(entries) == 1
        entry = pickle.loads(entries[0])
        # No AST substitution in these full-ranker captures; R8 excerpt independently checks the boundary.
        reports.append(export_case(args.output, 'extra_' + name, entry, expected, helper, excerpt,
                                    dict(kind='unchanged-full-ranker-additional-case', original_input_unchanged=True)))
    base = helper.base_frame(2); base['sort_index'] = [.9, .8]
    for column in ['chrom', 'pos', 'ref', 'alt', 'sort_index', 'ACMG_quant_score']:
        reports.append(export_case(args.output, 'failure_missing_' + column, base.drop(columns=column), None,
                       helper, excerpt, dict(kind='verbatim-R8-excerpt-failure-only')))
    write(args.output / 'oracle.json', dict(native_environment=helper.native_environment(), created_utc=helper.utc(),
        source_commit=helper.COMMIT, full_oracle=str(args.full_oracle), source_pins_sha256=FULL_PINS,
        full_capture_sha256=FULL_CAPTURE, harness_sha256=sha(Path(__file__)), excerpt=details, cases=reports,
        scope='R8 checkpoint boundary only; unchanged full ranker for successful cases; explicit source R8 excerpt failures',
        scientific_acceptance=False, production_backend_selected=False))
    write(args.output / 'manifest.json', {str(path.relative_to(args.output)): sha(path) for path in sorted(args.output.rglob('*')) if path.is_file()})
    print(json.dumps(dict(cases=len(reports), failures=sum(x['original_exception'] is not None for x in reports),
                          oracle_sha256=sha(args.output / 'oracle.json'), manifest_sha256=sha(args.output / 'manifest.json'),
                          native_environment=helper.native_environment())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest='command', required=True)
    capture_parser = command.add_parser('capture')
    capture_parser.add_argument('--full-oracle', type=Path, required=True)
    capture_parser.add_argument('--output', type=Path, required=True)
    capture(parser.parse_args())


if __name__ == '__main__':
    sys.dont_write_bytecode = True
    main()
