#!/usr/bin/env python3
"""Qualify the narrow exon metadata projection against the C/Python unpicklers.

compare --probe BINARY --raw-pickle CHECKPOINT --collector-source PYTHON
        --work-root NEW_DIRECTORY
Numerical payloads are tested for safe skipping, not decoded by this reader.
"""
import argparse
import ast
import io
import itertools
import json
import os
from pathlib import Path
import pickle
import random
import subprocess
import sys
import numpy as np
from das_exon_pickle_compare import SetTrace


def all_sets(value, path=()):
    if type(value) is dict:
        for key, item in value.items():
            if type(key) is str:
                if type(item) is set:
                    yield path + (key,), item
                elif type(item) is dict:
                    yield from all_sets(item, path + (key,))


def decode(line):
    words = line.split()
    count = int(words[0])
    text = lambda word: bytes.fromhex(word if word != '-' else '').decode()
    path = tuple(map(text, words[1:count+1]))
    size = int(words[count+1])
    items = list(map(text, words[count+2:]))
    if len(items) != size:
        raise ValueError('native item count differs')
    return path, items


def check(path, probe, root, label, failures):
    with path.open('rb') as stream:
        expected = pickle.load(stream)
    with path.open('rb') as stream:
        trace = SetTrace(stream)
        model = trace.load()
    c_sets = list(all_sets(expected))
    py_sets = list(all_sets(model))
    expected_rows = []
    for ordinal, (left, right) in enumerate(itertools.zip_longest(c_sets, py_sets)):
        if left is None or right is None or left[0] != right[0] or list(left[1]) != list(right[1]):
            failures.append(dict(stage='unpickler_order', case=label, ordinal=ordinal))
            continue
        expected_rows.append((left[0], trace.histories[id(right[1])][1]))
    output = root / (label + '.native')
    command = [str(probe), 'extract', '--input', str(path), '--output', str(output)]
    with (root / (label + '.stdout.log')).open('x') as out, (root / (label + '.stderr.log')).open('x') as err:
        status = subprocess.run(command, stdout=out, stderr=err, timeout=180).returncode
    if status:
        failures.append(dict(stage='native_exit', case=label, returncode=status))
        return dict(case=label, returncode=status, compared=0)
    actual_rows = [decode(line) for line in output.read_text().splitlines()]
    for ordinal, (wanted, actual) in enumerate(itertools.zip_longest(expected_rows, actual_rows)):
        if wanted != actual:
            failures.append(dict(stage='opcode_history', case=label, ordinal=ordinal, expected=wanted, actual=actual))
    statistics = json.loads((root / (label + '.stdout.log')).read_text())
    return dict(case=label, returncode=status, compared=len(expected_rows), statistics=statistics)


def compare(args):
    if sys.version_info[:3] != (3, 12, 3) or os.environ.get('PYTHONHASHSEED') != '0':
        raise RuntimeError('comparison requires CPython 3.12.3 / hash seed 0')
    root = args.work_root
    root.mkdir(parents=True, exist_ok=False)
    failures = []
    results = []
    shared = {'exons': {'1/2', '2/2'}, 'aa_min': 1, 'aa_max': 50}
    objects = [
        {},
        {'gene': {'db': {'dom': {'ENST1': shared, 'ENST2': shared, 'legacy': set()}}}},
        {'unicode': {'ENSTé': {'exons': {'', 'é', 'Ā', '😀', '\x00', '\n'}}}},
        {'large': {'exons': {f'{i}/3000' for i in range(2019)}}},
        {'ignored': [None, True, -2**200, float('nan'), b'bytes', bytearray(b'other')], 'exons': {'1'}},
    ]
    for dtype in ['f8', 'f4', 'i8', 'u8', 'b1', '>f8']:
        array = np.arange(10000).astype(dtype)
        objects.append({'array': array, 'copy': array, 'scalar': array[0], 'exons': {'3/4', '1/4', '2/4'}})
    for ordinal, value in enumerate(objects):
        for protocol in [4, 5]:
            label = f'fixture_{ordinal:03d}_p{protocol}'
            path = root / (label + '.pkl')
            with path.open('xb') as stream:
                pickle.dump(value, stream, protocol=protocol)
            results.append(check(path, args.probe, root, label, failures))
    # Plain dict serialization normally writes each key only once. Exercise
    # the stack machine's first-position/last-value rules independently, with
    # numeric overwrites and both SETITEM and SETITEMS forms.
    rng = random.Random(303)
    def text(value):
        raw = value.encode()
        return b'\x8c' + bytes([len(raw)]) + raw
    for case in range(100):
        payload = b'\x80\x05}'
        batch = case % 2 == 0
        if batch: payload += b'('
        for _ in range(rng.randrange(1, 25)):
            payload += text('key' + str(rng.randrange(8)))
            if rng.randrange(5) == 0:
                payload += b'N'
            else:
                payload += b'\x8f('
                for _ in range(rng.randrange(12)):
                    payload += text(rng.choice(['', 'é', 'Ā', '😀', '1/2', '2/2']))
                payload += b'\x90'
            if not batch: payload += b's'
        payload += (b'u' if batch else b'') + b'.'
        label = f'assignments_{case:03d}'
        path = root / (label + '.pkl'); path.write_bytes(payload)
        results.append(check(path, args.probe, root, label, failures))
    results.append(check(args.raw_pickle, args.probe, root, 'raw_checkpoint', failures))
    # The source finalizer converts score dicts into the actual NumPy arrays
    # emitted by the retained final model, without recollecting any VCF rows.
    body = [node for node in ast.parse(args.collector_source.read_text()).body
            if isinstance(node, ast.FunctionDef) and node.name == '_finalize_scores_dict']
    if len(body) != 1:
        raise RuntimeError('source finalizer absent')
    scope = dict(np=np)
    exec(compile(ast.fix_missing_locations(ast.Module(body=body, type_ignores=[])), str(args.collector_source), 'exec'), scope)
    with args.raw_pickle.open('rb') as stream:
        model = pickle.load(stream)
    scope['_finalize_scores_dict'](model)
    finalized = root / 'finalized_checkpoint.pkl'
    with finalized.open('xb') as stream:
        pickle.dump(model, stream, protocol=5)
    del model
    results.append(check(finalized, args.probe, root, 'finalized_checkpoint', failures))
    # These are explicitly outside this reader's projection profile. Ensure a
    # refusal is named and produces no partial output, without calling Python.
    refusals = {
        'protocol_three': pickle.dumps({'exons': {'1'}}, protocol=3),
        'nontext_set': pickle.dumps({'exons': {1}}, protocol=5),
        'nontext_metadata_key': pickle.dumps({1: {'exons': {'1'}}}, protocol=5),
        'trailing_bytes': pickle.dumps({}, protocol=5) + b'x',
        'truncated': pickle.dumps({'exons': {'1'}}, protocol=5)[:-1],
        'unset_memo': b'\x80\x05}h\x05.',
        'unknown_global': b'\x80\x05}\x8c\x08builtins\x8c\x04eval\x93.',
    }
    cyclic = {}; cyclic['cycle'] = cyclic
    refusals['cyclic_dictionary'] = pickle.dumps(cyclic, protocol=5)
    for label, payload in refusals.items():
        path = root / (label + '.pkl'); path.write_bytes(payload)
        output = root / (label + '.native')
        command = [str(args.probe), 'extract', '--input', str(path), '--output', str(output)]
        with (root / (label + '.stdout.log')).open('x') as out, (root / (label + '.stderr.log')).open('x') as err:
            status = subprocess.run(command, stdout=out, stderr=err, timeout=10).returncode
        if status != 1 or output.exists() or 'pickle byte' not in (root / (label + '.stderr.log')).read_text():
            failures.append(dict(stage='projection_refusal', case=label, returncode=status, output_exists=output.exists()))
    (root / 'mismatches.jsonl').write_text(''.join(json.dumps(row) + '\n' for row in failures))
    result = dict(passed=not failures, comparisons=sum(row['compared'] for row in results),
                  cases=len(results), refusals=len(refusals), mismatches=len(failures), results=results,
                  python=sys.version, numpy=np.__version__, numerical_values_decoded=False,
                  full_final_model_examined=False, full_das_qualified=False)
    (root / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key:result[key] for key in ['passed','comparisons','cases','refusals','mismatches']}), flush=True)
    if not result['passed']:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest='command', required=True)
    command = subcommands.add_parser('compare')
    for name in ['probe', 'raw-pickle', 'collector-source', 'work-root']:
        command.add_argument('--' + name, type=Path, required=True)
    command.set_defaults(function=compare)
    args = parser.parse_args()
    args.function(args)


if __name__ == '__main__':
    main()
