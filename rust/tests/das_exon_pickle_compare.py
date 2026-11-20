#!/usr/bin/env python3
"""Compare text-set opcode reconstruction with the pinned Python unpickler.

compare --probe BINARY --raw-pickle CHECKPOINT --work-root NEW_DIRECTORY
Run in a process with an explicit PYTHONHASHSEED. The checkpoint remains
unchanged; this reads its modern transcript exon sets, not a whole-model cache.
"""
import argparse
import io
import itertools
import json
import os
from pathlib import Path
import pickle
import random
import subprocess
import sys


class SetTrace(pickle._Unpickler):
    """Observe ADDITEMS order while keeping exact builtin sets and their layout.

    CPython 3.12's C load_additems passes a tuple to _PySet_Update. Insertion
    from that tuple is distinct from the optimized set-to-set update path.
    All remaining opcodes use the installed Python standard library unchanged.
    """
    dispatch = pickle._Unpickler.dispatch.copy()

    def __init__(self, stream):
        super().__init__(stream)
        self.histories = {}

    def load_empty_set(self):
        value = set()
        self.histories[id(value)] = (value, [])
        self.append(value)

    dispatch[pickle.EMPTY_SET[0]] = load_empty_set

    def load_additems(self):
        items = self.pop_mark()
        value = self.stack[-1]
        if type(value) is not set or id(value) not in self.histories:
            raise ValueError('expected an EMPTY_SET-owned ADDITEMS target')
        self.histories[id(value)][1].extend(items)
        value.update(tuple(items))

    dispatch[pickle.ADDITEMS[0]] = load_additems


def exon_sets(value, path=()):
    if type(value) is not dict:
        return
    if 'exons' in value and type(value['exons']) is set:
        yield path, value['exons']
    for key, child in value.items():
        if type(child) is dict:
            yield from exon_sets(child, path + (key,))


def compare(args):
    seed = os.environ.get('PYTHONHASHSEED')
    if seed is None or not seed.isdigit() or not 0 <= int(seed) <= 2**32 - 1:
        raise RuntimeError('explicit PYTHONHASHSEED is required')
    if sys.version_info[:3] != (3, 12, 3):
        raise RuntimeError('this comparison binds CPython 3.12.3')
    root = args.work_root
    root.mkdir(parents=True, exist_ok=False)
    with args.raw_pickle.open('rb') as stream:
        original = pickle.load(stream)
    originals = list(exon_sets(original))
    with args.raw_pickle.open('rb') as stream:
        trace = SetTrace(stream)
        traced = trace.load()
    traced_sets = list(exon_sets(traced))
    rows = []
    failures = []
    for ordinal, (left, right) in enumerate(itertools.zip_longest(originals, traced_sets)):
        if left is None or right is None or left[0] != right[0]:
            failures.append(dict(stage='trace_path', ordinal=ordinal))
            continue
        path, expected = left
        actual = right[1]
        if expected != actual or list(expected) != list(actual):
            failures.append(dict(stage='trace_set', path=path, python=list(expected), traced=list(actual)))
        history = trace.histories[id(actual)][1]
        if any(type(item) is not str for item in history):
            raise TypeError('the exon reconstruction profile requires exact strings')
        rows.append(dict(stage='checkpoint_unpickle', path=path, items=history, expected=list(expected)))
    checkpoint_count = len(rows)
    multiple = sum(len(value) > 1 for _, value in originals)
    # A second serialization is an observable extra reconstruction, not an
    # identity operation or a claim about an unknown historical process seed.
    roundtrip = pickle.loads(pickle.dumps([value for _, value in originals], protocol=5))
    changed = 0
    for (path, first), second in zip(originals, roundtrip, strict=True):
        changed += list(first) != list(second)
        rows.append(dict(stage='controlled_second_unpickle', path=path, items=list(first), expected=list(second)))
    del original, originals, trace, traced, traced_sets, roundtrip
    # Independent small fixtures cover empty/Unicode strings, duplicates,
    # table growth and ADDITEMS batches beyond the pickler's 1,000-item limit.
    fixtures = [[], [''], ['1/2', '2/2'], ['é', 'Ā', '😀', '', '\x00', '\n']]
    rng = random.Random(90210)
    for size in [1, 4, 5, 18, 19, 76, 77, 300, 1_000, 1_001, 2_019]:
        words = [f'{index}/3000' for index in range(size)]
        fixtures += [words, list(reversed(words)), rng.sample(words, len(words)) + words[:9]]
    fixture_count = 0
    for case, items in enumerate(fixtures):
        value = set(items)
        for protocol in [4, 5]:
            payload = pickle.dumps(value, protocol=protocol)
            traced_reader = SetTrace(io.BytesIO(payload))
            actual = traced_reader.load()
            expected = pickle.loads(payload)
            if list(actual) != list(expected):
                failures.append(dict(stage='fixture_trace', case=case, protocol=protocol))
            rows.append(dict(stage='fixture_unpickle', case=case, protocol=protocol,
                             items=traced_reader.histories[id(actual)][1], expected=list(expected)))
            fixture_count += 1
    incoming = root / 'sets.input'
    with incoming.open('x') as stream:
        for row in rows:
            stream.write(str(len(row['items'])) + ''.join(' ' + (value.encode().hex() or '-') for value in row['items']) + '\n')
    output = root / 'sets.native'
    command = [str(args.probe), 'batch', '--input', str(incoming), '--output', str(output), '--seed', seed]
    with (root / 'native.stdout.log').open('x') as stdout, (root / 'native.stderr.log').open('x') as stderr:
        status = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=120).returncode
    compared = 0
    with output.open() as stream:
        for ordinal, (row, line) in enumerate(itertools.zip_longest(rows, stream)):
            if row is None or line is None:
                failures.append(dict(stage='native_shape', ordinal=ordinal))
                continue
            words = line.split()
            actual = [bytes.fromhex(word if word != '-' else '').decode() for word in words[1:]]
            compared += 1
            if not words or int(words[0]) != len(actual) or actual != row['expected']:
                failures.append(dict(stage='native_order', ordinal=ordinal, reference=row, native=actual))
    with (root / 'mismatches.jsonl').open('x') as stream:
        for failure in failures:
            stream.write(json.dumps(failure) + '\n')
    # Retain opcode histories separately from public count-only metadata.
    with (root / 'cases.jsonl').open('x') as stream:
        for row in rows:
            stream.write(json.dumps(row) + '\n')
    result = dict(passed=status == 0 and not failures and compared == len(rows), hash_seed=int(seed),
                  checkpoint_sets=checkpoint_count, multiple_exon_sets=multiple, fixture_cases=fixture_count,
                  controlled_second_load_order_changes=changed, compared=compared, mismatches=len(failures),
                  returncode=status, python=sys.version, final_model_history_reconstructed=False,
                  unknown_historical_seed_equivalence=False, canonicalization_adopted=False)
    (root / 'report.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)
    if not result['passed']:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest='command', required=True)
    command = subcommands.add_parser('compare')
    for option in ['probe', 'raw-pickle', 'work-root']:
        command.add_argument('--' + option, type=Path, required=True)
    command.set_defaults(function=compare)
    args = parser.parse_args()
    args.function(args)


if __name__ == '__main__':
    main()
