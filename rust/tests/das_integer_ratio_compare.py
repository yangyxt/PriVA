#!/usr/bin/env python3
"""Compare exact integer division with pinned CPython, including binary64 edges.

check --binary PROBE --work-root NEW_DIRECTORY --report JSON
The oracle evaluates Python's / directly. Every result is compared by bits;
zero-division and overflow exceptions are compared by exact category.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import random
import struct
import subprocess
import sys


def bits(value):
    return struct.pack('>d', value).hex()


def cases():
    result = []
    def add(a, b, name):
        for x, y in [(a, b), (-a, b), (a, -b), (-a, -b)]:
            result.append((x, y, name))
    for a in [0, 1, 2, 3, 17, 2**53-1, 2**53+1, 2**63, 2**64, 10**4299]:
        for b in [0, 1, 2, 3, 7, 2**53-1, 2**53+1, 2**64, 10**4299]:
            add(a, b, 'decimal_and_integer_boundaries')
    for exponent in [-1100, -1076, -1075, -1074, -1073, -1023, -1022, -1021, -54, -53, -52, -1, 0, 1, 52, 53, 54, 1022, 1023, 1024, 1100]:
        if exponent < 0:
            for delta in [-2, -1, 0, 1, 2]:
                add(1, 2**(-exponent)+delta, 'power_denominator')
        else:
            for delta in [-2, -1, 0, 1, 2]:
                add(2**exponent+delta, 1, 'power_numerator')
    # Adjacent real binary64 values define their exact rational midpoint.
    # Denominator perturbations and large common factors exercise each side
    # without calculating expected results through the native algorithm.
    words = [0, 1, 2, 3, 2**52-2, 2**52-1, 2**52, 2**52+1,
             0x3feffffffffffffe, 0x3fefffffffffffff, 0x3ff0000000000000,
             0x3ff0000000000001, 0x4330000000000000,
             0x7feffffffffffffd, 0x7feffffffffffffe]
    rng = random.Random(0xD1A1DE)
    words += [rng.randrange(0x7fefffffffffffff) for _ in range(160)]
    for word in words:
        left = struct.unpack('>d', word.to_bytes(8, 'big'))[0]
        right = struct.unpack('>d', (word+1).to_bytes(8, 'big'))[0]
        midpoint = (Fraction(left)+Fraction(right))/2
        for scale in [1, 2**73+17, 10**330+19]:
            a, b = midpoint.numerator*scale, midpoint.denominator*scale
            for delta in [-1, 0, 1]:
                add(a+delta, b, 'adjacent_float_midpoint')
    overflow = (2**54-1)*2**970
    for delta in [-2, -1, 0, 1, 2]:
        add(overflow+delta, 1, 'overflow_rounding')
    for _ in range(5000):
        sizes = [0, 1, 2, 31, 53, 64, 65, 127, 512, 1023, 1074, 1138, 4096, 14000]
        a = rng.getrandbits(rng.choice(sizes))
        b = rng.getrandbits(rng.choice(sizes))
        add(a, b, 'deterministic_random')
    return result


def identity(path):
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.file_digest(path.open('rb'), 'sha256').hexdigest())


def check(args):
    if sys.version_info[:3] != (3, 12, 3):
        raise RuntimeError('this oracle is pinned to CPython 3.12.3')
    root = args.work_root
    root.mkdir(parents=True, exist_ok=False)
    pins = [identity(Path(__file__)), identity(args.binary)]
    records = cases()
    source = root / 'pairs.tsv'
    with source.open('x') as output:
        for a, b, _ in records:
            output.write(f'{a}\t{b}\n')
    with source.open() as input_stream, (root / 'native.txt').open('x') as output, (root / 'native.stderr.log').open('x') as error:
        status = subprocess.run([str(args.binary), 'divide'], stdin=input_stream,
                                stdout=output, stderr=error, timeout=180).returncode
    actual = (root / 'native.txt').read_text().splitlines()
    failures = []
    exceptions = {}
    counts = {}
    with (root / 'python.txt').open('x') as output:
        for ordinal, (a, b, name) in enumerate(records):
            counts[name] = counts.get(name, 0)+1
            try:
                expected = 'ok '+bits(a/b)
            except (OverflowError, ZeroDivisionError) as error:
                category = type(error).__name__
                exceptions[category] = exceptions.get(category, 0)+1
                expected = 'err '+category
            output.write(expected+'\n')
            observed = actual[ordinal] if ordinal < len(actual) else None
            if observed != expected:
                failures.append(dict(ordinal=ordinal, case=name, numerator=str(a), denominator=str(b),
                                     python=expected, native=observed))
    if len(actual) != len(records):
        failures.append(dict(case='result_count', expected=len(records), actual=len(actual)))
    (root / 'mismatches.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in failures))
    for pin in pins:
        if identity(Path(pin['path'])) != pin:
            raise RuntimeError('comparison source/binary changed')
    result = dict(passed=status == 0 and not failures, cases=len(records), mismatches=len(failures),
                  native_exit=status, exceptions=exceptions, groups=counts, inputs=pins,
                  python_version=sys.version, comparison='all binary64 bits or exact exception category')
    (root / 'report.json').write_text(json.dumps(result, indent=2)+'\n')
    args.report.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ['passed', 'cases', 'mismatches', 'exceptions']}), flush=True)
    if not result['passed']:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    command = commands.add_parser('check')
    for name in ['binary', 'work-root', 'report']:
        command.add_argument('--'+name, type=Path, required=True)
    command.set_defaults(function=check)
    args = parser.parse_args()
    args.function(args)


if __name__ == '__main__':
    main()
