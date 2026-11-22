#!/usr/bin/env python3
"""Compare exact NumPy seeds, sampled indices and every resulting PCG64 state.

python numpy_rng_compare.py compare --probe BIN --work-root NEW_DIR
Uses only installed NumPy 1.26.4 as the independent oracle. Covers the choice
algorithm boundary, rejection sampling, cached uint32/uint64 interleaving,
zero-sized and impossible draws, repeated sampling and chained operations.
"""
import argparse
import hashlib
import importlib
import itertools
import json
from pathlib import Path
import random
import subprocess

import numpy as np


def cases():
    generator = random.Random(40780228)
    entropy = [[], [0], [42], [0, 42], [42, 0], [2**32-1]*4,
               [3735928559, 195939070, 229505742, 305419896]]
    entropy += [[generator.getrandbits(32) for _ in range(length)] for length in [1, 2, 3, 4, 5, 8, 16, 100]]
    for raw in entropy:
        for count in [0, 1, 4, 8, 17, 129]:
            yield ['words', count, *raw]
    seeds = [0, 1, 42, 99, 2**32-1, 2**32, 2**64-1]
    for seed in seeds:
        yield ['seed', seed]
        yield ['raw', 100]
        yield ['uint32', 1]
        yield ['raw', 3]
        yield ['uint32', 3]
        for maximum in [0, 1, 2, 255, 256, 2**31, 2**31+1, 2**32-2, 2**32-1,
                        2**32, 2**63, 2**63+1, 2**64-2, 2**64-1]:
            yield ['bounded', maximum, 257]
        for population in [0, 1, 2, 7, 15, 16, 17, 9999, 10000, 10001, 28973, 200000]:
            sizes = sorted(set([0, 1, max(0, population//50-1), population//50,
                                population//50+1, population//2, population, population+1]))
            for size in sizes:
                for replace in [0, 1]:
                    # Keep generated fixtures bounded; whole 200000 populations
                    # are covered once while the other seeds use the boundary.
                    if population==200000 and size>5000 and seed!=42:
                        continue
                    yield ['choice', population, size, replace]
                yield ['raw', 1]
        for population in [2**32-1, 2**32, 2**32+1, 2**63-1]:
            for replace in [0, 1]:
                yield ['choice', population, 19, replace]
        yield ['uint32', 1]
        yield ['choice', 0, 1, 0]
        yield ['choice', 2, 3, 0]
        yield ['uint32', 1]
    for _ in range(250):
        yield ['seed', generator.getrandbits(64)]
        for _ in range(20):
            population = generator.choice([1, 13, 9999, 10000, 10001, 28973, 200000, 2**32+1])
            size = generator.randrange(0, min(population, 1000)+1)
            yield ['choice', population, size, generator.randrange(2)]
            yield ['uint32', generator.randrange(4)]
            yield ['raw', generator.randrange(4)]


def state(rng):
    value = rng.bit_generator.state
    return f"{value['state']['state']:032x}:{value['state']['inc']:032x}:{value['has_uint32']}:{value['uinteger']:08x}"


def compare(args):
    if np.__version__ != '1.26.4':
        raise RuntimeError('NumPy reference version differs')
    root = args.work_root.absolute()
    root.mkdir(parents=True, exist_ok=False)
    count = values = failures = words = states = 0
    rng = np.random.default_rng(0)
    with (root/'input.txt').open('x') as incoming, (root/'python.txt').open('x') as truth:
        for command in cases():
            count += 1
            incoming.write(' '.join(map(str, command))+'\n')
            name, *arguments = command
            if name=='words':
                actual = np.random.SeedSequence(arguments[1:]).generate_state(arguments[0])
                truth.write('words'+''.join(f' {int(value):016x}' for value in actual)+'\n')
                words += 1
                values += len(actual)
                continue
            try:
                if name=='seed':
                    rng = np.random.default_rng(arguments[0])
                    actual = []
                elif name=='raw':
                    actual = rng.bit_generator.random_raw(arguments[0])
                elif name=='uint32':
                    actual = rng.integers(0, 2**32, size=arguments[0], dtype=np.uint32)
                elif name=='bounded':
                    actual = rng.integers(0, arguments[0], endpoint=True, size=arguments[1], dtype=np.uint64)
                else:
                    actual = rng.choice(arguments[0], size=arguments[1], replace=bool(arguments[2]))
                values += len(actual)
                result = 'ok'+''.join(f' {int(value):016x}' for value in actual)
            except Exception as error:
                failures += 1
                result = type(error).__name__+' '+str(error)
            states += 1
            truth.write(result+'|'+state(rng)+'\n')
    with (root/'input.txt').open('rb') as incoming, (root/'native.txt').open('xb') as output, (root/'native.stderr.log').open('xb') as errors:
        status = subprocess.run([str(args.probe), 'compare', '--profile', 'numpy-1.26.4'], stdin=incoming, stdout=output, stderr=errors).returncode
    mismatches = compared = 0
    with (root/'input.txt').open() as inputs, (root/'python.txt').open() as truth, (root/'native.txt').open() as actual, (root/'mismatches.jsonl').open('x') as ledger:
        for line, (command, expected, observed) in enumerate(itertools.zip_longest(inputs, truth, actual), 1):
            compared += 1
            if expected != observed:
                mismatches += 1
                ledger.write(json.dumps(dict(line=line, command=command, python=expected, native=observed))+'\n')
    modules = {}
    for name in ['numpy.random.bit_generator', 'numpy.random._pcg64', 'numpy.random._generator', 'numpy.random._bounded_integers']:
        path = Path(importlib.import_module(name).__file__)
        modules[name] = dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    result = dict(passed=not mismatches and status==0 and compared==count, cases=count, sampled_values=values,
        seed_word_cases=words, generator_states=states, reference_failures=failures, mismatches=mismatches,
        native_returncode=status, numpy_version=np.__version__, modules=modules,
        probe_sha256=hashlib.sha256(args.probe.read_bytes()).hexdigest(),
        scope='unspawned pool4 SeedSequence, PCG64 and unweighted 1D choice with shuffle=True; full DAS remains pending')
    (root/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
    if not result['passed']:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('compare')
    p.add_argument('--probe', type=Path, required=True)
    p.add_argument('--work-root', type=Path, required=True)
    args = parser.parse_args()
    compare(args)


if __name__ == '__main__':
    main()
