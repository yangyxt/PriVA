#!/usr/bin/env python3
"""Exact DAS BH outputs against pinned statsmodels; no caller activation.

python statsmodels_fdr_compare.py compare --probe BIN --work-root NEW_DIR
Retains every input, expected/native line and mismatch. The adjusted IEEE-754
bits, rejection masks and defined empty-input exception are compared.
"""
import argparse
import hashlib
import inspect
import itertools
import json
import math
from pathlib import Path
import struct
import subprocess
import warnings

import numpy as np
import statsmodels
from statsmodels.stats.multitest import multipletests
from fdr_compare import bits, cases, number


def compare(args):
    root = args.work_root.absolute()
    root.mkdir(parents=True, exist_ok=False)
    alphas = [.01, .05, 0., 1., -.5, 1.5, math.nan, math.inf, -math.inf]
    count = failures = values = 0

    def extended():
        for index, (name, raw) in enumerate(cases(args.random_cases)):
            yield name, alphas[index % len(alphas)], raw
        payloads = [number(0x7ff8000000000001), number(0xfff8000000000042), math.nan]
        for length in [2, 3, 15, 16, 17, 31, 32, 33, 127, 128, 129, 4096]:
            for offset in range(3):
                yield f'nan-payloads-{length}-{offset}', .01, [payloads[(i+offset)%3] for i in range(length)]
            for alpha in alphas:
                yield f'alpha-{length}-{bits(alpha)}', alpha, [.01 if i%2 else .05 for i in range(length)]

    with (root/'input.bits').open('x') as incoming, (root/'python.bits').open('x') as truth, (root/'cases.txt').open('x') as catalogue:
        for name, alpha, raw in extended():
            count += 1
            values += len(raw)
            catalogue.write(name+'\n')
            incoming.write(' '.join(bits(value) for value in [alpha, *raw])+'\n')
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    rejected, adjusted, _, _ = multipletests(np.array(raw, dtype=np.float64), alpha=alpha, method='fdr_bh')
            except Exception as error:
                failures += 1
                truth.write(type(error).__name__+' '+str(error)+'\n')
            else:
                truth.write('ok'+''.join(f' {int(reject)}:{bits(value)}' for reject, value in zip(rejected, adjusted))+'\n')
    with (root/'input.bits').open('rb') as incoming, (root/'native.bits').open('xb') as output, (root/'native.stderr.log').open('xb') as errors:
        status = subprocess.run([str(args.probe), 'fdr', '--profile', 'statsmodels'], stdin=incoming, stdout=output, stderr=errors).returncode
    mismatches = compared = 0
    with (root/'cases.txt').open() as names, (root/'python.bits').open() as truth, (root/'native.bits').open() as actual, (root/'mismatches.jsonl').open('x') as ledger:
        for line, (name, expected, observed) in enumerate(itertools.zip_longest(names, truth, actual), 1):
            compared += 1
            if expected != observed:
                mismatches += 1
                ledger.write(json.dumps(dict(line=line, case=name, python=expected, native=observed))+'\n')
    source = Path(inspect.getsourcefile(multipletests))
    result = dict(passed=not mismatches and status==0 and compared==count, cases=count, values=values,
        reference_failures=failures, mismatches=mismatches, native_returncode=status,
        numpy_version=np.__version__, statsmodels_version=statsmodels.__version__, source=str(source),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(), probe_sha256=hashlib.sha256(args.probe.read_bytes()).hexdigest(),
        scope='1D float64 multipletests fdr_bh rejection mask and adjusted bits; unused alpha outputs excluded; full DAS remains pending')
    (root/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
    if not result['passed']:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    compare_parser = commands.add_parser('compare')
    compare_parser.add_argument('--probe', type=Path, required=True)
    compare_parser.add_argument('--work-root', type=Path, required=True)
    compare_parser.add_argument('--random-cases', type=int, default=20000)
    args = parser.parse_args()
    compare(args)


if __name__ == '__main__':
    main()
