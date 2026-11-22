#!/usr/bin/env python3
"""Compare one-dimensional float64 BH adjustment against installed SciPy.

python fdr_compare.py compare --probe BIN --work-root NEW_DIR [--rmc TSV ...]
Every adjusted IEEE-754 bit and failure class/message is compared. No tolerance,
threshold normalization or scientific caller switch is used.
"""
import argparse
import csv
import hashlib
import inspect
import itertools
import json
import math
from pathlib import Path
import random
import struct
import subprocess
import time

import numpy as np
import scipy
from scipy.stats import false_discovery_control


def bits(value):
    return struct.pack('>d', float(value)).hex()


def number(value):
    return struct.unpack('>d', value.to_bytes(8, 'big'))[0]


def cases(random_cases):
    yield 'empty', []
    edge = [0., -0., math.nextafter(0., 1.), math.nextafter(1., 0.), 1.,
            .01, .025, .05, .4, math.nextafter(.05, 0.), math.nextafter(.05, 1.)]
    for i, value in enumerate(edge):
        yield 'singleton-'+str(i), [value]
    for i, value in enumerate([math.nan, math.inf, -math.inf, -math.nextafter(0., 1.), math.nextafter(1., 2.)]):
        yield 'invalid-singleton-'+str(i), [value]
        yield 'invalid-with-valid-'+str(i), [.001, value, .5]
    for length in [2, 3, 7, 15, 16, 17, 31, 32, 33, 63, 64, 65, 255, 256, 257, 4096]:
        for i, value in enumerate(edge):
            yield f'tied-{length}-{i}', [value]*length
        yield 'mixed-zero-'+str(length), [(-0. if i % 2 else 0.) for i in range(length)]
        yield 'increasing-'+str(length), np.linspace(0.,1.,length).tolist()
        yield 'decreasing-'+str(length), np.linspace(1.,0.,length).tolist()
    for i, values in enumerate(itertools.permutations([0., .001, .025, .05, 1.])):
        yield 'permutation-'+str(i), list(values)
    rng=random.Random(40780228)
    for i in range(random_cases):
        length=rng.randrange(0,258)
        mode=i % 5
        if mode==0: values=[rng.random() for _ in range(length)]
        elif mode==1: values=[rng.choice(edge) for _ in range(length)]
        elif mode==2: values=[number(rng.randrange(0x3ff0000000000001)) for _ in range(length)]
        elif mode==3: values=[rng.randrange(0,101)/100. for _ in range(length)]
        else: values=[rng.choice([.05,.01])/rng.randrange(1,500) for _ in range(length)]
        yield 'seeded-'+str(i), values


def compare(args):
    root=args.work_root.absolute();root.mkdir(parents=True,exist_ok=False)
    inputs=root/'input.bits';expected=root/'python.bits';names=root/'cases.txt';native=root/'native.bits'
    count=values=failures=0;resources=[]
    start=time.monotonic()

    def add(name, raw, incoming, truth, catalogue):
        nonlocal count, values, failures
        count+=1;values+=len(raw);catalogue.write(name+'\n')
        incoming.write(' '.join(bits(v) for v in raw)+'\n')
        try:
            adjusted=false_discovery_control(np.array(raw,dtype=np.float64),method='bh')
        except ValueError as error:
            failures+=1;truth.write('ValueError '+str(error)+'\n')
        else:
            truth.write('ok'+''.join(' '+bits(v) for v in adjusted)+'\n')

    with inputs.open('x') as incoming, expected.open('x') as truth, names.open('x') as catalogue:
        for name, raw in cases(args.random_cases):add(name,raw,incoming,truth,catalogue)
        for path in args.rmc:
            before=hashlib.sha256(path.read_bytes()).hexdigest()
            with path.open() as stream:
                raw=[float(row['p']) for row in csv.DictReader(stream,delimiter='\t')]
            add('actual-'+path.name,raw,incoming,truth,catalogue)
            after=hashlib.sha256(path.read_bytes()).hexdigest()
            if after!=before:raise RuntimeError('RMC source changed')
            resources.append(dict(path=str(path),sha256=before,rows=len(raw)))
    reference_seconds=time.monotonic()-start
    start=time.monotonic()
    with inputs.open('rb') as incoming,native.open('xb') as output,(root/'native.stderr.log').open('xb') as errors:
        process=subprocess.run([str(args.probe)],stdin=incoming,stdout=output,stderr=errors)
    native_seconds=time.monotonic()-start
    mismatch_count=0;compared=0
    with expected.open() as truth,native.open() as actual,names.open() as catalogue,(root/'mismatches.jsonl').open('x') as errors:
        for i,(name,a,b) in enumerate(itertools.zip_longest(catalogue,truth,actual),1):
            compared+=1
            if a!=b:
                mismatch_count+=1
                errors.write(json.dumps(dict(line=i,case=name,python=a,native=b))+'\n')
    source=Path(inspect.getsourcefile(false_discovery_control))
    report=dict(passed=process.returncode==0 and mismatch_count==0 and compared==count,
        cases=count,values=values,reference_failures=failures,mismatches=mismatch_count,
        native_returncode=process.returncode,python_seconds=reference_seconds,native_seconds=native_seconds,
        scipy_version=scipy.__version__,numpy_version=np.__version__,scipy_source=str(source),
        scipy_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        probe_sha256=hashlib.sha256(args.probe.read_bytes()).hexdigest(),resources=resources,
        scope='1D float64 method=bh: exact adjusted bits and ValueError class/message; full PM1 producer not yet tested')
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    if not report['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('compare')
    p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--work-root',type=Path,required=True)
    p.add_argument('--rmc',type=Path,nargs='*',default=[])
    p.add_argument('--random-cases',type=int,default=20000)
    p.set_defaults(function=compare)
    args=parser.parse_args();args.function(args)


if __name__=='__main__':main()
