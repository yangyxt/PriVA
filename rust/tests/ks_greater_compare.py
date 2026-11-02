#!/usr/bin/env python3
"""Exact-bit comparison of the DAS one-sided KS primitive with SciPy 1.14.1.

python ks_greater_compare.py compare --probe BIN --work-root NEW_DIR
The recorded oracle inputs/results and every mismatch are retained. This
qualifies numerical primitives, not the complete DAS producer or caches.
"""
import argparse
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
import warnings

import numpy as np
import scipy
from scipy import special, stats


def bits(value):
    return struct.pack('>d',float(value)).hex()


def cases(count):
    for n in range(401):
        for k in range(n+2):
            yield f'binom-{n}-{k}', 'binom', (n,k)
    rng=random.Random(40780228)
    for i in range(count):
        n=rng.randrange(1,20001);k=rng.randrange(n+1)
        yield 'binom-random-'+str(i), 'binom', (n,k)
    edge=[[],[0.],[-0.],[1.],[-math.inf],[math.inf],[math.nan],
        [0.,1.],[0.,0.],[1.,1.],[1.,0.],[0.,math.nan],[math.nan,1.],
        [-math.inf,0.,math.inf],[-0.,0.,-0.,0.]]
    for i,(a,b) in enumerate(itertools.product(edge,repeat=2)):
        yield 'edge-'+str(i),'ks',(a,b)
    for n in [2,3,6,7,8,15,16,17,31,32,33,49,50,51,63,64,65,85,86,99,100,127,128,129,170,171,172,255,256,257,500,1000,10000,10001]:
        for shift in [0.,.001,.1,.3,1.]:
            a=np.linspace(0.,1.,n);b=a+shift
            yield f'equal-{n}-{shift}','ks',(a,b)
            yield f'equal-reverse-{n}-{shift}','ks',(b,a)
        for m in [n-1,n+1,50]:
            for shift in [0.,.1,.5,1.]:
                yield f'width-{n}-{m}-{shift}','ks',(np.linspace(0.,1.,n),np.linspace(0.,1.,m)+shift)
    for i in range(count):
        n=rng.randrange(1,151);m=n if i%3==0 else rng.randrange(1,151)
        if i%4==0:
            a=[rng.randrange(11)/10 for _ in range(n)];b=[rng.randrange(11)/10 for _ in range(m)]
        else:
            a=[rng.random() for _ in range(n)];b=[rng.random() for _ in range(m)]
        yield 'random-'+str(i),'ks',(a,b)


def compare(args):
    if scipy.__version__!='1.14.1' or np.__version__!='1.26.4':
        raise RuntimeError('reference package version changed')
    root=args.work_root;root.mkdir(parents=True,exist_ok=False)
    count=values=fallbacks=0;types={};start=time.monotonic()
    with (root/'input.bits').open('x') as incoming,(root/'python.bits').open('x') as truth,(root/'cases.txt').open('x') as catalogue:
        for name,kind,arguments in cases(args.random_cases):
            count+=1;types[kind]=types.get(kind,0)+1;catalogue.write(name+'\n')
            if kind=='binom':
                incoming.write(f'binom {arguments[0]} {arguments[1]}\n')
                truth.write('ok '+bits(special.binom(*arguments))+'\n')
            else:
                a,b=arguments;values+=len(a)+len(b)
                incoming.write('ks '+str(len(a))+''.join(' '+bits(x) for x in itertools.chain(a,b))+'\n')
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    value=stats.ks_2samp(np.asarray(a,dtype=np.float64),np.asarray(b,dtype=np.float64),alternative='greater')
                fallback=any('Exact calculation unsuccessful' in str(w.message) for w in caught)
                fallbacks+=fallback
                truth.write(f'ok {bits(value.statistic)} {bits(value.pvalue)} {int(fallback)}\n')
    reference_seconds=time.monotonic()-start;start=time.monotonic()
    with (root/'input.bits').open('rb') as incoming,(root/'native.bits').open('xb') as output,(root/'native.stderr.log').open('xb') as errors:
        status=subprocess.run([str(args.probe)],stdin=incoming,stdout=output,stderr=errors).returncode
    native_seconds=time.monotonic()-start;mismatches=0;compared=0;by_kind={}
    with (root/'cases.txt').open() as catalogue,(root/'python.bits').open() as truth,(root/'native.bits').open() as actual,(root/'mismatches.jsonl').open('x') as errors:
        for i,(name,a,b) in enumerate(itertools.zip_longest(catalogue,truth,actual),1):
            compared+=1
            if a!=b:
                mismatches+=1;kind='binom' if name and name.startswith('binom-') else 'ks';by_kind[kind]=by_kind.get(kind,0)+1
                errors.write(json.dumps(dict(line=i,case=name,python=a,native=b))+'\n')
    sources=[Path(stats._stats_py.__file__),Path(inspect.getsourcefile(stats.ks_2samp)),Path(special._ufuncs.__file__)]
    report=dict(passed=status==0 and mismatches==0 and compared==count,cases=count,types=types,values=values,
        exact_fallbacks=fallbacks,mismatches=mismatches,mismatches_by_kind=by_kind,native_returncode=status,
        python_seconds=reference_seconds,native_seconds=native_seconds,scipy_version=scipy.__version__,numpy_version=np.__version__,
        sources=[dict(path=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest()) for source in sources],
        probe_sha256=hashlib.sha256(args.probe.read_bytes()).hexdigest(),
        scope='1D float64 greater/auto: statistic/pvalue IEEE bits, exact-path fallback; positive integer binom <=20000; full DAS pending')
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
    if not report['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    command=commands.add_parser('compare')
    command.add_argument('--probe',type=Path,required=True)
    command.add_argument('--work-root',type=Path,required=True)
    command.add_argument('--random-cases',type=int,default=3000)
    args=parser.parse_args();compare(args)


if __name__=='__main__':main()
