#!/usr/bin/env python3
"""Independent exact-bit oracle for Hartigan's dip statistic.

python dip_stat_compare.py compare --probe BIN --work-root NEW_DIR
diptest 0.11.0 is loaded by this development-only Python comparator. It is
neither linked to the Rust executable nor copied into its implementation.
"""
import argparse
import hashlib
import itertools
import inspect
import json
import math
from pathlib import Path
import random
import struct
import subprocess
import sys
import time

import diptest
import numpy as np


def bits(value):
    return struct.pack('>d',float(value)).hex()


def cases(count):
    yield 'empty',[]
    edge=[0.,-0.,1.,-.5,math.nextafter(0.,1.),math.nextafter(1.,0.),math.nextafter(1.,2.)]
    for n in [1,2,3,4,5,6,7,8,15,16,17,31,32,33,63,64,65,127,128,129,255,256,257,1000,10000]:
        for i,value in enumerate(edge):
            yield f'constant-{n}-{i}',[value]*n
        for label,values in [('ramp',list(range(n))),('squared',[i*i for i in range(n)]),
            ('sqrt',[math.sqrt(i) for i in range(n)]),('two',[0.]*(n//2)+[1.]*(n-n//2)),
            ('two-skew',[0.]*(n//4)+[1.]*(n-n//4)),('three',[float(i%3) for i in range(n)])]:
            yield f'{label}-{n}',values
            yield f'{label}-mirrored-{n}',[-float(v) for v in values]
    for n in range(1,9):
        for i,values in enumerate(itertools.combinations_with_replacement([0.,1.,2.,3.],n)):
            yield f'grid-{n}-{i}',values
    rng=random.Random(40780228)
    for i in range(count):
        n=rng.randrange(1,513)
        if i%4==0:values=[rng.randrange(11)/10 for _ in range(n)]
        elif i%4==1:values=[rng.random() for _ in range(n)]
        elif i%4==2:values=[rng.gauss(0.,1.) for _ in range(n)]
        else:values=[rng.choice([-.5,.5])+rng.gauss(0.,.1) for _ in range(n)]
        yield 'random-'+str(i),values
    for grid_index,grid in enumerate([[0.,.1,.2,.3],[-.1,0.,.1,.2],
        [0.,math.nextafter(.1,0.),.1,math.nextafter(.1,1.)]]):
        for n in range(1,13):
            for i,values in enumerate(itertools.combinations_with_replacement(grid,n)):
                yield f'decimal-grid-{grid_index}-{n}-{i}',values
    reductions=[[0,1,0,1,1,1,0,1,1,1,0],
        [0,41,43,28,41,40,33,40,42,0,0],[12,16,11,11,10,13,8,14,14,0,0],
        [0,0,46,49,40,50,1,0,0,0,0],[0,0,5,6,2,6,6,1,6,6,0],
        [0,0,16,27,27,1,6,5,0,0,0]]
    for i,counts in enumerate(reductions):
        sample=[level/10 for level,count in enumerate(counts) for _ in range(count)]
        for exponent in [-1000,-500,-1,0,1,500,1000]:
            transformed=[math.ldexp(x,exponent) for x in sample]
            yield f'reduced-{i}-scaled-{exponent}',transformed
            yield f'reduced-{i}-reflected-{exponent}',[-x for x in reversed(transformed)]
        for j in range(10):
            rng.shuffle(sample);yield f'reduced-{i}-permutation-{j}',sample.copy()
    finite=[0.,-0.,math.nextafter(0.,1.),-math.nextafter(0.,1.),sys.float_info.min,
        -sys.float_info.min,sys.float_info.max,-sys.float_info.max]
    for n in [2,3,4,5,7,8,16,31,128]:
        for offset in range(len(finite)):
            yield f'finite-extremes-{n}-{offset}',[finite[(offset+i)%len(finite)] for i in range(n)]
    for i in range(1000):
        sample=[]
        for _ in range(rng.randrange(2,65)):
            word=rng.getrandbits(64)
            if (word>>52)&2047==2047:word^=1<<52
            sample.append(struct.unpack('>d',word.to_bytes(8,'big'))[0])
        yield f'finite-bits-{i}',sample
    for special in [float('inf'),float('-inf'),float('nan')]:
        for n in [1,2,3,4,10]:
            yield f'nonfinite-constant-{special}-{n}',[special]*n
        for position in range(5):
            sample=[0.,.1,.5,.7,1.];sample[position]=special
            yield f'nonfinite-mixed-{special}-{position}',sample
    yield 'nonfinite-both-infinities',[-float('inf'),0.,float('inf')]


def oracle():
    values=[struct.unpack('>d',bytes.fromhex(word))[0] for word in sys.stdin.read().split()]
    try:value=diptest.dipstat(np.asarray(values,dtype=np.float64),allow_zero=True)
    except Exception as error:print('error '+type(error).__name__)
    else:print('ok '+bits(value))


def compare(args):
    if diptest.__version__!='0.11.0' or np.__version__!='1.26.4':
        raise RuntimeError('reference version changed')
    root=args.work_root;root.mkdir(parents=True,exist_ok=False)
    count=values=failures=0;start=time.monotonic();boundaries=[];isolated=[]
    with (root/'input.bits').open('x') as incoming,(root/'python.bits').open('x') as truth,(root/'cases.txt').open('x') as catalogue:
        for name,sample in cases(args.random_cases):
            count+=1;values+=len(sample);catalogue.write(name+'\n')
            request=' '.join(bits(x) for x in sample)+'\n';incoming.write(request)
            if any(not math.isfinite(value) for value in sample):
                with (root/(name+'.oracle.stderr.log')).open('xb') as error_log:
                    try:
                        child=subprocess.run([sys.executable,__file__,'oracle'],input=request.encode(),stdout=subprocess.PIPE,stderr=error_log,timeout=3)
                    except subprocess.TimeoutExpired:
                        boundaries.append(dict(case=name,status='timeout',seconds=3));truth.write('unqualified timeout\n')
                    else:
                        isolated.append(dict(case=name,returncode=child.returncode))
                        if child.returncode:
                            boundaries.append(dict(case=name,status='exit',returncode=child.returncode));truth.write('unqualified exit '+str(child.returncode)+'\n')
                        else:
                            result=child.stdout.decode();truth.write(result);failures+=result.startswith('error ')
                continue
            try:value=diptest.dipstat(np.asarray(sample,dtype=np.float64),allow_zero=True)
            except Exception as error:
                failures+=1;truth.write('error '+type(error).__name__+'\n')
            else:truth.write('ok '+bits(value)+'\n')
    python_seconds=time.monotonic()-start;start=time.monotonic()
    with (root/'input.bits').open('rb') as incoming,(root/'native.bits').open('xb') as output,(root/'native.stderr.log').open('xb') as errors:
        status=subprocess.run([str(args.probe)],stdin=incoming,stdout=output,stderr=errors,timeout=300).returncode
    native_seconds=time.monotonic()-start;mismatches=compared=0;categories={};first=[]
    with (root/'cases.txt').open() as catalogue,(root/'python.bits').open() as truth,(root/'native.bits').open() as actual,(root/'mismatches.jsonl').open('x') as errors:
        for i,(name,a,b) in enumerate(itertools.zip_longest(catalogue,truth,actual),1):
            compared+=1
            if a!=b:
                mismatches+=1
                if a and b and a.startswith('ok ') and b.startswith('ok '):
                    left=int(a.split()[1],16);right=int(b.split()[1],16);distance=abs(left-right)
                    category='one_ulp' if distance==1 else ('zero_vs_nonzero' if not left or not right else 'larger')
                else:category='failure'
                categories[category]=categories.get(category,0)+1
                if len(first)<25:first.append(dict(case=name.strip() if name else None,category=category))
                errors.write(json.dumps(dict(line=i,case=name,python=a,native=b))+'\n')
    import diptest.lib._diptest_core as core
    sources=[Path(diptest.__file__),Path(inspect.getsourcefile(diptest.dipstat)),Path(core.__file__)]
    report=dict(passed=status==0 and mismatches==0 and compared==count and not boundaries,cases=count,values=values,
        reference_failures=failures,mismatches=mismatches,mismatch_categories=categories,first_mismatches=first,
        native_returncode=status,python_seconds=python_seconds,native_seconds=native_seconds,
        isolated_oracles=isolated,unqualified_resource_boundaries=boundaries,
        diptest_version=diptest.__version__,numpy_version=np.__version__,
        sources=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources],
        probe_sha256=hashlib.sha256(args.probe.read_bytes()).hexdigest(),scope='1D float64 dipstat allow_zero=True; full DAS pending')
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
    if not report['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    command=commands.add_parser('compare')
    command.add_argument('--probe',type=Path,required=True)
    command.add_argument('--work-root',type=Path,required=True)
    command.add_argument('--random-cases',type=int,default=10000)
    commands.add_parser('oracle')
    args=parser.parse_args()
    if args.command=='oracle':oracle()
    else:compare(args)


if __name__=='__main__':main()
