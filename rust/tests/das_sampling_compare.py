#!/usr/bin/env python3
"""Execute unchanged DAS sampling functions against the native procedures.

python das_sampling_compare.py compare --probe BIN --reference PY_SOURCE --work-root NEW_DIR
Only the two selected function ASTs are executed; no installer, I/O or worker
pool is invoked. Every result bit, failure category/message and PCG64 state is
compared, including state consumed before a later reference fails.
"""
import argparse
import ast
import hashlib
import inspect
import itertools
import json
from pathlib import Path
import random
import struct
import subprocess
import time
from types import SimpleNamespace
from typing import List
import warnings

import numpy as np
import diptest
import scipy
from scipy import stats


PROCEDURE_SOURCE_SHA='1b7a6100c11e7f013bb2ac858d7ce3d467ddfa690ad159959898348e6e688fe7'


def bits(value):
    return struct.pack('>d',float(value)).hex()


def state(rng):
    value=rng.bit_generator.state
    return f"{value['state']['state']} {value['state']['inc']} {value['has_uint32']} {value['uinteger']}"


def arrays(values):
    return ' '.join([str(len(values)),*(bits(value) for value in values)])


def cases():
    rng=random.Random(40780228)
    special=[0.,-0.,.1,.5,1.,float('inf'),float('-inf'),float('nan'),
        struct.unpack('>d',bytes.fromhex('7ff8000000001234'))[0]]
    references=[[],[[]],[[1.],[]],[[],[1.]],[[1.],[2.]],[[.1,.2,.3],[.4,.5],[.6]],
        [special,special[::-1]],[[float(i)/10000 for i in range(10001)]],
        [[float(i)/20000 for i in range(20001)],[1.,2.]]]
    for i,refs in enumerate(references):
        for seed in [0,42,99,2**32+17,2**64-1]:
            for target in [-7,-1,0,1,2,3,7,49,50,51,200,201,1000]:
                yield f'composite-edge-{i}-{seed}-{target}','composite',seed,target,refs
    for i in range(4000):
        refs=[[rng.random() for _ in range(rng.randrange(101))] for _ in range(rng.randrange(8))]
        yield f'composite-random-{i}','composite',rng.getrandbits(64),rng.randrange(-3,1001),refs
    for n in [0,1,2,3,7,8,49,50,51,99,100,9999,10000,10001]:
        query=[rng.random() for _ in range(n)]
        for m in [0,1,2,49,50,51,10000,10001]:
            composite=[rng.random() for _ in range(m)]
            for cap in [None,-1,0,1,49,50,10001]:
                yield f'ks-boundary-{n}-{m}-{cap}','ks',42,cap,[query,composite]
    for n in [2,3,7,8,9,15,16,17,127,128,129]:
        for i in range(len(special)):
            query=[special[(k+i)%len(special)] for k in range(n)]
            for cap in [None,0,3,50]:
                yield f'ks-special-{n}-{i}-{cap}','ks',99,cap,[query,special*17]
    for i in range(2000):
        n,m=rng.randrange(0,301),rng.randrange(0,801)
        query=[rng.randrange(11)/10 for _ in range(n)]
        composite=[rng.random() for _ in range(m)]
        cap=None if i%2 else rng.randrange(-2,251)
        yield f'ks-random-{i}','ks',rng.getrandbits(64),cap,[query,composite]
    for n in [0,1,2,3,4,7,16,32,100]:
        for m in [0,1,2,7,16,100,1000]:
            query=[rng.randrange(11)/10 for _ in range(n)]
            composite=[rng.random() for _ in range(m)]
            for permutations in [-1,0,1,2,10,500]:
                for seed in [0,42,2**32+17]:
                    yield f'dip-boundary-{n}-{m}-{permutations}-{seed}','dip',seed,permutations,[query,composite]
    for i in range(500):
        n,m=rng.randrange(1,257),rng.randrange(0,2049)
        query=[rng.random() for _ in range(n)]
        composite=[rng.randrange(101)/100 for _ in range(m)]
        yield f'dip-random-{i}','dip',rng.getrandbits(64),rng.choice([0,1,2,10,31,100,500]),[query,composite]


def compare(args):
    if np.__version__!='1.26.4' or scipy.__version__!='1.14.1' or diptest.__version__!='0.11.0':raise RuntimeError('reference package version changed')
    source=args.reference.read_bytes()
    source_sha=hashlib.sha256(source).hexdigest()
    wanted=['build_fixed_equal_weight_composite','ks_vs_fixed_composite','_dip_rescue_worker']
    nodes=[node for node in ast.parse(source).body if isinstance(node,ast.FunctionDef) and node.name in wanted]
    if [node.name for node in nodes]!=wanted:raise RuntimeError('reference procedure selection changed')
    ast_sha=hashlib.sha256(ast.dump(ast.Module(body=nodes,type_ignores=[]),include_attributes=False).encode()).hexdigest()
    procedure_source='\n\n'.join(ast.get_source_segment(source.decode(),node) for node in nodes)
    procedure_sha=hashlib.sha256(procedure_source.encode()).hexdigest()
    if procedure_sha!=PROCEDURE_SOURCE_SHA:raise RuntimeError('reference sampling procedure source changed')
    namespace={'np':np,'stats':stats,'List':List,'diptest':diptest}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(args.reference),'exec'),namespace)
    root=args.work_root;root.mkdir(parents=True,exist_ok=False)
    count=values=failures=0;by_kind={};start=time.monotonic()
    with (root/'input.bits').open('x') as incoming,(root/'python.bits').open('x') as truth,(root/'cases.txt').open('x') as catalogue:
        for name,kind,seed,parameter,samples in cases():
            count+=1;by_kind[kind]=by_kind.get(kind,0)+1;values+=sum(map(len,samples));catalogue.write(name+'\n')
            prefix=f'{kind} {seed} {parameter if parameter is not None else "none"}'
            if kind=='composite':prefix+=' '+str(len(samples))
            incoming.write(prefix+' '+' '.join(arrays(sample) for sample in samples)+'\n')
            rng=np.random.default_rng(seed);samples=[np.asarray(sample,dtype=np.float64) for sample in samples]
            def track_generator(value):
                nonlocal_rng=np.random.default_rng(value)
                observed.append(nonlocal_rng)
                return nonlocal_rng
            observed=[]
            class NumpyObserver:
                random=SimpleNamespace(default_rng=track_generator)
                def __getattr__(self,name):return getattr(np,name)
            namespace['np']=NumpyObserver()
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    if kind=='composite':
                        output=namespace[wanted[0]](samples,parameter,rng)
                        result='ok '+' '.join(bits(value) for value in output)
                    elif kind=='ks':
                        output=namespace[wanted[1]](samples[0],samples[1],rng,cap=parameter)
                        result='ok '+' '.join([bits(output['ks_stat']),bits(output['p_value']),str(output['n_query']),str(output['n_comp']),bits(output['median_query']),bits(output['median_comp'])])
                    else:
                        _,output=namespace[wanted[2]](('fixture',samples[0],samples[1],parameter,seed))
                        result='ok '+' '.join([bits(output['dip_query']),bits(output['dip_ref_mean']),bits(output['dip_ref_max']),bits(output['p_value']),str(output['n_query'])])
            except Exception as error:
                failures+=1;result=f'error {type(error).__name__} {error}'
            if kind=='dip':
                if len(observed)!=1:raise RuntimeError('rescue worker RNG observation changed')
                rng=observed[0]
            truth.write(result+' | '+state(rng)+'\n')
    python_seconds=time.monotonic()-start;start=time.monotonic()
    with (root/'input.bits').open('rb') as incoming,(root/'native.bits').open('xb') as output,(root/'native.stderr.log').open('xb') as errors:
        status=subprocess.run([str(args.probe)],stdin=incoming,stdout=output,stderr=errors,timeout=300).returncode
    native_seconds=time.monotonic()-start;mismatches=compared=state_mismatches=0;first=[]
    with (root/'cases.txt').open() as catalogue,(root/'python.bits').open() as truth,(root/'native.bits').open() as actual,(root/'mismatches.jsonl').open('x') as errors:
        for i,(name,a,b) in enumerate(itertools.zip_longest(catalogue,truth,actual),1):
            compared+=1
            if not a or not b or a.split(' | ')[-1]!=b.split(' | ')[-1]:state_mismatches+=1
            if a!=b:
                mismatches+=1
                if len(first)<20:first.append(name.strip() if name else None)
                errors.write(json.dumps(dict(line=i,case=name,python=a,native=b))+'\n')
    paths=[Path(np.__file__),Path(inspect.getsourcefile(stats.ks_2samp)),Path(stats._stats_py.__file__)]
    report=dict(passed=status==0 and mismatches==0 and compared==count,cases=count,by_procedure=by_kind,values=values,
        reference_failures=failures,matched_rng_states=count-state_mismatches,rng_state_mismatches=state_mismatches,mismatches=mismatches,first_mismatches=first,
        native_returncode=status,reference_sha256=source_sha,procedure_names=wanted,
        procedure_ast_sha256=ast_sha,procedure_source_sha256=procedure_sha,
        numpy_version=np.__version__,scipy_version=scipy.__version__,diptest_version=diptest.__version__,
        sources=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],
        probe_sha256=hashlib.sha256(args.probe.read_bytes()).hexdigest(),python_seconds=python_seconds,native_seconds=native_seconds,
        scope='1D float64 composite, KS and complete dip-rescue numerical procedure ASTs shared by original and sorted references; full DAS pending')
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
    if not report['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    command=commands.add_parser('compare');command.add_argument('--probe',type=Path,required=True)
    command.add_argument('--reference',type=Path,required=True);command.add_argument('--work-root',type=Path,required=True)
    args=parser.parse_args();compare(args)


if __name__=='__main__':main()
