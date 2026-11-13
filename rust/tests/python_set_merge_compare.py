#!/usr/bin/env python3
"""Compare native string-set copy/update and serialization stages with CPython.

compare --probe PROGRAM --work-root NEW_DIRECTORY
reference --work-root NEW_DIRECTORY

Seeds are explicit local oracle profiles. This does not recover an unknown
historical seed or qualify an entire Alpha-to-DAS installation.
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import pickle
import random
import subprocess
import sys

SEEDS=[0,1,42,4294967295]


def encode(value):return value.encode().hex() or '-'


def cases():
    sizes=[0,1,4,5,18,19,76,77,306,307,1228,1229]
    for left,right in itertools.product(sizes,repeat=2):
        for overlap in [False,True]:
            a=[f'e{i}/2000' for i in range(left)]
            b=[f'e{i+(left//2 if overlap else left)}/2000' for i in range(right)]
            yield f'size_{left}_{right}_{overlap}',a,b
    for size in [50000,50001,78643]:
        values=[f'x{i}' for i in range(size)]
        yield f'large_copy_{size}',values,[]
        yield f'large_overlap_{size}',values,values[::-1]
    chars=['','a','é','α','𐀀','😀','\0','\n','\r','\t','Aé','Aα','A😀']
    yield 'unicode',chars,chars[::-1]*3
    collision=[]
    for i in range(50000):
        text=f'collision{i}'
        if hash(text)&255==0:collision.append(text)
        if len(collision)==96:break
    if len(collision)!=96:raise RuntimeError('bounded collision fixture exhausted')
    yield 'collisions',collision[:77],collision[19:][::-1]
    rng=random.Random(16071)
    for ordinal in range(256):
        values=[rng.choice(chars)+str(rng.randrange(400)) for _ in range(rng.randrange(750))]
        left=rng.randrange(len(values)+1)
        yield f'random_{ordinal}',values[:left],values[left//2:][::-1]


def reference(args):
    if not (sys.version_info[:3]==(3,12,3) and sys.byteorder=='little'
            and sys.hash_info.algorithm=='siphash13' and sys.hash_info.width==64):
        raise RuntimeError('reference requires pinned CPython 3.12.3 64-bit SipHash13')
    root=args.work_root;root.mkdir(parents=True,exist_ok=False)
    rows=[]
    with (root/'operations.txt').open('x') as requests,(root/'python.order').open('x') as output:
        for ordinal,(name,left,right) in enumerate(cases()):
            requests.write('clear\n');sets={};stage=0
            def show(key):
                nonlocal stage
                label=f'{ordinal}:{stage}';stage+=1
                values=list(sets[key]);requests.write(f'order {label} {key}\n')
                output.write(f'{label} {len(values)} '+ ' '.join(map(encode,values))+'\n')
                rows.append(dict(case=name,stage=stage-1,members=len(values)))
            def new(key,values):
                requests.write(f'new {key}\n');sets[key]=set()
                for value in values:
                    requests.write(f'add {key} {encode(value)}\n');sets[key].add(value)
                show(key)
            def copy(dest,source):
                requests.write(f'copy {dest} {source}\n');sets[dest]=sets[source].copy();show(dest)
            def update(dest,source):
                requests.write(f'update {dest} {source}\n');sets[dest].update(sets[source]);show(dest)
            def reload(dest,source):
                requests.write(f'reload {dest} {source}\n')
                sets[dest]=pickle.loads(pickle.dumps(sets[source],protocol=5));show(dest)
            new('a',left);new('b',right);copy('original_a','a')
            update('a','b');update('a','b');show('b');copy('copied','a')
            update('b','original_a');new('empty',[]);update('empty','a')
            reload('worker','a');reload('final','worker');copy('final_copy','final')
            update('worker','original_a');reload('consumer','worker');copy('consumer_copy','consumer')
        requests.write('different_seed\n');output.write('different_seed refused\n')
    (root/'manifest.json').write_text(json.dumps(rows,indent=2)+'\n')


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def compare(args):
    root=args.work_root;root.mkdir(parents=True,exist_ok=False);results=[];mismatches=0;stages=0;members=0
    for seed in SEEDS:
        folder=root/str(seed)
        with (root/f'reference_{seed}.log').open('x') as log:
            subprocess.run([sys.executable,__file__,'reference','--work-root',str(folder)],
                           env=dict(os.environ,PYTHONHASHSEED=str(seed)),stdout=log,stderr=subprocess.STDOUT,
                           timeout=180,check=True)
        manifest=json.loads((folder/'manifest.json').read_text());stages+=len(manifest)
        members+=sum(row['members'] for row in manifest)
        with (folder/'native.stderr.log').open('x') as log:
            native=subprocess.run([str(args.probe),'batch','--input',str(folder/'operations.txt'),
                                   '--output',str(folder/'native.order'),'--seed',str(seed)],
                                  stdout=log,stderr=subprocess.STDOUT,timeout=180)
        delta=0
        with ((folder/'python.order').open() as expected,(folder/'native.order').open() as actual,
              (folder/'mismatches.jsonl').open('x') as ledger):
            for ordinal,(left,right) in enumerate(itertools.zip_longest(expected,actual)):
                if left!=right:
                    delta+=1;ledger.write(json.dumps(dict(ordinal=ordinal,
                        case=manifest[ordinal] if ordinal<len(manifest) else None,python=left,native=right))+'\n')
        mismatches+=delta
        results.append(dict(seed=seed,exit=native.returncode,stages=len(manifest),mismatches=delta,
                            operations_sha256=sha(folder/'operations.txt'),python_sha256=sha(folder/'python.order'),
                            native_sha256=sha(folder/'native.order')))
    report=dict(passed=not mismatches and all(row['exit']==0 for row in results),stages=stages,
                members=members,mismatches=mismatches,results=results,seeds=SEEDS,python=sys.version,
                binary_sha256=sha(args.probe),harness_sha256=sha(Path(__file__)),
                full_alpha_to_das_qualified=False,unknown_historical_seed_qualified=False)
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
    if not report['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    for name,function in [('reference',reference),('compare',compare)]:
        call=commands.add_parser(name);call.add_argument('--work-root',type=Path,required=True)
        if name=='compare':call.add_argument('--probe',type=Path,required=True)
        call.set_defaults(function=function)
    args=parser.parse_args();args.function(args)


if __name__=='__main__':main()
