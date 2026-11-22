#!/usr/bin/env python3
"""Compare ordered AF filtering to the unchanged Python source procedures.

compare --probe PROGRAM --reference SORTED_SOURCE --work-root NEW
NaN order, strict thresholds, global fallback and scalar types are exact. The
additional virtual-count cases use Python integer division/round directly, so
64-bit count boundaries do not require enormous allocated residue dictionaries.
"""
import argparse
import ast
from collections import Counter
import hashlib
import itertools
import json
import logging
import math
from pathlib import Path
import random
import struct
import subprocess
import numpy as np

PIN='1a2b4efc312c41c6fb20afa3c1d90566da7cd0772d2f261c2f7f09d094fae834'
def bits(value):return struct.pack('>d',value).hex()
def array(values):return '-1' if values is None else ' '.join([str(len(values)),*[bits(x) for x in values]])


def cases():
    values=[0.,-0.,.01,.05,math.nextafter(.01,0.),math.nextafter(.01,1.),
        math.nextafter(.05,0.),math.nextafter(.05,1.),1.,float('nan'),float('inf'),-float('inf')]
    lists=[[],*[list(pair) for pair in itertools.product(values,repeat=2)],
        [float('nan'),.1,.02],[.1,float('nan'),.02],[.1,.02,float('nan')]]
    thresholds=[(.05,.01,.05),(.05,.01,-1.),(.05,.01,0.),(.05,.01,1.),
        (float('nan'),.01,.05),(.05,float('nan'),.05),(.05,.01,float('nan'))]
    for residue,flat,threshold in itertools.product([None,[],*[ [v] for v in lists]],
            [None,[],[.01,.05],[math.nextafter(.05,1.),float('nan')]],thresholds):
        for global_present in ([False,True] if residue is None else [True]):
            yield global_present,residue,flat,*threshold
    rng=random.Random(0xAFDA5)
    for _ in range(20000):
        residue=None if rng.randrange(4)==0 else [[rng.choice(values) if rng.randrange(3)==0 else rng.random()*.1
            for _ in range(rng.randrange(8))] for _ in range(rng.randrange(30))]
        flat=None if rng.randrange(4)==0 else [rng.choice(values) for _ in range(rng.randrange(30))]
        yield residue is not None or bool(rng.randrange(2)),residue,flat,*rng.choice(thresholds)
    # Exactly one residue in twenty is retained; one in nineteen is removed.
    for total in [1,19,20,21,99,100,101,9999,10000,10001]:
        for count in [0,1,total//2,total]:
            yield True,[[.1,.02]]*count+[[]]*(total-count),None,.05,.01,.05


def encoded(summary,with_data):
    fraction=summary['frac_common']
    if type(fraction) not in [int,float]:raise AssertionError('unexpected fraction scalar')
    return ' '.join([str(summary['n_fully_common_residues']),str(summary['n_total_residues']),
        'int:'+str(fraction) if type(fraction) is int else 'float:'+bits(fraction),
        str(int(summary['removed'])),str(int(with_data))])


def compare(args):
    args.work_root.mkdir(parents=True,exist_ok=False)
    if hashlib.sha256(args.reference.read_bytes()).hexdigest()!=PIN:raise RuntimeError('reference source changed')
    names=['filter_common_missense','_collect_gnomad_residue_leaves','_collect_gnomad_leaves']
    functions=[node for node in ast.parse(args.reference.read_text()).body if isinstance(node,ast.FunctionDef) and node.name in names]
    if len(functions)!=3:raise RuntimeError('source functions absent')
    logger=logging.getLogger('common');logger.disabled=True;scope=dict(logger=logger,np=np)
    exec(compile(ast.fix_missing_locations(ast.Module(body=functions,type_ignores=[])),str(args.reference),'exec'),scope)
    count=0;paths=Counter();count_cases=0
    with (args.work_root/'input.bits').open('x') as incoming,(args.work_root/'python.bits').open('x') as expected:
        for global_present,residues,flat,common,rare,fraction in cases():
            count+=1;leaf={};model={'G':{'DB':{'ID':leaf}}}
            if residues is not None:leaf['gnomAD_residue_distribution']={i:values for i,values in enumerate(residues)}
            if flat is not None:leaf['gnomAD_distribution']=np.array(flat,dtype=np.float64)
            if global_present:model['OTHER']={'gnomAD_residue_distribution':{}}
            selected={'G:DB:ID'}
            kept,removed,summary=scope[names[0]](selected,model,common,rare,fraction)
            if kept|removed!=selected or kept&removed or (bool(removed)!=summary['G:DB:ID']['removed']):raise AssertionError('set result differs')
            if list(summary['G:DB:ID'])!=['n_fully_common_residues','n_total_residues','frac_common','removed']:raise AssertionError('field order differs')
            residue_map={};scope[names[1]](model,residue_map,[])
            use_residue=bool(residue_map.get('G:DB:ID'));use_flat=not residue_map and flat is not None
            paths['residue' if use_residue else 'flat' if use_flat else 'missing']+=1
            incoming.write(' '.join(['filter',str(int(global_present)),bits(common),bits(rare),bits(fraction),
                '-1' if residues is None else str(len(residues)),
                '' if residues is None else ' '.join(array(values) for values in residues),array(flat)])+'\n')
            expected.write(encoded(summary['G:DB:ID'],use_residue or use_flat)+'\n')
        rng=random.Random(0xC017)
        totals=[1,19,20,21,2**32-1,2**32,2**53-1,2**53,2**53+1,2**63,2**64-1]
        for _ in range(20000):
            total=rng.choice(totals) if rng.randrange(3)==0 else rng.randrange(1,2**64)
            common=rng.choice([0,1,total,total//20,total//2,rng.randrange(total+1)])
            fraction=common/total;threshold=rng.choice([.05,fraction,math.nextafter(fraction,-math.inf),math.nextafter(fraction,math.inf)])
            incoming.write(f'counts {total} {common} {bits(threshold)}\n');count_cases+=1
            expected.write(encoded(dict(n_fully_common_residues=common,n_total_residues=total,
                frac_common=round(fraction,4),removed=fraction>threshold),True)+'\n')
    with (args.work_root/'input.bits').open('rb') as incoming,(args.work_root/'native.bits').open('xb') as out,(args.work_root/'native.stderr.log').open('xb') as err:
        status=subprocess.run([str(args.probe),'filter'],stdin=incoming,stdout=out,stderr=err,timeout=300).returncode
    mismatch=0
    with (args.work_root/'python.bits').open() as left,(args.work_root/'native.bits').open() as right,(args.work_root/'mismatches.jsonl').open('x') as ledger:
        for ordinal,(a,b) in enumerate(itertools.zip_longest(left,right)):
            if a!=b:mismatch+=1;ledger.write(json.dumps(dict(ordinal=ordinal,python=a,native=b))+'\n')
    result=dict(passed=status==0 and not mismatch,cases=count+count_cases,source_filter_cases=count,virtual_count_cases=count_cases,
        source_paths=dict(paths),mismatches=mismatch,native_returncode=status,source_sha256=PIN,
        probe_sha256=hashlib.sha256(args.probe.read_bytes()).hexdigest(),
        scope='Ordered typed binary64 AF lists, global fallback, exact scalar bits/types and strict removal boundaries; full DAS separate')
    (args.work_root/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    if not result['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    cmd=commands.add_parser('compare')
    for name in ['probe','work-root']:cmd.add_argument('--'+name,type=Path,required=True)
    cmd.add_argument('--reference',type=Path,default=Path(__file__).parents[2]/'scripts/am_pick_intolerant_domains.py')
    cmd.set_defaults(function=compare);args=parser.parse_args();args.function(args)


if __name__=='__main__':main()
