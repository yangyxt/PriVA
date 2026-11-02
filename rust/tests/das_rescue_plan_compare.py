#!/usr/bin/env python3
"""Compare rescue planning/aggregation to whole unchanged Python procedures.

compare --probe PROGRAM --reference SORTED_SOURCE --work-root NEW
The process pool supplies prescribed numerical worker results in varying
completion order. Numerical worker equality is a separate test; no dip test is
mocked into a claim about the full DAS procedure. Array identities hash every
binary64 value in source order. Error categories and summary order are exact.
"""
import argparse
import ast
import hashlib
import itertools
import json
import logging
import math
from pathlib import Path
import random
import struct
import subprocess
import sys
from types import SimpleNamespace
import warnings
import numpy as np
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore',category=RuntimeWarning)

PIN='1a2b4efc312c41c6fb20afa3c1d90566da7cd0772d2f261c2f7f09d094fae834'
NAMES={'go':'rescue_bimodal_go_mf_domains','non':'rescue_non_gomf_domains'}
def bits(value):return struct.pack('>d',value).hex()
def number(value):return struct.unpack('>d',bytes.fromhex(value))[0]
def values(array):return np.array([number(x) for x in array],dtype=np.float64)
def digest(array):return hashlib.sha256(array.astype('>f8',copy=False).tobytes()).hexdigest()


def oracle(reference):
    if hashlib.sha256(reference.read_bytes()).hexdigest()!=PIN:raise RuntimeError('reference source changed')
    functions=[node for node in ast.parse(reference.read_text()).body if isinstance(node,ast.FunctionDef) and node.name in NAMES.values()]
    if len(functions)!=2:raise RuntimeError('rescue functions absent')
    logger=logging.getLogger('plan');logger.disabled=True
    scope=dict(np=np,logger=logger,multipletests=multipletests,_dip_rescue_worker=object())
    exec(compile(ast.fix_missing_locations(ast.Module(body=functions,type_ignores=[])),str(reference),'exec'),scope)
    return scope


def evaluate(scope,case):
    global_composite=values(case['global']);composites={kind:values(array) for kind,array in case['composites']}
    scores={};entries={};errors=set();probabilities={};tolerant=set();calls=[];completed=[]
    for row in case['domains']:
        path=row['domain'];name=path.split(':',1)[-1];tolerant.add(path)
        if row['scores'] is not None:scores[path]=values(row['scores'])
        entry=row['entry'];entry=None if entry is None else ('IPR',entry[0],'short','name',[('GO','molecular_function' if entry[1] else 'biological_process','description')],'')
        if name in entries and entries[name]!=entry:raise RuntimeError('fixture annotation resolves inconsistently')
        entries[name]=entry;probabilities[path]=number(row['pvalue'])
        if row.get('lookup_error',False):errors.add(name)
    class Normalizer:
        def query_interpro_entry_vep_anno(self,name,*unused):
            calls.append(name)
            if name in errors:raise KeyError(name)
            return dict(interpro_entries=[] if entries[name] is None else [entries[name]])
    class Pool:
        def __init__(self,workers):
            if workers<=0:raise ValueError('Number of processes must be at least 1')
        def __enter__(self):return self
        def __exit__(self,*unused):pass
        def imap_unordered(self,worker,tasks):
            order=list(range(len(tasks)));random.Random(case.get('completion_seed',0)).shuffle(order)
            for index in order:
                domain,query,composite,permutations,seed=tasks[index];completed.append(domain)
                yield domain,dict(dip_query=.25,dip_ref_mean=.125,dip_ref_max=.3,p_value=probabilities[domain],n_query=len(query))
    scope['mp']=SimpleNamespace(get_context=lambda unused:SimpleNamespace(Pool=Pool))
    fn=scope[NAMES[case['mode']]];captured={}
    def profile(frame,event,arg):
        if frame.f_code is fn.__code__ and event=='return':captured.update(frame.f_locals)
    previous=sys.getprofile();sys.setprofile(profile)
    try:
        rescued,summary=fn(tolerant,scores,global_composite,Normalizer(),{},type_composites=composites,
            threads=case['threads'],n_permutations=case['permutations'],fdr_threshold=number(case['alpha']),
            min_median_am=number(case['minimum']),seed=case['seed'])
    except (ValueError,KeyError) as error:
        return dict(error=type(error).__name__,calls=calls),dict(completed=len(completed),error=type(error).__name__)
    finally:sys.setprofile(previous)
    tasks=[]
    for domain,query,composite,permutations,seed in captured.get('tasks',[]):
        kind=next((name for name,array in composites.items() if array is composite),None)
        tasks.append([domain,len(query),digest(query),kind,len(composite),digest(composite),permutations,seed])
    result=[]
    for domain,row in summary.items():
        if list(row)!=['dip_query','dip_ref_mean','dip_ref_max','p_value','n_query','fdr_corrected_p_value','rescued']:
            raise AssertionError('summary field order changed')
        result.append([domain,bits(row['p_value']),bits(row['fdr_corrected_p_value']),row['rescued']])
    direct=sorted(captured.get('directly_rescued',set()))
    return dict(error=None,calls=calls,eligible=captured.get('n_go_mf',captured.get('n_eligible_type')),
        workers=captured.get('n_cores',0),tasks=tasks,direct=direct,rescued=sorted(rescued),summary=result),dict(
            tasks=len(tasks),direct=len(direct),tested=len(summary),completed=len(completed),
            completion_reordered=completed!=[task[0] for task in tasks])


def cases():
    def row(path,kind='Family',mf=True,scores=(.8,.9),p=.5):
        return dict(domain=path,scores=None if scores is None else list(map(bits,scores)),entry=None if kind is None else [kind,mf],pvalue=bits(p))
    def case(mode,rows,**changes):
        result=dict(mode=mode,domains=rows,**{'global':list(map(bits,[.1,.2,.8]))},composites=[['Domain',list(map(bits,[.2,.9]))]],
            minimum=bits(.564 if mode=='go' else .75),seed=42 if mode=='go' else 99,threads=2,permutations=500,alpha=bits(.01))
        result.update(changes);return result
    for mode in NAMES:
        yield case(mode,[])
        for threads in [-1,0,1,12]:
            yield case(mode,[row('G:Pfam:direct','Domain',mode=='go')],threads=threads)
            for p in [float('nan'),0.,.01,.5]:
                yield case(mode,[row('G:Pfam:direct','Domain',True),row('G:Pfam:test','Family' if mode=='go' else 'Repeat',mode=='go',p=p)],threads=threads)
    rng=random.Random(0xDA55E7)
    kinds=['Domain','Family','Repeat','Homologous_superfamily','Conserved_site','Binding_site','Active_site','PTM','Unknown',None,'δ','']
    dbs=['Gene3D','Superfamily','PANTHER','Pfam','PIRSF','SMART','PROSITE_profiles','Prints','Other','PROSITE_patterns']
    prefixes=['A','Z','é','α','𐀀','😀','\x00','\n','\t','']
    samples=[[],[.564],[math.nextafter(.564,0.)],[math.nextafter(.564,1.)],[.75],[math.nextafter(.75,0.)],
        [math.nextafter(.75,1.)],[float('nan')],[.9,float('nan')],[.1,.99],None]
    probs=[0.,-0.,.001,.01,.1,.5,1.,float('nan'),float('inf'),-float('inf'),-.1,1.1]
    for index in range(3000):
        mode=rng.choice(list(NAMES));rows=[]
        for i in range(rng.randrange(40)):
            path=f'{rng.choice(prefixes)}:{rng.choice(dbs)}:{i:03d}'
            values=rng.choice(samples) if rng.randrange(3)==0 else [rng.random() for _ in range(rng.randrange(2,16))]
            value=row(path,rng.choice(kinds),bool(rng.randrange(2)),values,rng.choice(probs))
            if rng.randrange(100)==0:value['lookup_error']=True
            rows.append(value)
        if rows and rng.randrange(4)==0:rows.append(rng.choice(rows).copy())
        rng.shuffle(rows)
        comps=[]
        for kind in kinds[:-3]:
            if kind is not None and rng.randrange(2):comps.append([kind,list(map(bits,[rng.random() for _ in range(rng.randrange(8))]))])
        yield case(mode,rows,threads=rng.choice([-1,0,1,2,12]),minimum=bits(rng.choice([.564,.75,float('nan')])),
            seed=rng.choice([-1,0,42,99,2**32-2,2**32-1,2**64-1]),composites=comps,
            alpha=bits(rng.choice([.01,0.,1.,float('nan')])),completion_seed=index)


def compare(args):
    args.work_root.mkdir(parents=True,exist_ok=False);scope=oracle(args.reference);counts=dict(cases=0,tasks=0,direct=0,tested=0,completion_reordered=0)
    failures={};successes=0
    with (args.work_root/'input.jsonl').open('x') as incoming,(args.work_root/'python.jsonl').open('x') as expected:
        for case in cases():
            counts['cases']+=1;incoming.write(json.dumps(case)+'\n');result,stats=evaluate(scope,case)
            expected.write(json.dumps(result)+'\n')
            if result['error']:failures[result['error']]=failures.get(result['error'],0)+1
            else:
                successes+=1
                for name in ['tasks','direct','tested','completion_reordered']:counts[name]+=stats[name]
    with (args.work_root/'input.jsonl').open('rb') as incoming,(args.work_root/'native.jsonl').open('xb') as out,(args.work_root/'native.stderr.log').open('xb') as error:
        status=subprocess.run([str(args.probe),'plan'],stdin=incoming,stdout=out,stderr=error,timeout=300).returncode
    mismatches=0
    with (args.work_root/'python.jsonl').open() as left,(args.work_root/'native.jsonl').open() as right,(args.work_root/'mismatches.jsonl').open('x') as ledger:
        for ordinal,(a,b) in enumerate(itertools.zip_longest(left,right)):
            if a is None or b is None or json.loads(a)!=json.loads(b):
                mismatches+=1;ledger.write(json.dumps(dict(ordinal=ordinal,python=a,native=b))+'\n')
    result=dict(passed=status==0 and not mismatches,counts=counts,successes=successes,defined_failure_categories=failures,
        mismatches=mismatches,native_returncode=status,source_sha256=PIN,
        probe_sha256=hashlib.sha256(args.probe.read_bytes()).hexdigest(),
        scope='Whole unchanged rescue procedures with prescribed numerical worker results; exact task arrays/types/seeds, candidate and summary order, FDR bits, direct early returns and failure categories. Numeric workers and full DAS are separate.')
    (args.work_root/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    if not result['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    cmd=commands.add_parser('compare')
    for name in ['probe','work-root']:cmd.add_argument('--'+name,type=Path,required=True)
    cmd.add_argument('--reference',type=Path,default=Path(__file__).parents[2]/'scripts/am_pick_intolerant_domains.py')
    cmd.set_defaults(function=compare);args=parser.parse_args();args.function(args)


if __name__=='__main__':main()
