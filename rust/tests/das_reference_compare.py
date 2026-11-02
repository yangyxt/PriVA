#!/usr/bin/env python3
"""DAS reference/Family threshold and binary64 comparisons against original ASTs.

compare --probe PROGRAM --reference SORTED_SOURCE --work-root NEW
Only producer-shaped signed-64-bit metadata is represented. Missing metadata
and explicit null protein length remain distinct; no scientific tolerance.
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

PIN='1a2b4efc312c41c6fb20afa3c1d90566da7cd0772d2f261c2f7f09d094fae834'
def bits(value):return struct.pack('>d',value).hex()
def text(value):return value.encode().hex() or '-'


def cases():
    kinds=['Family','Homologous_superfamily','Repeat','Domain','PTM','','family','δ']
    states=[('missing',1,270,300),('null',1,270,300),('some',1,270,300),('some',1,269,300),
        ('some',1,271,300),('some',1,225,250),('some',1,226,251),('some',1,0,-5),('some',1,1,0),
        ('some',-(2**63),2**63-1,2**63-1),('some',2**63-1,-2**63,-2**63),
        ('some',2**53+1,2**54+3,2**53+1),('some',-1,2**63-1,3),('some',5,1,7)]
    thresholds=[0.,.9,math.nextafter(.9,0.),math.nextafter(.9,1.),float('nan'),float('inf'),-float('inf')]
    for family,entry,kind,mf,meta,threshold,minimum in itertools.product(['reference','family'],[False,True],kinds,[False,True],states,thresholds,[-1,0,250,300]):
        yield (family,entry,kind,mf,*meta,threshold,minimum)
    rng=random.Random(0xDA510)
    special=[-2**63,-2**53-1,-2**53,-300,-1,0,1,250,300,2**53-1,2**53,2**53+1,2**63-1]
    for _ in range(50000):
        choose=lambda:rng.choice(special) if rng.randrange(3)==0 else rng.randrange(-2**63,2**63)
        yield (rng.choice(['reference','family']),True,rng.choice(kinds),rng.choice([False,True]),'some',
            choose(),choose(),choose(),rng.choice(thresholds),rng.choice(special))


def compare(args):
    args.work_root.mkdir(parents=True,exist_ok=False)
    if hashlib.sha256(args.reference.read_bytes()).hexdigest()!=PIN:raise RuntimeError('reference source changed')
    tree=ast.parse(args.reference.read_text());names=['identify_functional_domain','filter_whole_protein_family']
    functions=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in names]
    if len(functions)!=2:raise RuntimeError('source functions absent')
    logger=logging.getLogger('reference');logger.disabled=True;scope=dict(logger=logger)
    exec(compile(ast.fix_missing_locations(ast.Module(body=functions,type_ignores=[])),str(args.reference),'exec'),scope)
    class Normalizer:
        def __init__(self,entries):self.entries=entries
        def query_interpro_entry_vep_anno(self,*unused):return dict(interpro_entries=self.entries)
    count=0;failures=Counter();procedures=Counter()
    with (args.work_root/'input.bits').open('x') as incoming,(args.work_root/'python.bits').open('x') as expected:
        for family,present,kind,mf,state,start,end,length,threshold,minimum in cases():
            count+=1;procedures[family]+=1
            incoming.write(' '.join([family,str(int(present)),text(kind),str(int(mf)),state,str(start),str(end),str(length),bits(threshold),str(minimum)])+'\n')
            entry=('IPR',kind,'short','name',[('GO','molecular_function' if mf else 'biological_process','description')],'')
            normalizer=Normalizer([entry] if present else [])
            metadata={} if state=='missing' else {'G:DB:ID':dict(aa_start=start,aa_end=end,prot_length=None if state=='null' else length)}
            try:
                if family=='reference':
                    value=scope[names[0]]('G:DB:ID',normalizer,interpro_map_dict={},domain_metadata=metadata,
                        coverage_threshold=threshold,min_protein_length=minimum)
                    output='ok '+{None:'unknown','Non-functional':'non','Functional':'functional'}[value]
                else:
                    _,_,summary=scope[names[1]]({'G:DB:ID'},metadata,normalizer,{},coverage_threshold=threshold,min_protein_length=minimum)
                    value=summary['G:DB:ID']
                    if list(value)!=['interpro_type','coverage','prot_length','removed']:raise AssertionError('field order changed')
                    output=' '.join(['ok','*' if value['interpro_type'] is None else text(value['interpro_type']),
                        '*' if value['coverage'] is None else bits(value['coverage']),
                        '*' if value['prot_length'] is None else str(value['prot_length']),str(int(value['removed']))])
            except (ZeroDivisionError,TypeError) as error:
                failures[type(error).__name__]+=1;output='err '+type(error).__name__
            expected.write(output+'\n')
    with (args.work_root/'input.bits').open('rb') as incoming,(args.work_root/'native.bits').open('xb') as output,(args.work_root/'native.stderr.log').open('xb') as error:
        status=subprocess.run([str(args.probe),'classify'],stdin=incoming,stdout=output,stderr=error,timeout=300).returncode
    mismatch=0
    with (args.work_root/'python.bits').open() as left,(args.work_root/'native.bits').open() as right,(args.work_root/'mismatches.jsonl').open('x') as ledger:
        for ordinal,(a,b) in enumerate(itertools.zip_longest(left,right)):
            if a!=b:mismatch+=1;ledger.write(json.dumps(dict(ordinal=ordinal,python=a,native=b))+'\n')
    result=dict(passed=status==0 and not mismatch,cases=count,procedures=dict(procedures),defined_failures=dict(failures),
        mismatches=mismatch,native_returncode=status,source_sha256=PIN,probe_sha256=hashlib.sha256(args.probe.read_bytes()).hexdigest(),
        scope='Typed i64 metadata, exact binary64 coverage/rounded fields, threshold equalities and source exception categories; full DAS separate')
    (args.work_root/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    if not result['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    cmd=commands.add_parser('compare')
    for name in ['probe','work-root']:cmd.add_argument('--'+name,type=Path,required=True)
    cmd.add_argument('--reference',type=Path,default=Path(__file__).parents[2]/'scripts/am_pick_intolerant_domains.py')
    cmd.set_defaults(function=compare);args=parser.parse_args();args.function(args)


if __name__=='__main__':main()
