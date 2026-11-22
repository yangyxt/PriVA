#!/usr/bin/env python3
"""Compare the DAS-local InterPro lookup with the unchanged Python procedures.

compare --binary PRIVA --source PROTEIN_DOMAIN_MAPPING --normalization-json JSON --work-root NEW
reference --source SOURCE --normalization-json JSON --mapping PICKLE --domains JSONL --output JSONL
The source-produced schema-9 profile uses exact six-tuples and tuple GO terms.
Full DAS, constructor side effects and arbitrary nonproducer pickle shapes are
separate gates; this test preserves normalized-key and source-entry order.
"""
import argparse
import ast
import hashlib
import itertools
import json
import logging
import os
from pathlib import Path
import pickle
import random
import re
import subprocess
import sys
from typing import Optional,Tuple

PIN='477a515793cf034db30ba94aee08d1c0afb2508eb0520357ec1ab7ba8e28ac7b'
SEEDS=[0,1,42,4294967295]


def normalizer(source,document):
    if hashlib.sha256(source.read_bytes()).hexdigest()!=PIN:raise RuntimeError('Python lookup source changed')
    tree=ast.parse(source.read_text());cls=next(node for node in tree.body if isinstance(node,ast.ClassDef) and node.name=='DomainNormalizer')
    functions=[node for node in cls.body if isinstance(node,ast.FunctionDef) and node.name in ['normalize_vep_to_interpro','query_interpro_entry_vep_anno']]
    if len(functions)!=2:raise RuntimeError('original lookup procedures absent')
    scope=dict(re=re,Optional=Optional,Tuple=Tuple,FinalMapping=dict)
    declaration=ast.ClassDef(name='SourceNormalizer',bases=[],keywords=[],body=functions,decorator_list=[])
    exec(compile(ast.fix_missing_locations(ast.Module(body=[declaration],type_ignores=[])),str(source),'exec'),scope)
    subject=scope['SourceNormalizer']();subject.logger=logging.getLogger('lookup_oracle');subject.logger.disabled=True
    subject._vep_to_ipr_map=document['database_name_mappings']['vep_to_interpro']
    subject._vep_to_ipr_id_rules=document['accession_id_normalization']['vep_to_interpro_key']
    return subject


def reference(args):
    subject=normalizer(args.source,json.loads(args.normalization_json.read_text()))
    with args.mapping.open('rb') as stream:mapping=pickle.load(stream)
    with args.domains.open() as source,args.output.open('x') as output:
        for line in source:
            domain=json.loads(line);annotation=domain.split(':',1)[-1]
            result=subject.query_interpro_entry_vep_anno(annotation,mapping);entries=result.get('interpro_entries')
            entry=None
            if entries:
                value=entries[0]
                entry=dict(id=value[0],type=value[1],go_molecular_function=any(term[1]=='molecular_function' for term in value[4]))
            output.write(json.dumps(dict(domain=domain,entry=entry))+'\n')


def fixtures(subject):
    annotations=['Pfam:PF00001','Pfam:PF00002','Pfam:PF00003','PANTHER:PTHR12345:SF1','PANTHER:PTHR12345:SF2',
        'PANTHER:PTHR23456','PANTHER:bad','Gene3D:1.2.3.4','Gene3D:G3DSA:1.2.3.4','SUPFAM:SSF1','Superfamily:SSF1',
        'Superfamily:SSF2','SMART:SM1','Prints:PR1','PIRSF:PIRSF1','PROSITE_profiles:PS1','PROSITE_patterns:PS2',
        'TIGRFAM:NF1','NCBIFAM:NF1','SFLD:S1','CDD:C1','HAMAP:H1','Pfam:α','Pfam:😀','Pfam:a:b',
        'Pfam:','Unknown:X','NoColon','','Transmembrane_helices:TM1','Low_complexity_(Seg):L1']
    keys=[]
    for annotation in annotations:
        if ':' not in annotation:continue
        value=subject.normalize_vep_to_interpro(*annotation.split(':',1))
        if value is not None and value not in keys:keys.append(value)
    kinds=['Domain','Family','Repeat','Homologous_superfamily','Conserved_site','Binding_site','PTM','Active_site','','δ']
    mapping={}
    for index,key in enumerate(keys):
        terms=[('GO:1','biological_process','desc'),('GO:2','molecular_function' if index%2 else 'cellular_component','description')]
        first=(f'IPR{index}',kinds[index%len(kinds)],'short','full',terms,'abstract')
        second=(f'IPR_other{index}',kinds[(index+1)%len(kinds)],'short','full',[('GO:3','molecular_function','other')],'abstract')
        mapping[key]=[] if index%7==0 else [first,first,second]
    domains=['G:'+value for value in annotations]
    domains+=['G:'+a+'&'+b for a,b in itertools.product(annotations,repeat=2)]
    domains+=['G:'+a+'&'+b+'&'+a for a,b in itertools.product(annotations[::2],repeat=2)]
    domains+=['',':','G','Pfam:PF00001','α:Pfam:α','😀:Pfam:😀','G:Pfam:PF00001\n','G:Pfam:PF00001\t']
    rng=random.Random(0x1A7E)
    for _ in range(4000):domains.append('G:'+'&'.join(rng.choices(annotations,k=rng.randrange(1,14))))
    return mapping,domains


def compare(args):
    args.work_root.mkdir(parents=True,exist_ok=False);root=args.work_root
    subject=normalizer(args.source,json.loads(args.normalization_json.read_text()));mapping,domains=fixtures(subject)
    with (root/'mapping.pkl').open('xb') as stream:pickle.dump(mapping,stream,protocol=4)
    (root/'domains.jsonl').write_text(''.join(json.dumps(value)+'\n' for value in domains))
    with (root/'convert.stdout.log').open('xb') as out,(root/'convert.stderr.log').open('xb') as error:
        subprocess.run([str(args.binary),'cache-convert','interpro_mapping_pickle','--pickle',str(root/'mapping.pkl'),
            '--output',str(root/'mapping.rkyv')],stdout=out,stderr=error,check=True)
    rows=[];mismatches=[]
    for seed in SEEDS:
        expected=root/f'python_{seed}.jsonl';actual=root/f'native_{seed}.jsonl'
        with (root/f'python_{seed}.stderr.log').open('xb') as error:
            subprocess.run([sys.executable,__file__,'reference','--source',str(args.source),
                '--normalization-json',str(args.normalization_json),'--mapping',str(root/'mapping.pkl'),
                '--domains',str(root/'domains.jsonl'),'--output',str(expected)],
                env=dict(os.environ,PYTHONHASHSEED=str(seed)),stdout=subprocess.DEVNULL,stderr=error,check=True)
        with actual.open('xb') as out,(root/f'native_{seed}.stderr.log').open('xb') as error:
            status=subprocess.run([str(args.binary),'das-domains','lookup','--domains',str(root/'domains.jsonl'),
                '--interpro',str(root/'mapping.rkyv'),'--normalization-json',str(args.normalization_json),
                '--python-hash-seed',str(seed)],stdout=out,stderr=error).returncode
        count=0
        with expected.open() as left,actual.open() as right:
            for ordinal,(a,b) in enumerate(itertools.zip_longest(left,right)):
                count+=1
                if a!=b:mismatches.append(dict(seed=seed,ordinal=ordinal,python=a,native=b))
        rows.append(dict(seed=seed,returncode=status,rows=count))
    (root/'mismatches.json').write_text(json.dumps(mismatches,indent=2)+'\n')
    result=dict(passed=not mismatches and all(row['returncode']==0 and row['rows']==len(domains) for row in rows),
        cases=len(domains)*len(SEEDS),mapping_keys=len(mapping),mismatches=len(mismatches),results=rows,
        reference_sha256=PIN,normalization_sha256=hashlib.sha256(args.normalization_json.read_bytes()).hexdigest(),
        binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(),
        scope='Source-produced schema-9 tuples, exact ordered first entry and tuple GO-MF projection; full DAS remains separate')
    (root/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    if not result['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    for name,function in [('reference',reference),('compare',compare)]:
        cmd=commands.add_parser(name)
        fields=['source','normalization-json']+(['mapping','domains','output'] if name=='reference' else ['binary','work-root'])
        for field in fields:cmd.add_argument('--'+field,type=Path,required=True)
        cmd.set_defaults(function=function)
    args=parser.parse_args();args.function(args)


if __name__=='__main__':main()
