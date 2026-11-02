#!/usr/bin/env python3
"""Compare the source DAS benign filter's values and exact failure boundary.

check --binary PRIVA --source ORIGINAL.py --fixture-helper AA_FIXTURES.py
      --work-root NEW_DIRECTORY --report JSON
oracle invokes the unchanged source function; both arms emit requested domain
order. This tests a filter procedure, not the full DAS candidate-set history.
"""
import argparse
import ast
from collections import defaultdict
import hashlib
import importlib.util
import json
import logging
from pathlib import Path
import re
import shutil
import subprocess
import sys

FIELDS=['DOMAINS','Gene','Consequence','Feature_type','Feature','TAIL']
DEFAULTS=dict(DOMAINS='DB:A:B',Gene='ENSG',Consequence='missense_variant',Feature_type='Transcript',Feature='ENST_ONE',TAIL='end')
DOMAINS=['ENSG:DB:A','ENSG:DB:A:B','ENSG:DB:X','MISSING:DB:A','ENSG:DB:A','ENSG::A','ENSG:DB:','ENSG_α:DB:Δ']


def annotation(fields=FIELDS,**values):
    data=dict(DEFAULTS,**values);return '|'.join(data.get(key,'value') for key in fields)


def row(pos=10,csq=None,sig='Benign',review='reviewed_by_expert_panel',extra='',chrom='chr1'):
    return f'{chrom}\t{pos}\t.\tA\tT\t.\tPASS\tCLNSIG={sig};CLNREVSTAT={review};CSQ={annotation() if csq is None else csq}{extra}\n'


def cases():
    result=dict(normal={},empty=dict(rows=[]),three_occurrences=dict(rows=[row(),row(pos=11),row(pos=12)]),
        repeated_annotations=dict(rows=[row(csq=annotation()+','+annotation()+','+annotation(Feature='ENST_TWO'))]),
        repeated_prefixes=dict(rows=[row(csq=annotation(DOMAINS='DB:A:B&DB:A:B&DB:A'))]),
        overlapping_query_repeats={},empty_candidates=dict(domains=[]),minimum_zero=dict(minimum=0),minimum_negative=dict(minimum=-1),minimum_large=dict(minimum=99),
        plain=dict(compression='plain',index=False),gzip=dict(compression='gzip',index=False),no_index=dict(index=False),no_eof=dict(no_eof=True),
        empty_no_contigs=dict(contigs=[],rows=[],index=False),no_contigs=dict(contigs=[],index=False),
        unicode=dict(rows=[row(csq=annotation(Gene='ENSG_α',DOMAINS='DB:Δ',Feature='ENST_β'))]),
        empty_hierarchy=dict(rows=[row(csq=annotation(DOMAINS=':A&DB:&DB:A::B&DB&'))]),
        wrong_feature_type=dict(rows=[row(csq=annotation(Feature_type='RegulatoryFeature'))]),
        wrong_feature=dict(rows=[row(csq=annotation(Feature='NM_ONE'))]),
        no_missense=dict(rows=[row(csq=annotation(Consequence='synonymous_variant'))]),
        missense_substring=dict(rows=[row(csq=annotation(Consequence='prefix_missense_variant_suffix'))]),
        empty_gene=dict(rows=[row(csq=annotation(Gene=''))]),empty_domains=dict(rows=[row(csq=annotation(DOMAINS=''))]),
        duplicate_info_first=dict(rows=[row(extra=';CLNSIG=Pathogenic;CLNREVSTAT=other;CSQ=short')]),
        short_annotation=dict(rows=[row(csq='short')]),short_nontx=dict(rows=[row(csq='A|B|C|notTranscript')]),
        short_non_enst=dict(rows=[row(csq='A|B|C|Transcript|NM_1')]),
        short_then_valid=dict(rows=[row(csq='short,'+annotation())]),valid_then_missing=dict(rows=[row(csq=annotation()+',,'+annotation())]),
        no_format=dict(description='no format annotation'),
        other_last_format=dict(extra_header='##INFO=<ID=OTHER,Number=.,Type=String,Description="Format: Feature|TAIL">\n'),
        duplicate_header=dict(extra_header='##INFO=<ID=CSQ,Number=1,Type=Integer,Description="Format: OTHER|TAIL">\n'),
        duplicate_fields=dict(fields=FIELDS[:-1]+['Gene','Feature','TAIL']),
        polluted_last_required=dict(fields=FIELDS[:-1]),crlf=dict(crlf=True),
        no_clinical_info=dict(rows=['chr1\t10\t.\tA\tT\t.\tPASS\t.\n']),
        absent_csq=dict(rows=[row().split(';CSQ=')[0]+'\n']),
        invalid_unused_gt=dict(sample='GT',sample_value='bad',rows=[row(sig='Pathogenic')]),
        invalid_unused_format=dict(sample='DP',sample_value='bad',format_type='Integer',rows=[row(sig='Pathogenic')]),
        negative_pos=dict(mutate=('chr1\t10\t','chr1\t-1\t')),bad_pos=dict(mutate=('chr1\t10\t','chr1\tbad\t')),
        empty_pos=dict(mutate=('chr1\t10\t','chr1\t\t')),missing_file=dict(missing=True))
    for name in FIELDS[:-1]:result['missing_'+name]=dict(fields=[x for x in FIELDS if x!=name])
    for prefix in ['sig','review','csq']:
        for number,kind in [('1','String'),('.','Integer'),('1','Integer'),('.','Float'),('1','Float'),('0','Flag'),('G','String')]:
            result[f'{prefix}_{number}_{kind}']={prefix+'_number':number,prefix+'_type':kind}
        result[prefix+'_undefined']=dict(omit_headers=[dict(sig='CLNSIG',review='CLNREVSTAT',csq='CSQ')[prefix]])
    for sig in ['', '.', 'Likely_benign','Benign/Likely_benign','Benign,Benign','Pathogenic','Benign,']:
        result['significance_'+str(len(result))]=dict(rows=[row(sig=sig)])
    for review in ['', '.', 'practice_guideline','criteria_provided,_multiple_submitters,_no_conflicts','reviewed_by_expert_panel,practice_guideline','reviewed_by_expert_panel,']:
        result['review_'+str(len(result))]=dict(rows=[row(review=review)])
    for csq in ['', '.', '0']:
        for number in ['.','1']:result['csq_empty_'+str(len(result))]=dict(csq_number=number,rows=[row(csq=csq)])
    for size in [5,6,7]:result['short_record_'+str(size)]=dict(mutate=(row(),'\t'.join(row().rstrip().split('\t')[:size])+'\n'))
    return result


def original(source):
    import pysam
    body=[node for node in ast.parse(source.read_text()).body if isinstance(node,ast.FunctionDef) and node.name=='filter_benign_clinvar']
    if len(body)!=1:raise ValueError('original benign function absent')
    logger=logging.getLogger('benign_oracle');logger.disabled=True
    scope=dict(pysam=pysam,defaultdict=defaultdict,logger=logger)
    exec(compile(ast.fix_missing_locations(ast.Module(body=body,type_ignores=[])),str(source),'exec'),scope)
    return scope['filter_benign_clinvar']


def oracle(args):
    domains=[json.loads(line) for line in args.domains.read_text().splitlines()]
    _,_,summary=original(args.source)(set(domains),str(args.vcf),args.minimum)
    for domain in domains:print(json.dumps(dict(domain=domain,**summary[domain])),flush=True)


def identity(path):
    return dict(path=str(path),bytes=path.stat().st_size,sha256=hashlib.file_digest(path.open('rb'),'sha256').hexdigest())


def category(path):
    text=path.read_text(errors='replace')
    names=re.findall(r'\b([A-Za-z]+(?:Error|Exception)):',text)
    return names[-1] if names else None


def typed(value):
    if type(value) is dict:return ['dict',[[typed(key),typed(item)] for key,item in value.items()]]
    if type(value) is list:return ['list',[typed(item) for item in value]]
    return [type(value).__name__,value]


def check(args):
    root=args.work_root;root.mkdir(parents=True,exist_ok=False)
    pins=[identity(path) for path in [Path(__file__),args.binary,args.source,args.fixture_helper]]
    shutil.copy2(__file__,root/Path(__file__).name)
    spec=importlib.util.spec_from_file_location('clinvar_vcf_fixture_inputs',args.fixture_helper)
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    results=[];failures=[]
    for name,specification in cases().items():
        folder=root/name;folder.mkdir();specification=dict(specification)
        specification.setdefault('fields',FIELDS)
        specification.setdefault('rows',[row(csq=annotation(specification['fields']))])
        vcf=helper.make_case(folder,specification)
        if specification.get('missing'):vcf=folder/'missing.vcf'
        domains=folder/'domains.jsonl';domains.write_text(''.join(json.dumps(value)+'\n' for value in specification.get('domains',DOMAINS)))
        minimum=str(specification.get('minimum',3))
        commands=[('python',[sys.executable,__file__,'oracle','--source',args.source,'--vcf',vcf,'--domains',domains,'--minimum',minimum]),
            ('native',[args.binary,'das-domains','benign','--vcf',vcf,'--domains',domains,'--min-benign',minimum])]
        statuses=[]
        for label,command in commands:
            with (folder/(label+'.stdout.log')).open('x') as out,(folder/(label+'.stderr.log')).open('x') as err:
                statuses.append(subprocess.run(list(map(str,command)),stdout=out,stderr=err,timeout=60).returncode)
        categories=[category(folder/(label+'.stderr.log')) if code else None for (label,_),code in zip(commands,statuses)]
        passed=statuses[0]==statuses[1] and (not statuses[0] or categories[0] is not None and categories[0]==categories[1])
        differences=[];compared=0
        if statuses==[0,0]:
            expected=[json.loads(line) for line in (folder/'python.stdout.log').read_text().splitlines()]
            actual=[json.loads(line) for line in (folder/'native.stdout.log').read_text().splitlines()]
            if len(expected)!=len(actual):differences.append(dict(field='row_count',python=len(expected),native=len(actual)))
            for ordinal,(left,right) in enumerate(zip(expected,actual)):
                compared+=1
                if typed(left)!=typed(right):differences.append(dict(ordinal=ordinal,python=left,native=right))
            passed=passed and not differences
        elif (folder/'native.stdout.log').stat().st_size:
            passed=False;differences.append(dict(field='partial_native_output'))
        (folder/'mismatches.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in differences))
        result=dict(case=name,passed=passed,python_exit=statuses[0],native_exit=statuses[1],python_category=categories[0],native_category=categories[1],compared=compared,mismatches=len(differences))
        results.append(result)
        if not passed:failures.append(result)
    for pin in pins:
        if identity(Path(pin['path']))!=pin:raise RuntimeError('comparison input/source changed')
    result=dict(passed=not failures,cases=len(results),matched=sum(row['passed'] for row in results),
        successes=sum(row['python_exit']==0 for row in results),defined_failures=sum(row['python_exit']==1 for row in results),
        compared=sum(row['compared'] for row in results),failures=failures,results=results,inputs=pins,
        output_order='requested domain order; full DAS candidate-set history is separate',full_das_qualified=False)
    (root/'report.json').write_text(json.dumps(result,indent=2)+'\n');args.report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['passed','cases','matched','successes','defined_failures','compared','failures']}),flush=True)
    if not result['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    check_cmd=commands.add_parser('check')
    for name in ['binary','source','fixture-helper','work-root','report']:check_cmd.add_argument('--'+name,type=Path,required=True)
    check_cmd.set_defaults(function=check)
    oracle_cmd=commands.add_parser('oracle')
    for name in ['source','vcf','domains']:oracle_cmd.add_argument('--'+name,type=Path,required=True)
    oracle_cmd.add_argument('--minimum',type=int,required=True);oracle_cmd.set_defaults(function=oracle)
    args=parser.parse_args();args.function(args)


if __name__=='__main__':main()
