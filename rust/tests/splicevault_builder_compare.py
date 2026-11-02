#!/usr/bin/env python3
"""Differential SpliceVault fixtures using real pandas, pysam and CrossMap.

python splicevault_builder_compare.py fixtures --binary PATH --python-script PATH --work-root PATH --report PATH

Each arm has an isolated output/temp directory. The transparent CrossMap
launcher records argv and execs the real installed tool without changing its
arguments. Compare exact table bytes, exit status, called stages and retained
sidecars; normalize only the arm directory, UUID and upstream timestamp.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import time

HEADER='#chrom\tpos\tref\talt\tsplicevault_site_loc'
UUID=re.compile(r'[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}')

def sha(data):return hashlib.sha256(data).hexdigest()

def normalize(data,arm):
    data=UUID.sub('<UUID>',data.decode('utf-8')).replace(str(arm),'<ARM>')
    return re.sub(r'(?m)^##liftOverDate=.*$','##liftOverDate=<DATE>',data).encode()

def fixture_cases():
    row='chr1\t20\tA\tC\tchr1:21'
    second='chr1\t30\tA\tG\tchr1:31'
    third='chr2\t40\tG\tT\tchr2:41'
    cases=[]
    def add(name,text=None,chunk=2,**kwargs):
        cases.append(dict(name=name,text=text if text is not None else HEADER+'\n'+row+'\n',chunk=chunk,**kwargs))
    add('ordinary');add('three_chunks',HEADER+'\n'+row+'\n'+second+'\n'+third+'\n',chunk=1)
    add('reverse_strand',HEADER+'\n'+third+'\n')
    add('duplicate_variant_ids',HEADER+'\textra\n'+row+'\tfirst\n'+row.replace('chr1:21','chr1:31')+'\tlast\n')
    add('duplicate_ids_across_chunks',HEADER+'\textra\n'+row+'\tfirst\n'+row.replace('chr1:21','chr1:31')+'\tlast\n',chunk=1)
    add('vcf_unmapped',HEADER+'\n'+row.replace('\t20\t','\t1200\t')+'\n'+row+'\n')
    add('bed_unmapped',HEADER+'\n'+row.replace('chr1:21','chr1:1200')+'\n')
    add('all_unmapped',HEADER+'\n'+row.replace('\t20\t','\t1200\t')+'\n')
    add('header_only',HEADER+'\n');add('empty','');add('blank_only',' \n\n')
    add('empty_first_header','\n'+HEADER+'\n'+row+'\n')
    add('bom','\ufeff'+HEADER+'\n'+row+'\n')
    add('crlf',(HEADER+'\n'+row+'\n').replace('\n','\r\n'))
    add('bare_cr',(HEADER+'\n'+row+'\n').replace('\n','\r'))
    add('spaces_and_blank',HEADER+'\n \n\n'+row+'\n  \n')
    add('trailing_spaces_no_newline',HEADER+'\n'+row+'\n   ')
    add('header_trailing_tab',HEADER+'\t\n'+row+'\t007\n')
    add('duplicate_columns',HEADER+'\tpos\tvariant_id\tvariant_id\n'+row+'\t009\told\tkeep\n')
    add('implicit_duplicate_index',HEADER+'\nlabel\t'+row+'\nlabel\t'+second+'\n')
    add('implicit_multi_index',HEADER+'\na\tb\t'+row+'\na\tb\t'+second+'\n')
    add('implicit_unique_multi_index',HEADER+'\na\tb\t'+row+'\na\tc\t'+second+'\n')
    add('implicit_index_with_unmapped',HEADER+'\na\t'+row+'\na\t'+second.replace('\t30\t','\t1200\t')+'\n')
    for size in (1,2,3):
        add('extra_later_chunk'+str(size),HEADER+'\n'+row+'\n'+second+'\textra\n'+third+'\n',chunk=size)
        add('short_later_chunk'+str(size),HEADER+'\textra\n'+row+'\tx\n'+second+'\n'+third+'\ty\n',chunk=size)
        add('extra_first_row_of_chunk'+str(size),HEADER+'\n'+(row+'\n')*size+second+'\textra\n'+third+'\n',chunk=size)
    for name,value in [('float','20.9'),('zero','0'),('negative','-1'),('missing','NA'),('text','twenty'),('boolean','True'),('int32_boundary','2147483648'),('int64_overflow','9223372036854775808'),('infinity','inf')]:
        add('pos_'+name,HEADER+'\n'+row.replace('\t20\t','\t'+value+'\t')+'\n')
    for name,value in [('missing','NA'),('integer','123'),('boolean','False'),('dot','.'),('comma','A,G'),('symbolic','<DEL>')]:
        add('ref_'+name,HEADER+'\n'+row.replace('\tA\t','\t'+value+'\t')+'\n')
        add('alt_'+name,HEADER+'\n'+row.replace('\tC\t','\t'+value+'\t')+'\n')
    for name,value in [('bad_chrom','chrZ'),('numeric_chrom','1'),('leading_zero','chr01'),('mitochondrial','chrM')]:
        add(name,HEADER+'\n'+row.replace('chr1\t',value+'\t',1)+'\n')
    for name,value in [('zero','chr1:0'),('negative','chr1:-1'),('spaces','chr1: 2_1 '),('unicode','chr1:２１'),('many_colons','chr1:21:x'),('missing','NA'),('number','123'),('empty_chrom',':21'),('negative_large','chr1:-99999999999999999999999999999999999999999')]:
        add('site_'+name,HEADER+'\n'+row.replace('chr1:21',value)+'\n')
    for col in HEADER.split('\t'):
        names=HEADER.split('\t');values=row.split('\t');i=names.index(col);names.pop(i);values.pop(i)
        add('missing_column_'+col.replace('#',''), '\t'.join(names)+'\n'+'\t'.join(values)+'\n')
    for name,values in [('integer',['007','+9']),('float',['0.30000000000000004','1e-7']),('bool',['true','FALSE']),('missing',['NA','x']),('quoted',['"a\tb"','"line\nline"']),('quote',['"a""b"','quote"mid']),('unsigned',['9223372036854775808','18446744073709551615']),('unsigned_missing',['9223372036854775808','NA']),('overflow',['18446744073709551616','NA'])]:
        add('extra_'+name,HEADER+'\textra\n'+row+'\t'+values[0]+'\n'+second+'\t'+values[1]+'\n')
    add('infer_per_chunk',HEADER+'\textra\n'+row+'\t007\n'+second+'\ttext\n',chunk=1)
    add('malformed_quote',HEADER+'\textra\n'+row+'\t"unclosed\n')
    add('nul_extra',HEADER+'\textra\n'+row+'\ta\0b\n')
    add('invalid_utf8_read_ahead',raw_bytes=(HEADER+'\n'+row+'\n').encode()+b'\xff\n')
    add('later_failed_chunk',HEADER+'\n'+row+'\n'+second.replace('chr1:31','broken')+'\n',chunk=1)
    add('unmapped_symbolic',HEADER+'\n'+row.replace('\t20\t','\t1200\t').replace('\tC\t','\t<DEL>\t')+'\n')
    add('gzip_output',output_name='output.tsv.gz')
    add('gzip_three_chunks',HEADER+'\n'+row+'\n'+second+'\n'+third+'\n',chunk=1,output_name='output.tsv.gz')
    add('gzip_header_only',HEADER+'\n',output_name='output.tsv.gz')
    add('gzip_later_failed_chunk',HEADER+'\n'+row+'\n'+second.replace('chr1:31','broken')+'\n',chunk=1,output_name='output.tsv.gz')
    add('plain_existing_mode',output_mode=0o640)
    add('gzip_existing_mode',output_name='output.tsv.gz',output_mode=0o640)
    for suffix in ('.bz2','.xz','.zip','.tar.gz'):
        add('unsupported_output'+suffix,output_name='output.tsv'+suffix,unsupported=True)
    for value in ('0','-1','1.5','bad','1_0','１２'):
        add('chunk_option_'+value,chunk=value)
    add('missing_input',missing_input=True);add('temp_is_file',temp_is_file=True)
    add('output_is_directory',output_is_directory=True);add('input_is_output',input_is_output=True)
    add('missing_chain',missing_chain=True);add('missing_reference',missing_reference=True)
    add('help',argv=['--help']);add('missing_required',argv=[])
    add('unknown_option',extra_args=['--unknown'])
    add('invalid_type_then_help',extra_args=['--chunk_size','bad','--help'])
    return cases

def fixtures(args):
    import pysam
    work=args.work_root;work.mkdir(exist_ok=False,parents=True)
    resources=work/'resources';resources.mkdir()
    reference=resources/'reference with space.fa'
    reference.write_text('>chr1\n'+'A'*3000+'\n>chr2\n'+'C'*3000+'\n')
    pysam.faidx(str(reference))
    chain=resources/'chain with quote\'.chain'
    chain.write_text('chain 1 chr1 2000 + 0 1000 chr1 3000 + 100 1100 1\n1000\n\n'
                     'chain 1 chr2 2000 + 0 1000 chr2 3000 - 100 1100 2\n1000\n\n')
    crossmap=shutil.which('CrossMap');assert crossmap
    launch=work/'bin';launch.mkdir()
    wrapper=launch/'CrossMap'
    wrapper.write_text('#!'+sys.executable+'\nimport json,os,sys\n'
        'with open(os.environ["SPLICEVAULT_COMMAND_TRACE"],"a") as f: f.write(json.dumps(sys.argv[1:])+"\\n")\n'
        'os.execv('+repr(crossmap)+',['+repr(crossmap)+']+sys.argv[1:])\n')
    wrapper.chmod(0o755)
    env=dict(os.environ,PATH=str(launch)+':'+os.environ['PATH'],POLARS_MAX_THREADS='1')
    results={};started=time.monotonic()
    selected=fixture_cases()
    if args.case:
        names=set(args.case);selected=[case for case in selected if case['name'] in names]
        assert names=={case['name'] for case in selected},'Unknown fixture name'
    else:
        selected=[case for case in selected if not case.get('unsupported')]
    for case in selected:
        name=case['name'];folder=work/name;folder.mkdir()
        arms={}
        for label in ('python','native'):
            arm=folder/label;arm.mkdir();source=arm/'input.tsv'
            if 'raw_bytes' in case:source.write_bytes(case['raw_bytes'])
            else:source.write_text(case['text'],newline='')
            output=arm/case.get('output_name','output.tsv');temp=arm/'tmp'
            if case.get('missing_input'):source.unlink()
            if case.get('temp_is_file'):temp.write_text('retained')
            if case.get('output_is_directory'):output.mkdir()
            if case.get('input_is_output'):output=source
            if not case.get('input_is_output') and not case.get('output_is_directory'):output.write_bytes(b'previous-output\n')
            if 'output_mode' in case:output.chmod(case['output_mode'])
            command=[sys.executable,str(args.python_script)] if label=='python' else [str(args.binary),'splicevault-convert']
            options=['--input_tsv',str(source),'--output_tsv',str(output),'--chain_file',str(chain if not case.get('missing_chain') else resources/'missing.chain'),
                     '--reference_fasta',str(reference if not case.get('missing_reference') else resources/'missing.fa'),
                     '--temp_dir',str(temp),'--chunk_size',str(case['chunk'])]
            options=case.get('argv',options)+case.get('extra_args',[])
            trace=arm/'commands.jsonl';this_env=dict(env,SPLICEVAULT_COMMAND_TRACE=str(trace))
            with (arm/'run.log').open('wb') as log:
                completed=subprocess.run(command+options,env=this_env,stdout=log,stderr=subprocess.STDOUT,timeout=60)
            payload=output.read_bytes() if output.is_file() else None
            if payload is not None and payload[:2]==b'\x1f\x8b':payload=gzip.decompress(payload)
            commands=[json.loads(line) for line in trace.read_text().splitlines()] if trace.exists() else []
            normalized=normalize(json.dumps(commands,ensure_ascii=False).encode(),arm)
            sidecars=[]
            if temp.is_dir():
                for path in sorted(temp.iterdir()):
                    sidecars.append((UUID.sub('<UUID>',path.name),sha(normalize(path.read_bytes(),arm)) if path.is_file() else 'directory'))
            arms[label]=dict(returncode=completed.returncode,output_sha256=sha(payload) if payload is not None else None,
                output_bytes=len(payload) if payload is not None else None,command_sha256=sha(normalized),
                output_mode=stat.S_IMODE(output.stat().st_mode) if output.is_file() else None,
                calls=len(commands),sidecars=sorted(sidecars),temp_kind='directory' if temp.is_dir() else 'file' if temp.is_file() else 'absent')
        equal=arms['python']==arms['native']
        results[name]=dict(passed=equal,scope='unsupported compression mismatch' if case.get('unsupported') else 'supported table route',**arms)
    report=dict(kind='native-splicevault-mature-tool-fixtures-v1',passed=all(v['passed'] for v in results.values()),
        cases=results,case_count=len(results),failed=[n for n,v in results.items() if not v['passed']],
        binary=dict(path=str(args.binary),sha256=sha(args.binary.read_bytes())),
        python_script=dict(path=str(args.python_script),sha256=sha(args.python_script.read_bytes())),
        CrossMap=crossmap,wall_seconds=time.monotonic()-started)
    args.report.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='cases'},sort_keys=True))
    return 0 if report['passed'] else 1

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    procedures=parser.add_subparsers(dest='procedure',required=True)
    p=procedures.add_parser('fixtures')
    for name in ('binary','python-script','work-root','report'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--case',action='append',help='Run only named cases; unsupported formats remain explicit mismatches')
    p.set_defaults(run=fixtures)
    args=parser.parse_args();raise SystemExit(args.run(args))

if __name__=='__main__':main()
