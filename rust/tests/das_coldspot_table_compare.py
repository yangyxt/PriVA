#!/usr/bin/env python3
"""Compare the DAS coldspot TSV prefix with its unchanged Python statements.

compare --binary BINARY --source SORTED_REFERENCE --work-root NEW_DIRECTORY
oracle --source SORTED_REFERENCE --input TABLE --output JSONL
real --binary BINARY --source SORTED_REFERENCE --input TABLE --work-root NEW_DIRECTORY

Only parsed coldspot rows are produced. Interval filtering and complete DAS
publication are separate gates. Every mismatch and resource-limited oracle is
retained; a killed or timed-out Python process is never a matched exception.
"""
import argparse
import ast
import bz2
import gzip
import hashlib
import io
import json
import lzma
import os
from pathlib import Path
import re
import resource
import struct
import subprocess
import sys
import tarfile
import time
import zipfile

PIN='1a2b4efc312c41c6fb20afa3c1d90566da7cd0772d2f261c2f7f09d094fae834'
HEADER='type\tgene\taa_start_pos\taa_end_pos\n'


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def scalar(value):
    if isinstance(value,str):return dict(kind='str',value=value)
    if isinstance(value,bool):return dict(kind='bool',value=value)
    if isinstance(value,int):return dict(kind='int',value=str(value))
    if isinstance(value,float):return dict(kind='float',bits=struct.pack('>d',value).hex())
    raise TypeError(type(value).__name__)


def oracle(args):
    import pandas as pd
    if digest(args.source)!=PIN:raise RuntimeError('different producing source')
    tree=ast.parse(args.source.read_text())
    fn=next(node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name=='filter_coldspot_overlap')
    body=[]
    for node in fn.body:
        if isinstance(node,ast.Expr) and isinstance(node.value,ast.Call):break
        body.append(node)
    if not isinstance(body[-1],ast.For):raise RuntimeError('coldspot prefix changed')
    result=[]
    class Rows(list):
        def __init__(self,key):super().__init__();self.key=key
        def append(self,pair):
            result.append(dict(gene=scalar(self.key),start=str(pair[0]),end=str(pair[1])))
            super().append(pair)
    class Mapping(dict):
        def __missing__(self,key):
            value=Rows(key);self[key]=value;return value
    scope=dict(pd=pd,hcseeker_spots_tsv=args.input,defaultdict=lambda factory:Mapping())
    exec(compile(ast.fix_missing_locations(ast.Module(body=body,type_ignores=[])),str(args.source),'exec'),scope)
    with args.output.open('x') as out:
        for row in result:out.write(json.dumps(row,ensure_ascii=False,separators=(',',':'))+'\n')


def cases():
    normal=HEADER+'coldspot\tG\t1\t10\nhotspot\tH\t2\t3\ncoldspot\tG\t8\t4\n'
    result=[]
    def add(name,text=normal,**kw):result.append(dict(name=name,payload=text.encode(),**kw))
    add('normal');add('empty','');add('blank',' \n\r\n');add('header_only',HEADER)
    add('duplicates',normal+normal.split('\n',1)[1]);add('unicode',normal.replace('G','基因'))
    add('bom','\ufeff'+normal);add('crlf',normal.replace('\n','\r\n'));add('bare_cr',normal.replace('\n','\r'))
    add('spaces_blank',HEADER+'  \n\ncoldspot\tG\t1\t10\n  \n')
    add('bare_cr_spaced',normal.replace('coldspot',' coldspot').replace('\n','\r'))
    for count in [1,2,3,8]:
        for ending in ['\n','\r','\r\n']:
            add('leading_space_'+str(count)+'_'+ending.encode().hex(),
                HEADER.replace('\n',ending)+(' coldspot\tG\t1\t10'+ending)*count)
    add('mixed_cr_space',HEADER+'coldspot\tG\t1\t10\r coldspot\tG\t1\t10\r')
    add('bom_cr_space','\ufeff'+normal.replace('coldspot',' coldspot').replace('\n','\r'))
    add('quote',HEADER+'coldspot\t"G\tA\nB"\t"1"\t"10"\n')
    add('nul',HEADER+'coldspot\tG\x00ignored\t1\x00bad\t10\n')
    add('malformed_quote',HEADER+'coldspot\t"G\t1\t10\n')
    add('trailing_tabs',normal.replace('\n','\t\n'))
    add('implicit_index',HEADER+'i\tcoldspot\tG\t1\t10\ni\tcoldspot\tG\t2\t9\n')
    add('implicit_multi_index',HEADER+'i\tj\tcoldspot\tG\t1\t10\n')
    add('extra_later',normal+'coldspot\tG\t1\t10\textra\n')
    add('short_row',HEADER+'coldspot\tG\t1\n')
    add('short_ignored',HEADER+'hotspot\tG\t1\n')
    add('missing_input',missing=True);add('directory',directory=True)
    result.append(dict(name='invalid_utf8',payload=normal.encode()+b'\xff\n'))
    for field in ['type','gene','aa_start_pos','aa_end_pos']:
        names=HEADER.strip().split('\t');values=['coldspot','G','1','10'];i=names.index(field)
        names.pop(i);values.pop(i)
        for selected in [False,True]:
            add('missing_'+field+'_'+str(selected),'\t'.join(names)+'\n'+'\t'.join(values).replace('coldspot','coldspot' if selected else 'hotspot')+'\n')
    add('duplicate_names','type\tgene\tgene\taa_start_pos\taa_end_pos\n'+'coldspot\tG\tOTHER\t1\t10\n')
    for name,value in [('int','007'),('negative','-7'),('float','1.25'),('true','true'),('false','False'),
                       ('missing','NA'),('empty',''),('nan','NaN'),('infinity','inf'),('uint','18446744073709551615')]:
        add('gene_'+name,HEADER+'coldspot\t'+value+'\t1\t10\n')
    numbers=['007','+007','-0','2.9','-2.9','1e20','1e100','True','false','NA','nan','inf','-inf',
             ' 2_0 ','１２','-１２','\u00a0+١_٢\u00a0','1__2','2.5x','9223372036854775808',
             '18446744073709551615','18446744073709551616','9'*4300,'9'*4301]
    for i,value in enumerate(numbers):
        add('start_'+str(i),HEADER+'coldspot\tG\t'+value+'\t10\n')
        add('end_'+str(i),HEADER+'coldspot\tG\t1\t'+value+'\n')
    for extension in ['gz','bz2','xz','zst','zip','tar','tar.gz','tar.bz2','tar.xz','GZ','ZIP']:
        add('codec_'+extension,codec=extension)
    add('gzip_members',codec='gz',members=2)
    for ending in [b'\0',b'\0'*4096,b'x',b'xx',b'\x1f',b'\x1f\x8b']:
        add('gzip_trailing_'+str(len(ending))+'_'+ending[:2].hex(),codec='gz',trailing=ending)
    add('gzip_padded_members',codec='gz',members=2,padding=b'\0'*17)
    for codec in ['gz','bz2','xz','zst']:
        add('truncated_'+codec,codec=codec,truncate=5)
        add('invalid_'+codec,codec=codec,invalid=True)
        for keep in [0,1,2,3,4,5,8,12,16,24,32,48,64]:
            add('prefix_'+codec+'_'+str(keep),codec=codec,keep=keep)
        for missing in [1,2,3,4,6,7,8,12,16]:
            add('suffix_'+codec+'_'+str(missing),codec=codec,truncate=missing)
    add('zstd_incomplete_large',normal*10000,codec='zst',truncate=1)
    add('gzip_crc',codec='gz',crc=True)
    for codec in ['zip','tar']:
        for members in [0,2]:add(codec+'_members_'+str(members),codec=codec,members=members)
    # 32 columns give a 16,384-row C inference block. Exercise both conversion
    # and first-buffered-row surplus behavior without creating huge fixtures.
    wide=HEADER.rstrip()+'\t'+'\t'.join('x'+str(i) for i in range(28))+'\n'
    def row(gene='G',start='1',kind='coldspot'):
        return '\t'.join([kind,gene,start,'10']+['x']*28)+'\n'
    for name,left,right in [('text_int','007','G'),('text_uint','18446744073709551615','G'),
                            ('mixed_uint_signed','18446744073709551615','-1'),('bool_text','true','G')]:
        add('block_'+name,wide+row(gene=left)*16384+row(gene=right))
    add('block_coordinate',wide+row(start='1.9')*16384+row(start=' 2_0 '))
    add('block_type',wide+row(kind='1')*16384+row())
    for at in [16383,16384,16385,32768]:
        add('block_extra_'+str(at),wide+row()*at+row().rstrip('\n')+'\textra\n'+row())
    return result


def write_input(path,case):
    payload=case['payload'];codec=case.get('codec','').lower();members=case.get('members',1)
    if codec=='gz':
        payload=(gzip.compress(payload[:len(payload)//2],mtime=0)+case.get('padding',b'')+gzip.compress(payload[len(payload)//2:],mtime=0)) if members==2 else gzip.compress(payload,mtime=0)
    elif codec=='bz2':payload=bz2.compress(payload)
    elif codec=='xz':payload=lzma.compress(payload)
    elif codec=='zst':
        import zstandard
        payload=zstandard.ZstdCompressor().compress(payload)
    elif codec=='zip':
        out=io.BytesIO()
        with zipfile.ZipFile(out,'w',compression=zipfile.ZIP_DEFLATED) as archive:
            for i in range(members):archive.writestr(str(i)+'.tsv',payload)
        payload=out.getvalue()
    elif codec.startswith('tar'):
        out=io.BytesIO();mode={'tar':'w','tar.gz':'w:gz','tar.bz2':'w:bz2','tar.xz':'w:xz'}[codec]
        with tarfile.open(fileobj=out,mode=mode) as archive:
            for i in range(members):
                member=tarfile.TarInfo(str(i)+'.tsv');member.size=len(payload);archive.addfile(member,io.BytesIO(payload))
        payload=out.getvalue()
    if case.get('truncate'):payload=payload[:-case['truncate']]
    if 'keep' in case:payload=payload[:case['keep']]
    payload+=case.get('trailing',b'')
    if case.get('invalid'):payload=b'not a compressed file\n'
    if case.get('crc'):
        payload=bytearray(payload);payload[-8]^=1;payload=bytes(payload)
    if case.get('directory'):path.mkdir()
    elif not case.get('missing'):path.write_bytes(payload)


def run(command,folder,label,output=None):
    error=folder/(label+'.stderr.log');outpath=output or folder/(label+'.stdout.log')
    with outpath.open('xb') as out,error.open('xb') as err:
        def bounded():resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
        process=subprocess.Popen(list(map(str,command)),stdout=out,stderr=err,preexec_fn=bounded)
        limited=False
        try:code=process.wait(timeout=20)
        except subprocess.TimeoutExpired:
            limited=True;process.kill();code=process.wait()
    lines=error.read_text(errors='replace').splitlines()
    # Native runtime errors have the program prefix; it is not an exception
    # category. Retain the actual category, including zlib's lowercase error.
    category=next((m[1].rsplit('.',1)[-1] for line in reversed(lines)
                   if (m:=re.match(r'([A-Za-z_.]+):',line.removeprefix('priva: ')))),None)
    allocation_failure=any('out of memory' in line.lower() for line in lines)
    return dict(exit=code,category=category,
                resource_limited=limited or code<0 or category=='MemoryError' or allocation_failure)


def one(args,folder,path):
    expected=folder/'python.jsonl';actual=folder/'native.jsonl'
    py=run([sys.executable,__file__,'oracle','--source',args.source,'--input',path,'--output',expected],folder,'python')
    native=run([args.binary,'das-domains','coldspot-table','--spots',path],folder,'native',actual)
    mismatches=[];rows=0
    if py['resource_limited'] or native['resource_limited']:
        return dict(qualified=False,passed=False,python=py,native=native,scope='resource-dependent parser outcome')
    if py['exit']!=native['exit']:mismatches.append('exit code differs')
    elif py['exit']:
        if py['category'] is None or native['category'] is None:mismatches.append('unrecognized exception category')
        if py['category']!=native['category']:mismatches.append('exception category differs')
        if actual.stat().st_size:mismatches.append('partial output after failure')
    else:
        import itertools
        with expected.open() as left,actual.open() as right,(folder/'mismatches.jsonl').open('x') as ledger:
            for ordinal,(a,b) in enumerate(itertools.zip_longest(left,right)):
                rows+=1
                if a is None or b is None or json.loads(a)!=json.loads(b):
                    ledger.write(json.dumps(dict(ordinal=ordinal,python=a,native=b))+'\n');mismatches.append(ordinal)
    return dict(qualified=True,passed=not mismatches,python=py,native=native,rows=rows,mismatches=len(mismatches),
                failure_reasons=[v for v in mismatches if isinstance(v,str)],
                python_sha256=digest(expected) if expected.exists() else None,native_sha256=digest(actual))


def compare(args):
    root=args.work_root;root.mkdir(parents=True,exist_ok=False);start=time.monotonic();records=[]
    selected=cases()
    if args.case:
        wanted=set(args.case);selected=[v for v in selected if v['name'] in wanted]
        assert {v['name'] for v in selected}==wanted
    for case in selected:
        folder=root/case['name'];folder.mkdir();path=folder/('spots.tsv'+('.'+case['codec'] if case.get('codec') else ''))
        write_input(path,case);record=dict(case=case['name'],**one(args,folder,path));records.append(record)
        print(json.dumps({k:v for k,v in record.items() if k not in ['python_sha256','native_sha256']}),flush=True)
    result=dict(passed=all(v['passed'] for v in records),cases=len(records),matched=sum(v['passed'] for v in records),
        qualified=sum(v['qualified'] for v in records),resource_boundaries=[v for v in records if not v['qualified']],
        mismatches=[v for v in records if v['qualified'] and not v['passed']],results=records,
        binary_sha256=digest(args.binary),source_sha256=digest(args.source),harness_sha256=digest(__file__),
        full_das_qualified=False,wall_seconds=time.monotonic()-start)
    (root/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    if not result['passed']:raise SystemExit(1)


def real(args):
    root=args.work_root;root.mkdir(parents=True,exist_ok=False)
    before=digest(args.input);result=one(args,root,args.input)
    if before!=digest(args.input):raise RuntimeError('input changed')
    result.update(input=dict(path=str(args.input),bytes=args.input.stat().st_size,sha256=before),
                  binary_sha256=digest(args.binary),source_sha256=digest(args.source),harness_sha256=digest(__file__),full_das_qualified=False)
    (root/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='input'}),flush=True)
    if not result['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    call=commands.add_parser('oracle')
    for name in ['source','input','output']:call.add_argument('--'+name,type=Path,required=True)
    call.set_defaults(function=oracle)
    call=commands.add_parser('compare')
    for name in ['binary','source','work-root']:call.add_argument('--'+name,type=Path,required=True)
    call.add_argument('--case',action='append');call.set_defaults(function=compare)
    call=commands.add_parser('real')
    for name in ['binary','source','input','work-root']:call.add_argument('--'+name,type=Path,required=True)
    call.set_defaults(function=real)
    args=parser.parse_args();args.function(args)


if __name__=='__main__':main()
