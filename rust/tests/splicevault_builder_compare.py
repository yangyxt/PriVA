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
import bz2
import gzip
import hashlib
import io
import json
import lzma
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import time
import zipfile
import zlib

HEADER='#chrom\tpos\tref\talt\tsplicevault_site_loc'
UUID=re.compile(r'[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}')

def sha(data):return hashlib.sha256(data).hexdigest()

def normalize(data,arm):
    data=UUID.sub('<UUID>',data.decode('utf-8')).replace(str(arm),'<ARM>')
    return re.sub(r'(?m)^##liftOverDate=.*$','##liftOverDate=<DATE>',data).encode()

def output_signature(payload, started, ended):
    """Compare content, stream boundaries and archive metadata, not codec bytes.

    ZIP timestamps are checked against the actual producing invocation. ZIP
    deflate bytes and gzip timestamps can differ without any data difference.
    Logical members include duplicates in order; reading by filename alone
    would incorrectly discard earlier duplicate ZIP members.
    """
    if payload is None:return None
    if payload.startswith(b'PK'):
        entries=[]
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            for member in archive.infolist():
                value=archive.read(member)
                stamp=time.mktime((*member.date_time,0,0,-1))
                entries.append(dict(name=member.filename,bytes=len(value),sha256=sha(value),
                    method=member.compress_type,flag_bits=member.flag_bits,
                    external_attr=member.external_attr,internal_attr=member.internal_attr,
                    create_system=member.create_system,create_version=member.create_version,
                    extract_version=member.extract_version,extra_hex=member.extra.hex(),
                    comment_hex=member.comment.hex(),timestamp_in_run=started is None or started-2<=stamp<=ended+2))
            comment=archive.comment.hex()
        return dict(format='zip',members=entries,comment_hex=comment,
                    timestamps_valid=all(x['timestamp_in_run'] for x in entries))
    streams=[];decoded=payload;codec='plain'
    if payload.startswith((b'BZh',b'\xfd7zXZ\0',b'\x1f\x8b',b'\x28\xb5\x2f\xfd')):
        codec=('bz2' if payload.startswith(b'BZh') else 'xz' if payload.startswith(b'\xfd7zXZ\0')
               else 'gzip' if payload.startswith(b'\x1f\x8b') else 'zstd')
        remaining=payload;pieces=[]
        while remaining:
            if codec=='bz2':decoder=bz2.BZ2Decompressor()
            elif codec=='xz':decoder=lzma.LZMADecompressor()
            elif codec=='gzip':decoder=zlib.decompressobj(31)
            else:
                import zstandard
                decoder=zstandard.ZstdDecompressor().decompressobj()
            value=decoder.decompress(remaining)
            assert decoder.eof,'Incomplete '+codec+' stream'
            streams.append(dict(bytes=len(value),sha256=sha(value)))
            pieces.append(value);remaining=decoder.unused_data
        decoded=b''.join(pieces)
    members=[]
    try:
        with tarfile.open(fileobj=io.BytesIO(decoded),mode='r:') as archive:
            for member in archive.getmembers():
                value=archive.extractfile(member).read()
                members.append(dict(name=member.name,bytes=len(value),sha256=sha(value),
                    mode=member.mode,uid=member.uid,gid=member.gid,mtime=member.mtime,
                    type=member.type.decode(),uname=member.uname,gname=member.gname,
                    pax_headers=member.pax_headers))
        codec='tar+'+codec
    except tarfile.ReadError:pass
    return dict(format=codec,bytes=len(decoded),sha256=sha(decoded),streams=streams,members=members)

def publication_signature(arm, output, started, ended):
    records=[json.loads(line) for line in (arm/'publication.jsonl').read_text().splitlines()]
    states=[];last=None
    for record in records+[dict(file=output.name,mtime_ns=output.stat().st_mtime_ns)]:
        payload=(arm/record['file']).read_bytes()
        current=dict(sha256=sha(payload),mtime_ns=record['mtime_ns'])
        states.append(dict(output=output_signature(payload,started,ended),
            bytes_unchanged=None if last is None else current['sha256']==last['sha256'],
            mtime_unchanged=None if last is None else current['mtime_ns']==last['mtime_ns']))
        last=current
    return states

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
    # Retain the original case names so saved Python oracles remain reusable.
    for suffix in ('.bz2','.xz','.zip','.tar.gz'):
        add('unsupported_output'+suffix,output_name='output.tsv'+suffix,compression=True)
    for suffix in ('bz2','xz','zst','zip','tar','tar.gz','tar.bz2','tar.xz'):
        for profile,text,chunk in (
            ('one',HEADER+'\n'+row+'\n',2),
            ('three',HEADER+'\n'+row+'\n'+second+'\n'+third+'\n',1),
            ('empty_middle',HEADER+'\n'+row+'\n'+second.replace('\t30\t','\t1200\t')+'\n'+third+'\n',1),
            ('empty_last',HEADER+'\n'+row+'\n'+second.replace('\t30\t','\t1200\t')+'\n',1),
            ('later_failure',HEADER+'\n'+row+'\n'+second.replace('chr1:31','broken')+'\n',1),
            ('header_only',HEADER+'\n',2),
            ('all_unmapped',HEADER+'\n'+row.replace('\t20\t','\t1200\t')+'\n',2)):
            if profile=='one' and suffix in ('bz2','xz','zip','tar.gz'):continue
            add('compression_'+suffix+'_'+profile,text,chunk=chunk,output_name='output.tsv.'+suffix,
                output_mode=0o640,compression=True)
    for basename in ('output.tsv.BZ2','output.tsv.XZ','output.tsv.ZST','output.tsv.ZIP',
                     'output.tsv.TAR','output.tsv.TAR.GZ','output.tsv.tar.GZ','output.tsv.TAR.gz',
                     'output.tsv.TAR.BZ2','output.tsv.TAR.XZ','é.tsv.zip','é.tsv.tar.gz',
                     'a'*110+'.tar.gz','.zip','.tar'):
        add('compression_name_'+basename,HEADER+'\n'+row+'\n'+second+'\n',chunk=1,
            output_name=basename,compression=True)
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
    wrapper.write_text('#!'+sys.executable+'\nimport json,os,sys\nfrom pathlib import Path\n'
        'if os.environ.get("SPLICEVAULT_PUBLICATION_TRACE"):\n'
        '    output=Path(os.environ["SPLICEVAULT_PUBLICATION_TRACE"])\n'
        '    log=output.parent/"publication.jsonl"\n'
        '    index=len(log.read_text().splitlines()) if log.exists() else 0\n'
        '    snapshot=output.parent/("publication_%03d.bin"%index)\n'
        '    state=output.stat();snapshot.write_bytes(output.read_bytes())\n'
        '    with log.open("a") as f:f.write(json.dumps(dict(file=snapshot.name,mtime_ns=state.st_mtime_ns))+"\\n")\n'
        'with open(os.environ["SPLICEVAULT_COMMAND_TRACE"],"a") as f: f.write(json.dumps(sys.argv[1:])+"\\n")\n'
        'os.execv('+repr(crossmap)+',['+repr(crossmap)+']+sys.argv[1:])\n')
    wrapper.chmod(0o755)
    env=dict(os.environ,PATH=str(launch)+':'+os.environ['PATH'],POLARS_MAX_THREADS='1')
    results={};started=time.monotonic()
    reused={}
    if args.reuse_python_results:
        assert args.reuse_python_root,'A saved oracle also requires its unchanged fixture root'
        saved=json.loads(args.reuse_python_results.read_text())
        assert saved['python_script']['sha256']==sha(args.python_script.read_bytes()),'Different producing Python source'
        for resource in resources.iterdir():
            assert resource.read_bytes()==(args.reuse_python_root/'resources'/resource.name).read_bytes(),'Different fixture resource'
        reused=saved['cases']
    selected=fixture_cases()
    if args.case:
        names=set(args.case);selected=[case for case in selected if case['name'] in names]
        assert names=={case['name'] for case in selected},'Unknown fixture name'
    elif args.compression_only:
        selected=[case for case in selected if case.get('compression')]
    else:
        selected=[case for case in selected if not case.get('compression')]
    for case in selected:
        name=case['name'];folder=work/name;folder.mkdir()
        arms={};encoded={}
        for label in ('python','native'):
            reuse=label=='python' and name in reused
            arm=args.reuse_python_root/name/'python' if reuse else folder/label
            source=arm/'input.tsv'
            if reuse:
                assert source.read_bytes()==case.get('raw_bytes',case['text'].encode()),'Different saved fixture input'
            else:
                arm.mkdir()
                if 'raw_bytes' in case:source.write_bytes(case['raw_bytes'])
                else:source.write_text(case['text'],newline='')
            output=arm/case.get('output_name','output.tsv');temp=arm/'tmp'
            if not reuse:
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
            if args.publication_trace and name not in reused:this_env['SPLICEVAULT_PUBLICATION_TRACE']=str(output)
            if reuse:
                returncode=reused[name]['python']['returncode'];start=end=None
            else:
                with (arm/'run.log').open('wb') as log:
                    start=time.time()
                    completed=subprocess.run(command+options,env=this_env,stdout=log,stderr=subprocess.STDOUT,timeout=60)
                    end=time.time();returncode=completed.returncode
            payload=output.read_bytes() if output.is_file() else None
            if reuse:
                old_payload=gzip.decompress(payload) if payload and payload.startswith(b'\x1f\x8b') else payload
                assert sha(old_payload)==reused[name]['python']['output_sha256'],'Saved oracle output changed'
            signature=output_signature(payload,start,end)
            encoded[label]=dict(sha256=sha(payload) if payload is not None else None,
                bytes=len(payload) if payload is not None else None,started=start,ended=end,reused=reuse,
                observed_storage_mode=stat.S_IMODE(output.stat().st_mode) if output.is_file() else None)
            commands=[json.loads(line) for line in trace.read_text().splitlines()] if trace.exists() else []
            resource_root=args.reuse_python_root/'resources' if reuse else resources
            normalized=normalize(json.dumps(commands,ensure_ascii=False).encode(),arm).replace(str(resource_root).encode(),b'<RESOURCES>')
            sidecars=[]
            if temp.is_dir():
                for path in sorted(temp.iterdir()):
                    sidecars.append((UUID.sub('<UUID>',path.name),sha(normalize(path.read_bytes(),arm).replace(str(resource_root).encode(),b'<RESOURCES>')) if path.is_file() else 'directory'))
            arms[label]=dict(returncode=returncode,output=signature,command_sha256=sha(normalized),
                # An archive store may make a saved oracle read-only. Its
                # immutable producing receipt records publication permissions.
                output_mode=reused[name]['python']['output_mode'] if reuse else encoded[label]['observed_storage_mode'],
                calls=len(commands),sidecars=sorted(sidecars),temp_kind='directory' if temp.is_dir() else 'file' if temp.is_file() else 'absent')
            if args.publication_trace and name not in reused and commands:
                arms[label]['publication']=publication_signature(arm,output,start,end)
        equal=arms['python']==arms['native'] and all(a['output'] is None or a['output'].get('timestamps_valid',True) for a in arms.values())
        results[name]=dict(passed=equal,scope='compression/append route' if case.get('compression') else 'supported table route',
                           encoded=encoded,compression_bytes_equal=encoded['python']['sha256']==encoded['native']['sha256'],**arms)
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
    p.add_argument('--case',action='append',help='Run only named cases')
    p.add_argument('--compression-only',action='store_true',help='Run inferred compression and archive append cases')
    p.add_argument('--reuse-python-root',type=Path,help='Read-only fixture root with saved Python outputs')
    p.add_argument('--reuse-python-results',type=Path,help='Original report binding saved Python outputs and producing source')
    p.add_argument('--publication-trace',action='store_true',help='Compare output snapshots and relative mtimes before each real CrossMap call')
    p.set_defaults(run=fixtures)
    args=parser.parse_args();raise SystemExit(args.run(args))

if __name__=='__main__':main()
