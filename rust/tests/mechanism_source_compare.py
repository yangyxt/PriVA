#!/usr/bin/env python3
"""Differential acquisition-source checks; no downloads or installed writes.

check --producer SCRIPT --binary PRIVA --work-root NEW_DIRECTORY
oracle --producer SCRIPT --exception JSON -- original arguments...
Compare exact exits, categories, diagnostics, payload/receipt bytes and modes.
Only unpredictable temporary sibling names and the native diagnostic envelope
are normalized. Unsupported JSON representations remain measured mismatches.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import time
sys.dont_write_bytecode=True
URL='https://selected.invalid/source.tsv'
EREPO='ClinVar Variation Id\tApplied Evidence Codes (Met)\n123\tPS3\n'
DOSAGE='#Gene Symbol\tGenomic Location\tHaploinsufficiency Score\tHaploinsufficiency Description\tHaploinsufficiency Disease ID\tTriplosensitivity Score\nGENE\t1:1-2\t3\tdescription\tMONDO:1\t0\n'
G2P='gene symbol,disease name,molecular mechanism\nGENE,condition,mechanism\n'
XML='<JDBOR><DisorderList><Disorder/></DisorderList></JDBOR>'
BASE={'clingen_erepo':EREPO,'clingen_dosage':DOSAGE,'clingen_dosage_runtime':DOSAGE,'g2p_ddg2p':G2P,'orphadata':XML}

def pin(path):
    path=Path(path)
    with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    return dict(path=str(path),bytes=path.stat().st_size,sha256=digest)

def save(path,value):
    with Path(path).open('x') as out:json.dump(value,out,indent=2,sort_keys=True);out.write('\n')

def oracle(args):
    spec=importlib.util.spec_from_file_location('source_oracle',args.producer)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);last={}
    def observe(frame,event,arg):
        if event=='exception' and frame.f_code.co_filename==str(args.producer):
            category,exception,_=arg
            if category is not SystemExit:last.update(category=category.__name__,message=str(exception))
        return observe
    sys.settrace(observe)
    try:code=module.main(args.arguments)
    except Exception as exc:
        sys.settrace(None);last.update(category=type(exc).__name__,message=str(exc));raise
    finally:
        sys.settrace(None)
        if last:save(args.exception,last)
    raise SystemExit(code)

def cases():
    result={}
    def add(name,payload=b'payload\x00bytes\n',**kw):result[name]=dict(payload=payload.hex() if isinstance(payload,bytes) else payload.encode().hex(),**kw)
    def validate(name,kind,payload,**kw):add(name,payload,command='validate',kind=kind,**kw)
    for kind,payload in BASE.items():
        validate(kind,kind,payload)
        for name,spec in [('empty',dict(payload='')),('missing',dict(input_kind='missing')),('directory',dict(input_kind='directory')),
                          ('denied',dict(input_kind='denied')),('symlink',dict(input_kind='symlink')),('loop',dict(input_kind='loop'))]:
            validate(kind+'_'+name,kind,spec.pop('payload',payload),**spec)
        if kind=='orphadata':continue
        header,row=payload.split('\n',1);sep=',' if kind=='g2p_ddg2p' else '\t'
        for name,text in [('header_only',header+'\n'),('blank_rows',header+'\n\n'+sep*(len(header.split(sep))-1)+'\n'),
            ('short_row',payload+'short\n'),('long_row',payload+row.rstrip('\n')+sep+'extra\n'),
            ('wrong_header','bad'+sep+'header\n'+row),('bom','\ufeff'+payload),('comment','#comment\n'+payload),
            ('crlf',payload.replace('\n','\r\n')),('cr',payload.replace('\n','\r')),('no_final_newline',payload.rstrip('\n')),
            ('blank_header','\n'+payload),('quote_header',sep.join('"'+v+'"' for v in header.split(sep))+'\n'+row),
            ('quote_after_close',header+'\n"123"x'+sep+'X\n'),('unclosed_quote',header+'\n"123'),
            ('field_at_limit',header+'\n'+sep.join(['X']*(len(header.split(sep))-1)+['X'*131072])+'\n'),
            ('field_over_limit',header+'\n'+sep.join(['X']*(len(header.split(sep))-1)+['X'*131073])+'\n'),
            ('nul_value',payload.replace('PS3','P\x00S3').replace('description','d\x00escription').replace('mechanism\n','mech\x00anism\n'))]:
            validate(kind+'_'+name,kind,text)
        for name,raw in [('bad_utf8_early',b'\xff'+payload.encode()),('bad_utf8_after_row',payload.encode()+b'\xff'),
            ('bad_utf8_late',payload.encode()+b'\n'*8192+b'\xff'),('bom_bad_utf8',b'\xef\xbb\xbf'+payload.encode()+b'\xff'),
            ('truncated_utf8',payload.encode()+b'\xe2\x82')]:validate(kind+'_'+name,kind,raw)
    flexible='symbol\thi score\tts score\nGENE\t3\t0\n'
    validate('dosage_aliases','clingen_dosage',flexible)
    validate('dosage_runtime_aliases','clingen_dosage_runtime',flexible)
    validate('dosage_filtered_multiline','clingen_dosage','symbol\thi score\tts score\n"GENE\n#ignored\nstill gene"\t3\t0\n')
    validate('g2p_multiline','g2p_ddg2p','gene symbol,disease name,molecular mechanism\nG,"line1\r\nline2",X\n')
    validate('g2p_escaped_quote','g2p_ddg2p','gene symbol,disease name,molecular mechanism\nG,"a""b",X\n')
    validate('g2p_unquoted_quote','g2p_ddg2p','gene symbol,disease name,molecular mechanism\nG,a"b,X\n')
    validate('g2p_unicode_header','g2p_ddg2p','GENE SYMBOL,Disease Label,Mutation Consequence\nG,D,M\n')
    for value in ['', 'NaN','nan','NA','N/A','NULL','null','None','<NA>','#N/A','#N/A N/A','#NA','-1.#IND','-1.#QNAN','-NaN','-nan','1.#IND','1.#QNAN','n/a',
                  'NAN','Inf','infinity','-Infinity','+nan','0','-1','1.9','1e300','1e309','1e-999','1_000','1__0','_1','1_','١٢.٣','１２e２','²','\u00a0+1\u2003','1\u001c','1x','0x10']:
        label='numeric_'+str(len(result))
        validate(label,'clingen_erepo',EREPO.replace('123\t',value+'\t'))
    # First duplicate header wins at list.index, while membership sees all.
    validate('duplicate_numeric_header','clingen_erepo','ClinVar Variation Id\tClinVar Variation Id\tApplied Evidence Codes (Met)\ninvalid\t123\tPS3\n')
    for name,payload in [('wrong_root','<root><DisorderList/></root>'),('missing_list','<JDBOR/>'),('nested_list','<JDBOR><x><DisorderList/></x></JDBOR>'),
        ('default_namespace','<JDBOR xmlns="urn:x"><DisorderList/></JDBOR>'),('child_namespace','<JDBOR><DisorderList xmlns="urn:x"/></JDBOR>'),
        ('attribute_namespace','<JDBOR xmlns:p="urn:x" p:a="x"><DisorderList/></JDBOR>'),('namespace_then_malformed','<JDBOR xmlns="urn:x"><bad'),
        ('late_namespace','<JDBOR>'+' '*17000+'<DisorderList xmlns="u"/></JDBOR>'),('truncated','<JDBOR><DisorderList>'),
        ('junk_after',XML+'junk'),('mismatched','<JDBOR><DisorderList></JDBOR>'),('bad_entity','<JDBOR>&absent;<DisorderList/></JDBOR>'),
        ('internal_entity','<!DOCTYPE JDBOR [<!ENTITY x "<DisorderList/>">]><JDBOR>&x;</JDBOR>'),
        ('external_entity','<!DOCTYPE JDBOR [<!ENTITY x SYSTEM "file:///not-read.txt">]><JDBOR>&x;<DisorderList/></JDBOR>'),
        ('external_dtd','<!DOCTYPE JDBOR SYSTEM "file:///not-read.dtd">'+XML),('comment_pi','<?xml version="1.0"?>\n<!--comment--><?p x?>'+XML),
        ('wrongcase','<jdbor><DisorderList/></jdbor>'),('lower_list','<JDBOR><disorderlist/></JDBOR>'),('invalid_utf8',XML.encode()+b'\xff')]:
        validate('xml_'+name,'orphadata',payload)
    for encoding in ('utf-8','utf-16','utf-16-le','utf-16-be','iso-8859-1','windows-1252','utf-32','latin-1','unknown'):
        document='<?xml version="1.0" encoding="'+encoding+'"?>'+XML
        validate('xml_encoding_'+encoding,'orphadata',document.encode(encoding if encoding not in ('unknown',) else 'utf-8'))
    for url in [URL,'HTTPS://EXAMPLE.invalid/a','http://user:pass@host:badport/a','https://[::1]/','https://[fe80::1%scope]/','https://[v1.ab]/','https://é.invalid/α','https://h/a?b=c#d',
        'ftp://h/a','//h/a','https:///a','https://h/a b',' https://h/a','https://h/\t','https://[a]/','https://[127.0.0.1]/','https://[vX.a]/','https://[::1','https://h／x/a']:
        add('url_'+str(len(result)),command='record-source',url=url,argument=not url.startswith(('https://selected','HTTPS://','http://user','https://[::1]/','https://[fe80','https://[v1','https://é','https://h/a?')))
    digest=hashlib.sha256(b'payload\x00bytes\n').hexdigest()
    valid=dict(kind='priva_clingen_source_v1',url=URL,sha256=digest)
    receipts={'missing':None,'valid':valid,'extra':dict(valid,extra=2),'wrong_url':dict(valid,url='https://other.invalid'),
        'wrong_kind':dict(valid,kind='other'),'wrong_hash':dict(valid,sha256='0'*64),'upper_hash':dict(valid,sha256=digest.upper()),
        'null':None,'list':[],'number':1,'string':'value','bool':True,
        'download_null':dict(valid,download_sha256=None),'download_nested':dict(valid,download_sha256={'b':[1,True,None],'a':'α'}),
        'download_nonfinite':dict(valid,download_sha256=float('nan')),'download_number':dict(valid,download_sha256=9007199254740993)}
    for command in ('record-source','check-source','metadata'):
        for name,receipt in receipts.items():
            add(command+'_'+name,command=command,receipt=receipt,receipt_present=name!='missing')
        for name,raw in [('malformed','{bad'),('utf8',b'\xff'),('bom',b'\xef\xbb\xbf{}'),('duplicate',json.dumps(valid)[:-1]+',"url":"wrong"}'),
                         ('wide_integer',json.dumps(dict(valid,download_sha256=2**160))),('surrogate',json.dumps(dict(valid,download_sha256='\ud800'))),('deep','['*600+'0'+']'*600)]:
            add(command+'_json_'+name,command=command,receipt_raw=(raw if isinstance(raw,bytes) else raw.encode()).hex(),
                unqualified='Shared JSON representation boundary' if name in ('wide_integer','surrogate','deep') else None)
        for name in ('missing','directory','empty','denied','symlink','loop'):
            add(command+'_input_'+name,command=command,input_kind=name,receipt=valid,receipt_present=True)
        for name in ('directory','denied','symlink','broken_link'):
            add(command+'_sidecar_'+name,command=command,receipt=valid,receipt_present=True,sidecar_kind=name)
    add('record_preserves_download_hash',command='record-source',receipt=dict(valid,download_sha256='9'*64),receipt_present=True)
    add('record_changes_url_resets_download',command='record-source',receipt=dict(valid,url='https://old.invalid',download_sha256='9'*64),receipt_present=True)
    add('metadata_explicit_receipt',command='metadata',receipt=valid,receipt_present=True,explicit=True,input_kind='missing')
    add('check_explicit_receipt',command='check-source',receipt=valid,receipt_present=True,explicit=True)
    add('record_denied_parent',command='record-source',parent_denied=True)
    add('relative_source',command='record-source',relative=True)
    return result

def normalized(text,root):
    text=text.replace(str(root),'$CASE')
    return re.sub(r'(?<=\.source\.json\.)[A-Za-z0-9_]{6,8}(?![A-Za-z0-9_])','RANDOM',text)

def snapshot(folder):
    result={}
    for path in sorted(folder.rglob('*')):
        key=normalized(str(path.relative_to(folder)),folder);info=path.lstat();mode=stat.S_IMODE(info.st_mode)
        if path.is_symlink():result[key]=dict(kind='symlink',target=normalized(os.readlink(path),folder),mode=mode)
        elif path.is_file():
            old_mode=mode
            if not mode&0o400:path.chmod(mode|0o400)
            result[key]=dict(kind='file',bytes=info.st_size,sha256=pin(path)['sha256'],mode=old_mode)
            if not old_mode&0o400:path.chmod(old_mode)
        elif path.is_dir():result[key]=dict(kind='directory',mode=mode)
    return result

def prepare(folder,spec):
    folder.mkdir();path=folder/'payload';kind=spec.get('input_kind','file');payload=bytes.fromhex(spec['payload'])
    if kind=='directory':path.mkdir()
    elif kind=='loop':path.symlink_to('payload')
    elif kind=='symlink':(folder/'target').write_bytes(payload);path.symlink_to('target')
    elif kind!='missing':path.write_bytes(b'' if kind=='empty' else payload)
    if kind=='denied':path.chmod(0)
    sidecar=folder/('explicit.json' if spec.get('explicit') else 'payload.source.json')
    if spec.get('receipt_present') or 'receipt_raw' in spec:
        content=bytes.fromhex(spec['receipt_raw']) if 'receipt_raw' in spec else (json.dumps(spec['receipt'])+'\n').encode()
        sidecar.write_bytes(content)
    sk=spec.get('sidecar_kind')
    if sk=='directory':sidecar.unlink();sidecar.mkdir()
    elif sk=='denied':sidecar.chmod(0)
    elif sk=='symlink':sidecar.rename(folder/'receipt_target');sidecar.symlink_to('receipt_target')
    elif sk=='broken_link':sidecar.unlink();sidecar.symlink_to('absent_receipt')
    if spec.get('parent_denied'):folder.chmod(0o500)
    name='payload' if spec.get('relative') else str(path)
    command=spec['command'];args=[command if command!='metadata' else 'check-source','--input',name]
    if command=='validate':args+=['--kind',spec['kind']]
    else:args+=['--url',spec.get('url',URL)]
    if command=='metadata':args+=['--metadata-only']
    if spec.get('explicit'):args+=['--receipt',str(sidecar)]
    return args

def measure(command,folder,label,cwd):
    started=time.monotonic();stdout=folder/(label+'.stdout');stderr=folder/(label+'.stderr')
    with stdout.open('xb') as out,stderr.open('xb') as err:
        child=subprocess.Popen(list(map(str,command)),cwd=cwd,stdout=out,stderr=err,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONHASHSEED='0'))
        _,status,usage=os.wait4(child.pid,0);child.returncode=os.waitstatus_to_exitcode(status)
    return dict(command=list(map(str,command)),returncode=child.returncode,wall_seconds=time.monotonic()-started,peak_rss_kib=usage.ru_maxrss,stdout=pin(stdout),stderr=pin(stderr))

def outcome(measurement,folder,exception,native=False,argument=False):
    code=measurement['returncode'];stderr=Path(measurement['stderr']['path']).read_text();stdout=Path(measurement['stdout']['path']).read_bytes()
    if code==0:return dict(exit=0,stdout=stdout.hex(),category=None,message=None)
    if argument:return dict(exit=code,stdout=stdout.hex(),category='argument' if code==2 else 'unexpected',message=None)
    if not stderr:return dict(exit=code,stdout=stdout.hex(),category=None,message=None)
    if native:
        match=re.fullmatch(r'priva: ([A-Za-z]+): (.*)\n',stderr,re.S)
        category,message=match.groups() if match else ('unclassified',stderr)
    else:
        value=json.loads(exception.read_text()) if exception.exists() else dict(category='unclassified',message=stderr)
        category,message=value['category'],value['message']
    return dict(exit=code,stdout=stdout.hex(),category=category,message=normalized(message,folder))

def compare_case(root,name,spec,producer,binary):
    folder=root/name;folder.mkdir();measurements={};outcomes={};after={};preserved={};before={}
    for arm in ('python','native'):
        work=folder/arm;args=prepare(work,spec);before[arm]=snapshot(work);exception=folder/(arm+'.exception.json')
        command=[sys.executable,__file__,'oracle','--producer',str(producer),'--exception',str(exception),'--',*args] if arm=='python' else [str(binary),'mechanism-source',*args]
        measurements[arm]=measure(command,folder,arm,work)
        outcomes[arm]=outcome(measurements[arm],work,exception,arm=='native',spec.get('argument',False))
        if spec.get('parent_denied'):work.chmod(0o700)
        after[arm]=snapshot(work)
        preserved[arm]=all(after[arm].get(key)==value for key,value in before[arm].items() if key in ('payload','target'))
    exact=outcomes['python']==outcomes['native'] and after['python']==after['native'] and all(preserved.values())
    value=dict(spec=spec,measurements=measurements,outcomes=outcomes,before=before,after=after,source_preserved=preserved,exact=exact,
        stderr_bytes_equal=Path(measurements['python']['stderr']['path']).read_bytes()==Path(measurements['native']['stderr']['path']).read_bytes(),unqualified=spec.get('unqualified'))
    save(folder/'comparison.json',value);return value

def check(args):
    args.work_root.mkdir();results={}
    for name,spec in cases().items():results[name]=compare_case(args.work_root,name,spec,args.producer,args.binary)
    qualified={k:v for k,v in results.items() if not v['unqualified']};unqualified={k:v for k,v in results.items() if v['unqualified']}
    failed=[k for k,v in qualified.items() if not v['exact']]
    value=dict(producer=pin(args.producer),binary=pin(args.binary),comparator=pin(__file__),cases=results,
        summary=dict(qualified=len(qualified),successes=sum(v['outcomes']['python']['exit']==0 for v in qualified.values()),
            failed=failed,unqualified=len(unqualified),unqualified_mismatches=[k for k,v in unqualified.items() if not v['exact']],
            stderr_byte_differences=sum(not v['stderr_bytes_equal'] for v in qualified.values())),passed=not failed,
        limitations=['Unpredictable temporary basenames are normalized; bytes/mode/publication preservation remain exact.',
            'Native category envelope and Python CLI diagnostics are recorded separately, not called byte-identical.',
            'Shared JSON representation refusals remain unqualified, never matched failures.'])
    save(args.work_root/'comparison.json',value);print(json.dumps(value['summary']));return 0 if value['passed'] else 1

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='procedure',required=True)
    q=sub.add_parser('check')
    for key in ('producer','binary','work-root'):q.add_argument('--'+key,type=Path,required=True)
    q.set_defaults(function=check)
    q=sub.add_parser('oracle');q.add_argument('--producer',type=Path,required=True);q.add_argument('--exception',type=Path,required=True)
    q.add_argument('arguments',nargs=argparse.REMAINDER);q.set_defaults(function=oracle)
    args=p.parse_args()
    if args.procedure=='oracle' and args.arguments[:1]==['--']:args.arguments=args.arguments[1:]
    raise SystemExit(args.function(args))

if __name__=='__main__':main()
