#!/usr/bin/env python3
"""Compare the optional indexed-plugin helper with its unchanged Python source.

compare --producer SCRIPT --binary PRIVA --work-root NEW_DIRECTORY --output JSON
reference is the internal original-CLI launcher. bgzip, tabix and GNU sort must
already be on PATH. Fixtures use private directories and never installed data.
Saved Python arms may be reused only with identical source and fixture bodies;
their files and original output hashes are checked before use.
"""
from __future__ import annotations
import argparse
import ast
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import shlex
import stat
import struct
import subprocess
import sys
import time
sys.dont_write_bytecode = True


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def fixture_body(source, name):
    node = next(node for node in ast.parse(source).body
                if isinstance(node, ast.FunctionDef) and node.name == name)
    return ast.dump(node, include_attributes=False)


def saved_case(folder):
    """Follow recorded reuse links without relabeling the producing arm."""
    seen=set()
    while not (folder/'python.json').is_file():
        assert folder not in seen, 'Cyclic saved-oracle reference'
        seen.add(folder)
        reference=json.loads((folder/'python_reused.json').read_text())
        target=Path(reference['path'])
        assert hashlib.sha256(target.read_bytes()).hexdigest()==reference['sha256']
        folder=target.parent
    return folder


def reference(args):
    """Observe the original CLI's own catch boundary, without editing its code."""
    spec=importlib.util.spec_from_file_location('pinned_plugin',args.source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    exception={}
    def trace(frame,event,value):
        if frame.f_code is module.main.__code__ and event=='exception':
            kind,error,_=value
            exception.update(category=kind.__name__,message=str(error))
        return trace
    sys.settrace(trace)
    try:
        code=module.main(json.loads(args.arguments.read_text()))
    except BaseException as error:
        code=error.code if isinstance(error,SystemExit) else 1
        if not isinstance(error,SystemExit):print(str(error),file=sys.stderr)
    finally:
        sys.settrace(None)
        args.exception_output.write_text(json.dumps(exception))
    raise SystemExit(code)


def identity(path):
    value=path.stat()
    return (value.st_dev,value.st_ino,value.st_mode,value.st_size,value.st_mtime_ns,value.st_ctime_ns)


def files(root):
    result={}
    for path in sorted(root.rglob('*')):
        name=str(path.relative_to(root))
        info=path.lstat()
        if path.is_symlink():
            result[name]={'kind':'link','target':os.readlink(path),'mode':stat.S_IMODE(info.st_mode)}
        elif path.is_file():
            data=path.read_bytes()
            result[name]={'kind':'file','bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),
                          'mode':stat.S_IMODE(info.st_mode)}
        elif path.is_dir():result[name]={'kind':'directory','mode':stat.S_IMODE(info.st_mode)}
    return result


def setup(root,case,tools):
    root.mkdir()
    table=root/'input.tsv.bgz';cache=root/'cache'
    header=b'#chr\tpos\tend\tvalue\n'
    rows=b'1\t10\t12\tA\n1\t20\t22\tB\n2\t5\t8\tC\n'
    payload=case.get('payload',header+rows)
    if case.get('missing'):return table,cache
    if case.get('directory'):
        table.mkdir();return table,cache
    if case.get('plain'):
        table.write_bytes(payload)
    elif case.get('gzip'):
        table.write_bytes(gzip.compress(payload,mtime=0))
    else:
        with table.open('wb') as out:
            subprocess.run([tools['bgzip'],'-c'],input=payload,stdout=out,check=True)
    if case.get('empty'):table.write_bytes(b'')
    if case.get('corrupt_crc'):
        data=bytearray(table.read_bytes());data[-36]^=1;table.write_bytes(data)
    if case.get('index') in ('fresh','stale'):
        command=[tools['tabix'],'-s','1','-b','2','-e','3' if case.get('kind')=='mavedb' else '2']
        if case.get('zero_based'):command.append('-0')
        if 'meta' in case:command.extend(['-c',chr(case['meta']) if case['meta'] else ''])
        if 'skip' in case:command.extend(['-S',str(case['skip'])])
        subprocess.run(command+[str(table)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        if case['index']=='stale':os.utime(str(table)+'.tbi',ns=(1_000_000_000,1_000_000_000))
    if 'settings' in case:
        n,fmt,seq,begin,end,meta,skip,names=case['settings']
        Path(str(table)+'.tbi').write_bytes(gzip.compress(
            b'TBI\x01'+struct.pack('<8i',n,fmt,seq,begin,end,meta,skip,names),mtime=0))
    if 'broken_index' in case:Path(str(table)+'.tbi').write_bytes(case['broken_index'])
    if case.get('index_directory'):Path(str(table)+'.tbi').mkdir()
    if case.get('read_only'):table.chmod(0o444)
    os.utime(table,ns=(2_000_000_000,2_000_000_000))
    return table,cache


def fixture_cases():
    cases=[]
    for kind in ['mavedb','splicevault','primateai']:
        for index in ['missing','fresh','stale']:
            cases.append({'name':f'{kind}_{index}','kind':kind,'index':index})
    for name,extra in [
        ('missing',{'missing':True}),('directory',{'directory':True}),('empty',{'empty':True}),
        ('plain',{'plain':True}),('ordinary_gzip',{'gzip':True}),('crc',{'corrupt_crc':True}),
        ('read_only',{'read_only':True}),('index_directory',{'index_directory':True}),
        ('broken_index',{'broken_index':b'not gzip'}),('short_index',{'broken_index':gzip.compress(b'TBI\x01',mtime=0)}),
        ('bad_coordinates',{'payload':b'1\twrong\t12\n'}),
        ('negative',{'payload':b'1\t-1\t2\n'}),('zero',{'payload':b'1\t0\t1\n'}),
        ('no_records',{'payload':b'#header\n\n'}),
        ('unterminated',{'payload':b'1\t10\t12\tA'}),
        ('crlf',{'payload':b'#header\r\n1\t10\t12\r\n'}),
        ('cr_only',{'payload':b'#header\r1\t10\t12\r'}),
        ('utf8',{'payload':'#é\n1\t10\t12\té\n'.encode()}),
        ('invalid_utf8',{'payload':b'#\xff\n1\t10\t12\n'}),
        ('unicode_coordinate',{'payload':'1\t١٠\t١٢\n'.encode()}),
        ('underscore_coordinate',{'payload':b'1\t1_0\t1_2\n'}),
        ('huge_coordinate',{'payload':b'1\t'+b'9'*70+b'\t'+b'9'*70+b'\n'}),
        ('4096_headers',{'payload':b'#header\n'*4096+b'1\t10\t12\n'}),
        ('4095_headers',{'payload':b'#header\n'*4095+b'1\t10\t12\n'}),
        ('long_header',{'payload':b'#'+b'a'*(1024*1024)+b'\n1\t10\t12\n'}),
        ('exact_limit_header',{'payload':b'#'+b'a'*(1024*1024-2)+b'\n1\t10\t12\n'}),
    ]:cases.append(dict(name=name,**extra))
    for name,settings in [
        ('columns',(1,0,1,2,3,35,0,2)),('nongeneric',(1,2,1,2,2,35,0,2)),
        ('zero_based',(1,65536,1,2,2,35,0,2)),('bad_format',(1,-1,1,2,2,35,0,2)),
        ('zero_refs',(0,0,1,2,2,35,0,2)),('meta',(1,0,1,2,2,33,0,2)),
        ('skip',(1,0,1,2,2,35,1,2)),('bad_seq',(1,0,0,2,2,35,0,2)),
        ('bad_meta',(1,0,1,2,2,128,0,2)),('negative_skip',(1,0,1,2,2,35,-1,2)),
        ('no_names',(1,0,1,2,2,35,0,0)),
    ]:
        case={'name':'settings_'+name,'settings':settings}
        if name=='meta':case['payload']=b'!header\n1\t10\t12\n'
        if name=='skip':case['payload']=b'ignored noncomment\n1\t10\t12\n'
        cases.append(case)
    original=b''.join(f'#metadata {i}\n'.encode() for i in range(11))+b'chr\tpos\tref\talt\tscore\n'
    records=b'2\t2\tA\tG\t0.5\n1\t10\tG\tT\t0.7\n1\t2\tC\tA\t0.2\n1\t2\tA\tG\t0.1\n'
    for assembly in ['hg19','GRCh37','hg38','GRCh38']:
        for repeat in [False,True]:
            cases.append({'name':f'raw_{assembly}_{repeat}','action':'format-primate','assembly':assembly,
                          'payload':original+records,'repeat':repeat})
    for name,payload in [('raw_last_no_newline',original+records.rstrip(b'\n')),
                         ('raw_blank_rows',original+b'\n'+records+b'\n'),
                         ('raw_crlf',(original+records).replace(b'\n',b'\r\n')),
                         ('raw_bad_allele_text',original+records.replace(b'0.5',b'\xff')),
                         ('raw_only_header',original),('raw_empty',b''),
                         ('prepared',b'#chr\tpos\tref\n1\t10\tA\n')]:
        cases.append({'name':name,'action':'format-primate','assembly':'GRCh38','payload':payload})
    for threads in ['0','-1','2','1_0','١','invalid']:
        cases.append({'name':'threads_'+threads,'action':'format-primate','assembly':'hg38',
                      'payload':original+records,'threads':threads})
    for tool,mode in [('bgzip','exit'),('bgzip','missing'),('tabix','exit'),('tabix','missing'),
                      ('sort','exit'),('sort','missing'),('bgzip','compression_exit')]:
        case={'name':f'tool_{tool}_{mode}','tool':tool,'tool_mode':mode}
        if tool=='sort' or mode=='compression_exit':
            case.update(action='format-primate',assembly='hg38',payload=original+records)
        cases.append(case)
    for mutation in ['source','algorithm','output_hash','index_hash','invalid_json','list_receipt',
                     'missing_receipt','missing_output','index_directory']:
        cases.append({'name':'reuse_'+mutation,'action':'format-primate','assembly':'hg38',
                      'payload':original+records,'repeat':True,'mutate_receipt':mutation})
    for tokens in [['--threads','bad','--help'],['--help','--threads','bad'],
                   ['--threads','bad','--threads','1'],['--threads','1','--threads','bad'],
                   ['--threads=-0'],['--threads=-000'],['--threads=+١'],['--threads=\u001c1'],
                   ['--assembly','bad','--help'],['--help','--assembly','bad']]:
        cases.append({'name':'admission_'+str(len(cases)),'action':'format-primate','assembly':'hg38',
                      'payload':original+records,'extra_argv':tokens})
    return cases


def child_environment(case,root,tools):
    env=dict(os.environ)
    if not case.get('tool'):return env
    bindir=root/'tools';bindir.mkdir(exist_ok=True)
    # Only selected test tools are intercepted; these are local executable
    # fixtures, never substitutes on a scientific run's PATH.
    tool=case['tool'];mode=case['tool_mode']
    if mode=='missing':
        for name,source in tools.items():
            if name!=tool:(bindir/name).symlink_to(source)
        if tool!='sort':(bindir/'sort').symlink_to(shutil.which('sort'))
        env['PATH']=str(bindir)
        return env
    if mode=='compression_exit':
        body='if [ "$1" = "-@" ]; then /bin/cat >/dev/null; exit 23; fi\nexec '+shlex.quote(tools[tool])+' "$@"\n'
    elif tool=='sort':
        body='/bin/cat >/dev/null\nexit 17\n'
    else:body='exit 17\n'
    executable=bindir/tool;executable.write_text('#!/bin/sh\n'+body);executable.chmod(0o755)
    env['PATH']=str(bindir)+os.pathsep+env['PATH'];return env


def mutate_receipt(cache,mutation):
    ready=next(cache.glob('*/*/source.ready.json'))
    data=json.loads(ready.read_text())
    if mutation=='missing_receipt':ready.unlink();return
    if mutation=='invalid_json':ready.write_text('{invalid');return
    if mutation=='list_receipt':ready.write_text('[]');return
    output=next(ready.parent.glob('*.bgz'))
    if mutation=='missing_output':output.unlink();return
    if mutation=='index_directory':
        index=Path(str(output)+'.tbi');index.unlink();index.mkdir();return
    if mutation=='source':data['source']['assembly']='different'
    if mutation=='algorithm':data['algorithm']='different'
    if mutation=='output_hash':data['output_sha256']='0'*64
    if mutation=='index_hash':data['index_sha256']='0'*64
    ready.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')


def normalized(text,root):
    text=text.replace(str(root),'<CASE>')
    # Temporary names are generated independently by both standard libraries.
    text=re.sub(r'(\.index-|\.build-)[A-Za-z0-9_]{6,8}',r'\1<TMP>',text)
    return text


def compare(args):
    binary=args.binary.resolve();source=args.producer.resolve()
    root=args.work_root.absolute();root.mkdir()
    shutil.copy2(__file__,root/Path(__file__).name)
    tools={name:shutil.which(name) for name in ['bgzip','tabix']}
    assert all(tools.values())
    cases=fixture_cases();rows=[];started=time.monotonic()
    reused={}
    if args.reuse:
        prior=json.loads(args.reuse.read_text());prior_root=Path(prior['root'])
        assert prior['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
        previous=prior_root/prior.get('controller_name','Indexed_Plugin_Validation.py')
        old_text=previous.read_text();new_text=Path(__file__).read_text()
        # Compare definitions without importing an old host-specific controller.
        for name in ('fixture_cases','setup','child_environment','mutate_receipt','normalized'):
            assert fixture_body(old_text,name)==fixture_body(new_text,name), name
        for case in cases:reused[case['name']]=saved_case(prior_root/case['name'])
    for case in cases:
        case_root=root/case['name'];case_root.mkdir()
        active=reused.get(case['name'],case_root)/'work'
        assert not active.exists()
        outcomes={}
        for arm in ['python','native']:
            if arm=='python' and case['name'] in reused:
                prior_case=reused[case['name']]
                outcomes[arm]=json.loads((prior_case/'python.json').read_text())
                assert files(prior_case/'python')==outcomes[arm]['files']
                for attempt,observation in enumerate(outcomes[arm]['runs']):
                    assert normalized((prior_case/f'python.{attempt}.stdout').read_text(),active)==observation['stdout']
                write_json(case_root/'python_reused.json',{
                    'path':str(prior_case/'python.json'),
                    'sha256':hashlib.sha256((prior_case/'python.json').read_bytes()).hexdigest()})
                continue
            table,cache=setup(active,case,tools)
            env=child_environment(case,active,tools)
            argv=[case.get('action','ensure-index'),'--input',str(table)]
            if argv[0]=='ensure-index':argv+=['--kind',case.get('kind','primateai')]
            else:
                argv+=['--cache-dir',str(cache),'--assembly',case.get('assembly','hg38'),
                       '--threads',case.get('threads','1')]
            argv+=case.get('extra_argv',[])
            arguments=case_root/(arm+'.arguments.json');arguments.write_text(json.dumps(argv))
            error_output=case_root/(arm+'.exception.json')
            source_before=identity(table) if table.is_file() else None
            before=files(active)
            observations=[]
            for attempt in range(2 if case.get('repeat') else 1):
                if attempt and case.get('mutate_receipt'):mutate_receipt(cache,case['mutate_receipt'])
                command=([sys.executable,__file__,'reference','--source',str(source),
                          '--arguments',str(arguments),'--exception-output',str(error_output)]
                         if arm=='python' else [str(binary),'indexed-plugin-cache',*argv])
                begin=time.monotonic();run=subprocess.run(command,capture_output=True,env=env)
                stdout=run.stdout.decode('utf-8',errors='backslashreplace')
                stderr=run.stderr.decode('utf-8',errors='backslashreplace')
                category=None
                if run.returncode:
                    if arm=='python':
                        detail=json.loads(error_output.read_text())
                        category=detail.get('category')
                        message=detail.get('message','')
                    else:
                        final=stderr.splitlines()[-1] if stderr.splitlines() else ''
                        match=re.match(r'priva: ([A-Za-z_]+): (.*)',final)
                        category,message=match.groups() if match else ('ArgumentError',final.removeprefix('priva: '))
                    if run.returncode==2:
                        category='ArgumentError'
                        if arm=='python':
                            message=stderr.splitlines()[-1].split(': error: ',1)[-1]
                else:message=''
                observations.append({'exit':run.returncode,'stdout':normalized(stdout,active),
                    'category':category,'message':normalized(message,active),
                    'seconds':time.monotonic()-begin})
                (case_root/f'{arm}.{attempt}.stdout').write_bytes(run.stdout)
                (case_root/f'{arm}.{attempt}.stderr').write_bytes(run.stderr)
            unchanged=source_before is None or identity(table)==source_before
            outcome={'runs':observations,'files':files(active),'input_identity_unchanged':unchanged}
            outcomes[arm]=outcome
            shutil.copytree(active,case_root/arm,symlinks=True)
            shutil.rmtree(active)
            write_json(case_root/(arm+'.json'),outcome)
        a,b=outcomes['python'],outcomes['native']
        mismatch=[]
        for i,(left,right) in enumerate(zip(a['runs'],b['runs'])):
            for field in ['exit','stdout','category','message']:
                if left[field]!=right[field]:mismatch.append({'attempt':i,'field':field,'python':left[field],'native':right[field]})
        if a['files']!=b['files']:
            mismatch.append({'field':'artifact_inventory','different_paths':[
                name for name in sorted(set(a['files'])|set(b['files'])) if a['files'].get(name)!=b['files'].get(name)]})
        if not a['input_identity_unchanged'] or not b['input_identity_unchanged']:
            mismatch.append({'field':'selected_input_identity'})
        rows.append({'name':case['name'],'passed':not mismatch,'mismatches':mismatch})
    result={'root':str(root),'controller_name':Path(__file__).name,'producer':str(source),'binary':str(binary),
        'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'binary_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),
        'cases':len(rows),'passed':all(row['passed'] for row in rows),'results':rows,
        'reused_python_cases':len(reused),'new_python_cases':len(rows)-len(reused),
        'seconds':time.monotonic()-started,'scope':'synthetic actual-tool fixtures; no installed-cache or full-resource qualification'}
    write_json(root/'comparison.json',result);write_json(args.output,result)
    print(json.dumps({'passed':result['passed'],'cases':len(rows),
                      'failed':[r['name'] for r in rows if not r['passed']]}))
    if not result['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    subs=parser.add_subparsers(dest='procedure',required=True)
    cmd=subs.add_parser('reference');cmd.add_argument('--source',type=Path,required=True)
    cmd.add_argument('--arguments',type=Path,required=True);cmd.add_argument('--exception-output',type=Path,required=True)
    cmd.set_defaults(function=reference)
    cmd=subs.add_parser('compare');cmd.add_argument('--producer',type=Path,required=True)
    cmd.add_argument('--binary',type=Path,required=True);cmd.add_argument('--work-root',type=Path,required=True)
    cmd.add_argument('--output',type=Path,required=True)
    cmd.add_argument('--reuse',type=Path,help='Reuse immutable same-source/same-fixture Python arms')
    cmd.set_defaults(function=compare)
    args=parser.parse_args();args.function(args)


if __name__=='__main__':main()
