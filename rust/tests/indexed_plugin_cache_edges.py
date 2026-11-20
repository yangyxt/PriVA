#!/usr/bin/env python3
"""Bounded synthetic failure-order review, without installed resource writes.

check --producer SCRIPT --binary PRIVA --work-root NEW_DIRECTORY --output JSON
The actual tool commands are retained; explicit local wrappers inject one
permission or signal fault and record subprocess order for each arm.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent

def pin(path):
    path=Path(path)
    with path.open('rb') as stream:sha=hashlib.file_digest(stream,'sha256').hexdigest()
    return dict(path=str(path),bytes=path.stat().st_size,sha256=sha)

def save(path,value):
    with Path(path).open('x') as stream:json.dump(value,stream,indent=2,sort_keys=True);stream.write('\n')

def snapshot(root):
    result={}
    for folder,dirs,names in os.walk(root,followlinks=False):
        for name in dirs+names:
            path=Path(folder)/name;value=path.lstat();key=str(path.relative_to(root));mode=stat.S_IMODE(value.st_mode)
            item=dict(mode=mode)
            if stat.S_ISLNK(value.st_mode):item.update(kind='link',target=os.readlink(path))
            elif stat.S_ISREG(value.st_mode):item.update(kind='file',bytes=value.st_size,sha256=pin(path)['sha256'])
            else:
                item.update(kind='directory')
                if not mode&0o500:item['unreadable_directory']=True
            result[key]=item
    return result

def normalize(value,root):
    return re.sub(r'(\.index-|\.build-)[A-Za-z0-9_]{6,8}',r'\1<TMP>',value.replace(str(root),'<CASE>'))

def check(args):
    binary=args.binary.resolve();source=args.producer.resolve()
    root=args.work_root.absolute();root.mkdir(parents=True)
    frozen=root/'indexed_plugin_cache_compare.py';shutil.copyfile(P/'indexed_plugin_cache_compare.py',frozen)
    shutil.copyfile(__file__,root/Path(__file__).name)
    tools={key:shutil.which(key) for key in ('bgzip','tabix','sort')};assert all(tools.values())
    previous=json.loads(args.reuse_report.read_text()) if args.reuse_report else None
    if previous:assert previous['producer']['sha256']==pin(source)['sha256'] and previous.get('oracle_executed')
    reuse={row['case']:row for row in previous['cases']} if previous else {}
    results=[]
    for label in ('destination_permission','cleanup_permission','signal_usr1','cr_utf8_boundary','cleanup_permission_failure','cleanup_parent_permission','cr_eof','lf_bad_utf8','cr_lf_bad_utf8','cr_before_boundary','signal_term','signal_rt35','sort_pipe_crc','bgzip_pipe','sort_pipe_crc_header8192','sort_pipe_crc_total8192','sort_pipe_crc_total8193','sort_pipe_crc_header4096','sort_pipe_crc_total4096','sort_pipe_crc_total4097'):
        if label in reuse and previous['binary']['sha256']==pin(binary)['sha256']:
            results.append(reuse[label]);continue
        case=root/label;case.mkdir();work=Path(reuse[label]['arms']['python']['stdout_identity']['path']).parent/'work' if label in reuse else case/'work';outcomes={}
        for arm in ('python','native'):
            if arm=='python' and label in reuse:
                old=reuse[label]['arms']['python']
                for field in ('stdout_identity','stderr_identity'):assert pin(old[field]['path'])==old[field]
                old_root=Path(old['stdout_identity']['path']).parent/'python'
                for name,entry in old['files'].items():
                    path=old_root/name
                    if entry['kind']=='file':assert pin(path)['sha256']==entry['sha256'] and path.stat().st_size==entry['bytes']
                    elif entry['kind']=='link':assert path.is_symlink() and os.readlink(path)==entry['target']
                    else:assert path.is_dir()
                outcomes[arm]=old;continue
            work.mkdir();table=work/'input.tsv.bgz';cache=work/'cache';env=dict(os.environ)
            payload=b'1\t10\t10\tA\n1\t20\t20\tB\n'
            if label.startswith('sort_pipe_crc') or label=='bgzip_pipe':
                import random
                data=random.Random(124).randbytes(262140).hex().encode()
                payload=(b'#'+data+b'\n')*11+b'chr\tpos\tvalue\n1\t10\tA\n'
                if label in ('sort_pipe_crc_header8192','sort_pipe_crc_header4096'):
                    size=int(label.rsplit('header',1)[1]);payload=(b'#'+data+b'\n')*11+b'chr\tpos\t'+b'A'*(size-10)+b'\n'
                if label.startswith('sort_pipe_crc_total'):
                    size=int(label.rsplit('total',1)[1])
                    payload=(b'#'+data+b'\n')*11+b'chr\tpos\tvalue\n'+b'1\t10\t'+b'A'*(size-15-6)+b'\n'
            if label=='destination_permission':payload=b'# header\n'*11+b'chr\tpos\tvalue\n1\t10\tA\n'
            if label.startswith(('cr_','lf_')):
                endings={'cr_utf8_boundary':b'\r\xff\n','cr_eof':b'\r','lf_bad_utf8':b'\n\xff\n','cr_lf_bad_utf8':b'\r\n\xff\n','cr_before_boundary':b'\r\xff\n'}
                payload=b'1\t10\t10\t'+b'A'*((8190 if label=='cr_before_boundary' else 8191)-8)+endings[label]
            with table.open('wb') as output:subprocess.run([tools['bgzip'],'-c'],input=payload,stdout=output,check=True)
            if label.startswith('sort_pipe_crc'):
                data=bytearray(table.read_bytes());data[-36]^=1;table.write_bytes(data)
            if label=='destination_permission':
                identity=dict(path=str(table.absolute()),sha256=pin(table)['sha256'],size=table.stat().st_size,assembly='hg38')
                key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:20]
                parent=cache/'hg38';parent.mkdir(parents=True);blocked=work/'blocked';blocked.mkdir();(blocked/'target').mkdir();blocked.chmod(0)
                (parent/key).symlink_to(blocked/'target',target_is_directory=True)
                argv=['format-primate','--input',str(table),'--cache-dir',str(cache),'--assembly','hg38']
            elif label.startswith('sort_pipe_crc') or label=='bgzip_pipe':argv=['format-primate','--input',str(table),'--cache-dir',str(cache),'--assembly','hg38']
            else:argv=['ensure-index','--input',str(table),'--kind','primateai']
            if label.startswith(('cleanup_','signal_')) or label.startswith('sort_pipe_crc') or label in ('destination_permission','bgzip_pipe'):
                bindir=work/'tools';bindir.mkdir();env['PATH']=str(bindir)+os.pathsep+env['PATH']
                for tool in ('bgzip','tabix','sort'):
                    script='#!'+sys.executable+'\nimport json,os,pathlib,signal,subprocess,sys\n'
                    script+='with open('+repr(str(work/'toolcalls.jsonl'))+',"a") as stream:stream.write(json.dumps([sys.argv[0],*sys.argv[1:]])+"\\n")\n'
                    if label.startswith('signal_') and tool=='bgzip':script+='os.kill(os.getpid(),'+str({'signal_usr1':10,'signal_term':15,'signal_rt35':35}[label])+')\n'
                    else:
                        script+='code=subprocess.call(['+repr(tools[tool])+',*sys.argv[1:]])\n'
                        if label in ('cleanup_permission','cleanup_permission_failure') and tool=='tabix':
                            script+='if code==0 and "-s" in sys.argv:\n parent=pathlib.Path(sys.argv[-1]).parent/"blocked"\n parent.mkdir()\n (parent/"retained").write_text("synthetic fixture\\n")\n parent.chmod(0)\n'
                        if label=='cleanup_parent_permission' and tool=='tabix':script+='if code==0 and \"-s\" in sys.argv:pathlib.Path(sys.argv[-1]).parent.parent.chmod(0o500)\n'
                        if label=='cleanup_permission_failure' and tool=='tabix':script+='if \"-s\" in sys.argv:code=17\n'
                        script+='sys.exit(code)\n'
                    if label.startswith('sort_pipe_crc') and tool=='sort':script='#!/bin/sh\nexec 0<&-\nexit 17\n'
                    if label=='bgzip_pipe' and tool=='bgzip':script='#!/bin/sh\nif [ "$1" = "-@" ]; then exec 0<&-; exit 17; fi\nexec '+__import__('shlex').quote(tools['bgzip'])+' "$@"\n'
                    path=bindir/tool;path.write_text(script);path.chmod(0o755)
            argument_file=case/(arm+'.argv.json');save(argument_file,argv);exception=case/(arm+'.exception.json')
            command=[sys.executable,str(frozen),'reference','--source',str(source),'--arguments',str(argument_file),'--exception-output',str(exception)] if arm=='python' else [str(binary),'indexed-plugin-cache',*argv]
            before=pin(table)
            if label in reuse:
                original=reuse[label]['arms']['python']
                assert before['sha256']==original['files']['input.tsv.bgz']['sha256']
                saved_argv=json.loads((Path(original['stdout_identity']['path']).parent/'python.argv.json').read_text())
                assert argv==saved_argv
                for name,item in snapshot(work).items():
                    if name.startswith('tools/') and item['kind']=='file':
                        # Reuse the original fault-injection executable itself,
                        # including its absolute alias path and original mode.
                        # A later fixture spelling (SIGUSR1 versus integer 10)
                        # must not be relabeled as the old producing bytes.
                        prior=Path(original['stdout_identity']['path']).parent/'python'/name
                        assert pin(prior)['sha256']==original['files'][name]['sha256']
                        shutil.copyfile(prior,work/name)
                        (work/name).chmod(original['files'][name]['mode'])
                        assert snapshot(work)[name]==original['files'][name]
            run=subprocess.run(command,capture_output=True,env=env,timeout=60)
            stdout=case/(arm+'.stdout');stderr=case/(arm+'.stderr');stdout.write_bytes(run.stdout);stderr.write_bytes(run.stderr)
            observed=snapshot(work);assert pin(table)==before
            if arm=='python':assert exception.exists(), 'Python oracle did not execute the original helper'
            calls=[json.loads(line) for line in (work/'toolcalls.jsonl').read_text().splitlines()] if (work/'toolcalls.jsonl').exists() else []
            outcomes[arm]=dict(exit=run.returncode,stdout=normalize(run.stdout.decode(errors='backslashreplace'),work),stderr=normalize(run.stderr.decode(errors='backslashreplace'),work),
                exception=json.loads(exception.read_text()) if exception.exists() else None,files=observed,tool_calls=calls,
                selected_input_preserved=True,command=command,stdout_identity=pin(stdout),stderr_identity=pin(stderr))
            adjustments=[]
            work_mode=stat.S_IMODE(work.stat().st_mode)
            if work_mode!=0o755:
                work.chmod(0o755);adjustments.append(dict(path=str(work),producing_mode=work_mode,storage_mode=0o755))
            for folder,dirs,names in os.walk(work,followlinks=False):
                for name in dirs:
                    path=Path(folder)/name
                    if path.is_symlink():continue
                    mode=stat.S_IMODE(path.stat().st_mode)
                    if not mode&0o500:path.chmod(mode|0o700);adjustments.append(dict(path=str(path),producing_mode=mode,storage_mode=mode|0o700))
            outcomes[arm]['storage_adjustments']=adjustments
            work.rename(case/arm)
        row=dict(case=label,arms=outcomes);save(case/'comparison.json',row);results.append(row)
    def comparable(arm):
        category=None;message=''
        if arm['exit']:
            if arm['exception'] is not None:category=arm['exception'].get('category');message=arm['exception'].get('message','')
            else:
                match=re.search(r'(?:^|\n)priva: ([A-Za-z_]+): (.*)',arm['stderr'],re.S)
                category,message=match.groups() if match else ('unknown',arm['stderr'])
            active=Path(arm['stdout_identity']['path']).parent/'work'
            # A reused arm keeps its original absolute active fixture path.
            message=normalize(message,active).rstrip('\n')
        artifacts={normalize(k,Path('/unmatchable')):v for k,v in arm['files'].items() if k!='toolcalls.jsonl' and k!='tools' and not k.startswith('tools/')}
        if '--input' in arm['command']:
            arguments=arm['command']
        else:
            arguments=json.loads((Path(arm['stdout_identity']['path']).parent/'python.argv.json').read_text())
        active=Path(arguments[arguments.index('--input')+1]).parent
        calls=[[normalize(item,active) for item in call] for call in arm['tool_calls']]
        return dict(exit=arm['exit'],stdout=arm['stdout'],category=category,message=message,files=artifacts,tool_calls=calls)
    for row in results:
        left,right=map(comparable,[row['arms']['python'],row['arms']['native']]);row['differences']=[k for k in left if left[k]!=right[k]];row['passed']=not row['differences']
    value=dict(passed=all(row['passed'] for row in results),reference_report=pin(args.reuse_report) if args.reuse_report else None,binary=pin(binary),producer=pin(source),oracle=pin(frozen),controller=pin(root/Path(__file__).name),oracle_executed=True,tools={k:pin(v) for k,v in tools.items()},cases=results)
    save(root/'review.json',value);save(args.output,value)
    print(json.dumps(dict(passed=value['passed'],cases=len(results),failures=[dict(case=row['case'],differences=row['differences']) for row in results if not row['passed']])))
    if not value['passed'] and not args.expect_mismatches:raise RuntimeError('Edge comparison failed; all measurements retained')

def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='procedure',required=True)
    p=sub.add_parser('check');p.add_argument('--producer',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--reuse-report',type=Path);p.add_argument('--expect-mismatches',action='store_true');p.set_defaults(run=check)
    args=parser.parse_args();args.run(args)

if __name__=='__main__':main()
