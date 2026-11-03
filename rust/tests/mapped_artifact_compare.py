#!/usr/bin/env python3
"""Compare mapped publication procedures and filesystem effects with their source.

python mapped_artifact_compare.py check --repo CHECKOUT --binary PRIVA --work-root EMPTY

Scientific producer scripts here are labelled tiny byte fixtures. No scientific
calculation, public download or installed resource write occurs. Fault injection
interposes libc rename only in these child processes, never in production code.
"""
from __future__ import annotations
import argparse,contextlib,fcntl,hashlib,importlib.util,io,json,os,shutil,stat,subprocess,sys
from pathlib import Path
from types import SimpleNamespace


def write(path,value):
    Path(path).write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')

def load(path):
    spec=importlib.util.spec_from_file_location('mapped_fixture_'+str(abs(hash(str(path)))),path)
    value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

PRODUCER=r'''# Labelled byte producer fixture, never scientific calculation.
import json,os,signal,sys
from pathlib import Path
args=sys.argv[1:]
Path(os.environ['FIXTURE_ARGV']).write_text(json.dumps(args))
output=Path(args[args.index('--output')+1]);mode=os.environ.get('FIXTURE_MODE','success')
if mode=='missing':sys.exit(0)
if mode=='symlink':output.symlink_to(os.environ['FIXTURE_OUTSIDE']);sys.exit(0)
if mode=='hardlink':os.link(os.environ['FIXTURE_OUTSIDE'],output);sys.exit(0)
output.write_bytes(b'' if mode=='empty' else b'labelled fixture payload\x00with duplicates B A B\n')
os.utime(output,ns=(1700000000000000000,1700000000000000000))
if mode=='log':(output.parent/'producer.log').write_text('retained fixture log\n')
if mode in ('config','mapping','producer','prior'):
    p=Path(os.environ['FIXTURE_MUTATE']);p.write_bytes(p.read_bytes()+b'changed\n')
if mode=='exit':sys.exit(9)
if mode=='signal':os.kill(os.getpid(),signal.SIGTERM)
'''
SHIM=r'''#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int ends(const char *s,const char *tail){size_t a=strlen(s),b=strlen(tail);return a>=b&&!strcmp(s+a-b,tail);}
int rename(const char *a,const char *b){
 static int (*next)(const char*,const char*)=0; if(!next)next=dlsym(RTLD_NEXT,"rename");
 int rc=next(a,b);int saved=errno;const char *log=getenv("FIXTURE_RENAMES");
 if(log){FILE *f=fopen(log,"a");if(f){fprintf(f,"%d\t%s\t%s\n",rc,a,b);fclose(f);}}
 const char *hit=getenv("FIXTURE_FAIL_RENAME"),*kind=getenv("FIXTURE_FAIL_KIND");
 if(!rc&&((hit&&!strcmp(hit,b))||(kind&&!strcmp(kind,"staged_receipt")&&strstr(b,"/.mapped-artifact-")&&ends(b,".nm_to_enst.json")))){errno=EIO;return -1;}
 errno=saved;return rc;
}
'''

class Fixture:
    def __init__(self,root,repo,shim,kind='motifs',bound=False):
        self.root=root
        alpha=load(repo/'rust/tests/alpha_annotation_base_compare.py') if bound else None
        base=alpha.Fixture(root,repo) if bound else None
        if not bound:root.mkdir(parents=True)
        self.cache=base.cache if bound else root/'managed cache';self.output=self.cache/'plugins/AlphaMissense';self.output.mkdir(parents=True,exist_ok=True)
        self.source=base.source if bound else root/'source';self.scripts=self.source/'scripts';self.scripts.mkdir(parents=True,exist_ok=True)
        for name in ['mapped_artifact_publish.py','biomart_nm_to_enst.py','biomart_ftp_nm_to_enst.py']:
            shutil.copy2(repo/'scripts'/name,self.scripts/name)
        self.kind=kind;self.producer=self.scripts/('am_pick_missense_constrained_segments.py' if kind=='motifs' else 'prepare_pm1_regions.py');self.producer.write_text(PRODUCER)
        self.config=base.config if bound else root/'runtime config.json'
        if not bound:self.config.write_text('{"fixture":"actual runtime configuration"}\n')
        self.mapping=root/'selected mapping.tsv';self.mapping.write_text('deliberately unpinned fixture mapping\n')
        self.inputs=[]
        for name in ['gnomad.tsv','spots.tsv']:
            p=root/name;p.write_text('fixture anchor\n');self.inputs.append(p)
        self.target=self.output/'selected.pkl';self.outside=root/'outside.pkl';self.outside.write_bytes(b'unchanged outside fixture')
        self.env=dict(base.env if bound else os.environ,PYTHONDONTWRITEBYTECODE='1',LD_PRELOAD=str(shim),FIXTURE_RENAMES=str(root/'renames.tsv'),FIXTURE_ARGV=str(root/'argv.json'),FIXTURE_OUTSIDE=str(self.outside))
        self.module=load(self.scripts/'mapped_artifact_publish.py')
        self.kw=dict(artifact=str(self.target),mapping=str(self.mapping),input=list(map(str,self.inputs)),cache_root=str(self.cache),output_dir=str(self.output),assembly='hg38',config=str(self.config),producer=str(self.producer),kind=kind,annotation_base=None,native_lineage=None)
        self.arguments=['--literal',"spaces ' quotes ; $(never_exec) `never_exec`",'--output',str(self.target),'--skip_das']
        self.expected=json.dumps(self.module.BIOMART.mapping_identity(self.mapping,self.inputs))
        if bound:
            base.published();self.kw['annotation_base']=base.cfg['alphamissense_vep_base_vcf']
    def complete(self,target=None):
        target=target or self.target;target.write_bytes(b'previous complete fixture')
        self.module.BIOMART.atomic_json(Path(str(target)+'.nm_to_enst.json'),dict(mapping_identity=self.module.BIOMART.mapping_identity(self.mapping,self.inputs),artifact_bytes=target.stat().st_size,artifact_sha256=digest(target)))
        # Pin only freshly generated fixture metadata before observing it.
        # A filesystem may round a just-created mtime by one nanosecond later;
        # production publication retains its exact before/after identity checks.
        for p in [target,Path(str(target)+'.nm_to_enst.json')]:os.utime(p,ns=(1700000000000000000,1700000000000000000))
    def pair(self):
        values=[]
        for p in [self.target,Path(str(self.target)+'.nm_to_enst.json')]:
            if not os.path.lexists(p):values.append(None);continue
            m=p.lstat();v=dict(inode=m.st_ino,size=m.st_size,mtime_ns=m.st_mtime_ns,ctime_ns=m.st_ctime_ns,mode=m.st_mode,nlink=m.st_nlink,uid=m.st_uid)
            if stat.S_ISREG(m.st_mode):v['sha256']=digest(p)
            if p.is_symlink():v['link']=os.readlink(p)
            values.append(v)
        return values
    def lineage(self):
        source=self.root/'historical.pkl';self.complete(source)
        value=json.loads(Path(str(source)+'.nm_to_enst.json').read_text());value['receipt_sha256']=digest(Path(str(source)+'.nm_to_enst.json'))
        path=self.root/'lineage.json';write(path,value);return source,path
    def marker(self):return Path(str(self.target)+'.publish.pending.json')


def reference(args):
    request=json.loads(args.request.read_text());m=load(Path(request['helper'])/'scripts/mapped_artifact_publish.py')
    values=request['options']
    for key in ['artifact','mapping','cache_root','output_dir','config','producer','annotation_base','native_lineage']:
        if values.get(key) is not None:values[key]=Path(values[key])
    values['input']=list(map(Path,values['input']));values['arguments']=(['--'] if request['command']=='run' else [])+request['arguments'];values['expected_identity']=request['expected']
    selected=SimpleNamespace(**values)
    try:
        if request['command']=='adopt-native':result=m.adopt_native(selected,Path(request['source']),json.loads(Path(request['lineage']).read_text()))
        elif request['command']=='context':result=m.context(selected)
        else:result=getattr(m,request['command'])(selected)
        value=dict(success=True,value=result)
    except Exception as error:value=dict(success=False,category=type(error).__name__,message=str(error))
    write(args.output,value)
    return 0


class Harness:
    def __init__(self,args):
        self.args=args;self.root=args.work_root;self.root.mkdir();self.rows=[];self.counter=0
        c=self.root/'rename_fixture.c';c.write_text(SHIM);self.shim=self.root/'rename_fixture.so'
        subprocess.run(['cc','-shared','-fPIC','-O0',str(c),'-o',str(self.shim),'-ldl'],check=True)
    def fixture(self,label,arm,kind='motifs'):
        return Fixture(self.root/(f'{self.counter:03}_{label}_{arm}'),self.args.repo,self.shim,kind)
    def normalize(self,value,f):
        if not isinstance(value,str):return value
        value=value.replace(str(f.root),'<FIXTURE>')
        import re
        return re.sub(r'\.mapped-artifact-(motifs|pm1)-[^/\s\'"\]]+',r'.mapped-artifact-\1-<STAGE>',value)
    def invoke(self,f,command,arm,source=None,lineage=None):
        index=len(list(f.root.glob('request_*.json')));tag=str(index)
        request=dict(helper=str(f.source),options=f.kw,arguments=f.arguments,expected=f.expected,command=command,source=str(source) if source else None,lineage=str(lineage) if lineage else None)
        request_path=f.root/('request_'+tag+'.json');write(request_path,request)
        if arm=='reference':
            output=f.root/('outcome_'+tag+'.json')
            argv=[sys.executable,__file__,'reference','--request',str(request_path),'--output',str(output)]
        else:
            argv=[str(self.args.binary),'mapped-artifact',command,'--helper-root',str(f.source),'--python',sys.executable]
            for key,value in f.kw.items():
                if key=='input':
                    for path in value:argv.extend(['--input',path])
                elif value is not None:argv.extend(['--'+key.replace('_','-'),value])
            if command=='run':argv.extend(['--expected-identity',f.expected,'--',*f.arguments])
            if command=='adopt-native':argv.extend(['--source',str(source),'--lineage',str(lineage)])
        result=subprocess.run(argv,env=f.env,capture_output=True,text=True)
        (f.root/('stdout_'+tag+'.log')).write_text(result.stdout);(f.root/('stderr_'+tag+'.log')).write_text(result.stderr)
        if arm=='reference':
            if result.returncode:raise RuntimeError('reference observer failed: '+result.stderr)
            value=json.loads(output.read_text());value['returncode']=0 if value['success'] else 1
        else:
            value=dict(success=result.returncode==0,returncode=result.returncode)
            if result.returncode:
                category,_,message=result.stderr.strip().rsplit('priva: ',1)[-1].partition(': ');value.update(category=category,message=message)
            elif command=='context':value['value']=json.loads(result.stdout)
        value.update(argv=argv,stdout=result.stdout)
        return value
    def record(self,label,arms,*,success=None,message=True,extra=True):
        ref,native=arms
        passed=ref['success']==native['success'] and native['returncode']==ref['returncode'] and extra
        if success is not None:passed &= ref['success']==success and native['success']==success
        if not ref['success']:
            passed &= ref['category']==native.get('category')
            if message:passed &= ref['normalized_message']==native.get('normalized_message')
        if ref['success'] and 'value' in native:passed &= ref.get('value')==native['value']
        self.rows.append(dict(name=label,passed=bool(passed),reference=ref,native=native,expected_success=success,exact_message_checked=not ref['success'] and message))
        write(self.root/'progress.json',dict(cases=len(self.rows),passed=sum(r['passed'] for r in self.rows),rows=self.rows))
    def simple(self,label,setup=lambda f:None,*,command='run',kind='motifs',success=None,preserve=False,message=True):
        self.counter+=1;arms=[];contents=[]
        for arm in ['reference','native']:
            f=self.fixture(label,arm,kind);setup(f);before=f.pair();outside=dict(bytes=f.outside.read_bytes(),stat=f.outside.stat())
            value=self.invoke(f,command,arm)
            value['normalized_message']=self.normalize(value.get('message',''),f)
            value['prior_pair']=before;value['after_pair']=f.pair()
            value['preserved']=not preserve or before==value['after_pair']
            value['outside_unchanged']=outside['bytes']==f.outside.read_bytes() and (f.env.get('FIXTURE_MODE')=='hardlink' or outside['stat'].st_ctime_ns==f.outside.stat().st_ctime_ns)
            if value['success'] and command=='run':
                f.module.BIOMART.check_receipt(f.target,Path(f.kw['mapping']),list(map(Path,f.kw['input'])))
                contents.append((f.target.read_bytes(),Path(str(f.target)+'.nm_to_enst.json').read_bytes()))
                argv=json.loads((f.root/'argv.json').read_text());at=f.arguments.index('--output')+1;argv[at]=f.arguments[at]
                value['argv_unchanged_except_output']=argv==f.arguments
            arms.append(value)
        self.record(label,arms,success=success,message=message,extra=all(x['preserved'] and x['outside_unchanged'] and x.get('argv_unchanged_except_output',True) for x in arms) and (not contents or contents[0]==contents[1]))
    def basic(self):
        for kind in ['motifs','pm1']:
            self.simple('cold_'+kind,kind=kind,success=True)
            self.simple('replace_'+kind,lambda f:f.complete(),kind=kind,success=True)
            self.simple('recover_absent_'+kind,command='recover',kind=kind,success=True,preserve=True)
        for mode in ['exit','signal','empty','missing','symlink','hardlink','config','mapping','producer','prior','log']:
            def setup(f,mode=mode):
                f.complete();f.env['FIXTURE_MODE']=mode
                if mode in ['config','mapping','producer','prior']:f.env['FIXTURE_MUTATE']=str({'config':f.config,'mapping':f.mapping,'producer':f.producer,'prior':f.target}[mode])
            self.simple('producer_'+mode,setup,success=mode=='log',preserve=mode not in ['prior','log'],message=True)
        for label,arguments in [('missing_output',['literal']),('duplicate_output',['--output','x','--output','y']),('output_no_value',['--output']),('different_output',['--output','other'])]:
            self.simple(label,lambda f,a=arguments:setattr(f,'arguments',a),success=False,preserve=True)
        self.simple('wrong_identity',lambda f:setattr(f,'expected','{}'),success=False,preserve=True)
        self.simple('wide_ignored_identity',lambda f:setattr(f,'expected',f.expected[:-1]+',"extra":'+str(2**150)+'}'),success=False,preserve=True)
        self.simple('surrogate_ignored_identity',lambda f:setattr(f,'expected',f.expected[:-1]+',"extra":"\\ud800"}'),success=False,preserve=True)
        self.simple('unknown_existing',lambda f:f.target.write_bytes(b'unknown supplied bytes'),success=False,preserve=True)
        self.simple('orphan_receipt',lambda f:Path(str(f.target)+'.nm_to_enst.json').write_text('{}'),success=False,preserve=True)
        self.simple('recover_unknown_existing',lambda f:f.target.write_bytes(b'unknown supplied bytes'),command='recover',success=True,preserve=True)
        for which in ['producer','config']:
            self.simple('relative_'+which,lambda f,k=which:f.kw.update({k:'relative'}),success=False,preserve=True)
        self.simple('producer_wrong_name',lambda f:f.kw.update(producer=str(f.config)),success=False,preserve=True)
        self.simple('relative_target',lambda f:f.kw.update(artifact='relative.pkl'),success=False,preserve=True)
        self.simple('parent_target',lambda f:f.kw.update(artifact=str(f.output/'x/../selected.pkl')),success=False,preserve=True)
        self.simple('outside_target',lambda f:f.kw.update(artifact=str(f.root/'outside_target.pkl')),success=False,preserve=True)
        self.simple('relative_root',lambda f:f.kw.update(cache_root='relative'),success=False,preserve=True)
        self.simple('unknown_assembly_run',lambda f:f.kw.update(assembly='unvalidated'),success=True)
        self.simple('pm1_base_rejected',lambda f:f.kw.update(annotation_base=str(f.outside)),kind='pm1',success=False,preserve=True)
    def recovery(self):
        for phase in ['staged_receipt','marker','artifact','receipt']:
            for kind in ['motifs','pm1']:
                self.counter+=1;arms=[]
                for arm in ['reference','native']:
                    f=self.fixture('interrupt_'+phase,arm,kind);f.complete();before=f.pair()
                    if phase=='staged_receipt':f.env['FIXTURE_FAIL_KIND']=phase
                    else:f.env['FIXTURE_FAIL_RENAME']=str({'marker':f.marker(),'artifact':f.target,'receipt':Path(str(f.target)+'.nm_to_enst.json')}[phase])
                    failed=self.invoke(f,'run',arm);failed['normalized_message']=self.normalize(failed.get('message',''),f)
                    f.env.pop('FIXTURE_FAIL_KIND',None);f.env.pop('FIXTURE_FAIL_RENAME',None)
                    calls=(f.root/'argv.json').read_bytes();value=self.invoke(f,'recover',arm)
                    value['normalized_message']=self.normalize(value.get('message',''),f)
                    value['injected_failure']=failed;value['producer_not_repeated']=(f.root/'argv.json').read_bytes()==calls
                    value['marker_removed']=not f.marker().exists();value['prior_preserved']=phase!='staged_receipt' or f.pair()==before
                    arms.append(value)
                self.record('recover_'+kind+'_'+phase,arms,success=True,extra=all(not a['injected_failure']['success'] and a['injected_failure']['category']=='OSError' and a['producer_not_repeated'] and a['marker_removed'] and a['prior_preserved'] for a in arms))
    def paths_and_locks(self):
        def alias(f):
            p=f.root/'cache alias';p.symlink_to(f.cache,target_is_directory=True);f.kw['cache_root']=str(p)
        self.simple('root_alias',alias,success=True)
        self.simple('relative_output_context',lambda f:f.kw.update(output_dir='relative directory'),success=True)
        self.simple('empty_inputs',lambda f:(f.kw.update(input=[]),setattr(f,'expected',json.dumps(f.module.BIOMART.mapping_identity(f.mapping,[])))),success=True)
        self.simple('reordered_duplicate_inputs',lambda f:(f.kw.update(input=list(map(str,[f.inputs[1],f.inputs[0],f.inputs[1]]))),setattr(f,'expected',json.dumps(f.module.BIOMART.mapping_identity(f.mapping,[f.inputs[1],f.inputs[0],f.inputs[1]])))),success=True)
        self.simple('leading_literal_separator',lambda f:f.arguments.insert(0,'--'),success=True)
        self.simple('unknown_tool_option_literal',lambda f:f.arguments.extend(['--unknown-scientific-option','-2']),success=True)
        for which in ['producer','config','target','receipt']:
            def setup(f,which=which):
                f.complete();p={'producer':f.producer,'config':f.config,'target':f.target,'receipt':Path(str(f.target)+'.nm_to_enst.json')}[which]
                kept=f.root/('preserved_'+which);p.rename(kept);p.symlink_to(kept)
            self.simple('linked_'+which,setup,success=False,preserve=True)
        def shared(f):
            f.complete();os.link(f.target,f.root/'hardlinked_history')
        self.simple('shared_target',shared,success=False,preserve=True)
        self.simple('legacy_relative_root_recover',lambda f:(f.complete(),f.kw.update(cache_root='legacy')),command='recover',success=True,preserve=True)
        self.simple('missing_mapping_absent_recover',lambda f:f.mapping.unlink(),command='recover',success=True,preserve=True)
        for label,raw in [('unterminated','{'),('trailing_comma','[1,]'),('bom','\ufeff{}')]:
            self.simple('expected_json_'+label,lambda f,raw=raw:setattr(f,'expected',raw),success=False,preserve=True,message=False)
        self.counter+=1;arms=[]
        for arm in ['reference','native']:
            f=self.fixture('directory_lock',arm);f.complete();before=f.pair();descriptor=os.open(f.output,os.O_RDONLY|os.O_DIRECTORY)
            try:
                fcntl.flock(descriptor,fcntl.LOCK_EX|fcntl.LOCK_NB);value=self.invoke(f,'run',arm)
            finally:os.close(descriptor)
            value['normalized_message']=self.normalize(value.get('message',''),f);value['preserved']=before==f.pair();arms.append(value)
        self.record('nonblocking_directory_lock',arms,success=False,extra=all(v['preserved'] for v in arms))

    def pending_cases(self):
        mutations=['producer_job_bigint','producer_job_surrogate','producer_job_nested','ctime_bigint','ctime_surrogate','float_size','leading_whitespace','root_extra','schema','context_extra_bigint','context_extra_surrogate','prior_extra','prior_list','prior_unhashable','new_list','size_string','size_null','size_zero','size_float_negative','stage_integer','stage_escape','stage_mode','stage_setgid','target_escape','receipt_first','marker_symlink','marker_hardlink','marker_oversize','staged_symlink','receipt_extra_bigint','receipt_extra_surrogate','receipt_extra_nested','receipt_duplicate_last']
        valid={'producer_job_bigint','producer_job_surrogate','producer_job_nested','ctime_bigint','ctime_surrogate','float_size','leading_whitespace','stage_setgid','receipt_duplicate_last'}
        for mutation in mutations:
            self.counter+=1;arms=[]
            for arm in ['reference','native']:
                f=self.fixture(mutation,arm);f.complete();f.env['FIXTURE_FAIL_RENAME']=str(f.marker())
                stopped=self.invoke(f,'run','reference')
                if stopped['success']:raise RuntimeError('source marker fixture did not interrupt')
                f.env.pop('FIXTURE_FAIL_RENAME');value=json.loads(f.marker().read_text());stage=Path(value['stage']);receipt=Path(str(stage/f.target.name)+'.nm_to_enst.json')
                selected=None
                if mutation.startswith('producer_job_'):
                    value['producer_job_id']={'producer_job_bigint':2**150,'producer_job_surrogate':'\ud800','producer_job_nested':[[False,None,1.5]]}[mutation]
                    if mutation=='producer_job_nested':
                        for _ in range(599):value['producer_job_id']=[value['producer_job_id']]
                elif mutation.startswith('ctime_'):value['new']['artifact']['ctime_ns']=2**150 if mutation=='ctime_bigint' else '\ud800'
                elif mutation=='float_size':value['new']['artifact']['size']=float(value['new']['artifact']['size'])
                elif mutation=='root_extra':value['ignored']=None
                elif mutation=='schema':value['schema']='unknown'
                elif mutation.startswith('context_extra_'):value['context']['ignored']=2**150 if mutation.endswith('bigint') else '\ud800'
                elif mutation=='prior_extra':value['prior']['artifact']['ignored']=1
                elif mutation=='prior_list':value['prior']=['artifact','receipt']
                elif mutation=='prior_unhashable':value['prior']=['unknown',[]]
                elif mutation=='new_list':value['new']=['artifact','receipt']
                elif mutation.startswith('size_'):value['new']['artifact']['size']={'size_string':'1','size_null':None,'size_zero':0,'size_float_negative':-0.5}[mutation]
                elif mutation=='stage_integer':value['stage']=7
                elif mutation=='stage_escape':value['stage']=str(f.root)
                elif mutation in ['stage_mode','stage_setgid']:stage.chmod(0o750 if mutation=='stage_mode' else 0o2700)
                elif mutation=='target_escape':value['target']=str(f.outside)
                elif mutation=='receipt_first':os.replace(receipt,Path(str(f.target)+'.nm_to_enst.json'))
                elif mutation=='staged_symlink':
                    path=stage/f.target.name;path.unlink();path.symlink_to(f.outside)
                elif mutation.startswith('receipt_extra_') or mutation=='receipt_duplicate_last':
                    v=json.loads(receipt.read_text())
                    if mutation=='receipt_duplicate_last':raw='{"artifact_bytes":0,'+json.dumps(v)[1:]
                    else:
                        v['ignored']={'receipt_extra_bigint':2**150,'receipt_extra_surrogate':'\ud800','receipt_extra_nested':[False,None,1.5]}[mutation];raw=json.dumps(v)
                    receipt.write_text(raw);value['new']['receipt']=f.module.file_state(receipt)
                write(f.marker(),value)
                if mutation=='producer_job_nested':f.marker().write_text(json.dumps(value,separators=(',',':')))
                if mutation=='leading_whitespace':f.marker().write_text(' \n\t'+f.marker().read_text())
                if mutation in ['marker_symlink','marker_hardlink']:
                    selected=f.root/'preserved_marker.json';f.marker().rename(selected)
                    if mutation=='marker_symlink':f.marker().symlink_to(selected)
                    else:os.link(selected,f.marker())
                if mutation=='marker_oversize':f.marker().write_bytes(b' '*65537)
                before=f.pair();calls=(f.root/'argv.json').read_bytes();result=self.invoke(f,'recover',arm)
                result['normalized_message']=self.normalize(result.get('message',''),f)
                result['no_producer_call']=(f.root/'argv.json').read_bytes()==calls
                result['prior_preserved']=mutation in valid or before==f.pair()
                arms.append(result)
            self.record(mutation,arms,success=mutation in valid,extra=all(v['no_producer_call'] and v['prior_preserved'] for v in arms))

    def context_cases(self):
        self.counter+=1;f=self.fixture('context','same');arms=[self.invoke(f,'context',a) for a in ['reference','native']]
        self.record('context_exact_same_profile',arms,success=True)

    def bound_cases(self):
        alpha=load(self.args.repo/'rust/tests/alpha_annotation_base_compare.py')
        alpha.Harness(SimpleNamespace(repo=self.args.repo,binary=self.args.binary,work_root=self.root/'bound_support'))
        self.counter+=1;f=Fixture(self.root/(str(self.counter)+'_bound_context'),self.args.repo,self.shim,bound=True)
        self.record('bound_context_exact',[self.invoke(f,'context',a) for a in ['reference','native']],success=True)
        for missing in ['alphamissense_annotation_base.py','alphamissense_cache.py']:
            self.counter+=1;arms=[]
            for arm in ['reference','native']:
                f=Fixture(self.root/(str(self.counter)+'_missing_asset_'+arm),self.args.repo,self.shim,bound=True)
                (f.scripts/missing).unlink();value=self.invoke(f,'context',arm);value['normalized_message']=self.normalize(value.get('message',''),f);arms.append(value)
            self.record('bound_missing_'+missing,arms,success=False)
        for mutation in ['none','wide_provenance','false_provenance','schema','historical_artifact','historical_receipt','lineage_symlink']:
            self.counter+=1;arms=[]
            for arm in ['reference','native']:
                f=Fixture(self.root/(str(self.counter)+'_bound_'+mutation+'_'+arm),self.args.repo,self.shim,bound=True)
                source,lineage=f.lineage();v=json.loads(lineage.read_text());v.update(schema='priva-alphamissense-motif-input-v1',provenance={'fixture':'verified retained input'},artifact_path=str(source))
                args=SimpleNamespace(**{k:(Path(x) if k in ['artifact','mapping','cache_root','output_dir','config','producer','annotation_base'] and x is not None else x) for k,x in f.kw.items()});args.input=f.inputs
                with __import__('unittest').mock.patch.dict(os.environ,f.env,clear=True):v['scientific_input']=f.module.context(args)['scientific_input']
                if mutation=='wide_provenance':v['provenance']=2**150
                if mutation=='false_provenance':v['provenance']=False
                if mutation=='schema':v['schema']='wrong'
                write(lineage,v);f.kw['native_lineage']=str(lineage)
                f.env['FIXTURE_FAIL_RENAME']=str(f.target)
                stopped=self.invoke(f,'adopt-native',arm,source,lineage)
                if stopped['success']:raise RuntimeError('bound adoption did not interrupt')
                f.env.pop('FIXTURE_FAIL_RENAME');f.kw['native_lineage']=None
                if mutation=='historical_artifact':source.write_bytes(source.read_bytes()+b'changed')
                if mutation=='historical_receipt':Path(str(source)+'.nm_to_enst.json').write_text('{}')
                if mutation=='lineage_symlink':
                    kept=f.root/'lineage_saved.json';lineage.rename(kept);lineage.symlink_to(kept)
                before=f.pair();result=self.invoke(f,'recover',arm);result['normalized_message']=self.normalize(result.get('message',''),f)
                result['pair_unchanged']=result['success'] or before==f.pair();arms.append(result)
            self.record('bound_recover_'+mutation,arms,success=mutation in ['none','wide_provenance'],extra=all(v['pair_unchanged'] for v in arms))

    def changed_mapping(self):
        for mutation in ['recover','same_mapping','run_pending']:
            self.counter+=1;arms=[]
            for arm in ['reference','native']:
                f=self.fixture('namespace_'+mutation,arm);f.complete();old=f.pair();f.mapping.write_text('new selected mapping bytes\n')
                selected=f.module.BIOMART.mapping_identity(f.mapping,f.inputs)
                target=f.module.BIOMART.changed_mapping_target(identity=selected,cache_root=f.cache,output_dir=f.output,assembly='hg38',kind=f.kind)
                f.kw['artifact']=str(target);f.expected=json.dumps(selected);at=f.arguments.index('--output')+1;f.arguments[at]=str(target)
                f.env['FIXTURE_FAIL_RENAME']=str(target);stopped=self.invoke(f,'run','reference');f.env.pop('FIXTURE_FAIL_RENAME')
                if stopped['success']:raise RuntimeError('namespace preparation did not interrupt')
                calls=(f.root/'argv.json').read_bytes()
                if mutation=='same_mapping':f.complete();old=f.pair()
                if mutation=='run_pending':f.arguments=[]
                else:f.kw['artifact']=str(f.target)
                result=self.invoke(f,'run' if mutation=='run_pending' else 'recover',arm);result['normalized_message']=self.normalize(result.get('message',''),f)
                result['selected_preserved']=f.pair()==old;result['producer_not_repeated']=(f.root/'argv.json').read_bytes()==calls
                if result['success']:f.module.BIOMART.check_receipt(target,f.mapping,f.inputs)
                arms.append(result)
            self.record('namespace_'+mutation,arms,success=mutation!='same_mapping',extra=all(v['selected_preserved'] and v['producer_not_repeated'] for v in arms))

    def adoption_edges(self):
        for mutation in ['pm1','lineage_null','lineage_big_hash','lineage_bad_receipt','same_source','source_symlink','source_hardlink','warm_ignored','warm_big_hash']:
            self.counter+=1;arms=[]
            for arm in ['reference','native']:
                f=self.fixture('adopt_edge_'+mutation,arm,'pm1' if mutation=='pm1' else 'motifs');source,lineage=f.lineage()
                v=json.loads(lineage.read_text())
                if mutation=='lineage_null':v=None
                if mutation=='lineage_big_hash':v['artifact_sha256']=2**150
                if mutation=='lineage_bad_receipt':v['receipt_sha256']='0'*64
                write(lineage,v)
                if mutation=='same_source':f.kw['artifact']=str(source)
                if mutation=='source_symlink':
                    kept=f.root/'kept.pkl';source.rename(kept);source.symlink_to(kept)
                if mutation=='source_hardlink':os.link(source,f.root/'linked.pkl')
                if mutation.startswith('warm_'):
                    first=self.invoke(f,'adopt-native','reference',source,lineage)
                    if not first['success']:raise RuntimeError('warm source fixture did not publish')
                    if mutation=='warm_ignored':
                        receipt=Path(str(f.target)+'.nm_to_enst.json');r=json.loads(receipt.read_text());r['ignored']={'x':2**150,'s':'\ud800','other':[False,None,0.25]};write(receipt,r)
                    else:v['artifact_sha256']=2**150;write(lineage,v)
                before=f.pair();result=self.invoke(f,'adopt-native',arm,source,lineage);result['normalized_message']=self.normalize(result.get('message',''),f)
                result['selected_preserved']=before==f.pair();arms.append(result)
            self.record('adopt_edge_'+mutation,arms,success=mutation=='warm_ignored',extra=all(v['selected_preserved'] for v in arms))

    def adoption(self):
        for extra in ['none','bigint','surrogate','nested']:
            self.counter+=1;arms=[];payloads=[]
            for arm in ['reference','native']:
                f=self.fixture('adopt_'+extra,arm);source,lineage=f.lineage()
                if extra!='none':
                    p=Path(str(source)+'.nm_to_enst.json');value=json.loads(p.read_text());value['ignored']={'bigint':2**150,'surrogate':'\ud800','nested':{'x':[False,None,0.25]}}[extra];write(p,value)
                    v=json.loads(lineage.read_text());v['receipt_sha256']=digest(p);v['ignored']=value['ignored'];write(lineage,v)
                before=[f.module.file_state(source),f.module.file_state(Path(str(source)+'.nm_to_enst.json'))]
                value=self.invoke(f,'adopt-native',arm,source,lineage);value['normalized_message']=self.normalize(value.get('message',''),f)
                value['history_unchanged']=before==[f.module.file_state(source),f.module.file_state(Path(str(source)+'.nm_to_enst.json'))]
                if value['success']:
                    payloads.append((f.target.read_bytes(),Path(str(f.target)+'.nm_to_enst.json').read_bytes()))
                    warm=self.invoke(f,'adopt-native',arm,source,lineage);value['warm_success']=warm['success']
                arms.append(value)
            self.record('adopt_'+extra,arms,success=True,extra=all(v['history_unchanged'] and v.get('warm_success') for v in arms) and len(payloads)==2 and payloads[0]==payloads[1])
    def check(self):
        self.basic();self.recovery();self.paths_and_locks();self.pending_cases();self.context_cases();self.changed_mapping();self.adoption();self.adoption_edges();self.bound_cases()
        report=dict(cases=len(self.rows),passed=sum(r['passed'] for r in self.rows),successes=sum(r['reference']['success'] for r in self.rows),rows=self.rows)
        write(self.root/'complete.json',report)
        print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
        return 0 if all(r['passed'] for r in self.rows) else 1


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    check=sub.add_parser('check');check.add_argument('--repo',type=Path,required=True);check.add_argument('--binary',type=Path,required=True);check.add_argument('--work-root',type=Path,required=True)
    ref=sub.add_parser('reference');ref.add_argument('--request',type=Path,required=True);ref.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();return reference(args) if args.command=='reference' else Harness(args).check()
if __name__=='__main__':raise SystemExit(main())
