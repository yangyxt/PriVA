#!/usr/bin/env python3
"""Exact optional installed-environment comparison, using private fixtures.

fixtures --reference scripts/installed_environment.py --python PYTHON --work DIR
         --report JSON [--candidate PRIVA] [--reuse-report JSON]
recompare --comparison JSON --report JSON
Successful Python arms are immutable oracles; candidate replays verify their
source, fixture specifications, and saved artifacts before reusing them.
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
import time
import yaml

MANAGER = r'''import hashlib,json,os,subprocess,sys
from pathlib import Path
r=Path(os.environ['EXPORT_CASE_ROOT']);s=json.loads((r/'scenario.json').read_text())
with (r/'trace.jsonl').open('a') as f:f.write(json.dumps(['manager',*sys.argv[1:]])+'\n')
name=sys.argv[-1];prefix=None
if '--prefix' in sys.argv:prefix=Path(name);name=s.get('names',{}).get(str(prefix),'cadd-unknown')
action=s.get('action')
out=r/'cache/priva_acmg.installed.yml'
def git(*args):subprocess.run([s['git'],'-C',str(r/'cache'),*args],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
if action=='new-link':out.symlink_to(r/'foreign')
elif action=='replace-object':
 target=out.resolve();p=target.with_name('new');p.write_bytes(target.read_bytes());p.chmod(0o444);p.replace(target)
elif action=='change-object':
 target=out.resolve();target.chmod(0o644);target.write_bytes(b'changed');target.chmod(0o444)
elif action=='retarget':out.unlink();out.symlink_to(r/'foreign')
elif action=='owner':git('config','annex.uuid','22222222-2222-2222-2222-222222222222')
elif action=='unrelated-commit':
 (r/'cache/note').write_text('unrelated metadata\n');git('add','note');git('commit','-qm','unrelated')
elif action=='change-config':(r/'config.yml').write_text('cache_root: relative\n')
elif action=='remove-root':shutil=None;out.unlink(missing_ok=True);(r/'cache').rmdir()
document=s.get('document','name: priva_acmg\nchannels: [conda-forge]\ndependencies: [python=3.12.3=build]\n')
if prefix is not None:
 document='name: null\nprefix: '+json.dumps(str(prefix))+'\nchannels: [conda-forge]\ndependencies: [python=3.12.3=build]\n'
if s.get('fail_name')==name:sys.exit(19)
sys.stdout.buffer.write(document.encode().replace(b'$ROOT',str(r).encode()))
if s.get('stdout_hex') is not None:sys.stdout.buffer.write(bytes.fromhex(s['stdout_hex']))
if s.get('stderr_hex') is not None:sys.stderr.buffer.write(bytes.fromhex(s['stderr_hex']))
sys.exit(s.get('exit',0))
'''
GIT = r'''import json,os,subprocess,sys
from pathlib import Path
r=Path(os.environ['EXPORT_CASE_ROOT']);s=json.loads((r/'scenario.json').read_text())
with (r/'trace.jsonl').open('a') as f:f.write(json.dumps(['git',*sys.argv[1:]])+'\n')
if 'ls-tree' in sys.argv and s.get('git_action'):
 count=sum('ls-tree' in json.loads(line) for line in (r/'trace.jsonl').read_text().splitlines())
 if count==2:
  out=r/'cache/priva_acmg.installed.yml';action=s['git_action']
  if action=='replace-object':
   target=out.resolve();p=target.with_name('new');p.write_bytes(target.read_bytes());p.chmod(0o444);p.replace(target)
  elif action=='retarget':out.unlink();out.symlink_to(r/'foreign')
  elif action=='unrelated-commit':
   (r/'cache/note').write_text('unrelated metadata\n')
   for args in [['add','note'],['commit','-qm','unrelated']]:subprocess.run([s['git'],'-C',str(r/'cache'),*args],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
sys.exit(subprocess.run([s['git'],*sys.argv[1:]]).returncode)
'''
REFERENCE = r'''import json,os,runpy,sys
from pathlib import Path
source=sys.argv.pop(1);last=None
def trace(frame,event,arg):
 global last
 if event=='exception' and frame.f_code.co_filename==source and not isinstance(arg[1],SystemExit):last={'category':arg[0].__name__,'message':str(arg[1])}
 return trace
sys.settrace(trace)
try:runpy.run_path(source,run_name='__main__')
finally:
 sys.settrace(None)
 Path(os.environ['EXPORT_CASE_ROOT'],'caught.json').write_text(json.dumps(last))
'''
DEFAULT = dict(name='priva_acmg',channels=['conda-forge'],dependencies=['python=3.12.3=build'])

def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def pin(path):
    path=Path(path);return dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path))

def save(path,value):
    with Path(path).open('x') as stream:json.dump(value,stream,indent=2,sort_keys=True);stream.write('\n')

def specs():
    cases=[]
    def add(label,**kw):cases.append(dict(case=label,**kw))
    add('new');add('regular-identical',existing='identical');add('regular-changed',existing='changed')
    add('regular-foreign',existing='foreign');add('regular-directory',existing='directory')
    add('regular-invalid-yaml',existing='invalid-yaml');add('regular-wrong-name',existing='wrong-name')
    add('regular-extra-key',existing='extra');add('new-link-during-manager',action='new-link')
    for name in ['', '../bad','a/b','-bad','é','name space','ok.name-1_']:
        add('name-'+str(len(cases)),name=name)
    for name,config in [('none',''),('list','[one]'),('null','cache_root: null'),('relative','cache_root: relative'),('slash','cache_root: /'),('missing','cache_root: $ROOT/missing'),('regular','cache_root: $ROOT/foreign'),('mapping','cache_root: {a: b}'),('alias','root: &a $ROOT/cache\ncache_root: *a'),('double-key','cache_root: relative\ncache_root: $ROOT/cache'),('merge','base: &a {cache_root: $ROOT/cache}\n<<: *a'),('block','cache_root: >-\n  $ROOT/cache'),('dots','cache_root: $ROOT/cache/../cache')]:
        add('config-'+name,config=config)
    add('config-file-missing',config_missing=True);add('config-invalid-utf8',config_hex='ff')
    add('manager-failure',exit=9);add('manager-missing',manager='missing');add('manager-not-executable',manager='not-executable');add('manager-directory',manager='directory')
    add('manager-invalid-stdout',stdout_hex='ff');add('manager-invalid-stderr',stderr_hex='ff');add('manager-failed-invalid-stderr',exit=9,stderr_hex='ff')
    add('manager-empty',document='');add('manager-sequence',document='[one]');add('manager-wrong-name',document=yaml.safe_dump(dict(DEFAULT,name='other')))
    add('manager-null-name',document=yaml.safe_dump(dict(DEFAULT,name=None)));add('manager-no-name',document='dependencies: [python]\n')
    for i,value in enumerate([None,[],{},'package',[None],[''],['  '],['\u001c'],[7],[True],[{'pip':[]}],[{'pip':['a','b']}],[{'pip':'a'}],[{'pip':[3]}],[{'pip':None}],[{'pip':['']}],[{'pip':['a'],'other':'b'}]]):
        add(f'dependencies-{i}',document=yaml.safe_dump(dict(DEFAULT,dependencies=value),sort_keys=False))
    for i,value in enumerate([None,{},'conda-forge',[1],[None],[''],['conda-forge'],['https://example.test/channel'],['https://@example.test'],['https://:@example.test'],['https://user@example.test'],['https://:pass@example.test'],['https://example.test?q=1'],['https://example.test#x'],['https://example.test?'],['https://example.test#'],['file:foo'],['FILE:/foo'],['/local'],[' https://example.test '],['https://[::1]/channel'],['https://[bad]'],['https://[v1.x]/channel'],['https://[vZ.x]'],['https://[127.0.0.1]'],['https://example.test：80/path'],['https://[::1']]):
        add(f'channels-{i}',document=yaml.safe_dump(dict(DEFAULT,channels=value),sort_keys=False))
    for i,value in enumerate(['true','null','123','1.0e+3','2024-01-01','a: b','a # b','- a',"both '\"",'line\nbreak','trailing ','éβ😀','\t','\x00',' '.join(['package']*20)]):
        add(f'quoted-package-{i}',document=yaml.safe_dump(dict(DEFAULT,dependencies=[value,{'pip':[value]}]),sort_keys=False))
    for name,document in [('shared-list','channels: &a [package]\ndependencies: *a'),('shared-pip','dependencies: [&a {pip: [one]}, *a]'),('recursive-list','dependencies: &a [*a]'),('merge-pip','base: &a {pip: [one]}\ndependencies: [{<<: *a}]'),('duplicate-deps','dependencies: []\ndependencies: [one]')]:
        add(name,document=document)
    malformed=['[unterminated','a: [1,] trailing','a:\n\tb: x','a: {x:','a: "bad\\q"','*missing','a: &a x\nb: &a y','[x]: y','x: !!python/object x','a: x\n---\nb: y','a: [x}\n','a: !!int bad','a: !!str [x]','a: !!bool bad','a: !!timestamp bad','a: !!binary @bad']
    for i,value in enumerate(malformed):add(f'malformed-{i}',document=value)
    for name in ['good','relative','missing','directory-history','reported-wrong','reported-null','reported-symlink']:
        add('prefix-'+name,prefix=name)
    for name in ['identical-md5','identical-sha256','changed','manager-failure','uncommitted','foreign','dangling','corrupt','object-symlink','objects-symlink','no-owner','invalid-owner','invalid-schema','large','retarget','change-object','replace-object','owner','unrelated-commit']:
        add('annex-'+name,annex=name,exit=7 if name=='manager-failure' else 0,action=name if name in ['retarget','change-object','replace-object','owner','unrelated-commit'] else None)
    for name in ['none','null','nonstring','empty','missing','two','incomplete','late-failure','directory-recipe','invalid-name','change-config','symlink-base']:
        add('installation-'+name,installation=name,action='change-config' if name=='change-config' else None)
    for name,args in [('no-command',[]),('invalid-command',['other']),('missing-manager',['export','--config','$ROOT/config.yml']),('unknown-option',['export','--config','$ROOT/config.yml','--manager','$ROOT/tools/manager','--unknown']),('help',['--help']),('export-help',['export','--help']),('abbrev',['export','--con','$ROOT/config.yml','--man','$ROOT/tools/manager'])]:
        add('arguments-'+name,args=args)
    for name,document in [
        ('omap-empty-channels','channels: !!omap []\ndependencies: [one]'),
        ('pairs-empty-channels','channels: !!pairs []\ndependencies: [one]'),
        ('omap-dependencies','dependencies: !!omap [{a: b}]'),
        ('pairs-dependencies','dependencies: !!pairs [{a: b}]'),
        ('omap-wrong-node','ignored: !!omap x\nother: !!int bad'),
        ('seq-wrong-node','ignored: !!seq x\nother: !!int bad'),
        ('map-wrong-node','ignored: !!map x\nother: !!int bad'),
        ('set-wrong-node','ignored: !!set x\nother: !!int bad'),
        ('scalar-cycle','dependencies: [!!str &a {=: *a}]'),
        ('omap-wrong-entry','ignored: !!omap [x]'),
        ('omap-multiple-entry','ignored: !!omap [{x: a, y: b}]'),
        ('value-tag','dependencies: [!!str {=: one}]'),
        ('merge-wrong-node','dependencies: [{<<: x}]'),
        ('merge-wrong-sequence','dependencies: [{<<: [x]}]'),
        ('unknown-explicit','ignored: !unknown x\ndependencies: [one]'),
        ('yaml-set','ignored: !!set {a, b}\ndependencies: [one]'),
        ('yaml-set-key','? !!set {a, b}\n: x'),
        ('float-invalid','ignored: !!float bad'),
        ('float-empty','ignored: !!float ""'),
        ('int-empty','ignored: !!int ""'),
        ('binary-valid','ignored: !!binary YQ==\ndependencies: [one]'),
        ('binary-invalid','ignored: !!binary abc'),
        ('binary-invalid-length','ignored: !!binary a'),
        ('binary-unicode','ignored: !!binary é'),
        ('timestamp-valid','ignored: !!timestamp 2024-1-1\ndependencies: [one]'),
        ('timestamp-hour','ignored: 2024-01-01T25:00:00'),
        ('timestamp-month','ignored: 2024-13-01'),
        ('timestamp-day','ignored: 2023-02-29'),
        ('timestamp-offset','ignored: 2024-01-01T00:00:00+24:00'),
        ('timestamp-compound','ignored: 2024-13-01T00:00:00+24:00'),
    ]:add('yaml-edge-'+name,document=document)
    for name in ['key-trailing-dot','utf16','nul-payload','object-loop','staged-foreign','invalid-utf8','utf16-truncated','utf16-surrogate','utf16-low']:
        add('annex-edge-'+name,annex=name)
    for name in ['replace-object','retarget','unrelated-commit']:
        add('annex-git-race-'+name,annex='identical-md5',git_action=name)
    return cases


def identity(path):
    path=Path(path)
    try:s=path.lstat()
    except FileNotFoundError:return None
    result=dict(mode=stat.S_IMODE(s.st_mode),kind='link' if path.is_symlink() else 'file' if path.is_file() else 'directory',stat=[s.st_dev,s.st_ino,s.st_mode,s.st_size,s.st_mtime_ns,s.st_ctime_ns,s.st_nlink])
    if path.is_symlink():result['literal']=os.readlink(path)
    elif path.is_file():result.update(bytes=s.st_size,sha256=sha(path))
    return result

def files(root):
    result={}
    for folder,dirs,names in os.walk(root,followlinks=False):
        for name in dirs+names:
            path=Path(folder)/name;relative=str(path.relative_to(root))
            if relative.startswith('cache/') or relative=='foreign':
                item=identity(path)
                if item['kind']!='directory':result[relative]=item
    return result

def setup(folder,spec,python):
    folder.mkdir();(folder/'tools').mkdir();cache=folder/'cache';cache.mkdir();(folder/'foreign').write_text('foreign untouched\n')
    env=dict(os.environ,EXPORT_CASE_ROOT=str(folder),PYTHONDONTWRITEBYTECODE='1',GIT_CONFIG_NOSYSTEM='1',GIT_CONFIG_GLOBAL=os.devnull,
             GIT_AUTHOR_DATE='2001-01-01T00:00:00+00:00',GIT_COMMITTER_DATE='2001-01-01T00:00:00+00:00')
    for key in ('GIT_DIR','GIT_COMMON_DIR','GIT_WORK_TREE','GIT_INDEX_FILE'):env.pop(key,None)
    git=shutil.which('git');scenario=dict(spec,git=git)
    for name,program in [('manager',MANAGER),('git',GIT)]:
        path=folder/'tools'/name;path.write_text('#!'+str(python)+'\n'+program);path.chmod(0o755)
    env['PATH']=str(folder/'tools')+os.pathsep+env['PATH']
    config=spec.get('config','cache_root: $ROOT/cache\n').replace('$ROOT',str(folder))
    if not spec.get('config_missing'):(folder/'config.yml').write_bytes(bytes.fromhex(spec['config_hex']) if 'config_hex' in spec else config.encode())
    output=cache/'priva_acmg.installed.yml';payload=yaml.safe_dump(DEFAULT,sort_keys=False).encode()
    existing=spec.get('existing')
    if existing=='directory':output.mkdir()
    elif existing:output.write_bytes(dict(identical=payload,changed=yaml.safe_dump(dict(DEFAULT,dependencies=['old']),sort_keys=False).encode(),foreign=b'foreign',**{'invalid-yaml':b'[bad','wrong-name':yaml.safe_dump(dict(DEFAULT,name='other')).encode(),'extra':yaml.safe_dump(dict(DEFAULT,variables={})).encode()})[existing])
    args=['export','--config',str(folder/'config.yml'),'--manager',str(folder/'tools/manager')]
    if 'name' in spec:args+=['--name',spec['name']]
    if spec.get('manager')=='missing':args[-1]=str(folder/'missing')
    if spec.get('manager')=='not-executable':(folder/'tools/manager').chmod(0o644)
    if spec.get('manager')=='directory':args[-1]=str(cache)
    prefix_kind=spec.get('prefix')
    if prefix_kind:
        prefix=folder/'prefix';prefix.mkdir();(prefix/'conda-meta').mkdir()
        if prefix_kind=='directory-history':(prefix/'conda-meta/history').mkdir()
        elif prefix_kind!='missing':(prefix/'conda-meta/history').write_text('history')
        args+=['--prefix',str(prefix) if prefix_kind!='relative' else 'prefix']
        scenario['document']=yaml.safe_dump(dict(DEFAULT,prefix=str(prefix)),sort_keys=False)
        # Select custom explicit-prefix replies instead of automatic CADD replies.
        manager=(folder/'tools/manager').read_text().replace("if prefix is not None:","if prefix is not None and False:")
        (folder/'tools/manager').write_text(manager)
        if prefix_kind=='reported-wrong':scenario['document']=yaml.safe_dump(dict(DEFAULT,prefix=str(folder/'other')),sort_keys=False)
        if prefix_kind=='reported-null':scenario['document']=yaml.safe_dump(dict(DEFAULT,name=None,prefix=str(prefix)),sort_keys=False)
        if prefix_kind=='reported-symlink':(folder/'prefix-link').symlink_to(prefix);scenario['document']=yaml.safe_dump(dict(DEFAULT,prefix=str(folder/'prefix-link')),sort_keys=False)
    annex=spec.get('annex')
    if annex:
        def run(*args):subprocess.run([git,'-C',str(cache),*args],env=env,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
        run('init','-q');run('config','user.name','Export Fixture');run('config','user.email','fixture@example.invalid');run('config','core.hooksPath',os.devnull);run('config','commit.gpgSign','false');run('config','annex.uuid','11111111-1111-1111-1111-111111111111')
        if annex=='invalid-schema':payload=b'foreign: export\n'
        if annex=='utf16':payload=payload.decode().encode('utf-16')
        if annex=='invalid-utf8':payload=b'\xff'
        if annex=='utf16-truncated':payload=b'\xff\xfeA'
        if annex=='utf16-surrogate':payload=b'\xff\xfe\x00\xd8a\x00'
        if annex=='utf16-low':payload=b'\xff\xfe\x00\xdc'
        if annex=='nul-payload':payload=b'name: priva_acmg\x00' 
        if annex=='large':payload=b'X'*(1024**2+1)
        algorithm='SHA256E' if annex=='identical-sha256' else 'MD5E';digest=hashlib.sha256(payload).hexdigest() if algorithm=='SHA256E' else hashlib.md5(payload).hexdigest();key=f'{algorithm}-s{len(payload)}--{digest}.'+('' if annex=='key-trailing-dot' else 'yml')
        obj=cache/'.git/annex/objects/a/b'/key/key;obj.parent.mkdir(parents=True);obj.write_bytes(payload);obj.chmod(0o444)
        output.symlink_to(os.path.relpath(obj,cache));run('add',output.name);run('commit','-qm','export fixture')
        if annex=='changed':scenario['document']=yaml.safe_dump(dict(DEFAULT,dependencies=['different']),sort_keys=False)
        if annex=='uncommitted':run('rm','--cached','-q',output.name);run('commit','-qm','remove recorded pointer')
        if annex=='foreign':output.unlink();output.symlink_to(folder/'foreign')
        if annex=='dangling':obj.unlink()
        if annex=='corrupt':obj.chmod(0o644);obj.write_bytes(b'corrupt');obj.chmod(0o444)
        if annex=='object-symlink':obj.unlink();obj.symlink_to(folder/'foreign')
        if annex=='object-loop':obj.unlink();obj.symlink_to(obj.name)
        if annex=='staged-foreign':
            literal=os.readlink(output);output.unlink();output.symlink_to('uncommitted');run('add',output.name);output.unlink();output.symlink_to(literal)
        if annex=='objects-symlink':objects=cache/'.git/annex/objects';objects.rename(cache/'.git/annex/retained-objects');objects.symlink_to('retained-objects')
        if annex=='no-owner':run('config','--unset','annex.uuid')
        if annex=='invalid-owner':run('config','annex.uuid','invalid')
    installation=spec.get('installation')
    if installation:
        args[0]='export-installation';base=folder/'cadd';recipes=base/'envs'
        if installation not in ['none','null','nonstring','empty','missing']:
            recipes.mkdir(parents=True)
            for name in ['b','a']:
                if installation=='invalid-name' and name=='b':name='bad name'
                recipe=recipes/(name+'.yml');recipe.write_text('dependencies: [python]\n')
                digest=hashlib.md5(os.path.realpath(recipes/'conda').encode()+recipe.read_bytes()).hexdigest();prefix=recipes/'conda'/(digest+'_')
                # Recipes intentionally have distinct bytes for independent prefixes.
                if name=='b':recipe.write_text('dependencies: [perl]\n');digest=hashlib.md5(os.path.realpath(recipes/'conda').encode()+recipe.read_bytes()).hexdigest();prefix=recipes/'conda'/(digest+'_')
                (prefix/'conda-meta').mkdir(parents=True,exist_ok=True);(prefix/'conda-meta/history').write_text('history')
                if installation!='incomplete' or name=='a':Path(str(prefix)+'.env_setup_done').write_text('complete')
                scenario.setdefault('names',{})[str(prefix)]='cadd-'+name
            if installation=='late-failure':scenario['fail_name']='cadd-b'
            if installation=='directory-recipe':(recipes/'z.yml').mkdir()
            if installation=='symlink-base':(folder/'cadd-link').symlink_to(base);base=folder/'cadd-link'
        base_value={'none':None,'null':'null','nonstring':7,'empty':''}.get(installation,str(base))
        if installation!='none':(folder/'config.yml').write_text(config+'cadd_base_dir: '+json.dumps(base_value)+'\n')
    if 'args' in spec:args=[s.replace('$ROOT',str(folder)) for s in spec['args']]
    (folder/'scenario.json').write_text(json.dumps(scenario));return args,env

def normalize(value,folder):
    if isinstance(value,str):return value.replace(str(folder),'$ROOT')
    if isinstance(value,list):return [normalize(v,folder) for v in value]
    if isinstance(value,dict):return {k:normalize(v,folder) for k,v in value.items()}
    return value

def execute(folder,spec,args):
    argv,env=setup(folder,spec,args.python);before=files(folder)
    if args.arm=='python':
        wrapper=folder/'reference.py';wrapper.write_text(REFERENCE);command=[str(args.python),str(wrapper),str(args.reference),*argv]
    else:command=[str(args.candidate),'installed-environment',*argv]
    start=time.monotonic()
    with (folder/'stdout').open('xb') as out,(folder/'stderr').open('xb') as err:
        child=subprocess.Popen(command,cwd=folder,env=env,stdout=out,stderr=err);_,status,usage=os.wait4(child.pid,0);child.returncode=os.waitstatus_to_exitcode(status)
    after=files(folder);stderr=(folder/'stderr').read_text(errors='backslashreplace');diagnostic=None
    if args.arm=='python':diagnostic=json.loads((folder/'caught.json').read_text())
    elif child.returncode==1:
        match=re.match(r'priva: ([A-Za-z]+(?:Error|Expired)): ([\s\S]*)\n$',stderr)
        diagnostic=dict(category=match[1],message=match[2]) if match else dict(category='UNCLASSIFIED',message=stderr)
    # A caught internal exception is not a terminal failure on a success path.
    if child.returncode==0:diagnostic=None
    trace=[json.loads(v) for v in (folder/'trace.jsonl').read_text().splitlines()] if (folder/'trace.jsonl').exists() else []
    def comparable(items):return {k:{a:b for a,b in v.items() if a!='stat'} for k,v in items.items()}
    preserved={k:after.get(k)==v for k,v in before.items()}
    result=dict(returncode=child.returncode,diagnostic=diagnostic,stdout=(folder/'stdout').read_text(errors='backslashreplace'),trace=trace,files=comparable(after),preserved=preserved)
    result=normalize(result,folder)
    value=dict(case=spec['case'],spec=spec,arm=args.arm,command=command,result=result,before=before,after=after,
               stdout=pin(folder/'stdout'),stderr=pin(folder/'stderr'),wall_seconds=time.monotonic()-start,peak_rss_kib=usage.ru_maxrss)
    save(folder/'result.json',value);return value

def comparison_view(arm):
    """Bind setup-specific metadata before comparing observable results.

    Git indexes encode each fixture's inode/timestamps. Require the saved raw
    bytes, compare staged modes/blob IDs/path order, and retain each arm's
    independent raw preservation check. CADD names encode its absolute recipe
    root; verify that exact MD5 before comparing the named recipe invocation.
    """
    result=json.loads(json.dumps(arm['result']));root=Path(arm['stdout']['path']).parent;proofs={}
    if arm['arm']=='python' and result['returncode']==1 and result['diagnostic'] is None:
        # CPython disables a trace callback when its recursion limit is reached.
        # Recover only a complete, explicit terminal traceback diagnostic from
        # the preserved stderr; a missing category is never a matched failure.
        assert sha(arm['stderr']['path'])==arm['stderr']['sha256']
        stderr=Path(arm['stderr']['path']).read_text()
        terminal=re.search(r'(?:^|\n)([A-Za-z_][A-Za-z0-9_.]*): ([^\n]*)\n\Z',stderr)
        if 'Traceback (most recent call last):\n' in stderr and terminal:
            result['diagnostic']=normalize(dict(category=terminal[1],message=terminal[2]),root)
            proofs['uncaught_python_diagnostic']=dict(stderr_sha256=arm['stderr']['sha256'],terminal=terminal[0].strip('\n'))
    if result['returncode']==2:
        assert sha(arm['stderr']['path'])==arm['stderr']['sha256']
        stderr=Path(arm['stderr']['path']).read_text()
        if arm['arm']=='python':
            message=next((line.split(': error: ',1)[1] for line in stderr.splitlines() if ': error: ' in line),None)
        else:message=stderr.removeprefix('priva: ').rstrip('\n')
        result['diagnostic']=dict(category='ArgumentError',message=message)

    index_name='cache/.git/index'
    if index_name in result['files']:
        index=root/index_name;assert sha(index)==arm['after'][index_name]['sha256']
        scenario=json.loads((root/'scenario.json').read_text())
        env=dict(os.environ,GIT_OPTIONAL_LOCKS='0',GIT_CONFIG_NOSYSTEM='1',GIT_CONFIG_GLOBAL=os.devnull)
        for key in ('GIT_DIR','GIT_COMMON_DIR','GIT_WORK_TREE','GIT_INDEX_FILE'):env.pop(key,None)
        staged=subprocess.check_output([scenario['git'],'-C',str(root/'cache'),'--literal-pathspecs','ls-files','--stage','-z'],env=env)
        value=result['files'][index_name];value.pop('sha256');value['staged_entries']=staged.decode()
        proofs['git_index']=dict(raw_sha256=sha(index),staged_sha256=hashlib.sha256(staged).hexdigest(),raw_unchanged=arm['result']['preserved'][index_name])
    recipes=root/'cadd/envs';prefixes=[]
    if recipes.is_dir():
        mapping={}
        for recipe in recipes.glob('*.yml'):
            if not recipe.is_file():continue
            digest=hashlib.md5(os.path.realpath(recipes/'conda').encode()+recipe.read_bytes()).hexdigest()+'_'
            mapping[digest]=recipe.stem
        for row in result['trace']:
            if row[:3]==['manager','env','export'] and '--prefix' in row:
                selected=row[-1];name=Path(selected).name
                if name in mapping:
                    prefixes.append(dict(selected=selected,recipe=mapping[name],exact_root_and_recipe_md5=True))
                    row[-1]='$CADD_PREFIX['+mapping[name]+']'
        proofs['cadd_prefixes']=prefixes
    return result,proofs

def verify_saved_arm(arm):
    for stream in ['stdout','stderr']:
        assert sha(arm[stream]['path'])==arm[stream]['sha256']
    for name,meta in arm['after'].items():
        if meta.get('sha256'):
            assert sha(Path(arm['stdout']['path']).parent/name)==meta['sha256']

def recompare(args):
    saved=json.loads(args.comparison.read_text())
    assert saved['fixture_specs']==specs()
    for key in ['reference','python','candidate']:
        assert saved[key] and sha(saved[key]['path'])==saved[key]['sha256']
    rows=[]
    for old in saved['cases']:
        for arm in ['python','native']:verify_saved_arm(old[arm])
        left,left_proof=comparison_view(old['python']);right,right_proof=comparison_view(old['native'])
        differences=[key for key in left if left[key]!=right[key]]
        rows.append(dict(case=old['case'],python=old['python'],native=old['native'],
            comparison_proofs=dict(python=left_proof,native=right_proof),exact=not differences,differences=differences))
    result=dict(saved,cases=rows,exact=sum(row['exact'] for row in rows),strict_passed=all(row['exact'] for row in rows),
        comparator=pin(__file__),reinterpreted_saved_report=pin(args.comparison),producer_rerun=False)
    save(args.report,result)
    print(json.dumps(dict(count=result['count'],exact=result['exact'],mismatches=[dict(case=row['case'],fields=row['differences']) for row in rows if not row['exact']])))

def fixtures(args):
    args.work.mkdir();saved=json.loads(args.reuse_report.read_text()) if args.reuse_report else None
    old={}
    if saved:
        assert saved['reference']['sha256']==sha(args.reference)
        wanted={v['case']:v for v in specs()}
        assert all(wanted[v['case']]==v for v in saved['fixture_specs'])
        old={row['case']:row['python'] for row in saved['cases']}
    partial={}
    if args.reuse_partial:
        assert sha(args.reuse_partial.parent/'installed_environment.py')==sha(args.reference)
        partial={json.loads(p.read_text())['case']:json.loads(p.read_text()) for p in args.reuse_partial.glob('*/python/result.json')}
    assert len({s['case'] for s in specs()})==len(specs()),'Duplicate fixture name'
    rows=[]
    for spec in specs():
        folder=args.work/spec['case'];folder.mkdir();arms={}
        if spec['case'] in old or spec['case'] in partial:
            arms['python']=old[spec['case']] if spec['case'] in old else partial[spec['case']]
            assert arms['python']['spec']==spec
            verify_saved_arm(arms['python'])
        else:
            args.arm='python';arms['python']=execute(folder/'python',spec,args)
        if args.candidate:
            args.arm='native';arms['native']=execute(folder/'native',spec,args)
            left,left_proof=comparison_view(arms['python']);right,right_proof=comparison_view(arms['native'])
            differences=[key for key in left if left[key]!=right[key]]
            arms['comparison_proofs']=dict(python=left_proof,native=right_proof)
        else:differences=[]
        rows.append(dict(case=spec['case'],**arms,exact=not differences,differences=differences))
    result=dict(reference=pin(args.reference),python=pin(args.python),candidate=pin(args.candidate) if args.candidate else None,
                fixture_specs=specs(),cases=rows,exact=sum(row['exact'] for row in rows),count=len(rows),reused_python=bool(saved),
                strict_passed=all(row['exact'] for row in rows),runtime_stderr_policy='Compare exact terminal exception class/message; native prefix and Python tracebacks remain separately retained raw bytes.')
    save(args.report,result);print(json.dumps(dict(count=result['count'],exact=result['exact'],mismatches=[dict(case=row['case'],fields=row['differences']) for row in rows if not row['exact']])))

def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(required=True)
    p=sub.add_parser('fixtures');p.add_argument('--reference',type=Path,required=True);p.add_argument('--python',type=Path,required=True);p.add_argument('--candidate',type=Path);p.add_argument('--work',type=Path,required=True);p.add_argument('--report',type=Path,required=True);p.add_argument('--reuse-report',type=Path);p.add_argument('--reuse-partial',type=Path);p.set_defaults(function=fixtures)
    p=sub.add_parser('recompare');p.add_argument('--comparison',type=Path,required=True);p.add_argument('--report',type=Path,required=True);p.set_defaults(function=recompare)
    args=parser.parse_args();args.function(args)
if __name__=='__main__':main()
