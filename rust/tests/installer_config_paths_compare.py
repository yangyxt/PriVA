#!/usr/bin/env python3
"""Differential installer path plans and private byte-adoption fixtures.

python installer_config_paths_compare.py check --repo CHECKOUT --binary PRIVA --work-root EMPTY
python installer_config_paths_compare.py reference --request JSON --output JSON

Each arm uses an independent private tree. Scientific resources are opaque byte
fixtures except explicitly labelled bounded indexed-header fixtures. Never run
this harness against an installed cache/configuration. Full source outcomes and
all failed/interrupted fixtures are retained, not silently replaced.
"""
from __future__ import annotations
import argparse, contextlib, fcntl, hashlib, importlib.util, io, json, os, re, shutil, stat, subprocess, sys, time
from pathlib import Path
from types import SimpleNamespace

def write(path,value):Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=True)+'\n')
def load(path):
 spec=importlib.util.spec_from_file_location('config_paths_'+str(abs(hash(str(path)))),path)
 m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def export(value):
 if isinstance(value,Path):return str(value)
 if isinstance(value,dict):return {str(k):export(v) for k,v in value.items()}
 if isinstance(value,(tuple,list)):return list(map(export,value))
 if isinstance(value,set):return sorted(map(export,value))
 return value
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def xattrs(path):return {name:os.getxattr(path,name).hex() for name in sorted(os.listxattr(path))}
def reference(args):
 request=json.loads(args.request.read_text());m=load(Path(request['helper'])/'scripts/installer_config_paths.py');options=request['options'];command=request['command']
 output=io.StringIO()
 try:
  with contextlib.redirect_stdout(output):
   path=Path(options['config']) if 'config' in options else None
   root=Path(options['source-root']) if 'source-root' in options and 'no-source-root' not in options else None
   if command=='resolve':
    value=m.resolve_config(SimpleNamespace(config=path,source_root=root,check='check' in options,output=Path(options['output']) if 'output' in options else None,managed_updates_json=options.get('managed-updates-json')))
   elif command=='render':value=m.render_config(path.read_text(),json.loads(options['updates-json']));output.write(value);value=None
   elif command=='adopt-sources':value=m.adopt_sources(json.loads(options['plan-json']),Path(options['cache-root']))
   elif command=='hash-order':value=dict(literal=list(m.REPLACEMENT_TABLE_KEYS),difference=list(m.REPLACEMENT_TABLE_KEYS-{'hgnc_table'}))
   elif command in ['cache-identity','sha256']:value=getattr(m,'cache_identity' if command=='cache-identity' else 'sha256_file')(Path(options['path']))
   elif command=='source-format':value=m.spliceai_source_format(options['key'],options['path'])
   elif command=='reverse-chain-source':value=m.reverse_chain_source(json.loads(options['value-json']),root)
   elif command=='bundled-resource':value=m.bundled_resource(options['key'],json.loads(options['value-json']),root)
   elif command in ['gnomad-members','gnomad-validate']:
    g=m.gnomad_inspector()
    value=g.family_paths(Path(options['path'])) if command=='gnomad-members' else g.validate_family(chr_x=Path(options['path']),assembly=options['assembly'],reference=options.get('reference'))
   else:
    config=m.yaml.load(path.read_text(),Loader=m.UniqueLoader)
    if command=='managed-roots':value=m.configure_cache_paths(path,config,'--check' if 'check' in options else 'write',root)
    elif command=='portable-updates':value=m.portable_updates(config,root)
    elif command=='resolve-cache-paths':value=m.resolve_cache_paths(config,root)
    elif command=='defaults':value=m._portable_defaults(config,root)
    elif command=='resource-groups':value=m.resource_groups(config,json.loads(options['updates-json']))
    elif command=='adoption-plan':value=m.adoption_plan(config,json.loads(options['updates-json']) if 'updates-json' in options else m.portable_updates(config,root),root)
    elif command=='source-needed':value=m.source_is_needed(options['key'],config)
    else:raise AssertionError(command)
  outcome=dict(success=True,value=export(value),stdout=output.getvalue())
 except Exception as e:outcome=dict(success=False,category=type(e).__name__,message=str(e),stdout=output.getvalue())
 write(args.output,outcome)

SHIM=r'''#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/file.h>
int flock(int fd,int operation){
 static int (*next)(int,int)=0;if(!next)next=dlsym(RTLD_NEXT,"flock");
 const char *path=getenv("FIXTURE_FLOCK_ENTERED");
 if(path&&(operation&LOCK_EX)){FILE *f=fopen(path,"a");if(f){fputs("entered\n",f);fclose(f);}}
 return next(fd,operation);
}
int rename(const char *a,const char *b){
 static int (*next)(const char*,const char*)=0;if(!next)next=dlsym(RTLD_NEXT,"rename");
 const char *needle=getenv("FIXTURE_RENAME_FAIL");int before=getenv("FIXTURE_FAIL_BEFORE")!=0;
 const char *base=strrchr(b,'/');base=base?base+1:b;
 int hit=needle&&!strcmp(base,needle);int rc;
 if(hit&&before){errno=EIO;rc=-1;}else{rc=next(a,b);if(!rc&&hit){errno=EIO;rc=-1;}}
 int saved=errno;const char *log=getenv("FIXTURE_RENAME_LOG");if(log){FILE *f=fopen(log,"a");if(f){fprintf(f,"%d\t%s\t%s\n",rc,a,b);fclose(f);}}
 const char *change=getenv("FIXTURE_CHANGE_AFTER_RENAME");if(change&&!rc&&strstr(a,".adopting-")){FILE *f=fopen(change,"a");if(f){fputs("# concurrent fixture writer\n",f);fclose(f);}unsetenv("FIXTURE_CHANGE_AFTER_RENAME");}
 errno=saved;return rc;
}
'''

class Fixture:
 def __init__(self,root,repo,shim,seed=0):
  self.root=root;root.mkdir(parents=True);self.repo=repo;self.source=root/'source';(self.source/'scripts').mkdir(parents=True)
  for name in ['installer_config_paths.py','spliceai_cache_input.py','gnomad_cache.py']:
   shutil.copy2(repo/'scripts'/name,self.source/'scripts'/name)
  self.cache=root/'cache';self.config=root/'config.yaml';self.config_data=dict(cache_root=str(self.cache),tmp_dir=str(root/'tmp'),assembly='hg38',vep_installed_version='116')
  self.env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONHASHSEED=str(seed),LD_PRELOAD=str(shim),FIXTURE_RENAME_LOG=str(root/'renames.log'))
  for key in list(self.env):
   if key.startswith('PRIVA_FORCE_'):self.env.pop(key)
  self.options=dict(config=str(self.config),**{'source-root':str(self.source),'helper-root':str(self.source)})
  self.command='resolve';self.raw=None;self.sources=[];self.extra=[];self.hash_aliases={};self.protected=[];self.attempt=0;self.expect_source_unchanged=True
 def file(self,relative,data=b'opaque selected source\x00B A B\n',mode=0o640):
  p=self.root/relative;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data);p.chmod(mode);os.utime(p,ns=(1700000000000000000,1700000000000000000));return p
 def selected(self,key,relative='outside/selected.dat',data=b'opaque selected bytes\n'):
  p=self.file(relative,data);self.config_data[key]=str(p);self.sources.append(p);return p
 def complete_bundles(self):
  m=load(self.source/'scripts/installer_config_paths.py')
  for relative in [*m.BUNDLED_PATHS.values(),*m.BUNDLED_CACHES.values(),'data/loeuf/loeuf_dataset_hg38.tsv.gz']:
   p=self.source/relative
   if not p.exists():self.sources.append(self.file('source/'+relative,b'opaque bundled snapshot; existence/nonzero admission only\n'))
 def save(self):
  import yaml
  self.config.write_text(self.raw if self.raw is not None else yaml.safe_dump(self.config_data,sort_keys=False))
  self.config.chmod(0o640);os.utime(self.config,ns=(1700000000000000000,1700000000000000000))
 def normalize(self,value):
  if isinstance(value,str):
   for actual,alias in self.hash_aliases.items():value=value.replace(actual,alias)
   value=value.replace(str(self.root),'<FIXTURE>')
   value=re.sub(r'(?<=\.installer-retained/)\d{8}T\d{6}\.\d{6}Z','<STAMP>',value)
   value=re.sub(r'\.adopting-[^/\s\'"\],)]+','.adopting-<TEMP>',value)
   value=re.sub(r'\.adoptions-[^/\s\'"\],)]+','.adoptions-<TEMP>',value)
   value=re.sub(r'\.adoption-metadata-[^/\s\'"\],)]+','.adoption-metadata-<TEMP>',value)
   value=re.sub(r'\.config\.yaml\.[^/\s\'"\],)]+','.config.yaml.<TEMP>',value)
   value=re.sub(r'\.output\.yaml\.[^/\s\'"\],)]+','.output.yaml.<TEMP>',value)
   return value
  if isinstance(value,list):return [self.normalize(v) for v in value]
  if isinstance(value,dict):return {self.normalize(k):self.normalize(v) for k,v in value.items()}
  return value
 def source_state(self):
  values={}
  for p in self.sources:
   m=p.stat();values[str(p)]=dict(identity=[m.st_ino,m.st_size,m.st_mtime_ns,m.st_ctime_ns],sha256=digest(p),mode=stat.S_IMODE(m.st_mode),nlink=m.st_nlink,xattrs=xattrs(p))
  return values
 def products(self):
  values={};errors=[]
  if hasattr(self,'prepared') and (self.prepared/'adapter.ready.json').is_file():
   sha=digest(self.prepared/'adapter.ready.json');self.hash_aliases[sha]='<ORIGINAL_RECEIPT_SHA>';self.hash_aliases[sha[:16]]='<ORIGINAL_RECEIPT_PREFIX>'
   target=self.cache/'VEP_plugin_cache/SpliceAI2/GRCh38'
   if (target/'adapter.ready.json').is_file():
    try:
     selected=json.loads((target/'adapter.ready.json').read_text())
     if 'relocation' in selected:
      observed=digest(target/'adapter.ready.json');record=json.loads((target/'relocation.ready.json').read_text())
      if record['relocated_receipt_sha256']!=observed or record['original_receipt_sha256']!=sha:errors.append('prepared receipt hash binding differs')
      self.hash_aliases[observed]='<RELOCATED_RECEIPT_SHA>'
    except (ValueError,KeyError):pass
  roots=[self.cache,self.config,self.root/'output.yaml']
  def visit(p):
   if not os.path.lexists(p):return
   relative=self.normalize(str(p.relative_to(self.root)));m=p.lstat()
   if p.is_symlink():values[relative]=dict(kind='link',target=self.normalize(os.readlink(p)));return
   if p.is_dir():
    values[relative]=dict(kind='directory',mode=stat.S_IMODE(m.st_mode),xattrs=xattrs(p))
    for child in p.iterdir():visit(child)
    return
   data=p.read_bytes();record=dict(kind='file',mode=stat.S_IMODE(m.st_mode),xattrs=xattrs(p),data=self.normalize(data.decode('utf-8',errors='backslashreplace')))
   # A failed receipt publication retains its complete .adoptions-* temporary
   # receipt. Verify its identities exactly like the named receipt's.
   if p.name in ['installer-source-adoptions.json','relocation.ready.json'] or (p.parent==self.cache and p.name.startswith('.adoptions-')):
    try:v=json.loads(data)
    except (ValueError,UnicodeError):v=None
    items=v.get('copies',v.get('members',[])) if isinstance(v,dict) else []
    if not isinstance(items,list):items=[]
    identities=[]
    for item in items:
     if not isinstance(item,dict) or 'source_identity' not in item:continue
     source=Path(item['source']);now=source.stat();observed=[now.st_ino,now.st_size,now.st_mtime_ns,now.st_ctime_ns]
     if len(item['source_identity'])==5:observed.insert(0,now.st_dev)
     if item['source_identity']!=observed:errors.append('adoption source identity differs')
     if 'sha256' in item and (item['sha256']!=digest(source) or item['bytes']!=now.st_size):errors.append('adoption source bytes/hash differ')
     identities.append(item['source_identity'])
     if 'copied_identity' in item:
      target=Path(item['destination']);now=target.stat();observed=[now.st_dev,now.st_ino,now.st_size,now.st_mtime_ns,now.st_ctime_ns]
      if item['copied_identity']!=observed:errors.append('prepared copied identity differs')
      identities.append(item['copied_identity'])
    # Replace only independently checked inode/time arrays. Preserve every
    # other original output byte, including key order and float spellings.
    def identity_text(match):
     try:observed=json.loads(match.group(2))
     except ValueError:return match.group(0)
     return match.group(1)+'["<actual-source-identity>"]' if observed in identities else match.group(0)
    record['data']=self.normalize(re.sub(r'("(?:source_identity|copied_identity)"\s*:\s*)(\[[^\]]*\])',identity_text,data.decode('utf-8')))
   values[relative]=record
  for p in roots:visit(p)
  events=(self.root/'renames.log').read_text(errors='surrogateescape').splitlines() if (self.root/'renames.log').exists() else []
  return dict(paths=values,events=self.normalize(events),receipt_checks=errors)

 def prepared_pair(self,where='external',variant=None):
  import gzip
  m=load(self.source/'scripts/spliceai_cache_input.py');folder=self.root/where;folder.mkdir(parents=True,exist_ok=True);self.prepared=folder;outputs={}
  for role in ['snv','indel']:
   text='##fileformat=VCFv4.2\n##assembly=GRCh38\n##SpliceAI2Version=2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
   p=self.file(str((folder/f'SpliceAI2_GRCh38.{role}.vep.vcf.gz').relative_to(self.root)),gzip.compress(text.encode(),mtime=0));index=self.file(str(p.relative_to(self.root))+'.tbi',b'synthetic metadata-only index')
   self.sources.extend([p,index]);outputs[role]=dict(path=str(p),index_path=str(index),bytes=p.stat().st_size,sha256=digest(p),index_sha256=digest(index));self.config_data['spliceai_'+role+'_prescore']=str(p)
  context=self.file(str(folder.relative_to(self.root))+'/SpliceAI2_GRCh38.gene_context.tsv',b'chrom\tstart\tend\tstrand\tgene\nchr1\t1\t100\t+\tGENE\n');self.sources.append(context)
  data=dict(schema_version=1,resource='SpliceAI2',version='2.0',assembly='GRCh38',status='prepared_annotation_only',coverage='complete_published',scientific_qualified=False,source_revision=m.REVISION,source_manifest_sha256='a'*64,outputs=outputs,gene_context=dict(path=str(context),sha256=digest(context)))
  if variant=='incomplete':data['coverage']='validation_subset'
  receipt=self.file(str(folder.relative_to(self.root))+'/adapter.ready.json',json.dumps(data).encode());self.sources.append(receipt);self.file(str(folder.relative_to(self.root))+'/.install.lock',b'')
  if variant=='missing':receipt.unlink();self.sources.remove(receipt)

 def gnomad(self,assembly='hg38',variant=None):
  import pysam
  g=load(self.source/'scripts/gnomad_cache.py');family=g.family_paths(self.root/'gnomad/selected.chrX.vcf.bgz');lengths={c:100000+i for i,c in enumerate(g.CHROMOSOMES)};lengths['chr1']=g.CHR1_LENGTH[assembly]
  for chrom,path in family.items():
   text='##fileformat=VCFv4.2\n'+''.join(f'##contig=<ID={c},length={n+(1 if variant=="wrong_length" and c==chrom=="chr7" else 0)}>\n' for c,n in lengths.items())
   for field in g.REQUIRED_INFO:
    if variant=='missing_info' and chrom=='chr7' and field=='AF_joint':continue
    kind='Float' if field.startswith('AF_') else 'Integer'
    if variant=='wrong_type' and chrom=='chr7' and field=='AF_joint':kind='String'
    text+=f'##INFO=<ID={field},Number=1,Type={kind},Description="fixture">\n'
   text+='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
   if variant!='empty_contigs':text+=f'{chrom}\t1\t.\tA\tC\t.\tPASS\t.\n'
   plain=self.file('gnomad/plain.vcf',text.encode());path.parent.mkdir(exist_ok=True);pysam.tabix_compress(str(plain),str(path),force=True);pysam.tabix_index(str(path),preset='vcf',force=True,csi=variant=='csi');plain.unlink()
   for p in [path,Path(str(path)+('.csi' if variant=='csi' else '.tbi'))]:os.utime(p,ns=(1700000000000000000,1700000000000000000));self.sources.append(p)
  if variant in ['missing_member','missing_index','corrupt_index']:
   p=family['chr7'] if variant=='missing_member' else Path(str(family['chr7'])+'.tbi')
   if variant=='corrupt_index':p.write_bytes(b'broken');os.utime(p,ns=(1700000000000000000,1700000000000000000))
   else:p.unlink();self.sources.remove(p)
  self.config_data.update(assembly=assembly,gnomad_vcf_chrX=str(family['chrX']))
  self.options.update(path=str(family['chrX']),assembly=assembly)

class Harness:
 def __init__(self,args):
  self.args=args;self.root=args.work_root;self.root.mkdir();self.rows=[];self.n=0
  c=self.root/'rename_fixture.c';c.write_text(SHIM);self.shim=self.root/'rename_fixture.so';subprocess.run(['cc','-shared','-fPIC','-O0',str(c),'-o',str(self.shim),'-ldl'],check=True)
 def invoke(self,f,arm):
  attempt=f.attempt;f.attempt+=1
  req=f.root/f'request_{attempt}.json';outcome_path=f.root/f'outcome_{attempt}.json';write(req,dict(helper=str(f.source),command=f.command,options=f.options));before=f.source_state()
  if arm=='reference':
   argv=[sys.executable,__file__,'reference','--request',str(req),'--output',str(outcome_path)]
  else:
   argv=[str(self.args.binary),'installer-config-paths',f.command]
   for key,value in f.options.items():
    argv.append('--'+key)
    if key not in ['check','no-source-root']:argv.append(value)
  lock_check=None
  if getattr(f,'lock_case',False):
   marker=f.root/'lock-entered';f.env['FIXTURE_FLOCK_ENTERED']=str(marker)
   f.cache.mkdir(exist_ok=True);original=f.config.read_bytes()
   with (f.cache/'.installer-adoption.lock').open('a') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX)
    process=subprocess.Popen(argv,env=f.env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,cwd=f.root)
    deadline=time.monotonic()+10
    while not marker.exists() and process.poll() is None and time.monotonic()<deadline:time.sleep(.01)
    time.sleep(.05)
    lock_check=dict(reached_lock=marker.exists(),blocked=process.poll() is None,config_preserved=f.config.read_bytes()==original,no_copy=sorted(p.name for p in f.cache.iterdir())==['.installer-adoption.lock'])
    fcntl.flock(lock,fcntl.LOCK_UN)
   stdout,stderr=process.communicate(timeout=30);done=subprocess.CompletedProcess(argv,process.returncode,stdout,stderr)
  else:done=subprocess.run(argv,env=f.env,capture_output=True,text=True,cwd=f.root)
  (f.root/f'stdout_{attempt}.log').write_text(done.stdout);(f.root/f'stderr_{attempt}.log').write_text(done.stderr)
  if arm=='reference':
   if done.returncode:raise RuntimeError('reference observer failed: '+done.stderr)
   out=json.loads(outcome_path.read_text());out['returncode']=0 if out['success'] else 1
  else:
   out=dict(success=done.returncode==0,returncode=done.returncode,stdout=done.stdout,value=None)
   if done.returncode:
    category,_,message=done.stderr.strip().rsplit('priva: ',1)[-1].partition(': ');out.update(category=category,message=message)
   elif f.command not in ['resolve','managed-roots','adopt-sources','render']:
    out['value']=done.stdout.strip() if f.command=='sha256' else json.loads(done.stdout)
  # Library adapters return JSON; original functions return objects directly.
  if f.command not in ['resolve','managed-roots','adopt-sources','render']:out['stdout']=''
  compared={k:out.get(k) for k in ['value','stdout','category','message']}
  # Python's real stderr uses backslashreplace for unencodable code points.
  # Preserve the raw exception above, but compare its actual printable channel.
  if isinstance(compared['message'],str):compared['message']=compared['message'].encode('utf-8',errors='backslashreplace').decode()
  out['products']=f.products();out['normalized']=f.normalize(compared);out['source_unchanged']=before==f.source_state();out['source_guard_passed']=out['source_unchanged']==f.expect_source_unchanged;out['argv']=argv;out['lock_check']=lock_check
  target=f.env.get('FIXTURE_RENAME_FAIL');out['failure_target_verified']=target is None or any(event.split('\t')[0]=='-1' and Path(event.split('\t')[-1]).name==target for event in out['products']['events'])
  return out
 def case(self,name,setup=lambda f:None,success=None,message=True,seed=0,recover=False):
  self.n+=1;arms=[];fixtures=[]
  for arm in ['reference','native']:
   f=Fixture(self.root/(f'{self.n:04}_{name}_{arm:_<9}'),self.args.repo,self.shim,seed);setup(f);f.save();fixtures.append(f);arms.append(self.invoke(f,arm))
  a,b=arms;passed=a['success']==b['success'] and a['returncode']==b['returncode'] and a['source_guard_passed'] and b['source_guard_passed'] and a['failure_target_verified'] and b['failure_target_verified']
  if success is not None:passed &= a['success']==success and b['success']==success
  if a['success']:passed &= a['normalized']['value']==b['normalized']['value'] and a['normalized']['stdout']==b['normalized']['stdout']
  else:
   passed &= a.get('category')==b.get('category')
   if message:passed &= a['normalized']['message']==b['normalized']['message']
  passed &= a['products']==b['products'] and not a['products']['receipt_checks'] and not b['products']['receipt_checks']
  if a['lock_check'] is not None or b['lock_check'] is not None:passed &= a['lock_check']==b['lock_check'] and all(a['lock_check'].values())
  recovery=[]
  if recover:
   for f,arm in zip(fixtures,['reference','native']):
    for key in ['FIXTURE_RENAME_FAIL','FIXTURE_FAIL_BEFORE','FIXTURE_CHANGE_AFTER_RENAME']:f.env.pop(key,None)
    # The successful publication must reuse already copied data without inode
    # replacement. Receipts/configs have their own intentional publication.
    existing={p:(p.stat().st_ino,p.stat().st_mtime_ns,p.stat().st_ctime_ns,digest(p)) for p in f.cache.rglob('*') if p.is_file() and not p.is_symlink() and not p.name.startswith('.') and p.name!='installer-source-adoptions.json'} if f.cache.exists() else {}
    recovered=self.invoke(f,arm)
    recovered['prior_payloads_preserved']=all(p.exists() and (p.stat().st_ino,p.stat().st_mtime_ns,p.stat().st_ctime_ns,digest(p))==identity for p,identity in existing.items())
    recovery.append(recovered)
   x,y=recovery
   passed &= all(v['success'] and v['source_guard_passed'] and v['prior_payloads_preserved'] and not v['products']['receipt_checks'] for v in recovery) and x['normalized']==y['normalized'] and x['products']==y['products']
  row=dict(name=name,passed=bool(passed),expected_success=success,exact_message_checked=message and not a['success'],reference=a,native=b,recovery=recovery);self.rows.append(row)
  write(self.root/'progress.json',dict(cases=len(self.rows),passed=sum(r['passed'] for r in self.rows),rows=self.rows))
  print(json.dumps(dict(case=name,passed=bool(passed),reference=a['success'],native=b['success'])),flush=True)
 def run(self):
  def mode(command,check=False,update=None):
   def setup(f):
    f.command=command
    if check:f.options['check']=''
    if update:f.config_data.update(update)
   return setup
  for seed in [*range(24),101,997,65535,4294967295]:self.case('hash_order_'+str(seed),mode('hash-order'),True,seed=seed)
  for assembly in ['hg19','GRCh37','hg38','GRCh38']:
   for command in ['portable-updates','resolve-cache-paths','resolve','managed-roots']:
    self.case(assembly+'_'+command,mode(command,command in ['resolve','managed-roots'],dict(assembly=assembly)),True)
  self.case('resolve_write',mode('resolve'),True)
  self.case('managed_write',mode('managed-roots'),True)
  self.case('output_preflight',lambda f:f.options.update(output=str(f.root/'output.yaml')),True)
  invalid={
   'empty_doc':'','root_sequence':'[]','root_string':'text','root_integer':'7',
   'duplicate':'a: 1\na: 2\n','nested_duplicate':'ignore: {a: 1, a: 2}',
   'bool_int_duplicate':'true: x\n1: y','float_int_duplicate':'1.0: x\n1: y','nan_duplicate':'.nan: x\n.nan: y',
   'unhashable':'? [x]\n: y','merge':'defaults: &d {x: 1}\nother: {<<: *d}',
   'value_tag':'=: value','recursive_map':'a: &a {self: *a}','duplicate_before_bad_value':'a: 1\na: !bad x',
   'bad_tag':'a: !bad value','bad_yaml':'a: [1,','second_doc':'a: 1\n---\nb: 2',
  }
  for name,text in invalid.items():self.case(name,lambda f,text=text:setattr(f,'raw',text),False)
  for key,value in [('cache_root','/'),('cache_root','relative'),('cache_root','/with space'),('cache_root',None),('tmp_dir','/'),('assembly','GRCh39'),('assembly',[]),('vep_installed_version','1.0'),('vep_installed_version',True),('cache_folders',[]),('cache_folders',False),('cache_folders',dict(wrong='here')),('cache_folders',dict(gnomad='../x')),('cache_folders',dict(clinvar=1)),('vep_cache','flat'),('repeat_region_file_name','../bad'),('repeat_region_file_name','.'),('repeat_region_file_name',False),('hpo_release',17),('nm_to_enst_tsv',1),('hgnc_table',[]),('alphamissense_vcf',True)]:
   self.case('invalid_'+key+'_'+str(len(self.rows)),mode('resolve',True,{key:value}),False)
  for key,value in [('pext_enabled',True),('pext_enabled','TRUE'),('hpo_release','v2020'),('mondo_release',17),('cache_root','//cache-test'),('cache_folders',dict(vep_cache='a',clinvar='b'))]:self.case('valid_'+key+'_'+str(len(self.rows)),mode('portable-updates',False,{key:value}),True)
  def yaml_extra(f,extra):
   import yaml
   f.raw=yaml.safe_dump(f.config_data,sort_keys=False)+extra
  for name,extra in [('anchors','ignored: &x {a: true, b: null}\ncopy: *x\n'),('flow','ignored: {x: [a, b], y: 1}\n'),('recursive_sequence','ignored: &x [*x]\n'),('date','ignored: 2024-01-02\n'),('set','ignored: !!set {a: null, b: null}\n'),('binary','ignored: !!binary SGVsbG8=\n'),('wide_int','ignored: '+str(10**100)+'\n')]:self.case('yaml_'+name,lambda f,e=extra:yaml_extra(f,e),True)
  self.case('different_rounded_numeric_keys',lambda f:yaml_extra(f,'18014398509481985: a\n18014398509481984.0: b\n'),True)
  self.case('wide_hex_version',mode('resolve',True,{'vep_installed_version':2**160-1}),True)
  def raw_hex(f):
   import yaml
   f.raw=yaml.safe_dump(f.config_data,sort_keys=False).replace("vep_installed_version: '116'","vep_installed_version: 0x"+'f'*40)
   f.options['check']=''
  self.case('raw_hex_version',raw_hex,True)
  def json_config(f,value,check=False):
   f.raw=json.dumps(dict(f.config_data,ignored=value),ensure_ascii=True)
   if check:f.options['check']=''
  for name,value in [('ordinary',dict(unicode='é',truth=True,missing=None,small=1e-7,zero=-0.0)),('wide',10**100),('deep',[[None]]),('controls','\u0000\n\t'),('unicode','é𝄞')]:self.case('json_'+name,lambda f,v=value:json_config(f,v),True)
  self.case('json_surrogate_check',lambda f:json_config(f,'\ud800',True),True)
  self.case('json_surrogate_write',lambda f:json_config(f,'\ud800'),False)
  self.case('json_surrogates_write',lambda f:json_config(f,'\ud800\ud801'),False)
  def surrogate_adopt(f):f.selected('primateai_prescore');json_config(f,'\ud800');f.command='managed-roots'
  self.case('json_surrogate_managed_adoption',surrogate_adopt,False)
  def surrogate_missing_output(f):json_config(f,'\ud800');f.options['output']=str(f.root/'absent/output.yaml')
  self.case('json_surrogate_missing_output',surrogate_missing_output,False)
  for key in ['primateai_prescore','ref_genome','gnomad_vcf_chrX','nm_to_enst_tsv','hgnc_table','clinvar_gene_stat','mondo_obo']:
   self.case('missing_'+key,lambda f,k=key:f.config_data.update({k:str(f.root/'missing')}),False)
  for key in ['clinvar_vcf','alphamissense_pd_stat','gofcards_mechanism_review_tsv']:
   self.case('missing_allowed_'+key,lambda f,k=key:f.config_data.update({k:str(f.root/'missing')}),True)
  self.case('disabled_pext',lambda f:f.config_data.update(pext_enabled=False,pext_mean_bw=str(f.root/'missing')),True)
  self.case('grch37_mave_dormant',lambda f:f.config_data.update(assembly='hg19',mavedb_file=str(f.root/'missing')),True)
  bundles={
   'clinvar_gene_stat':'data/ClinVar/clinvar_2star_stats.pkl.gz','interpro_mapping_pickle':'data/InterPro/Interpro_entry_mapping.pkl.gz',
   'hgnc_table':'data/hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv','dispensable_gene_list':'data/dispensable_genes/dispensable_gene_list.txt',
   'gene_nonlof_mechanism_schema':'data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.schema.json','mondo_disease_scope_overrides':'data/mondo/disease_scope_overrides.tsv',
   'clingen_map':'data/clingen/clingen_map.hg38.pkl.gz','gene_dosage_sensitivity':'data/clingen/gene_dosage_sensitivity.hg38.tsv','loeuf_prescore':'data/loeuf/loeuf_dataset_hg38.tsv.gz',
   'gofcards_exact_gof_cache':'data/gofcards/gofcards_exact_gof.json.gz','gene_nonlof_mechanism_json':'data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.json.gz',
  }
  def bundle(f,key,kind):
   relative=bundles[key];p=f.file('source/'+relative,b'opaque bundle, no format checks\xff\n');f.sources.append(p)
   if kind=='relative':f.config_data[key]=relative
   elif kind=='absolute':f.config_data[key]=str(p)
   elif kind=='legacy':f.config_data[key]=p.name
   elif kind=='empty':p.write_bytes(b'');f.config_data[key]=relative
   elif kind=='external':
    p.unlink();q=f.file('outside/bundle_source',b'opaque external same basename\n');p.symlink_to(q);f.config_data[key]=str(p);f.sources.remove(p);f.sources.append(q)
  for key in bundles:
   self.case('bundle_'+key,lambda f,k=key:bundle(f,k,'relative'),True)
  for key in ['clinvar_gene_stat','interpro_mapping_pickle']:
   for kind in ['absolute','legacy','empty','external']:self.case('bundle_'+key+'_'+kind,lambda f,k=key,t=kind:bundle(f,k,t),True)
  def hgnc_alias(f):
   bundle(f,'hgnc_table','absolute');alias=f.source/'data/hgnc/non_alt_loci_set.tsv';alias.symlink_to(Path(bundles['hgnc_table']).name);f.config_data['hgnc_table']='data/hgnc/non_alt_loci_set.tsv'
  self.case('hgnc_approved_alias',hgnc_alias,True)
  def forced(f):
   f.file('source/'+bundles['gofcards_exact_gof_cache'],b'opaque bundle');f.env['PRIVA_FORCE_GOFCARDS_CACHE']='1';f.command='portable-updates'
  self.case('force_unset_bundle',forced,True)
  def forced_missing(f):f.config_data['gofcards_mechanism_review_tsv']=str(f.root/'missing');f.env['PRIVA_FORCE_GOFCARDS_CACHE']='1'
  self.case('force_missing_review',forced_missing,False)
  def source_file(f):f.selected('primateai_prescore')
  self.case('file_adoption',source_file,True)
  def idempotent(f):source_file(f);f.complete_bundles()
  self.case('idempotent_second_resolution',idempotent,True,recover=True)
  def locked(f):source_file(f);f.lock_case=True
  self.case('blocking_advisory_lock',locked,True)
  def companions(f):
   p=f.selected('alphamissense_prescore','outside/alpha.tsv.gz')
   for suffix in ['.tbi','.csi','.fai','.gzi','.provenance.json','.nm_to_enst.json','.ready.json','.source.json','.annotation-base.json']:
    q=f.file(str(p.relative_to(f.root))+suffix,b'opaque companion\n');f.sources.append(q)
  self.case('all_companions',companions,True)
  def existing(f,same=False,link=False):
   p=f.selected('primateai_prescore');target=f.cache/'VEP_plugin_cache/PrimateAI/supplied/primateai_prescore'/p.name
   if link:target.parent.mkdir(parents=True);q=f.file('outside/keep.dat',b'outside protected\n');f.sources.append(q);target.symlink_to(q)
   else:f.file(str(target.relative_to(f.root)),p.read_bytes() if same else b'old bytes\n',mode=0o600)
  self.case('retain_previous',lambda f:existing(f),True)
  self.case('identical_destination',lambda f:existing(f,True),True)
  self.case('retain_symlink',lambda f:existing(f,link=True),True)
  def directory(f):
   a=f.file('old_plugins/child/a.txt');b=f.file('old_plugins/child/b.txt');f.sources.extend([a,b]);f.config_data['vep_plugins_dir']=str(f.root/'old_plugins')
  self.case('directory_adoption',directory,True)
  def child(f):
   directory(f);f.config_data['loftee_parent_dir']=str(f.root/'old_plugins/child')
  self.case('directory_child_dedup',child,True)
  def dangling(f):
   p=f.root/'broken';p.symlink_to(f.root/'absent');f.config_data['primateai_prescore']=str(p)
  self.case('broken_source',dangling,False)
  def wrong_type(f):
   p=f.root/'directory';p.mkdir();f.config_data['primateai_prescore']=str(p)
  self.case('wrong_file_type',wrong_type,False)
  def escape(f):
   source_file(f);outside=f.root/'external';outside.mkdir();(f.cache/'VEP_plugin_cache').mkdir(parents=True);(f.cache/'VEP_plugin_cache/PrimateAI').symlink_to(outside,target_is_directory=True)
  self.case('managed_directory_escape',escape,False)
  def source_link(f):
   p=f.selected('primateai_prescore');link=f.root/'selected-link';link.symlink_to(p);f.config_data['primateai_prescore']=str(link)
  self.case('source_file_link',source_link,True)
  def existing_directory_link(f):
   directory(f);f.cache.mkdir();outside=f.root/'outside-directory';outside.mkdir();(f.cache/'VEP_plugin').symlink_to(outside,target_is_directory=True)
  self.case('existing_directory_link',existing_directory_link,False)
  def nested_directory_link(f):
   directory(f);dest=f.cache/'VEP_plugin';dest.mkdir(parents=True);outside=f.root/'outside-directory';outside.mkdir();(dest/'child').symlink_to(outside,target_is_directory=True)
  self.case('nested_directory_link',nested_directory_link,False)
  def xattrs(f):
   p=f.selected('primateai_prescore');os.setxattr(p,b'user.priva-fixture',b'preserve-xattr')
  self.case('copy_xattrs',xattrs,True)
  for name,ending in [('non_bmp','𝄞'),('private_transport','\U000f0000'),('private_escape','\U000f0800'),('literal_escape','\\ud800'),('byte_surrogate','\udc80')]:
   self.case('unicode_source_'+name,lambda f,s=ending:f.selected('primateai_prescore','outside/selected_'+s+'.dat'),True)
  def byte_directory(f):
   p=f.file('old_plugins/child/byte_\udcff.txt');f.sources.append(p);f.config_data['vep_plugins_dir']=str(f.root/'old_plugins')
  self.case('byte_directory_entry',byte_directory,True)
  def transport_metadata(f):
   import yaml
   f.raw=yaml.safe_dump(f.config_data,sort_keys=False)+'ignored: "\\U000f0000\\ud800\\U000f0800"\n'
  self.case('literal_private_and_surrogate_metadata',transport_metadata,True)
  def receipt(f,value):
   source_file(f);f.cache.mkdir();(f.cache/'installer-source-adoptions.json').write_text(value)
  for name,value in [('malformed','{'),('missing_copies','{}'),('wrong_copies','{"copies":1}'),('array','[]'),('null','null'),('string','"text"'),('integer','1'),('float','1.0'),('empty','{"copies":[]}'),('wide','{"copies":[],"extra":'+str(10**100)+'}'),('surrogate','{"copies":[],"extra":"\\ud800"}'),('nested','{"copies":[],"extra":'+'['*600+'null'+']'*600+'}'),('duplicate','{"extra":1,"copies":[],"extra":2}'),('mixed','{"copies":[],"extra":{"null":null,"true":true,"false":false,"float":-0.0,"unicode":"é"}}'),('float_forms','{"copies":[],"extra":[1e-7,1e30,1.25,1.0,1e-5,1e-6,NaN,Infinity,-Infinity]}')]:
   self.case('receipt_'+name,lambda f,v=value:receipt(f,v),name not in ['malformed','missing_copies','wrong_copies','array','null','string','integer','float'],message=name!='malformed')
  def fail(f,target,before=False):
   source_file(f);f.complete_bundles();f.env['FIXTURE_RENAME_FAIL']=target
   if before:f.env['FIXTURE_FAIL_BEFORE']='1'
  for target in ['selected.dat','installer-source-adoptions.json','config.yaml']:
   for before in [False,True]:self.case('rename_'+target+'_'+str(before),lambda f,t=target,b=before:fail(f,t,b),False,recover=True)
  def change_config(f):source_file(f);f.env['FIXTURE_CHANGE_AFTER_RENAME']=str(f.config)
  self.case('concurrent_config',change_config,False)
  def change_source(f):
   p=f.selected('primateai_prescore');f.env['FIXTURE_CHANGE_AFTER_RENAME']=str(p);f.expect_source_unchanged=False
  self.case('concurrent_source_after_copy',change_source,False)
  def prepared(f,where='external',variant=None,check=False):
   f.prepared_pair(where,variant)
   if check:f.options['check']=''
  for where in ['external','cache/selected']:
   self.case('prepared_'+where.replace('/','_'),lambda f,w=where:prepared(f,w),True)
  def root_alias(f,selected_alias):
   prepared(f,'cache/selected');alias=f.root/'cache_alias';alias.symlink_to(f.cache,target_is_directory=True);f.config_data['cache_root']=str(alias)
   if selected_alias:
    for key in ['spliceai_snv_prescore','spliceai_indel_prescore']:f.config_data[key]=str(alias/Path(f.config_data[key]).relative_to(f.cache))
  for aliased in [False,True]:self.case('prepared_root_alias_'+str(aliased),lambda f,a=aliased:root_alias(f,a),True)
  def prepared_collision(f,external=True):
   f.prepared_pair('cache/VEP_plugin_cache/SpliceAI2/GRCh38');f.prepared_pair('external' if external else 'cache/selected')
  self.case('prepared_collision',prepared_collision,False)
  self.case('prepared_managed_other_canonical',lambda f:prepared_collision(f,False),True)
  def prepared_external_link(f):
   prepared(f);f.cache.mkdir();link=f.cache/'outside_link';link.symlink_to(f.prepared,target_is_directory=True)
   for key in ['spliceai_snv_prescore','spliceai_indel_prescore']:f.config_data[key]=str(link/Path(f.config_data[key]).name)
  self.case('prepared_external_link',prepared_external_link,True)
  def native_splice(f):
   m=load(f.source/'scripts/spliceai_cache_input.py');header=('\t'.join(m.native_columns())+'\n').encode()
   for key,name in [('spliceai_snv_prescore','renamed.bin'),('spliceai_indel_prescore','renamed.vcf.gz')]:f.selected(key,'outside/'+name,header)
  self.case('native_spliceai_sources_preserved',native_splice,True)
  self.case('prepared_missing',lambda f:prepared(f,'external','missing'),False)
  self.case('prepared_incomplete',lambda f:prepared(f,'external','incomplete'),False)
  self.case('prepared_incomplete_plan',lambda f:prepared(f,'external','incomplete',True),True)
  def prepared_wrong_assembly(f):prepared(f,'cache/selected');f.config_data['assembly']='hg19'
  self.case('prepared_wrong_assembly',prepared_wrong_assembly,False)
  def mixed_splice(f):
   prepared(f,'cache/selected');path=Path(f.config_data['spliceai_indel_prescore']);path.write_bytes(b'##fileformat=VCFv4.2\n##reference=GRCh38\n##SpliceAIVersion=1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
  self.case('mixed_splice',mixed_splice,False)
  for variant in [None,'csi','empty_contigs','missing_member','missing_index','corrupt_index','missing_info','wrong_type','wrong_length']:
   def setup(f,v=variant):f.gnomad(variant=v);f.command='gnomad-validate'
   self.case('gnomad_'+str(variant),setup,variant in [None,'csi','empty_contigs'])
  self.case('gnomad_adopt',lambda f:f.gnomad(),True)
  def gnomad37(f):f.gnomad('hg19');f.command='gnomad-validate'
  self.case('gnomad37',gnomad37,True)
  for key in ['repeat_region_file_name','hpo_release','primateai_prescore','nm_to_enst_tsv','clinvar_vcf','vep_installed_version','assembly']:
   for check in [True,False]:
    def escaped_field(f,k=key,c=check):
     import yaml
     f.config_data.pop(k,None)
     f.raw=yaml.safe_dump(f.config_data,sort_keys=False)+k+': "\\ud800"\n'
     if c:f.options['check']=''
    self.case('escaped_consumed_'+key+'_'+str(check),escaped_field)
  result=dict(cases=len(self.rows),passed=sum(r['passed'] for r in self.rows),successes=sum(r['reference']['success'] for r in self.rows),rows=self.rows)
  write(self.root/'complete.json',result);return 0 if result['cases']==result['passed'] else 1

def main():
 parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
 child=commands.add_parser('check')
 for key in ['repo','binary','work-root']:child.add_argument('--'+key,type=Path,required=True)
 child=commands.add_parser('reference');child.add_argument('--request',type=Path,required=True);child.add_argument('--output',type=Path,required=True)
 args=parser.parse_args();return Harness(args).run() if args.command=='check' else reference(args)
if __name__=='__main__':raise SystemExit(main())
