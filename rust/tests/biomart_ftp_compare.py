#!/usr/bin/env python3
"""Exact bounded original/native BioMart FTP projection comparisons.

check --producer SCRIPT --binary PRIVA --work-root NEW --report JSON
oracle --producer SCRIPT --request JSON --result JSON
No network, installed writes, table filtering or scientific source substitution.
Each fixture has separate input/output files. Reuse accepts only identical
source/case definitions and preserved prior oracle bytes/modes/link structure.
"""
import argparse,copy,gzip,hashlib,importlib.util,json,os,re,shutil,sqlite3,stat,subprocess,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
DATASET='hsapiens_gene_ensembl'
PROJECTIONS={
 'transcripts':(DATASET+'__transcript__main',['gene_id_1020_key','transcript_id_1064_key','stable_id_1023','stable_id_1066','biotype_1064']),
 'main':(DATASET+'__translation__main',['gene_id_1020_key','transcript_id_1064_key','translation_id_1068_key','stable_id_1023','stable_id_1066','biotype_1064']),
 'hgnc':(DATASET+'__ox_hgnc__dm',['gene_id_1020_key','display_label_1074']),
 'swissprot':(DATASET+'__ox_uniprotswissprot__dm',['translation_id_1068_key','dbprimary_acc_1074']),
 'refseq':(DATASET+'__ox_refseq_mrna__dm',['transcript_id_1064_key','dbprimary_acc_1074'])}
ROWS={
 'transcripts':[[1,1,'ENSG1','ENST1','protein_coding'],[2,2,'ENSG2','ENST2','lncRNA'],[3,3,'ENSG3','ENST3','lncRNA']],
 'main':[[1,1,11,'ENSG1','ENST1','protein_coding'],[2,2,None,'ENSG2','ENST2','lncRNA'],[3,3,0,'ENSG3','ENST3','lncRNA']],
 'hgnc':[[1,'GENE1'],[1,'GENE1'],[2,''],[None,'NULL'],[0,'ZERO']],
 'swissprot':[[11,'P1'],[11,'P2'],[None,'NULL'],[0,'']],
 'refseq':[[1,'NM_1'],[1,'NM_2'],[2,None],[None,'NULL']]}
def pin(path):
 path=Path(path)
 with path.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
 return dict(path=str(path),bytes=path.stat().st_size,sha256=digest)
def save(path,value):
 with Path(path).open('x') as f:json.dump(value,f,indent=2,sort_keys=True);f.write('\n')
def module(path):
 spec=importlib.util.spec_from_file_location('ftp_source',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def oracle(a):
 m=module(a.producer);r=json.loads(a.request.read_text());o=r['options'];cmd=r['command'];p=lambda k:Path(o[k]);db=None
 try:
  if cmd=='ftp-schema':value=m.read_schema(p('schema'))
  elif cmd=='ftp-project':
   schema=m.read_schema(p('schema'));db=sqlite3.connect(p('db'));db.execute('PRAGMA cache_size=-131072');imports={}
   for alias,(table,fields) in m.PROJECTIONS.items():imports[alias]=m.load_projection(db,p('raw'),alias,table,fields,schema[table])
   audit=m.relationship_audit(db);export=m.export_candidate(db,p('output'),audit)
   m.require(db.execute('PRAGMA integrity_check').fetchone()==('ok',),'FTP projection SQLite integrity failed')
   value=dict(imports=imports,audit=audit,export=export)
  else:
   db=sqlite3.connect(p('db'))
   if cmd=='ftp-load':value=m.load_projection(db,p('raw'),o['alias'],o['table'],json.loads(o['columns']),json.loads(p('schema').read_text()))
   elif cmd=='ftp-audit':value=m.relationship_audit(db)
   elif cmd=='ftp-export':value=m.export_candidate(db,p('output'),json.loads(p('audit').read_text()))
   else:raise AssertionError(cmd)
 except Exception as e:
  save(a.result,dict(success=False,category=type(e).__name__,message=str(e)));raise SystemExit(1)
 finally:
  if db is not None:db.close()
 save(a.result,dict(success=True,value=value));print(json.dumps(value))
def schemas():
 result={table:[(name,'int(10),' if name.endswith('_key') else 'varchar(100),') for name in columns] for table,columns in PROJECTIONS.values()}
 result['meta_conf__dataset__main']=[(v,'varchar(10),') for v in ['dataset_id_key','dataset','display_name','description','type','visible','version','modified']]
 result['meta_conf__xml__dm']=[(v,'varchar(10),') for v in ['dataset_id_key','xml','compressed_xml','message_digest']]
 return result
def schema_text(schema):return ''.join('CREATE TABLE `'+table+'` (\n'+''.join('  `'+name+'` '+kind+'\n' for name,kind in columns)+') ENGINE=MyISAM CHARSET=latin1;\n' for table,columns in schema.items())
def encoded(value):
 if value is None:return b'\\N'
 if isinstance(value,bytes):return value
 value=str(value).encode('ascii');out=bytearray()
 escapes={0:b'\\0',8:b'\\b',10:b'\\n',13:b'\\r',9:b'\\t',26:b'\\Z',92:b'\\\\'}
 for b in value:out.extend(escapes.get(b,bytes([b])))
 return bytes(out)
def dump(rows):return b''.join(b'\t'.join(map(encoded,row))+b'\n' for row in rows)
def cases():
 result={}
 def add(name,command,**spec):result[name]=dict(command=command,**spec)
 base=schema_text(schemas())
 for kind,payload in {'ordinary':base,'crlf':base.replace('\n','\r\n'),'cr':base.replace('\n','\r'),'no_final_newline':base.rstrip('\n'),'latin1_comment':'-- café\n'+base,'unselected_bad_charset':'CREATE TABLE `ignored` (\n) CHARSET=utf8;\n'+base,'unselected_bad_column':'CREATE TABLE `ignored` (\n  `bad\n'+base,'missing_close':base.rsplit(') ENGINE',1)[0],'create_unquoted':'CREATE TABLE bad (\n'+base,'create_empty':'CREATE TABLE `` (\n'+base,'column_empty':base.replace('`gene_id_1020_key`','``',1),'column_unclosed':base.replace('`gene_id_1020_key`','`gene_id_1020_key',1),'column_one_space':base.replace('int(10),','',1),'column_two_spaces':base.replace('int(10),',' ',1),'column_nbsp':base.replace('` int(10),','`\xa0int(10),'),'charset':base.replace('CHARSET=latin1','CHARSET=utf8',1),'key_type':base.replace('int(10),','bigint(20),',1),'key_case':base.replace('int(10),','INT(10),',1),'leading_other':'-- anything\x00\xff\n'+base,'boundary_cr':'--'+'x'*8189+'\r\n'+base}.items():add('schema_'+kind,'ftp-schema',schema_hex=payload.encode('latin1').hex())
 for name in schemas():
  missing=schemas();del missing[name];add('schema_missing_'+name,'ftp-schema',schema_hex=schema_text(missing).encode().hex())
  duplicate=base+schema_text({name:schemas()[name]});add('schema_duplicate_'+name,'ftp-schema',schema_hex=duplicate.encode().hex())
 for kind in ['missing','directory','notgzip','truncated','badcrc','multimember']:
  add('schema_input_'+kind,'ftp-schema',schema_input=kind)
 columns=['k_key','stable_id_1023','stable_id_1066','label'];row=[1,'ENSG1','ENST1','A']
 for kind,rows in {'ordinary':[row,[2,'ENSG2','ENST2','']],'nulls':[[None,None,None,None],row],'zero':[[0,'ENSG0','ENST0','']],'max_i64':[[2**63-1,'ENSG1','ENST1','A']],'overflow':[[2**63,'ENSG1','ENST1','A']],'leading_zeros':[[b'0001','ENSG1','ENST1','A']],'zero_4300':[[b'0'*4300,'ENSG1','ENST1','A']],'zero_4301':[[b'0'*4301,'ENSG1','ENST1','A']],'stable_short':[[1,'ENSG','ENST1','A']],'stable_unicode':[[1,b'ENSG\xff','ENST1','A']],'stable_suffix':[[1,'ENSG1.2','ENST1','A']],'negative':[[-1,'ENSG1','ENST1','A']],'plus':[[b'+1','ENSG1','ENST1','A']],'key_space':[[b' 1','ENSG1','ENST1','A']],'key_decimal':[[b'1.0','ENSG1','ENST1','A']],'key_empty':[[b'','ENSG1','ENST1','A']],'key_nonascii':[[b'\xff','ENSG1','ENST1','A']],'short':[[1,'ENSG1']],'long':[row+['extra']],'empty':[],'null_then_blank':[[1,None,None,''],[None,'ENSG1','ENST1',None]],'overflow_then_bad':[[2**63,'ENSG1','ENST1','A'],[1,'bad','ENST1','A']]}.items():add('load_'+kind,'ftp-load',rows_hex=dump(rows).hex())
 for payload in [b'bad\\',b'\\N',b'\\\\N',b'\\q',b'\\0',b'\\b',b'\\n',b'\\r',b'\\t',b'\\Z',b'\\\\',b'\xff',b'\x00',b'\r',b'"quoted"']:
  add('load_escape_'+payload.hex(),'ftp-load',rows_hex=dump([[1,'ENSG1','ENST1',payload]]).hex())
 for count in [9999,10000,10001]:
  add('load_batch_'+str(count),'ftp-load',repeat=count)
  add('load_batch_bad_'+str(count),'ftp-load',repeat=count,tail_hex=dump([[1,'bad','ENST1','A']]).hex())
 for count in [9998,9999,10000]:add('load_deferred_overflow_'+str(count),'ftp-load',repeat=count,tail_hex=dump([[2**63,'ENSG1','ENST1','A'],[1,'bad','ENST1','A']]).hex())
 for kind in ['missing','directory','notgzip','truncated','badcrc','multimember','no_newline','crlf','trailing_blank']:
  add('load_input_'+kind,'ftp-load',raw_input=kind)
 for kind in ['missing_column','duplicate_column','wrong_alias','existing_table','existing_index','raw_is_db','db_directory','db_parent_missing']:
  add('load_setup_'+kind,'ftp-load',setup=kind)
 for kind in ['ordinary','duplicate_transcript','no_transcripts','null_gene','zero_gene','null_stable','missing_main','changed_main','hgnc_orphan','swissprot_orphan','refseq_orphan','no_untranslated','no_noncoding','untranslated_swissprot','biotype_null','biotype_null_both','missing_table','dimension_null','zero_orphan']:
  add('audit_'+kind,'ftp-audit',rows_kind=kind)
 for kind in ['ordinary','empty','missing_table','unicode','unicode_run','blob','number','float','first_unicode','nulls','all_blank','wrong_count','missing_count','audit_null','audit_list']:
  add('export_'+kind,'ftp-export',rows_kind=kind)
 for kind in ['file','directory','link','broken','loop','part_file','part_directory','part_link','part_broken','part_loop','parent_missing']:
  add('export_destination_'+kind,'ftp-export',output_kind=kind)
 for kind in ['ordinary','main_bad_row','swissprot_missing','hgnc_empty','bad_audit','existing_output','existing_db']:
  add('project_'+kind,'ftp-project',project_kind=kind)
 return result

def mutate(rows,kind):
 rows=copy.deepcopy(rows)
 if kind=='duplicate_transcript':rows['transcripts'].append(rows['transcripts'][0])
 if kind=='no_transcripts':rows['transcripts']=[]
 if kind=='null_gene':rows['transcripts'][0][0]=None
 if kind=='zero_gene':rows['transcripts'][0][0]=0
 if kind=='null_stable':rows['transcripts'][0][2]=None
 if kind=='missing_main':rows['main'].pop()
 if kind=='changed_main':rows['main'][0][3]='ENSG9'
 if kind.endswith('_orphan') and kind.split('_')[0] in ['hgnc','swissprot','refseq']:rows[kind.split('_')[0]].append([99,'ORPHAN'])
 if kind=='no_untranslated':
  rows['main'][1][2]=22;rows['main'][2][2]=33
 if kind=='no_noncoding':
  for alias in ['transcripts','main']:
   for row in rows[alias]:row[-1]='protein_coding'
 if kind=='untranslated_swissprot':rows['swissprot'].append([0,'BAD'])
 if kind in ['biotype_null','biotype_null_both']:
  rows['transcripts'][0][-1]=None
  if kind.endswith('_both'):rows['main'][0][-1]=None
 if kind=='dimension_null':rows['hgnc']=[[None,None]];rows['swissprot']=[[None,None]];rows['refseq']=[[None,None]]
 if kind=='zero_orphan':rows['hgnc'].append([0,'ORPHAN'])
 return rows

def setup(folder,spec,producer):
 folder.mkdir();raw=folder/'raw';raw.mkdir();source=folder/'source/scripts';source.mkdir(parents=True)
 for name in ['biomart_ftp_nm_to_enst.py','biomart_nm_to_enst.py']:shutil.copyfile(producer.with_name(name),source/name)
 cmd=spec['command'];options={};schema=schemas();data=bytes.fromhex(spec.get('schema_hex',schema_text(schema).encode().hex()));schema_path=folder/'schema.sql.gz';schema_path.write_bytes(gzip.compress(data,mtime=0))
 kind=spec.get('schema_input')
 if kind=='missing':schema_path.unlink()
 if kind=='directory':schema_path.unlink();schema_path.mkdir()
 if kind=='notgzip':schema_path.write_bytes(b'not gzip')
 if kind=='truncated':schema_path.write_bytes(schema_path.read_bytes()[:-5])
 if kind=='badcrc':d=bytearray(schema_path.read_bytes());d[-8]^=1;schema_path.write_bytes(d)
 if kind=='multimember':schema_path.write_bytes(gzip.compress(data[:100],mtime=0)+gzip.compress(data[100:],mtime=0))
 db=folder/'projection.sqlite';output=folder/'output.tsv'
 if cmd=='ftp-schema':options['schema']=str(schema_path)
 elif cmd=='ftp-load':
  columns=['k_key','stable_id_1023','stable_id_1066','label'];fields=[(v,'ignored') for v in columns]
  payload=bytes.fromhex(spec.get('rows_hex',dump([[1,'ENSG1','ENST1','A']]).hex()))
  if 'repeat' in spec:payload=dump([[1,'ENSG1','ENST1','A']])*spec['repeat']+bytes.fromhex(spec.get('tail_hex',''))
  path=raw/'table.txt.gz';path.write_bytes(gzip.compress(payload,mtime=0));kind=spec.get('raw_input')
  if kind=='missing':path.unlink()
  if kind=='directory':path.unlink();path.mkdir()
  if kind=='notgzip':path.write_bytes(b'not gzip')
  if kind=='truncated':path.write_bytes(path.read_bytes()[:-5])
  if kind=='badcrc':d=bytearray(path.read_bytes());d[-8]^=1;path.write_bytes(d)
  if kind=='multimember':path.write_bytes(gzip.compress(payload[:3],mtime=0)+gzip.compress(payload[3:],mtime=0))
  if kind=='no_newline':path.write_bytes(gzip.compress(payload.rstrip(b'\n'),mtime=0))
  if kind=='crlf':path.write_bytes(gzip.compress(payload.replace(b'\n',b'\r\n'),mtime=0))
  if kind=='trailing_blank':path.write_bytes(gzip.compress(payload+b'\n',mtime=0))
  alias='x';kind=spec.get('setup')
  if kind=='missing_column':fields=fields[1:]
  if kind=='duplicate_column':columns.append(columns[0])
  if kind=='wrong_alias':alias='x"invalid'
  if kind in ['existing_table','existing_index']:
   with sqlite3.connect(db) as conn:
    if kind=='existing_table':conn.execute('CREATE TABLE x(z)')
    else:conn.execute('CREATE TABLE old(z)');conn.execute('CREATE INDEX x_k_key ON old(z)')
  if kind=='raw_is_db':db=path
  if kind=='db_directory':db.mkdir()
  if kind=='db_parent_missing':db=folder/'missing/db'
  schema_path=folder/'columns.json';schema_path.write_text(json.dumps(fields))
  options=dict(db=str(db),raw=str(raw),alias=alias,table='table',columns=json.dumps(columns),schema=str(schema_path))
 else:
  rows=mutate(ROWS,spec.get('rows_kind','ordinary'));kind=spec.get('rows_kind')
  if cmd=='ftp-project':
   kind=spec.get('project_kind');options=dict(db=str(db),raw=str(raw),schema=str(schema_path),output=str(output))
   if kind=='main_bad_row':rows['main'][1][3]='bad'
   if kind=='hgnc_empty':rows['hgnc']=[]
   if kind=='bad_audit':rows['main'].pop()
   if kind=='existing_output':output.write_bytes(b'existing output')
   if kind=='existing_db':
    with sqlite3.connect(db) as conn:conn.execute('CREATE TABLE transcripts(x)')
   for alias,(table,columns) in PROJECTIONS.items():
    if kind=='swissprot_missing' and alias=='swissprot':continue
    (raw/(table+'.txt.gz')).write_bytes(gzip.compress(dump(rows[alias]),mtime=0))
  else:
   if cmd=='ftp-export':
    if kind=='empty':rows['main']=[]
    if kind in ['unicode','unicode_run','first_unicode']:rows['hgnc'][0 if kind=='first_unicode' else 1][1]='éΩ' if kind=='unicode_run' else 'é'
    if kind=='blob':rows['hgnc'][1][1]=b'BLOB'
    if kind in ['number','float']:rows['main'][1][3]=2 if kind=='number' else 2.5
    if kind=='nulls':rows['main'][1][3]=None
    if kind=='all_blank':rows['main'][1][3]='';rows['main'][1][4]=''
   with sqlite3.connect(db) as conn:
    for alias,(table,columns) in PROJECTIONS.items():
     if kind=='missing_table' and alias=='swissprot':continue
     conn.execute('CREATE TABLE '+alias+' (__row INTEGER PRIMARY KEY, '+','.join('"'+n+'"' for n in columns)+')')
     conn.executemany('INSERT INTO '+alias+' VALUES ('+','.join('?' for _ in range(len(columns)+1))+')',[(i+1,*row) for i,row in enumerate(rows[alias])])
   if cmd=='ftp-audit':options=dict(db=str(db))
   else:
    with sqlite3.connect(db) as conn:
     try:expected=conn.execute('SELECT COUNT(*) FROM main m LEFT JOIN hgnc h ON h.gene_id_1020_key=m.gene_id_1020_key LEFT JOIN swissprot u ON u.translation_id_1068_key=m.translation_id_1068_key LEFT JOIN refseq r ON r.transcript_id_1064_key=m.transcript_id_1064_key').fetchone()[0]
     except sqlite3.Error:expected=0
    audit={'expected_output_rows':expected+(kind=='wrong_count')}
    if kind=='missing_count':audit={}
    if kind=='audit_null':audit=None
    if kind=='audit_list':audit=[]
    path=folder/'audit.json';path.write_text(json.dumps(audit));options=dict(db=str(db),output=str(output),audit=str(path))
    outkind=spec.get('output_kind','');target=folder/'target';target.write_bytes(b'old target')
    part=Path(str(output)+'.part')
    dest=part if outkind.startswith('part_') else output;outkind=outkind.removeprefix('part_')
    if outkind=='file':dest.write_bytes(b'existing')
    if outkind=='directory':dest.mkdir()
    if outkind=='link':dest.symlink_to(target)
    if outkind=='broken':dest.symlink_to(folder/'missing-target')
    if outkind=='loop':dest.symlink_to(dest)
    if outkind=='parent_missing':options['output']=str(folder/'missing/output.tsv')
 for path in folder.rglob('*'):
  if path.is_file() and not path.is_symlink():os.utime(path,ns=(1700000000000000000,1700000000000000000))
 return dict(command=cmd,options=options,helper_root=str(source.parent))
def inventory(folder):
 result={}
 for path in sorted(folder.rglob('*')):
  if path.name in ['stdout','stderr','request.json','result.json','arm.json']:continue
  s=path.lstat();v=dict(mode=stat.S_IMODE(s.st_mode),kind=stat.S_IFMT(s.st_mode),inode=s.st_ino,device=s.st_dev,nlink=s.st_nlink,bytes=s.st_size,mtime_ns=s.st_mtime_ns,ctime_ns=s.st_ctime_ns)
  if path.is_symlink():v['link']=os.readlink(path)
  elif path.is_file():v['sha256']=pin(path)['sha256']
  result[str(path.relative_to(folder))]=v
 return result
def normalized(value,root):return json.loads(json.dumps(value).replace(str(root),'$ARM'))
def comparable(before,after,root):
 result={}
 for name,v in after.items():
  if v['kind']==stat.S_IFDIR:continue
  item={k:x for k,x in v.items() if k not in ['inode','device','mtime_ns','ctime_ns']}
  if 'link' in item:item['link']=item['link'].replace(str(root),'$ARM');item['bytes']=len(item['link'].encode())
  if name in before:item['unchanged_identity']=before[name]==v
  result[name]=item
 return result
def execute(root,spec,a,native):
 request=setup(root,spec,a.producer);save(root/'request.json',request);before=inventory(root)
 if native:
  command=[str(a.binary),'biomart',request['command'],'--helper-root',request['helper_root']]
  for k,v in request['options'].items():command+=['--'+k,str(v)]
 else:command=[sys.executable,__file__,'oracle','--producer',str(a.producer),'--request',str(root/'request.json'),'--result',str(root/'result.json')]
 start=time.monotonic();process=subprocess.run(command,cwd=root,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
 (root/'stdout').write_bytes(process.stdout);(root/'stderr').write_bytes(process.stderr)
 lines=process.stdout.decode(errors='replace').splitlines();progress=[v for v in lines if v.startswith('IMPORTED ')]
 if native:
  if not process.returncode:
   remaining='\n'.join(v for v in lines if not v.startswith('IMPORTED '));value=dict(success=True,value=json.loads(remaining))
  else:
   diagnostic=process.stderr.decode(errors='replace').rstrip('\n');match=re.fullmatch(r'priva: (\w+): ([\s\S]*)',diagnostic)
   value=dict(success=False,category=match[1] if match else 'NATIVE_BOUNDARY',message=match[2] if match else diagnostic)
  save(root/'result.json',value)
 else:value=json.loads((root/'result.json').read_text())
 return_order=key_order(json.loads('\n'.join(v for v in lines if not v.startswith('IMPORTED ')))) if not process.returncode else None
 after=inventory(root);record=dict(return_dictionary_order=return_order,command=command,returncode=process.returncode,result=value,progress=progress,before=before,after=after,comparable=comparable(before,after,root),wall_seconds=time.monotonic()-start,stdout=pin(root/'stdout'),stderr=pin(root/'stderr'))
 save(root/'arm.json',record);return record

def key_order(value):
 result=[];pending=[([],value)]
 while pending:
  path,item=pending.pop()
  if isinstance(item,dict):
   result.append([path,list(item)])
   pending.extend((path+[key],child) for key,child in reversed(list(item.items())))
  elif isinstance(item,list):pending.extend((path+[i],child) for i,child in reversed(list(enumerate(item))))
 return result
def exact_value(left,right):
 pending=[(left,right)]
 while pending:
  a,b=pending.pop()
  if type(a) is not type(b):return False
  if isinstance(a,dict):
   if a.keys()!=b.keys():return False
   pending.extend((v,b[k]) for k,v in a.items())
  elif isinstance(a,list):
   if len(a)!=len(b):return False
   pending.extend(zip(a,b))
  elif a!=b:return False
 return True

def readback(path):
 record=json.loads(path.read_text());now=inventory(path.parent);assert set(now)==set(record['after'])
 drift={}
 for name,old in record['after'].items():
  current=now[name];fields=['kind','mode'] if old['kind']==stat.S_IFDIR else ['kind','mode','nlink','bytes','sha256','link']
  assert all(old.get(k)==current.get(k) for k in fields),name
  if old!=current:drift[name]=dict(producing=old,storage=current)
 return record,dict(arm=pin(path),storage_metadata_differences=drift,bytes_modes_links_exact=True)
def check(a):
 a.work_root.mkdir();producer=pin(a.producer);binary=pin(a.binary) if a.binary else None;fixtures=cases();prior=None;old={};proofs={};results=[]
 if a.reuse_report:
  prior=json.loads(a.reuse_report.read_text());assert prior['producer']['sha256']==producer['sha256'];old={v['name']:v for v in prior['cases']}
 for name,spec in fixtures.items():
  folder=a.work_root/name;folder.mkdir()
  if name in old and prior['fixtures'][name]==spec:
   python_path=Path(old[name]['python']);python,proof=readback(python_path);proofs[name]=proof
  else:
   python_path=folder/'python/arm.json';python=execute(python_path.parent,spec,a,False)
  if a.oracle_only:
   results.append(dict(name=name,exact=False,differences=['native_unmeasured'],python=str(python_path),native=None,python_reused=name in proofs));continue
  native_path=folder/'native/arm.json';native=execute(native_path.parent,spec,a,True)
  differences=[]
  for key in ['returncode','result','progress','comparable','return_dictionary_order']:
   if not exact_value(normalized(python[key],python_path.parent),normalized(native[key],native_path.parent)):differences.append(key)
  record=dict(name=name,exact=not differences,differences=differences,python=str(python_path),native=str(native_path),python_reused=name in proofs)
  save(folder/'comparison.json',record);results.append(record)
 assert pin(a.producer)==producer and (not a.binary or pin(a.binary)==binary)
 result=dict(oracle_only=a.oracle_only,producer=producer,binary=binary,harness=pin(__file__),fixtures=fixtures,cases=results,total=len(results),exact=sum(v['exact'] for v in results),mismatches=[v for v in results if not v['exact']],reuse_report=pin(a.reuse_report) if a.reuse_report else None,reuse_readbacks=proofs,source_stable=True)
 save(a.report,result);print(json.dumps(dict(total=result['total'],exact=result['exact'],mismatches=result['mismatches'])))
 return 0 if a.oracle_only or not result['mismatches'] else 1

def main():
 p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(required=True)
 q=s.add_parser('oracle');q.add_argument('--producer',type=Path,required=True);q.add_argument('--request',type=Path,required=True);q.add_argument('--result',type=Path,required=True);q.set_defaults(run=oracle)
 q=s.add_parser('check');q.add_argument('--producer',type=Path,required=True);q.add_argument('--binary',type=Path);q.add_argument('--work-root',type=Path,required=True);q.add_argument('--report',type=Path,required=True);q.add_argument('--reuse-report',type=Path);q.add_argument('--oracle-only',action='store_true');q.set_defaults(run=check)
 a=p.parse_args();sys.exit(a.run(a) or 0)
if __name__=='__main__':main()
