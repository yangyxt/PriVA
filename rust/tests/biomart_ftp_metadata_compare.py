#!/usr/bin/env python3
"""Exact offline BioMart raw-checksum and publisher-metadata comparisons.

check --producer SCRIPT --binary PRIVA --xml-reference XML --dataset-row JSON
      --work-root NEW --report JSON [--reuse-report OLD] [--oracle-only]
The XML reference is a retained pinned publisher file. No network is used.
"""
import argparse,copy,gzip,hashlib,importlib.util,json,os,re,shutil,stat,subprocess,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('ftp_compare_support',Path(__file__).with_name('biomart_ftp_compare.py'))
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
pin,save,module=base.pin,base.save,base.module

def oracle(a):
 m=module(a.producer);r=json.loads(a.request.read_text());o=r['options']
 try:
  if r['command']=='ftp-verify-raw':
   if 'name' in o:p=m.PINS[o['name']]
   else:
    p=json.loads(o['pin'])
    if len(p)==3 and isinstance(p[1],list):p[1]=tuple(p[1])
   value=m.verify_raw(Path(o['path']),p)
  else:value=m.metadata(Path(o['raw']),json.loads(Path(o['schemas']).read_text()),Path(o['output-dir']))
 except Exception as e:
  save(a.result,dict(success=False,category=type(e).__name__,message=str(e)));raise SystemExit(1)
 save(a.result,dict(success=True,value=value));print(json.dumps(value))

def cases():
 result={}
 def add(name,command,**spec):result[name]=dict(command=command,**spec)
 for kind in ['ordinary','empty','no_sha','wrong_sha','wrong_size','sha_before_sum','size_before_sha','missing','directory','link','broken','loop','path_parent_file','short_pin','long_pin','unknown_name']:
  add('raw_'+kind,'ftp-verify-raw',kind=kind)
 for suffix in ['.gz','.gz.part','.GZ','.gz.other','']:
  for kind in ['good','empty','multimember','padding','notgzip','magic1','header_truncated','deflate_truncated','trailer_truncated','badcrc','badsize','method','trailing','deflate_invalid']:
   add('raw_gzip_'+suffix.replace('.','_')+'_'+kind,'ftp-verify-raw',kind='gzip',suffix=suffix,gzip=kind)
 for suffix in ['.gz','.gz.part']:
  for kind in ['deflate_lengths','deflate_distance','deflate_codes']:add('raw_deflate_extra_'+suffix.replace('.','_')+'_'+kind,'ftp-verify-raw',kind='gzip',suffix=suffix,gzip=kind)
 for kind in ['real','empty','one','extra','noninteger_first','noninteger_second','invalid_utf8','exit','exit_utf8','signal','missing','not_executable','signed','underscores','unicode','wide','long_integer','wrong','newline','stderr','mutate_raw','mutate_gzip','remove_gzip']:
  add('sum_'+kind,'ftp-verify-raw',kind='sum',sum=kind)
 for kind in ['ordinary','no_final_lf','crlf','empty','duplicate','short','other_short','other_empty','other_escape','null','nonascii','escape','absent','wrong_name','short_selected','bad_assembly','bad_display','extra_fields']:
  add('metadata_dataset_'+kind,'ftp-metadata',dataset=kind)
 for kind in ['ordinary','no_final_lf','crlf','empty','duplicate','absent','short','long','other_short','other_escape','escape','id_escape','null_blob','null_compressed','mismatch','message_null','message_nonascii','compressed_empty','compressed_notgzip','compressed_truncated','compressed_badcrc','compressed_badsize','compressed_method','compressed_multimember','compressed_padding','compressed_deflate_invalid','xml_invalid','xml_wrong_identity','xml_pin_changed']:
  add('metadata_xml_'+kind,'ftp-metadata',xml=kind)
 for kind in ['deflate_lengths','deflate_distance','deflate_codes']:add('metadata_compressed_extra_'+kind,'ftp-metadata',xml='compressed_'+kind)
 for side in ['dataset','xml']:
  for kind in ['missing','directory','notgzip','truncated','badcrc','multimember']:
   add('metadata_'+side+'_file_'+kind,'ftp-metadata',input_side=side,input_kind=kind)
 for kind in ['dataset_missing','dataset_reordered','dataset_short_pair','dataset_long_pair','xml_missing','xml_reordered','xml_short_pair','xml_long_pair','dataset_before_missing_file','decode_before_xml_schema']:
  add('metadata_schema_'+kind,'ftp-metadata',schema=kind)
 for kind in ['new','same','changed','directory','link_same','link_changed','broken','loop','parent_missing','parent_file','parent_link','part_file']:
  add('metadata_view_'+kind,'ftp-metadata',view=kind)
 return result

def gz_variant(blob,kind):
 good=gzip.compress(blob,mtime=0)
 if kind in ['deflate_lengths','deflate_distance','deflate_codes']:
  payload={'deflate_lengths':b'\x01\x01\x00\x00\x00X','deflate_distance':b'\x03\x02\x00','deflate_codes':b'\x05\x00\x00\x00'+b'\0'*16}[kind]
  return good[:10]+payload+good[-8:]
 if kind=='good':return good
 if kind=='empty':return b''
 if kind=='multimember':return gzip.compress(blob[:len(blob)//2],mtime=0)+gzip.compress(blob[len(blob)//2:],mtime=0)
 if kind=='padding':return good+b'\0'*8193
 if kind=='notgzip':return b'not gzip'
 if kind=='magic1':return b'\x1f'
 if kind=='header_truncated':return good[:5]
 if kind=='deflate_truncated':return good[:12]
 if kind=='trailer_truncated' or kind=='truncated':return good[:-5]
 if kind=='trailing':return good+b'bad'
 data=bytearray(good)
 if kind=='badcrc':data[-8]^=1
 elif kind=='badsize':data[-4]^=1
 elif kind=='method':data[2]=9
 elif kind=='deflate_invalid':data[10]=7
 else:raise AssertionError(kind)
 return bytes(data)

def mysql(value):
 if value is None:return b'\\N'
 if not isinstance(value,bytes):value=str(value).encode('ascii')
 return value.replace(b'\\',b'\\\\').replace(b'\0',b'\\0').replace(b'\b',b'\\b').replace(b'\n',b'\\n').replace(b'\r',b'\\r').replace(b'\t',b'\\t').replace(b'\x1a',b'\\Z')
def dump(rows):return b''.join(b'\t'.join(mysql(v) for v in row)+b'\n' for row in rows)

def setup(folder,spec,a):
 folder.mkdir();source=folder/'source/scripts';source.mkdir(parents=True)
 for n in ['biomart_ftp_nm_to_enst.py','biomart_nm_to_enst.py']:shutil.copyfile(a.producer.with_name(n),source/n)
 options={};env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1');cmd=spec['command']
 if cmd=='ftp-verify-raw':
  kind=spec['kind'];path=folder/('raw'+spec.get('suffix','.bin'));content=b'ordinary\0raw bytes\n'
  if kind=='sum' and spec['sum'].endswith('_gzip'):
   path=folder/'raw.gz';content=gzip.compress(content,mtime=0)
  if kind=='empty':content=b''
  if kind=='gzip':content=gz_variant(content,spec['gzip'])
  path.write_bytes(content);size=len(content);sha=hashlib.sha256(content).hexdigest();bsd=None
  if kind=='sum':
   bsd=list(map(int,subprocess.check_output(['sum','-r','--',str(path)],text=True).split()[:2]));which=spec['sum']
   if which!='real':
    tools=folder/'tools';tools.mkdir();program=tools/'sum';env['PATH']=str(tools)
    first,second=map(str,bsd)
    raw={'empty':b'','one':first.encode(),'extra':f'{first} {second} ignored other'.encode(),'noninteger_first':b'bad 1','noninteger_second':f'{first} bad'.encode(),'invalid_utf8':b'\xff 1','exit_utf8':b'\xff 1','signed':f'+{first} +{second}'.encode(),'underscores':f'{first[0]}_{first[1:]} {second}'.encode(),'unicode':''.join(chr(0x660+int(c)) if c.isdigit() else c for c in f'{first} {second}').encode(),'wide':b'9'*80+b' 1','long_integer':b'1'*4301+b' 1','wrong':b'0 0','newline':f'{first}\r\n{second}\r'.encode()}.get(which,f'{first} {second}'.encode())
    effect='os.write(2,b"sum diagnostic\\n")\n' if which=='stderr' else 'open(sys.argv[3],"wb").write(b"changed bytes")\n' if which in ['mutate_raw','mutate_gzip'] else 'os.unlink(sys.argv[3])\n' if which=='remove_gzip' else ''
    program.write_text('#!'+sys.executable+'\nimport os,sys,signal\nos.write(1,bytes.fromhex('+repr(raw.hex())+'))\n'+effect+('os.kill(os.getpid(),signal.SIGTERM)\n' if which=='signal' else 'sys.exit(7)\n' if which in ['exit','exit_utf8'] else ''))
    program.chmod(0o644 if which=='not_executable' else 0o755)
    if which=='missing':program.unlink()
  p=[size,bsd,sha]
  if kind=='no_sha':p[2]=None
  if kind in ['wrong_sha','sha_before_sum']:p[2]='0'*64
  if kind=='sha_before_sum':p[1]=[0,0]
  if kind in ['wrong_size','size_before_sha']:p[0]+=1;p[2]='0'*64
  if kind in ['link','broken','loop']:
   original=folder/'target';path.rename(original);path.symlink_to(original if kind=='link' else folder/'absent' if kind=='broken' else path)
  if kind=='missing':path.unlink()
  if kind=='directory':path.unlink();path.mkdir()
  if kind=='path_parent_file':path=path/'child'
  if kind=='short_pin':p=p[:2]
  if kind=='long_pin':p+=[1]
  options={'path':str(path),'pin':json.dumps(p)}
  if kind=='unknown_name':options.pop('pin');options['name']='unknown'
 else:
  raw=folder/'raw';raw.mkdir();out=folder/'view';out.mkdir();row=copy.deepcopy(a.reference_row)
  blob=a.reference_xml;compressed=gzip.compress(blob,mtime=0);fields=[row[0],blob,compressed,'unused'];rows=[row];xmlrows=[fields]
  kind=spec.get('dataset','ordinary')
  if kind=='empty':rows=[]
  if kind=='duplicate':rows=[row,row]
  if kind=='short':rows=[['short']]
  if kind=='other_short':rows=[['short'],row]
  if kind=='other_empty':rows=[[],row]
  if kind=='other_escape':rows=[['x','unselected','bad\\'],row]
  if kind=='null':row[3]=None
  if kind=='nonascii':row[3]=b'\xff'
  if kind=='absent':rows=[['x','other']]
  if kind=='wrong_name':row[1]='other'
  if kind=='short_selected':rows=[row[:3]]
  if kind=='bad_assembly':row[6]='GRCh37'
  if kind=='bad_display':row[2]='changed'
  if kind=='extra_fields':row+=['extra']
  dataset=dump(rows)
  if kind=='escape':dataset=dataset.replace(b'\t'+mysql(row[3])+b'\t',b'\tbad\\\t',1)
  if kind=='no_final_lf':dataset=dataset.rstrip(b'\n')
  if kind=='crlf':dataset=dataset.replace(b'\n',b'\r\n')
  kind=spec.get('xml','ordinary')
  if kind=='empty':xmlrows=[]
  if kind=='duplicate':xmlrows=[fields,fields]
  if kind=='absent':fields[0]='absent'
  if kind=='short':xmlrows=[fields[:3]]
  if kind=='long':fields+=['extra']
  if kind=='other_short':xmlrows=[['other'],fields]
  if kind=='other_escape':xmlrows=[['other',b'bad\\',None,None],fields]
  if kind=='null_blob':fields[1]=None
  if kind=='null_compressed':fields[2]=None
  if kind=='mismatch':fields[2]=gzip.compress(b'other',mtime=0)
  if kind=='message_null':fields[3]=None
  if kind=='message_nonascii':fields[3]=b'\xff'
  if kind.startswith('compressed_'):fields[2]=gz_variant(blob,kind.removeprefix('compressed_'))
  if kind in ['xml_invalid','xml_wrong_identity','xml_pin_changed']:
   fields[1]=b'<unclosed' if kind=='xml_invalid' else b'<Other/>' if kind=='xml_wrong_identity' else blob+b'\n';fields[2]=gzip.compress(fields[1],mtime=0)
  xml=dump(xmlrows)
  if kind=='escape':xml=xml[:-1]+b'\\\n'
  if kind=='id_escape':xml=xml.replace(mysql(row[0])+b'\t',b'\\'+mysql(row[0])+b'\t',1)
  if kind=='no_final_lf':xml=xml.rstrip(b'\n')
  if kind=='crlf':xml=xml.replace(b'\n',b'\r\n')
  schema=base.schemas();kind=spec.get('schema','')
  side='meta_conf__dataset__main' if kind.startswith('dataset_') else 'meta_conf__xml__dm'
  if kind.endswith('_missing'):del schema[side]
  if kind.endswith('_reordered'):schema[side].reverse()
  if kind.endswith('_short_pair'):schema[side][0]=['one']
  if kind.endswith('_long_pair'):schema[side][0]=['one','two','three']
  if kind=='dataset_before_missing_file':schema['meta_conf__dataset__main'].reverse()
  if kind=='decode_before_xml_schema':schema['meta_conf__xml__dm'].reverse();dataset=dataset.replace(mysql(row[3]),b'\xff',1)
  paths={k:raw/('meta_conf__'+k+'__'+('main' if k=='dataset' else 'dm')+'.txt.gz') for k in ['dataset','xml']}
  paths['dataset'].write_bytes(gzip.compress(dataset,mtime=0));paths['xml'].write_bytes(gzip.compress(xml,mtime=0))
  if spec.get('input_kind'):
   path=paths[spec['input_side']];kind=spec['input_kind'];payload=dataset if spec['input_side']=='dataset' else xml
   if kind=='missing':path.unlink()
   elif kind=='directory':path.unlink();path.mkdir()
   else:path.write_bytes(gz_variant(payload,kind))
  if spec.get('schema')=='dataset_before_missing_file':paths['dataset'].unlink()
  view=out/'hsapiens_gene_ensembl.dataset.xml';kind=spec.get('view','new')
  if kind in ['same','changed']:view.write_bytes(blob if kind=='same' else b'old')
  if kind=='directory':view.mkdir()
  if kind.startswith('link_'):
   target=folder/'target.xml';target.write_bytes(blob if kind=='link_same' else b'old');view.symlink_to(target)
  if kind=='broken':view.symlink_to(folder/'missing.xml')
  if kind=='loop':view.symlink_to(view)
  if kind=='parent_missing':out.rmdir()
  if kind=='parent_file':out.rmdir();out.write_bytes(b'file')
  if kind=='parent_link':out.rename(folder/'real_view');out.symlink_to(folder/'real_view',target_is_directory=True)
  if kind=='part_file':view.with_name(view.name+'.part').write_bytes(b'unrelated preserved')
  schema_path=folder/'schemas.json';save(schema_path,schema);options={'raw':str(raw),'schemas':str(schema_path),'output-dir':str(out)}
 return dict(command=cmd,options=options,helper_root=str(source.parent)),env

def execute(root,spec,a,native):
 request,env=setup(root,spec,a);save(root/'request.json',request);before=base.inventory(root)
 if native:
  command=[str(a.binary),'biomart',request['command'],'--helper-root',request['helper_root']]
  for k,v in request['options'].items():command+=['--'+k,str(v)]
 else:command=[sys.executable,__file__,'oracle','--producer',str(a.producer),'--request',str(root/'request.json'),'--result',str(root/'result.json')]
 started=time.monotonic();process=subprocess.run(command,cwd=root,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=90)
 (root/'stdout').write_bytes(process.stdout);(root/'stderr').write_bytes(process.stderr)
 if native:
  if not process.returncode:value=dict(success=True,value=json.loads(process.stdout))
  else:
   diagnostic=process.stderr.decode(errors='replace').rstrip('\n');match=re.search(r'(?:^|\n)priva: (\w+): ([\s\S]*)',diagnostic)
   value=dict(success=False,category=match[1] if match else 'NATIVE_BOUNDARY',message=match[2] if match else diagnostic)
  save(root/'result.json',value)
 else:value=json.loads((root/'result.json').read_text())
 after=base.inventory(root);result=dict(command=command,returncode=process.returncode,result=value,return_dictionary_order=base.key_order(json.loads(process.stdout)) if not process.returncode else None,before=before,after=after,comparable=base.comparable(before,after,root),stdout=pin(root/'stdout'),stderr=pin(root/'stderr'),wall_seconds=time.monotonic()-started)
 save(root/'arm.json',result);return result

def check(a):
 a.work_root.mkdir();producer=pin(a.producer);binary=pin(a.binary) if a.binary else None
 reference=pin(a.xml_reference);assert reference['sha256']=='a5afb502ad97bada17f9fd4cdc3886be18ba8b237f3b896568793eb0cdd74fa7'
 a.reference_xml=a.xml_reference.read_bytes();a.reference_row=json.loads(a.dataset_row);fixtures=cases();old={};prior=None;proofs={};results=[]
 if a.reuse_report:
  prior=json.loads(a.reuse_report.read_text());assert prior['producer']['sha256']==producer['sha256'];assert prior['xml_reference']==reference and prior['dataset_row']==a.reference_row
  old={v['name']:v for v in prior['cases']}
 for name,spec in fixtures.items():
  folder=a.work_root/name;folder.mkdir()
  if name in old and prior['fixtures'][name]==spec:
   original_path=Path(old[name]['python']);assert pin(original_path)==old[name]['python_record'];original,proof=base.readback(original_path)
   for stream in ['stdout','stderr']:assert pin(original[stream]['path'])==original[stream]
   proofs[name]=proof
  else:
   original_path=folder/'python/arm.json';original=execute(original_path.parent,spec,a,False)
  if a.oracle_only:results.append(dict(name=name,exact=False,differences=['native_unmeasured'],python=str(original_path),native=None,python_reused=name in proofs,python_record=pin(original_path)));continue
  native_path=folder/'native/arm.json';native=execute(native_path.parent,spec,a,True);differences=[]
  original_external=Path(original['stderr']['path']).read_bytes()
  native_external=Path(native['stderr']['path']).read_bytes().split(b'priva: ',1)[0]
  if original_external.replace(os.fsencode(original_path.parent),b'$ARM')!=native_external.replace(os.fsencode(native_path.parent),b'$ARM'):differences.append('external_stderr')
  for key in ['returncode','result','comparable','return_dictionary_order']:
   if not base.exact_value(base.normalized(original[key],original_path.parent),base.normalized(native[key],native_path.parent)):differences.append(key)
  record=dict(name=name,exact=not differences,differences=differences,python=str(original_path),native=str(native_path),python_reused=name in proofs,python_record=pin(original_path),native_record=pin(native_path));save(folder/'comparison.json',record);results.append(record)
 assert pin(a.producer)==producer and (not a.binary or pin(a.binary)==binary) and pin(a.xml_reference)==reference
 value=dict(producer=producer,binary=binary,xml_reference=reference,dataset_row=a.reference_row,harness=pin(__file__),support=pin(Path(base.__file__)),fixtures=fixtures,cases=results,total=len(results),exact=sum(v['exact'] for v in results),mismatches=[v for v in results if not v['exact']],reuse_readbacks=proofs,source_stable=True,oracle_only=a.oracle_only)
 save(a.report,value);print(json.dumps(dict(total=value['total'],exact=value['exact'],mismatches=value['mismatches'])))
 return 0 if a.oracle_only or not value['mismatches'] else 1

def main():
 p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(required=True)
 q=s.add_parser('oracle');q.add_argument('--producer',type=Path,required=True);q.add_argument('--request',type=Path,required=True);q.add_argument('--result',type=Path,required=True);q.set_defaults(run=oracle)
 q=s.add_parser('check');q.add_argument('--producer',type=Path,required=True);q.add_argument('--binary',type=Path);q.add_argument('--xml-reference',type=Path,required=True);q.add_argument('--dataset-row',required=True);q.add_argument('--work-root',type=Path,required=True);q.add_argument('--report',type=Path,required=True);q.add_argument('--reuse-report',type=Path);q.add_argument('--oracle-only',action='store_true');q.set_defaults(run=check)
 a=p.parse_args();sys.exit(a.run(a) or 0)
if __name__=='__main__':main()
