#!/usr/bin/env python3
"""Source differential checks for optional BioMart receipt publication.

check --producer SCRIPT --binary PRIVA --work-root NEW --report JSON
      [--reuse-report OLD] [--oracle-only]
All file writes are private fixtures. No network or installed cache is used.
Random eight-character tempfile tokens and source-root paths have explicit
per-arm bindings; receipt contents, modes and preservation checks stay exact.
"""
import argparse,hashlib,importlib.util,json,math,os,re,shutil,stat,subprocess,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('publication_support',Path(__file__).with_name('biomart_ftp_compare.py'))
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
pin=base.pin
def save(path,value):
    with Path(path).open('x') as f:f.write(json.dumps(value,ensure_ascii=True,separators=(',',':'))+'\n')
def module(path):
    sys.path.insert(0,str(path.parent));s=importlib.util.spec_from_file_location('publication_original',path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def oracle(a):
    r=json.loads(a.request.read_text());o=r['options'];m=module(Path(r['helper_root'])/'scripts/biomart_nm_to_enst.py')
    method=m.FTP if o.get('method')=='ftp' else m
    try:
        record=json.loads(Path(o['record']).read_text()) if 'record' in o else None
        if r['command']=='atomic-json':value=method.atomic_json(Path(o['output']),record)
        elif r['command']=='runtime-receipt':value=str(method.runtime_receipt(record,Path(o['workspace']),o['basename']))
        elif r['command']=='external-path':
            if o['method']=='service':value=m.external_acquisition_path(Path(o['path']))
            else:
                # The FTP function's original first check, without creating a
                # receipt solely to inspect the pure path admission step.
                value=m.FTP.require(not Path(o['path']).resolve().is_relative_to(Path(m.FTP.__file__).resolve().parents[1]),'BioMart acquisition workspace must be outside the checkout')
        elif r['command']=='ftp-portable-provenance':value=m.FTP.portable_provenance(record)
        elif r['command']=='ftp-make-provenance':value=m.FTP.make_provenance(Path(o['mapping']),m.validate_table,o['operation'])
        else:raise AssertionError(r['command'])
    except Exception as e:
        save(a.result,dict(success=False,category=type(e).__name__,message=str(e)));raise SystemExit(1)
    save(a.result,dict(success=True,value=value));print(json.dumps(value))
def cases():
    result={}
    def add(name,command,**spec):result[name]=dict(command=command,**spec)
    records=['ordinary','empty','array','null','bool','int','float','floats','wide','digits4300','surrogates','surrogate_keys','duplicates','unicode_order','negative_zero','invalid','utf8','bom','digits4301']
    records+=['depth'+str(n) for n in [399,500,980,993,994,995,1200]]
    outputs=['new','file','directory','link','broken','loop','hardlink','part_file','part_directory','part_link','part_broken','part_loop','missing_parent','parent_file','readonly_parent','readonly_file','dot','empty','root']
    for method in ['service','ftp']:
        for kind in records:add(method+'_record_'+kind,'atomic-json',method=method,record=kind)
        for kind in outputs:add(method+'_output_'+kind,'atomic-json',method=method,output=kind)
        for kind in ['outside','missing','file','inside','inside_missing','link_outside','link_inside','broken','loop','parent_file','dots_outside','dots_inside','readonly']:
            add(method+'_external_'+kind,'external-path',method=method,path=kind)
            add(method+'_runtime_'+kind,'runtime-receipt',method=method,path=kind)
        for kind in ['empty','space','unicode','long','limit','nested','nested_missing','dots','absolute']:
            add(method+'_basename_'+kind,'runtime-receipt',method=method,basename=kind)
        for kind in ['wide','surrogates','surrogate_keys','depth500','depth993','depth994','depth995','depth1200']:
            add(method+'_runtime_record_'+kind,'runtime-receipt',method=method,record=kind)
    for kind in records:add('portable_'+kind,'ftp-portable-provenance',record=kind)
    for kind in ['null','bool','int','float','array','string']:add('portable_root_'+kind,'ftp-portable-provenance',root=kind)
    for kind in ['missing','empty','directory','broken','loop','small','wrong_hash','parent_file']:
        add('provenance_'+kind,'ftp-make-provenance',mapping=kind)
    return result
def record_text(kind):
    values={'ordinary':'{"z": [true, null, -3], "a": "hello\\nworld", "empty": {}}','empty':'{}','array':'[0, {}, []]','null':'null','bool':'true','int':'7','float':'1.5','floats':'[1e-5, 1e16, -0.0, 0.0001, 1.2345678901234567, 1e999, -1e999, NaN]','wide':'[1461501637330902918203684832716283019655932542976, -1461501637330902918203684832716283019655932542976]','surrogates':'["\\ud800", "\\udfff", "\\ud800x", "\\ud800\\ud800", "\\ud83d\\ude00"]','surrogate_keys':'{"\\ud800": 1, "\\u0000": 2, "\\udfff": 3, "\\ud800": 4}','duplicates':'{"z":1,"a":{"x":1,"x":2},"z":3}','unicode_order':'{"😀": 0, "\\ue000": 1, "\\ud800": 2, "é": 3}','negative_zero':'[-0, -0.0, -0e0]','invalid':'{"valid": 1, "x": [}','utf8':b'{"x":"\xff"}','bom':'\ufeff{}'}
    if kind.startswith('depth'):
        n=int(kind[5:]);return '['*n+'0'+']'*n
    if kind.startswith('digits'):return '9'*int(kind[6:])
    return values[kind]
def setup(folder,spec,a):
    folder.mkdir();source=folder/'source/scripts';source.mkdir(parents=True)
    for name in ['biomart_nm_to_enst.py','biomart_ftp_nm_to_enst.py']:shutil.copyfile(a.producer.with_name(name),source/name)
    cmd=spec['command'];options={};env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',SLURM_JOB_ID='12345')
    if 'method' in spec:options['method']=spec['method']
    if cmd in ['atomic-json','runtime-receipt','ftp-portable-provenance']:
        text=record_text(spec.get('record','ordinary'))
        if cmd=='ftp-portable-provenance':
            # All source-selected values are retained without interpreting
            # them; unrelated valid JSON is deliberately ignored.
            fields=dict(a.facts['ftp_specification'],bytes=28785905,sha256='snapshot',statistics=None,extra='ignored',metadata_kind='wrong')
            encoded=json.dumps(fields);text=encoded.replace('"statistics": null','"statistics": '+(text.decode('latin1') if isinstance(text,bytes) else text))
            if spec.get('record') in ['invalid','utf8','bom','digits4301']:text=record_text(spec['record'])
            if 'root' in spec:text={'null':'null','bool':'true','int':'3','float':'1.0','array':'[]','string':'"text"'}[spec['root']]
        record=folder/'record.json';record.write_bytes(text if isinstance(text,bytes) else text.encode());options['record']=str(record)
    if cmd=='atomic-json':
        output=folder/'result-data.json';target=folder/'prior';target.write_bytes(b'previous target\n');kind=spec.get('output','new');dest=Path(str(output)+'.part') if kind.startswith('part_') else output
        subtype=kind.removeprefix('part_')
        if subtype in ['file','readonly_file']:dest.write_bytes(b'existing\n')
        if subtype=='readonly_file':dest.chmod(0o400)
        if subtype=='directory':dest.mkdir()
        if subtype in ['link','broken','loop']:dest.symlink_to(target if subtype=='link' else folder/'missing' if subtype=='broken' else dest)
        if subtype=='hardlink':os.link(target,dest)
        if kind=='missing_parent':output=folder/'absent/output.json'
        if kind=='parent_file':output=target/'output.json'
        if kind=='readonly_parent':parent=folder/'readonly';parent.mkdir();parent.chmod(0o500);output=parent/'output.json'
        if kind=='dot':output=Path('.')
        if kind=='empty':output=Path('')
        # Absolute filesystem root is never a publication destination. A
        # private path ending in / is the safe directory-equivalent probe.
        if kind=='root':output=folder
        options['output']=str(output)
    if cmd in ['runtime-receipt','external-path']:
        work=folder/'workspace';work.mkdir();inside=source.parent/'workspace';inside.mkdir();kind=spec.get('path','outside')
        if kind=='missing':work=folder/'absent'
        if kind=='file':work=folder/'file';work.write_text('file')
        if kind=='inside':work=inside
        if kind=='inside_missing':work=inside/'absent'
        if kind in ['link_outside','link_inside','broken','loop']:
            alias=folder/'alias';alias.symlink_to(work if kind=='link_outside' else inside if kind=='link_inside' else folder/'absent' if kind=='broken' else alias);work=alias
        if kind=='parent_file':parent=folder/'file';parent.write_text('file');work=parent/'absent'
        if kind=='dots_outside':work=inside/'../../workspace'
        if kind=='dots_inside':work=folder/'workspace/../source/workspace'
        if kind=='readonly':work.chmod(0o500)
        if cmd=='external-path':options['path']=str(work)
        else:
            basename='mapping.tsv';b=spec.get('basename','ordinary')
            if b=='empty':basename=''
            if b=='space':basename='mapping file'
            if b=='unicode':basename='表-😀'
            if b=='long':basename='x'*240
            if b=='limit':basename='x'*224
            if b in ['nested','nested_missing','dots']:
                if b!='nested_missing':(folder/'workspace/sub').mkdir()
                basename='sub/../name' if b=='dots' else 'sub/name'
            if b=='absolute':basename=str(folder/'absolute')
            options.update(workspace=str(work),basename=basename)
    if cmd=='ftp-make-provenance':
        mapping=folder/'mapping.tsv';kind=spec['mapping']
        if kind in ['small','empty','wrong_hash']:mapping.write_bytes(b'x' if kind=='small' else b'')
        if kind=='wrong_hash':
            with mapping.open('wb') as f:f.truncate(28785905)
        if kind=='directory':mapping.mkdir()
        if kind in ['broken','loop']:mapping.symlink_to(folder/'absent' if kind=='broken' else mapping)
        if kind=='parent_file':mapping.write_text('file');mapping=mapping/'child'
        options.update(mapping=str(mapping),operation='registered_existing_exact_snapshot')
    for path in folder.rglob('*'):
        if path.is_file() and not path.is_symlink():os.utime(path,ns=(1700000000000000000,1700000000000000000))
    return dict(command=cmd,options=options,helper_root=str(source.parent)),env
TOKEN=re.compile(r'(?<=\.acquisition\.)[a-z0-9_]{8}(?=\.json(?:\.part)?)')
def normalized(value,root):return json.loads(TOKEN.sub('$TOKEN',json.dumps(value,ensure_ascii=True).replace(str(root),'$ARM')))
def equal(left,right):
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
        elif isinstance(a,float) and math.isnan(a) and math.isnan(b):continue
        elif a!=b:return False
    return True
def execute(root,spec,a,native):
    request,env=setup(root,spec,a);save(root/'request.json',request);before=base.inventory(root)
    if native:
        command=[str(a.binary),'biomart',request['command'],'--helper-root',request['helper_root']]
        for k,v in request['options'].items():command+=['--'+k,str(v)]
    else:command=[sys.executable,__file__,'oracle','--producer',str(a.producer),'--request',str(root/'request.json'),'--result',str(root/'result.json')]
    started=time.monotonic();process=subprocess.run(command,cwd=root,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=60,env=env)
    (root/'stdout').write_bytes(process.stdout);(root/'stderr').write_bytes(process.stderr)
    if native:
        if not process.returncode:value=dict(success=True,value=json.loads(process.stdout) if process.stdout else None)
        else:
            diagnostic=process.stderr.decode(errors='replace').rstrip('\n');match=re.fullmatch(r'priva: (\w+): ([\s\S]*)',diagnostic)
            value=dict(success=False,category=match[1] if match else 'NATIVE_BOUNDARY',message=match[2] if match else diagnostic)
        save(root/'result.json',value)
    else:value=json.loads((root/'result.json').read_text())
    order=base.key_order(json.loads(process.stdout)) if process.stdout and not process.returncode else None
    if value.get('success') and value.get('value') is None:order=None
    after=base.inventory(root);comp=base.comparable(before,after,root)
    record=dict(returncode=process.returncode,result=value,return_dictionary_order=order,before=before,after=after,comparable=comp,command=command,stdout=pin(root/'stdout'),stderr=pin(root/'stderr'),wall_seconds=time.monotonic()-started)
    save(root/'arm.json',record);return record
def check(a):
    a.work_root.mkdir();a.facts=json.loads((a.producer.parent.parent/'rust/src/biomart/facts.json').read_text());fixtures=cases();results=[];old={};proofs={}
    producer=pin(a.producer);binary=pin(a.binary) if a.binary else None
    if a.reuse_report:
        prior=json.loads(a.reuse_report.read_text());assert prior['producer']['sha256']==producer['sha256'];old={v['name']:v for v in prior['cases']}
    for name,spec in fixtures.items():
        if a.only and name not in a.only:continue
        folder=a.work_root/name;folder.mkdir()
        if name in old and prior['fixtures'][name]==spec:
            original_path=Path(old[name]['python']);assert pin(original_path)==old[name]['python_record'];original,proof=base.readback(original_path)
            for stream in ['stdout','stderr']:assert pin(original[stream]['path'])==original[stream]
            proofs[name]=proof
        else:original_path=folder/'python/arm.json';original=execute(original_path.parent,spec,a,False)
        if a.oracle_only:results.append(dict(name=name,python=str(original_path),python_record=pin(original_path)));continue
        native_path=folder/'native/arm.json';native=execute(native_path.parent,spec,a,True);differences=[]
        for key in ['returncode','result','comparable','return_dictionary_order']:
            if not equal(normalized(original[key],original_path.parent),normalized(native[key],native_path.parent)):differences.append(key)
        # json.dump's depth limit is CPython's call stack. Python stops with
        # RecursionError and a staging prefix; the native writer has no such limit.
        stated=bool(differences) and original['result'].get('category')=='RecursionError' and native['result'].get('success')
        record=dict(name=name,exact=not differences,stated_difference='python_recursion_limit' if stated else None,differences=differences,python=str(original_path),native=str(native_path),python_record=pin(original_path),native_record=pin(native_path),python_reused=name in proofs)
        save(folder/'comparison.json',record);results.append(record)
    assert pin(a.producer)==producer and (not a.binary or pin(a.binary)==binary)
    value=dict(producer=producer,binary=binary,harness=pin(__file__),support=pin(Path(base.__file__)),fixtures=fixtures,cases=results,total=len(results),exact=sum(v.get('exact',False) for v in results),stated=[v['name'] for v in results if v.get('stated_difference')],mismatches=[v for v in results if not v.get('exact',False) and not v.get('stated_difference')],reuse_readbacks=proofs,source_stable=True,oracle_only=a.oracle_only)
    save(a.report,value);print(json.dumps(dict(total=value['total'],exact=value['exact'],stated=value['stated'],mismatches=value['mismatches'])))
    return 0 if a.oracle_only or not value['mismatches'] else 1
def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(required=True)
    q=s.add_parser('oracle');q.add_argument('--producer',type=Path,required=True);q.add_argument('--request',type=Path,required=True);q.add_argument('--result',type=Path,required=True);q.set_defaults(run=oracle)
    q=s.add_parser('check');q.add_argument('--producer',type=Path,required=True);q.add_argument('--binary',type=Path);q.add_argument('--work-root',type=Path,required=True);q.add_argument('--report',type=Path,required=True);q.add_argument('--reuse-report',type=Path);q.add_argument('--oracle-only',action='store_true');q.add_argument('--only',action='append');q.set_defaults(run=check)
    a=p.parse_args();sys.exit(a.run(a) or 0)
if __name__=='__main__':main()
