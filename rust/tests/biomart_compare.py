#!/usr/bin/env python3
"""Source-bound BioMart procedure comparison; no network or installed writes.

check --producer SCRIPT --binary PRIVA --work-root NEW --report JSON
oracle --producer SCRIPT --request JSON --result JSON
The native CLI uses named arguments; the oracle calls the same original
procedure. Compare exits/categories/diagnostics, complete JSON/ordered arrays,
file bytes/modes/link targets and mutation order visible in partial outputs.
Raw per-arm metadata is retained; no storage-mode or timestamp substitution.
"""
from __future__ import annotations
import argparse, copy, hashlib, importlib.util, json, math, os, re, shutil, stat, struct, subprocess, sys
from pathlib import Path
sys.dont_write_bytecode=True
HEAD='Gene stable ID\tTranscript stable ID\tHGNC symbol\tUniProtKB/Swiss-Prot ID\tRefSeq mRNA ID\n'
ROW='ENSG1\tENST1\tGENE\tP1\tNM_1\n'
BODY=HEAD+ROW
GR37='ensembl-grch37-mart116'; FTP='ensembl-grch38-mart116-ftp'
BASENAME='ensembl_GRCh37_mart116_nm_to_enst.tsv'

def pin(path):
    path=Path(path)
    with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    return dict(path=str(path),bytes=path.stat().st_size,sha256=digest)
def save(path,value):
    with Path(path).open('x') as stream:json.dump(value,stream,indent=2,sort_keys=True,default=str);stream.write('\n')
def module(path):
    spec=importlib.util.spec_from_file_location('biomart_oracle',path);value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value
def encode(value):return value.hex() if isinstance(value,bytes) else value.encode().hex()
def oracle(args):
    b=module(args.producer);r=json.loads(args.request.read_text());o=r['options'];command=r['command']
    p=lambda k:Path(o[k]);inputs=[Path(v) for v in r.get('inputs',[])];value=None
    try:
        if command=='identity':value=b.mapping_identity(p('mapping'),inputs)
        elif command=='receipt-check':b.check_receipt(p('artifact'),p('mapping'),inputs)
        elif command=='receipt-write':
            artifact=p('artifact')
            if not artifact.is_file() or not artifact.stat().st_size:raise ValueError(f'Scientific command did not create a nonempty output; no mapping receipt written: {artifact}')
            expected=json.loads(o['expected-identity'])
            if expected!=b.mapping_identity(p('mapping'),inputs):raise ValueError('Mapping content/provenance changed during scientific generation; output has no valid mapping receipt')
            b.atomic_json(Path(str(artifact)+'.nm_to_enst.json'),dict(mapping_identity=expected,artifact_bytes=artifact.stat().st_size,artifact_sha256=b.sha256(artifact)))
        elif command=='artifact-receipt':value=b._artifact_receipt(p('artifact'))
        elif command=='artifact-route':
            value={'artifact':str(b.artifact_route(artifact=p('artifact'),mapping=p('mapping'),cache_root=p('cache-root'),output_dir=p('output-dir'),assembly=o['assembly'],kind=o['kind'],inputs=inputs))}
        elif command=='changed-target':value=str(b.changed_mapping_target(identity=json.loads(o['identity']),cache_root=p('cache-root'),output_dir=p('output-dir'),assembly=o['assembly'],kind=o['kind']))
        elif command=='verify':value=b.verify_pinned(p('mapping'),o.get('assembly'))
        elif command=='validate-table':value=b.validate_table(p('table'),raw=r.get('raw',False),output=p('output') if 'output' in o else None)
        elif command=='parse-metadata':value=b.parse_metadata(p('registry').read_text(),p('datasets').read_text(),o.get('source',GR37))
        elif command=='source':
            if 'assembly' in o:value=b.source_for_assembly(o['assembly'],o.get('vep-release'))
            elif 'basename' in o:value=b.source_for_basename(o['basename'])
            else:value=b._source(o.get('source',GR37))
        elif command=='snapshot-identity':value=b.snapshot_identity(p('mapping'),o.get('source',GR37))
        elif command=='is-bundled':value=b.is_bundled_mapping(p('mapping'))
        elif command=='ftp-specification':value=b.FTP.specification()
        elif command=='ftp-table-identity':value=b.FTP.table_identity(p('mapping'),b.validate_table)
        elif command=='ftp-verify-record':value=b.FTP.verify_record(p('mapping'),json.loads(p('record').read_text()),b.validate_table,o.get('assembly'))
        elif command=='ftp-parse-xml':value=b.FTP.parse_xml(p('xml').read_bytes(),json.loads(o['dataset-row']))
        else:raise AssertionError(command)
    except Exception as exc:
        save(args.result,dict(success=False,category=type(exc).__name__,message=str(exc)));raise SystemExit(1)
    save(args.result,dict(success=True,value=value))
    if value is not None:print(json.dumps(value,sort_keys=True,default=str))

def cases():
    result={}
    def add(name,command,**kw):result[name]=dict(command=command,**kw)
    tables={'normal':BODY,'duplicates':BODY+ROW+ROW.replace('ENST1','ENST2'),
        'no_nm':BODY.replace('NM_1',''),'other_refseq':BODY+ROW.replace('NM_1','NC_1.2'),
        'empty':'','header_only':HEAD,'wrong_header':'wrong\n'+ROW,'short':HEAD+'ENSG1\tENST1\n',
        'long':BODY+ROW.rstrip('\n')+'\textra\n','blank_row':BODY+'\n','blank_before':'\n'+BODY,
        'nm_version':BODY.replace('NM_1','NM_1.2'),'nm_letters':BODY.replace('NM_1','NM_x'),
        'xm_version':BODY.replace('NM_1','XM_1.2')+ROW,'bad_refseq':BODY.replace('NM_1','bad'),
        'gene_prefix':BODY.replace('ENSG1','ENSG'),'transcript_prefix':BODY.replace('ENST1','ENSTx'),
        'unicode_decimal':BODY.replace('ENSG1','ENSG١').replace('ENST1','ENST１').replace('NM_1','NM_٢'),
        'unicode_numeric':BODY.replace('ENSG1','ENSG²'),'nul_symbol':BODY.replace('GENE','G\x00ENE'),
        'unicode_symbol':BODY.replace('GENE','αβ'),'crlf':BODY.replace('\n','\r\n'),'cr':BODY.replace('\n','\r'),
        'no_final_newline':BODY.rstrip('\n'),'cr_tail':BODY+'\r\n','utf8_header':b'\xff'+BODY.encode(),
        'utf8_row':HEAD.encode()+b'\xff\n','utf8_after':BODY.encode()+b'\xff\n','utf8_partial':BODY.encode()+b'\xe2\x82',
        'bom':'\ufeff'+BODY,'quote':'"'+HEAD[:-1]+'"\n'+ROW}
    for name,payload in tables.items():add('table_'+name,'validate-table',payload=encode(payload),output=True)
    for name,suffix in [('complete','[success]\n'),('crlf','[success]\r\n'),('whitespace','[success]\n \t\r\n\v\f'),('duplicate','[success]\n[success]\n'),('trailing_row','[success]\n'+ROW),('no_stamp',''),('space',' [success]\n'),('junk','[success]\nX')]:
        add('raw_'+name,'validate-table',payload=encode(BODY+suffix),output=True,raw=True)
    for kind in ('missing','directory','symlink','broken','loop','same'):
        add('table_input_'+kind,'validate-table',payload=encode(BODY),output=True,input_kind=kind)
    for kind in ('directory','symlink','broken','loop','parent_missing'):
        add('table_output_'+kind,'validate-table',payload=encode(BODY),output=True,output_kind=kind)
    for command in ('identity','is-bundled'):
        for kind in ('normal','empty','missing','directory','symlink','broken','loop'):
            add(command+'_'+kind,command,input_kind=kind)
        for kind in ('normal','empty','missing','symlink','broken','root_alias','child_link','bad_receipt','receipt_link','receipt_broken'):
            add(command+'_bundle_'+kind,command,bundle=True,bundle_kind=kind)
        for kind in ('canonical_mapping','aliased_mapping'):
            add(command+'_helper_alias_'+kind,command,bundle=True,helper_alias=kind)
    for kind in ('normal','duplicate','reverse','missing','directory','unicode'):
        add('anchor_'+kind,'identity',anchor_kind=kind)
    for source in (GR37,'ensembl-grch38-mart113',FTP,'unknown'):
        add('source_'+source,'source',options={'source':source})
    for assembly in ('hg19','GRCh37','GRCh37.p13','hg38','GRCh38','GRCh38.p14','unknown'):
        for release in (None,'113','116','117','116.0','True'):
            add('assembly_'+assembly+'_'+str(release),'source',options=dict(assembly=assembly,**({'vep-release':release} if release is not None else {})))
    for name in (BASENAME,'ensembl_GRCh38_mart113_nm_to_enst.tsv','ensembl_GRCh38_mart116_ftp_candidate.tsv','unknown'):
        add('basename_'+name,'source',options={'basename':name})
    for command in ('receipt-check','artifact-receipt','receipt-write','artifact-route'):
        for kind in ('valid','missing','invalid','utf8','bom','null','list','extra','wrong_digest','upper_digest','bool_bytes','float_bytes','negative_bytes','missing_identity_key','extra_identity_key','anchor_missing','anchor_bad','artifact_missing','artifact_empty','artifact_changed','identity_changed','metadata_changed','duplicate_keys','wide_extra','surrogate_extra','deep_extra'):
            add(command+'_'+kind,command,receipt_kind=kind)
    # Required-success ignored fields and exact source shape/equality failures.
    for command in ('receipt-check','artifact-receipt','artifact-route'):
        for kind in ('wide_mapping_size','negative_wide_mapping_size','surrogate_mapping_digest',
                     'surrogate_provenance','wide_artifact_size','surrogate_artifact_digest',
                     'wide_anchor_size','surrogate_anchor_name','surrogate_anchor_digest',
                     'empty_anchor_name','escaped_anchor_name','wrong_anchor_type',
                     'nonfinite_extra','duplicate_ignored','duplicate_mapping_valid',
                     'surrogate_key_extra','malformed_ignored','wide_root','surrogate_root','deep_root'):
            add(command+'_'+kind,command,receipt_kind=kind)
    for kind in ('valid','existing','part_existing','part_symlink','part_broken','part_directory','receipt_symlink','receipt_directory','expected_invalid'):
        add('write_'+kind,'receipt-write',write_kind=kind)
    for kind in ('same_bad_options','changed','anchor_changed','root_alias','relative_root','missing_root','file_root','filesystem_root','outside','output_link','parent_file','target_link','target_hardlink','receipt_link','target_unknown','target_good','orphan_receipt'):
        add('route_'+kind,'artifact-route',route_kind=kind)
    registry='<r><MartURLLocation name="ENSEMBL_MART_ENSEMBL" z="extra" database="ensembl_mart_116" displayName="Ensembl Genes 116" host="grch37.ensembl.org" path="/biomart/martservice" serverVirtualSchema="default"/></r>'
    datasets='TableSet\thsapiens_gene_ensembl\tHuman genes (GRCh37.p13)\t1\tGRCh37.p13\textra\n'
    variations={'valid':registry,'empty':'','missing':'<r/>','bad':'<r>','junk':registry+'junk','duplicate':registry.replace('</r>',registry[3:-4]+'</r>'),
        'wrong_database':registry.replace('database="ensembl_mart_116"','database="other"'),'wrong_order':registry.replace('name="ENSEMBL_MART_ENSEMBL"','unexpected="value" name="ENSEMBL_MART_ENSEMBL"').replace('host="grch37.ensembl.org"','host="bad"'),
        'namespace':registry.replace('<r>','<r xmlns="urn:test">'),'entity':'<!DOCTYPE r [<!ENTITY x "value">]>'+registry.replace('z="extra"','z="&x;"'),
        'external':'<!DOCTYPE r [<!ENTITY x SYSTEM "file:///must-not-read">]>'+registry.replace('</r>','&x;</r>'),'utf16_declaration':'<?xml version="1.0" encoding="utf-16"?>'+registry}
    for name,value in variations.items():add('metadata_'+name,'parse-metadata',registry=encode(value),datasets=encode(datasets))
    for name,value in [('missing',''),('wrong',datasets.replace('GRCh37.p13','GRCh38.p14')),('duplicate',datasets*2),('short','x\thsapiens_gene_ensembl\n'),('crlf',datasets.replace('\n','\r\n')),('unicode_split','ignored\u0085'+datasets),('wrong_then_badxml','')]:
        add('datasets_'+name,'parse-metadata',registry=encode('<r>' if name=='wrong_then_badxml' else registry),datasets=encode(value))
    for name in ('valid','short','wrong_assembly','wrong_root','wrong_attribute','duplicate_attribute','wrong_main','nested_main','empty_key','malformed','namespace','utf16'):
        add('ftp_xml_'+name,'ftp-parse-xml',xml_kind=name)
    add('ftp_specification','ftp-specification')
    for name in ('normal','missing','directory','symlink'):
        add('snapshot_'+name,'snapshot-identity',input_kind=name)
    for source in (GR37,'ensembl-grch38-mart113',FTP):
        for change in ('valid','assembly','null','list','source','filename','bytes','sha256','statistics','extra_statistic','extra','flag_zero','metadata_xml','metadata_record','metadata_after','query','portable_schema','forbidden','operation'):
            add('verify_'+source+'_'+change,'verify',verify_source=source,verify_mutation=change)
    return result

def setup(folder,spec,producer,resources=None):
    folder.mkdir();source=folder/'source';(source/'scripts').mkdir(parents=True)
    for name in ('biomart_nm_to_enst.py','biomart_ftp_nm_to_enst.py'):shutil.copyfile(producer.parent/name,source/'scripts'/name)
    selected=folder/'mapping.tsv';artifact=folder/'artifact.pkl';anchor=folder/'anchor.tsv';second=folder/'second.tsv';output=folder/'output.tsv'
    selected.write_bytes(bytes.fromhex(spec.get('payload',encode(BODY))));artifact.write_bytes(b'original artifact\n');anchor.write_bytes(b'anchor\n');second.write_bytes(b'second\n')
    command=spec['command'];options=dict(spec.get('options',{}));inputs=[]
    kind=spec.get('input_kind','normal')
    def replace_kind(path,kind):
        if kind=='normal':return
        if path.exists() or path.is_symlink():path.unlink()
        if kind=='empty':path.touch()
        elif kind=='directory':path.mkdir()
        elif kind in ('symlink','broken','loop'):
            target=path.with_name(path.name+'.target')
            if kind=='symlink':target.write_bytes(BODY.encode())
            path.symlink_to(path.name if kind=='loop' else target.name)
    if kind!='same':replace_kind(selected,kind)
    if spec.get('bundle'):
        selected=source/'data/hgnc'/BASENAME;selected.parent.mkdir(parents=True);selected.write_bytes(b'arbitrary bundled bytes\n')
        bk=spec.get('bundle_kind');replace_kind(selected,bk)
        side=Path(str(selected)+'.provenance.json')
        if bk=='bad_receipt':side.write_text('not JSON\n')
        if bk in ('receipt_link','receipt_broken'):
            target=folder/'bundle-receipt-target.json'
            if bk=='receipt_link':target.write_text('ignored provenance content\n')
            side.symlink_to(target)
        if bk=='root_alias':
            alias=folder/'source-alias';alias.symlink_to(source,target_is_directory=True);selected=alias/'data/hgnc'/BASENAME
        if bk=='child_link':
            origin=source/'data/hgnc';destination=folder/'outside-hgnc';origin.rename(destination);origin.symlink_to(destination,target_is_directory=True)
    if 'anchor_kind' in spec:
        ak=spec['anchor_kind'];inputs=[anchor,second]
        if ak=='duplicate':inputs=[anchor,anchor,second]
        if ak=='reverse':inputs.reverse()
        if ak=='missing':anchor.unlink()
        if ak=='directory':anchor.unlink();anchor.mkdir()
        if ak=='unicode':target=folder/'α.tsv';anchor.rename(target);inputs=[target,second]
    receipt=Path(str(artifact)+'.nm_to_enst.json');base_identity=dict(mapping_sha256=hashlib.sha256(selected.read_bytes()).hexdigest(),mapping_bytes=selected.stat().st_size,provenance_sha256=None,anchor_inputs=[]) if command.startswith('receipt') or command in ('artifact-receipt','artifact-route') else None
    if base_identity is not None:
        record=dict(mapping_identity=base_identity,artifact_bytes=artifact.stat().st_size,artifact_sha256=hashlib.sha256(artifact.read_bytes()).hexdigest());rk=spec.get('receipt_kind','valid')
        if rk=='extra':record['ignored']={'x':[1,True,None]}
        if rk=='wrong_digest':record['mapping_identity']['mapping_sha256']='0'*64
        if rk=='upper_digest':record['mapping_identity']['mapping_sha256']=record['mapping_identity']['mapping_sha256'].upper()
        if rk=='bool_bytes':record['mapping_identity']['mapping_bytes']=True
        if rk=='float_bytes':record['mapping_identity']['mapping_bytes']=float(record['mapping_identity']['mapping_bytes'])
        if rk=='negative_bytes':record['mapping_identity']['mapping_bytes']=-1
        if rk=='missing_identity_key':record['mapping_identity'].pop('provenance_sha256')
        if rk=='extra_identity_key':record['mapping_identity']['extra']=None
        if rk=='anchor_missing':record['mapping_identity']['anchor_inputs']=[{}]
        if rk=='anchor_bad':record['mapping_identity']['anchor_inputs']=[dict(name='',bytes=1,sha256='0'*64)]
        if rk=='artifact_missing':artifact.unlink()
        if rk=='artifact_empty':artifact.write_bytes(b'')
        if rk=='artifact_changed':artifact.write_bytes(b'changed artifact')
        if rk=='identity_changed':selected.write_bytes(BODY.replace('ENST1','ENST2').encode())
        if rk=='metadata_changed':Path(str(selected)+'.provenance.json').write_text('{}')
        if rk=='wide_extra':record['ignored']=2**160
        if rk=='surrogate_extra':record['ignored']='\ud800'
        if rk in ('wide_mapping_size','negative_wide_mapping_size'):record['mapping_identity']['mapping_bytes']=(-1 if rk.startswith('negative') else 1)*2**160
        if rk=='surrogate_mapping_digest':record['mapping_identity']['mapping_sha256']='\ud800'
        if rk=='surrogate_provenance':record['mapping_identity']['provenance_sha256']='\ud800'
        if rk=='wide_artifact_size':record['artifact_bytes']=2**160
        if rk=='surrogate_artifact_digest':record['artifact_sha256']='\ud800'
        if rk in ('wide_anchor_size','surrogate_anchor_name','surrogate_anchor_digest','empty_anchor_name','escaped_anchor_name','wrong_anchor_type'):
            item=dict(name='anchor.tsv',bytes=7,sha256='0'*64)
            if rk=='wide_anchor_size':item['bytes']=2**160
            if rk=='surrogate_anchor_name':item['name']='\ud800'
            if rk=='surrogate_anchor_digest':item['sha256']='\ud800'
            if rk=='empty_anchor_name':item['name']=''
            if rk=='escaped_anchor_name':item['name']='a,]\\\"b'
            record['mapping_identity']['anchor_inputs']=[None if rk=='wrong_anchor_type' else item]
        if rk=='nonfinite_extra':record['ignored']=[float('nan'),float('inf'),-float('inf'),-0.0]
        if rk=='surrogate_key_extra':record['\ud800']={'\ud801':2**160}
        raw=json.dumps(record)
        if rk=='invalid':raw='{bad'
        if rk=='utf8':raw=b'\xff'
        if rk=='bom':raw='\ufeff'+raw
        if rk=='null':raw='null'
        if rk=='list':raw='[]'
        if rk=='duplicate_keys':raw=raw[:-1]+',"artifact_bytes":1}'
        if rk=='deep_extra':raw=raw[:-1]+',"ignored":'+'['*600+'0'+']'*600+'}'
        if rk=='duplicate_ignored':raw=raw[:-1]+',"ignored":1,"ignored":{"b":2,"a":3,"b":4}}'
        if rk=='duplicate_mapping_valid':raw='{"mapping_identity":'+str(2**160)+','+raw[1:]
        if rk=='malformed_ignored':raw=raw[:-1]+',"ignored":[1,]}'
        if rk=='wide_root':raw=str(2**160)
        if rk=='surrogate_root':raw='"\\ud800"'
        if rk=='deep_root':raw='['*600+'0'+']'*600
        if rk!='missing':receipt.write_bytes(raw if isinstance(raw,bytes) else raw.encode())
        if command=='receipt-write':options['expected-identity']=json.dumps(base_identity)
    if command in ('identity','is-bundled','verify','snapshot-identity','ftp-table-identity'):options['mapping']=str(selected)
    elif command in ('receipt-check','receipt-write'):options.update(artifact=str(artifact),mapping=str(selected))
    elif command=='artifact-receipt':options['artifact']=str(artifact)
    elif command=='validate-table':
        options['table']=str(selected)
        if spec.get('output'):
            output.write_bytes(b'previous output\n');ok=spec.get('output_kind')
            if ok:replace_kind(output,ok)
            if ok=='parent_missing':output=folder/'missing-parent/output.tsv'
            options['output']=str(selected if kind=='same' else output)
    elif command=='parse-metadata':
        for key in ('registry','datasets'):
            path=folder/(key+'.txt');path.write_bytes(bytes.fromhex(spec[key]));options[key]=str(path)
    elif command=='artifact-route':
        cache=folder/'cache';cache.mkdir();dest=cache/'plugins/AlphaMissense';route=spec.get('route_kind','same')
        if route not in ('same','same_bad_options'):selected.write_bytes(BODY.replace('ENST1','ENST9').encode())
        options.update(artifact=str(artifact),mapping=str(selected),**{'cache-root':str(cache),'output-dir':str(dest),'assembly':'hg38','kind':'motifs'})
        if route=='same_bad_options':options.update(assembly='bad',kind='bad',**{'cache-root':'relative','output-dir':'relative'})
        if route=='anchor_changed':inputs=[anchor]
        if route=='root_alias':alias=folder/'alias';alias.symlink_to(cache,target_is_directory=True);options.update(**{'cache-root':str(alias),'output-dir':str(alias/'plugins/AlphaMissense')})
        if route=='relative_root':options['cache-root']='relative'
        if route=='missing_root':options['cache-root']=str(folder/'missing')
        if route=='file_root':options['cache-root']=str(artifact)
        if route=='filesystem_root':options['cache-root']='/'
        if route=='outside':options['output-dir']=str(folder/'outside')
        if route=='output_link':dest.parent.mkdir();dest.symlink_to(folder,target_is_directory=True)
        if route=='parent_file':(cache/'plugins').write_bytes(b'file')
        if route.startswith(('target_','receipt_','orphan_')):
            b=module(source/'scripts/biomart_nm_to_enst.py');identity=b.mapping_identity(selected,inputs);target=b.changed_mapping_target(identity=identity,cache_root=cache,output_dir=dest,assembly='hg38',kind='motifs');target.parent.mkdir(parents=True);side=Path(str(target)+'.nm_to_enst.json')
            if route=='target_link':target.symlink_to(artifact)
            elif route=='target_hardlink':os.link(artifact,target)
            elif route=='receipt_link':side.symlink_to(receipt)
            elif route in ('target_unknown','target_good','orphan_receipt'):
                target.write_bytes(b'new artifact');side.write_text(json.dumps(dict(mapping_identity=identity,artifact_bytes=target.stat().st_size,artifact_sha256=hashlib.sha256(target.read_bytes()).hexdigest())))
                if route=='target_unknown':side.unlink()
                if route=='orphan_receipt':target.unlink()
    elif command=='ftp-parse-xml':
        import xml.etree.ElementTree as ET
        b=module(source/'scripts/biomart_nm_to_enst.py');f=b.FTP
        row=['95',f.DATASET,'Human genes (GRCh38.p14)','Ensembl genes','TableSet','1',f.ASSEMBLY,'2026-04-30']
        root=ET.Element('DatasetConfig',dataset=f.DATASET,datasetID='95',version=f.ASSEMBLY,displayName=row[2],interfaces='default')
        for name,(key,table,field,display) in f.ATTRIBUTES.items():ET.SubElement(root,'AttributeDescription',internalName=name,key=key,tableConstraint=table,field=field,displayName=display)
        for level,key in zip(('gene','transcript','translation'),('gene_id_1020_key','transcript_id_1064_key','translation_id_1068_key')):
            ET.SubElement(root,'MainTable').text=f.DATASET+'__'+level+'__main';ET.SubElement(root,'Key').text=key
        xk=spec['xml_kind']
        if xk=='short':row=row[:3]
        if xk=='wrong_assembly':row[6]='wrong'
        if xk=='wrong_root':root.tag='Wrong'
        if xk=='wrong_attribute':root.find('AttributeDescription').set('key','other')
        if xk=='duplicate_attribute':ET.SubElement(root,'AttributeDescription',**root.find('AttributeDescription').attrib)
        if xk=='wrong_main':root.find('MainTable').text='other'
        if xk=='nested_main':node=root.find('MainTable');root.remove(node);ET.SubElement(root,'container').append(node)
        if xk=='empty_key':root.find('Key').text=None
        if xk=='namespace':root.set('xmlns','urn:test')
        raw=ET.tostring(root,encoding='utf-16' if xk=='utf16' else 'utf-8')
        if xk=='malformed':raw+=b'junk'
        path=folder/'dataset.xml';path.write_bytes(raw);options.update(xml=str(path),**{'dataset-row':json.dumps(row)})
    wk=spec.get('write_kind')
    if wk:
        part=Path(str(receipt)+'.part')
        if wk=='part_existing':part.write_text('retained staging')
        if wk in ('part_symlink','part_broken'):part.symlink_to(artifact if wk=='part_symlink' else folder/'absent')
        if wk=='part_directory':part.mkdir()
        if wk=='receipt_symlink':receipt.unlink();receipt.symlink_to(folder/'replaced-link-target')
        if wk=='receipt_directory':receipt.unlink();receipt.mkdir()
        if wk=='expected_invalid':options['expected-identity']='{bad'
    if 'verify_source' in spec:
        b=module(source/'scripts/biomart_nm_to_enst.py');selected_source=b._source(spec['verify_source']);reference=(resources or producer.parent.parent/'data/hgnc')/selected_source['basename'];record=json.loads(Path(str(reference)+'.provenance.json').read_text());change=spec['verify_mutation']
        selected.unlink();selected.symlink_to(reference)
        if change=='assembly':options['assembly']='wrong'
        if change=='source':record['source_id']='unknown'
        if change=='filename':selected=folder/(b.GRCH38_BASENAME if selected_source['basename']==b.BASENAME else b.BASENAME);selected.symlink_to(reference);options['mapping']=str(selected)
        if change=='bytes':record['bytes']+=1
        if change=='sha256':record['sha256']='0'*64
        if change=='statistics':record['statistics']['rows']+=1
        if change=='extra_statistic':record['statistics']['extra']=0
        if change=='extra':record['ignored']={'unused':[1,True,None]}
        if change=='flag_zero':record['scientific_qualified']=0;record['completion_stamp_verified']=0
        if change=='metadata_xml':record.setdefault('metadata_before',{})['registry_xml']='<broken>'
        if change=='metadata_record':record.setdefault('metadata_before',{})['registry_record']={}
        if change=='metadata_after':record.setdefault('metadata_after',{})['dataset_record']=[]
        if change=='query':record['query_sha256']='0'*64
        if change=='portable_schema':record['metadata_schema_version']=True
        if change=='forbidden':record['query_xml']='forbidden'
        if change=='operation':record.pop('metadata_kind',None);record['operation']='unknown'
        if change=='null':record=None
        if change=='list':record=[]
        Path(str(selected)+'.provenance.json').write_text(json.dumps(record))
    # Fixed setup timestamps make preservation observations independently exact.
    for path in folder.rglob('*'):
        if path.is_file() and not path.is_symlink():os.utime(path,ns=(1700000000000000000,1700000000000000000))
    helper_root=source
    if spec.get('helper_alias'):
        helper_root=folder/'helper-alias';helper_root.symlink_to(source,target_is_directory=True)
        if spec['helper_alias']=='aliased_mapping':options['mapping']=str(helper_root/'data/hgnc'/BASENAME)
    return dict(command=command,options=options,inputs=list(map(str,inputs)),raw=spec.get('raw',False),helper_root=str(helper_root))

def inventory(folder):
    result={}
    for path in sorted(folder.rglob('*')):
        if path.name in ('request.json','result.json','stdout','stderr','arm.json'):continue
        s=path.lstat();record=dict(mode=stat.S_IMODE(s.st_mode),kind=stat.S_IFMT(s.st_mode),inode=s.st_ino,device=s.st_dev,nlink=s.st_nlink,bytes=s.st_size,mtime_ns=s.st_mtime_ns,ctime_ns=s.st_ctime_ns)
        if path.is_symlink():record['link']=os.readlink(path)
        elif path.is_file():record['sha256']=pin(path)['sha256']
        result[str(path.relative_to(folder))]=record
    return result
def comparable(before,after,folder):
    result={}
    for name,item in after.items():
        if item['kind']==stat.S_IFDIR:continue
        result[name]={k:v for k,v in item.items() if k not in ('inode','device','mtime_ns','ctime_ns')}
        if 'link' in item:
            result[name]['link']=item['link'].replace(str(folder),'$ARM')
            result[name]['bytes']=len(result[name]['link'].encode())
        if name in before:
            result[name]['unchanged_identity']=before[name]==item
    return result
def normalize(value,folder):
    return json.loads(json.dumps(value,default=str).replace(str(folder),'$ARM'))
def exact_value(left,right):
    # Preserve JSON types and IEEE float values; NaN is a source-supported JSON
    # constant. Iteration avoids imposing the harness's recursion limit.
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
        elif isinstance(a,float):
            if not (math.isnan(a) and math.isnan(b)) and struct.pack('!d',a)!=struct.pack('!d',b):return False
        elif a!=b:return False
    return True

def key_order(value):
    result=[];pending=[((),value)]
    while pending:
        path,item=pending.pop()
        if isinstance(item,dict):
            result.append((path,list(item)))
            pending.extend((path+(key,),child) for key,child in reversed(list(item.items())))
        elif isinstance(item,list):pending.extend((path+(i,),child) for i,child in reversed(list(enumerate(item))))
    return result

def readback(arm):
    value=json.loads(arm.read_text());observed=inventory(arm.parent);saved=value['after']
    assert set(observed)==set(saved), 'Saved oracle entries changed'
    drift={}
    for key,record in saved.items():
        now=observed[key]
        # Device IDs and directory byte counts are node-local observations.
        # Times are retained verbatim, never rounded or substituted. Oracle
        # reuse binds original payload bytes, publication mode and link shape;
        # it does not assert a cross-node stat timestamp is unchanged.
        fields=('kind','mode') if record['kind']==stat.S_IFDIR else ('kind','mode','nlink','bytes','sha256','link')
        assert all(now.get(k)==record.get(k) for k in fields), 'Saved oracle bytes/mode/link changed: '+key
        if now!=record:drift[key]={k:[record.get(k),now.get(k)] for k in set(record)|set(now) if record.get(k)!=now.get(k)}
    return value,dict(arm=pin(arm),observed=observed,raw_metadata_differences=drift,original_bytes_modes_links_exact=True)
def run_arm(folder,spec,producer,binary,native,resources=None):
    request=setup(folder,spec,producer,resources);save(folder/'request.json',request);before=inventory(folder)
    if native:
        command=[str(binary),'biomart',request['command'],'--helper-root',request['helper_root']]
        for key,value in request['options'].items():command+=['--'+key,str(value)]
        for path in request['inputs']:command+=['--input',path]
        if request['raw']:command+=['--raw']
    else:command=[sys.executable,str(Path(__file__).resolve()),'oracle','--producer',request['helper_root']+'/scripts/biomart_nm_to_enst.py','--request',str(folder/'request.json'),'--result',str(folder/'result.json')]
    process=subprocess.run(command,cwd=folder,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=30,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
    (folder/'stdout').write_bytes(process.stdout);(folder/'stderr').write_bytes(process.stderr)
    if native:
        if process.returncode==0:result=dict(success=True,value=json.loads(process.stdout) if process.stdout else None)
        else:
            message=process.stderr.decode(errors='replace').rstrip('\n');match=re.fullmatch(r'priva: (\w+): ([\s\S]*)',message)
            result=dict(success=False,category=match[1] if match else 'NATIVE_BOUNDARY',message=match[2] if match else message)
        save(folder/'result.json',result)
    else:result=json.loads((folder/'result.json').read_text())
    after=inventory(folder);value=dict(command=command,returncode=process.returncode,result=result,before=before,after=after,comparable=comparable(before,after,folder),stdout_sha256=hashlib.sha256(process.stdout).hexdigest(),stderr_sha256=hashlib.sha256(process.stderr).hexdigest())
    save(folder/'arm.json',value);return value
def check(args):
    args.work_root.mkdir();producer_pin=pin(args.producer);ftp_pin=pin(args.producer.with_name('biomart_ftp_nm_to_enst.py'));binary_pin=pin(args.binary) if args.binary else None
    assert args.binary or args.oracle_only
    fixtures=cases();results=[];reuse=None;reuse_cases={};readbacks={}
    rerun_native=set(getattr(args,'rerun_native_case',[]) or [])
    if rerun_native-set(fixtures):raise ValueError('Unknown native replay case')
    if rerun_native and not args.reuse_native_root:raise ValueError('Native case replay requires --reuse-native-root')
    resources=args.resources or args.producer.parent.parent/'data/hgnc'
    resource_pins={str(p):pin(p) for p in sorted(resources.glob('ensembl*.tsv*')) if p.is_file()}
    if args.reuse_native_root:
        assert args.reuse_native_harness
        spec=importlib.util.spec_from_file_location('prior_biomart_cases',args.reuse_native_harness);old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
        assert old.cases()==fixtures
    if args.reuse_report:
        reuse=json.loads(args.reuse_report.read_text())
        assert reuse['producer']['sha256']==producer_pin['sha256'] and reuse['ftp_producer']['sha256']==ftp_pin['sha256']
        reuse_cases={v['name']:v for v in reuse['cases']}
    for name,spec in fixtures.items():
        root=args.work_root/name;root.mkdir()
        if name in reuse_cases and reuse['fixtures'][name]==spec:
            arm=Path(reuse_cases[name]['python']);python,proof=readback(arm);python_root=arm.parent
            readbacks[name]={'python':proof}
        else:
            python_root=root/'python';python=run_arm(python_root,spec,args.producer,args.binary,False,resources)
        if args.oracle_only:
            record=dict(name=name,exact=False,differences=['native_unmeasured'],python=str(python_root/'arm.json'),native=None,python_reused=python_root!=root/'python')
            results.append(record);print(json.dumps(record),flush=True);continue
        saved_native=args.reuse_native_root/name/'native/arm.json' if args.reuse_native_root and name not in rerun_native else None
        if saved_native and saved_native.is_file():
            native,proof=readback(saved_native);native_root=saved_native.parent
            assert pin(native['command'][0])['sha256']==binary_pin['sha256']
            readbacks.setdefault(name,{})['native']=proof
        else:native_root=root/'native';native=run_arm(native_root,spec,args.producer,args.binary,True,resources)
        # Re-evaluate saved observed files under the same explicit per-arm path
        # binding. Original raw values remain in arm.json and reuse readbacks.
        python['comparable']=comparable(python['before'],python['after'],python_root)
        native['comparable']=comparable(native['before'],native['after'],native_root)
        reasons=[]
        for field in ('returncode','result','comparable'):
            if not exact_value(normalize(python[field],python_root),normalize(native[field],native_root)):reasons.append(field)
        order=None
        if spec['command']=='artifact-receipt' and python['result']['success'] and native['result']['success']:
            # The original returns this parsed object unchanged. Saved oracle
            # result.json used sorted keys, so recover producing order from its
            # source-bound, read-back receipt rather than storage metadata.
            py_request=json.loads((python_root/'request.json').read_text())
            original=json.loads(Path(py_request['options']['artifact']+'.nm_to_enst.json').read_text())
            assert exact_value(original,python['result']['value'])
            observed=json.loads((native_root/'stdout').read_bytes())
            order=key_order(original)==key_order(observed)
            if not order:reasons.append('return_dictionary_order')
        record=dict(name=name,exact=not reasons,differences=reasons,python=str(python_root/'arm.json'),native=str(native_root/'arm.json'),python_reused=python_root!=root/'python',native_reused=native_root!=root/'native',return_dictionary_order_exact=order)
        results.append(record);print(json.dumps(record),flush=True)
    assert pin(args.producer)==producer_pin and (not args.binary or pin(args.binary)==binary_pin) and pin(args.producer.with_name('biomart_ftp_nm_to_enst.py'))==ftp_pin
    assert all(pin(path)==identity for path,identity in resource_pins.items())
    report=dict(full_return_contract='artifact-receipt preserves complete JSON values and dictionary insertion order; lexical JSON encoding is not claimed byte-identical',producer=producer_pin,ftp_producer=ftp_pin,binary=binary_pin,harness=pin(__file__),fixtures=fixtures,cases=results,total=len(results),exact=sum(v['exact'] for v in results),mismatches=[v for v in results if not v['exact']],source_stable=True,reuse_report=pin(args.reuse_report) if args.reuse_report else None,reuse_readbacks=readbacks,resource_pins=resource_pins,resources_preserved=True,reuse_native_harness=pin(args.reuse_native_harness) if args.reuse_native_harness else None)
    report['rerun_native_cases']=sorted(rerun_native)
    save(args.report,report)
def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(required=True)
    q=sub.add_parser('oracle');q.add_argument('--producer',type=Path,required=True);q.add_argument('--request',type=Path,required=True);q.add_argument('--result',type=Path,required=True);q.set_defaults(run=oracle)
    q=sub.add_parser('check');q.add_argument('--producer',type=Path,required=True);q.add_argument('--binary',type=Path);q.add_argument('--work-root',type=Path,required=True);q.add_argument('--report',type=Path,required=True);q.add_argument('--reuse-report',type=Path);q.add_argument('--reuse-native-root',type=Path);q.add_argument('--reuse-native-harness',type=Path);q.add_argument('--resources',type=Path);q.add_argument('--oracle-only',action='store_true');q.set_defaults(run=check)
    q.add_argument('--rerun-native-case',action='append',default=[],help='Replay this named native case while reusing other source-bound native outcomes')
    args=parser.parse_args();args.run(args)
if __name__=='__main__':main()
