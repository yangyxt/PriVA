#!/usr/bin/env python3
"""Compare optional prepared condition-evidence procedures with unchanged Python.

check --producer SCRIPT --binary PRIVA --work-root NEW_DIRECTORY
oracle --producer SCRIPT --request JSON
All fixture inputs, measured commands, exact outputs, diagnostics and modes are
retained. Runtime category/message comparison is separate from raw stderr bytes.
This harness does not download sources or call the installer.
"""
from __future__ import annotations
import argparse
import bz2
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import lzma
import os
from pathlib import Path
import resource
import stat
import subprocess
import sys
import time
import types
from xml.sax.saxutils import escape
sys.dont_write_bytecode=True
G2P_URL='https://www.ebi.ac.uk/gene2phenotype/api/panel/all/download'
ORPHA_URL='https://www.orphadata.com/data/xml/en_product6.xml'
G2P_HEADER=['gene symbol','disease name','molecular mechanism','confidence','allelic requirement','cross cutting modifier','phenotypes','pmids','g2p id','disease mim','disease MONDO']
REGISTRY_HEADER=['disease_id','mondo_id','disease_scope','priva_scope','scope_review_status']

def pin(path):
    path=Path(path)
    with path.open('rb') as handle:digest=hashlib.file_digest(handle,'sha256').hexdigest()
    return dict(path=str(path),bytes=path.stat().st_size,sha256=digest)

def save(path,value):
    with Path(path).open('x') as handle:json.dump(value,handle,indent=2,sort_keys=True);handle.write('\n')

def encode(names,rows,delimiter=','):
    stream=io.StringIO(newline='');writer=csv.writer(stream,delimiter=delimiter,lineterminator='\n');writer.writerow(names);writer.writerows(rows);return stream.getvalue().encode()

def g2p(**values):
    base=dict(zip(G2P_HEADER,['GENE','Condition','loss of function','definitive','monoallelic','complete penetrance','HP:0034950','123','G2P:1','1','MONDO:1']))
    base.update(values);return encode(G2P_HEADER,[[base[k] for k in G2P_HEADER]])

def association(gene='GENE',status='Assessed',relationship='Disease-causing germline mutation(s) (loss of function)',validation='123[PMID] 456 [PMID] 123[PMID]'):
    return '<DisorderGeneAssociation><DisorderGeneAssociationStatus><Name>'+escape(status)+'</Name></DisorderGeneAssociationStatus><DisorderGeneAssociationType><Name>'+escape(relationship)+'</Name></DisorderGeneAssociationType><Gene><Symbol>'+escape(gene)+'</Symbol></Gene><SourceOfValidation>'+escape(validation)+'</SourceOfValidation></DisorderGeneAssociation>'

def disorder(associations=None,code='1',disease='Condition'):
    return '<Disorder><OrphaCode>'+escape(code)+'</OrphaCode><Name>'+escape(disease)+'</Name><DisorderGeneAssociationList>'+(association() if associations is None else associations)+'</DisorderGeneAssociationList></Disorder>'

def xml(content=None):return ('<JDBOR><DisorderList>'+(disorder() if content is None else content)+'</DisorderList></JDBOR>').encode()

def cases():
    result={}
    def add(name,command='parse',kind='g2p',payload=None,**options):
        result[name]=dict(command=command,kind=kind,payload=(g2p() if payload is None else payload).hex(),**options)
    add('g2p_normal')
    normal=g2p()
    for label,raw in [
        ('empty',b''),('blank',b' \n\r\n'),('header_only',normal.split(b'\n')[0]+b'\n'),
        ('wrong_header',b'other,fields\na,b\n'),('bom',b'\xef\xbb\xbf'+normal),('double_bom',b'\xef\xbb\xbf'*2+normal),
        ('crlf',normal.replace(b'\n',b'\r\n')),('cr',normal.replace(b'\n',b'\r')),('unterminated_final',normal.rstrip(b'\n')),
        ('quoted_multiline',g2p(**{'disease name':'a\nb\rc\r\nd'})),('unclosed',b'gene symbol,disease name\nG,"bad'),
        ('unclosed_after_blank',b'gene symbol,disease name\n\n\nG,"bad'),
        ('unclosed_embedded_line',b'gene symbol,disease name\nG,"bad\nmore\n'),
        ('quote_follower',b'gene symbol,disease name\nG,"a"b\n'),('nul_header',b'gene symbol\0more,disease name\nG,D\n'),
        ('nul_value',g2p(**{'disease name':'alpha\0beta'})),('no_gene_header',b'gene,disease name\nG,D\n'),
        ('space_line',b'gene symbol,disease name\n   \nG,D\n'),('tab_line',b'gene symbol,disease name\n\t\nG,D\n'),
        ('long_first',b'gene symbol,disease name\nindex,G,D\nindex,H,E\n'),
        ('multi_index',b'gene symbol,disease name\nindex,level,G,D\ni,l,H,E\n'),
        ('long_later',b'gene symbol,disease name\nG,D\nindex,H,E\nJ,F\n'),
        ('longer_later',b'gene symbol,disease name\ni,G,D\nx,y,H,E\nz,J,F\n'),
        ('short_later',b'gene symbol,disease name\ni,G,D\nJ\n'),
        ('short_first',b'gene symbol,disease name,molecular mechanism\nG\nH,D,LOF\n'),
        ('duplicate_header',b'gene symbol,gene symbol,disease name\nFIRST,SECOND,D\n'),
        ('normalized_duplicate',b'gene symbol,gene_symbol,disease name\nFIRST,SECOND,D\n'),
        ('candidate_precedence',b'gene symbol,disease label,disease name\nG,LABEL,NAME\n'),
        ('unicode_header',g2p().replace(b'gene symbol','GENE SYMBOL'.encode())),
        ('wide_utf8',g2p(**{'disease name':'α'*140000})),
        ('bad_utf8',normal+b'\xff'),('bad_utf8_late',normal+b'\n'*262144+b'\xff'),
        ('bad_bom',b'\xef\xbb\xbf'+normal+b'\xff'),('truncated_utf8',normal+b'\xe2\x82'),
    ]:add('csv_'+label,payload=raw)
    for value in ['', 'NA','NaN','nan','N/A','NULL','None','null','<NA>','#N/A','#N/A N/A','#NA','-1.#IND','-1.#QNAN','-NaN','-nan','1.#IND','1.#QNAN','n/a','NAN','0007','0','True',' nan ','\u001cGENE\u001f']:
        add('missing_'+str(len(result)),payload=g2p(**{'gene symbol':value,'disease name':value,'molecular mechanism':value}))
    for value in ['dominant negative','dominant\u001cnegative','xdominant negative','dominant negativex','dominant_negative','dominant-negative',
        'gain of function','gain function','activating','activated','activation','increased activity','LOF','absent gene product','absence of gene product',
        'reduced function','decreased function','partial loss','haploinsufficiency','haplo insufficient','triplosensitivity','triplo sensitivity',
        'increased gene dosage','decreased dosage','altered gene product','altered protein','all missense','inframe','in-frame','unmapped','\u2003',
        'loss of function; gain of function; dominant negative; partial loss']:
        add('mechanism_'+str(len(result)),payload=g2p(**{'molecular mechanism':value}))
    for value in ['1','2','3',' 3 ','03','definitive','sufficient','green','strong','moderate','emerging','limited','amber','supportive','refuted','disputed','red','unlikely','unknown','not sufficient, disputed','red strong','strong moderate','limited red']:
        add('confidence_'+str(len(result)),payload=g2p(**{'confidence':value}))
    for raw,terms in [('complete penetrance','HP:0003829'),('incomplete penetrance','HP:0034950'),('age dependent penetrance',''),
        ('variable expressivity',''),('late onset',''),('high penetrance','HP:4000160;HP:0034950;HP:4000160'),('', 'hp:0034950,HP:4000158'),('', 'HP:0034950|HP:0003829')]:
        add('penetrance_'+str(len(result)),payload=g2p(**{'cross cutting modifier':raw,'phenotypes':terms}))
    for suffix in ('.gz','.bz2','.xz','.zst','.zip','.tar','.tar.gz'):
        add('compression_'+suffix.replace('.','_'),suffix=suffix,compress=suffix)
    for suffix in ('.gz','.bz2','.xz','.zst'):
        for label,payload in [('invalid',normal+b'\xff'),('late',normal+b'\n'*262144+b'\xff'),('terminal',normal+b'\xe2\x82')]:
            add('compressed_utf8_'+suffix[1:]+'_'+label,suffix=suffix,compress=suffix,payload=payload)
    for label,payload in [
        ('one_before',b'A'*(262144-1)+b'\xff'),('at_boundary',b'A'*262144+b'\xff'),
        ('bom_boundary',b'\xef\xbb\xbf'+b'A'*(262144-3)+b'\xff'),
        ('multibyte_boundary',b'A'*(262144-1)+'α'.encode()+b'\xff'),
        ('unicode_read_ratio',('α'*140000).encode()+b'\xff'),
        ('emoji_read_ratio',('😀'*70000).encode()+b'\xff'),
        ('mixed_read_ratio',('α'*140000+'A'*150000).encode()+b'\xff'),
    ]:add('decode_'+label,payload=payload)
    for suffix in ('.zip','.tar','.tar.gz'):
        for layout in ('empty','multiple','directory','duplicate','unicode'):
            add('archive_'+suffix.replace('.','_')+'_'+layout,suffix=suffix,compress=suffix,archive_layout=layout)
    add('orpha_normal',kind='orphadata',payload=xml())
    for label,content in [
        ('no_disorder',''),('missing_code',disorder(code='')),('empty_gene',disorder(association(gene=' '))),
        ('duplicate',disorder()+disorder()),('case_status',disorder(association(status='assessed'))),
        ('unassessed',disorder(association(status='Not yet assessed'))),('somatic',disorder(association(relationship='Disease-causing somatic mutation(s) (loss of function)'))),
        ('no_mechanism',disorder(association(relationship='Disease-causing germline mutation(s)'))),
        ('gof',disorder(association(relationship='Disease-causing germline mutation(s) (gain of function)'))),
        ('both',disorder(association(relationship='Disease-causing germline mutation(s) (gain of function and loss of function)'))),
        ('activating',disorder(association(relationship='Disease-causing germline activating'))),
        ('haplo',disorder(association(relationship='Disease-causing germline haploinsufficiency'))),
        ('unicode_pmids',disorder(association(validation='١٢٣\u2003[PMID] 123[PMID] ١٢٣[PMID] ²[PMID]'))),
        ('first_code','<Disorder><OrphaCode>1</OrphaCode><OrphaCode>2</OrphaCode><DisorderGeneAssociationList>'+association()+'</DisorderGeneAssociationList></Disorder>'),
        ('nested_disorder',disorder(disorder()+association())),
        ('namespaced_disorder',disorder().replace('<Disorder>','<Disorder xmlns="urn:test">')),
        ('namespaced_child',disorder().replace('<Name>Condition','<Name xmlns="urn:test">Condition')),
        ('nested_text',disorder(disease='START').replace('START','a<b>b</b>c')),
        ('nested_association','<Disorder><OrphaCode>1</OrphaCode><other><DisorderGeneAssociationList>'+association()+'</DisorderGeneAssociationList></other></Disorder>'),
        ('comments_pi',disorder().replace('Condition','a<!--comment--><?pi text?>b')),
    ]:add('xml_'+label,kind='orphadata',payload=xml(content))
    for label,raw in [
        ('wrong_root',xml().replace(b'JDBOR',b'root')),('truncated',xml()[:-10]),('empty',b''),('junk_after',xml()+b'junk'),
        ('external',b'<!DOCTYPE JDBOR [<!ENTITY e SYSTEM "file:///not-read">]>'+xml().replace(b'Condition',b'&e;')),
        ('internal',b'<!DOCTYPE JDBOR [<!ENTITY e "Text">]>'+xml().replace(b'Condition',b'&e;')),
        ('bad_entity',xml().replace(b'Condition',b'&absent;')),('invalid_byte',xml()+b'\xff'),
        ('after_good_block',xml(disorder()+(' '*17000)+'<bad>').replace(b'</DisorderList>',b'')),
        ('nested_clear_boundary',xml(disorder(associations=disorder()+(' '*16384)+association()))),
    ]:add('xml_'+label,kind='orphadata',payload=raw)
    for encoding in ('utf-8','utf-16','utf-16-le','utf-16-be','iso-8859-1','windows-1252','utf-32','unknown'):
        raw=('<?xml version="1.0" encoding="'+encoding+'"?>'+xml().decode()).encode('utf-8' if encoding=='unknown' else encoding)
        add('xml_encoding_'+encoding,kind='orphadata',payload=raw)
    registry=encode(REGISTRY_HEADER,[['OMIM:1','MONDO:2','germline','include','reviewed'],['ORPHA:1','MONDO:3','germline','include','reviewed']],delimiter='\t')
    for label,options in [
        ('both',{}),('missing_both',dict(missing_sources=True)),('default_bad_g2p',dict(g2p_payload=b'gene symbol\n"bad'.hex())),
        ('custom_bad_g2p',dict(g2p_payload=b'gene symbol\n"bad'.hex(),g2p_url='https://selected.invalid/source')),
        ('default_bad_xml',dict(xml_payload=b'<bad>'.hex())),('custom_bad_xml',dict(xml_payload=b'<bad>'.hex(),orpha_url='https://selected.invalid/source')),
        ('both_default_bad',dict(g2p_payload=b'"bad'.hex(),xml_payload=b'<bad>'.hex())),
        ('custom_good',dict(g2p_url='https://selected.invalid/g2p',orpha_url='https://selected.invalid/xml')),
        ('no_registry',dict(no_registry=True)),('missing_registry',dict(registry_kind='missing')),
        ('empty_registry',dict(registry='')),('bad_registry',dict(registry=b'disease_id\nOMIM:1\n'.hex())),
        ('registry_directory',dict(registry_kind='directory')),('registry_denied',dict(registry_kind='denied')),
        ('fill_missing_mondo',dict(g2p_payload=g2p(**{'disease MONDO':''}).hex())),
        ('registry_first_duplicate',dict(registry=(registry+registry.split(b'\n')[1]+b'\n').hex())),
        ('empty_skips_bad_registry',dict(missing_sources=True,registry=b'"bad'.hex())),
    ]:
        options.setdefault('registry',registry.hex());add('prepared_'+label,command='prepared',**options)
    for label,names,rows,reg in [
        ('fill',['gene_symbol','source_condition_id'],[['A','OMIM:1'],['B','ORPHA:1'],['A','absent']],registry),
        ('source_wins',['source_condition_id','mondo_id','priva_scope'],[['OMIM:1','MONDO:999','exclude'],['OMIM:1','','']],registry),
        ('duplicate_left',['source_condition_id'],[['OMIM:1'],['OMIM:1'],['ORPHA:1']],registry),
        ('empty_keys',['source_condition_id','mondo_id'],[['',''],['NA',''],['OMIM:1','NA']],registry+b'\tMONDO:empty\tscope\tinclude\treviewed\n'),
        ('first_duplicate',['source_condition_id'],[['OMIM:1']],registry+b'OMIM:1\tMONDO:last\tlast\texclude\tlast\n'),
        ('empty_skips',['source_condition_id'],[],b'broken\n'),
        ('missing_source_column',['gene_symbol'],[['A']],registry),
        ('missing_registry_columns',['source_condition_id'],[['OMIM:1']],b'disease_id\tmondo_id\nOMIM:1\tM\n'),
        ('bad_registry_row',['source_condition_id'],[['OMIM:1']],registry+b'A\tB\tC\tD\tE\textra\n'),
        ('implicit_registry_index',['source_condition_id'],[['OMIM:1']],encode(REGISTRY_HEADER,[['index','OMIM:1','M','scope','include','reviewed']],delimiter='\t')),
        ('early_width_late_decode',['source_condition_id'],[['OMIM:1']],registry+b'A\tB\tC\tD\tE\textra\n'+b'\n'*262144+b'\xff'),
    ]:
        add('attach_'+label,command='attach',payload=encode(names,rows,delimiter='\t'),registry=reg.hex())
    for kind in ('g2p','orphadata'):
        for file_kind in ('missing','directory','denied','symlink','dangling','loop'):
            add(kind+'_input_'+file_kind,kind=kind,payload=g2p() if kind=='g2p' else xml(),input_kind=file_kind)
        for output_kind in ('existing_different','existing_identical','existing_directory','denied_parent','temporary_directory','temporary_file','symlink','dangling','hardlink'):
            add(kind+'_output_'+output_kind,kind=kind,payload=g2p() if kind=='g2p' else xml(),output_kind=output_kind,
                unqualified='One historical FileExistsError versus IsADirectoryError remains unexplained; same-path producer/open/OS/native repetitions matched. Preserve this diagnostic boundary.' if kind=='g2p' and output_kind=='temporary_directory' else None)
        for suffix in ('.tsv','.tsv.gz','.zip','.tar.gz'):
            add(kind+'_output_suffix_'+suffix.replace('.','_'),kind=kind,payload=g2p() if kind=='g2p' else xml(),output_suffix=suffix)
    for code in (0x345,0x5b0,0x93e,0x2160,0xb2,0x200d,0x301,0x1e4d0,0x105c0):
        for side in ('prefix','suffix'):
            term=chr(code)+'dominant negative' if side=='prefix' else 'dominant negative'+chr(code)
            add('word_'+str(code)+'_'+side,payload=g2p(**{'molecular mechanism':term}))
    boundaries=set();previous=False
    for code in range(0x110001):
        current=code<0x110000 and chr(code).isalnum()
        if current!=previous:
            boundaries.update(c for c in (code-1,code) if 0<=c<0x110000 and not 0xd800<=c<=0xdfff)
        previous=current
    add('word_all_unicode_category_boundaries',payload=encode(['gene symbol','molecular mechanism'],
        [[str(code)+'_'+side,chr(code)+'dominant negative' if side=='prefix' else 'dominant negative'+chr(code)] for code in sorted(boundaries) for side in ('prefix','suffix')]))
    # A decoder delivers a complete read before the tokenizer can see its rows.
    # Competing failures test that order, including unfinished prefix records.
    for label,raw in [
        ('width_same_read',registry+b'A\tB\tC\tD\tE\textra\n'+b'\xff'),
        ('quoted_prefix',registry+b'A\tB\t"'+b'x'*262144+b'\xff'),
        ('partial_wide_prefix',registry+b'A\tB\tC\tD\tE\t'+b'x'*262144+b'\xff'),
        ('blank_prefix',b'\n'*262144+b'\xff'),
        ('valid_prefix',registry+b'\n'*262144+b'\xff'),
        ('width_multibyte_late',registry+b'A\tB\tC\tD\tE\textra\n'+('α'*270000).encode()+b'\xff'),
        ('width_multibyte_same_read',registry+b'A\tB\tC\tD\tE\textra\n'+('α'*140000).encode()+b'\xff'),
        ('width_crlf_late',registry.replace(b'\n',b'\r\n')+b'A\tB\tC\tD\tE\textra\r\n'+b'\r\n'*140000+b'\xff'),
    ]:
        add('failure_order_'+label,command='attach',payload=encode(['source_condition_id'],[['OMIM:1']],delimiter='\t'),registry=raw.hex())
    return result

def module_at(path):
    sys.path.insert(0,str(path.parent));spec=importlib.util.spec_from_file_location('mechanism_builder_oracle',path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module);return module

def oracle(args):
    request=json.loads(args.request.read_text());module=module_at(args.producer);paths=request['paths'];spec=request['spec']
    try:
        if spec['command']=='parse':frame=(module.parse_g2p if spec['kind']=='g2p' else module.parse_orphadata)(Path(paths['input']));module.write_tsv(Path(paths['output']),module.UNIFIED_COLUMNS,frame)
        elif spec['command']=='prepared':
            module.configure_source_urls(types.SimpleNamespace(g2p_ddg2p_url=spec.get('g2p_url'),orphadata_gene_disease_url=spec.get('orpha_url')))
            frame=module.parse_all_sources(Path(paths['raw-dir']),Path(paths['disease-scope-registry']) if 'disease-scope-registry' in paths else None)
            module.write_gene_pathogenic_mechanism_evidence(Path(paths['output']),frame)
        else:
            frame=module.pd.read_csv(paths['input'],sep='\t',dtype=str,low_memory=False).fillna('')
            frame=module.attach_mondo_condition_identity(frame,Path(paths['disease-scope-registry']));module.write_tsv(Path(paths['output']),list(frame.columns),frame)
    except Exception as exc:
        save(Path(request['exception']),dict(category=type(exc).__name__,message=str(exc)))
        print(type(exc).__name__+': '+str(exc),file=sys.stderr);return 1
    return 0

def admission_oracle(args):
    module=module_at(args.producer)
    parser=argparse.ArgumentParser();sub=parser.add_subparsers(dest='procedure',required=True)
    parse=sub.add_parser('parse');parse.add_argument('--kind',choices=['g2p','orphadata'],required=True)
    parse.add_argument('--input',required=True);parse.add_argument('--output',required=True)
    prepared=sub.add_parser('prepared');prepared.add_argument('--raw-dir',required=True);prepared.add_argument('--output',required=True)
    prepared.add_argument('--disease-scope-registry')
    for name in ('g2p-ddg2p-url','orphadata-gene-disease-url'):prepared.add_argument('--'+name,type=module.source_url)
    attach=sub.add_parser('attach')
    for name in ('input','output','disease-scope-registry'):attach.add_argument('--'+name,required=True)
    parser.parse_args(args.arguments[1:] if args.arguments[:1]==['--'] else args.arguments)
    return 0

def admission(args,env):
    folder=args.work_root/'admission';folder.mkdir();base=['prepared','--raw-dir','missing','--output','result.tsv'];cases=[
        [],['parse'],['parse','--kind','invalid'],['parse','--kind','invalid','--kind','g2p','--input','x','--output','result.tsv'],
        ['parse','--kind','invalid','--help'],['parse','--help','--kind','invalid'],['prepared','--help'],
        base+['--unknown','value'],['attach','--input','x'],
    ]
    for url in ('','relative','ftp://example.test/source','https://has space.test/source','https://[invalid]/source'):
        cases.append(base+['--g2p-ddg2p-url',url])
    cases.extend([base+['--g2p-ddg2p-url','invalid','--g2p-ddg2p-url',G2P_URL],
                  base+['--g2p-ddg2p-url','invalid','--help'],base+['--help','--g2p-ddg2p-url','invalid']])
    rows=[]
    for index,arguments in enumerate(cases):
        case=folder/str(index);case.mkdir();measurements={}
        for arm in ('python','native'):
            root=case/arm;root.mkdir()
            command=[sys.executable,Path(__file__).resolve(),'admission-oracle','--producer',args.producer,'--',*arguments] if arm=='python' else [args.binary,'gene-mechanism',*arguments]
            measurements[arm]=measure(command,root,env)
        rows.append(dict(arguments=arguments,measurements=measurements,passed=measurements['python']['exit']==measurements['native']['exit'] and not list(case.glob('*/result.tsv*'))))
    save(folder/'comparison.json',rows);return rows

def compressed(payload,suffix,layout='single'):
    if suffix=='.gz':return gzip.compress(payload,mtime=0)
    if suffix=='.bz2':return bz2.compress(payload)
    if suffix=='.xz':return lzma.compress(payload)
    if suffix=='.zst':
        import zstandard;return zstandard.ZstdCompressor().compress(payload)
    if suffix=='.zip':
        import zipfile
        out=io.BytesIO()
        names=[] if layout=='empty' else ['a.csv','b.csv'] if layout=='multiple' else ['dir/'] if layout=='directory' else ['a.csv','a.csv'] if layout=='duplicate' else ['α.csv'] if layout=='unicode' else ['input.csv']
        with zipfile.ZipFile(out,'w') as archive:
            for name in names:archive.writestr(zipfile.ZipInfo(name,(2000,1,1,0,0,0)),b'' if layout=='directory' else payload)
        return out.getvalue()
    import tarfile
    out=io.BytesIO()
    with tarfile.open(fileobj=out,mode='w') as archive:
        names=[] if layout=='empty' else ['a.csv','b.csv'] if layout=='multiple' else ['dir/'] if layout=='directory' else ['a.csv','a.csv'] if layout=='duplicate' else ['α.csv'] if layout=='unicode' else ['input.csv']
        for name in names:
            item=tarfile.TarInfo(name)
            if layout=='directory':item.type=tarfile.DIRTYPE;archive.addfile(item)
            else:item.size=len(payload);archive.addfile(item,io.BytesIO(payload))
    return gzip.compress(out.getvalue(),mtime=0) if suffix=='.tar.gz' else out.getvalue()

def snapshot(root):
    result={}
    for path in sorted(root.rglob('*')):
        if path.name in ('request.json','exception.json','stdout','stderr'):continue
        st=path.lstat();entry=dict(mode=stat.S_IMODE(st.st_mode),kind='link' if path.is_symlink() else 'directory' if path.is_dir() else 'file')
        if path.is_symlink():entry['target']=os.readlink(path)
        elif path.is_file():entry.update(bytes=st.st_size,sha256=pin(path)['sha256'] if st.st_mode & 0o400 else None,retained_mtime=st.st_mtime_ns==1700000000000000000)
        result[str(path.relative_to(root))]=entry
    return result

def prepare(root,spec,identical):
    root.mkdir();paths={};raw=bytes.fromhex(spec['payload']);input_path=root/('input'+spec.get('suffix','.csv' if spec['kind']=='g2p' else '.xml'))
    if spec.get('compress'):raw=compressed(raw,spec['compress'],spec.get('archive_layout','single'))
    file_kind=spec.get('input_kind','file')
    if file_kind=='directory':input_path.mkdir()
    elif file_kind=='loop':input_path.symlink_to(input_path.name)
    elif file_kind!='missing':
        if file_kind in ('symlink','dangling'):
            target=root/'input_target'
            if file_kind=='symlink':target.write_bytes(raw)
            input_path.symlink_to(target.name)
        else:
            input_path.write_bytes(raw)
            if file_kind=='denied':input_path.chmod(0)
    paths['input']=str(input_path)
    if spec['command']=='prepared':
        source_root=root/'raw';source_root.mkdir();(source_root/'orphadata').mkdir()
        if not spec.get('missing_sources'):
            (source_root/'AllG2P.csv').write_bytes(bytes.fromhex(spec.get('g2p_payload',g2p().hex())))
            (source_root/'orphadata/en_product6.xml').write_bytes(bytes.fromhex(spec.get('xml_payload',xml().hex())))
        paths={'raw-dir':str(source_root)}
    if spec['command'] in ('prepared','attach') and not spec.get('no_registry'):
        registry=root/'registry.tsv';kind=spec.get('registry_kind','file')
        if kind=='directory':registry.mkdir()
        elif kind!='missing':
            registry.write_bytes(bytes.fromhex(spec['registry']))
            if kind=='denied':registry.chmod(0)
        paths['disease-scope-registry']=str(registry)
    output=root/('output'+spec.get('output_suffix','.tsv'));kind=spec.get('output_kind','new')
    if kind=='denied_parent':
        parent=root/'blocked';parent.mkdir();parent.chmod(0o500);output=parent/output.name
    elif kind=='existing_directory':output.mkdir()
    elif kind in ('symlink','dangling'):
        target=root/'output_target'
        if kind=='symlink':target.write_bytes(b'preserved target\n')
        output.symlink_to(target.name)
    elif kind=='hardlink':
        target=root/'output_target';target.write_bytes(b'old linked\n');os.link(target,output)
    elif kind in ('existing_different','existing_identical'):
        output.write_bytes(identical if kind=='existing_identical' else b'old output\n');output.chmod(0o640);os.utime(output,ns=(1700000000000000000,1700000000000000000))
    temporary=Path(str(output)+'.tmp')
    if kind=='temporary_directory':temporary.mkdir()
    if kind=='temporary_file':temporary.write_bytes(b'old temp\n');temporary.chmod(0o600)
    paths['output']=str(output);return paths

def measure(command,root,env):
    started=time.monotonic()
    with (root/'stdout').open('xb') as out,(root/'stderr').open('xb') as err:
        child=subprocess.Popen(list(map(str,command)),cwd=root,env=env,stdout=out,stderr=err)
        _,status,usage=os.wait4(child.pid,0);child.returncode=os.waitstatus_to_exitcode(status)
    return dict(command=list(map(str,command)),exit=child.returncode,wall_seconds=time.monotonic()-started,peak_rss_kib=usage.ru_maxrss,
        stdout=pin(root/'stdout'),stderr=pin(root/'stderr'))

def check(args):
    args.work_root.mkdir(parents=True);env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONHASHSEED='0',POLARS_MAX_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    module=module_at(args.producer);identical={}
    seed=args.work_root/'identical-seeds';seed.mkdir()
    for kind,payload,fn in [('g2p',g2p(),module.parse_g2p),('orphadata',xml(),module.parse_orphadata)]:
        path=seed/kind;path.write_bytes(payload);target=seed/(kind+'.tsv');module.write_tsv(target,module.UNIFIED_COLUMNS,fn(path));identical[kind]=target.read_bytes()
    reuse={};reused=0
    if args.reuse_report:
        previous=json.loads(args.reuse_report.read_text())
        if previous['producer']['sha256']!=pin(args.producer)['sha256']:raise RuntimeError('Saved Python producer differs')
        if pin(previous['comparator']['path'])!=previous['comparator']:raise RuntimeError('Saved comparison source differs')
        for row in previous['cases']:
            if pin(row['receipt']['path'])!=row['receipt']:raise RuntimeError('Saved case receipt differs')
            value=json.loads(Path(row['receipt']['path']).read_text());reuse[row['case']]=value
    rows=[];failures=[];unqualified=[];raw_stderr_differences=0
    for name,spec in cases().items():
        folder=args.work_root/name;folder.mkdir();arms={}
        for arm in ('python','native'):
            if arm=='python' and name in reuse:
                old=reuse[name];old_spec={k:v for k,v in old['spec'].items() if k!='unqualified'};new_spec={k:v for k,v in spec.items() if k!='unqualified'}
                if old_spec!=new_spec:raise RuntimeError('Saved fixture differs: '+name)
                py=old['arms']['python'];old_root=Path(py['measurement']['stdout']['path']).parent
                if snapshot(old_root)!=py['after']:raise RuntimeError('Saved fixture/output changed: '+name)
                for stream in ('stdout','stderr'):
                    if pin(py['measurement'][stream]['path'])!=py['measurement'][stream]:raise RuntimeError('Saved measurement changed: '+name)
                arms[arm]=py;reused+=1;continue
            root=folder/arm;paths=prepare(root,spec,identical[spec['kind']]);request=root/'request.json';save(request,dict(paths=paths,spec=spec,exception=str(root/'exception.json')))
            before=snapshot(root)
            if arm=='native' and name in reuse and before!=reuse[name]['arms']['python']['before']:raise RuntimeError('New fixture differs from saved Python input: '+name)
            if arm=='python':command=[sys.executable,Path(__file__).resolve(),'oracle','--producer',args.producer,'--request',request]
            else:
                command=[args.binary,'gene-mechanism',spec['command']]
                if spec['command']=='parse':command+=['--kind',spec['kind']]
                for key,value in paths.items():command+=['--'+key,value]
                for key,flag in [('g2p_url','g2p-ddg2p-url'),('orpha_url','orphadata-gene-disease-url')]:
                    if key in spec:command+=['--'+flag,spec[key]]
            measured=measure(command,root,env);after=snapshot(root)
            stderr=(root/'stderr').read_text(errors='backslashreplace').replace(str(root),'<CASE>')
            normalized=stderr.removeprefix('priva: ') if stderr.startswith('priva: ') else stderr.replace('\npriva: ','\n')
            arms[arm]=dict(measurement=measured,before=before,after=after,stderr=stderr,normalized_stderr=normalized,
                exception=json.loads((root/'exception.json').read_text()) if (root/'exception.json').exists() else None)
        py,rs=arms['python'],arms['native'];differences=[]
        if py['measurement']['exit']!=rs['measurement']['exit']:differences.append('exit')
        if py['normalized_stderr']!=rs['normalized_stderr']:differences.append('stderr_category_message')
        if py['measurement']['stdout']['sha256']!=rs['measurement']['stdout']['sha256']:differences.append('stdout')
        if py['after']!=rs['after']:differences.append('output_bytes_modes_links_preservation')
        if py['stderr']!=rs['stderr']:raw_stderr_differences+=1
        row=dict(case=name,spec=spec,arms=arms,differences=differences,passed=not differences)
        save(folder/'comparison.json',row);rows.append(dict(case=name,passed=not differences,differences=differences,
            python_exit=py['measurement']['exit'],native_exit=rs['measurement']['exit'],receipt=pin(folder/'comparison.json')))
        if spec.get('unsupported') or spec.get('unqualified'):unqualified.append(rows[-1])
        elif differences:failures.append(rows[-1])
    admitted=admission(args,env)
    summary=dict(cases=len(rows),qualified_cases=len(rows)-len(unqualified),successes=sum(v['python_exit']==0 for v in rows if v not in unqualified),
        exact_failures=sum(v['python_exit']==1 and v['passed'] and v not in unqualified for v in rows),failures=len(failures),unqualified=len(unqualified),raw_stderr_differences=raw_stderr_differences,
        reused_python_cases=reused,admission_cases=len(admitted),admission_failed=sum(not v['passed'] for v in admitted))
    report=dict(kind='gene-mechanism-prepared-differential-v1',producer=pin(args.producer),comparator=pin(Path(__file__)),binary=pin(args.binary),summary=summary,
        cases=rows,failures=failures,unqualified=unqualified,admission=admitted,reference_report=pin(args.reuse_report) if args.reuse_report else None,
        passed=not failures and all(v['passed'] for v in admitted),stderr_envelope='Original function exceptions are captured without traceback; native priva prefix and per-arm fixture roots are normalized separately from raw-byte recording. Argument admission compares exact exits; usage wording is recorded without claiming byte equality.')
    save(args.work_root/'comparison.json',report);print(json.dumps(summary));return int(not report['passed'])

def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('oracle');p.add_argument('--producer',type=Path,required=True);p.add_argument('--request',type=Path,required=True);p.set_defaults(run=oracle)
    p=sub.add_parser('admission-oracle');p.add_argument('--producer',type=Path,required=True);p.add_argument('arguments',nargs=argparse.REMAINDER);p.set_defaults(run=admission_oracle)
    p=sub.add_parser('check');p.add_argument('--producer',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True);p.add_argument('--reuse-report',type=Path);p.set_defaults(run=check)
    args=parser.parse_args();return args.run(args)

if __name__=='__main__':raise SystemExit(main())
