#!/usr/bin/env python3
"""Compare the unchanged Python annotation contract with native procedures.

  python clinvar_contract_compare.py check --source CONTRACT.py --binary PRIVA --work-root NEW_DIRECTORY

Fixtures are synthetic. Cache conversion admission is reported separately from
shape equivalence: the Python shape checker deliberately accepts more value
types than the existing native archive schemas. No producer is executed.
"""
from __future__ import annotations
import argparse
import copy
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import re
import stat
import subprocess
import sys


def load(path):
    spec=importlib.util.spec_from_file_location('contract_reference',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def captured(function):
    try: return {'success':True,'result':function()}
    except Exception as error: return {'success':False,'category':type(error).__name__,'message':str(error)}


def run(command, folder, label):
    result=subprocess.run(list(map(str,command)),capture_output=True)
    (folder/(label+'.stdout.log')).write_bytes(result.stdout)
    (folder/(label+'.stderr.log')).write_bytes(result.stderr)
    match=re.match(r'priva: ([A-Za-z]+): ',result.stderr.decode(errors='replace'))
    value={'exit':result.returncode,'success':result.returncode==0}
    if match:
        value['category']=match[1]
        value['message']=result.stderr.decode(errors='replace')[match.end():].rstrip('\n')
    if result.returncode==0 and result.stdout.strip():value['result']=json.loads(result.stdout)
    return value


def paired(args,module,folder,name,operation,options,function,original_cli=None,outputs=()):
    before={p:(p.read_bytes(),stat.S_IMODE(p.stat().st_mode)) if p.is_file() and not p.is_symlink() else None for p in outputs}
    expected=captured(function)
    if original_cli is not None:
        cli=run([sys.executable,args.source,*original_cli],folder,name+'.python')
        expected['exit']=cli['exit']
        if cli['success']!=expected['success']:raise RuntimeError('reference function and CLI disagree: '+name)
    else:expected['exit']=0 if expected['success'] else 1
    expected_outputs={str(p): (p.read_bytes().hex(),stat.S_IMODE(p.stat().st_mode)) if p.is_file() and not p.is_symlink() else None for p in outputs}
    for p,item in before.items():
        if item is None:
            if p.is_file() and not p.is_symlink():p.unlink()
        else:
            p.write_bytes(item[0]);p.chmod(item[1])
    command=[args.binary,'clinvar-contract',operation]
    for key,value in options.items():command += ['--'+key,str(value)]
    actual=run(command,folder,name+'.native')
    actual_outputs={str(p):(p.read_bytes().hex(),stat.S_IMODE(p.stat().st_mode)) if p.is_file() and not p.is_symlink() else None for p in outputs}
    ok=(expected['exit']==actual['exit'] and expected['success']==actual['success'])
    if expected['success']:
        ok &= expected.get('result')==actual.get('result') and expected_outputs==actual_outputs
    else:
        ok &= expected['category']==actual.get('category') and expected_outputs==actual_outputs
    check_message=not expected['success'] and expected['message'].startswith(('ClinVar ','Incomplete ClinVar ',
        'Missing/non-file/empty ClinVar ','pair-capture requires'))
    if check_message:
        expected_message=expected['message']
        if operation=='cache':expected_message=expected_message.replace(original_cli[0],str(options['input']))
        ok &= expected_message==actual.get('message')
    value={'name':name,'scope':operation,'python':expected,'native':actual,'output_state_equal':expected_outputs==actual_outputs,'passed':bool(ok)}
    value['source_message_checked']=check_message
    # Fixture payloads stay beside their logs; the report does not duplicate them.
    (folder/(name+'.comparison.json')).write_text(json.dumps(value,indent=2)+'\n')
    return value


def headers(args,module,root):
    import pysam
    folder=root/'headers';folder.mkdir()
    csq=sorted(module.AA_CSQ_FIELDS|module.GENE_CSQ_FIELDS|module.AF_CSQ_FIELDS)
    infos=sorted(module.TABLE_INFO_FIELDS|module.AF_INFO_FIELDS|{'CLNSIG','CLNREVSTAT','ALLELEID'})
    prefix='##fileformat=VCFv4.2\n##contig=<ID=2>\n##contig=<ID=1>\n'
    info=''.join(f'##INFO=<ID={n},Number=.,Type=String,Description="fixture">\n' for n in infos)
    description='##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: '+'|'.join(csq)+'">\n'
    columns='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    text=prefix+info+description+columns
    cases={'valid':text,'bad_body':text+'unparsed body\n','empty':'','no_fileformat':text.split('\n',1)[1],
        'no_columns':text[:-len(columns)],'crlf':text.replace('\n','\r\n'),'no_final_newline':text.rstrip('\n'),
        'duplicate_contig':text.replace('##contig=<ID=1>','##contig=<ID=2>\n##contig=<ID=1>'),
        'no_description':text.replace(description,'##INFO=<ID=CSQ,Number=.,Type=String>\n'),
        'unquoted_description':text.replace(description,description.replace('"','')),
        'no_format':text.replace('Format: ','Format:'),'trailing_field':text.replace(description,description.replace('">','|">')),
        'duplicate_field':text.replace(description,description.replace('">','|Gene">')),
        'double_format':text.replace(description,description.replace('">','|Format: second">')),
        'escaped_suffix':text.replace(description,description.replace('">','|escaped\\"quote\\\\tail">')),
        'leading_space_attribute':text.replace(description,description.replace('Description',' Description')),
        'duplicate_description_attribute':text.replace(description,description.replace('Description=', 'Description="ignored",Description=')),
        'duplicate_info_valid_first':text.replace(description,description+'##INFO=<ID=CSQ,Number=.,Type=String,Description="wrong">\n'),
        'duplicate_info_absent_description_first':text.replace(description,'##INFO=<ID=CSQ,Number=.,Type=String>\n'+description),
        'missing_number':text.replace('ID=CSQ,Number=.,','ID=CSQ,'),
        'missing_type':text.replace('ID=CSQ,Number=.,Type=String,','ID=CSQ,Number=.,'),
        'format_no_sample':text.replace(columns,columns.rstrip('\n')+'\tFORMAT\n'),
        'empty_sample':text.replace(columns,columns.rstrip('\n')+'\tFORMAT\t\n'),
        'duplicate_sample':text.replace(columns,columns.rstrip('\n')+'\tFORMAT\tS\tS\n')}
    for name in infos+['CSQ']:
        cases['missing_info_'+name]=''.join(line+'\n' for line in text.splitlines() if not line.startswith('##INFO=<ID='+name+','))
    for name in csq:cases['missing_csq_'+name]=text.replace(description,description.replace('Format: '+'|'.join(csq),'Format: '+'|'.join(v for v in csq if v!=name)))
    rows=[]
    def check(name,path,stage='complete'):
        rows.append(paired(args,module,folder,name,'header',{'vcf':path,'stage':stage},lambda:module.validate_header(str(path),stage),[str(path),'--stage',stage]))
    for name,content in cases.items():
        path=folder/(name+'.vcf');path.write_text(content);check(name,path)
    original=folder/'valid.vcf'
    for stage in ['aa','gene','af','vep']:check('stage_'+stage,original,stage)
    for stage in ['aa','gene','af','vep','complete']:
        # A minimal stage-specific header detects accidental requirement widening.
        req_info={'CSQ','CLNSIG','CLNREVSTAT'};req_csq=set()
        if stage in ['aa','vep','complete']:req_csq.update(module.AA_CSQ_FIELDS)
        if stage in ['gene','vep','complete']:req_csq.update(module.GENE_CSQ_FIELDS);req_info.add('ALLELEID')
        if stage in ['af','complete']:req_csq.update(module.AF_CSQ_FIELDS);req_info.update(module.AF_INFO_FIELDS)
        if stage in ['gene','complete']:req_info.update(module.TABLE_INFO_FIELDS)
        content=prefix+''.join(f'##INFO=<ID={n},Number=.,Type=String,Description="fixture">\n' for n in sorted(req_info-{'CSQ'}))
        content+='##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: '+'|'.join(sorted(req_csq))+'">\n'+columns
        path=folder/('minimal_'+stage+'.vcf');path.write_text(content);check('minimal_'+stage,path,stage)
    bgzf=folder/'valid.bgz';pysam.tabix_compress(str(original),str(bgzf));check('bgzf',bgzf)
    path=folder/'missing_eof.bgz';path.write_bytes(bgzf.read_bytes()[:-28]);check('missing_eof',path)
    path=folder/'ordinary.gzip';path.write_bytes(gzip.compress(text.encode()));check('ordinary_gzip',path)
    bcf=folder/'valid.bcf'
    with pysam.VariantFile(str(original)) as source:
        with pysam.VariantFile(str(bcf),'wb',header=source.header):pass
    check('bcf',bcf)
    path=folder/'raw.bcf';path.write_bytes(gzip.decompress(bcf.read_bytes()));check('raw_bcf',path)
    check('missing',folder/'absent');check('directory',folder)
    return rows


SCHEMAS={'aa':'clinvar_aa_stat','splice':'clinvar_splice_stat','gene-af':'clinvar_patho_af_stat','exon-af':'clinvar_patho_exon_af_stat','gene':'clinvar_gene_stat'}
def caches(args,module,root):
    folder=root/'caches';folder.mkdir()
    sets={k:set() for k in ['clinvar_variants','clinvar_pathogenic','clinvar_not_pathogenic','putative_nmd_variants','small_aachange','large_aachange','no_aachange']}
    gene=dict(gene_symbol='S',ensembl_id='G',**sets,transcripts={},domains={})
    splice=dict(chrom='1',pos=1,ref='A',alt='T',exon='1/2',intron='',hgvsc='x',consequence='splice',clinvar_sig='P',clinvar_review='R',splice_ai={'SpliceAI_pred_DS_AG':'0.2'})
    valid={'aa':{'T':{'1':{'p.Ala1Val':{'CLNSIG':[],'CLNREVSTAT':['P']}}}},'splice':{'T':[splice]},
        'gene-af':{'G':{'chrom':'1','pos':1,'af':0.0,'acc_patho_af':-0.0},'N':None},
        'exon-af':{'T':{'1/2':(0.0,1.0,2.0,3.0,(4.0,5.0))}},'gene':{'G':gene}}
    specs=[]
    for kind,data in valid.items():specs += [(kind+'_valid',kind,data),(kind+'_empty',kind,{})]
    specs += [('aa_empty_positions','aa',{'T':{}}),('aa_empty_changes','aa',{'T':{'1':{}}}),
        ('splice_empty_records','splice',{'T':[]}),('splice_nullable','splice',{'T':[dict(splice,alt=None,exon=None,consequence=None,splice_ai={})]}),
        ('gene_af_empty_record','gene-af',{'G':{}}),('gene_af_missing_field','gene-af',{'G':{'af':1.0}}),
        ('gene_af_none_only','gene-af',{'G':None}),('gene_af_boolean','gene-af',{'G':dict(valid['gene-af']['G'],af=True)}),
        ('exon_empty','exon-af',{'T':{}}),('gene_nontext','gene',{'G':dict(gene,gene_symbol=3.0)})]
    rows=[];admission=[]
    for name,kind,data in specs:
        source=folder/(name+'.pkl');source.write_bytes(pickle.dumps(data,protocol=4))
        archive=folder/(name+'.rkyv')
        converted=run([args.binary,'cache-convert',SCHEMAS[kind],'--pickle',source,'--output',archive],folder,name+'.convert')
        if not converted['success']:raise RuntimeError('representable cache conversion failed: '+name)
        def oracle(source=source,kind=kind,archive=archive):
            value=module.validate_cache(source,kind);value['path']=str(archive);return value
        rows.append(paired(args,module,folder,name,'cache',{'input':archive,'kind':kind},oracle,[str(source),'--cache-kind',kind]))
    # Existing archive representation boundaries, not helper shape mismatches.
    for name,kind,data in [('aa_untyped_values','aa',{'T':{'1':{'x':{'CLNSIG':1,'CLNREVSTAT':None}}}}),
        ('aa_extra_field','aa',{'T':{'1':{'x':{'CLNSIG':[],'CLNREVSTAT':[],'extra':1}}}}),
        ('splice_minimum_keys','splice',{'T':[{'chrom':None,'pos':None,'clinvar_sig':None,'splice_ai':None}]}),
        ('gene_af_none_scalar','gene-af',{'G':dict(valid['gene-af']['G'],af=None)}),
        ('gene_plain_integer_symbol','gene',{'G':dict(gene,gene_symbol=3)})]:
        source=folder/(name+'.pkl');source.write_bytes(pickle.dumps(data,protocol=4))
        archive=folder/(name+'.rkyv')
        expected=captured(lambda:module.validate_cache(source,kind))
        actual=run([args.binary,'cache-convert',SCHEMAS[kind],'--pickle',source,'--output',archive],folder,name+'.convert')
        admission.append(dict(name=name,python_shape_success=expected['success'],native_conversion_exit=actual['exit'],
            scope='existing native schema admission; no native archive produced',passed=expected['success'] and actual['exit']==2 and not archive.exists()))
    return rows,admission


def metadata(args,module,root):
    folder=root/'metadata';folder.mkdir();rows=[]
    def setup(name):
        work=folder/name;work.mkdir()
        options={k:str(work/k) for k in ['vcf','producer','aa','splice','receipt','captured','config']}
        options['contract-source']=str(args.source)
        for name in ['vcf','producer','aa','splice','config']:Path(options[name]).write_text(name+' bytes\n')
        Path(options['vcf']+'.tbi').write_text('tbi\n')
        return work,options
    def python_options(options):return [str(x) for pair in options.items() if pair[0]!='contract-source' for x in ['--'+pair[0],pair[1]]]
    def pair(name,operation,change=None,prepare=False):
        work,options=setup(name)
        if prepare:
            module.pair_command(['pair-capture',*python_options(options)])
            module.pair_command(['pair-complete',*python_options(options)])
        if change:change(work,options)
        outputs=[Path(options[k]) for k in ['receipt','captured'] if k in options]
        argv=[operation,*python_options(options)]
        rows.append(paired(args,module,work,name,operation,options,lambda:module.pair_command(argv),argv,outputs))
    pair('capture','pair-capture')
    pair("capture space ' Ω 😀",'pair-capture')
    pair('complete','pair-complete',prepare=True)
    pair('check','pair-check',prepare=True)
    pair('capture_missing_option','pair-capture',lambda w,o:o.pop('captured'))
    pair('capture_no_index_first','pair-capture',lambda w,o:(o.pop('captured'),Path(o['vcf']+'.tbi').unlink()))
    pair('capture_no_outputs_needed','pair-capture',lambda w,o:(Path(o['aa']).unlink(),Path(o['splice']).unlink()))
    pair('complete_missing_option','pair-complete',lambda w,o:o.pop('captured'))
    pair('complete_missing_aa_first','pair-complete',lambda w,o:(o.pop('captured'),Path(o['aa']).unlink()))
    pair('check_changed_producer','pair-check',lambda w,o:Path(o['producer']).write_text('changed'),prepare=True)
    pair('check_new_csi','pair-check',lambda w,o:Path(o['vcf']+'.csi').write_text('csi'),prepare=True)
    pair('complete_changed_input','pair-complete',lambda w,o:Path(o['vcf']).write_text('changed input'),prepare=True)
    pair('capture_both_indexes','pair-capture',lambda w,o:Path(o['vcf']+'.csi').write_text('csi'))
    for label,payload in [('extra','{"extra": 3}'),('extra_bigint','{"extra": '+str(10**100)+'}'),
        ('extra_surrogate','{"extra": "\\ud800"}'),('extra_surrogate_nonpair','{"extra": "\\ud800\\u0041"}'),
        ('extra_surrogate_pair','{"extra": "\\ud83d\\ude00"}'),('extra_surrogate_key','{"\\ud800": 1}'),
        ('extra_digit_limit','{"extra": '+'9'*4300+'}'),('extra_over_digit_limit','{"extra": '+'9'*4301+'}'),
        ('extra_invalid','{"extra": [}'),('extra_deep','{"extra":'+'['*450+'0'+']'*450+'}'),
        ('extra_very_deep','{"extra":'+'['*2000+'0'+']'*2000+'}'),
        *[('extra_depth_'+str(n),'{"extra":'+'['*n+'0'+']'*n+'}') for n in [9000,9995,9996,9997,9998,11000]]]:
        def edit(w,o,payload=payload):
            target=Path(o['receipt']);text=target.read_text().rstrip();target.write_text(text[:-1]+','+payload[1:])
        pair('check_'+label,'pair-check',edit,prepare=True)
    for label,payload in [('true','true'),('float','1.0'),('wrong','2'),('nan','NaN'),('bigint','9'*100),('surrogate','"\\ud800"')]:
        def edit(w,o,payload=payload):
            target=Path(o['receipt']);target.write_text(target.read_text().replace('"schema": 1','"schema": '+payload))
        pair('check_schema_'+label,'pair-check',edit,prepare=True)
    for label,payload in [('array','[]'),('null','null'),('integer','3'),('invalid','{'),('invalid_utf8',b'\xff')]:
        def edit(w,o,payload=payload):
            target=Path(o['receipt']);target.write_bytes(payload if isinstance(payload,bytes) else payload.encode())
        pair('check_json_'+label,'pair-check',edit,prepare=True)
    for operation in ['pair-capture','pair-check','pair-complete']:
        for key,target in [('aa','vcf'),('splice','aa'),('aa','config'),('aa','producer'),('aa','contract-source'),('receipt','vcf'),('captured','receipt')]:
            pair(operation+'_alias_'+key+'_'+target,operation,lambda w,o,key=key,target=target:o.__setitem__(key,o[target]))
        for label,target in [('live','vcf'),('broken','absent')]:
            pair(operation+'_metadata_symlink_'+label,operation,lambda w,o,target=target:Path(o['receipt']).symlink_to(w/target))
        pair(operation+'_hardlink',operation,lambda w,o:(Path(o['aa']).unlink(),os.link(o['vcf'],o['aa'])))
        pair(operation+'_duplicate_before_metadata',operation,lambda w,o:(o.__setitem__('aa',o['vcf']),o.__setitem__('captured',o['receipt'])))
        pair(operation+'_metadata_before_hardlink',operation,lambda w,o:(Path(o['aa']).unlink(),os.link(o['vcf'],o['aa']),o.__setitem__('captured',o['receipt'])))
        pair(operation+'_missing_parent',operation,lambda w,o:o.__setitem__('captured',str(w/'missing'/'capture.json')))
    for label in ['ordinary','empty','directory','missing','broken','loop','relative','dotdot','double_slash']:
        work,options=setup('metadata_'+label);path=Path(options['vcf'])
        if label=='empty':path.write_bytes(b'')
        if label=='directory':path=work
        if label=='missing':path=work/'absent'
        if label in ['broken','loop']:
            path=work/label;path.symlink_to('absent' if label=='broken' else 'loop')
        value=str(path)
        if label=='relative':value=os.path.relpath(path)
        if label=='dotdot':value=str(work)+'/../'+work.name+'/vcf'
        if label=='double_slash':value='/'+str(path)
        rows.append(paired(args,module,work,'metadata_'+label,'pair-metadata',{'input':value},lambda value=value:module.pair_metadata(value)))
    # AF's preflight does not require any file to exist; all six names matter.
    work,options=setup('af_paths')
    af={k:options[k] for k in ['vcf','config']};af.update(gene=str(work/'gene'),exon=str(work/'exon'))
    argv=['af-paths-check',*[str(x) for k,v in af.items() for x in ['--'+k,v]]]
    rows.append(paired(args,module,work,'af_paths','af-paths-check',af,lambda:module.af_paths_command(argv[1:]),argv))
    return rows


def check(args):
    args.work_root.mkdir(parents=True)
    module=load(args.source)
    rows=headers(args,module,args.work_root)
    cache_rows,admission=caches(args,module,args.work_root);rows+=cache_rows
    rows+=metadata(args,module,args.work_root)
    report={'source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(),
        'binary_sha256':hashlib.sha256(args.binary.read_bytes()).hexdigest(),'cases':rows,'schema_admission':admission,
        'passed':all(r['passed'] for r in rows+admission)}
    (args.work_root/'complete.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':report['passed'],'cases':len(rows),'mismatches':[r['name'] for r in rows if not r['passed']],
                      'schema_admission_controls':len(admission)}))
    if not report['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    sub=commands.add_parser('check')
    for option in ['source','binary','work-root']:sub.add_argument('--'+option,type=Path,required=True)
    sub.set_defaults(function=check)
    args=parser.parse_args();args.function(args)


if __name__=='__main__':main()
