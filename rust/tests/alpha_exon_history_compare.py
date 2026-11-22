#!/usr/bin/env python3
"""Compare native exon construction history to the unchanged Alpha collector.

compare --binary PRIVA --source-dir CHECKOUT --work-root NEW_DIRECTORY
oracle executes the actual serial/Pool.map source and one final pickle load.
Generated model pickles are test oracles only, never native producer inputs.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import struct
import subprocess
import sys

sys.dont_write_bytecode=True
MARKERS={'distribution','gnomAD_distribution','gnomAD_residue_distribution','min_distribution','max_distribution'}
PIN='cc7b2d60665e8052e6362107e359cc133e31c398f243589ee63ad7e8f90341c4'

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    value=importlib.util.module_from_spec(spec);sys.modules[name]=value;spec.loader.exec_module(value);return value

def scalar(value):
    if type(value) is int:return ['int',str(value)]
    if type(value) is float:return ['float',struct.pack('>d',value).hex()]
    raise TypeError(type(value).__name__)

def flatten(model):
    result=[]
    def walk(gene,ordinal,path,node):
        for order,(key,value) in enumerate(node.items()):
            base=dict(gene=gene,gene_order=ordinal,path=path,member_order=order,key=key)
            if key in MARKERS:
                result.append(dict(base,kind='history-marker'))
            elif isinstance(value,dict) and 'exons' in value:
                result.append(dict(base,kind='history-transcript',exons=list(value['exons']),
                    minimum=scalar(value['aa_min']),maximum=scalar(value['aa_max'])))
            elif isinstance(value,dict):
                branch=[*path,key];result.append(dict(base,path=branch,kind='history-branch'))
                walk(gene,ordinal,branch,value)
            else:raise TypeError('unexpected original history layout')
    for ordinal,(gene,node) in enumerate(model.items()):
        result.append(dict(gene=gene,gene_order=ordinal,path=[],member_order=0,key=gene,kind='history-gene'))
        walk(gene,ordinal,[],node)
    return result

def oracle(args):
    source=args.source_dir/'scripts/stat_protein_domain_amscores.py'
    if hashlib.sha256(source.read_bytes()).hexdigest()!=PIN:raise ValueError('collector source differs')
    original=module('alpha_exon_history_original',source)
    model=original.main(str(args.vcf),str(args.work_root/'original.pkl'),threads=args.threads)
    (args.work_root/'expected.json').write_text(json.dumps(flatten(model),ensure_ascii=False)+'\n')
    # Actual final serialization, followed by the exact original mapper. This
    # load is intentionally one additional stage after the worker transport.
    with (args.work_root/'original.pkl').open('rb') as stream:loaded=pickle.load(stream)
    mapper=module('alpha_exon_history_mapper',args.source_dir/'rust/tests/das_transcript_compare.py')
    mapping=mapper.mapping_functions(args.source_dir/'scripts/am_pick_intolerant_domains.py')['generate_transcript_pp_domain_map'](loaded)
    with (args.work_root/'mapping.pkl').open('xb') as stream:pickle.dump(mapping,stream,protocol=5)

def native_rows(path):
    import pyarrow.parquet as pq
    result=[]
    for batch in pq.ParquetFile(path).iter_batches():
        for row in batch.to_pylist():
            if row['item_order']!=0 or row['values'] or row['counts'] or row['typed_keys'] or row['auxiliary_int'] is not None:
                raise ValueError('history contains a numerical payload')
            selected={key:row[key] for key in ['gene','gene_order','path','member_order','key','kind']}
            if row['kind']=='history-transcript':
                if (row['minimum_int'] is None)==(row['minimum_float'] is None):raise ValueError('minimum scalar type absent/ambiguous')
                selected.update(exons=row['keys'],minimum=scalar(row['minimum_int'] if row['minimum_int'] is not None else row['minimum_float']),maximum=scalar(row['maximum']))
            elif row['keys'] or row['minimum_int'] is not None or row['minimum_float'] is not None or row['maximum'] is not None:
                raise ValueError('non-transcript row has transcript payload')
            result.append(selected)
    return result

def specs(alpha):
    chosen=['normal','empty','unicode','parent_multiplicity','duplicates','transcript_order','empty_protein_position',
        'empty_exon','hierarchy_empty_ids','skipped_sources','first_score_max_af','header_order','nonadjacent_gene',
        'variant_tuple_layout','tuple_nan_parent_multiplicity','tuple_nan_chromosomes']
    values={name:dict(alpha.cases()[name]) for name in chosen}
    values['gene_pending_order']=dict(contigs=['chr1','chr2','chr3'],rows=[
        alpha.row(csq=alpha.annotation(Gene='ENSG_A',EXON='0/30')),
        alpha.row(pos=11,csq=alpha.annotation(Gene='ENSG_B',EXON='1/30')),
        alpha.row(chrom='chr2',csq=alpha.annotation(Gene='ENSG_C',EXON='2/30')),
        alpha.row(chrom='chr3',csq=alpha.annotation(Gene='ENSG_A',EXON='3/30'))])
    values['overlapping_sets']=dict(contigs=['chr3','chr1','chr2','empty'],rows=[
        alpha.row(chrom=chrom,pos=i+1,variant=f'A{i}T',csq=alpha.annotation(EXON=f'{i}/100',Protein_position=f'{i}/100',
                  DOMAINS='DB:A:B:C&DB:Q&DB:A',Feature='ENST_B' if i%2 else 'ENST_A'))
        for chrom,start,end in [('chr1',0,37),('chr2',17,60),('chr3',35,78)] for i in range(start,end)])
    for size in [4,5,18,19,76,77,306,307,1228,1229]:
        values[f'set_boundary_{size}']=dict(rows=[alpha.row(pos=i+1,variant=f'A{i}T',
            csq=alpha.annotation(EXON=f'{i}/2000',Protein_position=f'{i+1}/2000')) for i in range(size)])
    return values

def failure_specs(alpha):
    selected=['missing_file','short_csq','score_tuple','af_scalar_float','variant_integer','csq_integer',
        'csq_none_after_prefix','no_format','last_other_format','negative_pos','no_eof','crc','truncated',
        'plain_indexed_route','no_index_parallel','empty_no_header_contigs_parallel']
    return {name:alpha.cases()[name] for name in selected}

def category(text):
    import re
    kinds=re.findall(r'(?:^|\n)(?:priva: )?([A-Za-z]+(?:Error|Layout)): ',text)
    return kinds[-1] if kinds else None

def history_text(rows):
    encode=lambda value:value.encode().hex() or '-'
    lines=[]
    for row in rows:
        if row['kind']!='history-transcript':continue
        path=[row['gene'],*row['path'],row['key'],'exons']
        words=[str(len(path)),*map(encode,path),str(len(row['exons'])),*map(encode,row['exons'])]
        lines.append(' '.join(words)+'\n')
    return ''.join(lines)

def compare(args):
    import pyarrow.parquet as pq
    root=args.work_root;root.mkdir(parents=True,exist_ok=False)
    alpha=module('alpha_history_cases',args.source_dir/'rust/tests/alpha_domain_compare.py')
    cache=module('alpha_history_cache_compare',args.source_dir/'rust/tests/cache_values_compare.py')
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',POLARS_MAX_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    differences=[];results=[];failures=[];boundaries=[]
    def run(folder,label,command,seed=0):
        with (folder/(label+'.stdout.log')).open('x') as out,(folder/(label+'.stderr.log')).open('x') as err:
            return subprocess.run(list(map(str,command)),env=dict(env,PYTHONHASHSEED=str(seed)),stdout=out,stderr=err,timeout=120).returncode
    all_specs=specs(alpha)
    if args.case:all_specs={name:all_specs[name] for name in args.case}
    for name,spec in all_specs.items():
        case=root/name;case.mkdir();vcf=alpha.make_case(case,spec)
        contigs=spec.get('contigs',['chr1']);parts=[]
        for ordinal,chrom in enumerate(contigs):
            part=case/f'raw{ordinal}.parquet';parts.append(part)
            status=run(case,f'raw{ordinal}',[args.binary,'alpha-exon-history','chromosome','--vcf',vcf,'--chrom',chrom,'--output',part])
            if status:raise RuntimeError(f'raw fixture failed: {name}/{chrom}: {status}')
        for route,threads in [('serial',1),('pool',2)]:
            statistics=case/(route+'.statistics.parquet')
            status=run(case,route+'_statistics',[args.binary,'alphamissense-domain','build','--vcf',vcf,'--threads',threads,'--output',statistics])
            if status:raise RuntimeError(f'native statistics fixture failed: {name}/{route}: {status}')
            source=pq.ParquetFile(statistics).metadata.metadata[b'priva.builder_source_sha256'].decode()
            for seed in args.seed or [0,1,42,4294967295]:
                folder=case/f'{route}_{seed}';folder.mkdir();native=folder/'history.parquet';export=folder/'history.native'
                statuses=[run(folder,'original',[sys.executable,__file__,'oracle','--source-dir',args.source_dir,
                    '--vcf',vcf,'--threads',threads,'--work-root',folder],seed)]
                command=[args.binary,'alpha-exon-history', 'build' if route=='serial' else 'merge']
                command+=['--vcf',vcf] if route=='serial' else [word for path in parts for word in ['--input',path]]
                command+=['--hash-seed',seed,'--output',native]
                statuses.append(run(folder,'native',command,seed))
                if any(statuses):failures.append(dict(case=name,route=route,seed=seed,stage='producer',statuses=statuses));continue
                wanted=json.loads((folder/'expected.json').read_text());actual=native_rows(native)
                delta=[]
                for i in range(max(len(wanted),len(actual))):
                    a=wanted[i] if i<len(wanted) else None;b=actual[i] if i<len(actual) else None
                    if a!=b:delta.append(dict(row=i,original=a,native=b))
                (folder/'history_mismatches.jsonl').write_text(''.join(json.dumps(row,ensure_ascii=False)+'\n' for row in delta))
                if delta:differences.append(dict(case=name,route=route,seed=seed,rows=len(delta)))
                status=run(folder,'export',[args.binary,'alpha-exon-history','export','--input',native,'--output',export],seed)
                exported=export.read_text().split('\n',1) if export.exists() else []
                if status or len(exported)!=2 or exported[1]!=history_text(wanted):
                    failures.append(dict(case=name,route=route,seed=seed,stage='export',status=status));continue
                output=folder/'mapping.rkyv'
                command=[args.binary,'das-domains','transcripts','--statistics',statistics,'--exon-history',export,
                    '--source-builder-sha256',source,'--python-hash-seed',seed,'--output',output]
                status=run(folder,'transcripts',command,seed)
                if status:failures.append(dict(case=name,route=route,seed=seed,stage='transcripts',status=status));continue
                comparison=cache.compare(argparse.Namespace(binary=args.binary,cache='alphamissense_tranx_domain_map',pickle=folder/'mapping.pkl',
                    archive=output,report=folder/'mapping_comparison.json',mismatches=folder/'mapping_mismatches.jsonl'))
                results.append(dict(case=name,route=route,seed=seed,rows=len(wanted),sets=sum(row['kind']=='history-transcript' for row in wanted),
                    history_matched=not delta,mapping_comparison=comparison))
                if not comparison['passed']:failures.append(dict(case=name,route=route,seed=seed,stage='mapping',mismatches=comparison['mismatches']))
    if not args.case:
        for name,spec in failure_specs(alpha).items():
            folder=root/('failure_'+name);folder.mkdir();vcf=alpha.make_case(folder,spec)
            threads=int(spec.get('threads','1'));output=folder/'history.parquet';sentinel=b'previous complete publication\n';output.write_bytes(sentinel)
            original=run(folder,'original',[sys.executable,__file__,'oracle','--source-dir',args.source_dir,'--vcf',vcf,
                '--threads',threads,'--work-root',folder])
            command=[args.binary,'alpha-exon-history','build','--vcf',vcf,'--hash-seed',0,'--output',output]
            if threads>1:
                # Admission to an indexed chromosome must retain the original
                # missing-index/plain-input failures. Zero contigs is tested by
                # the original native Alpha command (no chromosome is valid).
                command=[args.binary,'alpha-exon-history','chromosome','--vcf',vcf,'--chrom','chr1','--output',output]
            native=run(folder,'native',command)
            categories=[category((folder/(label+'.stderr.log')).read_text()) for label in ['original','native']]
            matched=original==native==1 and categories[0] is not None and categories[0]==categories[1] and output.read_bytes()==sentinel
            boundaries.append(dict(case=name,original_status=original,native_status=native,categories=categories,matched=matched,previous_output_preserved=output.read_bytes()==sentinel))
            if not matched:failures.append(dict(case=name,stage='failure_boundary',original_status=original,native_status=native,categories=categories))
    expected=len(all_specs)*2*len(args.seed or [0,1,42,4294967295])
    result=dict(passed=not differences and not failures and len(results)==expected,cases=len(results),expected=expected,
        mismatches=differences,failures=failures,boundaries=boundaries,results=results,full_real_history_qualified=False,full_das_qualified=False)
    (root/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['passed','cases','expected','mismatches','failures']}),flush=True)
    if not result['passed']:raise SystemExit(1)

def guards(args):
    """Reject incomplete/corrupt native companions without touching good files."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    root=args.work_root;root.mkdir(parents=True,exist_ok=False)
    case=args.fixtures/'normal';final=case/'serial_0/history.parquet';raw=case/'raw0.parquet'
    history=case/'serial_0/history.native';statistics=case/'serial.statistics.parquet'
    source=pq.ParquetFile(statistics).metadata.metadata[b'priva.builder_source_sha256'].decode()
    # Polars stores its profile in Parquet footer key/value entries rather than
    # Arrow's optional serialized schema metadata. Preserve that exact footer
    # when constructing corrupt-row fixtures, so row validation is reached.
    table=pq.read_table(final).replace_schema_metadata(pq.ParquetFile(final).metadata.metadata);results=[]
    env=dict(os.environ,POLARS_MAX_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    def run(name,command,expected_code=1,expected_category='ValueError',output=True):
        target=root/(name+'.output');sentinel=b'previous complete publication\n'
        if output:target.write_bytes(sentinel);command=[*command,'--output',target]
        with (root/(name+'.stdout.log')).open('x') as out,(root/(name+'.stderr.log')).open('x') as err:
            status=subprocess.run(list(map(str,[args.binary,*command])),env=env,stdout=out,stderr=err,timeout=30).returncode
        observed=category((root/(name+'.stderr.log')).read_text())
        preserved=not output or target.read_bytes()==sentinel
        passed=status==expected_code and (expected_category is None or observed==expected_category) and preserved
        results.append(dict(case=name,passed=passed,status=status,expected_status=expected_code,category=observed,
            expected_category=expected_category,previous_output_preserved=preserved))
    for name,field,value in [('wrong_schema',b'priva.history.schema',b'unknown'),('wrong_source',b'priva.producer',b'0'*64),
        ('wrong_builder',b'priva.builder_source_sha256',b'0'*64),('wrong_stage',b'priva.stage',b'raw'),
        ('wrong_profile',b'priva.history.profile',b'unknown'),('invalid_seed',b'priva.history.seed',b'-1')]:
        metadata=dict(table.schema.metadata);metadata[field]=value;altered=root/(name+'.parquet')
        pq.write_table(table.replace_schema_metadata(metadata),altered,compression='zstd')
        run(name,['alpha-exon-history','export','--input',altered])
    for name in ['duplicate_member','unknown_kind','gene_ordinal','minimum_type','extra_payload','missing_parent','member_order']:
        rows=table.to_pylist();transcript=next(row for row in rows if row['kind']=='history-transcript')
        if name=='duplicate_member':transcript['keys']*=2
        elif name=='unknown_kind':transcript['kind']='unknown'
        elif name=='gene_ordinal':rows[0]['gene_order']=1
        elif name=='minimum_type':transcript['minimum_float']=1.0
        elif name=='extra_payload':transcript['values']=[1.0]
        elif name=='missing_parent':transcript['path']=['missing','parent']
        elif name=='member_order':transcript['member_order']+=1
        altered=root/(name+'.parquet');pq.write_table(pa.Table.from_pylist(rows,schema=table.schema),altered,compression='zstd')
        run(name,['alpha-exon-history','export','--input',altered],
            expected_category='UnsupportedLayout' if name in ['duplicate_member','missing_parent'] else 'ValueError')
    run('raw_export',['alpha-exon-history','export','--input',raw])
    multi=args.fixtures/'overlapping_sets';parts=[multi/f'raw{i}.parquet' for i in range(4)]
    for name,selected in [('missing_shard',parts[:-1]),('reordered_shards',list(reversed(parts))),('duplicate_shards',[parts[0]]*4)]:
        run(name,['alpha-exon-history','merge',*[word for path in selected for word in ['--input',path]],'--hash-seed','0'])
    for name in ['wrong_input','wrong_history_seed','wrong_history_profile','wrong_history_source','missing_line','extra_line','changed_member']:
        contents=history.read_text();header,body=contents.split('\n',1);metadata=json.loads(header[2:])
        if name=='wrong_input':metadata['priva.input']='{}'
        elif name=='wrong_history_seed':metadata['priva.history.seed']='-1'
        elif name=='wrong_history_profile':metadata['priva.history.profile']='unknown'
        elif name=='wrong_history_source':metadata['priva.builder_source_sha256']='0'*64
        elif name=='missing_line':body=''
        elif name=='extra_line':body+=body.splitlines(keepends=True)[0]
        elif name=='changed_member':body=body.rsplit(' ',1)[0]+' 74616d7065726564\n'
        altered=root/(name+'.native');altered.write_text('# '+json.dumps(metadata)+'\n'+body)
        run(name,['das-domains','transcripts','--statistics',statistics,'--exon-history',altered,
            '--source-builder-sha256',source,'--python-hash-seed','0'])
    run('admission_seed',['alpha-exon-history','build','--vcf',case/'input.vcf.gz','--hash-seed','-1'],2,None)
    result=dict(passed=all(row['passed'] for row in results),cases=len(results),results=results,
        scope='native companion integrity/admission, separate from original Python failures')
    (root/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    if not result['passed']:raise SystemExit(1)

def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    check=commands.add_parser('compare')
    for name in ['binary','source-dir','work-root']:check.add_argument('--'+name,type=Path,required=True)
    check.add_argument('--case',action='append');check.add_argument('--seed',type=int,action='append');check.set_defaults(function=compare)
    reference=commands.add_parser('oracle')
    for name in ['source-dir','vcf','work-root']:reference.add_argument('--'+name,type=Path,required=True)
    reference.add_argument('--threads',type=int,required=True);reference.set_defaults(function=oracle)
    guard=commands.add_parser('guards')
    for name in ['binary','fixtures','work-root']:guard.add_argument('--'+name,type=Path,required=True)
    guard.set_defaults(function=guards)
    args=parser.parse_args();args.function(args)

if __name__=='__main__':main()
