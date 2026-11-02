#!/usr/bin/env python3
"""Independent retained-input DAS transcript-cache comparisons.

compare --binary PRIVA --projection-probe PROBE --source-dir CHECKOUT
        --work-root NEW_DIRECTORY
oracle runs only the unchanged Alpha collector and transcript mapping functions.
No complete DAS computation or installed cache is invoked.
"""
import argparse
import ast
import importlib.util
import json
import logging
import os
from pathlib import Path
import pickle
import struct
import subprocess
import sys
from typing import List


def module(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    value=importlib.util.module_from_spec(spec);sys.modules[name]=value;spec.loader.exec_module(value);return value


def mapping_functions(path):
    names={'process_domain_dict_for_mapping','generate_transcript_pp_domain_map'}
    body=[node for node in ast.parse(path.read_text()).body if isinstance(node,ast.FunctionDef) and node.name in names]
    if len(body)!=2:raise RuntimeError('original transcript functions absent')
    logger=logging.getLogger('transcript_oracle');logger.disabled=True
    scope=dict(List=List,logger=logger)
    exec(compile(ast.fix_missing_locations(ast.Module(body=body,type_ignores=[])),str(path),'exec'),scope)
    return scope


def oracle(args):
    source=module('alpha_transcript_fixture_source',args.collector)
    source.main(str(args.vcf),str(args.statistics),threads=1)
    with args.statistics.open('rb') as stream:model=pickle.load(stream)
    mapping=mapping_functions(args.das_source)['generate_transcript_pp_domain_map'](model)
    with args.mapping.open('xb') as stream:pickle.dump(mapping,stream,protocol=5)


def metadata_genes(statistics, exon_history, counts):
    """Independent Arrow projection of only fields read by the original mapper.

    Distribution contents are never inspected by these two original functions.
    Their dictionary positions and presence are retained, and every transcript
    gets the actual C-unpickler's set from the saved ADDITEMS operand sequence.
    One gene is held at a time; no numerical samples are expanded or loaded.
    """
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    file=pq.ParquetFile(statistics);metadata=file.metadata.metadata
    if metadata[b'priva.schema']!=b'priva-alpha-domain-parquet-v2' or metadata[b'priva.stage']!=b'final':
        raise ValueError('expected final typed Alpha model')
    genes=json.loads(metadata[b'priva.keys']);current=None
    columns=['gene','gene_order','path','member_order','item_order','key','kind','keys','minimum_int','minimum_float','maximum']
    with exon_history.open() as stream:
        for batch in file.iter_batches(batch_size=65536,columns=columns):
            counts['rows']+=len(batch)
            selected=batch.filter(pc.not_equal(batch.column(batch.schema.get_field_index('kind')),'residue'))
            for row in selected.to_pylist():
                gene,kind,key=row['gene'],row['kind'],row['key']
                if kind=='gene':
                    if current is not None:yield current[0],current[2]
                    if row['gene_order']!=counts['genes'] or gene!=genes[counts['genes']] or key!=gene or row['path']:
                        raise ValueError('gene source order differs')
                    current=(gene,row['gene_order'],{});counts['genes']+=1;continue
                if current is None or current[:2]!=(gene,row['gene_order']):raise ValueError('noncontiguous model gene')
                node=current[2];path=row['path'][:-1] if kind in ['source','domain'] else row['path']
                for part in path:node=node[part]
                if key in node or row['member_order']!=len(node) or row['item_order']!=0:
                    raise ValueError('member identity or order differs')
                if kind in ['source','domain','residues']:node[key]={}
                elif kind=='array':node[key]=None
                elif kind=='transcript':
                    words=stream.readline().split()
                    if not words:raise ValueError('exon history ended early')
                    decode=lambda word:bytes.fromhex(word if word!='-' else '').decode()
                    n=int(words[0]);wanted=[gene,*row['path'],key,'exons']
                    if n!=len(wanted) or [decode(word) for word in words[1:n+1]]!=wanted:
                        raise ValueError('exon path or mutation order differs')
                    size=int(words[n+1]);items=[decode(word) for word in words[n+2:]]
                    if len(items)!=size or len(set(items))!=size or sorted(items)!=row['keys']:
                        raise ValueError('exon membership differs')
                    # Do not use pickle.dumps(set(items)): that adds an extra
                    # observable set-iteration/serialization layer.
                    payload=bytearray(b'\x80\x05\x8f(')
                    for item in items:
                        encoded=item.encode();payload+=b'X'+struct.pack('<I',len(encoded))+encoded
                    payload+=b'\x90.'
                    exons=pickle.loads(payload)
                    node[key]=dict(exons=exons,aa_min=row['minimum_int'] if row['minimum_int'] is not None else row['minimum_float'],aa_max=row['maximum'])
                    counts['sets']+=1;counts['members']+=size
                else:raise ValueError('unexpected typed metadata kind: '+kind)
        if stream.readline():raise ValueError('extra exon history')
    if current is not None:yield current[0],current[2]
    if counts['genes']!=len(genes) or counts['rows']!=file.metadata.num_rows:
        raise ValueError('incomplete gene or row coverage')


def project(args):
    """Run the unchanged mapper over a streaming metadata-only dictionary view."""
    seed=os.environ.get('PYTHONHASHSEED')
    if seed is None or not seed.isdigit() or not 0<=int(seed)<=2**32-1:
        raise ValueError('explicit Python hash seed required')
    counts=dict(rows=0,genes=0,sets=0,members=0)
    class Genes:
        def items(self):return metadata_genes(args.statistics,args.exon_history,counts)
    mapping=mapping_functions(args.das_source)['generate_transcript_pp_domain_map'](Genes())
    with args.mapping.open('xb') as stream:pickle.dump(mapping,stream,protocol=5)
    result=dict(hash_seed=int(seed),counts=counts,transcripts=len(mapping),numerical_samples_loaded=False,
                fresh_native_install_qualified=False,full_das_qualified=False)
    args.report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


def mapping_differences(comparator, original, projected):
    with original.open('rb') as stream:left=pickle.load(stream)
    with projected.open('rb') as stream:right=pickle.load(stream)
    for key in sorted(left.keys()|right.keys()):
        if key not in left or key not in right:
            yield dict(key=key,field='root_key',original=key in left,projected=key in right);continue
        a=comparator.reference_value('alphamissense_tranx_domain_map',key,left[key])
        b=comparator.reference_value('alphamissense_tranx_domain_map',key,right[key])
        for path,wanted,actual in comparator.differences(a,b):
            yield dict(key=key,path=path,original=wanted,projected=actual)


def compare(args):
    import pyarrow.parquet as pq
    root=args.work_root;root.mkdir(parents=True,exist_ok=False)
    alpha=module('alpha_transcript_fixture_inputs',args.source_dir/'rust/tests/alpha_domain_compare.py')
    comparator=module('transcript_cache_independent_compare',args.source_dir/'rust/tests/cache_values_compare.py')
    specs={key:value for key,value in alpha.cases().items() if key in {
        'normal','empty','parent_multiplicity','duplicates','transcript_order','empty_protein_position',
        'empty_exon','hierarchy_empty_ids','unicode','residue_min_max'}}
    specs['pp_stable_ties']=dict(rows=[
        alpha.row(csq=alpha.annotation(DOMAINS='DB:B',EXON='3/9',Protein_position='5/9')),
        alpha.row(pos=11,csq=alpha.annotation(DOMAINS='DB:A',EXON='1/9',Protein_position='5/9')),
        alpha.row(pos=12,csq=alpha.annotation(DOMAINS='DB:B',EXON='2/9',Protein_position='8/9'))])
    specs['leaf_stops_nested']=dict(rows=[alpha.row(csq=alpha.annotation(DOMAINS='DB:A:B:C&DB:A&DB:X:Y',EXON='4/5'))])
    specs['joined_path_collisions']=dict(rows=[
        alpha.row(csq=alpha.annotation(Gene='ENSG:G',DOMAINS='DB:A',Protein_position='8/9')),
        alpha.row(pos=11,csq=alpha.annotation(Gene='ENSG',DOMAINS='G:DB:A',Protein_position='2/9'))])
    results=[];failures=[];boundaries=[];projection_cases=0
    def run(label,command,folder,env):
        with (folder/(label+'.stdout.log')).open('x') as out,(folder/(label+'.stderr.log')).open('x') as err:
            return subprocess.run(list(map(str,command)),env=env,stdout=out,stderr=err,timeout=120).returncode
    for case,spec in specs.items():
        for seed in [0,1,42,4294967295]:
            folder=root/f'{case}_{seed}';folder.mkdir()
            vcf=alpha.make_case(folder,spec);model=folder/'original.pkl';mapping=folder/'mapping.pkl'
            native=folder/'statistics.parquet';history=folder/'exons.native';output=folder/'mapping.rkyv'
            env=dict(os.environ,PYTHONHASHSEED=str(seed),POLARS_MAX_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
            command=[sys.executable,__file__,'oracle','--collector',args.source_dir/'scripts/stat_protein_domain_amscores.py',
                '--das-source',args.source_dir/'scripts/am_pick_intolerant_domains.py','--vcf',vcf,'--statistics',model,'--mapping',mapping]
            statuses=[run('original',command,folder,env)]
            statuses.append(run('alpha',[args.binary,'alphamissense-domain','build','--vcf',vcf,'--output',native],folder,env))
            statuses.append(run('history',[args.projection_probe,'extract','--input',model,'--output',history],folder,env))
            if any(statuses):
                failures.append(dict(case=case,seed=seed,stage='prerequisites',statuses=statuses));continue
            projected=folder/'projected.pkl'
            status=run('projected',[sys.executable,__file__,'project','--statistics',native,'--exon-history',history,
                '--das-source',args.source_dir/'scripts/am_pick_intolerant_domains.py','--mapping',projected,
                '--report',folder/'projection.json'],folder,env)
            if status:
                failures.append(dict(case=case,seed=seed,stage='oracle_projection',status=status));continue
            differences=list(mapping_differences(comparator,mapping,projected))
            (folder/'projection_mismatches.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in differences))
            if differences:failures.append(dict(case=case,seed=seed,stage='oracle_projection',mismatches=len(differences)))
            projection_cases+=1
            source=pq.ParquetFile(native).metadata.metadata[b'priva.builder_source_sha256'].decode()
            command=[args.binary,'das-domains','transcripts','--statistics',native,'--exon-history',history,
                '--source-builder-sha256',source,'--python-hash-seed',str(seed),'--output',output]
            status=run('native',command,folder,env)
            if status:
                failures.append(dict(case=case,seed=seed,stage='native',status=status));continue
            result=comparator.compare(argparse.Namespace(binary=args.binary,cache='alphamissense_tranx_domain_map',pickle=mapping,
                archive=output,report=folder/'comparison.json',mismatches=folder/'mismatches.jsonl'))
            results.append(dict(case=case,seed=seed,passed=result['passed'],comparison=result))
            if not result['passed']:failures.append(dict(case=case,seed=seed,stage='values',mismatches=result['mismatches']))
            if case=='normal' and seed==0:
                variants=[('source_identity','--source-builder-sha256','0'*64,1,'ValueError'),
                    ('invalid_seed','--python-hash-seed','-1',2,None),
                    ('missing_history','--exon-history',folder/'missing.native',1,'FileNotFoundError')]
                contents=history.read_text();lines=contents.splitlines(keepends=True)
                for label,text in [('empty_history',''),('extra_history',contents+lines[-1]),
                    ('wrong_member',contents.rstrip().rsplit(' ',1)[0]+' 74616d7065726564\n')]:
                    altered=folder/(label+'.native');altered.write_text(text)
                    variants.append((label,'--exon-history',altered,1,'ValueError'))
                for label,option,value,expected_code,category in variants:
                    failed_output=folder/(label+'.rkyv');sentinel=b'previous publication\n';failed_output.write_bytes(sentinel)
                    changed=list(command);changed[changed.index(option)+1]=value;changed[changed.index('--output')+1]=failed_output
                    code=run(label,changed,folder,env);error=(folder/(label+'.stderr.log')).read_text()
                    passed=code==expected_code and failed_output.read_bytes()==sentinel and (category is None or category+':' in error)
                    boundaries.append(dict(case=label,passed=passed,returncode=code,expected_code=expected_code,category=category,previous_output_preserved=failed_output.read_bytes()==sentinel))
                    if not passed:failures.append(dict(case=label,stage='boundary',returncode=code,expected_code=expected_code,category=category,error=error))
    (root/'failures.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in failures))
    result=dict(passed=not failures and len(results)==len(specs)*4,cases=len(results),projection_cases=projection_cases,expected_cases=len(specs)*4,boundaries=boundaries,
        failures=len(failures),results=results,fresh_native_install_qualified=False,full_das_qualified=False)
    (root/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['passed','cases','expected_cases','failures']}),flush=True)
    if not result['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    command=commands.add_parser('compare')
    for name in ['binary','projection-probe','source-dir','work-root']:command.add_argument('--'+name,type=Path,required=True)
    command.set_defaults(function=compare)
    oracle_cmd=commands.add_parser('oracle')
    for name in ['collector','das-source','vcf','statistics','mapping']:oracle_cmd.add_argument('--'+name,type=Path,required=True)
    oracle_cmd.set_defaults(function=oracle)
    project_cmd=commands.add_parser('project')
    for name in ['statistics','exon-history','das-source','mapping','report']:project_cmd.add_argument('--'+name,type=Path,required=True)
    project_cmd.set_defaults(function=project)
    args=parser.parse_args();args.function(args)


if __name__=='__main__':main()
