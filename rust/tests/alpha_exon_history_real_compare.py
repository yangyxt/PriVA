#!/usr/bin/env python3
"""Independent CPython operations over complete native insertion histories.

compare --source-dir CHECKOUT --raw-part PARQUET [--raw-part PARQUET ...]
        --final-parquet PARQUET --serialized-history LINES
        --work-root NEW_DIRECTORY

The unchanged source merge and transcript mapper consume metadata only. Actual
CPython set.add, pickle worker transport, set.update and final unpickle retain
their distinct stages. This oracle does not read the legacy score pickle or
execute the VCF collector; raw history membership/ranges need their own gate.
"""
import argparse
import hashlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import pickle
import struct
import sys
sys.dont_write_bytecode=True

PIN='cc7b2d60665e8052e6362107e359cc133e31c398f243589ee63ad7e8f90341c4'

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);value=importlib.util.module_from_spec(spec)
    sys.modules[name]=value;spec.loader.exec_module(value);return value

def rows(path):
    import pyarrow.parquet as pq
    for batch in pq.ParquetFile(path).iter_batches(batch_size=16384):
        yield from batch.to_pylist()

def raw_genes(path):
    import pyarrow.parquet as pq
    metadata=pq.ParquetFile(path).metadata.metadata;keys=json.loads(metadata[b'priva.keys'])
    current=None;seen=0
    for row in rows(path):
        if row['kind']=='history-gene':
            if current is not None:yield current[0],current[2]
            if row['gene_order']!=seen or row['gene']!=keys[seen] or row['key']!=row['gene'] or row['path']:
                raise ValueError('raw history gene order differs')
            current=(row['gene'],row['gene_order'],{});seen+=1;continue
        if current is None or current[:2]!=(row['gene'],row['gene_order']):raise ValueError('history gene rows are noncontiguous')
        branch=row['kind']=='history-branch';path=row['path'][:-1] if branch else row['path'];node=current[2]
        for key in path:node=node[key]
        key=row['key']
        if row['member_order']!=len(node) or key in node or row['item_order']!=0:raise ValueError('raw history member order differs')
        if row['kind'] in ['history-branch','history-marker']:node[key]={}
        elif row['kind']=='history-transcript':
            sequence=row['keys'];exons=set()
            for exon in sequence:
                if exon in exons:raise ValueError('duplicate first-insertion event')
                exons.add(exon)
            node[key]=dict(exons=exons,aa_min=row['minimum_int'] if row['minimum_int'] is not None else row['minimum_float'],aa_max=row['maximum'])
        else:raise ValueError('unexpected raw history kind')
    if current is not None:yield current[0],current[2]
    if seen!=len(keys):raise ValueError('incomplete raw history gene inventory')

def normalized(row):
    if row['item_order'] or row['values'] or row['counts'] or row['typed_keys'] or row['auxiliary_int'] is not None:
        raise ValueError('final history contains numeric payload')
    selected={key:row[key] for key in ['gene','gene_order','path','member_order','key','kind']}
    if row['kind']=='history-transcript':
        if (row['minimum_int'] is None)==(row['minimum_float'] is None):raise ValueError('ambiguous minimum type')
        minimum=['int',str(row['minimum_int'])] if row['minimum_int'] is not None else ['float',struct.pack('>d',row['minimum_float']).hex()]
        selected.update(exons=row['keys'],minimum=minimum,maximum=['int',str(row['maximum'])])
    elif row['keys'] or row['minimum_int'] is not None or row['minimum_float'] is not None or row['maximum'] is not None:
        raise ValueError('unexpected transcript fields')
    return selected

def compare(args):
    import pyarrow.parquet as pq
    root=args.work_root;root.mkdir(parents=True,exist_ok=False)
    seed=os.environ.get('PYTHONHASHSEED')
    if seed is None or not seed.isdigit() or not 0<=int(seed)<=2**32-1:raise ValueError('explicit Python hash seed required')
    collector=args.source_dir/'scripts/stat_protein_domain_amscores.py'
    if hashlib.sha256(collector.read_bytes()).hexdigest()!=PIN:raise ValueError('collector source changed')
    original=module('alpha_real_history_original',collector)
    fixture=module('alpha_real_history_fixture',args.source_dir/'rust/tests/alpha_exon_history_compare.py')
    mapper=module('alpha_real_history_mapper',args.source_dir/'rust/tests/das_transcript_compare.py')
    headers=[pq.ParquetFile(path).metadata.metadata for path in args.raw_part]
    first=headers[0];contigs=json.loads(first[b'priva.contigs'])
    if len(contigs)!=len(headers):raise ValueError('every raw header chromosome is required')
    ordered=[];last={}
    for ordinal,header in enumerate(headers):
        if header[b'priva.stage']!=b'raw' or header[b'priva.ordinal']!=str(ordinal).encode() or header[b'priva.history.profile']!=b'first-insertion-v1':
            raise ValueError('raw history stage/order differs')
        for key in [b'priva.input',b'priva.producer',b'priva.contigs',b'priva.builder_source_sha256',b'priva.history.schema']:
            if header[key]!=first[key]:raise ValueError('raw history source/profile differs')
        for gene in json.loads(header[b'priva.keys']):
            if gene not in last:ordered.append(gene)
            last[gene]=ordinal
    final=pq.ParquetFile(args.final_parquet).metadata.metadata
    if final[b'priva.stage']!=b'final' or final[b'priva.history.profile']!=b'cpython312-pool-v1' or final[b'priva.history.seed']!=seed.encode():
        raise ValueError('wrong final history seed/serialization route')
    if json.loads(final[b'priva.keys'])!=ordered:raise ValueError('final history gene inventory differs')
    for key in [b'priva.input',b'priva.producer',b'priva.contigs',b'priva.builder_source_sha256',b'priva.history.schema']:
        if final[key]!=first[key]:raise ValueError('final history source/profile differs')
    native_rows=iter(rows(args.final_parquet));counts=dict(genes=0,rows=0,sets=0,members=0,worker_genes=0);mismatches=0
    with args.serialized_history.open() as native, (root/'expected_serialized.native').open('x') as expected, (root/'mismatches.jsonl').open('x') as ledger:
        header=native.readline()
        if not header.startswith('# '):raise ValueError('native history export header absent')
        exported=json.loads(header[2:])
        for key in [b'priva.input',b'priva.producer',b'priva.history.schema',b'priva.history.profile',b'priva.history.seed',b'priva.builder_source_sha256']:
            if exported[key.decode()]!=final[key].decode():raise ValueError('native export profile differs')
        def mismatch(kind,wanted,actual):
            nonlocal mismatches
            mismatches+=1;ledger.write(json.dumps(dict(kind=kind,at=dict(counts),expected=wanted,native=actual),ensure_ascii=False)+'\n')
        def completed_genes():
            pending={};emitted=0
            for ordinal,path in enumerate(args.raw_part):
                for gene,node in raw_genes(path):
                    # multiprocessing's ForkingPickler uses default protocol 4
                    # in CPython 3.12; different genes do not share exon sets.
                    worker=pickle.loads(pickle.dumps({gene:node},protocol=4));counts['worker_genes']+=1
                    original._merge_scores(pending,worker)
                while emitted<len(ordered) and last[ordered[emitted]]<=ordinal:
                    gene=ordered[emitted];node=pending.pop(gene);emitted+=1
                    original._finalize_scores_dict({gene:node})
                    plain=original.to_plain_dict({gene:node})
                    projected=fixture.flatten(plain)
                    for wanted in projected:
                        wanted['gene_order']=counts['genes']
                        actual=next(native_rows,None);actual=normalized(actual) if actual is not None else None
                        if wanted!=actual:mismatch('typed_history_row',wanted,actual)
                        counts['rows']+=1
                        if wanted['kind']=='history-transcript':
                            line=fixture.history_text([wanted]);observed=native.readline();expected.write(line)
                            if line!=observed:mismatch('serialized_exon_sequence',line,observed)
                            counts['sets']+=1;counts['members']+=len(wanted['exons'])
                    counts['genes']+=1
                    # One final consumer load, using the actual C unpickler.
                    loaded=pickle.loads(pickle.dumps(plain,protocol=pickle.HIGHEST_PROTOCOL))
                    yield gene,loaded[gene]
            if pending or emitted!=len(ordered):raise ValueError('incomplete original merge')
        class Genes:
            def items(self):return completed_genes()
        mapping=mapper.mapping_functions(args.source_dir/'scripts/am_pick_intolerant_domains.py')['generate_transcript_pp_domain_map'](Genes())
        for extra in native:mismatch('extra_serialized_history',None,extra)
        for extra in native_rows:mismatch('extra_typed_row',None,normalized(extra))
    with (root/'mapping.pkl').open('xb') as stream:pickle.dump(mapping,stream,protocol=5)
    result=dict(passed=mismatches==0,python_hash_seed=int(seed),counts=counts,transcripts=len(mapping),mismatches=mismatches,
        source_sha256=PIN,raw_history_vcf_order_qualification='separate source-bound collector/metadata comparison',
        original_full_score_model_loaded=False,original_vcf_collector_rerun=False,legacy_pickle_input=False,full_das_qualified=False)
    (root/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
    if not result['passed']:raise SystemExit(1)

def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    check=commands.add_parser('compare')
    for name in ['source-dir','final-parquet','serialized-history','work-root']:check.add_argument('--'+name,type=Path,required=True)
    check.add_argument('--raw-part',type=Path,required=True,action='append');check.set_defaults(function=compare)
    args=parser.parse_args();args.function(args)

if __name__=='__main__':main()
