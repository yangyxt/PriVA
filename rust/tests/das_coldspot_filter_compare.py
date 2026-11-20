#!/usr/bin/env python3
"""Compare the assembled native coldspot procedure with unchanged Python.

check --binary PRIVA --source ORIGINAL.py --lookup-source protein_domain_mapping.py
      --normalization-json JSON --work-root NEW_DIRECTORY --report JSON
oracle --source ORIGINAL.py --lookup-source SOURCE --normalization-json JSON
       --mapping PICKLE --domains JSONL --metadata JSON --spots TSV --threshold FLOAT

The independent oracle executes the complete source function and real pandas
reader. Native inputs use the existing typed metadata and schema-9 archives.
Every requested summary value/type/float bit/order and defined failure category
is compared. Enclosing DAS candidate-set construction remains a separate gate.
"""
import argparse
import ast
from collections import defaultdict
import gzip
import hashlib
import importlib.util
import json
import logging
import os
from pathlib import Path
import pickle
import re
import shutil
import struct
import subprocess
import sys

sys.dont_write_bytecode = True
HEADER = 'type\tgene\taa_start_pos\taa_end_pos\n'
DOMAINS = ['G:Pfam:PF1', 'G:Pfam:PF2', 'G:SMART:SM1', 'G:Pfam:PF3',
           'G:Unknown:X', 'G:Pfam:PF1&SMART:SM1', 'G:Pfam:PF4', 'G:Pfam:PF5',
           'G:Pfam:PF1', 'MISSING:Pfam:PF1', 'H:Pfam:PF1']


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def original(path):
    import pandas as pd
    if hashlib.sha256(path.read_bytes()).hexdigest() not in [
            '37502dedf58f44eb67237c462a9a0c19a84bc696764ecc870edabcaa60b5e92c',
            '1a2b4efc312c41c6fb20afa3c1d90566da7cd0772d2f261c2f7f09d094fae834']:
        raise RuntimeError('original DAS source changed')
    body = [node for node in ast.parse(path.read_text()).body
            if isinstance(node, ast.FunctionDef) and node.name == 'filter_coldspot_overlap']
    if len(body) != 1:
        raise RuntimeError('original coldspot function absent')
    logger = logging.getLogger('coldspot_original')
    logger.disabled = True
    scope = dict(pd=pd, defaultdict=defaultdict, logger=logger)
    exec(compile(ast.fix_missing_locations(ast.Module(body=body, type_ignores=[])), str(path), 'exec'), scope)
    return scope['filter_coldspot_overlap']


def oracle(args):
    helper = load(Path(__file__).with_name('das_domain_lookup_compare.py'), 'source_lookup')
    normalizer = helper.normalizer(args.lookup_source, json.loads(args.normalization_json.read_text()))
    with args.mapping.open('rb') as stream:
        mapping = pickle.load(stream)
    domains = [json.loads(line) for line in args.domains.read_text().splitlines()]
    metadata = json.loads(args.metadata.read_text())
    filtered, excluded, summary = original(args.source)(set(domains), metadata, str(args.spots),
                                                       normalizer, mapping, args.threshold)
    if filtered != set(domains)-excluded or any(summary[d]['removed'] != (d in excluded) for d in domains):
        raise RuntimeError('source filter result membership is inconsistent')
    for domain in domains:
        print(json.dumps(dict(domain=domain, **summary[domain])), flush=True)


def cases():
    normal = HEADER+'coldspot\tG\t1\t10\nhotspot\tG\t1\t999\ncoldspot\tH\t40\t50\n'
    out = dict(normal={}, empty_candidates=dict(domains=[]), header_only=dict(text=HEADER),
               absent_gene=dict(symbol='ABSENT'), numeric_gene=dict(symbol='1', text=HEADER+'coldspot\t1\t1\t10\n'),
               numeric_text_gene=dict(symbol='1', text=HEADER+'coldspot\t1\t1\t10\ncoldspot\tG\t1\t10\n'),
               boolean_gene=dict(symbol='True', text=HEADER+'coldspot\tTrue\t1\t10\n'),
               unicode=dict(symbol='α', text=HEADER+'coldspot\tα\t1\t10\n'),
               tab_newline_symbol=dict(symbol='G\tA\nB', text=HEADER+'coldspot\t"G\tA\nB"\t1\t10\n'),
               repeated_spots=dict(text=normal+normal.split('\n', 1)[1]),
               reversed_spot=dict(text=HEADER+'coldspot\tG\t10\t1\n'),
               zero_length_spot=dict(text=HEADER+'coldspot\tG\t10\t9\n'),
               reversed_domain=dict(bounds=[9, 1]), negative_domain=dict(bounds=[-20, -1], text=HEADER+'coldspot\tG\t-10\t-1\n'),
               all_integer_endpoints=dict(bounds=[-2**63, 2**63-1], text=HEADER+f'coldspot\tG\t{-2**63}\t{2**63-1}\n'),
               one_residue=dict(bounds=[7, 7], text=HEADER+'coldspot\tG\t7\t7\n'),
               touching=dict(bounds=[10, 20]), disjoint=dict(bounds=[11, 20]),
               missing_file=dict(missing=True), empty_file=dict(text=''),
               bad_coordinate=dict(text=HEADER+'coldspot\tG\tbad\t10\n'),
               late_bad_coordinate=dict(text=normal+'coldspot\tOTHER\tbad\t10\n'),
               empty_candidates_bad_coordinate=dict(domains=[], text=HEADER+'coldspot\tG\tbad\t10\n'),
               exempt_bad_coordinate=dict(domains=['G:SMART:SM1'], text=HEADER+'coldspot\tG\tbad\t10\n'),
               nonselected_bad_coordinate=dict(text=normal+'hotspot\tG\tbad\tbad\n'),
               gzip=dict(compression='gzip'), gzip_truncated=dict(compression='truncated'),
               crlf=dict(text=normal.replace('\n', '\r\n')), bare_cr=dict(text=normal.replace('\n', '\r')),
               comment_is_data=dict(text=normal+'#coldspot\tG\tbad\tbad\n'),
               implicit_index=dict(text=HEADER+'idx\tcoldspot\tG\t1\t10\n'),
               surplus_later=dict(text=normal+'coldspot\tG\t1\t10\textra\n'))
    for value in ['0', '-0.0', '-1', '0.5', '0.6', '1', '1.1', 'NaN', 'inf', '-inf']:
        out['threshold_'+str(len(out))] = dict(threshold=value)
    for value in ['5', '10', '1.9', '-1.9', 'NaN', 'inf', '-inf', '1_0', '١٠', '9'*4300, '9'*4301]:
        out['coordinate_'+str(len(out))] = dict(text=HEADER+'coldspot\tG\t1\t'+value+'\n')
    for exponent in [64, 65, 128, 1022, 1023, 1074, 1075, 1076, 4096, 14000]:
        for delta in [-1, 0, 1]:
            out[f'large_denominator_{exponent}_{delta}'] = dict(bounds=[1, 1], text=HEADER+f'coldspot\tG\t1\t{2**exponent+delta}\n')
    for numerator in [1, 3, 5, 7, 2**53+1, 2**63-1]:
        out['subnormal_numerator_'+str(numerator)] = dict(bounds=[1, numerator], text=HEADER+f'coldspot\tG\t1\t{2**1075}\n')
    for exponent in [64, 1023, 1075, 14000]:
        out['signed_huge_'+str(exponent)] = dict(bounds=[-5, 5], text=HEADER+f'coldspot\tG\t{-2**exponent}\t{2**exponent}\n')
    for name in ['type', 'gene', 'aa_start_pos', 'aa_end_pos']:
        columns = HEADER.rstrip().split('\t')
        values = ['coldspot', 'G', '1', '10']
        at = columns.index(name)
        del columns[at], values[at]
        out['missing_'+name] = dict(text='\t'.join(columns)+'\n'+'\t'.join(values)+'\n')
        out['missing_'+name+'_unselected'] = dict(text='\t'.join(columns)+'\n'+'\t'.join(values).replace('coldspot', 'hotspot')+'\n')
    for value in out.values():
        value.setdefault('text', normal)
    return out


def category(path):
    text = path.read_text(errors='replace')
    matches = re.findall(r'^priva: ([A-Za-z_][A-Za-z_0-9]*(?:Error|Exception)):', text, re.M)
    if not matches:
        matches = re.findall(r'^(?:[A-Za-z_][A-Za-z_0-9]*\.)*([A-Za-z_][A-Za-z_0-9]*(?:Error|Exception)):', text, re.M)
    return matches[-1] if matches else None


def typed(value):
    if type(value) is dict:
        return ['dict', [[typed(key), typed(item)] for key, item in value.items()]]
    if type(value) is float:
        return ['float', struct.pack('>d', value).hex()]
    if type(value) is list:
        return ['list', [typed(item) for item in value]]
    return [type(value).__name__, value]


def identity(path):
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.file_digest(path.open('rb'), 'sha256').hexdigest())


def run(command, folder, label, env=None):
    with (folder/(label+'.stdout.log')).open('x') as output, (folder/(label+'.stderr.log')).open('x') as error:
        return subprocess.run(list(map(str, command)), stdout=output, stderr=error, env=env, timeout=60).returncode


def fixture_metadata(template, metadata, output):
    import pyarrow as pa
    import pyarrow.parquet as pq
    table = pq.read_table(template)
    rows = []
    original_row = table.to_pylist()[0]
    for ordinal, (domain, value) in enumerate(metadata.items()):
        row = dict(original_row, gene=domain, gene_order=ordinal, key=value['gene_symbol'],
                   minimum_int=value['aa_start'], maximum=value['aa_end'], auxiliary_int=value['prot_length'])
        rows.append(row)
    result = pa.Table.from_pylist(rows, schema=table.schema).replace_schema_metadata(
        {b'priva.schema': b'priva-alpha-domain-parquet-v2', b'priva.kind': b'metadata'})
    pq.write_table(result, output, compression='zstd')


def check(args):
    root = args.work_root
    root.mkdir(parents=True, exist_ok=False)
    source_dir = Path(__file__).parent
    pins = [identity(path) for path in [Path(__file__), args.binary, args.source, args.lookup_source,
            args.normalization_json, source_dir/'das_domain_lookup_compare.py', source_dir/'alpha_domain_compare.py']]
    helper = load(source_dir/'das_domain_lookup_compare.py', 'source_lookup')
    subject = helper.normalizer(args.lookup_source, json.loads(args.normalization_json.read_text()))
    mapping = {}
    for database, accession, kinds in [('Pfam', 'PF1', ['Family']), ('Pfam', 'PF2', ['Homologous_superfamily']),
                                      ('SMART', 'SM1', ['Domain']), ('Pfam', 'PF3', []),
                                      ('Pfam', 'PF4', ['Domain', 'Family']), ('Pfam', 'PF5', ['Family', 'Domain'])]:
        key = subject.normalize_vep_to_interpro(database, accession)
        if key is None:
            raise RuntimeError('fixture normalization profile changed')
        mapping[key] = [(f'IPR{i}_{accession}', kind, 'short', 'full', [('GO:1', 'molecular_function', 'desc')], 'abstract')
                        for i, kind in enumerate(kinds)]
    with (root/'mapping.pkl').open('xb') as stream:
        pickle.dump(mapping, stream, protocol=4)
    if run([args.binary, 'cache-convert', 'interpro_mapping_pickle', '--pickle', root/'mapping.pkl', '--output', root/'mapping.rkyv'], root, 'mapping'):
        raise RuntimeError('native mapping fixture conversion failed')
    alpha = load(source_dir/'alpha_domain_compare.py', 'alpha_fixture')
    fixture = root/'metadata_template'
    fixture.mkdir()
    vcf = alpha.make_case(fixture, {})
    if run([args.binary, 'domain-metadata', 'build', '--vcf', vcf, '--output', fixture/'metadata.parquet'], fixture, 'native'):
        raise RuntimeError('native metadata template construction failed')
    results = []
    for name, specification in cases().items():
        seeds = [0, 1, 42, 4294967295] if name in ['normal', 'threshold_32', 'repeated_spots'] else [0]
        for seed in seeds:
            folder = root/(name+'_'+str(seed))
            folder.mkdir()
            domains = specification.get('domains', DOMAINS)
            (folder/'domains.jsonl').write_text(''.join(json.dumps(domain)+'\n' for domain in domains))
            start, end = specification.get('bounds', [1, 10])
            metadata = {domain: dict(gene_symbol=specification.get('symbol', 'H' if domain.startswith('H:') else 'G'),
                                     aa_start=start, aa_end=end, prot_length=100)
                        for domain in dict.fromkeys(domains) if not domain.startswith('MISSING:')}
            (folder/'metadata.json').write_text(json.dumps(metadata)+'\n')
            fixture_metadata(fixture/'metadata.parquet', metadata, folder/'metadata.parquet')
            spots = folder/'spots.tsv'
            if specification.get('compression'):
                spots = folder/'spots.tsv.gz'
                payload = gzip.compress(specification['text'].encode(), mtime=0)
                if specification['compression'] == 'truncated':
                    payload = payload[:-5]
                spots.write_bytes(payload)
            elif not specification.get('missing'):
                spots.write_bytes(specification['text'].encode())
            threshold = specification.get('threshold', '0.5')
            common = ['--domains', folder/'domains.jsonl', '--spots', spots, '--normalization-json', args.normalization_json]
            python = [sys.executable, __file__, 'oracle', '--source', args.source, '--lookup-source', args.lookup_source,
                      '--mapping', root/'mapping.pkl', '--metadata', folder/'metadata.json', '--threshold='+threshold, *common]
            native = [args.binary, 'das-domains', 'coldspot', '--metadata', folder/'metadata.parquet',
                      '--interpro', root/'mapping.rkyv', '--python-hash-seed', str(seed), '--overlap-threshold', threshold, *common]
            statuses = [run(python, folder, 'python', dict(os.environ, PYTHONHASHSEED=str(seed))), run(native, folder, 'native')]
            categories = [category(folder/(label+'.stderr.log')) if code else None for label, code in zip(['python', 'native'], statuses)]
            differences = []
            passed = statuses[0] == statuses[1] and (not statuses[0] or categories[0] is not None and categories[0] == categories[1])
            compared = 0
            if statuses == [0, 0]:
                expected = [json.loads(line) for line in (folder/'python.stdout.log').read_text().splitlines()]
                actual = [json.loads(line) for line in (folder/'native.stdout.log').read_text().splitlines()]
                if len(expected) != len(actual):
                    differences.append(dict(field='row_count', python=len(expected), native=len(actual)))
                for ordinal, (left, right) in enumerate(zip(expected, actual)):
                    compared += 1
                    if typed(left) != typed(right):
                        differences.append(dict(ordinal=ordinal, python=left, native=right))
            elif (folder/'native.stdout.log').stat().st_size:
                differences.append(dict(field='partial_native_output'))
            passed = passed and not differences
            (folder/'mismatches.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in differences))
            results.append(dict(case=name, seed=seed, passed=passed, python_exit=statuses[0], native_exit=statuses[1],
                                python_category=categories[0], native_category=categories[1], compared=compared,
                                mismatches=len(differences)))
    for pin in pins:
        if identity(Path(pin['path'])) != pin:
            raise RuntimeError('comparison input/source changed')
    result = dict(passed=all(row['passed'] for row in results), cases=len(results), matched=sum(row['passed'] for row in results),
                  compared=sum(row['compared'] for row in results), defined_failures=sum(row['python_exit'] == 1 for row in results),
                  failures=[row for row in results if not row['passed']], results=results, inputs=pins,
                  full_das_qualified=False, output_order='requested domain order; candidate-set history separate')
    (root/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    args.report.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ['passed', 'cases', 'matched', 'compared', 'defined_failures', 'failures']}), flush=True)
    if not result['passed']:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for name, function in [('check', check), ('oracle', oracle)]:
        command = commands.add_parser(name)
        fields = ['source', 'lookup-source', 'normalization-json']
        fields += ['binary', 'work-root', 'report'] if name == 'check' else ['mapping', 'domains', 'metadata', 'spots']
        for field in fields:
            command.add_argument('--'+field, type=Path, required=True)
        if name == 'oracle':
            command.add_argument('--threshold', type=float, required=True)
        command.set_defaults(function=function)
    args = parser.parse_args()
    args.function(args)


if __name__ == '__main__':
    main()
