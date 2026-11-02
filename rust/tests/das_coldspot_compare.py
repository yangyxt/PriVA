#!/usr/bin/env python3
"""Compare the DAS inclusive coldspot loop with its unchanged Python source.

compare --probe BINARY --work-root NEW_DIRECTORY
Uses actual pandas rows but injects the parsed table; file parsing/CLI admission
and the later complete DAS output remain separate qualification gates.
"""
import argparse
import ast
from collections import defaultdict
import itertools
import json
import math
from pathlib import Path
import random
import struct
import subprocess
import pandas as pd


def bits(value):
    return struct.pack('>d', value).hex()


class Logger:
    def __init__(self): self.rows = []
    def info(self, value): self.rows.append(value)


class Lookup:
    def __init__(self, kind): self.kind = kind
    def query_interpro_entry_vep_anno(self, *unused):
        return None if self.kind is None else {'interpro_entries': [('id', self.kind)]}


def cases():
    thresholds = [-float('inf'), -1., -0., 0., .5, math.nextafter(.5, 1.), .6, 1., float('inf'), float('nan')]
    intervals = [(0, 0), (1, 10), (6, 15), (10, 1), (-5, 4), (-(2**63), 2**63-1),
                 (2**63-10, 2**63-1), (-(2**63), -(2**63)+9)]
    for kind, meta, present, domain, spot, threshold in itertools.product(range(5), [True, False], [True, False], intervals, intervals, thresholds):
        yield dict(kind=kind, meta=meta, present=present, domain=domain, spots=[spot], threshold=threshold)
    rng = random.Random(1693)
    for _ in range(12_000):
        anchor = rng.randrange(-(2**63), 2**63)
        coordinates = [max(-(2**63), min(2**63-1, anchor+rng.randrange(-1000, 1001))) for _ in range(14)]
        yield dict(kind=rng.randrange(5), meta=rng.randrange(5)!=0, present=rng.randrange(4)!=0,
                   domain=coordinates[:2], spots=list(zip(coordinates[2::2], coordinates[3::2])),
                   threshold=rng.choice(thresholds) if rng.randrange(2) else rng.random())
    for threshold in thresholds:
        yield dict(kind=1, meta=True, present=False, domain=(1, 10), spots=[], threshold=threshold)


def compare(args):
    root = args.work_root; root.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parents[2] / 'scripts/am_pick_intolerant_domains.py'
    body = [node for node in ast.parse(source.read_text()).body if isinstance(node, ast.FunctionDef) and node.name == 'filter_coldspot_overlap']
    if len(body) != 1: raise RuntimeError('original coldspot procedure absent')
    logger = Logger(); scope = dict(pd=pd, defaultdict=defaultdict, logger=logger)
    exec(compile(ast.fix_missing_locations(ast.Module(body=body, type_ignores=[])), str(source), 'exec'), scope)
    inputs = list(cases()); expected = []; path = 'gene:db:domain'; names = [None, 'Family', 'Homologous_superfamily', 'Domain', 'family']
    original_reader = pd.read_csv
    try:
        with (root/'requests.txt').open('x') as incoming:
            for row in inputs:
                rows = [{'type':'coldspot','gene':'G' if row['present'] else 'OTHER','aa_start_pos':a,'aa_end_pos':b} for a,b in row['spots']]
                # A non-coldspot row must never enter the interval loop.
                rows.append(dict(type='hotspot',gene='G',aa_start_pos=-(2**63),aa_end_pos=2**63-1))
                table = pd.DataFrame(rows)
                pd.read_csv = lambda *args, **kwargs: table
                logger.rows.clear()
                meta = {path:dict(gene_symbol='G',aa_start=row['domain'][0],aa_end=row['domain'][1])} if row['meta'] else {}
                filtered, excluded, summary = scope['filter_coldspot_overlap']({path},meta,'unused.tsv',Lookup(names[row['kind']]),{},row['threshold'])
                result = summary[path]
                evaluated = int(next(value.rsplit(': ',1)[1] for value in logger.rows if 'Domains evaluated (' in value))
                exempt = int(next(value.rsplit(': ',1)[1] for value in logger.rows if 'Exempt (not ' in value))
                if (path in excluded) != result['removed'] or (path in filtered) == result['removed']:
                    raise AssertionError('original summary/set disagreement')
                expected.append([bits(result['coldspot_max_frac_of_coldspot']),bits(result['coldspot_max_frac_of_domain']),
                                 str(int(result['removed'])),str(evaluated),str(exempt)])
                words = [str(row['kind']),str(int(row['meta'])),*map(str,row['domain']),str(int(row['present'])),str(len(row['spots'])),bits(row['threshold'])]
                for a,b in row['spots']: words += [str(a),str(b)]
                incoming.write(' '.join(words)+'\n')
    finally:
        pd.read_csv = original_reader
    with (root/'requests.txt').open('rb') as incoming,(root/'native.txt').open('xb') as out,(root/'native.stderr.log').open('xb') as err:
        status = subprocess.run([str(args.probe),'batch'],stdin=incoming,stdout=out,stderr=err,timeout=180).returncode
    failures = []
    with (root/'native.txt').open() as actual:
        for ordinal,(wanted,line) in enumerate(itertools.zip_longest(expected,actual)):
            if line is None or wanted is None or line.split()!=wanted:
                failures.append(dict(ordinal=ordinal,input=inputs[ordinal] if ordinal<len(inputs) else None,python=wanted,native=line))
    (root/'mismatches.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in failures))
    result = dict(passed=status==0 and not failures,cases=len(inputs),mismatches=len(failures),returncode=status,
                  exempt=sum(row[4]=='1' for row in expected),evaluated=sum(row[3]=='1' for row in expected),
                  removed=sum(row[2]=='1' for row in expected),pandas=pd.__version__,
                  file_parser_qualified=False,full_das_qualified=False)
    (root/'report.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result),flush=True)
    if not result['passed']: raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__); subcommands=parser.add_subparsers(dest='command',required=True)
    command=subcommands.add_parser('compare')
    for name in ['probe','work-root']: command.add_argument('--'+name,type=Path,required=True)
    command.set_defaults(function=compare); args=parser.parse_args(); args.function(args)


if __name__=='__main__': main()
