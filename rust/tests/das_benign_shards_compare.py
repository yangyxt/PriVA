#!/usr/bin/env python3
"""Qualify explicit indexed benign counts and complete ordered aggregation.

check --binary PRIVA --source ORIGINAL.py --fixture-helper AA_FIXTURES.py
      --work-root NEW_DIRECTORY --report JSON
oracle executes the unchanged original function. With --chrom it only replaces
its VariantFile iterator with that file's pysam.fetch(chrom) iterator. Successful
complete shard sums are compared independently with the ordinary sequential
function. Indexed parse failures remain separately identified from sequential
failures; they do not establish unrestricted replacement of that original route.
"""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

sys.dont_write_bytecode = True


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def base():
    return load(Path(__file__).with_name('das_benign_compare.py'), 'benign_source_oracle')


def oracle(args):
    source = base()
    function = source.original(args.source)
    counts = dict(records=0)
    if args.chrom is not None or args.count_output is not None:
        import pysam
        factory = pysam.VariantFile
        class Chromosome:
            def __init__(self, *values, **options):
                self.file = factory(*values, **options)
                self.header = self.file.header
            def __iter__(self):
                iterator = self.file if args.chrom is None else self.file.fetch(args.chrom)
                for record in iterator:
                    counts['records'] += 1
                    yield record
            def close(self):
                self.file.close()
        function.__globals__['pysam'] = SimpleNamespace(VariantFile=Chromosome)
    domains = [json.loads(line) for line in args.domains.read_text().splitlines()]
    _, _, summary = function(set(domains), str(args.vcf), args.minimum)
    for domain in domains:
        print(json.dumps(dict(domain=domain, **summary[domain])), flush=True)
    if args.count_output is not None:
        with args.count_output.open('x') as output:json.dump(counts,output)


def check(args):
    import pysam
    common = base()
    root = args.work_root
    root.mkdir(parents=True, exist_ok=False)
    pins = [common.identity(path) for path in [Path(__file__), Path(__file__).with_name('das_benign_compare.py'),
        args.binary, args.source, args.fixture_helper]]
    inputs = load(args.fixture_helper, 'benign_shard_vcf_inputs')
    env = dict(__import__('os').environ, PYTHONDONTWRITEBYTECODE='1')
    results, guards, route_boundaries = [], [], []
    def run(folder, label, command):
        with (folder/(label+'.stdout.log')).open('x') as out, (folder/(label+'.stderr.log')).open('x') as err:
            status = subprocess.run(list(map(str, command)), env=env, stdout=out, stderr=err, timeout=120).returncode
        return status, common.category(folder/(label+'.stderr.log')) if status else None
    def python_command(vcf, domains, minimum, chrom=None):
        command = [sys.executable, __file__, 'oracle', '--source', args.source, '--vcf', vcf,
                   '--domains', domains, '--minimum', str(minimum)]
        if chrom is not None:command += ['--chrom', chrom]
        return command
    def read(path):
        return [json.loads(line) for line in path.read_text().splitlines()]
    specifications = {name:copy.deepcopy(common.cases()[name]) for name in [
        'normal','empty','three_occurrences','repeated_annotations','repeated_prefixes','unicode',
        'duplicate_info_first','short_annotation','short_nontx','valid_then_missing','no_format',
        'other_last_format','negative_pos','bad_pos','empty_pos','invalid_unused_gt','invalid_unused_format','csq_1_Integer',
        'short_record_5','short_record_6','short_record_7',
        'sig_._Integer','empty_candidates','minimum_negative','minimum_huge','minimum_negative_huge',
        'minimum_unicode','minimum_underscores']}
    specifications.update(
        csi=dict(csi=True),
        header_order=dict(contigs=['chr2','chrEmpty','chr1'],rows=[common.row(),common.row(pos=11),common.row(chrom='chr2')]),
        exact_threshold_across_chromosomes=dict(contigs=['chr2','chrEmpty','chr1'],rows=[common.row(),common.row(pos=11),common.row(chrom='chr2')]),
        below_threshold=dict(contigs=['chr2','chrEmpty','chr1'],rows=[common.row(),common.row(chrom='chr2')]),
        all_empty=dict(contigs=['chr2','chrEmpty','chr1'],rows=[]),
        repeated_domains_across_chromosomes=dict(contigs=['chr2','chrEmpty','chr1'],
            rows=[common.row(csq=common.annotation(DOMAINS='DB:A:B&DB:A')),common.row(chrom='chr2',csq=common.annotation()+','+common.annotation())]))
    pipeline = None
    for name, spec in specifications.items():
        folder = root/name
        folder.mkdir()
        spec.setdefault('fields',common.FIELDS)
        spec.setdefault('rows',[common.row(csq=common.annotation(spec['fields']))])
        vcf = inputs.make_case(folder,spec)
        requested = spec.get('domains',common.DOMAINS)
        domains = folder/'domains.jsonl'
        domains.write_text(''.join(json.dumps(value)+'\n' for value in requested))
        minimum = spec.get('minimum',3)
        whole_python = run(folder,'whole_python',python_command(vcf,domains,minimum))
        whole_native = run(folder,'whole_native',[args.binary,'das-domains','benign','--domains',domains,'--vcf',vcf,'--min-benign',str(minimum)])
        whole_match = whole_python == whole_native and (whole_python[0] == 0 or whole_python[1] is not None)
        if whole_python[0] == 0:
            expected = read(folder/'whole_python.stdout.log')
            whole_match &= common.typed(expected) == common.typed(read(folder/'whole_native.stdout.log'))
        shards, indexed, totals = [], [], dict.fromkeys(requested,0)
        with pysam.VariantFile(str(vcf)) as handle:contigs = list(handle.header.contigs)
        for ordinal, chrom in enumerate(contigs):
            shard = folder/f'counts_{ordinal}.json'
            count_path = folder/f'python_{ordinal}.records.json'
            python = run(folder,f'python_{ordinal}',python_command(vcf,domains,minimum,chrom)+['--count-output',count_path])
            native = run(folder,f'native_{ordinal}',[args.binary,'das-domains','benign-chromosome',
                '--domains',domains,'--vcf',vcf,'--chrom',chrom,'--output',shard])
            matched = python == native and (python[0] == 0 or python[1] is not None)
            if python[0] == native[0] == 0:
                truth = read(folder/f'python_{ordinal}.stdout.log')
                payload = json.loads(shard.read_text())
                count_rows = payload['counts']
                want = [[domain,next(row['n_benign_missense'] for row in truth if row['domain']==domain)] for domain in dict.fromkeys(requested)]
                want = [row for row in want if row[1] > 0]
                matched &= common.typed(count_rows) == common.typed(want)
                matched &= payload['contigs'] == contigs and payload['chromosome'] == chrom and payload['domain_rows'] == len(requested)
                matched &= payload['records'] == json.loads(count_path.read_text())['records']
                for domain, count in count_rows:totals[domain] += count
                shards.append(shard)
            else:
                matched &= not shard.exists() and (folder/f'native_{ordinal}.stdout.log').stat().st_size == 0
                if whole_python != python:
                    route_boundaries.append(dict(case=name,chromosome=chrom,sequential=list(whole_python),indexed=list(python),
                        claim='indexed procedure has a distinct parse boundary; no general sequential replacement claim'))
            indexed.append(dict(chromosome=chrom,python=list(python),native=list(native),passed=bool(matched)))
        merged = None
        sequential_sum_equal = None
        if len(shards) == len(contigs) and whole_python[0] == 0:
            merged = run(folder,'merged',[args.binary,'das-domains','benign-merge','--domains',domains,
                '--min-benign',str(minimum),*[word for path in shards for word in ['--input',path]]])
            summed = [dict(domain=domain,n_benign_missense=totals[domain],removed=totals[domain]>=int(minimum)) for domain in requested]
            sequential_sum_equal = common.typed(summed) == common.typed(expected)
            merged_match = merged == (0,None) and common.typed(read(folder/'merged.stdout.log')) == common.typed(expected)
            if name == 'header_order':pipeline=(folder,domains,shards)
        else:merged_match=True
        results.append(dict(case=name,passed=bool(whole_match and all(row['passed'] for row in indexed) and merged_match and sequential_sum_equal is not False),
            sequential_python=list(whole_python),sequential_native=list(whole_native),indexed=indexed,
            successful_whole_function_aggregation=sequential_sum_equal,merge=list(merged) if merged else None))
    folder,domains,shards=pipeline
    def guard(name, selected=None, payload=None, expected=1):
        case=root/('guard_'+name);case.mkdir()
        paths=list(shards if selected is None else selected)
        if payload is not None:
            changed=case/'mutated.json';changed.write_text(json.dumps(payload)+'\n');paths[0]=changed
        status=run(case,'native',[args.binary,'das-domains','benign-merge','--domains',domains,'--min-benign','3',
            *[word for path in paths for word in ['--input',path]]])
        guards.append(dict(case=name,status=list(status),passed=status==(expected,'ValueError') and (case/'native.stdout.log').stat().st_size==0))
    guard('missing',shards[:-1]);guard('reversed',list(reversed(shards)));guard('duplicate',[shards[0],shards[0],shards[2]])
    original=json.loads(shards[0].read_text())
    for name,key,value in [('format','format','wrong'),('source','source_sha256','0'*64),('domain_hash','domains_sha256','0'*64),
        ('domain_rows','domain_rows',0),('chromosome','chromosome','chrX'),('catalogue','contigs',['chr2','chr1']),
        ('zero_count','counts',[[common.DOMAINS[0],0]]),('negative_count','counts',[[common.DOMAINS[0],-1]]),
        ('float_count','counts',[[common.DOMAINS[0],1.0]]),('duplicate_count','counts',[[common.DOMAINS[0],1],[common.DOMAINS[0],1]]),
        ('unrequested','counts',[['not-requested',1]])]:
        changed=copy.deepcopy(original);changed[key]=value;guard(name,payload=changed)
    case=root/'guard_existing_output';case.mkdir();destination=case/'counts.json';destination.write_bytes(b'previous output\n')
    status=run(case,'native',[args.binary,'das-domains','benign-chromosome','--domains',domains,'--vcf',folder/'input.vcf.gz',
        '--chrom','chr2','--output',destination])
    guards.append(dict(case='existing_output',status=list(status),passed=status==(1,'FileExistsError') and destination.read_bytes()==b'previous output\n'))
    for pin in pins:
        if common.identity(Path(pin['path']))!=pin:raise RuntimeError('source/input changed during comparison')
    result=dict(passed=all(row['passed'] for row in results+guards),cases=len(results),indexed_chromosomes=sum(len(row['indexed']) for row in results),
        successful_whole_function_aggregations=sum(row['successful_whole_function_aggregation'] is True for row in results),
        guards=guards,guard_cases=len(guards),results=results,route_boundaries=route_boundaries,
        failures=[row for row in results+guards if not row['passed']],inputs=pins,
        full_das_qualified=False,unrestricted_sequential_replacement_qualified=False)
    (root/'report.json').write_text(json.dumps(result,indent=2)+'\n');args.report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['passed','cases','indexed_chromosomes','successful_whole_function_aggregations','guard_cases','route_boundaries','failures']}),flush=True)
    if not result['passed']:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest='command',required=True)
    command=commands.add_parser('check')
    for name in ['binary','source','fixture-helper','work-root','report']:command.add_argument('--'+name,type=Path,required=True)
    command.set_defaults(function=check)
    command=commands.add_parser('oracle')
    for name in ['source','vcf','domains']:command.add_argument('--'+name,type=Path,required=True)
    command.add_argument('--minimum',type=int,required=True);command.add_argument('--chrom')
    command.add_argument('--count-output',type=Path);command.set_defaults(function=oracle)
    args=parser.parse_args();args.function(args)

if __name__=='__main__':main()
