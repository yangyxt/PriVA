#!/usr/bin/env python3
"""replay/throughput: complete synthetic pedigree qualification, no family data.

Replay compares a new native candidate to SAVED ACTUAL Python outcomes. It does
not execute either old oracle worker or reuse its disproved expectations.
All candidate fixtures/outputs/logs belong to the supplied fresh data-side root.
"""
import argparse
from collections import Counter
import contextlib
import hashlib
import importlib.util
import json
import logging
import os
from pathlib import Path
import shutil
import resource
import subprocess
import sys
import time
import traceback
from types import SimpleNamespace
import warnings

SOURCE_SHA = '06d243890da53d172ac1bc4e5d18f4c6897d3eba2b6c0a2d0404564705dd0c5b'
SAVED_OUTCOMES = {
    '0': 'e7b550a57c8e02521ba67e93838f51a2401c4056796bd5d8e302b8fd54f11623',
    '1': '83d279cfbba684ae1092fd6b6aad350d3de2c294f989f82387c5ebd07d22a179',
}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def typed(value):
    if value is None or type(value) in (str, bool, int):
        return dict(type=type(value).__name__, value=value)
    if type(value) is float:
        return dict(type='float', repr=repr(value))
    if type(value) in (tuple, list):
        return dict(type=type(value).__name__, items=[typed(item) for item in value])
    return dict(type=type(value).__module__ + '.' + type(value).__qualname__, repr=repr(value))


def untyped(value):
    kind = value['type']
    if kind in ('NoneType', 'str', 'bool', 'int'):
        return value['value']
    if kind in ('tuple', 'list'):
        result = [untyped(item) for item in value['items']]
        return tuple(result) if kind == 'tuple' else result
    raise RuntimeError('Unqualified saved plain type: ' + kind)


def proxy(chrom='1', samples=None, alts=('C', 'G', 'T')):
    return SimpleNamespace(chrom=chrom, samples=samples or {'p': {'GT': (0, 1)}},
                           alts=alts, pos=1, id='unit', ref='A')


class Capture(logging.Handler):
    def __init__(self):
        super().__init__()
        self.rows = []

    def emit(self, record):
        self.rows.append(dict(level=record.levelname, logger=record.name, function=record.funcName,
                             line=record.lineno, message=record.getMessage(), created=record.created))


def modules(candidate):
    sys.path.insert(0, str(candidate / 'scripts'))
    # Match the explicitly recorded imported namespace; canonical CLI is
    # independently executed below and remains a distinct scope.
    spec = importlib.util.spec_from_file_location('frozen_pedigree', candidate / 'scripts/pedigree_filtration_per_fam.py')
    reference = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reference)
    import priva_rust_pedigree_filtration as bridge
    import pysam
    bridge._admit(vars(reference))
    return bridge, reference, pysam


def output(path, pysam, preexisting=None):
    result = dict(exists=path.exists(), regular_file=path.is_file(), records=[])
    if not path.is_file():
        result['unchanged_from_preexisting'] = False
        return result
    result.update(sha256=sha(path), bytes=path.stat().st_size)
    try:
        with pysam.VariantFile(path) as source:
            result['header_text'] = str(source.header)
            for row in source:
                result['records'].append(dict(id=row.id, text=str(row), source=row.info.get('VARIANT_SOURCE'),
                    GT={sample: typed(row.samples[sample].get('GT')) for sample in row.samples}))
    except Exception as error:
        result['readback_exception'] = dict(type=type(error).__name__, message=str(error))
    result['unchanged_from_preexisting'] = preexisting is not None and result['sha256'] == preexisting
    return result


def normalize(value, root):
    if isinstance(value, str):
        return value.replace(str(root), '${SEED_ROOT}')
    if isinstance(value, list):
        return [normalize(item, root) for item in value]
    if isinstance(value, dict):
        return {key: normalize(item, root) for key, item in value.items()}
    return value


def projection(row, root):
    error = None if row['error'] is None else {name: row['error'][name] for name in ('type', 'qualified_type', 'message', 'args')}
    result = dict(id=row['id'], value=row['value'], error=normalize(error, root))
    if 'output' in row:
        result['output'] = {name: row['output'].get(name) for name in ('exists', 'regular_file', 'header_text',
            'records', 'unchanged_from_preexisting')}
    result['logs'] = [normalize({key: item[key] for key in ('level', 'logger', 'function', 'line', 'message')}, root)
                      for item in row['logs']]
    # Physical warning filenames/lines are retained, not falsely presented as
    # identical stacks after porting; category/message/order must be exact.
    result['warnings'] = [{key: item[key] for key in ('category', 'message')} for item in row['warnings']]
    return result


def operation(row, bridge, namespace, client, pysam, header, source_root, target_root, binary):
    name, scope, data = row['id'], row['scope'], row['input']
    def instance(record=None, sex='1', father='f', mother='m', proband='p'):
        return bridge.VariantSource(record if record is not None else proxy(), proband, father, mother, sex,
                                    session=client, namespace=namespace)
    if scope == 'raw-Python-genotype-domain':
        predicate = name.rsplit('-', 1)[1]
        return lambda: bridge._predicate(predicate, untyped(data['genotype']), session=client, namespace=namespace), None
    if scope == 'explicit-getter-proxy':
        record = proxy(samples={'p': {'GT': untyped(data['record_samples']['p']['GT'])}})
        return lambda: instance(record)._get_genotype(untyped(data['requested_sample'])), None
    if scope == 'chromosome-string-proxy':
        return lambda: instance(proxy(name[len('chrom-'):])).chrom_type, None
    if scope == 'malformed-proxy-error':
        if name == 'chrom-none':
            return lambda: instance(proxy(None)), None
        alts = None if name == 'eager-debug-alts-none' else ()
        return lambda: instance(proxy(alts=alts), father=None, mother=None).determine_source(), None
    if scope == 'explicit-mutable-class-proxy':
        def mutated():
            value = instance(father=None, mother=None)
            value.chrom_type = 'EXPLICIT_UNSUPPORTED_PROXY'
            return value.determine_source()
        return mutated, None
    if scope == 'explicit-uncoerced-class-sex':
        sex = 1 if name == 'raw-sex-1' else None
        return lambda: instance(proxy('X'), sex=sex, father=None, mother=None).determine_source(), None
    if scope in ('actual-pysam-parent-grid', 'actual-pysam-proband-grid'):
        def actual():
            record = header.new_record(contig=data['chrom'], start=0, alleles=('A', 'C', 'G', 'T'), id='synthetic')
            for sample, gt in data['actual_pysam_GT'].items():
                record.samples[sample]['GT'] = untyped(gt)
            return instance(record, sex=untyped(data['sex']), father=data['father'], mother=data['mother'], proband=data['proband']).determine_source()
        return actual, None
    if scope == 'actual-pysam-absent-sample':
        def absent():
            record = header.new_record(contig='1', start=0, alleles=('A', 'C', 'G', 'T'), id='synthetic')
            for sample, gt in (('p', (0, 1)), ('f', (0, 1)), ('m', (0, 0))):
                record.samples[sample]['GT'] = gt
            return instance(record, sex='2', proband='absent').determine_source()
        return absent, None
    if scope == 'actual-pysam-missing-GT':
        return lambda: instance(header.new_record(contig='1', start=0, alleles=('A', 'C', 'G', 'T'), id='no-format'), sex='2').determine_source(), None
    if scope == 'direct-source-method-proxy':
        method = data['function'].split('.')[-1]
        return lambda: getattr(instance(), method)(*untyped(data['operands'])), None
    require(scope == 'actual-pandas-pysam-whole-file', 'Unknown saved reference scope: ' + scope)
    case = data['case_spec']
    folder = target_root / 'cases' / case['id']
    folder.mkdir(parents=True)
    original = source_root / 'cases' / case['id']
    for filename in ('input.vcf', 'input.ped'):
        shutil.copyfile(original / filename, folder / filename)
    destination = folder / 'output.vcf'
    if case.get('preexisting'):
        destination.write_bytes(b'preexisting synthetic sentinel\n')
    if case.get('output_directory'):
        destination.mkdir()
    before = sha(destination) if destination.is_file() else None
    input_path = folder / ('absent.vcf' if case.get('input_missing') else 'input.vcf')
    statistics = {}
    def whole():
        return bridge.pedigree_filter(str(input_path), str(folder / 'input.ped'), untyped(data['target']), str(destination),
                                      binary=binary, namespace=namespace, statistics=statistics)
    return whole, dict(path=destination, before=before, statistics=statistics)


def replay(args):
    require(os.environ.get('PYTHONHASHSEED') == args.seed, 'Use actual fresh process hash seed, not a runtime assignment')
    require(not args.output.exists(), 'Preserve all candidate attempts')
    seed_root = args.saved_oracle / ('seed-' + args.seed)
    require(sha(seed_root / 'result.json') == SAVED_OUTCOMES[args.seed], 'Saved actual result pin changed')
    saved_result = json.loads((seed_root / 'result.json').read_bytes())
    require(sha(saved_result['outcome_pin']['path']) == saved_result['outcome_pin']['sha256'], 'Saved actual outcomes changed')
    saved = json.loads(Path(saved_result['outcome_pin']['path']).read_bytes())
    require(len(saved) == 652 and len({row['id'] for row in saved}) == 652, 'Saved exact denominator changed')
    args.output.mkdir(parents=True)
    bridge, reference, pysam = modules(args.candidate)
    namespace = vars(reference)
    capture = Capture()
    reference.logger.addHandler(capture)
    header = pysam.VariantHeader()
    for name in ('1', 'X', 'Y'):
        header.contigs.add(name, length=1000)
    header.formats.add('GT', 1, 'String', 'Genotype')
    for name in ('p', 'f', 'm'):
        header.add_sample(name)
    observed, differences = [], []
    actual_original_calls = Counter()
    per_case_calls = Counter()
    reference_codes = {value.__code__: name for name, value in namespace.items() if callable(value) and hasattr(value, '__code__')}
    reference_codes.update({method.__code__: 'VariantSource.' + name for name, method in vars(reference.VariantSource).items()
                            if hasattr(method, '__code__')})
    def profile(frame, event, _arg):
        if event == 'call' and frame.f_code in reference_codes:
            name = reference_codes[frame.f_code]
            actual_original_calls[name] += 1
            per_case_calls[name] += 1
    with bridge.NativeSession(args.binary) as client:
        for prior in saved:
            begin, requests = len(capture.rows), client.sequence
            completed_before = Counter(client.completed_operations)
            per_case_calls.clear()
            function, details = operation(prior, bridge, namespace, client, pysam, header, seed_root, args.output, args.binary)
            with warnings.catch_warnings(record=True) as recorded:
                sys.setprofile(profile)
                try:
                    value, error = typed(function()), None
                except Exception as failure:
                    value = None
                    error = dict(type=type(failure).__name__, qualified_type=type(failure).__module__ + '.' + type(failure).__qualname__,
                        message=str(failure), args=typed(failure.args), traceback=[dict(filename=row.filename, line=row.lineno, function=row.name)
                            for row in traceback.extract_tb(failure.__traceback__)])
                finally:
                    sys.setprofile(None)
            row = dict(id=prior['id'], scope=prior['scope'], value=value, error=error, logs=capture.rows[begin:],
                warnings=[dict(category=item.category.__name__, message=str(item.message), filename=item.filename, line=item.lineno) for item in recorded],
                native_requests=client.sequence - requests, actual_reference_body_calls=dict(per_case_calls))
            row['completed_native_operations'] = dict(Counter(client.completed_operations) - completed_before)
            if details is not None:
                row['output'] = output(details['path'], pysam, details['before'])
                row['whole_native_statistics'] = details['statistics']
                stats = details['statistics']
                native_errors = []
                if stats.get('backend') != 'rust' or stats.get('native_processes') != 1 or stats.get('reference_fallback') is not False:
                    native_errors.append('admitted whole-file case did not use exactly one native process without fallback')
                completed_ops = stats.get('completed_native_operations', {})
                if completed_ops.get('hello') != 1 or completed_ops.get('quit') != 1:
                    native_errors.append('native handshake/closed session proof missing')
                if stats.get('observed_native_operations') != completed_ops:
                    native_errors.append('a native request lacked a validated native response')
                if stats.get('written_records') != len(row['output']['records']):
                    native_errors.append('actual written prefix count differs from parsed output')
                if per_case_calls.get('pedigree_filter', 0):
                    native_errors.append('original complete Python filter was called')
                row['genuine_native_errors'] = native_errors
                if native_errors:
                    differences.append(dict(id=row['id'], native_execution_errors=native_errors, actual_statistics=stats))
            expected = projection(prior, seed_root)
            actual = projection(row, args.output)
            if expected != actual:
                differences.append(dict(id=row['id'], expected=expected, actual=actual))
            observed.append(row)
    bridge.close_sessions()
    cli = []
    for prior in saved_result['CLI_controls']:
        name = prior['id']
        folder = args.output / 'cli' / name
        folder.mkdir(parents=True)
        destination = folder / 'output.vcf'
        argv = [sys.executable, '-B', str(args.candidate / 'scripts/pedigree_filtration_per_fam.py'), '--pedigree-binary', str(args.binary)]
        if name == 'missing-cli-arguments':
            argv += ['-v', 'not-opened-without-required-arguments']
        else:
            source = args.output / 'cases' / name
            argv += ['-v', str(source / ('absent.vcf' if name == 'unknown-family-before-input-open' else 'input.vcf')),
                     '-p', str(source / 'input.ped'), '-f', 'absent' if name == 'unknown-family-before-input-open' else 'F', '-o', str(destination)]
        with (folder / 'stdout.log').open('xb') as stdout, (folder / 'stderr.log').open('xb') as stderr:
            done = subprocess.run(argv, stdout=stdout, stderr=stderr, timeout=60)
        actual_output = output(destination, pysam)
        fields = ('exists', 'regular_file', 'header_text', 'records')
        mismatch = done.returncode != prior['returncode'] or {key: actual_output.get(key) for key in fields} != {key: prior['output'].get(key) for key in fields}
        cli.append(dict(id=name, actual_argv=argv, returncode=done.returncode, saved_returncode=prior['returncode'], output=actual_output,
                        difference=mismatch, stdout_sha256=sha(folder / 'stdout.log'), stderr_sha256=sha(folder / 'stderr.log'),
                        stderr_scope='actual new optional-argument usage/adapter stack retained; not claimed byte-identical to old script stack'))
        if mismatch:
            differences.append(dict(id='CLI-' + name, actual=cli[-1], expected=prior))
    save(args.output / 'outcomes.json', observed)
    save(args.output / 'differences.json', differences)
    save(args.output / 'CLI.json', cli)
    result = dict(kind='pedigree-complete-saved-oracle-replay-v1', completed=not differences,
        actual_native_step=os.environ['SLURM_JOB_ID'] + '.' + os.environ['SLURM_STEP_ID'], seed=args.seed,
        saved_result_sha256=sha(seed_root / 'result.json'), saved_actual_outcomes_sha256=saved_result['outcome_pin']['sha256'],
        original_expectations_used=False, original_seed_worker_rerun=False,
        captured_calls=len(observed), actual_direct_CLI_controls=len(cli), differences=len(differences),
        actual_reference_body_calls=dict(actual_original_calls),
        reference_body_semantics='Actual getter/log/error/library boundaries and explicit unsupported fallback calls are counted, not relabelled native.',
        outcomes_sha256=sha(args.output / 'outcomes.json'), differences_sha256=sha(args.output / 'differences.json'),
        no_family_model_or_accuracy_acceptance=True)
    save(args.output / 'result.json', result)
    print(json.dumps(dict(completed=result['completed'], calls=len(observed), CLI=len(cli), differences=len(differences),
                         result_sha256=sha(args.output / 'result.json')), sort_keys=True))
    require(not differences, 'Preserve every candidate mismatch')


def measure(args):
    """One fresh synthetic file, full operation, no saved-oracle rerun."""
    bridge, reference, pysam = modules(args.candidate)
    reference.logger.handlers[:] = [logging.NullHandler()]
    reference.logger.propagate = False
    reference.logger.setLevel(logging.INFO)
    stats = {}
    original_calls = Counter()
    original_code = reference.pedigree_filter.__code__
    require(sys.version_info[:2] == (3, 12), 'Qualified code-local observer requires Python3.12')
    monitor = sys.monitoring
    tool_id = next((value for value in range(6) if monitor.get_tool(value) is None), None)
    require(tool_id is not None, 'Preserve occupied monitoring tools')
    def source_entry(code, *_position):
        if code is original_code:
            original_calls['pedigree_filter'] += 1
    monitor.use_tool_id(tool_id, 'pedigree-throughput-source-entry')
    monitor.register_callback(tool_id, monitor.events.PY_START, source_entry)
    monitor.set_local_events(tool_id, original_code, monitor.events.PY_START)
    started = time.perf_counter()
    before_self, before_child = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    try:
        if args.backend == 'python':
            reference.pedigree_filter(str(args.vcf), str(args.ped), 'F', str(args.output))
        else:
            bridge.pedigree_filter(str(args.vcf), str(args.ped), 'F', str(args.output),
                                   binary=args.binary, namespace=vars(reference), statistics=stats)
    finally:
        monitor.set_local_events(tool_id, original_code, 0)
        monitor.register_callback(tool_id, monitor.events.PY_START, None)
        monitor.free_tool_id(tool_id)
    elapsed = time.perf_counter() - started
    after_self, after_child = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    with pysam.VariantFile(args.output) as source:
        written = sum(1 for _ in source)
    result = dict(kind='pedigree-synthetic-throughput-measurement-v1', backend=args.backend,
        input_rows=args.rows, written_rows=written, full_operation_wall_seconds=elapsed,
        process_cpu_seconds=(after_self.ru_utime + after_self.ru_stime - before_self.ru_utime - before_self.ru_stime),
        children_cpu_seconds=(after_child.ru_utime + after_child.ru_stime - before_child.ru_utime - before_child.ru_stime),
        Python_process_peak_RSS_KiB=after_self.ru_maxrss, native_child_peak_RSS_KiB=after_child.ru_maxrss,
        RSS_scope='Separate rusage maxima; not a simultaneous aggregate memory measurement.',
        rows_per_second=args.rows / elapsed, actual_original_filter_calls=dict(original_calls),
        statistics=stats, output_sha256=sha(args.output), output_bytes=args.output.stat().st_size,
        observer_scope='Python3.12 code-local PY_START on the original whole-filter code only; no broad frame profiler.',
        timing_scope='Imports excluded; full call includes source/ELF admission, session start/close and all PED/VCF I/O.')
    save(args.output.with_suffix('.measurement.json'), result)


def throughput(args):
    """Compare complete ordered outputs and measure small/large synthetic files."""
    require(not args.output.exists(), 'Preserve throughput attempts')
    args.output.mkdir(parents=True)
    groups = []
    # Fresh deterministic fixtures are not the saved seeds or family data.
    rows = [('1', '0/1', '0/0', '0/0'), ('1', '0/2', '0/0', '0/0'),
            ('1', '0/0', '0/0', '0/0'), ('1', '1/1', '1/1', '0/0'),
            ('X', '0/1', '0/1', '0/0'), ('Y', '0/1', '0/1', '0/0'),
            ('1', './1', '0/0', '0/0'), ('1', '0/1', './.', '0/1'),
            ('1', '0/3', '0/0', '0/0'), ('1', '0/1', '0/1', '0/1')]
    for count in args.count:
        require(count > 0, 'Positive synthetic input denominator required')
        folder = args.output / str(count)
        folder.mkdir()
        ped = folder / 'input.ped'
        ped.write_text('F\tp\tf\tm\t2\t2\nF\tf\t0\t0\t1\t1\nF\tm\t0\t0\t2\t1\n')
        vcf = folder / 'input.vcf'
        with vcf.open('x') as stream:
            stream.write('##fileformat=VCFv4.2\n')
            for chrom in ('1', 'X', 'Y'):
                stream.write(f'##contig=<ID={chrom},length={count + 10}>\n')
            stream.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tp\tf\tm\n')
            for index in range(count):
                chrom, p, f, m = rows[index % len(rows)]
                stream.write(f'{chrom}\t{index + 1}\tbench-{index}\tA\tC,G,T\t.\tPASS\t.\tGT\t{p}\t{f}\t{m}\n')
        measurements = {}
        for backend in ('python', 'rust'):
            destination = folder / (backend + '.vcf')
            argv = [sys.executable, '-B', str(Path(__file__)), 'measure', '--candidate', str(args.candidate),
                    '--binary', str(args.binary), '--vcf', str(vcf), '--ped', str(ped),
                    '--rows', str(count), '--output', str(destination), '--backend', backend]
            began = time.perf_counter()
            with (folder / (backend + '.log')).open('x') as stream:
                done = subprocess.run(argv, stdout=stream, stderr=subprocess.STDOUT, timeout=600)
            require(done.returncode == 0, 'Retain failed throughput child: ' + backend)
            report = json.loads(destination.with_suffix('.measurement.json').read_bytes())
            report['fresh_process_wall_seconds_including_imports'] = time.perf_counter() - began
            measurements[backend] = report
        left, right = measurements['python'], measurements['rust']
        stats = right['statistics']
        checks = dict(byte_exact_ordered_output=(left['output_sha256'], left['output_bytes'], left['written_rows']) ==
                      (right['output_sha256'], right['output_bytes'], right['written_rows']),
            one_native_process=stats.get('native_processes') == 1,
            native_requested_and_completed=stats.get('observed_native_operations') == stats.get('completed_native_operations'),
            exact_native_record_denominator=stats.get('completed_native_operations', {}).get('filter-begin') == count,
            parsed_and_written_denominators=stats.get('parsed_records') == count and stats.get('written_records') == right['written_rows'],
            no_original_filter_fallback=stats.get('reference_fallback') is False and not right['actual_original_filter_calls'],
            original_called_once=left['actual_original_filter_calls'] == {'pedigree_filter': 1})
        groups.append(dict(input_rows=count, checks=checks, measurements=measurements,
                           native_over_python_wall_ratio=right['full_operation_wall_seconds'] / left['full_operation_wall_seconds']))
        save(folder / 'comparison.json', groups[-1])
    result = dict(kind='pedigree-synthetic-throughput-v1', completed=all(all(item['checks'].values()) for item in groups),
        native_step=os.environ['SLURM_JOB_ID'] + '.' + os.environ['SLURM_STEP_ID'], groups=groups,
        saved_oracle_worker_rerun=False, qualification_scope='Synthetic complete ordered-file parity and measured cost; no family/biological or general speedup claim.',
        measurement_limit='Single sequential Python-then-Rust pair per size; report any slowdown explicitly.')
    save(args.output / 'result.json', result)
    print(json.dumps(dict(completed=result['completed'], measurements=[dict(rows=item['input_rows'],
        native_over_python_wall_ratio=item['native_over_python_wall_ratio']) for item in groups],
        result_sha256=sha(args.output / 'result.json')), sort_keys=True))
    require(result['completed'], 'Retain exact throughput/check mismatch')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='procedure', required=True)
    compare = commands.add_parser('replay')
    compare.add_argument('--candidate', type=Path, required=True)
    compare.add_argument('--saved-oracle', type=Path, required=True)
    compare.add_argument('--binary', type=Path, required=True)
    compare.add_argument('--output', type=Path, required=True)
    compare.add_argument('--seed', choices=('0', '1'), required=True)
    compare.set_defaults(function=replay)
    benchmark = commands.add_parser('throughput')
    for name in ('candidate', 'binary', 'output'):
        benchmark.add_argument('--' + name, type=Path, required=True)
    benchmark.add_argument('--count', type=int, action='append', required=True)
    benchmark.set_defaults(function=throughput)
    worker = commands.add_parser('measure')
    for name in ('candidate', 'binary', 'vcf', 'ped', 'output'):
        worker.add_argument('--' + name, type=Path, required=True)
    worker.add_argument('--rows', type=int, required=True)
    worker.add_argument('--backend', choices=('python', 'rust'), required=True)
    worker.set_defaults(function=measure)
    args = parser.parse_args()
    args.function(args)


if __name__ == '__main__':
    main()
