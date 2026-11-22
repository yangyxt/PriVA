#!/usr/bin/env python3
"""Exact typed/float-bit/error comparison with the previously captured Python oracle.

compare --binary FILE --binary-sha256 SHA --oracle DIR
        --oracle-manifest-sha256 SHA --output NEW_DIR
Every mismatch, expected failure, subprocess stream and transport guard is retained.
Unsupported Python diagnostic objects are reported, never counted as native parity.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(args):
    assert sha(args.binary) == args.binary_sha256
    assert sha(args.oracle / 'manifest.json') == args.oracle_manifest_sha256
    oracle = json.loads((args.oracle / 'manifest.json').read_text())
    for name, wanted in oracle['files'].items(): assert sha(args.oracle / name) == wanted, name
    args.output.mkdir(parents=True, exist_ok=False)
    results, mismatches = [], []

    def invoke(name, source, output, *, source_error=None, guard=False, extra=None):
        before = sha(source)
        previous = output.read_bytes() if output.exists() else None
        argv = [str(args.binary), 'coerce-pext', '--input', str(source), '--output', str(output), *(extra or [])]
        process = subprocess.run(argv, capture_output=True)
        (args.output / (name + '.stdout.log')).write_bytes(process.stdout)
        (args.output / (name + '.stderr.log')).write_bytes(process.stderr)
        assert sha(source) == before and sha(args.binary) == args.binary_sha256
        item = dict(name=name, argv=argv, exit_code=process.returncode, input_unchanged=True,
                    stdout_sha256=hashlib.sha256(process.stdout).hexdigest(), stderr_sha256=hashlib.sha256(process.stderr).hexdigest())
        if source_error or guard:
            item['expected_source_error'] = source_error
            item['transport_guard_only'] = guard
            valid = process.returncode == 2 and (not output.exists() if previous is None else output.read_bytes() == previous)
            if source_error:
                valid = valid and process.stderr == ('priva: ' + source_error['error_type'] + ': ' + source_error['error_message'] + '\n').encode()
        else:
            valid = process.returncode == 0 and output.is_file() and not process.stderr
        item['passed'] = valid
        if not valid: mismatches.append(item)
        results.append(item)
        return process

    result = args.output / 'values.actual.csv'
    invoke('supported_values', args.oracle / 'values.input.csv', result)
    if result.exists():
        with result.open(newline='') as stream: actual = list(csv.reader(stream))
        with (args.oracle / 'values.expected.csv').open(newline='') as stream: expected = list(csv.reader(stream))
        if actual != expected:
            mismatches.append(dict(name='exact_rows_and_float_bits', expected_rows=len(expected), actual_rows=len(actual),
                                   first_differences=[dict(row=i, expected=e, actual=a) for i, (e, a) in enumerate(zip(expected, actual)) if e != a][:50]))
    for i, failure in enumerate(oracle['source_failures']):
        invoke(f'source_failure_{i}', args.oracle / failure['input'], args.output / f'failure_{i}.csv', source_error=failure)
    for kind, value in [('bytes', '00'), ('array', '[]'), ('complex', '1+2j'), ('float16', '3800'),
                        ('longdouble', '0.1'), ('bigint', '1_000'), ('str', 'ff'), ('numpy_bool', '2')]:
        name = 'guard_' + kind
        source = args.output / (name + '.input.csv')
        with source.open('x', newline='') as stream:
            writer = csv.writer(stream); writer.writerow(['index_kind', 'index_value', 'value_kind', 'value'])
            writer.writerow(['int', '1', kind, value])
        invoke(name, source, args.output / (name + '.csv'), guard=True)
    occupied = args.output / 'occupied.csv'; occupied.write_bytes(b'preserve existing destination\n')
    invoke('overwrite', args.oracle / 'values.input.csv', occupied, guard=True)
    invoke('repeated_flag', args.oracle / 'values.input.csv', args.output / 'repeat.csv', guard=True,
           extra=['--input', str(args.oracle / 'values.input.csv')])
    empty = args.output / 'empty.input.csv'
    empty.write_text('index_kind,index_value,value_kind,value\n')
    empty_out = args.output / 'empty.actual.csv'
    invoke('empty', empty, empty_out)
    if empty_out.exists() and list(csv.reader(empty_out.open())) != [['index_kind', 'index_value', 'value_kind', 'value']]:
        mismatches.append(dict(name='empty_header'))
    labels = args.output / 'labels.input.csv'
    with labels.open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['index_kind', 'index_value', 'value_kind', 'value'])
        for kind, value in [('str', ''), ('str', 'ceb1'), ('null', ''), ('str', 'ceb1'), ('float', '8000000000000000')]:
            writer.writerow([kind, value, 'str', '2e35'])
    label_out = args.output / 'labels.actual.csv'
    invoke('opaque_duplicate_labels', labels, label_out)
    if label_out.exists():
        original = list(csv.reader(labels.open()))[1:]; current = list(csv.reader(label_out.open()))[1:]
        if [x[:2] for x in original] != [x[:2] for x in current]: mismatches.append(dict(name='opaque_labels'))
    report = dict(schema='priva.coerce-pext-native-comparison.v1', binary_sha256=args.binary_sha256,
        oracle_manifest_sha256=args.oracle_manifest_sha256, successful_scalar_rows=len(oracle['successful_cases']),
        expected_source_errors=len(oracle['source_failures']), unsupported_diagnostics_not_parity=len(oracle['unsupported_diagnostics']),
        native=dict(job_id=os.environ.get('SLURM_JOB_ID'), step_id=os.environ.get('SLURM_STEP_ID')),
        results=results, mismatches=mismatches, successful=not mismatches,
        comparison='exact typed rows/order/index tokens/f64 bits and supported source exception category/message; no tolerances',
        scientific_acceptance=False, production_backend_selected=False)
    (args.output / 'comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items() if key not in ('results', 'mismatches')}))
    if mismatches:
        print(json.dumps(mismatches))
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest='command', required=True).add_parser('compare')
    for name in ('binary', 'oracle', 'output'): command.add_argument('--' + name, type=Path, required=True)
    for name in ('binary-sha256', 'oracle-manifest-sha256'): command.add_argument('--' + name, required=True)
    compare(parser.parse_args())


if __name__ == '__main__':
    main()
