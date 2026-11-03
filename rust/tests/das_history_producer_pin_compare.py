#!/usr/bin/env python3
"""Check explicit history-producer pins and unchanged default archive bytes.

check --binary FILE --reference-binary FILE --source-dir DIR --work-root DIR
Uses tiny generated inputs and native builders only. The artificial producer
hash tests admission mechanics; real historical-profile acceptance is separate.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check(args):
    import pyarrow.parquet as pq
    root = args.work_root
    root.mkdir(parents=True, exist_ok=False)
    source = args.source_dir / 'rust/tests/alpha_domain_compare.py'
    spec = importlib.util.spec_from_file_location('history_pin_cases', source)
    factory = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(factory)
    vcf = factory.make_case(root, factory.cases()['normal'])
    rows = []

    def run(label, command):
        result = subprocess.run(list(map(str, command)), capture_output=True)
        (root / (label + '.stdout.log')).write_bytes(result.stdout)
        (root / (label + '.stderr.log')).write_bytes(result.stderr)
        return result

    statistics, history, exported = [root / name for name in ['statistics.parquet', 'history.parquet', 'history.native']]
    commands = [
        ['alphamissense-domain', 'build', '--vcf', vcf, '--output', statistics],
        ['alpha-exon-history', 'build', '--vcf', vcf, '--hash-seed', '0', '--output', history],
        ['alpha-exon-history', 'export', '--input', history, '--output', exported],
    ]
    for index, command in enumerate(commands):
        assert run('prepare_' + str(index), [args.binary, *command]).returncode == 0
    source_builder = pq.ParquetFile(statistics).metadata.metadata[b'priva.builder_source_sha256'].decode()
    header, body = exported.read_text().split('\n', 1)
    original_header = json.loads(header.removeprefix('# '))
    producer = original_header['priva.builder_source_sha256']
    assert len(producer) == 64
    inputs = {str(p): (p.stat().st_size, p.stat().st_mtime_ns, digest(p)) for p in [vcf, statistics, history, exported]}

    def project(binary, label, selected=exported, producer_pin=None, succeeds=True,
                failure_exit=1, failure_message=b'ValueError'):
        output = root / (label + '.rkyv')
        sentinel = b'previous complete publication\n'
        if not succeeds:
            output.write_bytes(sentinel)
        command = [binary, 'das-domains', 'transcripts', '--statistics', statistics,
                   '--exon-history', selected, '--source-builder-sha256', source_builder,
                   '--python-hash-seed', '0', '--output', output]
        if producer_pin is not None:
            command += ['--history-builder-sha256', producer_pin]
        result = run(label, command)
        passed = result.returncode == (0 if succeeds else failure_exit)
        if not succeeds:
            passed &= output.read_bytes() == sentinel
            passed &= failure_message in result.stderr
        rows.append(dict(name=label, passed=passed, returncode=result.returncode,
                         output_preserved=None if succeeds else output.read_bytes() == sentinel))
        assert passed, label
        return output

    baseline = project(args.reference_binary, 'reference_default')
    current = project(args.binary, 'current_default')
    assert baseline.read_bytes() == current.read_bytes()
    rows.append(dict(name='default_archive_bytes', passed=True))

    def provenance(label, output):
        value = run(label, [args.binary, 'das-domains', 'transcript-provenance', '--input', output])
        assert value.returncode == 0
        return json.loads(value.stdout)

    default_provenance = provenance('default_provenance', current)
    consumer_source = digest(args.source_dir / 'rust/src/alpha_domain/das_transcripts.rs')
    assert dict(default_provenance)['history_builder_sha256'] == producer
    assert 'transcript_consumer_source_sha256' not in dict(default_provenance)
    pinned = project(args.binary, 'explicit_current', producer_pin=producer)
    pinned_provenance = provenance('pinned_provenance', pinned)
    assert pinned_provenance == default_provenance[:-2] + [
        ['transcript_consumer_source_sha256', consumer_source], *default_provenance[-2:]]

    # This is deliberately synthetic, not a relabeled scientific artifact.
    artificial = 'ab' * 32
    altered = dict(original_header, **{'priva.builder_source_sha256': artificial})
    historical = root / 'synthetic_older_header.native'
    historical.write_text('# ' + json.dumps(altered) + '\n' + body)
    project(args.binary, 'default_refuses_older', selected=historical, succeeds=False)
    project(args.binary, 'wrong_producer', selected=historical, producer_pin='cd' * 32, succeeds=False)
    accepted = project(args.binary, 'explicit_older', selected=historical, producer_pin=artificial)
    declared = dict(provenance('older_provenance', accepted))
    assert declared['history_builder_sha256'] == artificial
    assert declared['transcript_consumer_source_sha256'] == consumer_source
    rows.append(dict(name='producer_and_consumer_remain_distinct', passed=True))
    exports = []
    for label, output in [('default_values', current), ('pinned_values', pinned), ('older_values', accepted)]:
        result = run(label, [args.binary, 'cache-export', 'alphamissense_tranx_domain_map', '--input', output, '--output', '-'])
        assert result.returncode == 0
        exports.append(result.stdout)
    assert exports[0] == exports[1] == exports[2]
    rows.append(dict(name='complete_ordered_values', passed=True))

    for field, value in [('priva.history.schema', 'wrong'), ('priva.stage', 'raw'),
                         ('priva.history.profile', 'unknown'), ('priva.history.seed', '4294967296'),
                         ('priva.input', 'wrong input')]:
        path = root / (field + '.native')
        path.write_text('# ' + json.dumps(dict(altered, **{field: value})) + '\n' + body)
        project(args.binary, 'reject_' + field, selected=path, producer_pin=artificial, succeeds=False)
    for name, text in [('legacy', body), ('truncated', '# {\n'),
                       ('missing_rows', '# ' + json.dumps(altered) + '\n'),
                       ('extra_rows', '# ' + json.dumps(altered) + '\n' + body + body)]:
        path = root / (name + '.native')
        path.write_text(text)
        project(args.binary, 'reject_' + name, selected=path, producer_pin=artificial, succeeds=False)
    for name, value in [('short_pin', 'ab'), ('nonhex_pin', 'z' * 64)]:
        project(args.binary, name, selected=historical, producer_pin=value, succeeds=False,
                failure_exit=2, failure_message=b'--history-builder-sha256 requires SHA256')
    assert all((Path(p).stat().st_size, Path(p).stat().st_mtime_ns, digest(p)) == expected
               for p, expected in inputs.items())
    result = dict(passed=all(row['passed'] for row in rows), cases=len(rows), rows=rows,
                  default_archive_bytes_exact=True, complete_ordered_values_exact=True,
                  original_inputs_unchanged=True, artificial_producer_fixture_only=True,
                  full_real_historical_profile_qualified=False)
    (root / 'complete.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(required=True)
    command = sub.add_parser('check')
    for name in ('binary', 'reference-binary', 'source-dir', 'work-root'):
        command.add_argument('--' + name, type=Path, required=True)
    command.set_defaults(function=check)
    args = parser.parse_args()
    args.function(args)


if __name__ == '__main__':
    main()
