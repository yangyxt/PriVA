#!/usr/bin/env python3
"""Install selected-release CADD transcript candidates; never activate scoring.

python cadd_transcript_resources.py plan --config CONFIG
python cadd_transcript_resources.py install --config CONFIG [--verify-only]
python cadd_transcript_resources.py adopt --config CONFIG --component peptide \
    --raw-receipt RAW.json --raw-receipt-sha256 SHA \
    --view-receipt VIEW.json --view-receipt-sha256 SHA

Raw official Ensembl archives and unchanged plain views live under cache_root.
The selected assembly/release determines URLs; GRCh37 retains annotation-base75
lineage. Acquisition/explicit adoption records provenance, not feature coverage,
transcript compatibility, model inference or scientific/public qualification.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location('cadd_candidate_references', HERE / 'reference_resources.py')
R = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(R)
ENSEMBL = 'https://ftp.ensembl.org/pub'
FLAGS = {'scientific_qualified': False, 'feature_semantics_accepted': False,
         'cadd_enabled': False, 'model_inference_performed': False}


def read_json(path):
    path = Path(path)
    R.require(path.stat().st_size <= 1024 * 1024, 'Oversized candidate metadata')
    before = R.identity(path)
    value = json.loads(path.read_text())
    R.require(R.identity(path) == before, 'Candidate metadata changed during reading')
    return value


def publish_json(path, record):
    if os.path.lexists(path):
        R.require(read_json(path) == record, f'Existing candidate receipt differs; retained: {path}')
        return
    with tempfile.NamedTemporaryFile(mode='w', prefix='.candidate-', dir=path.parent, delete=False) as stream:
        json.dump(record, stream, indent=2, sort_keys=True); stream.write('\n')
        temporary = Path(stream.name)
    try: R.publish_new(temporary, path)
    finally: temporary.unlink(missing_ok=True)


def plan(config):
    assembly = {'hg19': 'GRCh37', 'GRCh37': 'GRCh37', 'hg38': 'GRCh38', 'GRCh38': 'GRCh38'}.get(config.get('assembly'))
    release = str(config.get('vep_installed_version', ''))
    R.require(assembly and re.fullmatch(r'[1-9][0-9]*', release), 'Select assembly and a numeric VEP release')
    root = R.checked_path(config.get('cache_root', ''))
    R.require(root.is_absolute() and root != Path(root.anchor), 'Set a dedicated absolute cache_root')
    directory = root / 'references' / ('Ensembl' + release) / assembly
    origin = ENSEMBL + ('/grch37' if assembly == 'GRCh37' else '') + '/release-' + release
    sources = {'peptide': ('fasta/homo_sapiens/pep', f'Homo_sapiens.{assembly}.pep.all.fa.gz')}
    if assembly == 'GRCh38': sources['gtf'] = ('gtf/homo_sapiens', f'Homo_sapiens.{assembly}.{release}.gtf.gz')
    return {'schema_version': 1, 'assembly': assembly, 'vep_release': release,
            'annotation_base_release': '75' if assembly == 'GRCh37' else release,
            'cache_root': str(root), 'directory': str(directory), 'resources': {
                name: {'url': origin + '/' + folder + '/' + filename,
                       'checksums_url': origin + '/' + folder + '/CHECKSUMS',
                       'path': str(directory / filename), 'view': str(directory / 'views' / filename[:-3])}
                for name, (folder, filename) in sources.items()}, **FLAGS}


@contextmanager
def locked(selected):
    root, directory = Path(selected['cache_root']), Path(selected['directory'])
    R.require(directory.resolve().is_relative_to(root.resolve()), 'Candidate directory escapes cache_root')
    for part in directory.relative_to(root).parts:
        root /= part
        R.require(not root.is_symlink(), f'Refuse a linked candidate directory: {root}')
    directory.mkdir(parents=True, exist_ok=True)
    with R.staging(directory / '.install.lock').open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield directory


def source_pin(resource):
    with urlopen(resource['checksums_url'], timeout=60) as response:
        text = response.read(1024 * 1024 + 1)
    R.require(len(text) <= 1024 * 1024, 'Oversized publisher CHECKSUMS')
    rows = [line.split() for line in text.decode('ascii').splitlines()]
    rows = [row for row in rows if len(row) == 3 and row[2] == Path(resource['path']).name]
    R.require(len(rows) == 1 and all(value.isdigit() for value in rows[0][:2]), 'Missing/ambiguous publisher checksum')
    with urlopen(Request(resource['url'], method='HEAD'), timeout=60) as response:
        length = response.headers.get('Content-Length', '')
        R.require(length.isdigit() and int(length) > 0, 'Publisher archive length is unavailable')
        return {'url': resource['url'], 'checksums_url': resource['checksums_url'],
                'publisher_checksums_sha256': hashlib.sha256(text).hexdigest(),
                'publisher_bsd_sum': list(map(int, rows[0][:2])), 'bytes': int(length),
                'etag': response.headers.get('ETag'), 'last_modified': response.headers.get('Last-Modified')}


def verify_archive(path, source):
    before = R.identity(path)
    R.require(path.stat().st_size == source['bytes'], 'Candidate archive length differs; retained')
    measured = list(map(int, subprocess.check_output(['sum', '-r', str(path)], text=True).split()[:2]))
    R.require(measured == source['publisher_bsd_sum'], 'Candidate archive publisher checksum differs; retained')
    checksum = R.digest(path)
    R.require(R.identity(path) == before, 'Candidate archive changed during verification')
    return checksum


def acquire(source, destination):
    """Exact resumable range into owned staging; publish only after checksum."""
    partial = R.staging(str(destination) + '.partial')
    offset = partial.stat().st_size if partial.exists() else 0
    R.require(offset <= source['bytes'], 'Candidate partial exceeds source length; retained')
    if offset < source['bytes']:
        headers = {'Range': f'bytes={offset}-'} if offset else {}
        if source.get('etag'): headers['If-Match'] = source['etag']
        elif source.get('last_modified'): headers['If-Unmodified-Since'] = source['last_modified']
        with urlopen(Request(source['url'], headers=headers), timeout=60) as response:
            R.require(response.status in (200, 206), 'Unexpected candidate HTTP status')
            if offset:
                R.require(response.status == 206 and response.headers.get('Content-Range') ==
                          f"bytes {offset}-{source['bytes'] - 1}/{source['bytes']}", 'Publisher refused exact resume; partial retained')
            if source.get('etag') and response.headers.get('ETag'):
                R.require(response.headers['ETag'] == source['etag'], 'Candidate ETag changed; partial retained')
            with partial.open('ab' if offset else 'wb') as stream:
                shutil.copyfileobj(response, stream, R.BLOCK)
    checksum = verify_archive(partial, source)
    R.publish_new(partial, destination)
    return checksum


def metadata(selected, component, stage, path, checksum, **extra):
    return {'schema_version': 1, 'stage': stage, 'component': component,
            **{key: selected[key] for key in ('assembly', 'vep_release', 'annotation_base_release')},
            'path': str(path), 'sha256': checksum, 'identity': R.identity(path), **extra, **FLAGS}


def verified(record, selected, component, stage, path):
    R.require(all(record.get(key) == value for key, value in {
        'schema_version': 1, 'stage': stage, 'component': component, 'path': str(path),
        **{key: selected[key] for key in ('assembly', 'vep_release', 'annotation_base_release')}, **FLAGS}.items()),
        'Candidate receipt scope differs; existing bytes retained')
    R.require(not path.is_symlink(), f'Refuse a linked managed candidate leaf: {path}')
    R.require(re.fullmatch('[0-9a-f]{64}', str(record.get('sha256', ''))) and path.is_file() and
              (R.identity(path) == record.get('identity') or R.digest(path) == record['sha256']),
              f'Candidate resource changed; retained: {path}')
    return record


def validate_view(path, component, assembly):
    """Bounded assembly/format check for downloaded external candidates only."""
    before = R.identity(path)
    with path.open('rb') as stream: prefix = stream.read(65536).decode('ascii')
    lines = prefix.splitlines()
    if component == 'peptide':
        R.require(lines and lines[0].startswith('>') and ':' + assembly + ':' in lines[0],
                  'Peptide FASTA header does not identify the selected assembly')
    else:
        rows = [line for line in lines if line and not line.startswith('#')]
        R.require(any(line.startswith('#!genome-build ' + assembly) or
                      line == '#!genome-version ' + assembly for line in lines) and
                  rows and len(rows[0].split('\t')) == 9, 'GTF header/format does not identify the selected assembly')
    R.require(R.identity(path) == before, 'Candidate view changed during bounded header inspection')


def resource_install(selected, name, resource, verify_only=False):
    directory = Path(selected['directory']); path = Path(resource['path']); view = Path(resource['view'])
    R.require(not view.parent.is_symlink(), 'Refuse a linked candidate view directory')
    R.require(not path.is_symlink() and not view.is_symlink(), 'Refuse a linked managed candidate leaf')
    raw_receipt, view_receipt = directory / (name + '.raw.ready.json'), directory / (name + '.view.ready.json')
    if raw_receipt.is_file():
        raw = verified(read_json(raw_receipt), selected, name, 'raw', path)
        R.require(raw['source']['url'] == resource['url'] and raw['source']['checksums_url'] == resource['checksums_url'],
                  'Candidate source/release changed; retained')
    else:
        R.require(not verify_only, 'Selected transcript archive is not installed')
        source_receipt = directory / (name + '.download.json')
        if source_receipt.is_file(): source = read_json(source_receipt)
        else:
            R.require(not os.path.lexists(path) and not os.path.lexists(str(path) + '.partial'),
                      'Unidentified candidate archive/partial retained; use explicit verified adoption')
            source = source_pin(resource); publish_json(source_receipt, source)
        R.require(source['url'] == resource['url'] and source['checksums_url'] == resource['checksums_url'],
                  'Staged candidate source differs; retained')
        checksum = verify_archive(path, source) if path.is_file() else acquire(source, path)
        raw = metadata(selected, name, 'raw', path, checksum, source=source)
        publish_json(raw_receipt, raw)
    if view_receipt.is_file():
        prepared = verified(read_json(view_receipt), selected, name, 'view', view)
        R.require(prepared.get('raw_sha256') == raw['sha256'], 'Candidate view/archive identity differs')
    else:
        R.require(not verify_only, 'Selected transcript plain view is not installed')
        view.parent.mkdir(exist_ok=True)
        R.require(not view.parent.is_symlink(), 'Refuse a linked candidate view directory')
        if os.path.lexists(view):
            checksum = R.source_digest(path)
            R.require(view.is_file() and R.digest(view) == checksum, 'Existing candidate view differs; retained')
        else:
            # A unique partial retains interrupted bytes; no existing file is truncated.
            try:
                with tempfile.NamedTemporaryFile(prefix=view.name + '.', suffix='.partial', dir=view.parent, delete=False) as stream:
                    partial = Path(stream.name); checksum = R.source_digest(path, stream)
            except (OSError, EOFError) as error:
                raise ValueError('Candidate gzip decompression failed; staged bytes retained') from error
            validate_view(partial, name, selected['assembly'])
            R.publish_new(partial, view)
        validate_view(view, name, selected['assembly'])
        prepared = metadata(selected, name, 'view', view, checksum, raw_sha256=raw['sha256'],
                            operation='gzip_decompression_preserving_all_original_bytes_and_record_order')
        publish_json(view_receipt, prepared)
    return {'raw': raw, 'view': prepared, 'raw_receipt': str(raw_receipt), 'view_receipt': str(view_receipt)}


def persistent_reference_identity(path):
    """Keep stable reference fields in receipts; st_dev is mount-local across nodes."""
    return dict(zip(('st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns'), R.identity(path)[1:]))


def install(config, verify_only=False):
    selected = plan(config)
    with locked(selected) as directory:
        resources = {name: resource_install(selected, name, resource, verify_only)
                     for name, resource in selected['resources'].items()}
        reference = None
        if config.get('ref_genome'):
            fasta = R.checked_path(config['ref_genome']); index = Path(str(fasta) + '.fai')
            R.require(fasta.is_file() and index.is_file(), 'Selected reference/index is not installed')
            reference = {'path': str(fasta), 'identity': persistent_reference_identity(fasta),
                         'index': str(index), 'index_identity': persistent_reference_identity(index),
                         'identity_policy': 'st_ino,st_size,st_mtime_ns,st_ctime_ns; node_local_st_dev_excluded'}
        record = {'schema_version': 1, 'stage': 'cadd-selected-transcript-candidate',
                  'plan': selected, 'resources': resources, 'selected_reference': reference,
                  'integrity': 'publisher_BSD_sum_and_length_plus_local_SHA256; publisher_cryptographic_digest_not_claimed',
                  'environment_candidates': {'CADD_ESM_TRANSCRIPTS': selected['resources']['peptide']['view'],
                      **({'CADD_MMSPLICE_GTF': selected['resources']['gtf']['view']} if 'gtf' in selected['resources'] else {}),
                      **({'CADD_REFERENCE_FASTA': reference['path']} if reference else {})}, **FLAGS}
        fingerprint = hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()
        receipt = directory / ('candidate.' + fingerprint[:16] + '.ready.json')
        if not verify_only: publish_json(receipt, record)
        else: R.require(read_json(receipt) == record, 'Selected transcript candidate manifest differs')
    return {'path': str(receipt), 'sha256': R.digest(receipt), 'assembly': selected['assembly'],
            'vep_release': selected['vep_release'], 'annotation_base_release': selected['annotation_base_release'],
            'environment_candidates': record['environment_candidates'], **FLAGS}


def pinned_json(path, checksum):
    R.require(R.digest(path) == checksum, 'Supplied receipt checksum differs')
    record = read_json(path)
    R.require(R.digest(path) == checksum, 'Supplied receipt changed while reading')
    return record


def adopt(config, component, raw_receipt, raw_sha256, view_receipt, view_sha256):
    """Copy proven complete bytes while retaining supplied original identities."""
    selected = plan(config)
    R.require(component in selected['resources'], 'Component does not belong to the selected assembly')
    resource = selected['resources'][component]
    raw, view = pinned_json(raw_receipt, raw_sha256), pinned_json(view_receipt, view_sha256)
    expected_component = 'ensembl' + selected['vep_release'] + '_' + component
    for record in (raw, view):
        R.require(record.get('assembly') == selected['assembly'] and record.get('component') == expected_component and
                  record.get('scientific_qualified') is False and record.get('cadd_enabled') is False,
                  'Supplied candidate assembly/release/scope differs')
    R.require(raw.get('stage') == 'raw_acquisition' and raw.get('status') == 'source_acquired_not_scientifically_qualified' and
              len(raw.get('files', [])) == 1 and raw['source']['url'] == resource['url'] and
              raw['source']['publisher_checksums'] == resource['checksums_url'] and
              view.get('stage') == 'uncompressed_view' and view.get('status') == 'view_prepared_not_scientifically_qualified' and
              view.get('raw_receipt') == str(raw_receipt) and view.get('raw_receipt_sha256') == raw_sha256 and
              view.get('operation') == 'gzip_decompression_preserving_all_original_bytes_and_record_order',
              'Supplied candidate acquisition/view provenance differs')
    origin = {'raw_receipt': str(raw_receipt), 'raw_receipt_sha256': raw_sha256,
              'view_receipt': str(view_receipt), 'view_receipt_sha256': view_sha256,
              'originals_retained': True}
    members = [raw['files'][0], view]
    validated = {}
    for member in members:
        path = Path(member['path'])
        R.require(path.is_file() and re.fullmatch('[0-9a-f]{64}', member.get('sha256', '')) and
                  (R.identity(path) == member.get('identity') or R.digest(path) == member['sha256']),
                  'Supplied candidate original changed; retained')
        validated[str(path)] = R.identity(path)
    validate_view(Path(view['path']), component, selected['assembly'])
    with locked(selected) as directory:
        for member, destination, stage in ((members[0], Path(resource['path']), 'raw'),
                                            (members[1], Path(resource['view']), 'view')):
            source = Path(member['path']); before = R.identity(source)
            R.require(before == validated[str(source)], 'Supplied original changed before adoption')
            destination.parent.mkdir(exist_ok=True)
            R.require(not destination.parent.is_symlink(), 'Refuse a linked candidate parent')
            R.require(not destination.is_symlink(), 'Refuse a linked managed candidate leaf')
            ready = directory / (component + '.' + stage + '.ready.json')
            if ready.is_file():
                previous = verified(read_json(ready), selected, component, stage, destination)
                R.require(previous['sha256'] == member['sha256'], 'Managed candidate differs from supplied original')
                continue
            if os.path.lexists(destination):
                R.require(destination.is_file() and not destination.samefile(source) and R.digest(destination) == member['sha256'],
                          'Existing managed candidate differs; retained')
            else:
                with tempfile.NamedTemporaryFile(prefix='.adopt-', suffix='.partial', dir=destination.parent, delete=False) as stream:
                    partial = Path(stream.name)
                    with source.open('rb') as input_stream: shutil.copyfileobj(input_stream, stream, R.BLOCK)
                R.require(R.identity(source) == before and R.digest(partial) == member['sha256'], 'Original changed during copy')
                R.publish_new(partial, destination)
            R.require(R.identity(source) == before, 'Original changed during adoption')
            extra = {'origin': origin}
            if stage == 'raw':
                extra['source'] = {'url': resource['url'], 'checksums_url': resource['checksums_url'],
                    'publisher_checksums_sha256': raw['source']['publisher_checksum_manifest_sha256'],
                    'publisher_bsd_sum': raw['source']['publisher_bsd_sum'], 'bytes': raw['source']['bytes']}
            else: extra.update(raw_sha256=members[0]['sha256'], operation=view['operation'])
            publish_json(ready, metadata(selected, component, stage, destination, member['sha256'], **extra))
    return {'component': component, 'managed_directory': selected['directory'], 'originals_retained': True, **FLAGS}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='procedure', required=True)
    for name in ('plan', 'install', 'adopt'):
        child = commands.add_parser(name)
        inputs = child.add_mutually_exclusive_group(required=True)
        inputs.add_argument('--config', type=Path); inputs.add_argument('--config-json')
        if name == 'install': child.add_argument('--verify-only', action='store_true')
        if name == 'adopt':
            child.add_argument('--component', choices=('peptide', 'gtf'), required=True)
            for flag in ('raw-receipt', 'view-receipt'): child.add_argument('--' + flag, type=Path, required=True)
            for flag in ('raw-receipt-sha256', 'view-receipt-sha256'): child.add_argument('--' + flag, required=True)
    args = parser.parse_args(argv)
    try:
        if args.config_json: config = json.loads(args.config_json)
        else:
            import yaml
            config = yaml.safe_load(args.config.read_text())
        if args.procedure == 'plan': result = plan(config)
        elif args.procedure == 'install': result = install(config, args.verify_only)
        else: result = adopt(config, args.component, args.raw_receipt, args.raw_receipt_sha256,
                             args.view_receipt, args.view_receipt_sha256)
        print(json.dumps(result, sort_keys=True)); return 0
    except (ValueError, OSError, KeyError, TypeError) as error:
        print('CADD transcript candidate: ' + str(error), file=sys.stderr); return 1


if __name__ == '__main__': raise SystemExit(main())
