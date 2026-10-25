#!/usr/bin/env python3
"""Reuse a completed *original* AlphaMissense domain-statistics producer.

2026-10-24 CT: an original-producer output must not impersonate a chromosome
builder generation. This helper only copies/publishes existing bytes; it does
not load pickle, rerun science, or qualify a replacement calculation.

Commands (all options named)::

  python alphamissense_original_cache.py snapshot --source-script SCRIPT --input-vcf VCF
  python alphamissense_original_cache.py import --source-script SCRIPT --input-vcf VCF \
      --original-pickle SAVED --execution-evidence RUN.json --output-pickle DEST --cache-root ROOT
  python alphamissense_original_cache.py check --source-script SCRIPT --input-vcf VCF \
      --output-pickle DEST --cache-root ROOT

RUN.json is an explicit, externally reviewed attestation, not proof manufactured
by this helper: schema RUN_SCHEMA, complete true, recipe from snapshot, output
{path,size,sha256}, native_execution {state:"COMPLETED",exit_code:"0:0",records:
[{path,sha256},...]}, and source_binding/input_binding {basis,description}.
Each basis is "recorded-at-run" or "retrospective-attestation". Retrospective
attestations must include a nonempty "limitation" describing missing historical
proof. Snapshot observes NOW and must never be described as a pre-run capture.

Import streams a private copy, fsyncs it, publishes without replacing existing
files, then closes its separate receipt. Interrupted copies remain retained;
retrying the same request recovers owned publications. Warm checks compare
unchanged input/output metadata and producer CONTENT, avoiding a repeated large
payload hash. Metadata excludes node-local st_dev. Source paths can differ
between matching PriVA/PriVA-nextflow checkouts. Bundled data is not inspected.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import stat
import sys
import uuid

PRODUCER_SHA256 = 'cc7b2d60665e8052e6362107e359cc133e31c398f243589ee63ad7e8f90341c4'
SCHEMA = 'priva-alphamissense-original-publication-v1'
RUN_SCHEMA = 'priva-alphamissense-original-execution-v1'
READY_SUFFIX = '.original.complete.json'
PENDING_SUFFIX = '.original.pending.json'


def require(value, message):
    if not value:
        raise ValueError(message)


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n').encode()


def identity(path):
    info = Path(path).stat()
    require(stat.S_ISREG(info.st_mode) and info.st_size > 0,
            'Expected a regular nonempty file: ' + str(path))
    return {name: str(getattr(info, name)) for name in
            ('st_ino', 'st_mode', 'st_size', 'st_mtime_ns', 'st_ctime_ns', 'st_nlink')}


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def entry(path, *, hashed=False):
    path = Path(path).absolute()
    before = identity(path)
    result = dict(path=str(path), resolved=str(path.resolve(strict=True)), identity=before)
    if hashed:
        result['sha256'] = digest(path)
    require(identity(path) == before, 'File changed while reading: ' + str(path))
    return result


def snapshot(*, source_script, input_vcf):
    source = entry(source_script, hashed=True)
    require(source['sha256'] == PRODUCER_SHA256, 'Scientific producer differs from pinned original content')
    vcf = Path(input_vcf).absolute()
    members = [vcf, *(Path(str(vcf) + suffix) for suffix in ('.tbi', '.csi')
                      if os.path.lexists(str(vcf) + suffix))]
    require(len(members) > 1, 'An indexed input VCF is required')
    return dict(producer_sha256=source['sha256'], inputs=[entry(path) for path in members])


def read_json(path):
    require(not Path(path).is_symlink(), 'Preserve linked metadata: ' + str(path))
    before = identity(path)
    with Path(path).open() as stream:
        result = json.load(stream)
    require(identity(path) == before and isinstance(result, dict), 'Metadata changed or is not an object')
    return result


def managed(path, cache_root):
    path, root = Path(path).absolute(), Path(cache_root).absolute()
    require('..' not in path.parts and '..' not in root.parts and path != root and path.is_relative_to(root),
            'Destination leaves configured cache_root: ' + str(path))
    require(path.resolve().is_relative_to(root.resolve()), 'Destination aliases outside cache_root: ' + str(path))
    for item in (path, *path.parents):
        require(not item.is_symlink(), 'Preserve linked destination: ' + str(item))
        if item == root:
            break
    return path


def no_alias(path, protected):
    path = Path(path)
    for name in protected:
        other = Path(name)
        require(path.resolve() != other.resolve() and not
                (path.exists() and other.exists() and os.path.samefile(path, other)),
                'Publication aliases selected input/source/evidence: ' + str(path))


def sync_directory(folder):
    fd = os.open(folder, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def write_json(path, value, *, replace=False):
    """Exclusive final publication; replacement is only for our owned journal."""
    path = Path(path)
    temporary = path.with_name(path.name + '.partial-' + uuid.uuid4().hex)
    with temporary.open('xb') as stream:
        payload = encoded(value)
        require(stream.write(payload) == len(payload), 'Short publication metadata write')
        stream.flush(); os.fsync(stream.fileno())
    if replace:
        os.replace(temporary, path)
    else:
        os.link(temporary, path)  # FileExistsError preserves nonowned metadata.
        temporary.unlink()
    sync_directory(path.parent)


@contextmanager
def locked(path):
    fd = os.open(path, os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
    try:
        require(stat.S_ISREG(os.fstat(fd).st_mode), 'Publication lock is not regular')
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield
    finally:
        os.close(fd)


def validate_attestation(value, recipe):
    require(isinstance(value, dict), 'Original execution attestation must be an object')
    require(value.get('schema') == RUN_SCHEMA and value.get('complete') is True,
            'Missing original-producer execution completion attestation')
    require(value.get('recipe') == recipe, 'Execution source/input/index attestation differs from selection')
    for field in ('source_binding', 'input_binding'):
        binding = value.get(field, {})
        require(isinstance(binding, dict) and binding.get('basis') in ('recorded-at-run', 'retrospective-attestation') and
                isinstance(binding.get('description'), str) and binding['description'].strip(),
                'Missing explicit historical binding: ' + field)
        if binding['basis'] == 'retrospective-attestation':
            require(isinstance(binding.get('limitation'), str) and binding['limitation'].strip(),
                    'Retrospective binding must disclose missing historical proof: ' + field)
    native = value.get('native_execution', {})
    require(isinstance(native, dict) and native.get('state') == 'COMPLETED' and native.get('exit_code') == '0:0',
            'Native producer execution did not complete successfully')
    records = native.get('records')
    require(isinstance(records, list) and records, 'Native execution evidence records are required')
    return records


def evidence(*, execution_evidence, original_pickle, recipe):
    evidence_file = entry(execution_evidence, hashed=True)
    value = read_json(execution_evidence)
    records = validate_attestation(value, recipe)
    for record in records:
        current = entry(record['path'], hashed=True)
        require(current['sha256'] == record.get('sha256'), 'Native evidence checksum differs')
    original = entry(original_pickle)
    output = value.get('output', {})
    require(output.get('path') == original['path'] and output.get('size') == original['identity']['st_size'] and
            isinstance(output.get('sha256'), str) and len(output['sha256']) == 64 and
            all(char in '0123456789abcdef' for char in output['sha256']), 'Original output attestation differs')
    return dict(file=evidence_file, attestation=value, original=original)


def publication_paths(*, output_pickle, cache_root, protected=()):
    output = managed(output_pickle, cache_root)
    paths = [output, *(managed(str(output) + suffix, cache_root) for suffix in
                       (READY_SUFFIX, PENDING_SUFFIX, '.original.lock'))]
    for path in paths:
        no_alias(path, protected)
    return paths


def validate_receipt(value, *, output, recipe):
    require(value.get('schema') == SCHEMA and value.get('complete') is True and
            value.get('output_path') == str(output) and value.get('recipe') == recipe,
            'Original-producer receipt differs from selected source/input/index')
    require(value.get('output', {}) == entry(output), 'Published original output changed')
    execution = value.get('execution', {})
    attestation = execution.get('attestation', {})
    validate_attestation(attestation, recipe)
    require(value.get('sha256') == attestation.get('output', {}).get('sha256') and
            value['output']['identity']['st_size'] == attestation.get('output', {}).get('size'),
            'Original output digest differs from execution attestation')
    require(execution.get('file') == entry(execution['file']['path'], hashed=True),
            'Original execution evidence changed')
    require(read_json(execution['file']['path']) == attestation,
            'Embedded execution attestation differs from its retained evidence')
    for record in attestation['native_execution']['records']:
        require(entry(record['path'], hashed=True)['sha256'] == record.get('sha256'),
                'Retained native execution evidence changed')


def check(*, source_script, input_vcf, output_pickle, cache_root):
    ready, pending = (Path(str(output_pickle) + suffix) for suffix in (READY_SUFFIX, PENDING_SUFFIX))
    if not os.path.lexists(ready) and not os.path.lexists(pending):
        return None  # CLI exit 3 means this producer does not own the selection.
    recipe = snapshot(source_script=source_script, input_vcf=input_vcf)
    protected = [source_script, *(item['path'] for item in recipe['inputs'])]
    output, ready, pending, _ = publication_paths(output_pickle=output_pickle, cache_root=cache_root, protected=protected)
    require(not os.path.lexists(pending), 'Original publication is pending; recover with the same import command')
    before = entry(ready, hashed=True)
    value = read_json(ready)
    validate_receipt(value, output=output, recipe=recipe)
    require(snapshot(source_script=source_script, input_vcf=input_vcf) == recipe and entry(ready, hashed=True) == before,
            'Source/input/index/receipt changed during original-cache check')
    return dict(output=str(output), receipt=str(ready), reused=True)


def copy_original(original, stage, expected):
    before = entry(original)
    result = hashlib.sha256()
    with Path(original).open('rb') as source, stage.open('xb') as target:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            require(target.write(block) == len(block), 'Short original-cache publication write')
            result.update(block)
        target.flush(); os.fsync(target.fileno())
    require(entry(original) == before, 'Original producer output changed during copy')
    require(identity(stage)['st_size'] == expected['size'] and result.hexdigest() == expected['sha256'],
            'Incomplete original output or publication checksum mismatch')
    return entry(stage)


def import_original(*, source_script, input_vcf, original_pickle, execution_evidence, output_pickle, cache_root):
    recipe = snapshot(source_script=source_script, input_vcf=input_vcf)
    execution = evidence(execution_evidence=execution_evidence, original_pickle=original_pickle, recipe=recipe)
    protected = [source_script, original_pickle, execution_evidence,
                 *(item['path'] for item in recipe['inputs']),
                 *(item['path'] for item in execution['attestation']['native_execution']['records'])]
    output, ready, pending, lock_path = publication_paths(output_pickle=output_pickle, cache_root=cache_root, protected=protected)
    request = dict(recipe=recipe, output_path=str(output), execution=execution)
    output.parent.mkdir(parents=True, exist_ok=True)

    def stable():
        require(snapshot(source_script=source_script, input_vcf=input_vcf) == recipe and
                evidence(execution_evidence=execution_evidence, original_pickle=original_pickle, recipe=recipe) == execution,
                'Source/input/index/original/execution evidence changed during publication')

    with locked(lock_path):
        for suffix in ('.complete.json', '.publish.pending.json'):
            require(not os.path.lexists(str(output) + suffix),
                    'Preserve chromosome-builder publication; original import must use a separate destination')
        if os.path.lexists(pending):
            journal = read_json(pending)
            require(journal.get('schema') == SCHEMA and journal.get('request') == request,
                    'Preserve nonowned or changed original publication journal')
        else:
            if os.path.lexists(ready):
                value = read_json(ready)
                validate_receipt(value, output=output, recipe=recipe)
                require(value.get('execution') == execution, 'Preserve another original execution receipt')
                stable()
                return dict(output=str(output), receipt=str(ready), reused=True)
            require(not os.path.lexists(output), 'Preserve existing output without original-producer ownership receipt')
            journal = dict(schema=SCHEMA, request=request,
                           stage=str(output.with_name(output.name + '.original-copy-' + uuid.uuid4().hex)), retained=[])
            write_json(pending, journal)
        stage = managed(journal['stage'], cache_root)
        require(stage.parent == output.parent and stage.name.startswith(output.name + '.original-copy-'), 'Invalid journal stage')
        no_alias(stage, protected)
        if os.path.lexists(ready):
            value = read_json(ready)
            validate_receipt(value, output=output, recipe=recipe)
            require(value.get('execution') == execution, 'Preserve another original execution receipt')
        else:
            if os.path.lexists(output):
                completed = journal.get('stage_ready', {})
                observed = entry(output, hashed=True)
                require(observed['identity']['st_ino'] == completed.get('identity', {}).get('st_ino') and
                        observed['identity']['st_size'] == execution['attestation']['output']['size'] and
                        observed['sha256'] == execution['attestation']['output']['sha256'],
                        'Preserve output not owned by this interrupted publication')
                if os.path.lexists(stage):
                    require(os.path.samefile(stage, output), 'Preserve changed recovery stage')
                    stage.unlink()
            else:
                if 'stage_ready' in journal:
                    require(entry(stage) == journal['stage_ready'] and digest(stage) == execution['attestation']['output']['sha256'],
                            'Preserve corrupted completed staging output')
                else:
                    if os.path.lexists(stage):
                        journal['retained'].append(str(stage))
                        stage = output.with_name(output.name + '.original-copy-' + uuid.uuid4().hex)
                        journal['stage'] = str(stage)
                        write_json(pending, journal, replace=True)
                    journal['stage_ready'] = copy_original(original_pickle, stage, execution['attestation']['output'])
                    stable()
                    write_json(pending, journal, replace=True)
                stable()
                os.link(stage, output)  # Atomic no-clobber publication; never os.replace(output).
                stage.unlink()
                sync_directory(output.parent)
            stable()
            published = entry(output)
            require(all(published['identity'][field] == journal['stage_ready']['identity'][field]
                        for field in ('st_ino', 'st_mode', 'st_size', 'st_mtime_ns')),
                    'Published payload changed before its receipt closed')
            value = dict(schema=SCHEMA, complete=True, **request, output=published,
                         sha256=execution['attestation']['output']['sha256'],
                         finished_utc=datetime.now(timezone.utc).isoformat(), installation_accepted=False)
            write_json(ready, value)
        stable()
        validate_receipt(value, output=output, recipe=recipe)
        pending.unlink(); sync_directory(output.parent)
        return dict(output=str(output), receipt=str(ready), reused=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('snapshot', 'check', 'import'):
        command = commands.add_parser(name)
        command.add_argument('--source-script', required=True)
        command.add_argument('--input-vcf', required=True)
        if name != 'snapshot':
            command.add_argument('--output-pickle', required=True)
            command.add_argument('--cache-root', required=True)
        if name == 'import':
            command.add_argument('--original-pickle', required=True)
            command.add_argument('--execution-evidence', required=True)
    args = vars(parser.parse_args())
    name = args.pop('command')
    try:
        result = (import_original if name == 'import' else globals()[name])(**args)
        if result is None:
            return 3
        print(json.dumps(result, sort_keys=True))
        return 0
    except (OSError, ValueError, KeyError, TypeError) as error:
        print('Original AlphaMissense cache: ' + str(error), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
