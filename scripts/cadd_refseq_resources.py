#!/usr/bin/env python3
"""Acquire pinned RefSeq candidates separately from CADD scientific acceptance.

python cadd_refseq_resources.py plan --config CONFIG
python cadd_refseq_resources.py acquire --config CONFIG
python cadd_refseq_resources.py check --config CONFIG

Only the reviewed GRCh38/merged-VEP116 source release is supported. ``plan``
reads local cache metadata; ``check`` is local and read-only. ``acquire`` saves
four original NCBI files beneath config.cache_root, with owned resumable staging
and no replacement of existing files. Downloads do not validate transcript,
sequence, coordinate or feature compatibility and do not activate a profile.
No Ensembl resource, upstream scorer, installer or existing result is changed.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import sys
import tempfile
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location('refseq_candidate_common', HERE / 'cadd_transcript_resources.py')
T = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(T)
R = T.R

RELEASE = 'GCF_000001405.40-RS_2025_08'
PREFIX = 'GCF_000001405.40_GRCh38.p14'
BASE = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/' + RELEASE + '/'
# Reviewed publisher MD5s and exact HTTP entity lengths, not a mutable "latest".
FILES = {
    'gff': (PREFIX + '_genomic.gff.gz', 78190483, '24b731562b9d4cae9e37b23404b5be16'),
    'protein': (PREFIX + '_protein.faa.gz', 28489875, '25a880236608ad046d7f9b65e895eb4b'),
    'rna': (PREFIX + '_rna.fna.gz', 135516534, 'b4a2ce202c90c0f24f22850c6bc7d774'),
    'assembly_report': (PREFIX + '_assembly_report.txt', 80454, '21f3ac4aa8245a99eb874082051b9dde'),
}
FLAGS = {**T.FLAGS, 'sequence_coordinate_validation_performed': False,
         'transcript_coverage_validated': False, 'raw_scoring': False}
MAX_METADATA = 1024 * 1024


def identity(path):
    """Persistent fields only: st_dev is node-local on this cluster."""
    path = R.checked_path(str(path)); stat = path.stat()
    return dict(path=str(path), target=str(path.resolve()), inode=stat.st_ino,
                size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                ctime_ns=stat.st_ctime_ns, mode=stat.st_mode)


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def read_json(path):
    return T.read_json(R.staging(str(path)))


def publish_json(path, value):
    T.publish_json(R.staging(str(path)), value)


def cache_source(config):
    """Pin the exact configured merged-cache info; never infer a new release."""
    R.require(isinstance(config, dict), 'Configuration must be a mapping')
    R.require(config.get('assembly') in ('hg38', 'GRCh38') and
              str(config.get('vep_installed_version')) == '116',
              'RefSeq candidate requires configured GRCh38 and VEP116')
    root_value = config.get('cache_root', '')
    R.require(isinstance(root_value, str) and Path(root_value).is_absolute(), 'Set an absolute cache_root')
    root = R.checked_path(root_value)
    R.require(root != Path(root.anchor), 'cache_root must be a dedicated directory')
    folders = config.get('cache_folders') or {}
    R.require(isinstance(folders, dict), 'cache_folders must be a mapping')
    folder = folders.get('vep_cache', 'VEP_cache')
    R.require(isinstance(folder, str) and folder not in ('', '.', '..') and Path(folder).name == folder,
              'vep_cache folder must be one directory name')
    cache = R.checked_path(config.get('vep_cache_dir') or str(root / folder))
    path = R.checked_path(str(cache / 'homo_sapiens_merged/116_GRCh38/info.txt'))
    before = identity(path)
    R.require(0 < before['size'] <= 65536, 'Missing/oversized merged-cache metadata')
    data = path.read_bytes(); values = {}
    for line in data.decode('utf-8').splitlines():
        if not line or line.startswith('#'): continue
        fields = line.split('\t', 1)
        R.require(len(fields) == 2 and fields[0] not in values, 'Ambiguous merged-cache metadata')
        values[fields[0]] = fields[1]
    expected = dict(assembly='GRCh38', source_assembly='GRCh38.p14', species='homo_sapiens',
                    source_refseq=RELEASE + ' - ' + PREFIX + '_genomic.gff',
                    bam=PREFIX + '_knownrefseq_alns.bam')
    R.require(all(values.get(k) == v for k, v in expected.items()),
              'Configured merged VEP116 cache source differs from the reviewed RefSeq release; no fallback')
    R.require(identity(path) == before, 'Merged-cache metadata changed while reading')
    return root, dict(path=str(path), identity=before, sha256=hashlib.sha256(data).hexdigest(), fields=expected)


def plan(config):
    root, source = cache_source(config)
    directory = root / 'references/RefSeq/GRCh38' / RELEASE
    cursor = root
    for part in directory.relative_to(root).parts:
        cursor /= part
        R.require(not cursor.is_symlink(), f'Refuse a linked candidate directory: {cursor}')
    return dict(schema_version=1, kind='cadd-refseq-acquisition-plan-v1', assembly='GRCh38',
                vep_release='116', refseq_release=RELEASE, cache_root=str(root),
                directory=str(directory), merged_cache=source,
                manifest_url=BASE + 'md5checksums.txt', resources={
                    name: dict(path=str(directory / filename), source_url=BASE + filename,
                               bytes=size, publisher_md5=md5)
                    for name, (filename, size, md5) in FILES.items()}, **FLAGS)


def http_metadata(response, url):
    R.require(response.geturl() == url and response.status == 200,
              'Publisher response redirected or changed status; no source fallback')
    encoding = response.headers.get('Content-Encoding')
    R.require(encoding in (None, '', 'identity'), 'Publisher ignored Accept-Encoding: identity')
    length = response.headers.get('Content-Length', '')
    R.require(length.isdigit() and int(length) > 0, 'Publisher entity length is unavailable')
    return dict(url=url, requested_accept_encoding='identity', content_encoding=encoding,
                content_length=int(length), etag=response.headers.get('ETag'),
                last_modified=response.headers.get('Last-Modified'),
                entity_bytes_policy='Raw HTTP entity bytes; no content decoding or decompression')


def publisher_pins(selected, timeout=60):
    """Small metadata GET plus HEADs only; the manifest SHA is observed, not assumed."""
    url = selected['manifest_url']
    headers = {'Accept-Encoding': 'identity', 'User-Agent': 'PriVA-RefSeq-candidate/1.0'}
    with urlopen(Request(url, headers=headers), timeout=timeout) as response:
        http = http_metadata(response, url)
        R.require(http['content_length'] <= MAX_METADATA, 'Oversized publisher MD5 manifest')
        data = response.read(MAX_METADATA + 1)
    R.require(len(data) == http['content_length'], 'Publisher MD5 manifest length differs')
    rows = [line.split() for line in data.decode('ascii').splitlines() if line.strip()]
    heads = {}
    for name, resource in selected['resources'].items():
        filename = Path(resource['path']).name
        matches = [r for r in rows if len(r) == 2 and r[1] == './' + filename]
        R.require(len(matches) == 1 and matches[0][0] == resource['publisher_md5'] and
                  re.fullmatch('[0-9a-f]{32}', matches[0][0]),
                  f'Publisher MD5 differs or is missing/ambiguous: {filename}')
        with urlopen(Request(resource['source_url'], headers=headers, method='HEAD'), timeout=timeout) as response:
            head = http_metadata(response, resource['source_url'])
        R.require(head['content_length'] == resource['bytes'], f'Publisher length differs: {filename}')
        heads[name] = head
    return data, dict(sha256=hashlib.sha256(data).hexdigest(), bytes=len(data), http=http), heads


def verify_file(path, resource, expected_sha=None):
    R.require(not path.is_symlink() and path.is_file(), f'Not a regular candidate file: {path}')
    before = identity(path)
    R.require(before['size'] == resource['bytes'], f'Candidate length differs; retained: {path}')
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(R.BLOCK), b''):
            md5.update(block); sha.update(block)
    R.require(identity(path) == before, f'Candidate changed while verifying: {path}')
    R.require(md5.hexdigest() == resource['publisher_md5'] and
              (expected_sha is None or sha.hexdigest() == expected_sha),
              f'Candidate checksum differs; retained: {path}')
    return sha.hexdigest(), before


def publish_bytes(path, data):
    if os.path.lexists(path):
        R.require(not path.is_symlink() and path.is_file() and path.stat().st_size == len(data) and path.read_bytes() == data,
                  f'Existing publisher manifest differs; retained: {path}')
        return
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix='.refseq-manifest-', delete=False) as stream:
        temporary = Path(stream.name); stream.write(data)
    try: R.publish_new(temporary, path)
    finally: temporary.unlink(missing_ok=True)


def source_record(selected, directory, timeout):
    receipt = directory / 'source.ready.json'
    if receipt.exists():
        source = read_json(receipt)
        R.require(source['plan'] == selected, 'RefSeq acquisition source/cache identity changed; retained')
        manifest = source['publisher_md5_manifest']; path = Path(manifest['path'])
        R.require(path == directory / 'md5checksums.txt' and not path.is_symlink() and
                  path.stat().st_size == manifest['bytes'] and R.digest(path) == manifest['sha256'],
                  'Frozen publisher MD5 manifest changed; retained')
        return source
    data, manifest, heads = publisher_pins(selected, timeout)
    R.require(cache_source_from_plan(selected), 'Merged-cache source changed before publication')
    manifest['path'] = str(directory / 'md5checksums.txt')
    publish_bytes(Path(manifest['path']), data)
    source = dict(schema_version=1, kind='cadd-refseq-acquisition-source-v1', plan=selected,
                  publisher_md5_manifest=manifest, heads=heads, **FLAGS)
    publish_json(receipt, source)
    return source


def cache_source_from_plan(selected):
    record = selected['merged_cache']; path = Path(record['path'])
    return identity(path) == record['identity'] and R.digest(path) == record['sha256']


def prefix_sha(path, size):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        while size:
            block = stream.read(min(size, R.BLOCK))
            R.require(block, 'Owned partial shrank; retained')
            h.update(block); size -= len(block)
    return h.hexdigest()


def replace_state(path, value, previous):
    """Only this locked acquisition's owned checkpoint may be replaced."""
    R.require(not path.is_symlink() and
              (identity(path) if path.exists() else None) == previous, 'Download checkpoint changed; retained')
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, prefix='.refseq-state-', delete=False) as stream:
        temporary = Path(stream.name)
        json.dump(value, stream, sort_keys=True, indent=2); stream.write('\n')
        stream.flush(); os.fsync(stream.fileno())
    try:
        if previous is None: R.publish_new(temporary, path)
        else:
            R.require(identity(path) == previous, 'Download checkpoint changed before publication')
            os.replace(temporary, path)
    finally: temporary.unlink(missing_ok=True)


def owned_partial(path, state_path, scope, size):
    """Accept exact checkpoints or same-inode crash growth with a verified prefix."""
    if state_path.exists():
        state = read_json(state_path)
        R.require(state.get('scope') == scope and state.get('kind') == 'cadd-refseq-owned-partial-v1',
                  'Partial source ownership differs; retained')
        R.require(path.exists() and not path.is_symlink(), 'Owned partial missing/replaced; retained')
        current, old = identity(path), state['identity']
        same = all(current[k] == old[k] for k in ('path', 'target', 'inode', 'mode'))
        growth = old['size'] < current['size'] <= size and current['mtime_ns'] >= old['mtime_ns'] and current['ctime_ns'] >= old['ctime_ns']
        R.require(same and (current == old or growth), 'Owned partial identity changed; retained')
        R.require(prefix_sha(path, old['size']) == state['prefix_sha256'] and identity(path) == current,
                  'Owned partial prefix changed; retained')
        return state
    R.require(not os.path.lexists(path) and not os.path.lexists(state_path), 'Unowned partial/checkpoint retained')
    with path.open('xb'): pass
    state = dict(kind='cadd-refseq-owned-partial-v1', scope=scope, identity=identity(path),
                 prefix_sha256=hashlib.sha256(b'').hexdigest())
    replace_state(state_path, state, None)
    return state


def checkpoint(path, state_path, state):
    previous = identity(state_path); before = identity(path)
    updated = {**state, 'identity': before, 'prefix_sha256': prefix_sha(path, before['size'])}
    R.require(identity(path) == before, 'Partial changed while checkpointing; retained')
    replace_state(state_path, updated, previous)
    return updated


def transfer(partial, state_path, state, resource, head, timeout):
    offset = identity(partial)['size']; remaining = resource['bytes'] - offset
    R.require(remaining >= 0, 'Owned partial exceeds pinned length; retained')
    if not remaining: return
    headers = {'Accept-Encoding': 'identity', 'User-Agent': 'PriVA-RefSeq-candidate/1.0'}
    if offset: headers['Range'] = f'bytes={offset}-'
    if head.get('etag'): headers['If-Match'] = head['etag']
    elif head.get('last_modified'): headers['If-Unmodified-Since'] = head['last_modified']
    with urlopen(Request(resource['source_url'], headers=headers), timeout=timeout) as response:
        R.require(response.geturl() == resource['source_url'], 'Payload redirected; no source fallback')
        R.require(response.headers.get('Content-Encoding') in (None, '', 'identity'), 'Payload content encoding differs')
        if offset or response.status == 206:
            R.require(response.status == 206 and response.headers.get('Content-Range') ==
                      f"bytes {offset}-{resource['bytes'] - 1}/{resource['bytes']}", 'Publisher refused exact resume; partial retained')
        else: R.require(response.status == 200, 'Unexpected payload HTTP status')
        R.require(response.headers.get('Content-Length') == str(remaining), 'Payload response length differs')
        for field, key in (('ETag', 'etag'), ('Last-Modified', 'last_modified')):
            value = response.headers.get(field)
            R.require(not (value and head.get(key) and value != head[key]), 'Payload source identity changed')
        before = identity(partial)
        fd = os.open(partial, os.O_WRONLY | os.O_APPEND | os.O_NOFOLLOW)
        try:
            stat = os.fstat(fd)
            R.require(stat.st_ino == before['inode'] and stat.st_size == offset, 'Partial changed before append')
            with os.fdopen(fd, 'ab', closefd=False) as stream:
                try:
                    while remaining:
                        block = response.read(min(R.BLOCK, remaining + 1))
                        if not block: break
                        R.require(len(block) <= remaining, 'Payload exceeded pinned length; partial retained')
                        stream.write(block); remaining -= len(block)
                finally:
                    stream.flush(); os.fsync(stream.fileno())
        finally:
            os.close(fd)
            R.require(identity(partial)['inode'] == before['inode'], 'Partial replaced during transfer')
            checkpoint(partial, state_path, state)
    R.require(remaining == 0, 'Payload ended before pinned length; owned partial retained')


def resource_acquire(selected, name, source, timeout):
    resource = selected['resources'][name]; path = Path(resource['path']); directory = path.parent
    receipt = directory / (name + '.ready.json'); source_sha = fingerprint(source)
    if receipt.exists(): return resource_check(selected, name, source)
    state_path = directory / (name + '.partial.json'); partial = R.staging(str(path) + '.partial')
    scope = dict(source_fingerprint=source_sha, component=name, resource=resource)
    if os.path.lexists(path):
        # Recover publication interrupted before its ready receipt, only with
        # the same owned inode. Identical but unidentified final bytes are held.
        R.require(not path.is_symlink() and state_path.is_file(), f'Unidentified existing final file retained: {path}')
        state = read_json(state_path)
        R.require(state.get('scope') == scope and state.get('kind') == 'cadd-refseq-owned-partial-v1' and
                  identity(path)['inode'] == state['identity']['inode'], 'Final-file ownership differs; retained')
        if partial.exists():
            R.require(not partial.is_symlink() and partial.samefile(path), 'Published staging ownership differs; retained')
            partial.unlink()
    else:
        state = owned_partial(partial, state_path, scope, resource['bytes'])
        transfer(partial, state_path, state, resource, source['heads'][name], timeout)
        verify_file(partial, resource)
        R.publish_new(partial, path)
    sha, current = verify_file(path, resource)
    record = dict(schema_version=1, kind='cadd-refseq-acquired-resource-v1', component=name,
                  **resource, sha256=sha, identity=current, source_fingerprint=source_sha, **FLAGS)
    publish_json(receipt, record)
    return record


def resource_check(selected, name, source):
    resource = selected['resources'][name]; path = Path(resource['path'])
    record = read_json(path.parent / (name + '.ready.json'))
    expected = dict(schema_version=1, kind='cadd-refseq-acquired-resource-v1', component=name,
                    **resource, source_fingerprint=fingerprint(source), **FLAGS)
    R.require(all(record.get(k) == v for k, v in expected.items()), 'Resource receipt scope differs; retained')
    sha, current = verify_file(path, resource, record['sha256'])
    R.require(current == record['identity'], 'Acquired resource identity changed; retained')
    return record


def acquisition_record(selected, source, resources):
    return dict(schema_version=1, kind='cadd-refseq-acquisition-v1', plan=selected,
                source_fingerprint=fingerprint(source), resources=resources,
                acceptance='Source acquisition only; sequence/coordinate/feature validation is separate', **FLAGS)


def acquire(config, timeout=60):
    selected = plan(config)
    with T.locked(selected) as directory:
        source = source_record(selected, directory, timeout)
        resources = {name: resource_acquire(selected, name, source, timeout) for name in selected['resources']}
        R.require(cache_source_from_plan(selected), 'Merged-cache source changed during acquisition')
        result = acquisition_record(selected, source, resources)
        publish_json(directory / 'acquisition.ready.json', result)
    return result


def check(config):
    selected = plan(config); directory = Path(selected['directory'])
    # No locked()/mkdir, metadata rewrite, network access or repair in check.
    source = read_json(directory / 'source.ready.json')
    R.require(source['plan'] == selected, 'RefSeq acquisition source/cache identity changed; retained')
    manifest = source['publisher_md5_manifest']; path = Path(manifest['path'])
    R.require(path == directory / 'md5checksums.txt' and not path.is_symlink() and
              path.stat().st_size == manifest['bytes'] and R.digest(path) == manifest['sha256'],
              'Frozen publisher MD5 manifest changed; retained')
    resources = {name: resource_check(selected, name, source) for name in selected['resources']}
    result = acquisition_record(selected, source, resources)
    R.require(read_json(directory / 'acquisition.ready.json') == result, 'Aggregate acquisition receipt differs')
    R.require(cache_source_from_plan(selected), 'Merged-cache source changed during check')
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='procedure', required=True)
    for name in ('plan', 'acquire', 'check'):
        child = commands.add_parser(name)
        inputs = child.add_mutually_exclusive_group(required=True)
        inputs.add_argument('--config', type=Path); inputs.add_argument('--config-json')
        if name == 'acquire': child.add_argument('--timeout-seconds', type=int, default=60)
    args = parser.parse_args(argv)
    try:
        if args.config_json: config = json.loads(args.config_json)
        else:
            import yaml
            config = yaml.safe_load(args.config.read_text())
        if args.procedure == 'acquire':
            R.require(0 < args.timeout_seconds <= 600, 'Timeout must be from 1 to 600 seconds')
            result = acquire(config, args.timeout_seconds)
        else: result = globals()[args.procedure](config)
        print(json.dumps(result, sort_keys=True)); return 0
    except (ValueError, OSError, KeyError, TypeError) as error:
        print('CADD RefSeq candidate: ' + str(error), file=sys.stderr); return 1


if __name__ == '__main__': raise SystemExit(main())
