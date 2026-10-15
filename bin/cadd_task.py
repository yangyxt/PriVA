#!/usr/bin/env python3
"""Offline CADD task contracts; scientific calculations stay in upstream Python.

The JSON profiles are trusted installation configuration, not untrusted inputs.
Resource locking is a separate installation action. Analysis never installs data.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
PIN = '1fead0c1e02d3f0a8599bdb04e6b3a2dd70b4866'
ALIASES = {'hg19': 'GRCh37', 'GRCh37': 'GRCh37', 'hg38': 'GRCh38', 'GRCh38': 'GRCh38'}


class ContractError(ValueError):
    pass


def digest(path):
    value = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def read_json(path):
    with open(path) as handle:
        return json.load(handle)


def write_json(path, value):
    with open(path, 'x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write('\n')


def source_identity():
    paths = ['bin/cadd_task.py', 'bin/cadd/upstream_annotate.sh',
             'bin/cadd/upstream_calculate.sh', 'cadd.nf',
             'modules/local/cadd_upstream.nf', 'conf/cadd.config', 'nextflow.config',
             'conf/cadd/GRCh37.reference.json', 'conf/cadd/GRCh38.reference.json']
    return {name: digest(ROOT / name) for name in paths}


def load_profile(path, assembly=None, stub=False, locking=False):
    p = read_json(path)
    if p.get('schema_version') != 1 or p.get('upstream_commit') != PIN:
        raise ContractError('Profile must pin the audited CADD v1.7.2 source revision')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', p.get('profile_id', '')):
        raise ContractError('A safe nonempty profile_id is required')
    if p.get('profile_kind') not in ['reference', 'adapted']:
        raise ContractError('profile_kind must be reference or adapted')
    if not p.get('software_id') or (not stub and p['software_id'].startswith('REPLACE_')):
        raise ContractError('software_id must identify the installed pinned environment/image')
    if bool(p.get('stub', False)) != bool(stub):
        raise ContractError('Stub profiles and scientific execution must not be mixed')
    c = p['contract']
    if c['assembly'] not in ['GRCh37', 'GRCh38']:
        raise ContractError('Profile assembly must be GRCh37 or GRCh38')
    if assembly and ALIASES[assembly] != c['assembly']:
        raise ContractError('Requested assembly differs from the profile assembly')
    canonical = read_json(ROOT / 'conf/cadd' / (c['assembly'] + '.reference.json'))['contract']['columns']
    if c.get('columns') != canonical:
        raise ContractError('Complete ordered CADD 1.7 annotation columns are required')
    if not c.get('schema_id'):
        raise ContractError('Annotation schema identity is required')
    if not locking and not re.fullmatch(r'[0-9a-f]{64}', c.get('resource_manifest_sha256', '')):
        raise ContractError('A locked resource manifest SHA256 is required')
    if p['profile_kind'] == 'reference' and 'annotation' in p and p['annotation']['env'].get('CADD_ENSEMBL_DB') != '110':
        raise ContractError('Updated VEP annotations require an adapted profile, not a reference label')
    return p


def record(path, checksum=False):
    path = Path(path).resolve(strict=True)
    stat = path.stat()
    if not path.is_file():
        raise ContractError('Resource is not a regular file: ' + str(path))
    item = {'path': str(path), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}
    if checksum:
        item['sha256'] = digest(path)
    return item


def inventory(roots, checksum=False):
    if not roots:
        raise ContractError('Each enabled stage needs explicit nonempty resource_roots')
    paths = set()
    for raw in roots:
        if not Path(raw).is_absolute():
            raise ContractError('Resource paths must be absolute')
        root = Path(raw).resolve(strict=True)
        if root.is_file():
            paths.add(root)
        elif root.is_dir():
            for directory, subdirs, files in os.walk(root):
                for name in subdirs:
                    child = Path(directory) / name
                    if child.is_symlink():
                        raise ContractError('Declare linked resource directories as separate resolved roots: ' + str(child))
                for name in files:
                    paths.add((Path(directory) / name).resolve(strict=True))
        else:
            raise ContractError('Invalid resource root: ' + str(root))
    if not paths:
        raise ContractError('Resource inventory is empty')
    return {str(path): record(path, checksum) for path in sorted(paths)}


def lock_resources(args):
    p = load_profile(args.profile, locking=True)
    resources = {}
    for name in ['annotation', 'calculation']:
        if name in p:
            resources[name] = inventory(p[name].get('resource_roots'), args.checksum)
    if not resources:
        raise ContractError('No enabled stages in profile')
    write_json(args.output, {'schema_version': 1, 'assembly': p['contract']['assembly'],
                            'resources': resources})
    print(json.dumps({'resource_manifest_path': str(Path(args.output).resolve()),
                      'resource_manifest_sha256': digest(args.output)}))


def verify_resources(p, stage):
    path = Path(p['resource_manifest_path'])
    if not path.is_absolute() or digest(path) != p['contract']['resource_manifest_sha256']:
        raise ContractError('Locked resource manifest checksum differs')
    m = read_json(path)
    if m.get('schema_version') != 1 or m.get('assembly') != p['contract']['assembly']:
        raise ContractError('Resource manifest schema/assembly differs')
    if stage not in m.get('resources', {}):
        raise ContractError('Resource manifest lacks stage: ' + stage)
    expected = m['resources'][stage]
    current = inventory(p[stage].get('resource_roots'))
    if set(current) != set(expected):
        raise ContractError('Resource inventory membership changed; relock during installation')
    for name, actual in current.items():
        wanted = expected[name]
        if actual != {k: v for k, v in wanted.items() if k != 'sha256'}:
            raise ContractError('Resource changed since installation: ' + name)
        if 'sha256' in wanted and digest(name) != wanted['sha256']:
            raise ContractError('Resource content changed since installation: ' + name)
    env = p[stage]['env']
    # Every directly configured data path must be covered by the locked inventory.
    for key, value in env.items():
        if key in ['CADD_SOURCE_ROOT', 'CADD_ASSEMBLY', 'CADD_ENSEMBL_DB', 'CADD_ESM_MODELS',
                   'CADD_ESM_BATCH_SIZE', 'CADD_PYTHON2', 'CADD_VEP', 'CADD_ESM_PYTHON',
                   'CADD_REGSEQ_PYTHON', 'CADD_MMSPLICE_PYTHON', 'CADD_BGZIP', 'CADD_TABIX']:
            continue
        candidate = Path(value).resolve(strict=True)
        if candidate.is_file() and str(candidate) not in current:
            raise ContractError('Configured data file is outside locked resources: ' + key)
        if candidate.is_dir() and not any(candidate in Path(name).parents for name in current):
            raise ContractError('Configured data directory is outside locked resources: ' + key)
    return expected


def verify_source(env):
    source = Path(env['CADD_SOURCE_ROOT'])
    if not source.is_absolute():
        raise ContractError('CADD_SOURCE_ROOT must be absolute')
    result = subprocess.run(['git', '-C', str(source), 'rev-parse', 'HEAD'],
                            check=True, capture_output=True, text=True)
    if result.stdout.strip() != PIN:
        raise ContractError('Installed CADD source differs from pinned v1.7.2')
    subprocess.run(['git', '-C', str(source), 'diff', '--quiet', 'HEAD', '--', 'src', 'config'], check=True)
    # Untracked modules can change Python imports even when tracked files are clean.
    untracked = subprocess.run(['git', '-C', str(source), 'ls-files', '--others', '--exclude-standard', '--', 'src', 'config'],
                               check=True, capture_output=True, text=True).stdout
    if untracked.strip():
        raise ContractError('Untracked files in CADD src/config prevent source qualification')


def stage_preflight(p, stage, stub=False):
    if stub:
        return
    if stage not in p:
        raise ContractError('Profile does not enable ' + stage)
    cfg = p[stage]
    env = cfg['env']
    if any(not key.startswith('CADD_') or not isinstance(value, str) for key, value in env.items()):
        raise ContractError('Backend environment must contain only string CADD_* entries')
    verify_source(env)
    verify_resources(p, stage)
    verified = set()
    for asset in cfg.get('verify_files', []):
        path = Path(asset['path'])
        if not path.is_absolute() or not re.fullmatch(r'[0-9a-f]{64}', asset['sha256']) or digest(path) != asset['sha256']:
            raise ContractError('Immutable asset checksum differs: ' + str(path))
        verified.add(path.resolve())
    if stage == 'calculation':
        for key in ['CADD_IMPUTE_CONFIG', 'CADD_MODEL', 'CADD_PHRED_TABLE']:
            if Path(env[key]).resolve() not in verified:
                raise ContractError(key + ' must be bound to a verified immutable asset')
    for key in ['CADD_PYTHON2', 'CADD_VEP', 'CADD_ESM_PYTHON', 'CADD_REGSEQ_PYTHON',
                'CADD_MMSPLICE_PYTHON', 'CADD_BGZIP', 'CADD_TABIX']:
        if key in env and (not Path(env[key]).is_absolute() or not os.access(env[key], os.X_OK)):
            raise ContractError('Configured executable must be an installed absolute path: ' + key)


def open_table(path):
    with open(path, 'rb') as handle:
        compressed = handle.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt', newline='') if compressed else open(path, newline='')


def rows(path, expected):
    with open_table(path) as handle:
        seen = False
        for line in handle:
            text = line.rstrip('\r\n')
            if not seen and text.startswith('##'):
                continue
            if not seen:
                if text.lstrip('#').split('\t') != expected:
                    raise ContractError('Table header differs from complete ordered annotation contract')
                seen = True
            else:
                fields = text.split('\t')
                if text.startswith('#') or len(fields) != len(expected):
                    raise ContractError('Invalid table row width or embedded header')
                yield fields
        if not seen:
            raise ContractError('Table has no header')


def canonicalize(source, destination, columns):
    count = 0
    # Remove only pre-header metadata; upstream predictSKmodel misinterprets a
    # metadata line as its primary header. Annotation values/order are unchanged.
    with gzip.open(destination, 'xt', newline='') as output:
        output.write('#' + '\t'.join(columns) + '\n')
        for fields in rows(source, columns):
            output.write('\t'.join(fields) + '\n')
            count += 1
    return count


def verify_bundle(bundle, p, stub=False):
    bundle = Path(bundle)
    receipt = read_json(bundle / 'manifest.json')
    if receipt.get('schema_version') != 1 or receipt.get('stage') not in ['annotation', 'import-annotations']:
        raise ContractError('Expected a completed annotation bundle')
    if bool(receipt.get('stub')) != bool(stub):
        raise ContractError('Stub and scientific annotation artifacts must not be mixed')
    for key, value in [('contract', p['contract']), ('profile_id', p['profile_id']),
                       ('software_id', p['software_id']), ('profile_kind', p['profile_kind']), ('upstream_commit', PIN)]:
        if receipt.get(key) != value:
            raise ContractError('Annotation bundle/profile mismatch: ' + key)
    table = bundle / 'annotations.tsv.gz'
    if digest(table) != receipt['output_sha256']:
        raise ContractError('Annotation bundle content changed')
    count = sum(1 for _ in rows(table, p['contract']['columns']))
    if count != receipt['row_count']:
        raise ContractError('Annotation row count differs from receipt')
    return table, receipt


def run_backend(p, stage, source, output, work):
    cfg = p[stage]
    env = {k: v for k, v in os.environ.items() if not k.startswith('CADD_')}
    env.update(cfg['env'])
    env.pop('PYTHONPATH', None)
    env.pop('PYTHONHOME', None)
    env.update(CADD_ASSEMBLY=p['contract']['assembly'], PYTHONDONTWRITEBYTECODE='1', PYTHONNOUSERSITE='1')
    script = ROOT / 'bin/cadd' / ('upstream_' + ('annotate' if stage == 'annotation' else 'calculate') + '.sh')
    subprocess.run(['bash', str(script), str(Path(source).resolve()), str(Path(output).resolve()),
                    str(Path(work).resolve())], env=env, check=True)


def receipt(p, profile, stage, source, output, count, stub, before):
    if source_identity() != before:
        raise ContractError('CADD adapter source changed during task')
    return {'schema_version': 1, 'stage': stage, 'stub': bool(stub),
            'profile_id': p['profile_id'], 'profile_kind': p['profile_kind'],
            'profile_sha256': digest(profile), 'contract': p['contract'],
            'software_id': p['software_id'], 'upstream_commit': PIN,
            'adapter_source': before, 'input_sha256': digest(source),
            'output_sha256': digest(output), 'row_count': count,
            'implementation': ('all-four-calculation-stages-upstream-python' if stage == 'calculation' else
                               'upstream-python-annotation' if stage == 'annotation' else 'declared-complete-table-import')}


def execute(args):
    p = load_profile(args.profile, args.assembly, args.stub)
    before = source_identity()
    profile_before = digest(args.profile)
    stage = 'calculation' if args.command == 'calculate' else 'annotation'
    stage_preflight(p, stage, args.stub)
    out = Path(args.output)
    if out.exists():
        raise ContractError('Output exists; use a fresh task directory')
    out.mkdir(parents=True)
    work = out.parent / (out.name + '.work')
    work.mkdir()
    columns = p['contract']['columns']
    source = Path(args.input)
    if args.command == 'calculate':
        source, _ = verify_bundle(source, p, args.stub)
        input_before = digest(source)
        scored = out / 'scores.tsv'
        if args.stub:
            with open(scored, 'x') as handle:
                handle.write('#' + '\t'.join(columns + ['RawScore', 'PHRED']) + '\n')
                for fields in rows(source, columns):
                    handle.write('\t'.join(fields + ['0.000000', '0']) + '\n')
        else:
            run_backend(p, stage, source, scored, work)
        count = 0
        original_rows = rows(source, columns)
        for fields in rows(scored, columns + ['RawScore', 'PHRED']):
            original = next(original_rows, None)
            if original != fields[:-2]:
                raise ContractError('Calculation changed annotation values or row order')
            count += 1
        if next(original_rows, None) is not None:
            raise ContractError('Calculation dropped annotation rows')
        result = receipt(p, args.profile, stage, source, scored, count, args.stub, before)
        result['input_manifest_sha256'] = digest(Path(args.input) / 'manifest.json')
    else:
        input_before = digest(source)
        raw = source
        if args.command == 'annotate':
            raw = work / 'raw.anno.tsv.gz'
            if args.stub:
                fields = ['NA'] * len(columns)
                for name, value in [('Chrom', '1'), ('Pos', '100'), ('Ref', 'A'), ('Alt', 'G'),
                                    ('Type', 'SNV'), ('Consequence', 'NON_SYNONYMOUS')]:
                    fields[columns.index(name)] = value
                with gzip.open(raw, 'wt') as handle:
                    handle.write('#' + '\t'.join(columns) + '\n' + '\t'.join(fields) + '\n')
            else:
                run_backend(p, stage, source, raw, work)
        output = out / 'annotations.tsv.gz'
        count = canonicalize(raw, output, columns)
        result = receipt(p, args.profile, stage if args.command == 'annotate' else args.command,
                         source, output, count, args.stub, before)
    stage_preflight(p, stage, args.stub)
    if digest(source) != input_before or digest(args.profile) != profile_before:
        raise ContractError('Input or profile changed during execution')
    # A manifest is the completion marker; failed work retains diagnostics only.
    write_json(out / 'manifest.json', result)


def preflight(args):
    p = load_profile(args.profile, args.assembly, args.stub)
    stages = ['annotation', 'calculation'] if args.mode == 'full' else [
        'annotation' if args.mode == 'annotate' else 'calculation']
    for stage in stages:
        stage_preflight(p, stage, args.stub)
    source = Path(args.input)
    if args.mode == 'calculate':
        source, _ = verify_bundle(source, p, args.stub)
    identity = {'schema_version': 1, 'mode': args.mode, 'stub': bool(args.stub),
                'input_sha256': digest(source), 'profile_sha256': digest(args.profile),
                'contract': p['contract'], 'software_id': p['software_id'],
                'adapter_source': source_identity()}
    if args.mode == 'calculate':
        identity['input_manifest_sha256'] = digest(Path(args.input) / 'manifest.json')
    root = Path(args.output)
    if root.exists() and (root / 'launch.json').is_file():
        if read_json(root / 'launch.json') != identity:
            raise ContractError('Output belongs to different inputs/profile/source; use a new --outdir')
    else:
        root.mkdir(parents=True, exist_ok=True)
        # Nextflow creates its report directory before the workflow preflight.
        if any(child.name != 'pipeline_info' for child in root.iterdir()):
            raise ContractError('Output directory is not empty and has no matching launch receipt')
        write_json(root / 'launch.json', identity)
    print(json.dumps(identity, sort_keys=True))


def parser():
    root = argparse.ArgumentParser(description=__doc__)
    sub = root.add_subparsers(dest='command', required=True)
    lock = sub.add_parser('lock-resources', help='Installation only: inventory prepared resource roots; never download')
    lock.add_argument('--profile', required=True)
    lock.add_argument('--output', required=True)
    lock.add_argument('--checksum', action='store_true', help='Hash every resource; potentially TB of I/O')
    lock.set_defaults(func=lock_resources)
    for name in ['annotate', 'import-annotations', 'calculate', 'preflight']:
        cmd = sub.add_parser(name)
        cmd.add_argument('--profile', required=True)
        cmd.add_argument('--assembly', required=True, choices=ALIASES)
        cmd.add_argument('--input', required=True, help='VCF, complete table, or annotation bundle according to command')
        cmd.add_argument('--output', required=True)
        cmd.add_argument('--stub', action='store_true')
        if name == 'preflight':
            cmd.add_argument('--mode', choices=['annotate', 'calculate', 'full'], required=True)
        cmd.set_defaults(func=preflight if name == 'preflight' else execute)
    return root


def main():
    args = parser().parse_args()
    try:
        args.func(args)
    except (ContractError, OSError, KeyError, TypeError, json.JSONDecodeError, subprocess.CalledProcessError) as error:
        print('PriVA CADD contract: ' + str(error), file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
