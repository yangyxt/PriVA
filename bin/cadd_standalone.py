#!/usr/bin/env python3
"""Explicit, hash-pinned bridge to the standalone original-Python CADD interface.

The two artifact receipts are intentionally not interchanged. A native import
retains identical annotations; this bridge records both receipt identities.
"""
import json
from pathlib import Path
import shutil
import subprocess

from cadd_shared_vep import sha256, json_file

IMPLEMENTATION = 'configured complete Python calculation; native postprocessing disabled'


class StandaloneError(ValueError):
    pass


def verify(settings):
    for key in ['program', 'sha256', 'source_root', 'revision']:
        if not settings.get(key):
            raise StandaloneError('Standalone CADD pin lacks ' + key)
    program = Path(settings['program'])
    source = Path(settings['source_root'])
    if not program.is_absolute() or not source.is_absolute() or sha256(program) != settings['sha256']:
        raise StandaloneError('Standalone CADD executable checksum/path differs')
    revision = subprocess.run(['git', '-C', str(source), 'rev-parse', 'HEAD'],
                              capture_output=True, text=True, check=True).stdout.strip()
    if revision != settings['revision']:
        raise StandaloneError('Standalone CADD source revision differs')
    subprocess.run(['git', '-C', str(source), 'diff', '--quiet', 'HEAD', '--', '.'], check=True)
    return program, source


def projected_profile(profile, source):
    calculation = profile['calculation']
    env = calculation['env']
    phred = Path(env['CADD_PHRED_TABLE']).resolve(strict=True)
    matches = [resource for resource in calculation['verify_files']
               if Path(resource['path']).resolve(strict=True) == phred]
    if len(matches) != 1:
        raise StandaloneError('Exactly one verified PHRED resource is required')
    return {'schema_version': 1, 'profile_id': profile['profile_id'], 'contract': profile['contract'],
            'calculation': {'backend': {'program': str(source / 'reference/calculate.sh'),
                            'env': env, 'verify_files': calculation['verify_files']},
                            'phred_table': matches[0]}}


def execute(profile, annotation_bundle, output, work):
    settings = profile['calculation']['standalone']
    program, source = verify(settings)
    annotation_bundle = Path(annotation_bundle).resolve(strict=True)
    original_manifest = annotation_bundle / 'manifest.json'
    original_table = annotation_bundle / 'annotations.tsv.gz'
    original_hash = sha256(original_manifest)
    original_table_hash = sha256(original_table)
    profile_path = Path(work) / 'standalone-profile.json'
    with open(profile_path, 'x') as handle:
        json.dump(projected_profile(profile, source), handle, indent=2, sort_keys=True)
        handle.write('\n')
    native_annotations = Path(work) / 'standalone-annotations'
    native_scores = Path(work) / 'standalone-scores'
    subprocess.run([str(program), 'import-annotations', '--profile', str(profile_path),
                    '--input', str(original_table), '--output', str(native_annotations)], check=True)
    imported = json_file(native_annotations / 'manifest.json')
    native_table = native_annotations / 'annotations.tsv.gz'
    if imported.get('contract') != profile['contract'] or imported.get('input_sha256') != original_table_hash:
        raise StandaloneError('Standalone import lost the verified input contract')
    if imported.get('output_sha256') != sha256(native_table) or imported.get('engine_sha256') != settings['sha256']:
        raise StandaloneError('Standalone import provenance does not match the pinned executable/artifact')
    subprocess.run([str(program), 'calculate', '--profile', str(profile_path),
                    '--annotations', str(native_annotations), '--output', str(native_scores)], check=True)
    score_path = native_scores / 'scores.tsv'
    scored = json_file(native_scores / 'manifest.json')
    expected = {'stage': 'calculation', 'implementation': IMPLEMENTATION,
                'engine_sha256': settings['sha256'], 'profile_sha256': sha256(profile_path),
                'profile_id': profile['profile_id'], 'contract': profile['contract'],
                'input_sha256': sha256(native_table),
                'input_manifest_sha256': sha256(native_annotations / 'manifest.json'),
                'output_sha256': sha256(score_path)}
    if any(scored.get(key) != value for key, value in expected.items()):
        raise StandaloneError('Standalone calculation receipt differs from the pinned all-Python contract')
    verify(settings)
    if sha256(original_manifest) != original_hash or sha256(original_table) != original_table_hash:
        raise StandaloneError('Original annotation bundle changed during standalone execution')
    with open(score_path, 'rb') as source_file, open(output, 'xb') as destination:
        shutil.copyfileobj(source_file, destination)
    return {'program_sha256': settings['sha256'], 'revision': settings['revision'],
            'projected_profile_sha256': sha256(profile_path),
            'original_annotation_manifest_sha256': original_hash,
            'native_annotation_manifest_sha256': sha256(native_annotations / 'manifest.json'),
            'native_calculation_manifest_sha256': sha256(native_scores / 'manifest.json'),
            'implementation': IMPLEMENTATION, 'profiling': scored.get('profiling')}
