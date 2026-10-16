#!/usr/bin/env python3
"""Recover a qualified CADD view from native ordinals in PriVA's one VEP run.

This selects existing CSQ rows in the order emitted by the installed official
VEP per_gene method. It never recreates ranking rules or reruns VEP.
"""
import argparse
import gzip
import json
from pathlib import Path
import re
import shutil
import tempfile

import cadd_shared_vep as S

ADAPTER = 'official-vep-native-per-gene-v1'


def select_rows(source, destination, context_sha256):
    columns = None
    variants = 0
    with S.text_file(source) as handle, gzip.open(destination, 'wt') as output:
        for line in handle:
            if line.startswith('##INFO=<ID=CSQ,'):
                match = re.search(r'Format: ([^"\n]+)', line)
                if not match:
                    raise S.SharedVepError('CSQ header lacks Format declaration')
                columns = match.group(1).rstrip('>').split('|')
                required = S.REQUIRED_CSQ + ['CADD_PICK_ORDINAL', 'CADD_PICK_COUNT', 'CADD_PICK_CONTEXT']
                if len(set(columns)) != len(columns) or set(required) - set(columns):
                    raise S.SharedVepError('Native selection fields are missing or duplicated')
            if line.startswith('#'):
                output.write(line)
                continue
            if columns is None:
                raise S.SharedVepError('Shared VEP artifact lacks CSQ schema')
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 8:
                raise S.SharedVepError('Shared VEP variant has fewer than eight fields')
            info = fields[7].split(';')
            indices = [i for i, item in enumerate(info) if item.startswith('CSQ=')]
            if len(indices) != 1:
                raise S.SharedVepError('Every variant must contain one complete CSQ field')
            selected = {}
            counts = set()
            for annotation in info[indices[0]][4:].split(','):
                values = annotation.split('|')
                if len(values) != len(columns):
                    raise S.SharedVepError('CSQ row width differs from header')
                row = dict(zip(columns, values))
                if row['CADD_PICK_CONTEXT'] != context_sha256:
                    raise S.SharedVepError('Native selection runtime/source context differs')
                try:
                    ordinal, count = int(row['CADD_PICK_ORDINAL']), int(row['CADD_PICK_COUNT'])
                except ValueError:
                    raise S.SharedVepError('Missing or invalid native ordinal/count; inspect VEP plugin warnings')
                if count < 1 or ordinal < 0 or ordinal > count:
                    raise S.SharedVepError('Native ordinal/count is out of range')
                counts.add(count)
                if ordinal:
                    if ordinal in selected:
                        raise S.SharedVepError('Duplicate native selection ordinal')
                    selected[ordinal] = annotation
            if len(counts) != 1 or set(selected) != set(range(1, next(iter(counts)) + 1)):
                raise S.SharedVepError('Native per_gene row set is incomplete; a plugin may have filtered rows')
            info[indices[0]] = 'CSQ=' + ','.join(selected[i] for i in range(1, len(selected) + 1))
            fields[7] = ';'.join(info)
            output.write('\t'.join(fields) + '\n')
            variants += 1
    return variants


def build(source, profile, output, stub=False):
    source, output = Path(source).resolve(strict=True), Path(output).absolute()
    settings = profile['annotation']['shared_vep']
    if settings.get('adapter_id') != ADAPTER or not settings.get('qualification'):
        raise S.SharedVepError('Native ordinal selection requires explicit profile qualification')
    context_sha = S.native_context(profile)
    qualification = settings['qualification']
    evidence = S.json_file(qualification['path'])
    if S.sha256(qualification['path']) != qualification['sha256'] or evidence.get('selection_context_sha256') != context_sha:
        raise S.SharedVepError('Native selection context is not covered by the qualification')
    if output.exists():
        raise S.SharedVepError('View output already exists')
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=output.name + '.partial-', dir=output.parent))
    try:
        source_sha = S.sha256(source)
        vcf = temporary / 'vep.vcf.gz'
        variants = select_rows(source, vcf, context_sha)
        receipt = {'schema_version': 1, 'stage': 'cadd-shared-vep-view',
                   'assembly': profile['contract']['assembly'], 'profile_id': profile['profile_id'],
                   'resource_manifest_sha256': profile['contract']['resource_manifest_sha256'],
                   'vep_release': settings['vep_release'], 'transcript_release': settings['transcript_release'],
                   'selection': 'cadd-per-gene', 'adapter_id': ADAPTER, 'stub': bool(stub),
                   'qualification_sha256': qualification['sha256'],
                   'selection_context_sha256': context_sha, 'transcript_policy': settings['transcript_policy'],
                   'source_vep': {'path': str(source), 'sha256': source_sha},
                   'output_sha256': S.sha256(vcf), 'variant_count': variants,
                   'expected_annotation_rows': sum(1 for _ in S.expected_features(vcf))}
        (temporary / 'manifest.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
        S.verify(temporary, profile, stub)
        temporary.rename(output)
    except Exception:
        shutil.rmtree(temporary)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--profile', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    build(args.input, S.json_file(args.profile), args.output)


if __name__ == '__main__':
    main()
