#!/usr/bin/env python3
"""Fail-closed contract for a CADD-compatible view of PriVA's single VEP run.

This validates a previously qualified selection adapter. It does not invent a
replacement for VEP --per_gene or silently select/drop PriVA transcript rows.
"""
import gzip
import hashlib
import json
from pathlib import Path
import re

REQUIRED_CSQ = ['Allele', 'Consequence', 'Feature_type', 'Gene', 'Feature', 'SYMBOL',
                'CCDS', 'INTRON', 'EXON', 'cDNA_position', 'CDS_position', 'Protein_position',
                'Amino_acids', 'DOMAINS', 'SIFT', 'PolyPhen', 'MOTIF_NAME', 'HIGH_INF_POS',
                'MOTIF_SCORE_CHANGE']


class SharedVepError(ValueError):
    pass


def sha256(path):
    result = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def json_file(path):
    with open(path) as handle:
        return json.load(handle)


def text_file(path):
    with open(path, 'rb') as handle:
        compressed = handle.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if compressed else open(path)


def expected_features(path):
    """Yield exact variant+FeatureID order after upstream motif-row aggregation."""
    csq_columns = None
    header = False
    with text_file(path) as handle:
        for line in handle:
            if line.startswith('##INFO=<ID=CSQ,'):
                match = re.search(r'Format: ([^"\n]+)', line)
                if not match:
                    raise SharedVepError('CSQ header has no Format declaration')
                csq_columns = match.group(1).rstrip('>').split('|')
                missing = set(REQUIRED_CSQ) - set(csq_columns)
                if missing or len(csq_columns) != len(set(csq_columns)):
                    raise SharedVepError('Missing or duplicate required CSQ fields: ' + ','.join(sorted(missing)))
            elif line.startswith('#CHROM\t'):
                if csq_columns is None:
                    raise SharedVepError('Shared VEP artifact lacks CSQ schema')
                header = True
            elif line.startswith('#'):
                continue
            else:
                if not header:
                    raise SharedVepError('Shared VEP artifact has no VCF column header')
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) < 8 or ',' in fields[4] or len(fields[0]) > 2:
                    raise SharedVepError('CADD view requires explicit supported-contig, biallelic variant representation')
                try:
                    if int(fields[1]) < 1:
                        raise ValueError()
                except ValueError:
                    raise SharedVepError('Invalid VCF coordinate')
                info = dict(item.split('=', 1) for item in fields[7].split(';') if '=' in item)
                consequences = []
                for annotation in info.get('CSQ', '').split(','):
                    values = annotation.split('|')
                    if len(values) != len(csq_columns):
                        raise SharedVepError('CSQ row width differs from declared schema')
                    consequences.append(dict(zip(csq_columns, values)))
                # annotateVEPvcf.py transfers motif information to other rows,
                # then removes motif-only rows when any non-motif row exists.
                mixed_motifs = any(c['Feature_type'] != 'MotifFeature' for c in consequences)
                for consequence in consequences:
                    if mixed_motifs and consequence['Feature_type'] == 'MotifFeature':
                        continue
                    yield (fields[0], fields[1], fields[3], fields[4], consequence['Feature'] or 'NA')
    if not header:
        raise SharedVepError('Shared VEP artifact has no VCF header')


def verify(bundle, profile, stub=False):
    settings = profile.get('annotation', {}).get('shared_vep')
    if not settings or not settings.get('qualification'):
        raise SharedVepError('Shared VEP selection is unqualified: qualify the CADD --per_gene view, transcript release and FeatureID joins first')
    release = settings.get('transcript_release')
    if not release or release != settings.get('vep_release') or release != settings.get('esm_transcript_release'):
        raise SharedVepError('Shared VEP and ESM transcript releases differ')
    if profile['contract']['assembly'] == 'GRCh38' and release != settings.get('mmsplice_transcript_release'):
        raise SharedVepError('Shared VEP and MMSplice transcript releases differ')
    qualification = settings['qualification']
    if sha256(qualification['path']) != qualification['sha256']:
        raise SharedVepError('Shared VEP qualification checksum differs')
    evidence = json_file(qualification['path'])
    required = {'schema_version': 1, 'status': 'qualified', 'assembly': profile['contract']['assembly'],
                'vep_release': settings['vep_release'], 'transcript_release': release,
                'selection': 'cadd-per-gene', 'adapter_id': settings['adapter_id'], 'stub': bool(stub)}
    if any(evidence.get(key) != value for key, value in required.items()):
        raise SharedVepError('Shared VEP selection qualification does not match the installed profile')
    if evidence.get('checks') != {'per_gene_selection': 'exact', 'feature_id_joins': 'exact',
                                 'variant_representation': 'exact', 'missingness': 'preserved'}:
        raise SharedVepError('Shared VEP qualification lacks required exactness checks')
    proof = evidence.get('evidence', {})
    if not proof.get('path') or sha256(proof['path']) != proof.get('sha256'):
        raise SharedVepError('Shared VEP qualification evidence is unavailable or changed')
    receipt_path = Path(bundle) / 'manifest.json'
    receipt = json_file(receipt_path)
    for key, value in {'schema_version': 1, 'stage': 'cadd-shared-vep-view',
                       'assembly': profile['contract']['assembly'], 'profile_id': profile['profile_id'],
                       'resource_manifest_sha256': profile['contract']['resource_manifest_sha256'],
                       'vep_release': settings['vep_release'], 'transcript_release': release,
                       'selection': 'cadd-per-gene', 'adapter_id': settings['adapter_id'],
                       'qualification_sha256': qualification['sha256'], 'stub': bool(stub)}.items():
        if receipt.get(key) != value:
            raise SharedVepError('Shared VEP receipt/profile mismatch: ' + key)
    source = receipt['source_vep']
    if not Path(source['path']).is_absolute() or sha256(source['path']) != source['sha256']:
        raise SharedVepError('Original single-run PriVA VEP artifact changed')
    vcf = Path(bundle) / 'vep.vcf.gz'
    if sha256(vcf) != receipt.get('output_sha256'):
        raise SharedVepError('Shared VEP CADD view changed')
    count = sum(1 for _ in expected_features(vcf))
    if count != receipt.get('expected_annotation_rows'):
        raise SharedVepError('Shared VEP expected feature-row count differs')
    return vcf, {'input_manifest_sha256': sha256(receipt_path),
                 'shared_vep_source_sha256': source['sha256'],
                 'shared_vep_view_sha256': sha256(vcf),
                 'selection_qualification_sha256': qualification['sha256']}


def verify_feature_join(vcf, annotation_rows, columns):
    feature_index = columns.index('FeatureID')
    expected = expected_features(vcf)
    for row in annotation_rows:
        actual = (row[0], row[1], row[2], row[3], row[feature_index])
        if actual != next(expected, None):
            raise SharedVepError('CADD annotation changed variant/FeatureID join or source order')
    if next(expected, None) is not None:
        raise SharedVepError('CADD annotation lost shared VEP feature rows')
