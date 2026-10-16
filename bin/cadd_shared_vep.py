#!/usr/bin/env python3
"""Fail-closed contract for a CADD-compatible view of PriVA's single VEP run.

This validates a previously qualified selection adapter. It does not invent a
replacement for VEP --per_gene or silently select/drop PriVA transcript rows.
"""
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import re

REQUIRED_CSQ = ['Allele', 'Consequence', 'Feature_type', 'Gene', 'Feature', 'SYMBOL',
                'CCDS', 'INTRON', 'EXON', 'cDNA_position', 'CDS_position', 'Protein_position',
                'Amino_acids', 'DOMAINS', 'SIFT', 'PolyPhen', 'MOTIF_NAME', 'HIGH_INF_POS',
                'MOTIF_SCORE_CHANGE']
NATIVE_ADAPTER = 'official-vep-native-per-gene-v1'
SHARED_TRANSCRIPT_POLICY = {
    'cache_mode': 'merged', 'priva_transcript_sources': ['Ensembl', 'RefSeq'],
    'use_transcript_ref': True, 'selection': 'official-per-gene-over-explicit-source-pool',
    'annotation_interpretation': 'CADD-derived',
}
NATIVE_REQUIRED_OPTIONS = {'merged': '1', 'use_transcript_ref': '1', 'regulatory': '1',
                           'sift': 'b', 'polyphen': 'b', 'ccds': '1', 'domains': '1',
                           'numbers': '1', 'canonical': '1', 'total_length': '1'}


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


def variant_keys(path):
    with text_file(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 8:
                raise SharedVepError('Source/view has an incomplete VCF record')
            yield (fields[0], fields[1], fields[3], fields[4])


def native_context(profile):
    settings = profile['annotation']['shared_vep']
    policy = settings.get('transcript_policy', {})
    if (profile.get('profile_kind') != 'adapted'
            or any(policy.get(key) != value for key, value in SHARED_TRANSCRIPT_POLICY.items())
            or policy.get('cadd_transcript_sources') not in [['Ensembl'], ['Ensembl', 'RefSeq']]):
        raise SharedVepError('Shared merged/transcript-reference annotation requires the explicit adapted CADD-derived policy')
    context = settings.get('selection_context', {})
    options = context.get('annotation_options', {})
    if (context.get('adapter_id') != NATIVE_ADAPTER or context.get('vep_release') != settings['vep_release']
            or context.get('perl_perturb_keys') != '0' or not context.get('perl_hash_seed')
            or context.get('cadd_transcript_sources') != policy['cadd_transcript_sources']
            or not context.get('pick_order') or not context.get('perl_version')
            or not re.fullmatch('[0-9a-f]{64}', context.get('plugin_sha256', ''))
            or not re.fullmatch('[0-9a-f]{64}', context.get('output_factory_sha256', ''))):
        raise SharedVepError('Native selection Perl/runtime context is incomplete or incompatible')
    if any(options.get(key) != value for key, value in NATIVE_REQUIRED_OPTIONS.items()) or options.get('assembly') != profile['contract']['assembly']:
        raise SharedVepError('Native selection context conflicts with the fixed PriVA merged/transcript-reference profile')
    return hashlib.sha256(json.dumps(context, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


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
    context_sha = None
    if settings['adapter_id'] == NATIVE_ADAPTER:
        context_sha = native_context(profile)
        if evidence.get('selection_context_sha256') != context_sha or evidence.get('transcript_policy') != settings['transcript_policy']:
            raise SharedVepError('Native selection context is not independently qualified')
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
    if context_sha and receipt.get('selection_context_sha256') != context_sha:
        raise SharedVepError('Native selection receipt context differs')
    if context_sha and receipt.get('transcript_policy') != settings['transcript_policy']:
        raise SharedVepError('Native selection receipt transcript policy differs')
    for original, selected in itertools.zip_longest(variant_keys(source['path']), variant_keys(vcf)):
        if original != selected:
            raise SharedVepError('Shared VEP view changed variant coverage, representation or record order')
    count = sum(1 for _ in expected_features(vcf))
    if count != receipt.get('expected_annotation_rows'):
        raise SharedVepError('Shared VEP expected feature-row count differs')
    return vcf, {'input_manifest_sha256': sha256(receipt_path),
                 'shared_vep_source_sha256': source['sha256'],
                 'shared_vep_view_sha256': sha256(vcf),
                 'selection_qualification_sha256': qualification['sha256'],
                 'selection_context_sha256': context_sha}


def verify_feature_join(vcf, annotation_rows, columns):
    feature_index = columns.index('FeatureID')
    expected = expected_features(vcf)
    for row in annotation_rows:
        actual = (row[0], row[1], row[2], row[3], row[feature_index])
        if actual != next(expected, None):
            raise SharedVepError('CADD annotation changed variant/FeatureID join or source order')
    if next(expected, None) is not None:
        raise SharedVepError('CADD annotation lost shared VEP feature rows')
