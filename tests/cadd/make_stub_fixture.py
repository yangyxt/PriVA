#!/usr/bin/env python3
"""Synthetic full-schema wiring fixtures; never use for scientific acceptance."""
import argparse
import hashlib
import sys
import gzip
import json
from pathlib import Path

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / 'bin'))
from cadd_shared_vep import REQUIRED_CSQ

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', required=True)
args = parser.parse_args()
out = Path(args.output)
out.mkdir(parents=True, exist_ok=False)
proof = out / 'synthetic-selection-proof.json'
proof.write_text(json.dumps({'scope': 'synthetic-wiring-only', 'cases': 1}))
original = out / 'single-priva-vep.vcf'
values = {key: '' for key in REQUIRED_CSQ}
values.update(Allele='G', Consequence='missense_variant', Feature_type='Transcript', Feature='STUB_FEATURE')
original.write_text('##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: ' + '|'.join(REQUIRED_CSQ) + '">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tG\t.\t.\tCSQ=' + '|'.join(values[key] for key in REQUIRED_CSQ) + '\n')
for assembly in ['GRCh37', 'GRCh38']:
    profile = json.loads((root / 'conf/cadd' / (assembly + '.reference.json')).read_text())
    profile['stub'] = True
    profile['software_id'] = 'synthetic-wiring-only'
    profile['contract']['resource_manifest_sha256'] = '0' * 64
    qualification = out / (assembly + '.qualification.json')
    qualification.write_text(json.dumps({'schema_version': 1, 'status': 'qualified', 'assembly': assembly,
        'vep_release': '110', 'transcript_release': '110', 'selection': 'cadd-per-gene', 'adapter_id': 'stub-only',
        'stub': True, 'checks': {'per_gene_selection': 'exact', 'feature_id_joins': 'exact', 'variant_representation': 'exact', 'missingness': 'preserved'},
        'evidence': {'path': str(proof), 'sha256': sha(proof)}}))
    profile['annotation']['shared_vep'].update(adapter_id='stub-only', qualification={'path': str(qualification), 'sha256': sha(qualification)})
    shared = out / ('shared-' + assembly)
    shared.mkdir()
    with gzip.open(shared / 'vep.vcf.gz', 'wb') as handle:
        handle.write(original.read_bytes())
    (shared / 'manifest.json').write_text(json.dumps({'schema_version': 1, 'stage': 'cadd-shared-vep-view',
        'assembly': assembly, 'profile_id': profile['profile_id'], 'resource_manifest_sha256': '0' * 64,
        'vep_release': '110', 'transcript_release': '110', 'selection': 'cadd-per-gene', 'adapter_id': 'stub-only',
        'qualification_sha256': sha(qualification), 'stub': True,
        'source_vep': {'path': str(original), 'sha256': sha(original)},
        'output_sha256': sha(shared / 'vep.vcf.gz'), 'expected_annotation_rows': 1}))
    (out / (assembly + '.json')).write_text(json.dumps(profile, indent=2) + '\n')
with gzip.open(out / 'input.vcf.gz', 'wt') as handle:
    handle.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tG\t.\t.\t.\n')
