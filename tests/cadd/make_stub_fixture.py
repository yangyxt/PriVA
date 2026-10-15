#!/usr/bin/env python3
"""Synthetic full-schema wiring fixtures; never use for scientific acceptance."""
import argparse
import gzip
import json
from pathlib import Path

root = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', required=True)
args = parser.parse_args()
out = Path(args.output)
out.mkdir(parents=True, exist_ok=False)
for assembly in ['GRCh37', 'GRCh38']:
    profile = json.loads((root / 'conf/cadd' / (assembly + '.reference.json')).read_text())
    profile['stub'] = True
    profile['software_id'] = 'synthetic-wiring-only'
    profile['contract']['resource_manifest_sha256'] = '0' * 64
    (out / (assembly + '.json')).write_text(json.dumps(profile, indent=2) + '\n')
with gzip.open(out / 'input.vcf.gz', 'wt') as handle:
    handle.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tG\t.\t.\t.\n')
