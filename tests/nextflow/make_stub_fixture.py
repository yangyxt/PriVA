#!/usr/bin/env python3
"""Create tiny synthetic Nextflow wiring fixtures; these are not biological truth."""
import argparse
import gzip
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()
root = Path(args.outdir)
root.mkdir(parents=True, exist_ok=True)
if any(root.iterdir()):
    raise SystemExit("Fixture output directory must be empty")
with gzip.open(root / "cohort.vcf.gz", "wt") as handle:
    handle.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tchild\tdad\tmom\tother\n1\t10\t.\tA\tG\t100\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0\n")
(root / "cohort.vcf.gz.tbi").write_text("STUB INDEX: NOT A VALID TABIX INDEX\n")
(root / "family.ped").write_text("OTHER\tother\t0\t0\t2\t2\nFAM\tchild\tdad\tmom\t1\t2\nFAM\tdad\t0\t0\t1\t1\nFAM\tmom\t0\t0\t2\t1\n")
(root / "legacy.yaml").write_text("assembly: hg19\naf_cutoff: 0.05\napply_clingen_override: true\n")
(root / "resources.json").write_text('{"bundle_id":"stub-only","stub":true}\n')
(root / "cohort.anno.cadd.tsv").write_text("#Chrom\tPos\tRef\tAlt\tPHRED\tFeatureID\n1\t10\tA\tG\t10\tSTUB\n")
