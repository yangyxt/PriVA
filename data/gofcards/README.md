# GoFCards exact GOF lookup

`gofcards_exact_gof_hgvsp.tsv.gz` is a compact PriVA runtime lookup derived
from the standalone `gofcards_hg38_normalizer` workflow.

PriVA uses this file only for exact variant-level GOF matching by normalized
HGNC symbol plus protein change. It must not be interpreted as evidence that
all variants in a listed gene are gain-of-function.

Current cache:

- refreshed: 2026-07-07
- source workbook: `/paedyl01/disk1/yangyxt/gofcards_hg38_normalizer/work_20260707_full/gofcards_hg38_normalized_workbook_20260707.xlsx`
- exported TSV: `/paedyl01/disk1/yangyxt/gofcards_hg38_normalizer/work_20260707_full/gofcards_priva_exact_gof_hgvsp_20260707.tsv.gz`
- rows: 7,687 transcript-level records
- unique exact-match keys: 3,935 normalized `symbol+hgvsp_key` pairs

The table includes hg19 and hg38 genomic coordinates/ref/alt, VEP-derived
HGVSc/HGVSp when available, and GoFCards audit metadata. The stable PriVA input
path remains `data/gofcards/gofcards_exact_gof_hgvsp.tsv.gz`.
