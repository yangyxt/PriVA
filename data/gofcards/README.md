# GoFCards exact GOF lookup

`gofcards_exact_gof_hgvsp.tsv.gz` is a compact PriVA runtime lookup derived
from the standalone `gofcards_hg38_normalizer` workflow.

PriVA uses this file only for exact variant-level GOF matching. Two evidence
routes are supported by `scripts/gene_mechanism_hub.py`:

- `match_gofcards_variant_gof(...)`: normalized HGNC symbol plus exact protein
  change/HGVSp.
- `match_gofcards_variant_gof_by_genomic(...)`: normalized HGNC symbol plus
  exact hg19/hg38 genomic allele key. `key_type=vcf` uses VCF-padded indel
  alleles; `key_type=genomic` uses the sparse GoFCards source fields;
  `key_type=auto` tries both.

Neither route may be interpreted as evidence that all variants in a listed gene
are gain-of-function.

Current cache:

- refreshed: 2026-07-07
- source workbook: `/paedyl01/disk1/yangyxt/gofcards_hg38_normalizer/work_20260707_full/gofcards_hg38_normalized_workbook_20260707.xlsx`
- exported TSV: `/paedyl01/disk1/yangyxt/gofcards_hg38_normalizer/work_20260707_full/gofcards_priva_exact_gof_hgvsp_20260707.tsv.gz`
- rows: 14,587 transcript/genomic records
- unique exact HGVSp-match keys: 4,005 normalized `symbol+hgvsp_key` pairs
- rows without HGVSp retained for future genomic matching: 6,813

The table includes hg19 and hg38 genomic coordinates/ref/alt, VEP-derived
HGVSc/HGVSp when available, GoFCards audit metadata, HGNC-resolved symbols, and
GoFCards-to-VEP HGVSc/HGVSp concordance flags. It also includes VCF-padded
`hg19_vcf_*` and `hg38_vcf_*` keys for deletion/insertion records whose source
REF or ALT was blank. The stable PriVA input path remains
`data/gofcards/gofcards_exact_gof_hgvsp.tsv.gz`.

Validated examples:

- `CFTR` + `p.Phe508del` matches by HGVSp.
- `CFTR` + hg38 VCF key `7|117559592|CTTT|C` matches by genomic allele.
- `TP53` + hg38 key `17|7669058|A|G` matches by genomic allele even though
  the GoFCards row has no usable HGVSp.
