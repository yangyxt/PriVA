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
- rows: 3,923 VEP-calibrated runtime records
- columns: 39
- unique exact HGVSp-match keys: 2,097 normalized `HGNC_Symbol+hgvsp_key` pairs
- rows without HGVSp retained for genomic matching: 51
- rows with missing hg19/hg38 padded allele keys: 0

The table is VEP-centric and compact: one normalized `HGNC_Symbol`, VEP-derived
`HGVSc`/`HGVSp` for concordant transcript matches, raw GoFCards HGVS for trace,
and nonblank hg19/hg38 padded genomic allele keys. Coordinate-only rows are
kept only when GoFCards lacks parseable HGVS but exact genomic matching remains
possible. The stable PriVA input path remains
`data/gofcards/gofcards_exact_gof_hgvsp.tsv.gz`.

Validated examples:

- `CFTR` + `p.Phe508del` matches by HGVSp.
- `CFTR` + hg38 VCF key `7|117559592|CTTT|C` matches by genomic allele.
- `TP53` + hg38 key `17|7669058|A|G` matches by genomic allele even though
  the GoFCards row has no usable HGVSp.
