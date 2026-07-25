# GoFCards exact GOF upstream cache

`gofcards_exact_gof_hgvsp.tsv.gz` is a compact build-time product from the
standalone `gofcards_hg38_normalizer` workflow. It is not an independent PriVA
runtime authority.

The canonical builder validates this table and embeds runtime-eligible alleles
in `gene_nonlof_mechanism_curated_assertions.json`. PriVA's exact matcher reads
that canonical JSON. Its two exact-match routes are:

- `match_gofcards_variant_gof(...)`: normalized HGNC symbol plus exact protein
  change/HGVSp.
- `match_gofcards_variant_gof_by_genomic(...)`: normalized HGNC symbol plus
  exact hg19/hg38 genomic allele key. `key_type=vcf` uses VCF-padded indel
  alleles; `key_type=genomic` uses the sparse GoFCards source fields;
  `key_type=auto` tries both.

Neither route may be interpreted as evidence that all variants in a listed gene
are gain-of-function.

Current compact cache:

- refreshed: 2026-07-25
- workflow: `/paedyl01/disk1/yangyxt/gofcards_hg38_normalizer/bin/gofcards_workflow.sh`
- SHA256: `039ef4b151d1afbf44d13f03619c87c786d35e544638e34c8a5ea8f6567b7f5a`
- rows: 4,009 transcript-level records across 2,033 unique source alleles
- eligible rows: 3,967, including 3,906 with HGVSp and 61 genomic-only rows
- quarantined rows: 42 transcript records across 32 genuinely discordant source alleles

The table preserves the curated `GoFCards_HGNC_Symbol` as `HGNC_Symbol` and
stores `VEP_HGNC_Symbol` separately. Current HGNC aliases are resolved before
the two genes are compared. A genuine disagreement is retained with
`match_eligibility=quarantined_gene_discordance` for audit, but it cannot enter
canonical exact matching, ClinVar linking, or the HPO condition cache.

`*_genomic_key` preserves the source allele representation; `*_vcf_key` stores
the normalized VCF representation used for exact caller-style matching. The
stable upstream input path remains
`data/gofcards/gofcards_exact_gof_hgvsp.tsv.gz`.
