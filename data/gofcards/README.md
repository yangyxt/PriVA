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

- refreshed: 2026-07-26
- workflow: `/paedyl01/disk1/yangyxt/gofcards_hg38_normalizer/bin/gofcards_workflow.sh`
- SHA256: `9b9f2d914192753b517b8a3fbbd3821c68299b1c0d10c6e332eb1606869e0843`
- rows: 4,009 transcript-level records across 2,033 unique source alleles
- eligible rows: 3,960, including 3,900 with HGVSp and 60 genomic-only rows
- quarantined rows: 49 transcript records across 32 genuinely discordant source alleles
- allele-level propagation: 7 otherwise eligible sibling rows quarantined

The table preserves the curated `GoFCards_HGNC_Symbol` as `HGNC_Symbol` and
stores `VEP_HGNC_Symbol` separately. Current HGNC aliases are resolved before
the two genes are compared. A genuine disagreement is retained with
`match_eligibility=quarantined_gene_discordance` for audit, but it cannot enter
canonical exact matching, ClinVar linking, or the HPO condition cache. Every
other row for the same source allele is marked
`quarantined_allele_gene_discordance`, so a concordant sibling cannot bypass
the allele-level quarantine.

`*_genomic_key` preserves the source allele representation; `*_vcf_key` stores
the normalized VCF representation used for exact caller-style matching. The
stable upstream input path remains
`data/gofcards/gofcards_exact_gof_hgvsp.tsv.gz`.
