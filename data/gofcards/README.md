# GoFCards exact GOF upstream cache

`gofcards_exact_gof_hgvsp.tsv.gz` is a build-time product. It is not an
independent PriVA runtime authority.

## How it is built

- Builder: `/paedyl01/disk1/yangyxt/PriVA/scripts/build_gofcards_exact_gof_cache.py`
- Deployment: `bash scripts/install_utils.sh gofcards_exact_gof_cache_install config.yaml`
- Methods and limitations: `/paedyl01/disk1/yangyxt/PriVA/docs/GOFCARDS_NORMALIZATION_METHODS.md`
- Polarity reviews: `/paedyl01/disk1/yangyxt/PriVA/data/gofcards/gofcards_mechanism_reviews.tsv`

Seven steps: fetch the public workbook and both backend tables, build a GRCh37
VCF and check every reference allele against the genome, normalize, annotate with
VEP on GRCh37, lift over to GRCh38 with CrossMap and normalize again, annotate
with VEP on GRCh38, then join and export. GoFCards publishes GRCh37 coordinates
only, so GRCh37 is the sole ground truth and everything else is derived from it.

## Match routes

The canonical builder validates this table and embeds runtime-eligible alleles in
`gene_nonlof_mechanism_curated_assertions.json`. PriVA's exact matcher reads that
canonical JSON. In order of preference:

1. `match_gofcards_variant_gof(...)`: normalized HGNC symbol plus exact protein
   change (HGVSp). This is the default route.
2. The same route on HGVSc, which distinguishes different coding changes that
   produce the same protein change.
3. `match_gofcards_variant_gof_by_genomic(...)`: normalized HGNC symbol plus an
   exact genomic allele key for the assembly the query variant is mapped to.
   `key_type=vcf` uses the normalized VCF representation; `key_type=genomic` uses
   the sparse GoFCards source fields; `key_type=auto` tries both. Coordinates are
   an alternative route, not the default.

Neither route may be interpreted as evidence that all variants in a listed gene
are gain-of-function.

## Why every transcript is stored

Across this dataset, 66% of variants have transcripts that disagree with each
other about the protein change. Storing one transcript would make a match depend
on PriVA having picked the same one. All transcripts from both VEP runs are kept,
each with its version, so a match succeeds whichever transcript the query used.

The two VEP runs exist because the GRCh37 and GRCh38 caches ship different
transcript catalogues (GENCODE 19 and GENCODE 47). Both runs use the same VEP
cache directory that `config.yaml` points at, so transcript versions align with
PriVA's own annotation by construction.

## Eligibility and quarantine

A row is eligible only if its reference allele matched GRCh37, its resolved
curated gene equals its resolved VEP gene, and it carries no reviewed polarity
decision marking it non-GOF.

The table preserves the curated `GoFCards_HGNC_Symbol` as `HGNC_Symbol` and
stores `VEP_HGNC_Symbol` separately; current HGNC aliases and withdrawn symbols
are resolved on both sides before they are compared. A genuine gene disagreement
is retained with `match_eligibility=quarantined_gene_discordance` for audit but
cannot enter canonical exact matching, ClinVar linking, or the HPO condition
cache. Every other row for the same source allele is marked
`quarantined_allele_gene_discordance`, so a concordant sibling cannot bypass the
allele-level quarantine. Reviewed non-GOF alleles are marked
`quarantined_mechanism_polarity`, with the reason in `reject_reason`.

Reference-check failures never reach this table. They are written to
`gofcards_reference_rejects.tsv` in the build working directory, because an
allele whose reference does not match the genome describes a variant that does
not exist. Liftover failures do reach the table: they keep their HGVS and lose
only their GRCh38 coordinate key, flagged by `liftover_status=unmapped`.

## Evidence

Source records are aggregated per allele rather than reduced to one
representative, so an allele carrying many independent literature records keeps
all of their PMIDs, P scores, and function text. `evidence_record_count` records
how many source rows were combined.
