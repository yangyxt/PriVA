# GoFCards exact GOF lookup

`gofcards_exact_gof_hgvsp.tsv.gz` is a compact PriVA runtime lookup derived
from the local GoFCards GOF extraction workflow on 2026-07-03.

PriVA uses this file only for exact variant-level GOF matching by normalized
HGNC symbol plus protein change. It must not be interpreted as evidence that
all variants in a listed gene are gain-of-function.

The table is intentionally smaller than the raw GoFCards download and contains
only fields needed for reproducible matching and audit output.
