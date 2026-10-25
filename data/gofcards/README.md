# GoFCards resources

`gofcards_exact_gof.json.gz` contains bundled upstream exact-variant records.
`gofcards_mechanism_reviews.tsv` supplies curated mechanism-polarity decisions.

| Installer function in `scripts/install_utils.sh` | Builder/consumer | Relationship |
| --- | --- | --- |
| `gofcards_exact_gof_cache_install` | `scripts/build_gofcards_exact_gof_cache.py` | Acquire and normalize the upstream exact-variant cache |
| `gene_nonlof_mechanism_cache_install` | `scripts/build_gene_nonlof_mechanism_cache.py` | Include eligible records in canonical mechanism assertions |
| Runtime mechanism annotation | `scripts/gene_mechanism_variants.py` | Read canonical assertions for exact variant matching |

Runtime authority is the canonical mechanism cache under `data/patho_mechanism`
or its configured managed replacement. A gene's inclusion does not classify
all its variants as gain-of-function. Transcript identities, assembly-specific
alleles and polarity/quarantine decisions remain part of matching.
See [cache schema](../../docs/GOFCARDS_CACHE_SCHEMA.md),
[normalization methods](../../docs/GOFCARDS_NORMALIZATION_METHODS.md) and the
[user manual](../../README.md).
