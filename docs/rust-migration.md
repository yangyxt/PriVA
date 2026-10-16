# Rust migration and exact scientific parity

The migration has two independent acceptance gates: first reproduce the frozen
Python/shell workflow under Nextflow, then replace Python modules one at a time
with Rust implementations. Mature-tool shell wrappers remain in place. Python
remains the production scientific backend until the supplied family passes both
gates using the same resource bytes, tool versions, parameters and input files.
An orchestration test or a synthetic Rust test does not establish family parity.

The baseline is the source snapshot commit `09223ef` on
`feature/nextflow-rust-migration`. Original `main` remains unchanged in its own
worktree. At this baseline there are 56 Python files under `scripts/`, of which
40 import pandas, and four Python files under `data/`, of which two import
pandas. The 19 existing Python test files are test infrastructure, not production
modules to rewrite. New orchestration and independent acceptance utilities may
remain small Python programs; all production pandas table operations are in the
Rust/Polars migration scope, including resource preparation.

The baseline also imports NumPy in 30 scripts, pysam in 15, SciPy in seven,
scikit-learn in one (three import statements), statsmodels and pybedtools in one
each. These are separate compatibility requirements beyond dataframe syntax.
`cadd_hub_util.py` already uses Python Polars; its cache-table operations can move
to the Rust crate, but this does not replace the official CADD scoring engine.

## Implemented first boundary

`rust/` contains `priva-tabular final-order`, an opt-in terminal sorting module.
It uses Polars grouping/reduction and multicolumn sorting, then native Rust
first-appearance factorization. It consumes **already computed** `sort_index`
and `ACMG_quant_score` after the legacy protein-coding filter and all ranking
penalties; it does not apply that filter, compute either value or assign any criterion.
It preserves the original TSV cell strings and columns and only changes row order
and `variant_rank`. It is not wired into production shell wrappers.

The order reproduces the last block of `sort_and_rank_variants`:

1. Maximum `sort_index` for each `(chrom, pos, ref, alt)` group, descending.
2. Row `ACMG_quant_score`, descending.
3. `chrom`, integer `pos`, `ref`, `alt`, ascending.
4. Original row position for a complete tie.
5. One-based `variant_rank` from first appearance of `variant_id` in that order.

Do not substitute dense rank of the score or assume every transcript of a
variant will be contiguous. The legacy secondary score can interleave variants
whose group maxima tie. A row's `sort_index` is not the secondary sort key.

The first implementation accepts nonmissing keys and finite Float64 scores,
positive integer positions and a present `variant_id`. Unsupported missing keys,
NaN/infinite scores or the legacy no-`variant_id` fallback fail explicitly.
`--chrom-type string|integer` is mandatory: pandas' inferred chromosome dtype can
change `10` versus `2` ordering. Select the dtype of the actual legacy dataframe;
do not guess from the genome build. String comparisons preserve original contig
spelling. No assembly, contig, allele or transcript normalization is added here.

Toolchain: Rust **1.90.0**, `polars-core = 0.51.0`, `csv = 1.3.1`, exact direct
versions and transitive `Cargo.lock`. The limited Polars feature set avoids
building unused readers and query engines for this first boundary. Subsequent
modules can add pinned features after review. Build with `--locked`; deployment
should package the tested binary in the digest-pinned pipeline image, record its
SHA256 and use portable target features. A local debug build is development
validation, not a production release image.

```bash
export CARGO_TARGET_DIR=/path/to/build/rust
cargo +1.90.0 build --manifest-path rust/Cargo.toml --locked --jobs 1
POLARS_MAX_THREADS=1 "$CARGO_TARGET_DIR/debug/priva-tabular" final-order \
  preordered.tsv new-ordered.tsv --chrom-type string
python rust/tests/differential_final_order.py "$CARGO_TARGET_DIR/debug/priva-tabular"
```

The Rust tool refuses to overwrite an existing output. Temporary/build files and
validation logs belong in the associated data/output tree. The bounded development
build disables debug symbols. Set `POLARS_MAX_THREADS` from the Nextflow task CPU
allocation; do not let every table process use all node CPUs.

The explicit `+1.90.0` selects the pinned toolchain when running from the repository
root; `rust/rust-toolchain.toml` also selects it when running inside `rust/`.

The differential test extracts and executes the **actual terminal AST block** of
the frozen `scripts/acmg_scoring.py`, using pandas as the reference. It covers
empty/header-only tables, ties, repeat variant IDs, chromosome dtype, numeric
position order, quoted TSV strings, existing rank replacement and finite score
differences. It does not import or test the earlier SciPy/control-frequency,
inheritance, PEXT or gene-list penalty calculations.

## Exact acceptance comparator

`bin/priva-compare` is independent of pandas and Polars. It requires the final
ranked TSV **and** the pre-ranking `.acmg.tsv` matrix from both runs:

```bash
bin/priva-compare \
  --reference-final legacy/family.anno.filtered.tsv \
  --candidate-final nextflow/family.anno.filtered.tsv \
  --reference-acmg legacy/family.anno.filtered.acmg.tsv \
  --candidate-acmg nextflow/family.anno.filtered.acmg.tsv \
  --report parity.json
```

Exit 0 means all parsed cells, header order and row order match. Exit 1 is a
difference; exit 2 is an invalid/missing input. The report records SHA256 for
each file and separately reports byte equality; LF versus CRLF or equivalent CSV
quoting does not change parsed TSV values. No other normalization is applied:

- `sort_index`, `variant_rank`, `ACMG_quant_score`, `ACMG_class` and
  `ACMG_criteria` are required in the final table; all additional columns are
  checked, including expression, inheritance/mechanism and control penalties.
- Every one of the 28 matrix criteria and all audit columns must agree. Strength
  `1` is not equivalent to `2`; summary token ordering and suffix spelling are
  checked exactly. Blank, `NaN`, `NA`, `.`, `0`, `False` and whitespace remain
  distinct strings. Even `0.950` versus `0.95` fails; there is no numeric
  tolerance or hidden rounding. Resolve serialization deliberately rather than
  relax this acceptance gate.
- Final diagnostic row keys include variant ID, coordinates, transcript
  `Feature` when present and duplicate occurrence. The legacy ACMG matrix lacks
  transcript identity, so its key is variant ID plus occurrence in **original
  matrix order**. The comparator never deduplicates repeated IDs or joins the
  matrix to the final sorted table by variant ID alone.
- Empty results fail the single-family acceptance by default. `--allow-empty`
  is only for a known empty test fixture; zero-byte legacy output requires all
  four files to be zero-byte. Header-only and zero-byte forms are not silently
  equated. Duplicate headers, missing required columns and ragged rows fail.

The comparator establishes output agreement, not resource identity or clinical
correctness. The run manifest must separately prove the same input hashes,
assembly, pedigree, controls, resource bytes and tool/model versions. Do not
refresh cache files during either run. Use isolated writable task directories
around the existing immutable resources. VCF hub state must be disabled for both
runs or frozen and qualified for the identical cohort/control input.

## Source-faithful summary and posterior milestone

The same crate now exposes `priva score-criteria`, with function modules for
`fix_clingen_term`, `translate_strength_level`, `create_criteria_summary` and
`calculate_posterior_probability`. It consumes all 28 **already assigned** integer
criteria (0–5). ClinGen pickle loading and evidence-map resolution remain upstream;
the caller either explicitly disables override or supplies resolved per-row strings.
Only the source defaults (prior 0.1, exponential base 2, odds 350) are qualified.

```text
Validate typed input → summarize all rows → replace one criterion across all rows
                     → repeat criteria in source order → posterior → serialize
```

The column-outer override order preserves which malformed row fails first.
Integer BA1 zero versus a mixed floating column's `0.0` and the original empty
DataFrame failure are part of the contract. No criterion assignment, ranking
penalties or full family interpretation are claimed by this interface.

The installed reference environment (Python 3.12.3, NumPy 1.26.4, pandas 2.2.3)
passed **19,754 exact parsed-cell rows**, **16 failure/overwrite cases**, **9 Rust
tests**, and the existing **60 terminal-order cases + 5 rejections**. The scoring
oracle checks source SHA256
`9111a76e3fc94a064e14d97b535bb1dd2132f7c61836f8e8892477d4c5e69646`.
Use `rust/tests/differential_scoring.py` with the built `priva` executable.
These synthetic tests do not qualify complete ACMG assignment or family results.

## Port order and boundaries

| Stage | Python scope | Rust/Polars replacement and gate |
| --- | --- | --- |
| 1 | Terminal block of `acmg_scoring.py` | Implemented experimental group max, ordering and factorization. Synthetic differential tests precede the family test. |
| 2 | `generate_varID_for_vcf.py`, `collapse_HPO_anno.py`, `combine_annotations.py` | Stable variant/transcript IDs, typed table schema, HPO fields, CADD transcript lookup, duplicate and join behavior. Compare intermediate tables before downstream scoring. |
| 3 | `pedigree_filtration_per_fam.py`, `find_cosegregation_vars.py`, `determine_phase.py` | Explicit PED/genotype/phase model, family/control membership and ordered transmission results. Keep mature VCF tools and wrappers. |
| 4 | `acmg_consequence.py`, `acmg_inheritance.py`, splicing/UTR helpers and mechanism modules | Typed evidence and condition/mechanism assertions; preserve all source provenance, exact-match precedence and unknown values. |
| 5 | `acmg_pvs1_*`, `acmg_ps*`, `acmg_pm*`, `acmg_pp*`, `acmg_bs*`, `acmg_variant_mechanism.py`, `acmg_criteria_assign.py` | One evidence family at a time, strengths and audit fields before classification. Freeze ClinGen override behavior. |
| 6 | Remaining `acmg_scoring.py`, `resort_step3_only.py` | Bayesian calculation, rounding, classification, SciPy binomial/Beta behavior and ranking penalties. Require exact published numerical fields and all classification boundaries. |
| 7 | `build_*`, `prepare_*`, `stat_*`, AlphaMissense/ClinVar/MaveDB conversion helpers | Resource builders and pandas table preparation; preserve source data identities and validate converted caches against the same bytes before creating a new immutable resource release. |

The inventory scope is all production scripts, not only files named `acmg_*`.
Pure Python graph/ontology/sequence/statistical logic also needs an explicit Rust
module, a maintained mature tool adapter, or a documented temporary compatibility
boundary. Small tests and orchestration validators are not scientific backends.

## Semantics that must be specified before a port is enabled

- **Nulls and dtypes:** pandas `read_table` inference and its default missing
  vocabulary differ from explicit Polars nulls. Preserve integer PED IDs versus
  strings, nullable booleans, float32 narrowing during annotation conversion,
  float64 scoring, string `nan` created by `astype(str)` and `errors='coerce'`.
  Do not equate missing evidence with strength zero without the legacy rule.
- **Join cardinality/order:** preserve one-to-many transcript rows, row alignment,
  first-match lookups, group order, `drop_duplicates(keep='first')` and any
  index-aligned assignment. `combine_annotations.py` keeps the first CADD row
  per `(variant_id, Feature)`. Rust hash-map iteration is not output order.
- **Grouping and ties:** pandas groupby drops null group keys by default;
  Polars groups nulls. `factorize(sort=False)` uses first appearance and maps
  missing IDs to code -1 (rank zero after +1). The initial Rust scope rejects
  those cases; a future broader module must reproduce them or document an
  explicitly reviewed behavior change. Multikey ties retain source order.
- **Numbers:** keep operation order as well as Float64 precision. NumPy power,
  Python four-place rounding and Rust `round` are not interchangeable promises.
  SciPy `binom.cdf` and `beta.ppf` need differential boundary fixtures; a generic
  Rust statistics replacement is not sufficient evidence of exact agreement.
  Do not recompute or round `sort_index` merely to make comparisons pass.
- **Criteria:** preserve criterion order and strength suffixes, ClinGen's
  replacement override (including missing curated evidence), frequency conflict
  gates, inheritance/mechanism compatibility, PP1/BS4/BP2/PM3 relationships and
  exact-class thresholds. Compare the complete matrix, not just final labels.
- **Pickle resources:** Rust cannot safely assume arbitrary Python pickle is a
  portable cache format. Initially use the unchanged Python reader boundary.
  Later export a typed Arrow/Parquet/JSON representation once, record source
  pickle SHA256 and converter identity and prove key/value/ordering parity.
  Keep the original cache bytes; conversion is a new derived resource, not an
  in-place cache replacement. Do not rebuild a cache from a newer public source.

Enable a Rust boundary only after intermediate differential fixtures, the exact
single-family acceptance and resource/tool provenance all pass. Keep a per-module
backend switch for investigation and rollback until the full scientific path is
validated. Full Rust migration, complete cache-format conversion, offline release
packaging and the user-supplied family acceptance remain separate deliverables.

API references used for the first module: [Polars 0.51 DataFrame](https://docs.rs/polars-core/0.51.0/polars_core/frame/struct.DataFrame.html)
and [Polars 0.51 GroupBy](https://docs.rs/polars-core/0.51.0/polars_core/frame/group_by/struct.GroupBy.html).
