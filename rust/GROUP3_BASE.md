# Group 3 base: what the seven ACMG threads build on

Group 3 converts the ACMG decision step (`scripts/acmg_criteria_assign.py` and the scripts it
imports) to Rust. This note records the shared base laid down before the seven threads start
(2026-10-30): the table layer, the module skeleton, the comparison harness, the probes of the real
resource files and the decisions they led to. The group 3 plan
(`PriVA_upgrade/PriVA_Rust_Group3_Plan_20261030.md`, sections 5 and 6) is the longer background;
where this note differs, this note holds (the owner's decisions of 2026-10-30 came later).

## 1. Decisions applied

1. PP1: `--pp1_vcf` and `--pp1_ped` are accepted and forwarded to PP1_criteria (fixed on main in
   c529805); `types::StepArguments` carries both.
2. `determine_phase.py`: per-gene results merge in the genes' groupby order (as main's Python now
   does with `pool.imap`).
3. `beta.ppf` may differ from SciPy in the last digit; the other SciPy results must give the same
   decisions (`acmg_step::stats`).
4. Ordinary values are reproduced exactly; odd values are refused loudly (exit 1, a message naming
   the value or column, outputs untouched), as in groups 1 and 2.
5. polars is used for a pandas call whenever it gives the identical result on the real column
   types; own code remains for reading and writing tables as pandas does and for pandas `object`
   columns that mix kinds of values (section 3).

## 2. How Rust holds the table (decision D1, revised)

The table is `pandas_table::Table`: column names in order, one `Column` per name, each keeping
its pandas dtype. The reader produces the first five variants; the step's own pandas code makes
the last three, so they exist for the threads:

| `Column` variant | pandas dtype | polars Series (`pandas_polars::to_series`) |
| --- | --- | --- |
| `Int(Vec<i64>)` | `int64` | `Int64` |
| `Float(Vec<f64>)` (NaN is missing) | `float64` | `Float64`, NaN as null |
| `Bool(Vec<bool>)` | `bool` | `Boolean` |
| `BoolOrMissing(Vec<Option<bool>>)` | `object` of True/False/NaN | `Boolean` with nulls |
| `Text(Vec<Option<String>>)` | `object` of str/NaN | `String` with nulls |
| `Object(Vec<ObjectValue>)` | `object` mixing text, numbers, booleans, NaN, None | refused: own code |
| `NullableInt(Vec<Option<i64>>)` | `Int64` | `Int64` with nulls |
| `NullableBool(Vec<Option<bool>>)` | `boolean` | `Boolean` with nulls (Kleene logic matches) |

`pandas_polars::from_series(series, PandasDtype)` brings a polars result back as the pandas type
the Python code would have (a null becomes NaN, `None` or `pd.NA` as that type holds it, and a
null where the type cannot hold one is an error). The conversions copy each column once; numbers
are moved as plain vectors. Why not a polars `DataFrame` for the whole table: polars has no
`object` dtype and keeps null apart from NaN, and the step's pandas code depends on both (the
`fillna("")` of a float column, `~` on an object column, `x in [0, np.nan]`).

`Table` helpers: `column(name)` (the `KeyError` pandas raises), `has_column`, `set_column`
(`df[name] = ...`, keeps an existing column's place), `take_rows`, `filter_rows`,
`Column::cell`, `missing_mask` (`isna`), `dtype_name`, `take`.

Reading (`pandas_table::read_table_with(path, &ReadOptions)`) covers every `read_table` call of
the step: the default inference, `dtype=str`, `keep_default_na=False` (padding of short rows
stays empty text), `header=None` with `names`, `usecols` by position (the headerless PED, longer
rows allowed) and by name (file order kept; a missing name gives pandas' `ValueError` text), and
`comment="#"`. A one-column table with a second field (pandas would make it the row index) is
refused. Writing (`write_table`) now also writes the three later variants as `to_csv` does and
quotes a lone missing value in a one-column table (`""`), as `csv.writer` does.

## 3. pandas calls of group 3: polars or own code

Counts are the occurrences of the call text in the 26 scripts (`grep -o`, 2026-10-30); the plan's
figures counted slightly differently. "polars" names the operation in `pandas_polars` (or the
polars method) that gives the identical result; each listed equivalence was checked against
pandas 2.2.3 and polars 0.51.0 on the probe values (`pandas_polars` unit tests carry them).

| pandas call (count) | Argument shapes in group 3 | Decision |
| --- | --- | --- |
| `fillna` (149) | `False` 46, `0` 46, `""` 39, `"."` 9, others 9 | polars `fill_null` when the column is numeric (`0`, numbers), boolean (`False`/`True` on `Boolean`, `BoolOrMissing`, `NullableBool` or a `str_contains` result; pandas 2.2 turns an object True/False/NaN column into `bool` here) or text (`""`, `"."` on `Text`). **Own code** when the fill changes the kind of values: `fillna("")` or `fillna(".")` on a float column gives an object column of floats and text, `fillna(0)` on a text column gives text and integers (`Column::Object`), because polars has no object dtype. |
| `astype` (125) | `str` 45, `bool` 39, `int` 33, `float` 2, `Int64` 4, `boolean` 2 | `astype(str)`: **own code**, Python `str()` of each value (`'1e+16'`, `'nan'`, `'None'`, `'True'`); polars prints numbers differently. `astype(bool)`: no-op on `bool`; `float_astype_bool` (polars, NaN counts true) on floats; **own code** (Python truthiness: empty text false, NaN true) on object columns. `astype(int)`: polars cast from `bool`; `float_astype_int` (polars, truncation; NaN or inf raise pandas' `IntCastingNaNError`, values past 64 bits refused); **own code** from object. `astype(float)`: polars cast from integers; **own code** from text (Python `float()`). `Int64`/`boolean`: `NullableInt`/`NullableBool`, polars nulls. |
| `.str.*` (115) | `contains` 82, `lower` 9, `split` 7, `strip` 5, `upper` 3, `len` 3, `replace` 2, `count` 2, `removeprefix` 1, `match` 1 | `contains`: `str_contains` (polars regex; missing stays missing as pandas' `na=None`; refuses a `$` pattern on values with a line break, where Python's `re` also matches before the final `\n`). `lower`/`upper`: polars `to_lowercase`/`to_uppercase` (same on the step's text; the caches hold no non-ASCII character). `len`: polars `str_len_chars`, returned as pandas does (`float64` when a value is missing). **Own code**: `strip` (Python's whitespace includes U+001C to U+001F; polars' does not: `python_text::python_strip`), `split` (lists, then indexing), `replace` and `count` (pandas' regex defaults per call), `removeprefix`; `match` is polars `contains` with the pattern anchored by `^(?:...)`. |
| `.loc[...]` (94) | row filters and `df.loc[mask, col] = value` | Row filter: `Table::filter_rows` / `take_rows`. Masked assignment: polars `zip_with` (`where_mask`) when the value has the column's type and the mask has no missing values; **own code** when pandas changes the column's dtype on assignment (an int column given a float, a new column assigned on part of the rows, which pandas fills with NaN). |
| `.map` (47) | dicts (gene to score, symbol to code) and functions | **Own code**: a hash map lookup per row; the result type follows pandas' inference from the values (integer values with a missing key give `float64` with NaN; text values give `object`). |
| `isin` (27) | lists or sets of text, a few numbers | `isin_text`, `isin_numbers` (polars `is_in`, missing values answered false as pandas does; NaN in a number list matches the missing values). A text list holding `None` or NaN is refused (pandas matches NaN only by NaN, not by `None`). |
| `np.where` (22) | masks with numbers or text | `where_mask` (polars `zip_with`) for masks without missing values; refused otherwise (NaN is true in numpy, polars takes the second value). The result type is numpy's: mixing `np.nan` with text gives the text `'nan'`; mind it per call. |
| `groupby` (13) | `sort=False` on genes and variant keys, default `sort=True`, `dropna=False` once | `first_appearance_groups(table, keys, dropna)` (polars `group_by_stable`; pandas drops groups with a missing key unless `dropna=False`). `sort=True`: the same groups sorted by key (polars sort; Python and Rust order text by code point). Aggregations: polars where they skip missing values as pandas does (`max`, `sum`); `apply` per group is own code. |
| `drop_duplicates` (5) | whole table on read, `subset=[...]` | `drop_duplicates_first(table, subset)` (polars `unique_stable`, keep first; `-0.0` equals `0.0`, missing equals missing; object columns are keyed by Python equality, `1 == 1.0 == True`). |
| `sort_values` (1) | several keys, mixed directions | `stable_sort` (polars, `maintain_order`, nulls last). pandas' single-key sort is an unstable quicksort: not used in group 3. |
| comparisons (many) | `> x`, `== x`, `!= x` on floats with NaN | `compare_number`, `compare_numbers`, `compare_text` (NaN as null, filled as pandas answers: false, and true for `!=`). Float arithmetic in polars makes NaN again (`0/0`): `nan_to_null` before comparing; integer `/` must be cast to float first (polars divides integers with floor division). |
| `to_numeric(errors="coerce")` (48), `apply(axis=1)` (10), `iterrows` (4) | per-row parsing and row functions | **Own code**. `apply(axis=1)` hands each row as a Series whose values pandas may upcast; read the values from `Cell`s and check the upcasting per call. `to_numeric` parsing must be probed per call (it is not the `read_table` parser). |

## 4. Module layout

`rust/src/acmg_step/`, one module per Python script, named like the script; each thread edits
only its own files (and the one row per script in `CONVERSION.md`). Registered once: `lib.rs`
(`acmg_step`, `pandas_polars`), `cli.rs` (`priva acmg-assign`, `priva acmg-step`),
`acmg_step/mod.rs`, and the module table of `dev_harness.rs`.

| Thread | Modules |
| --- | --- |
| shared | `types` (fixed shapes: `StepOutput`, `NamedMasks`, `AcmgError`, `StepArguments`, `EffectInputs`/`QueryEffect`, `SymbolResolver`, `HarnessCall`), `stats`, `frame_file`, `dev_harness`, `resource_summary`; `crate::pandas_table`, `crate::pandas_polars` |
| T1 | `acmg_criteria_assign` (and `priva acmg-assign`), `acmg_consequence` |
| T2 | `gene_mechanism_common`, `gene_mechanism_conditions`, `gene_mechanism_resources`, `gene_mechanism_annotation`, `gene_mechanism_hub`, `clinvar_vcv` |
| T3 | `gene_mechanism_variants`, `acmg_variant_mechanism`, `acmg_scoring` |
| T4 | `acmg_pvs1_null_variant`, `alternative_start_codon`, `acmg_ps1_pm5_known_variant` |
| T5 | `acmg_pm2_bs1_ba1_frequency`, `acmg_bs2_bp5_observation`, `stats` |
| T6 | `acmg_pm1_pp2_bp1_missense`, `acmg_pp3_bp4_bp7_insilico`, `acmg_pm4_bp3_protein_length`, `acmg_ps3_bs3_pp5_bp6_reported`, `mavedb_interpreter` |
| T7 | `acmg_pp1_bs4_bp2_pm3_family`, `determine_phase`, `find_cosegregation_vars`, `acmg_ps2_pm6_denovo` |

Cross-thread signatures fixed now (stub bodies until their owner lands them):
`stats::{binom_sf, binom_cdf, beta_ppf, chi2_contingency_2x2, fisher_exact_greater}` (T5, used by
T3 and T6); `gene_mechanism_variants::infer_query_variant_effect(&EffectInputs) -> QueryEffect`
(T3, used by T2); `types::SymbolResolver` (implemented by T2's hub, used by T3);
`acmg_variant_mechanism::{variant_condition_masks, variant_mechanism_masks}(&Table) -> NamedMasks`
(T3's adapters to the existing mask modules, used by T4, T5, T7). Group 2's
`protein_domain_mapping::DomainNormalizer`, `interpro_mapping_io::load_interpro_mapping` and
`splicing_var_analysis::parse_hgvsc_splice_position` and group 1's `spliceai_cutoffs::cutoffs_for_config`
and `hpo_penetrance` are already in the crate and are called directly. Data files the Python code
finds under `scripts/..` come from `protein_domain_mapping::data_directory()` (`PRIVA_DATA_DIR`,
else the source tree's `data/`), as in group 2; no `--data-dir` option is added.

A criterion function returns `StepOutput`: its criteria (`Column::Int` strengths or `Column::Bool`,
as the Python array is; the matrix prints them as they are), the columns it writes into the table in
the Python order, and values passed on to later criteria. The orchestrator (T1) applies them.

## 5. Comparison harness

- `priva acmg-step MODULE FUNCTION --input IN.frame --output OUT.frame [--options JSON]` runs one
  converted function; every module has a `run_harness` entry (now: "not converted").
- The frame file (`acmg_step/frame_file.rs`) holds a pandas frame exactly: dtypes per column and
  every cell typed (floats by their bits, text as hex, NaN, `None` and `pd.NA` apart).
- `rust/tests/acmg_step_oracle_common.py`: `write_frame`/`read_frame`, `capture` (runs the real
  Python function, writes `expected.frame` or `expected_error.json`), `run_rust`, `compare_frames`
  (cell by cell, no tolerance), `manifest` (refuses other than Python 3.12.3, pandas 2.2.3,
  numpy 1.26.4).
- `priva acmg-step table read|write|drop-duplicates` and `priva acmg-step inspect pickle|json` test
  the table layer and the readers.
- Cluster folder `PriVA_upgrade/validation/rust_conversion_20261029/group3_base`:
  `run_group3_base_parity.sh` (template for the threads' harnesses), `table_layer_cases.py`,
  `probe_group3_resources.py`, `probe_pandas_semantics.py`, `acmg_python_run.py` (runs the Python step
  from a configuration as the shell does and prints column types), `probe_numpy_pickle.py`.

## 6. Probe results (2026-10-30, HGSC; counts and types only)

Cluster folder `PriVA_upgrade/validation/rust_conversion_20261029/group3_base`; jobs 1529197,
1529198, 1529209, 1529210, 1529212. Configuration: the GRCh38 full run
`fullrun_20261029_grch38_spliceai2_calibrated/config.yaml`.

**Resource files (plan probes P1 and P2).** Every pickle and JSON cache the configuration names was
read by Python (a plain-data unpickler that refuses every global) and by the crate's readers
(`priva acmg-step inspect`); the type summaries (type counts at every position) are compared.

| Resource | Size (unpacked) | Content | Crate reader |
| --- | ---: | --- | --- |
| `clinvar_aa_stat` (`aa_change.stats.pkl`) | 217 MB | dict transcript -> dict HGVSp -> dict -> dict {`CLNSIG`, `CLNREVSTAT`} -> list of str; 6.2 M strings | identical; 24 s, 5.7 GB (Python 2.5 GB) |
| `clinvar_splice_stat` (`splice_change.stats.pkl`) | 1,156 MB | dict transcript -> list of 5.0 M dicts {`chrom`, `pos` (int), `ref`, `alt`, `hgvsc`, `exon`, `intron`, `consequence`, `clinvar_sig`, `clinvar_review`, `splice_ai` (dict of 9 SpliceAI fields)} | identical; 104 s, **19.6 GB** (Python 7.6 GB, 248 s) |
| `clinvar_patho_af_stat` (`gene-af.pkl`) | 1.3 MB | dict gene -> dict of 12 fields or None | identical |
| `clinvar_patho_exon_af_stat` (`exon-af.pkl`) | 84 MB | dict transcript (119,344) -> dict exon -> tuple (4 floats, **numpy array**); 657,577 arrays (float64 or int64, 1-D) | **refused**: `STACK_GLOBAL` of `numpy.core.multiarray._reconstruct`, `numpy.ndarray`, `numpy.dtype` |
| `clinvar_gene_stat` (`clinvar_2star_stats.pkl.gz`) | 120 MB | dict gene -> dict of 11 fields, sets of str | identical; 12 s, 2.6 GB |
| `alphamissense_tranx_domain_map` | 131 MB | dict transcript -> {`exon_domains`, `pp_domains`} with lists and (int, int, str) tuples | identical; 14 s, 2.4 GB |
| `all_intolerant_domains` | 4.7 MB | set of 127,742 str | identical |
| `pm1_regions_pkl` | 5.6 MB | dict of 4 keys, sets of str and int | identical |
| `interpro_mapping_pickle` | 40 MB | dict (db, accession) -> list of 6-tuples | identical |
| `clingen_map` | 0.4 MB | dict variant id -> str (9,405) or float NaN (241) | identical |
| CDS cache next to the FASTA (`.cds_cache.pkl.gz`) | 490 MB | {`cds_sequences`, `transcript_to_key`}: 764,856 str -> str | identical; 2.6 s, 1.1 GB |
| `hpo_condition_mechanism_json` | 22 MB | `_meta`, `genes` (5,569 genes); no NaN/Infinity, no repeated keys, ASCII only | identical; 1.2 s |
| `gene_nonlof_mechanism_json` | 101 MB | 804 genes; 3,161 floats, booleans, nulls; ASCII only | identical; 2.0 s |

The two ClinVar dictionaries are far larger than the plan assumed (44 and 211 MiB). The exon keys
of `exon-af.pkl` are text shaped `N/N` (656,117) or `N-N/N` (1,460), as VEP's `EXON` column; its
outer keys are Ensembl and RefSeq transcript IDs.

**Column types on the real table (plan probe P3).** The BH14233 chr19 trio slice's combined table
was taken through the Python splicing and UTR steps (111 s) to the pre-ACMG table: 7,363 rows, 137
columns (object 49, int64 7, float64 71, bool 10). No object column mixes kinds of values. The table
layer reads it into a frame identical to pandas', `to_csv` and the Rust writer give byte-identical
files (identical to the input), and `drop_duplicates()` keeps the same 7,363 rows. The Python step took
218 s; its criteria matrix has 7,363 rows by 34 columns: the 28 criteria are **object columns holding
Python ints** (not int64), `variant_id` text, the 5 audit columns int64; the final table has 4,670
rows by 189 columns (object 70, all text or missing; int64 17; float64 83; bool 19). The column list
with types is `group3_base/real/acmg_python_types.txt` (sample columns masked).

**Table layer against pandas on synthetic tables** (`table_layer_cases.py`, 18 hand-written tables
and 300 random ones, the step's `read_table` option sets): 361 read cases: 349 identical frames (and for
each of them `to_csv` byte-identical and `drop_duplicates` identical on all columns and on the first
column, 698 cases); 5 that fail in both (an unclosed quote, a `usecols` name the file lacks, with
pandas' message); 7 refused on purpose (2 with an integer past the signed 64-bit range, which pandas
reads as unsigned or keeps as text; 5 with a row longer than the names, which pandas turns into the
row index); none different.

**pandas and polars behaviour** (`pandas_semantics_probe.txt`, `polars_probe_1529199.log`): the
facts of section 3. polars 0.51.0 orders NaN above every number and counts `NaN == NaN` true, divides
integers with floor division and gives null for a zero divisor, takes the second value of `zip_with`
where the mask is null, and its regex `$` does not match before a final line break. pandas 2.2.3 turns
an object True/False/NaN column into `bool` on `fillna(False)`, makes `~` of an object boolean column
`-2`/`-1`, treats NaN as false in `&`/`|` on object columns, and keeps the text of `dtype=str` values
(`007`, `1.50`).

**Dependency and compile time.** `polars-ops =0.51.0` (features `is_in`, `strings`) is added for
`is_in` and the string namespace; it brings 3 crates (`polars-ops`, `argminmax`, `libm`; `regex` was
already locked). The CI commands (`cargo build --locked --jobs 1` then `cargo test --locked --jobs 1`,
fresh target folder, HGSC node c81n-35): 296 s before, 333 s after (build 239 s to 269 s, tests
57 s to 64 s); 316 tests before, 333 after, all passing.

## 7. What each thread does first

Every thread: read sections 2 to 5, copy `run_group3_base_parity.sh` as its own harness in its own
cluster folder, capture its Python functions with `acmg_step_oracle_common.capture` on the cases of
plan section 7.2, and compare with `priva acmg-step <module> <function>`. Reuse the `pandas_polars`
helpers rather than writing new missing-value rules; a pandas call not in section 3 gets a pandas
probe before it gets code.

- **T1** (orchestrator): start with the readers the step needs that do not depend on any criterion:
  the PED (`ReadOptions { all_text: true, keep_default_na: false }`, then the headerless fallback),
  the AlphaMissense mean scores, the read and `drop_duplicates` of the table (both already identical to
  pandas on the trio slice, section 6), then `vep_consq_interpret`, `get_variant_type`,
  `truncate_fraction` against captures. `python_args` handles only options with values: add the
  `store_true` flag (`--ignore_clingen_override`) there. Owner of `pickle_reader`: the numpy-array
  extension of section 8 lands here if the owner chooses it. The orchestrator wiring comes last.
- **T2** (hub): `gene_mechanism_common`, then the condition cache path of `gene_mechanism_conditions`
  on the real 22 MB cache (the crate's JSON reader reads it identically, section 6) and the HGNC
  resolver (`read_table_with` with `all_text`; `.fillna("")` is `unwrap_or_default`). Test
  `annotate_gene_mechanism_categories` with a stub `infer_query_variant_effect` until T3 lands it.
- **T3** (exact alleles, scoring tail): first the two mask adapters in `acmg_variant_mechanism`
  (T4, T5 and T7 wait on them), then `infer_query_variant_effect` (T2 waits on it), then the matcher on
  the 101 MB non-LOF cache, then `summarize_acmg_criteria`, `compute_control_common` (with T5's
  `stats`) and the ranking glue around the existing `crate::acmg_scoring` modules.
- **T4** (PVS1, start codon, PS1/PM5): the two ClinVar dictionaries are far larger than the plan
  assumed (section 6): decode them with a visitor or into the shape the code needs rather than as a
  generic `pickle_reader::Value` tree, and measure memory first. PVS1's exon look-ups use integer and
  digit-only keys while the exon dictionary's keys are `"4/36"`-style text, so they never match in
  Python: reproduce that (section 8). The CDS cache pickle exists next to the FASTA (plan 8.8).
- **T5** (frequency, carriers, stats): `stats` first (T3 and T6 call it), with the decision tests of
  plan D6. PM2 and BS1 read the exon allele-frequency pickle, which holds numpy arrays: blocked on
  the owner's decision in section 8 for real-data runs (synthetic captures can start).
- **T6** (missense, predictors, length, reported): PP3/BP4/BP7 first (plain columns), then PM1 (plain
  `pm1_regions` pickle), PP2/BP1 (the 120 MB `clinvar_2star_stats` pickle of dicts and sets reads
  identically) once `stats` lands, then PM4/BP3 and the MaveDB rules.
- **T7** (family): PS2/PM6 and `identify_fam_members` on the PED reader, then `determine_phase` with the
  fixed merge order (`first_appearance_groups` over `Gene`), then BS4 and BP2/PM3; PP1 and
  `find_cosegregation_vars` last, with `--pp1_vcf`/`--pp1_ped` forwarded.

## 8. Open

1. **numpy arrays in `clinvar_patho_exon_af_stat`** (owner's decision, plan 8.7). The step reads
   only positions 0 and 2 of each tuple, never the arrays. Options: (a) extend `pickle_reader` to
   accept exactly the three numpy globals with the `_reconstruct`/`BUILD` pattern of this file
   (1-D, C-ordered, little-endian float64 and int64; the array becomes an opaque value), about 150
   lines with tests against this file; (b) rebuild the cache with plain tuples (changes a cache
   identity). Recommendation: (a), done by T1 as owner of the reader, before T5 runs PM2/BS1 on real
   data (T4's PVS1 reads the same file).
2. **Memory of the generic pickle reader on the ClinVar splice dictionary** (19.6 GB against
   Python's 7.6 GB). T4 should decode the two ClinVar dictionaries straight into the shape PS1/PM5
   and `summarize_clinvar_gene_pathogenicity` need (a visitor on the opcode stream) and keep the job's
   memory request in mind.
3. **PVS1 exon look-ups never match in Python** (`downstream_exon_patho_af` and the start-lost and
   splicing branches look up integer or digit-only keys; the keys are `N/N` text), so they always see
   the default `(0,)`. Reproduced faithfully unless the owner wants a fix in Python and Rust together.
4. The criteria matrix's criteria are object columns of Python ints; they print as integers, and
   T3/T1 must keep that type where the matrix is used again (posterior, `isin([0..5])`).
5. `rust/source_functions.tsv` still pins `acmg_scoring.py` at an old digest (plan 4.4); re-pin it
   with T3's first commit that touches `acmg_scoring`.
6. Probe P4 (how many variants get different phase calls under two genes) was not run: the merge
   order is now fixed by decision, so the count no longer changes the design. Probes P5 to P7 belong
   to T5 (stats, `bcftools`) and T6 (`bedtools`); the reader parts of P7 were done in groups 1 and 2.
7. Single-column tables: the writer now quotes a lone missing value (`""`) as `csv.writer` does
   (checked on the synthetic tables); earlier converted steps never write one-column tables.
