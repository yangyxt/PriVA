# Converting PriVA's Python scripts to Rust

Goal (decided 2026-10-29): every customized PriVA Python script becomes part of the `priva` Rust
binary built from this folder, except PriVA's replacements inside upstream CADD and the one-off
research tools and examples, which stay in Python. This file lists every script, the order of
conversion and its status.

## How each script is converted, in time order

1. Write a `priva` subcommand that takes the same arguments and files as the Python script and
   reproduces what the script actually does, including quirks of the libraries it uses (for
   example pysam never reporting INFO END).
2. Run both on the same inputs: small synthetic files in `cargo test`, then public data (ClinVar)
   and the BH14233 trio slice on HGSC. Outputs must be identical after normalising VCFs through
   htslib; tables are compared field by field with `bin/priva-compare`.
3. Switch the shell or Nextflow caller to the subcommand only after step 2 passes.
4. Delete the Python script once nothing calls it and its group is complete.

Groups are converted in pipeline order, so each converted step can be checked on the real
intermediate files of the step before it.

## Scripts

### 1. Steps that prepare each sample's variants (variant IDs, VCF to table, CADD hub cache, family filtering, SpliceAI version check)

7 files, 2,883 lines.

| Script | Lines | Status |
| --- | ---: | --- |
| `scripts/cadd_hub_util.py` | 533 | Rust `priva cadd-hub`: identical files and exit codes on 20 synthetic cases covering every match, update and merge path, on 2,032,185 public CADD GRCh38 v1.7 rows (chr19:1.0-1.3 Mb) and on the BH14233 trio slice's CADD table (2026-10-29). Unlike the Python script's polars, it also runs on HGSC nodes without AVX2. The shell still calls Python until group 1 is complete. |
| `scripts/combine_annotations.py` | 936 | Rust `priva combine-annotations`: identical tables and exit codes on 9 synthetic cases (22,956 rows over two chunks with columns that change type, three HPO file layouts, missing sex fields, a header without the INFO keys the script reads, and the error exits), the BH14233 trio slice (7,363 rows by 124 columns, byte-identical with one worker and cell for cell identical to the script's 4-worker run as sorted rows) and public ClinVar GRCh38 chr22 with made-up genotypes (3,586,969 rows by 121 columns; 460 s against 1,792 s) (2026-10-30). With more than one worker the script returns records in a varying order; the Rust version keeps the VCF's order. The shell still calls Python until group 1 is complete. |
| `scripts/generate_varID_for_vcf.py` | 150 | Rust `priva variant-ids`: byte-identical output, encode and decode, on pysam edge cases, 4,597,985 public ClinVar records and the BH14233 trio slice (2026-10-29). The shell still calls Python until group 1 is complete. |
| `scripts/hpo_penetrance.py` | 143 | Rust module `hpo_penetrance`, used by `priva combine-annotations` and checked with it (2026-10-30). The Python module stays until its callers are Rust. |
| `scripts/pedigree_filtration_per_fam.py` | 302 | Rust `priva pedigree-filter`: identical outputs and exit codes on 20 synthetic cases (every GT combination of a four-person family on chr1, chrX and chrY under 13 pedigree layouts, alleles beyond a record's ALT alleles, the error cases and a rerun on its own output), the repository's TEST_FAM, the BH14233 trio slice and 4,597,985 public ClinVar GRCh38 records given made-up genotypes (2026-10-29). Like the Python script under pysam, it stops when the pedigree names a family member who is not a VCF sample. The shell still calls Python until group 1 is complete. |
| `scripts/spliceai_cache_input.py` | 703 | Rust `priva spliceai-cache-input` (`install`, `adopt-prepared`; also `inspect_header` and `canonical_assembly` as library functions for the ACMG group). Compared with the Python script on 309 synthetic cases (header classification, YAML configurations, SpliceAI v1 pairs with real tabix indexes, completed and relocated SpliceAI2 folders, native tables with a stand-in helper script, `adopt-prepared` with resumed, occupied and refused copies, JSON receipts with odd numbers, command-line forms) by exit code, printed JSON, error text and every file left behind (bytes, permissions, times): 299 identical, 6 identical except for PyYAML's own wording of its syntax errors, 4 deliberate refusals, 0 different; on public ClinVar slices standing in for SpliceAI v1 caches (bgzip and tabix) the registered configuration and the report are identical (2026-10-30). The command line is the script's argparse parser reproduced (options, abbreviations, usage, help, errors) under the program name `priva spliceai-cache-input`. Differences by design: YAML that PyYAML reads but the flat reader does not (anchors, aliases, tags, block scalars, multi-line scalars, keys with spaces, an indented top-level mapping) is refused with exit code 1 and a message saying so; the real `config.yaml` files hold none; the folder holding `install_utils.sh` comes from the environment variable `PRIVA_SCRIPT_DIR` (no extra option); copied files do not get extended attributes; zlib's reason for a corrupt deflate stream is generic. The shell still calls Python until group 1 is complete. |
| `scripts/spliceai_cutoffs.py` | 116 | Rust `priva spliceai-cutoffs CONFIG [SNV_VCF INDEL_VCF]` (and `cutoffs_for_config` as a library function). Compared with the Python script on 366 synthetic cases (every header form with each assembly setting, configuration file forms including byte order mark, `---` line, quoted keys and tabs, score files that are missing, mixed or damaged, command-line forms) by exit code, printed text and error text: 344 identical, 14 identical except for PyYAML's or zlib's own wording of a diagnostic, 8 deliberate refusals, 0 different; identical on the repository `config.yaml`, on the cluster run's `config.yaml` with its real SpliceAI2 score files (2026-10-30). Deliberate refusals, exit code 1 with a message: the YAML constructs listed for `priva spliceai-cache-input`. The shell still calls Python until group 1 is complete. |

### 2. Steps that interpret splicing and UTR effects for each sample

4 files, 2,252 lines.

| Script | Lines | Status |
| --- | ---: | --- |
| `scripts/interpro_mapping_io.py` | 54 | Rust `priva interpro-mapping inspect` and the pickle reader that `priva splicing-analysis` uses: identical output and exit codes on 15 of 16 cases (the bundled mapping, a synthetic one, plain and gzip files, a list, a set, a dict of transcripts, option forms, missing and damaged files); the 16th, a pickle that holds a Python class object, is refused with exit 1 where Python loads it; the help text differs by design (2026-10-30). The reader accepts plain data only and names the opcode it refuses. The shell still calls Python until group 2 is complete. |
| `scripts/protein_domain_mapping.py` | 764 | Python for its installer command line (`--vcf --xml --report --mapping_output`, which builds the InterPro pickle). Its lookup class `DomainNormalizer` (VEP domain to InterPro entry, functional or not) is Rust library code in `rust/src/protein_domain_mapping.rs`, used by `priva splicing-analysis` and checked through it (4,000-row tables against the bundled InterPro mapping and a synthetic one, the BH14233 trio slice) and by unit tests (2026-10-30). |
| `scripts/splicing_var_analysis.py` | 1,156 | Rust `priva splicing-analysis`: identical tables and exit codes on 67 synthetic cases (4,000 rows against the bundled InterPro mapping, awkward pandas columns, every error exit, argument forms, SpliceAI v1 and v2 configurations) and 10 command-line forms, and a byte-identical table on the BH14233 chr19 trio slice (7,363 rows, 133 columns; 16 s against 112 s) (2026-10-30). Differences: `splicing_affected_exons` is written in ascending order (Python's set order varies from run to run; the comparison treats it as a set), and two inputs Python accepts are refused with exit 1 (an InterPro mapping that is a plain list; a transcript domain map that holds Python class objects). It reads and writes the table the way pandas does. The shell still calls Python until group 2 is complete. |
| `scripts/utr_anno_interpret.py` | 278 | Rust `priva utr-anno-interpret`, on the crate's shared pandas-style table reader and plain-data pickle reader: identical rewritten tables and exit codes on 270 synthetic cases (a zoo of column types and odd numbers including 2,700 random number columns, table layouts, every effect and label pairing with more than twenty shapes of domain pickle, bad annotations and `DistanceToCDS` spellings, unreadable pickles and tables), 20 command line forms, the repository example's real transcript domain map, and the BH14233 trio slice's combined table (7,363 rows by 124 columns, byte-identical, and identical again when run on its own output; none of its rows has a UTR effect that scores true, so the true branches are checked on the synthetic and example tables only) (2026-10-30). Refused on purpose with exit 1 and the table left untouched, where the script succeeds (23 of the 293 cases): a table named like a compressed file, a data row longer than the header, numbers past the signed 64-bit range that pandas reads as unsigned integers or keeps as text, pickles that need more than plain data (protocol 2 sets, `defaultdict`), and a domain map that is not a dict. Making the output match needed fixes in the shared readers (booleans in any capitalisation, 17-digit exponents that wrap like pandas' 32-bit integer, the smallest int64 as pandas' missing marker, exact ties in float printing). The shell still calls Python until group 2 is complete. |

### 3. ACMG decision steps (criterion assignment, scoring and ranking for each sample)

26 files, 13,704 lines.

Shared base laid down on 2026-10-30 before seven threads convert these scripts in parallel:
the module skeleton `rust/src/acmg_step/` (one module per script, `priva acmg-assign` registered,
`priva acmg-step` to compare one function with Python on exact pandas frames), the table layer
(`pandas_table` read options of the step, the `object`, `Int64` and `boolean` column types, and
`pandas_polars`, which moves typed columns to polars where polars gives pandas' result), the fixed
cross-thread signatures and the probes of the real resource files. `rust/GROUP3_BASE.md` records the
design, the call-by-call choice between polars and own code, the probe results and where each
thread starts.

| Script | Lines | Status |
| --- | ---: | --- |
| `scripts/acmg_bs2_bp5_observation.py` | 460 | Python |
| `scripts/acmg_consequence.py` | 196 | Python |
| `scripts/acmg_criteria_assign.py` | 711 | Python |
| `scripts/acmg_inheritance.py` | 190 | Python |
| `scripts/acmg_pm1_pp2_bp1_missense.py` | 307 | Rust module `acmg_pm1_pp2_bp1_missense` (`pm1_criteria`, `pp2_bp1_criteria`, `analyze_bp1_pp2`), compared with the Python functions through `priva acmg-step` (2026-10-30): PM1 on 89 synthetic cases (85 identical, 4 failing in both; plain and gzip region pickles, sets, frozensets and lists of residues, odd `Protein_position` texts and types); PP2/BP1 on 50 synthetic cases, all identical; `analyze_bp1_pp2` for all 10,016 genes of the installed ClinVar gene statistics, identical; the BH14233 trio slice identical for both (PM1 18, PP2 1, BP1 0). Waits for `stats.rs` (thread T5): the 2x2 chi-square and Fisher p-values (598 distinct tables over the installed statistics, 1 gene on the trio slice) were supplied from SciPy in these comparisons; the step itself calls `stats::chi2_contingency_2x2` and `stats::fisher_exact_greater`. Refused: two genes whose statistics name one `ensembl_id` with different results (Python's answer depends on process completion order), non-ASCII digits in `Protein_position`. |
| `scripts/acmg_pm2_bs1_ba1_frequency.py` | 799 | Python |
| `scripts/acmg_pm4_bp3_protein_length.py` | 224 | Rust module `acmg_pm4_bp3_protein_length` (2026-10-30): the repeat-region test is native (regions sorted per chromosome with a running maximum of ends) instead of pybedtools; identical to `bedtools intersect -wa -f 1 -u` 2.31.1 on 20,000 intervals placed around region edges of the installed GRCh38 repeat BED. PM4/BP3 on 90 synthetic cases: 81 identical, 9 failing in both; the BH14233 trio slice identical (PM4 17, BP3 0, 3,242 variants in repeats). Refused: empty regions (bedtools widens them), positions below 1, an empty or missing `ref`, names `to_csv` would quote. |
| `scripts/acmg_pp1_bs4_bp2_pm3_family.py` | 602 | Python |
| `scripts/acmg_pp3_bp4_bp7_insilico.py` | 139 | Rust module `acmg_pp3_bp4_bp7_insilico` (2026-10-30), with the boolean-mask rules the T6 modules share (pandas' `&`/`|` on `object` masks holding NaN, `np.logical_not`, `fillna` downcasting). PP3/BP4 on 178 synthetic cases: 166 identical, 11 failing in both, 1 refused (a text `PrimateAI` column, which `pd.to_numeric` would parse); BP7 on 179: 174 identical, 5 failing in both; the BH14233 trio slice identical (PP3 1, BP4 587, BP7 11). Refused: number columns read through `to_numeric` that are text, `~` of an `object` boolean Series. |
| `scripts/acmg_ps1_pm5_known_variant.py` | 288 | Python |
| `scripts/acmg_ps2_pm6_denovo.py` | 215 | Python |
| `scripts/acmg_ps3_bs3_pp5_bp6_reported.py` | 258 | Rust module `acmg_ps3_bs3_pp5_bp6_reported` (2026-10-30), including the frame `df.apply(axis=1)` returns (its `infer_objects` re-typing is reproduced). PS3/BS3 on 128 synthetic cases with the shipped and a synthetic MaveDB metadata table: 114 identical, 12 failing in both, 2 refused (an assay family that is missing, which Python groups by NaN identity); PP5/BP6 on 120, all identical; the BH14233 trio slice identical for both (the slice has no MaveDB URNs; BP6 93). Refused: empty tables and tables without an `object` column on the MaveDB path, mixed or nullable columns there. |
| `scripts/acmg_pvs1_null_variant.py` | 721 | Python |
| `scripts/acmg_scoring.py` | 836 | Partly in Rust for the ACMG step (`rust/src/acmg_step/acmg_scoring.rs`, group 3 thread T3): `summarize_acmg_criteria` (matrix, validation, ClinGen map, summary column and override, with the column types pandas leaves, which change to `object` in the ways pandas changes them) and the orchestrator's posterior step are identical to Python on 150 and 300 synthetic cases (109 identical frames and 41 identical errors; 300 identical) and on the BH14233 chr19 trio slice (7,363 rows, the Python step's own criteria and the real ClinGen map) (2026-10-30). `compute_control_common` is identical on the 37 synthetic cases that need no binomial or beta function (including the IndexError Python raises when the frequency-only route covers some rows and flags one) and waits for thread T5's `stats` for the rest. `sort_and_rank_variants` (R1 to R8) is not converted yet. Refused: criteria that are neither int64 nor bool, empty tables, a ClinGen map path that `repr()` would escape. The step still runs in Python. |
| `scripts/acmg_variant_mechanism.py` | 591 | In Rust for the ACMG step (`rust/src/acmg_step/acmg_variant_mechanism.rs`, group 3 thread T3): the two mask adapters (over the earlier mask modules) are identical to Python on 305 + 305 synthetic frames (identical masks or the same error) and on the trio slice after the hub; `annotate_exact_nonlof_variants` is identical on 100 synthetic caches with tables and on the trio slice with the real non-LOF cache (7,363 rows; no row of the slice matches a curated allele) (2026-10-30). Refused: booleans in object or nullable score columns, numbers past 2**53 in a score column, an unreadable or malformed mechanism JSON (Python skips every lookup then). The step still runs in Python. |
| `scripts/alternative_start_codon.py` | 405 | Python |
| `scripts/clinvar_vcv.py` | 1,061 | Python |
| `scripts/determine_phase.py` | 518 | Python |
| `scripts/find_cosegregation_vars.py` | 228 | Python |
| `scripts/gene_mechanism_annotation.py` | 411 | Python |
| `scripts/gene_mechanism_common.py` | 431 | Python |
| `scripts/gene_mechanism_conditions.py` | 1,060 | Python |
| `scripts/gene_mechanism_hub.py` | 541 | Python |
| `scripts/gene_mechanism_resources.py` | 1,028 | Python |
| `scripts/gene_mechanism_variants.py` | 642 | In Rust for the ACMG step (`rust/src/acmg_step/gene_mechanism_variants.rs`, group 3 thread T3), except the two compatibility views `match_gofcards_variant_gof*` that the step never calls: `infer_query_variant_effect` is identical on 304 synthetic frames (302 identical, 2 the same OverflowError) and on the trio slice's 7,363 rows; `match_curated_nonlof_variant` is identical on 100 synthetic caches and on all 34,314 keys of the real non-LOF cache's three indexes (protein change, coding change, genomic allele) (2026-10-30). Refused: score or position texts that `float()` might read through non-ASCII digits, cache values of kinds the cache never holds (lists where text is expected, non-text PMIDs). The step still runs in Python. |
| `scripts/mavedb_interpreter.py` | 842 | Rust module `mavedb_interpreter` (`pick_interpreter_from_metadata` and the 13 family interpreters; the free-text `notes` are not built): call, standardized score, reason, PS3/BS3 flags and raised errors identical for 61,561 queries over every row of the shipped metadata (scores at and next to each band edge and percentile, zeros, infinities, NaN) and 2,873 queries over 122 synthetic rows with odd bands, flags, families and enum names (2026-10-30). |

### 4. Scripts that build reference data once, at installation

31 files, 20,058 lines.

| Script | Lines | Status |
| --- | ---: | --- |
| `scripts/SpliceVault_tsv_vcf_conversion.py` | 192 | Python |
| `scripts/alphamissense_annotation_base.py` | 611 | Python |
| `scripts/alphamissense_cache.py` | 363 | Python |
| `scripts/alphamissense_domain_cache.py` | 620 | Python |
| `scripts/alphamissense_original_cache.py` | 368 | Python |
| `scripts/alphmis_tsv2vcf.py` | 85 | Python |
| `scripts/am_pick_intolerant_domains.py` | 1,748 | Python |
| `scripts/am_pick_missense_constrained_segments.py` | 1,275 | Python |
| `scripts/audit_clinvar_gofcards_review_tiers.py` | 154 | Python |
| `scripts/audit_condition_scope_inventory.py` | 868 | Python |
| `scripts/audit_gofcards_source_coverage.py` | 356 | Python |
| `scripts/biomart_ftp_nm_to_enst.py` | 556 | Python |
| `scripts/biomart_nm_to_enst.py` | 628 | Python |
| `scripts/build_gene_nonlof_mechanism_cache.py` | 2,430 | Python |
| `scripts/build_gene_pathogenic_mechanism_cache.py` | 1,868 | Python |
| `scripts/build_gofcards_exact_gof_cache.py` | 1,061 | Python |
| `scripts/build_hpo_condition_mechanism_cache.py` | 2,102 | Python |
| `scripts/build_mondo_disease_scope.py` | 529 | Python |
| `scripts/clinvar_annotation_contract.py` | 244 | Python |
| `scripts/clinvar_gene_refresh.py` | 1,201 | Python |
| `scripts/clinvar_stat_variants.py` | 373 | Python |
| `scripts/collapse_HPO_anno.py` | 268 | Python |
| `scripts/domain_metadata_cache.py` | 168 | Python |
| `scripts/inject_clinvar_into_gofcards.py` | 232 | Python |
| `scripts/prepare_clingen_map.py` | 64 | Python |
| `scripts/prepare_pm1_regions.py` | 397 | Python |
| `scripts/resort_step3_only.py` | 76 | Python |
| `scripts/stat_aachange_clinvar.py` | 371 | Python |
| `scripts/stat_exon_patho_afs.py` | 238 | Python |
| `scripts/stat_gene_patho_afs.py` | 216 | Python |
| `scripts/stat_protein_domain_amscores.py` | 396 | Python |

### 5. Installer and cluster deployment helpers

22 files, 6,211 lines.

| Script | Lines | Status |
| --- | ---: | --- |
| `scripts/cache_resources.py` | 159 | Python |
| `scripts/check_cluster_environment.py` | 87 | Python |
| `scripts/check_cluster_readiness.py` | 133 | Python |
| `scripts/conservation_resource_check.py` | 212 | Python |
| `scripts/das_input_continuity.py` | 212 | Python |
| `scripts/das_install_state.py` | 140 | Python |
| `scripts/deploy_cluster_gnomad.py` | 35 | Python |
| `scripts/deploy_cluster_plugin_inputs.py` | 67 | Python |
| `scripts/deploy_cluster_resources.py` | 190 | Python |
| `scripts/gnomad_cache.py` | 174 | Python |
| `scripts/indexed_plugin_cache.py` | 306 | Python |
| `scripts/installed_environment.py` | 242 | Python |
| `scripts/installer_config_paths.py` | 899 | Python |
| `scripts/installer_prerequisites.py` | 532 | Python |
| `scripts/mapped_artifact_publish.py` | 370 | Python |
| `scripts/mechanism_source_input.py` | 155 | Python |
| `scripts/prepare_cluster_config.py` | 69 | Python |
| `scripts/reference_resources.py` | 470 | Python |
| `scripts/spliceai2_to_vep_vcf.py` | 544 | Python |
| `scripts/vep_chromosome_cache.py` | 450 | Python |
| `scripts/vep_chromosome_jobs.py` | 375 | Python |
| `scripts/vep_plugins_release.py` | 390 | Python |

### 6. Nextflow task adapters

2 files, 786 lines.

| Script | Lines | Status |
| --- | ---: | --- |
| `bin/priva_task.py` | 673 | Python |
| `bin/priva_task_boundary.py` | 113 | Python |

### Removed when their callers are Rust: Python-to-Rust bridges

10 files, 2,687 lines.

| Script | Lines | Status |
| --- | ---: | --- |
| `scripts/bootstrap_priva_rust.py` | 662 | Python |
| `scripts/priva_rust_control_multiply.py` | 117 | Python |
| `scripts/priva_rust_criteria_summary.py` | 278 | Python |
| `scripts/priva_rust_gene_weights.py` | 257 | Python |
| `scripts/priva_rust_pedigree_context.py` | 184 | Python |
| `scripts/priva_rust_pext.py` | 217 | Python |
| `scripts/priva_rust_posterior.py` | 381 | Python |
| `scripts/priva_rust_rank_setup.py` | 150 | Python |
| `scripts/priva_rust_ranking_compatibility.py` | 232 | Python |
| `scripts/priva_rust_terminal_order.py` | 209 | Python |

### Kept in Python: PriVA's replacements inside upstream CADD

4 files, 795 lines.

| Script | Lines | Status |
| --- | ---: | --- |
| `scripts/cadd_prescore_index.py` | 172 | Stays in Python |
| `scripts/cadd_scripts/MMSplice.py` | 256 | Stays in Python |
| `scripts/cadd_scripts/extract_scored.py` | 272 | Stays in Python |
| `scripts/cadd_scripts_install.py` | 95 | Stays in Python |

### Kept in Python: one-off research tools and examples

5 files, 1,241 lines.

| Script | Lines | Status |
| --- | ---: | --- |
| `data/InterPro/gemini_interpro_func_classifier.py` | 425 | Stays in Python |
| `data/MaveDB/cal_mavedb_neutral_bands.py` | 272 | Stays in Python |
| `data/MaveDB/gpt_mavedb_assay_mapper.py` | 282 | Stays in Python |
| `data/MaveDB/pretify_json.py` | 55 | Stays in Python |
| `example/smoke_test.py` | 207 | Stays in Python |

## Kept in Python (decided 2026-10-29)

These two groups stay in Python; every other script above is converted.

- `scripts/cadd_scripts/MMSplice.py` runs the MMSplice neural network (TensorFlow) inside upstream
  CADD's own environment, and `extract_scored.py` replaces a script inside upstream CADD's workflow.
  Converting them means re-implementing parts of CADD, which runs as upstream `CADD.sh`.
- `data/InterPro/gemini_interpro_func_classifier.py` and `data/MaveDB/gpt_mavedb_assay_mapper.py` call
  external AI services by hand to produce committed data files; the other files in that group are
  one-off data tools and the smoke-test example.
