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
| `scripts/spliceai_cache_input.py` | 703 | Python |
| `scripts/spliceai_cutoffs.py` | 116 | Python |

### 2. Steps that interpret splicing and UTR effects for each sample

4 files, 2,252 lines.

| Script | Lines | Status |
| --- | ---: | --- |
| `scripts/interpro_mapping_io.py` | 54 | Python |
| `scripts/protein_domain_mapping.py` | 764 | Python |
| `scripts/splicing_var_analysis.py` | 1,156 | Python |
| `scripts/utr_anno_interpret.py` | 278 | Python |

### 3. ACMG decision steps (criterion assignment, scoring and ranking for each sample)

26 files, 13,704 lines.

| Script | Lines | Status |
| --- | ---: | --- |
| `scripts/acmg_bs2_bp5_observation.py` | 460 | Python |
| `scripts/acmg_consequence.py` | 196 | Python |
| `scripts/acmg_criteria_assign.py` | 711 | Python |
| `scripts/acmg_inheritance.py` | 190 | Python |
| `scripts/acmg_pm1_pp2_bp1_missense.py` | 307 | Python |
| `scripts/acmg_pm2_bs1_ba1_frequency.py` | 799 | Python |
| `scripts/acmg_pm4_bp3_protein_length.py` | 224 | Python |
| `scripts/acmg_pp1_bs4_bp2_pm3_family.py` | 602 | Python |
| `scripts/acmg_pp3_bp4_bp7_insilico.py` | 139 | Python |
| `scripts/acmg_ps1_pm5_known_variant.py` | 288 | Python |
| `scripts/acmg_ps2_pm6_denovo.py` | 215 | Python |
| `scripts/acmg_ps3_bs3_pp5_bp6_reported.py` | 258 | Python |
| `scripts/acmg_pvs1_null_variant.py` | 721 | Python |
| `scripts/acmg_scoring.py` | 836 | Python |
| `scripts/acmg_variant_mechanism.py` | 591 | Python |
| `scripts/alternative_start_codon.py` | 405 | Python |
| `scripts/clinvar_vcv.py` | 1,061 | Python |
| `scripts/determine_phase.py` | 518 | Python |
| `scripts/find_cosegregation_vars.py` | 228 | Python |
| `scripts/gene_mechanism_annotation.py` | 411 | Python |
| `scripts/gene_mechanism_common.py` | 431 | Python |
| `scripts/gene_mechanism_conditions.py` | 1,060 | Python |
| `scripts/gene_mechanism_hub.py` | 541 | Python |
| `scripts/gene_mechanism_resources.py` | 1,028 | Python |
| `scripts/gene_mechanism_variants.py` | 642 | Python |
| `scripts/mavedb_interpreter.py` | 842 | Python |

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
