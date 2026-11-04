# PriVA native procedures and source map

`priva` groups the pipeline's native procedures in one executable. The shell
orchestration still invokes mature external tools, including official VEP and
upstream `CADD.sh` with its Python scorer. Native cache builders and production
Nextflow migration are separate from the converted runtime procedures.

```text
Shell / task adapter
    ├── priva: preparation → splicing / UTR → ACMG → ranking
    ├── official VEP and other external tools
    └── upstream CADD.sh and its Python environments

Installer
    ├── retained Python builders and optional native builders
    └── cache-install → checked native caches → runtime readers
```

## Invocation and formats

`scripts/common_bash_utils.sh` resolves `PRIVA_BINARY`, then the YAML
`priva_binary`, then `PATH`, and sets the resource directories for the executable.
Use each command's `--help` for its supported input profile and named options.
The shell runtime calls native procedures for preparation, interpretation and
ACMG assignment. An optional builder listed below does not imply that its
installer caller has switched from Python.

Keyed caches use checked `PRIVARKV`/rkyv archives. Builder tables use
Polars/Parquet-zstd; naturally shared text formats remain text. See
[CACHE_FORMATS.md](CACHE_FORMATS.md) for schema versions and conversion commands.
The bundled ClinVar gene pickle is converted under configured `cache_root`
during installation; its larger native archive is not a bundled Git resource.
There is no native pickle writer.

Preserve scientific execution and mutation order, types, missingness, joins and
numerical/failure behavior within the selected profile. Unsupported input
profiles must be explicit. Reusable comparisons live in `tests/` and
`rust/tests/`; run outputs belong outside the source checkout.

## Runtime source map

The tables identify procedure ownership, not whole-pipeline qualification.
Native scientific differences, including residue-aware PM5 and frameshift
handling, remain profile-specific. An available module does not make its
standalone diagnostic or cache-building entry point native.

### Steps that prepare each sample's variants (variant IDs, VCF to table, CADD hub cache, family filtering, SpliceAI version check)

| Python source | Native procedure or current role |
| --- | --- |
| [`scripts/cadd_hub_util.py`](../scripts/cadd_hub_util.py) | `priva cadd-hub`. |
| [`scripts/combine_annotations.py`](../scripts/combine_annotations.py) | `priva combine-annotations`. |
| [`scripts/generate_varID_for_vcf.py`](../scripts/generate_varID_for_vcf.py) | `priva variant-ids`. |
| [`scripts/hpo_penetrance.py`](../scripts/hpo_penetrance.py) | Native `hpo_penetrance`, used by `combine-annotations`. |
| [`scripts/pedigree_filtration_per_fam.py`](../scripts/pedigree_filtration_per_fam.py) | `priva pedigree-filter`. |
| [`scripts/spliceai_cache_input.py`](../scripts/spliceai_cache_input.py) | `priva spliceai-cache-input`. |
| [`scripts/spliceai_cutoffs.py`](../scripts/spliceai_cutoffs.py) | `priva spliceai-cutoffs`. |

### Steps that interpret splicing and UTR effects for each sample

| Python source | Native procedure or current role |
| --- | --- |
| [`scripts/interpro_mapping_io.py`](../scripts/interpro_mapping_io.py) | `priva interpro-mapping inspect`. |
| [`scripts/protein_domain_mapping.py`](../scripts/protein_domain_mapping.py) | `priva interpro-mapping build`. XML-to-mapping builder only; standalone VCF/report caller remains Python. |
| [`scripts/splicing_var_analysis.py`](../scripts/splicing_var_analysis.py) | `priva splicing-analysis`. |
| [`scripts/utr_anno_interpret.py`](../scripts/utr_anno_interpret.py) | `priva utr-anno-interpret`. |

### ACMG decision steps (criterion assignment, scoring and ranking for each sample)

| Python source | Native procedure or current role |
| --- | --- |
| [`scripts/acmg_bs2_bp5_observation.py`](../scripts/acmg_bs2_bp5_observation.py) | Native runtime `acmg_step::acmg_bs2_bp5_observation` through `acmg-assign`. |
| [`scripts/acmg_consequence.py`](../scripts/acmg_consequence.py) | Native runtime `acmg_step::acmg_consequence` through `acmg-assign`. |
| [`scripts/acmg_criteria_assign.py`](../scripts/acmg_criteria_assign.py) | `priva acmg-assign`; runtime orchestration. |
| [`scripts/acmg_inheritance.py`](../scripts/acmg_inheritance.py) | Native runtime `acmg_step::acmg_inheritance` through `acmg-assign`. |
| [`scripts/acmg_pm1_pp2_bp1_missense.py`](../scripts/acmg_pm1_pp2_bp1_missense.py) | Native runtime `acmg_step::acmg_pm1_pp2_bp1_missense` through `acmg-assign`. |
| [`scripts/acmg_pm2_bs1_ba1_frequency.py`](../scripts/acmg_pm2_bs1_ba1_frequency.py) | Native runtime `acmg_step::acmg_pm2_bs1_ba1_frequency` through `acmg-assign`. |
| [`scripts/acmg_pm4_bp3_protein_length.py`](../scripts/acmg_pm4_bp3_protein_length.py) | Native runtime `acmg_step::acmg_pm4_bp3_protein_length` through `acmg-assign`. |
| [`scripts/acmg_pp1_bs4_bp2_pm3_family.py`](../scripts/acmg_pp1_bs4_bp2_pm3_family.py) | Native runtime `acmg_step::acmg_pp1_bs4_bp2_pm3_family` through `acmg-assign`. |
| [`scripts/acmg_pp3_bp4_bp7_insilico.py`](../scripts/acmg_pp3_bp4_bp7_insilico.py) | Native runtime `acmg_step::acmg_pp3_bp4_bp7_insilico` through `acmg-assign`. |
| [`scripts/acmg_ps1_pm5_known_variant.py`](../scripts/acmg_ps1_pm5_known_variant.py) | Native runtime `acmg_step::acmg_ps1_pm5_known_variant` through `acmg-assign`. |
| [`scripts/acmg_ps2_pm6_denovo.py`](../scripts/acmg_ps2_pm6_denovo.py) | Native runtime `acmg_step::acmg_ps2_pm6_denovo` through `acmg-assign`. |
| [`scripts/acmg_ps3_bs3_pp5_bp6_reported.py`](../scripts/acmg_ps3_bs3_pp5_bp6_reported.py) | Native runtime `acmg_step::acmg_ps3_bs3_pp5_bp6_reported` through `acmg-assign`. |
| [`scripts/acmg_pvs1_null_variant.py`](../scripts/acmg_pvs1_null_variant.py) | Native runtime `acmg_step::acmg_pvs1_null_variant` through `acmg-assign`. |
| [`scripts/acmg_scoring.py`](../scripts/acmg_scoring.py) | Native runtime `acmg_step::acmg_scoring` through `acmg-assign`. |
| [`scripts/acmg_variant_mechanism.py`](../scripts/acmg_variant_mechanism.py) | Native runtime `acmg_step::acmg_variant_mechanism` through `acmg-assign`. |
| [`scripts/alternative_start_codon.py`](../scripts/alternative_start_codon.py) | Native `acmg_step::alternative_start_codon` runtime detector; standalone CDS builder is outside conversion scope. |
| [`scripts/clinvar_vcv.py`](../scripts/clinvar_vcv.py) | Native `acmg_step::clinvar_vcv` text input; installer XML parsing remains Python. |
| [`scripts/determine_phase.py`](../scripts/determine_phase.py) | Native runtime `acmg_step::determine_phase` through `acmg-assign`. |
| [`scripts/find_cosegregation_vars.py`](../scripts/find_cosegregation_vars.py) | Native runtime `acmg_step::find_cosegregation_vars` through `acmg-assign`. |
| [`scripts/gene_mechanism_annotation.py`](../scripts/gene_mechanism_annotation.py) | Native runtime `acmg_step::gene_mechanism_annotation` through `acmg-assign`. |
| [`scripts/gene_mechanism_common.py`](../scripts/gene_mechanism_common.py) | Native runtime `acmg_step::gene_mechanism_common` through `acmg-assign`. |
| [`scripts/gene_mechanism_conditions.py`](../scripts/gene_mechanism_conditions.py) | Native runtime `acmg_step::gene_mechanism_conditions`; installation and diagnostic helpers remain Python. |
| [`scripts/gene_mechanism_hub.py`](../scripts/gene_mechanism_hub.py) | Native runtime facade over mechanism modules; standalone diagnostics remain Python. |
| [`scripts/gene_mechanism_resources.py`](../scripts/gene_mechanism_resources.py) | Native runtime `acmg_step::gene_mechanism_resources`; installation and diagnostic helpers remain Python. |
| [`scripts/gene_mechanism_variants.py`](../scripts/gene_mechanism_variants.py) | Native runtime `acmg_step::gene_mechanism_variants` through `acmg-assign`. |
| [`scripts/mavedb_interpreter.py`](../scripts/mavedb_interpreter.py) | Native runtime `acmg_step::mavedb_interpreter` through `acmg-assign`. |

## Installation builders

Installer callers remain Python unless their row says the installer calls
`priva`, or they are routed through the native cache conversion interface. Optional native procedures can be invoked directly.

### Scripts that build reference data once, at installation

| Python source | Native procedure or current role |
| --- | --- |
| [`scripts/SpliceVault_tsv_vcf_conversion.py`](../scripts/SpliceVault_tsv_vcf_conversion.py) | Optional `priva splicevault-convert`. Installer calls `priva splicevault-convert`; CrossMap remains required. Installed BGZF route; general archive-byte/ZIP64 behavior is outside its qualified profile. |
| [`scripts/alphamissense_annotation_base.py`](../scripts/alphamissense_annotation_base.py) | Optional `priva alpha-annotation-base`: profile/header checks, base pair publication/recovery and completed motif adoption. `motif-runtime` explicitly observes the original Python package environment; motif adoption uses the native BioMart identity/receipt API. No annotation/model generation or installer caller switch. |
| [`scripts/alphamissense_cache.py`](../scripts/alphamissense_cache.py) | Optional `priva alphamissense-cache {install,validate,check-annotation}` with the script's options; the subcommand comes first. The pinned publisher size and MD5, the bounded schema probe, the install lock, adoption of a supplied or symlinked source, the resumed `curl` download, retained replaced artifacts, the rebuilt `tabix` index with its first-record probe and the ready receipt follow the script; `bgzip -t`, `tabix` and `curl` remain external tools. `check-annotation` is the native header check `alpha-annotation-base` already uses. Installer calls `priva alphamissense-cache install` in `AlphaMissense_install` and `check-annotation` in `AlphaMissense_anno`. Compared by `rust/tests/alphamissense_cache_compare.py` (success paths on the real hg19 publisher file; no case reaches the network). |
| [`scripts/alphamissense_domain_cache.py`](../scripts/alphamissense_domain_cache.py) | Replaced in the installer: `AlphaMissense_stat` calls `priva alphamissense-domain chromosome` for every VCF header contig (completed chromosomes are kept for a resumed run), then `merge`, and reuses a complete Parquet through `check --vcf --builder-sha256` with its recorded identity. Retained for manual use. |
| [`scripts/alphamissense_original_cache.py`](../scripts/alphamissense_original_cache.py) | Legacy/manual helper, outside native conversion scope. The installer no longer reads its receipts; completed Python statistics stay in use only through a valid Python DAS receipt. |
| [`scripts/alphmis_tsv2vcf.py`](../scripts/alphmis_tsv2vcf.py) | Optional `priva alphamissense-vcf`. Installer caller remains Python. |
| [`scripts/am_pick_intolerant_domains.py`](../scripts/am_pick_intolerant_domains.py) | Optional `priva das-domains lookup/transcripts/benign/coldspot-table/coldspot/scores/common/analyze` component procedures. Transcript projection accepts the explicit native Alpha exon-history companion. `scores` and `common` project the final typed statistics into the ordered domain score arrays and common-missense summaries (`das_projection_compare.py`). `analyze` runs the whole step in `analyze_domain_data` order and writes `domain_tolerance_analysis`, `domain_filter_results`, the native `all_intolerant_domains` set, and the curated-reference and excluded sets as JSON; KS rows follow query order rather than worker completion order (`das_assembly_compare.py`, four generated cases against the sorted-rescue Python reference). Installer calls `priva das-domains transcripts` and `priva das-domains analyze` in `AlphaMissense_pick_intolerant_domains` with `--python-hash-seed` from `priva_python_hash_seed`; `--metadata/--spots` and `--clinvar-vcf` are passed when the Python step would have received `--am_vcf` and `--clinvar_vcf`. The transcript map and `all_intolerant_domains` are written as direct `.rkyv` archives; completed Python DAS outputs with a valid Python receipt are kept for one transition. See the rescue-order and exon-history profiles below. |
| [`scripts/am_pick_missense_constrained_segments.py`](../scripts/am_pick_missense_constrained_segments.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/audit_clinvar_gofcards_review_tiers.py`](../scripts/audit_clinvar_gofcards_review_tiers.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/audit_condition_scope_inventory.py`](../scripts/audit_condition_scope_inventory.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/audit_gofcards_source_coverage.py`](../scripts/audit_gofcards_source_coverage.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/biomart_ftp_nm_to_enst.py`](../scripts/biomart_ftp_nm_to_enst.py) | Optional `priva biomart` FTP specification, metadata/XML and prepared-table identity/verification procedures. FTP acquisition and SQLite projection remain Python. |
| [`scripts/biomart_nm_to_enst.py`](../scripts/biomart_nm_to_enst.py) | Optional `priva biomart` identity, receipt, pinned-table verification and artifact routing with explicit `--helper-root`. Acquisition and installer caller migration remain pending. |
| [`scripts/build_gene_nonlof_mechanism_cache.py`](../scripts/build_gene_nonlof_mechanism_cache.py) | Retained Python; native replacement pending. |
| [`scripts/build_gene_pathogenic_mechanism_cache.py`](../scripts/build_gene_pathogenic_mechanism_cache.py) | Optional `priva gene-mechanism parse/build/attach` for prepared inputs. Acquisition, refresh, locking and installer orchestration remain Python; see the manual for the supported procedures. |
| [`scripts/build_gofcards_exact_gof_cache.py`](../scripts/build_gofcards_exact_gof_cache.py) | Retained Python; native replacement pending. |
| [`scripts/build_hpo_condition_mechanism_cache.py`](../scripts/build_hpo_condition_mechanism_cache.py) | Retained Python; native replacement pending. |
| [`scripts/build_mondo_disease_scope.py`](../scripts/build_mondo_disease_scope.py) | Optional `priva mondo-scope`. Installer calls `priva mondo-scope build`. |
| [`scripts/clinvar_annotation_contract.py`](../scripts/clinvar_annotation_contract.py) | Retained Python; native replacement pending. |
| [`scripts/clinvar_gene_refresh.py`](../scripts/clinvar_gene_refresh.py) | Retained Python; native replacement pending. |
| [`scripts/clinvar_stat_variants.py`](../scripts/clinvar_stat_variants.py) | Optional `priva clinvar-gene-stats build/merge` with an explicit Python hash seed and recorded completion order for reference comparisons. Installer calls `priva clinvar-gene-stats build` in `ClinVar_Gene_stat_regenerate` with `--python-hash-seed` from `priva_python_hash_seed` (default 0) and writes the refresh candidate `clinvar_2star_stats.rkyv`; the bundled pickle is still converted by `cache-install publish`. |
| [`scripts/collapse_HPO_anno.py`](../scripts/collapse_HPO_anno.py) | Optional `priva hpo-assertions build`. Installer calls `priva hpo-assertions build`. |
| [`scripts/domain_metadata_cache.py`](../scripts/domain_metadata_cache.py) | Optional `priva domain-metadata`. Installer calls `priva domain-metadata build --assembly` (reused through its `.source.json`) for the DAS coldspot filter and writes `domain_metadata.{asm}.parquet`. |
| [`scripts/inject_clinvar_into_gofcards.py`](../scripts/inject_clinvar_into_gofcards.py) | Retained Python; native replacement pending. |
| [`scripts/prepare_clingen_map.py`](../scripts/prepare_clingen_map.py) | Optional `priva clingen-map`. Installer calls `priva clingen-map build` and writes `cache_root/ClinGen/clingen_map.{asm}.rkyv`; bundled maps are reused unchanged. |
| [`scripts/prepare_pm1_regions.py`](../scripts/prepare_pm1_regions.py) | Optional `priva pm1-regions`. Installer calls `priva mapped-artifact run --native-producer pm1-regions` (`pm1-regions build`) in `prepare_pm1_regions` and writes `pm1_regions.{asm}.rkyv` (the `.rkyv` sibling of a configured pickle value); `--intolerant-domains` is `priva cache-install resolve --key all_intolerant_domains`, so a Python DAS pickle or a native DAS archive is read. Freshness follows the priva binary. Consumes prepared DAS; does not regenerate it. |
| [`scripts/resort_step3_only.py`](../scripts/resort_step3_only.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/stat_aachange_clinvar.py`](../scripts/stat_aachange_clinvar.py) | Optional `priva clinvar-aa-splice`. Installer calls `priva clinvar-aa-splice build` and writes the `.rkyv` siblings of configured pickle values; the pair receipts hash `priva builder-identity clinvar-aa-splice` as the producer. Completed Python pickle pairs stay reusable for one transition. |
| [`scripts/stat_exon_patho_afs.py`](../scripts/stat_exon_patho_afs.py) | Optional `priva clinvar-exon-af`. Installer calls `priva clinvar-exon-af build` and writes the `.rkyv` sibling of a configured pickle value. |
| [`scripts/stat_gene_patho_afs.py`](../scripts/stat_gene_patho_afs.py) | Optional `priva clinvar-gene-af`. Installer calls `priva clinvar-gene-af build` and writes the `.rkyv` sibling of a configured pickle value; an absent archive (no CSQ format) fails the stage. |
| [`scripts/stat_protein_domain_amscores.py`](../scripts/stat_protein_domain_amscores.py) | Optional `priva alphamissense-domain`. Installer calls `priva alphamissense-domain chromosome/merge/check/identity` in `AlphaMissense_stat` and writes the `.parquet` sibling of a configured pickle value. `check` prints the recorded and current `builder_source_sha256`; `--builder-sha256 SHA256` accepts an earlier recorded builder, other checks unchanged. `identity --input PARQUET` prints the recorded value for `das-domains --source-builder-sha256`. |

## Installer helpers

### Installer and cluster deployment helpers

| Python source | Native procedure or current role |
| --- | --- |
| [`scripts/cache_resources.py`](../scripts/cache_resources.py) | Optional `priva cache-resources list/show/refresh-plan/validate` with an explicit registry and Python hash seed. Existing callers remain Python; unusual JSON/set boundaries are documented in the user manual. |
| [`scripts/check_cluster_environment.py`](../scripts/check_cluster_environment.py) | Retained Python environment inspection. |
| [`scripts/check_cluster_readiness.py`](../scripts/check_cluster_readiness.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/clinvar_annotation_contract.py`](../scripts/clinvar_annotation_contract.py) | Optional `priva clinvar-contract` header, native cache-shape, AF path and AA/splice completion-association procedures. Original source assets remain explicit. Installer calls `priva clinvar-contract header`, `cache` and `af-paths-check` in the AA/splice, AF, gene-statistics regeneration and PM5 stages; ClinVar VCF deployment, the legacy pickle checks and the AA/splice pair receipts remain Python. |
| [`scripts/conservation_resource_check.py`](../scripts/conservation_resource_check.py) | Optional `priva conservation-resource`. Installer calls `priva conservation-resource validate/publish`. |
| [`scripts/das_input_continuity.py`](../scripts/das_input_continuity.py) | Optional `priva das-input-continuity check/record`; verifies prepared ordered-consumer equivalence proofs and preserves original DAS outputs. Requires the pinned checker source asset, never a Python process. Installer caller remains Python. |
| [`scripts/das_install_state.py`](../scripts/das_install_state.py) | Optional `priva das-install-state capture/check/record`. Preserves completion receipt v3, source-role hashes, metadata identity and bundled existence/nonzero checks. Installer calls `priva das-install-state` for native DAS builds, with `priva builder-identity` output as the source roles; the legacy Python DAS receipt check remains Python. |
| [`scripts/deploy_cluster_gnomad.py`](../scripts/deploy_cluster_gnomad.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/deploy_cluster_plugin_inputs.py`](../scripts/deploy_cluster_plugin_inputs.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/deploy_cluster_resources.py`](../scripts/deploy_cluster_resources.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/gnomad_cache.py`](../scripts/gnomad_cache.py) | Optional `priva gnomad-cache {members,validate,default-path}` with the script's options; the subcommand comes first. Uses the family naming, header, assembly and first-indexed-record checks that `installer-config-paths` adoption already reproduces, read by the native BGZF and tabix/CSI readers instead of htslib. Installer calls `priva gnomad-cache default-path/validate` in `gnomAD_install`. Compared by `rust/tests/gnomad_cache_compare.py` on generated families (no gnomAD data is read). |
| [`scripts/indexed_plugin_cache.py`](../scripts/indexed_plugin_cache.py) | Optional `priva indexed-plugin-cache ensure-index/format-primate`. Retains bgzip, tabix and GNU sort; installer calls `priva indexed-plugin-cache`. |
| [`scripts/installed_environment.py`](../scripts/installed_environment.py) | Optional `priva installed-environment export/export-installation` for package-only YAML exports. Retains the installed package manager and statically linked libyaml; installer calls `priva installed-environment export-installation`. |
| [`scripts/installer_config_paths.py`](../scripts/installer_config_paths.py) | Optional `priva installer-config-paths {resolve,managed-roots} --config CONFIG --source-root CHECKOUT` with `--check`, `--output` and `--managed-updates-json`, plus the library procedures (`portable-updates`, `resolve-cache-paths`, `defaults`, `resource-groups`, `adoption-plan`, `adopt-sources`, `render` and the path, identity and gnomAD checks) for comparison. Planning, defaults, byte-preserving adoption with companions, retained earlier files, the adoption lock and receipt, and the atomic configuration rewrite follow the script; configuration YAML is read with the script's duplicate-key refusal, and escaped surrogate code points survive. It needs an explicit `PYTHONHASHSEED` (the script's replacement-table set order) and refuses a checkout whose `spliceai_cache_input.py` or `gnomad_cache.py` differs from the built copy. `CADD_install` calls `managed-roots` with `PYTHONHASHSEED=0`; the configuration steps of `priva_prepare_install_config` run before `PriVA_rust_install` records `priva_binary` and remain Python. Compared by `rust/tests/installer_config_paths_compare.py`. |
| [`scripts/installer_prerequisites.py`](../scripts/installer_prerequisites.py) | Optional `priva installer-prerequisites CONFIG --source-root CHECKOUT` with every `--stage` and option; `--source-root` is required. Reads YAML through `yaml_config` (refuses anchors and other unsupported constructs) and reads only plain-data ClinVar pickles; `--allow-pinned-reference-defaults` plans and checks the references with the native `reference-resources` procedures. `--stage manual-inputs --prompt` at a terminal collects unset paths and rewrites the configuration through the native `render_config` of `installer_config_paths.py`. Installer calls `priva installer-prerequisites` for the mapping check of the AlphaMissense motif and PM1 producers, `ClinVar_Gene_stat` and the GoFCards build checks (including the build route's pinned-reference preflight). The main, manual-input, snapshot, NM-to-ENST installation and GoFCards main-preflight calls run before `PriVA_rust_install` records `priva_binary` and remain Python. Compared by `rust/tests/installer_prerequisites_compare.py`. |
| [`scripts/mapped_artifact_publish.py`](../scripts/mapped_artifact_publish.py) | Optional `priva mapped-artifact run/recover/adopt-native/context`: native BioMart checks, literal external producer arguments, pending pair publication and preserved-artifact adoption. Scientific producers remain external. Installer calls `priva mapped-artifact recover/run --native-producer pm1-regions` in `prepare_pm1_regions`; the AlphaMissense motif caller remains Python. `--native-producer pm1-regions` instead runs this executable's `pm1-regions build` (exact `--name value` pairs bound to the selected assembly, mapping and HCSeeker/RMC anchors) and publishes `pm1_regions.{asm}.rkyv`; the context and receipt record `native_producer`. `biomart artifact-route/changed-target --artifact-format rkyv` names that archive. |
| [`scripts/mechanism_source_input.py`](../scripts/mechanism_source_input.py) | Optional `priva mechanism-source validate/record-source/check-source` for acquired-source validation and receipts. Retains distinct acquisition/runtime checks and bundled reuse policy. Existing callers remain Python. |
| [`scripts/prepare_cluster_config.py`](../scripts/prepare_cluster_config.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/reference_resources.py`](../scripts/reference_resources.py) | Optional `priva reference-resources {sources,preflight,install} [CONFIG] --data-dir DATA` with `--scope`, `--use-cache-root` and `--write-config`; `--data-dir` is required for preflight and install. Planning (managed defaults, legacy-default redirection, cache_root routing, the inspection-root policy), FASTA, index and forward-chain checks, pinned-source receipts, per-resource locks, link-based publication, `samtools faidx` index repair and the configuration rewrite follow the script; YAML is read through `yaml_config` (refuses anchors and flow/JSON documents). The download of a missing pinned UCSC archive runs `curl` with the same URL, Range, status and Content-Range rules; transport failures carry curl's message. Installer calls `priva reference-resources install` in `reference_genome_install`, `assembly_references_install` and `gofcards_reference_install`; the selected-reference preflight in `priva_install_prerequisites` runs before `PriVA_rust_install` records `priva_binary` and remains Python. Compared by `rust/tests/reference_resources_compare.py` (no case reaches the network). |
| [`scripts/spliceai2_to_vep_vcf.py`](../scripts/spliceai2_to_vep_vcf.py) | Optional `priva spliceai2-prepare {prepare,gene-context}` with the script's options; the subcommand comes first. Manifest, input and tool verification, the GENCODE48 gene context with the VEP release's gene names, the gzip, projection and bgzip pipeline, the tabix index with its first-record probe, the receipt and the publication order follow the script. By decision (4 Nov 2026) the projection runs as this executable's `spliceai2-vcf` instead of a converter compiled from `spliceai2_converter.rs` with rustc at installation: `build-converter` is not provided, and the receipt's `converter` (`role` `priva_builtin`) and `driver_sha256` name this executable. `--converter-bin` still runs an explicit standalone converter. Installer calls `priva spliceai2-prepare prepare` in `SpliceAI2_install`. Rust-only `--jobs N` checks and converts the shards N at a time and joins the BGZF parts; `rust/tests/spliceai2_prepare_jobs_compare.py` checks it against one job on real shards. Compared by `rust/tests/spliceai2_prepare_compare.py`. |
| [`scripts/vep_chromosome_cache.py`](../scripts/vep_chromosome_cache.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/vep_chromosome_jobs.py`](../scripts/vep_chromosome_jobs.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/vep_plugins_release.py`](../scripts/vep_plugins_release.py) | Optional `priva vep-plugins {check,install,transition,clear-transition,ensure-utr-reporting} --release RELEASE` with the same options; the mode comes first. Receipt checks, local-patch recognition, staged publication with the release-switch and partial-retry rules, the transition record, and the locked UTRAnnotator 116 reporting transaction (pinned upstream and installed bytes, original and patch backups, the pending recovery receipt) follow the script; the patch is the same `difflib.unified_diff` text. The bundled reporting module is read from `PRIVA_DATA_DIR`. Installer calls `priva vep-plugins check/install/ensure-utr-reporting/clear-transition` in `vep_cache_api_install`; the `transition` call in `conda_install_vep` runs before `PriVA_rust_install` records `priva_binary` and remains Python. Compared by `rust/tests/vep_plugins_compare.py` (no case reaches the network). |

## Task adapters and retained bridges

Production Nextflow integration is incomplete. Existing adapters and bootstrap
bridges remain until their caller replacements are qualified.

### Nextflow task adapters

| Python source | Native procedure or current role |
| --- | --- |
| [`bin/priva_task.py`](../bin/priva_task.py) | Retained Python; native replacement pending. |
| [`bin/priva_task_boundary.py`](../bin/priva_task_boundary.py) | Retained Python; native replacement pending. |

### Bootstrap and compatibility bridges

| Python source | Native procedure or current role |
| --- | --- |
| [`scripts/bootstrap_priva_rust.py`](../scripts/bootstrap_priva_rust.py) | Retained bootstrap/compatibility bridge; remove only after all callers migrate. |
| [`scripts/priva_rust_control_multiply.py`](../scripts/priva_rust_control_multiply.py) | Retained bootstrap/compatibility bridge; remove only after all callers migrate. |
| [`scripts/priva_rust_criteria_summary.py`](../scripts/priva_rust_criteria_summary.py) | Retained bootstrap/compatibility bridge; remove only after all callers migrate. |
| [`scripts/priva_rust_gene_weights.py`](../scripts/priva_rust_gene_weights.py) | Retained bootstrap/compatibility bridge; remove only after all callers migrate. |
| [`scripts/priva_rust_pedigree_context.py`](../scripts/priva_rust_pedigree_context.py) | Retained bootstrap/compatibility bridge; remove only after all callers migrate. |
| [`scripts/priva_rust_pext.py`](../scripts/priva_rust_pext.py) | Retained bootstrap/compatibility bridge; remove only after all callers migrate. |
| [`scripts/priva_rust_posterior.py`](../scripts/priva_rust_posterior.py) | Retained bootstrap/compatibility bridge; remove only after all callers migrate. |
| [`scripts/priva_rust_rank_setup.py`](../scripts/priva_rust_rank_setup.py) | Retained bootstrap/compatibility bridge; remove only after all callers migrate. |
| [`scripts/priva_rust_ranking_compatibility.py`](../scripts/priva_rust_ranking_compatibility.py) | Retained bootstrap/compatibility bridge; remove only after all callers migrate. |
| [`scripts/priva_rust_terminal_order.py`](../scripts/priva_rust_terminal_order.py) | Retained bootstrap/compatibility bridge; remove only after all callers migrate. |

## Retained Python tools

CADD scientific packages and model assets remain upstream. The auxiliary data
tools below are separate from the native runtime.

### Kept in Python: PriVA's replacements inside upstream CADD

| Python source | Native procedure or current role |
| --- | --- |
| [`scripts/cadd_prescore_index.py`](../scripts/cadd_prescore_index.py) | Retained Python. |
| [`scripts/cadd_scripts/MMSplice.py`](../scripts/cadd_scripts/MMSplice.py) | Retained Python. |
| [`scripts/cadd_scripts/extract_scored.py`](../scripts/cadd_scripts/extract_scored.py) | Retained Python. |
| [`scripts/cadd_scripts_install.py`](../scripts/cadd_scripts_install.py) | Retained Python. |

### Kept in Python: one-off research tools and examples

| Python source | Native procedure or current role |
| --- | --- |
| [`data/InterPro/gemini_interpro_func_classifier.py`](../data/InterPro/gemini_interpro_func_classifier.py) | Retained Python. |
| [`data/MaveDB/cal_mavedb_neutral_bands.py`](../data/MaveDB/cal_mavedb_neutral_bands.py) | Retained Python. |
| [`data/MaveDB/gpt_mavedb_assay_mapper.py`](../data/MaveDB/gpt_mavedb_assay_mapper.py) | Retained Python. |
| [`data/MaveDB/pretify_json.py`](../data/MaveDB/pretify_json.py) | Retained Python. |
| [`example/smoke_test.py`](../example/smoke_test.py) | Retained Python. |

## DAS rescue-order profile

The explicitly generated Python reference uses `domain-sorted-rescue-v1`: sort the
tolerant domain identifiers before selecting candidates in both rescue passes,
then retain the original `42 + i` / `99 + i` seed assignment, 32-bit wraparound,
permutation count, candidate filters, numerical routines and FDR thresholds.
This makes each domain's seed independent of Python string hashing and KS worker
completion order. It is a deliberate profile difference from the original
unsorted Python behavior. Other DAS stage ordering is unchanged.

Existing caches retain their producing profile and receipts. This change does
not rewrite them or qualify them as outputs of the new profile. A future native
DAS implementation must match the explicitly pinned sorted Python reference.
Generate that reference with `python rust/tests/das_sorted_reference.py render
--source scripts/am_pick_intolerant_domains.py --output /work/sorted_reference.py`.
The generated source is exclusively created; its whole AST differs only at the
two approved iterators. Keep its reported source/output hashes with the selected
profile. The still-active Python installer script remains byte-identical to the
original profile so its source-hash reuse gate does not trigger a new calculation.

## Native Alpha exon-history companion

```text
Annotated Alpha VCF -> per-chromosome first-insertion history
                    -> header-ordered merge with an explicit Python hash seed
                    -> serialized-set history -> DAS transcript mapping
```

`alpha-exon-history chromosome` captures the original observation order without
scores or a Python pickle. Merge all header contigs in their original order,
including empty contigs, with `alpha-exon-history merge --input FIRST --input
SECOND ... --hash-seed UINT32 --output HISTORY.parquet`. The final companion
reproduces CPython 3.12 set construction, worker transport and ordered updates
for that explicit seed. `build --vcf VCF --hash-seed UINT32 --output PARQUET`
uses the original sequential route instead.

Export with `alpha-exon-history export --input HISTORY.parquet --output HISTORY`;
pass that file to `das-domains transcripts --exon-history HISTORY` alongside the
typed statistics and the remaining required options. This retains exon traversal
order that the statistics model's set membership alone cannot recover. It does
not recover an unknown historical hash seed or require the legacy score pickle.
Source/profile guards remain active: do not retag an older companion when a
guard rejects it. The installer (`AlphaMissense_pick_intolerant_domains`) builds
`<statistics>.exon_history.seed<S>.parquet` this way with the
`priva_python_hash_seed` value, keeps completed chromosomes for a resumed run,
and exports the `.lines` file beside it.

Native statistics and history files keep the builder identity recorded when
they were built; rebuild them only on explicit request. `alphamissense-domain
check` and `alpha-exon-history check` print the recorded and current
`builder_source_sha256` and, by default, require the current builder. With
`--builder-sha256 SHA256` they accept exactly that recorded builder; VCF,
producer, extractor and profile checks stay strict. Merge and export still
require the current builder. `identity --input PARQUET` in either family prints
the recorded value without qualifying the file, for `das-domains transcripts`
`--source-builder-sha256` and `--history-builder-sha256`. Installer receipts use
`priva builder-identity NAME`: a SHA-256 over that builder's own embedded source
modules and Cargo.lock, listed by `priva builder-identity --help`. Its
`alphamissense-domain`, `domain-metadata` and `alpha-exon-history` values equal
the identities those builders record in their footers.
