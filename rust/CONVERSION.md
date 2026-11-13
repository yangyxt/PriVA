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
| [`scripts/clinvar_vcv.py`](../scripts/clinvar_vcv.py) | Native `acmg_step::clinvar_vcv` text input; the installer's XML parsing is `priva gofcards-clinvar-injection`. |
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
| [`scripts/alphamissense_annotation_base.py`](../scripts/alphamissense_annotation_base.py) | Optional `priva alpha-annotation-base`: profile/header checks, base pair publication/recovery and completed motif adoption. `motif-runtime` explicitly observes the original Python package environment; motif adoption uses the native BioMart identity/receipt API. No annotation/model generation. `AlphaMissense_anno` calls it for path, recover, check, runtime-associated, restore, capture and publish (November 2026); on the G9 installation it captured the same record as the Python helper. |
| [`scripts/alphamissense_cache.py`](../scripts/alphamissense_cache.py) | Optional `priva alphamissense-cache {install,validate,check-annotation}` with the script's options; the subcommand comes first. The pinned publisher size and MD5, the bounded schema probe, the install lock, adoption of a supplied or symlinked source, the resumed `curl` download, retained replaced artifacts, the rebuilt `tabix` index with its first-record probe and the ready receipt follow the script; `bgzip -t`, `tabix` and `curl` remain external tools. `check-annotation` is the native header check `alpha-annotation-base` already uses. Installer calls `priva alphamissense-cache install` in `AlphaMissense_install` and `check-annotation` in `AlphaMissense_anno`. Compared by `rust/tests/alphamissense_cache_compare.py` (success paths on the real hg19 publisher file; no case reaches the network). |
| [`scripts/alphamissense_domain_cache.py`](../scripts/alphamissense_domain_cache.py) | Replaced in the installer: `AlphaMissense_stat` calls `priva alphamissense-domain chromosome` for every VCF header contig (completed chromosomes are kept for a resumed run), then `merge`, and reuses a complete Parquet through `check --vcf --builder-sha256` with its recorded identity. Retained for manual use. |
| [`scripts/alphamissense_original_cache.py`](../scripts/alphamissense_original_cache.py) | Legacy/manual helper, outside native conversion scope. The installer no longer reads its receipts; completed Python statistics stay in use only through a valid Python DAS receipt. |
| [`scripts/alphmis_tsv2vcf.py`](../scripts/alphmis_tsv2vcf.py) | `priva alphamissense-vcf build --input TSV --output VCF`, called by `AlphaMissense_anno` (November 2026). The annotation recipe fingerprint counts this line as the original `python "$convert_py"` line, so caches annotated before the switch keep their fingerprint (same recipe and converter fingerprints on the G9 installation). |
| [`scripts/am_pick_intolerant_domains.py`](../scripts/am_pick_intolerant_domains.py) | Optional `priva das-domains lookup/transcripts/benign/coldspot-table/coldspot/scores/common/analyze` component procedures. Transcript projection accepts the explicit native Alpha exon-history companion. `scores` and `common` project the final typed statistics into the ordered domain score arrays and common-missense summaries (`das_projection_compare.py`). `analyze` runs the whole step in `analyze_domain_data` order and writes `domain_tolerance_analysis`, `domain_filter_results`, the native `all_intolerant_domains` set, and the curated-reference and excluded sets as JSON; KS rows follow query order rather than worker completion order (`das_assembly_compare.py`, four generated cases against the sorted-rescue Python reference). Installer calls `priva das-domains transcripts` and `priva das-domains analyze` in `AlphaMissense_pick_intolerant_domains` with `--python-hash-seed` from `priva_python_hash_seed`; `--metadata/--spots` and `--clinvar-vcf` are passed when the Python step would have received `--am_vcf` and `--clinvar_vcf`. The transcript map and `all_intolerant_domains` are written as direct `.rkyv` archives; completed Python DAS outputs with a valid Python receipt are kept for one transition. See the rescue-order and exon-history profiles below. |
| [`scripts/am_pick_missense_constrained_segments.py`](../scripts/am_pick_missense_constrained_segments.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/audit_clinvar_gofcards_review_tiers.py`](../scripts/audit_clinvar_gofcards_review_tiers.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/audit_condition_scope_inventory.py`](../scripts/audit_condition_scope_inventory.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/audit_gofcards_source_coverage.py`](../scripts/audit_gofcards_source_coverage.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/biomart_ftp_nm_to_enst.py`](../scripts/biomart_ftp_nm_to_enst.py) | Optional `priva biomart` FTP specification, metadata/XML, projection, identity/verification, `ftp-make-provenance`, `ftp-portable-provenance`, and the FTP path of `atomic-json`, `runtime-receipt`, `external-path` and `fetch`. `fetch --source ensembl-grch38-mart116-ftp` is `install`: the same directory flock, pinned dumps, resume by Range through `curl`, projection and receipts. Installer caller remains Python; see the row below. |
| [`scripts/biomart_nm_to_enst.py`](../scripts/biomart_nm_to_enst.py) | Optional `priva biomart` identity, receipt, pinned-table verification, artifact routing, receipt publication and `fetch` with explicit `--helper-root`. `fetch` takes the script's options and sources, requests through `curl` and records `priva biomart fetch` with the `biomart-fetch` builder identity as receipt producer. Installer calls `priva biomart identity` and `artifact-route` in the AlphaMissense motif and PM1 steps; `nm_to_enst_mapping_install` stays on Python `fetch` and `verify` because `main_install` reaches it before `priva_binary` is recorded. Comparison: `rust/tests/biomart_fetch_compare.py` through a local HTTPS publisher. Python's `RecursionError` on JSON nested about 994 deep has no native counterpart, by design. |
| [`scripts/build_gene_nonlof_mechanism_cache.py`](../scripts/build_gene_nonlof_mechanism_cache.py) | Optional `priva gene-nonlof-mechanism --source-root CHECKOUT` with the script's options. Source freshness and receipts, the selected-URL checks, the G2P, ClinGen dosage, PanelApp and GoFCards readers (the workbook read as openpyxl and pandas read it), the HGNC mapping, the canonical JSON with its Draft 2020-12 schema check, the build lock, the weekly ClinVar VCV acquisition (curl) and the explicit VCV XML check (Expat) follow the script; `build_inputs` keeps the SHA-256 of the checkout's Python builder and `clinvar_vcv.py`, so either one reuses the other's cache. urllib transfers run through `curl --http1.1`; ISO date workbook cells and schema patterns beyond anchors, classes and simple quantifiers are refused. Installer calls `priva gene-nonlof-mechanism` when `gene_nonlof_mechanism_builder_script` is the checkout's script; another configured script still runs with Python. `rust/tests/gene_nonlof_mechanism_compare.py` compares them. |
| [`scripts/build_gene_pathogenic_mechanism_cache.py`](../scripts/build_gene_pathogenic_mechanism_cache.py) | Optional `priva gene-mechanism parse/build/attach` for prepared inputs. Acquisition, refresh, locking and installer orchestration remain Python; see the manual for the supported procedures. |
| [`scripts/build_gofcards_exact_gof_cache.py`](../scripts/build_gofcards_exact_gof_cache.py) | Optional `priva gofcards-exact-gof` with the script's options. The workbook (pandas' default column typing), the summary cache and its retries, the GRCh37 reference check, the bcftools, VEP and CrossMap runs, the endpoint recovery, the HGNC and review rules, the export, the rejects table and the stats file follow the script; the gzip is the same zlib stream, with `gzip.open`'s name and time in the header. Requests run through `curl`, so transport failures carry curl's message; `--summary-endpoint` replaces the script's constant for tests; a workbook with two columns mapping to one source field and a compressed FASTA are refused. Installer calls it in `gofcards_exact_gof_cache_install` when `gofcards_builder_script` is the checkout's own script; a custom script still runs with Python. Compared by `rust/tests/gofcards_exact_gof_compare.py` (no case reaches a publisher). |
| [`scripts/build_hpo_condition_mechanism_cache.py`](../scripts/build_hpo_condition_mechanism_cache.py) | Optional `priva hpo-condition-mechanism --source-root CHECKOUT` with the script's options. The HPO frame and its MONDO merge, the G2P/Orphadata mechanisms and penetrance, ClinGen haploinsufficiency, the GoFCards exact alleles with their ClinVar condition links, the deduced recessive LOF, the scope gate, the payload and provenance checks and the atomic publication follow the script; the default ClinGen table comes from the checkout, and the cache records the checkout's Python builder as its builder source, so either one accepts the other's cache. Input tables are decoded completely before their rows are read. Installer calls `priva hpo-condition-mechanism` when `hpo_condition_mechanism_builder_script` is the checkout's script; another configured script still runs with Python. `rust/tests/hpo_condition_mechanism_compare.py` compares them. |
| [`scripts/build_mondo_disease_scope.py`](../scripts/build_mondo_disease_scope.py) | Optional `priva mondo-scope`. Installer calls `priva mondo-scope build`. |
| [`scripts/clinvar_annotation_contract.py`](../scripts/clinvar_annotation_contract.py) | Replaced in the installer by `priva clinvar-contract`; see the installer helpers table. |
| [`scripts/clinvar_gene_refresh.py`](../scripts/clinvar_gene_refresh.py) | Optional `priva clinvar-gene-refresh {plan,install,deploy,retire-owned,retire-shards,check-header,acquire,file-date,candidate-route,candidate-begin,candidate-complete}` with the script's options; the subcommand comes first. The registry policy, the bounded fileDate checks, the NCBI MD5 and 1 MiB header probes (`curl`), the scientific profile and its signature, completed-receipt reuse, the locked deployment with retained generations and recovery journals, the candidate receipts and the configuration rewrite follow the script; the original ClinVar stages still run through `install_utils.sh`. The profile observes the Python runtime with one `python -c` so that signatures and workspaces stay those of the script and `vep_chromosome_jobs.py`. Installer calls `priva clinvar-gene-refresh` in `ClinVar_install`, `ClinVar_VCF_deploy` and `_ClinVar_VCF_deploy`; `annotation_vcf.sh` calls it for `retire-shards` and the candidate phases. Retained for `vep_chromosome_jobs.py` and the Python tests. `rust/tests/clinvar_gene_refresh_compare.py` compares it with the script. |
| [`scripts/clinvar_stat_variants.py`](../scripts/clinvar_stat_variants.py) | Optional `priva clinvar-gene-stats build/merge` with an explicit Python hash seed and recorded completion order for reference comparisons. Installer calls `priva clinvar-gene-stats build` in `ClinVar_Gene_stat_regenerate` with `--python-hash-seed` from `priva_python_hash_seed` (default 0) and writes the refresh candidate `clinvar_2star_stats.rkyv`; the bundled pickle is still converted by `cache-install publish`. |
| [`scripts/collapse_HPO_anno.py`](../scripts/collapse_HPO_anno.py) | Optional `priva hpo-assertions build`. Installer calls `priva hpo-assertions build`. |
| [`scripts/domain_metadata_cache.py`](../scripts/domain_metadata_cache.py) | Optional `priva domain-metadata`. Installer calls `priva domain-metadata build --assembly` (reused through its `.source.json`) for the DAS coldspot filter and writes `domain_metadata.{asm}.parquet`. |
| [`scripts/inject_clinvar_into_gofcards.py`](../scripts/inject_clinvar_into_gofcards.py) | Optional `priva gofcards-clinvar-injection --source-root CHECKOUT` with the script's options. The streamed VCV reader and matching rules of `clinvar_vcv.py` (the same Expat feed and clear timing as `iterparse`), the injection, the atomic gzip output (Python's bytes), the audit table, the statistics and the log lines follow the script. The cache records the checkout's Python injector as `clinvar.injector_sha256`, so either one accepts the other's cache. Installer calls `priva gofcards-clinvar-injection` when `gofcards_clinvar_injector_script` is the checkout's script; another configured script still runs with Python. Retained for the Python tests. `rust/tests/gofcards_clinvar_injection_compare.py` compares it with the script. |
| [`scripts/validate_gofcards_cache.py`](../scripts/validate_gofcards_cache.py) | Optional `priva gofcards-cache-validate CACHE [KEY=VALUE ...]`. Installer calls it in `validate_gofcards_exact_gof_cache`; the main preflight, which runs before `priva_binary` is recorded, stays on the script. `rust/tests/gofcards_clinvar_injection_compare.py` compares them. |
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
| [`scripts/clinvar_annotation_contract.py`](../scripts/clinvar_annotation_contract.py) | Optional `priva clinvar-contract` header, native cache-shape, AF path and AA/splice completion-association procedures. Original source assets remain explicit. `legacy-cache` applies the same shape checks to the Python pickles of an existing installation through the restricted pickle reader; a pickle it cannot read fails the check, so the stage rebuilds. Installer calls `priva clinvar-contract header`, `cache`, `legacy-cache`, `af-paths-check` and `pair-capture`/`pair-check`/`pair-complete` (with `--contract-source` naming the Python contract, whose SHA-256 the pair receipts keep) in the ClinVar VCF deployment, AA/splice, AF, gene-statistics regeneration and PM5 stages. The Python module stays for `clinvar_stat_variants.py` and the tests. |
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
| [`scripts/mechanism_source_input.py`](../scripts/mechanism_source_input.py) | Optional `priva mechanism-source validate/record-source/check-source` for acquired-source validation and receipts. Retains distinct acquisition/runtime checks and bundled reuse policy. `ClinGen_deploy` and `priva_clingen_download` call it (November 2026). |
| [`scripts/prepare_cluster_config.py`](../scripts/prepare_cluster_config.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/reference_resources.py`](../scripts/reference_resources.py) | Optional `priva reference-resources {sources,preflight,install} [CONFIG] --data-dir DATA` with `--scope`, `--use-cache-root` and `--write-config`; `--data-dir` is required for preflight and install. Planning (managed defaults, legacy-default redirection, cache_root routing, the inspection-root policy), FASTA, index and forward-chain checks, pinned-source receipts, per-resource locks, link-based publication, `samtools faidx` index repair and the configuration rewrite follow the script; YAML is read through `yaml_config` (refuses anchors and flow/JSON documents). The download of a missing pinned UCSC archive runs `curl` with the same URL, Range, status and Content-Range rules; transport failures carry curl's message. Installer calls `priva reference-resources install` in `reference_genome_install`, `assembly_references_install` and `gofcards_reference_install`; the selected-reference preflight in `priva_install_prerequisites` runs before `PriVA_rust_install` records `priva_binary` and remains Python. Compared by `rust/tests/reference_resources_compare.py` (no case reaches the network). |
| [`scripts/spliceai2_to_vep_vcf.py`](../scripts/spliceai2_to_vep_vcf.py) | Optional `priva spliceai2-prepare {prepare,gene-context}` with the script's options; the subcommand comes first. Manifest, input and tool verification, the GENCODE48 gene context with the VEP release's gene names, the gzip, projection and bgzip pipeline, the tabix index with its first-record probe, the receipt and the publication order follow the script. By decision (4 Nov 2026) the projection runs as this executable's `spliceai2-vcf` instead of a converter compiled from `spliceai2_converter.rs` with rustc at installation: `build-converter` is not provided, and the receipt's `converter` (`role` `priva_builtin`) and `driver_sha256` name this executable. `--converter-bin` still runs an explicit standalone converter. Installer calls `priva spliceai2-prepare prepare` in `SpliceAI2_install`. Rust-only `--jobs N` checks and converts the shards N at a time and joins the BGZF parts; `rust/tests/spliceai2_prepare_jobs_compare.py` checks it against one job on real shards. Compared by `rust/tests/spliceai2_prepare_compare.py`. |
| [`scripts/vep_chromosome_cache.py`](../scripts/vep_chromosome_cache.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/vep_chromosome_jobs.py`](../scripts/vep_chromosome_jobs.py) | Legacy/manual helper, outside native conversion scope. |
| [`scripts/vep_plugins_release.py`](../scripts/vep_plugins_release.py) | Optional `priva vep-plugins {check,install,transition,clear-transition,ensure-utr-reporting} --release RELEASE` with the same options; the mode comes first. Receipt checks, local-patch recognition, staged publication with the release-switch and partial-retry rules, the transition record, and the locked UTRAnnotator 116 reporting transaction (pinned upstream and installed bytes, original and patch backups, the pending recovery receipt) follow the script; the patch is the same `difflib.unified_diff` text. The bundled reporting module is read from `PRIVA_DATA_DIR`. Installer calls `priva vep-plugins check/install/ensure-utr-reporting/clear-transition` in `vep_cache_api_install`; the `transition` call in `conda_install_vep` runs before `PriVA_rust_install` records `priva_binary` and remains Python. Compared by `rust/tests/vep_plugins_compare.py` (no case reaches the network). |

## Task adapters and bridges

The Nextflow task adapters are converted, and the Python-to-Rust bridges are
removed (G10, November 2026).

| Former Python source | Native procedure or current role |
| --- | --- |
| `bin/priva_task.py` | `priva task {manifest,validate,launch,run}` ([`rust/src/nextflow_task/`](src/nextflow_task/mod.rs)). Same checks, refusal messages, JSON records and exit codes; the per-task configuration is written by editing the original text instead of re-dumping it, with PyYAML-only spellings (such as `yes`) rewritten so `yq` reads the values PyYAML loaded. Compared with the Python adapter by `rust/tests/nextflow_task_compare.py` (38 of 38 synthetic cases, check job 1538049; 37 cases since the raw refusal case was dropped with the raw gate); its `--source` must be a revision that still has the Python adapters, such as `918c87f`. |
| `bin/priva_task_boundary.py` | `priva task-boundary {check,exec}`, same Bubblewrap arguments and exit code 125. |
| `scripts/priva_rust_*.py` (nine ranking and criteria bridges) | Removed with their tests. They were the optional `--ranking_binary` path inside the Python ACMG step; the workflow runs `priva acmg-assign`, which accepts and ignores `--ranking_binary`. `scripts/acmg_scoring.py` and `scripts/acmg_criteria_assign.py` match `main` again. |

[`scripts/bootstrap_priva_rust.py`](../scripts/bootstrap_priva_rust.py) stays: the
installer runs it to build `priva` itself, so it cannot be replaced by `priva`.

## Raw mode and the task layout

Raw mode was refused by `priva task` until its annotation phase (VEP, CADD
and the other annotation steps) was compared with the shell workflow. The
refusal was removed after the BH14233 chr19 1–3.5 Mb trio passed (job 1538402
at `a50f82a`, 8 November 2026):

| Run | Compared with the shell workflow |
| --- | --- |
| Local executor | Annotated VCF records and header (dated lines aside), annotation CADD table, filtered VCF records, filtered CADD table, final table and ACMG matrix: all identical |
| Same launch with `-resume` | All four tasks reused |
| Slurm profile, each task its own Slurm job on three nodes | All identical, with 12 annotation CPUs |

Fixes found on the way:

| Problem | Fix |
| --- | --- |
| `--slurm_cluster_options --constraint=X` reached sbatch as `#SBATCH true` | Nextflow reads a value after a space as a separate flag; the launch now refuses it, and the README shows the `--name='value'` form |
| A failed launch check left no message in `.nextflow.log` | Its output joins the Nextflow error |
| Each resource check read 171,914 files one stat and realpath at a time: about 3.5 minutes, eight times per run | One lstat per entry, eight entries at a time; the manifest is unchanged and a pass takes about 55 seconds (raw local run 2,179 s → 1,325 s) |
| `-resume` re-ran three of four raw tasks; one launch was refused with "Resource changed: cadd_base_dir" | This network file system can report a folder's size differently between reads. The bundle passes between tasks as its files, sorted by name, and folder records no longer compare their own size; the refusal now names each changed field |

CADD's RegSeq feature columns differed by 1e-06 in 24 of 13,825 rows when
the Slurm annotation task had 10 CPUs instead of the 12 the shell and local
runs saw; RawScore, PHRED and every family output were unchanged, and with 12
CPUs the table was identical. Keep `--annotation_cpus` fixed to compare runs.

The four tasks stay as they are. The annotation phase is a chain in which
each step reads the previous step's output (the gnomAD allele-frequency
filter decides which variants reach VEP and CADD), so finer tasks would run
one after another anyway. Filtration and prioritization take seconds; the
time a task adds is spent in the resource checks, not in the task layout.
Families run in parallel as their own tasks (next section). Splitting would also
mean breaking `annotation_vcf.sh` into separately callable steps, which the
shell workflow on `main` shares.

Not covered by the trio: raw mode inside the Apptainer image (annotated mode
passed there) and the `task_readonly` profile (the cluster has no Bubblewrap).

## Whole cohorts, and the Snakemake workflow removed

The `Snakefile` ran the phase scripts for a whole cohort: one annotation,
then filtration and prioritization for every PED family, or once without a
family when no PED was given. Nextflow now does the same, and the
`Snakefile` (with the example's `smoke.Snakefile`) is removed. The
`snakemake` package stays in the environments, because upstream `CADD.sh`
is itself a Snakemake workflow.

| Step | Tasks |
| --- | --- |
| Validation | One: the VCF, every PED family (or the one `--family` names), the resources |
| Annotation | One, over the cohort VCF with the full PED |
| Filtration, prioritization | One each per family, in parallel; with no PED, one each without a family |

The task contract records `families` and `family_samples` per family;
`priva task run --family ID` names the family of a family task, and refuses
a family the contract does not list. Each family's first PED row must be its
affected proband. Without a PED, prioritization reads the whole annotation
CADD table, as the Snakefile's no-pedigree rule did.

Compared with the Snakefile at `02d16b6` on the BH14233 trio slice, with its
three samples as three one-person families and again without a PED (job
1538573 at `1a07641`, 9 November 2026): the annotated VCF, the annotation
CADD table, and every family's filtered VCF, filtered CADD table, final table
and ACMG matrix were identical; `-resume` reused all eight cohort tasks.

## Annotation as one task per step

In raw mode with both hubs off, `annotation_vcf.sh`'s functions run as their own
Nextflow tasks, called with the arguments `main_workflow` gives them and without
changes to the script: `priva task run --phase annotation-preprocess`, then
`annotation-gnomad`, `-clinvar`, `-control`, `-vep` and `-cadd` side by side on
copies of the preprocessed VCF, then `annotation-merge`. The merge
(`rust/src/nextflow_task/info_merge.rs`) reads the preprocessed VCF and the four
step VCFs record by record and appends each step's new INFO fields and header
lines in the script's order; it refuses a step that kept a different record or
changed QUAL, FILTER, genotypes or an existing INFO field. Each step checks only
the resources it reads (`phase_reads`), so a step does not re-read the 148,287
files of the CADD installation. The gnomAD step removes no variants (`af_cutoff`
is applied later, in filtration), so every step sees the same variants.

Compared with the Snakefile at `02d16b6` (job 1538935 at `d74eef8`, 12 CPUs as the
reference, 10 November 2026): the annotated VCF, the annotation CADD table and
every family's outputs were identical for three families and without a PED;
`-resume` reused all 14 tasks. In a 24-CPU job (job 1538842) only CADD's RegSeq
columns and 3 RawScore values moved by 1e-06, from the CPU count CADD sees.

Thread counts are left out of the AlphaMissense annotation and DAS completion
records (`d9eac64`): the G9 AlphaMissense cache built with 22 threads is now
reused at 8. VEP with 2 and 10 forks differed only in `HGNC_ID`, which PriVA does
not read (`fork_1538819`).

## Prioritization as five tasks, and the annotation hub

Prioritization runs as `PRIORITIZE_COMBINE`, `PRIORITIZE_SPLICING` and `PRIORITIZE_UTR` side by side, `PRIORITIZE_JOIN` (`rust/src/nextflow_task/table_join.rs`) and `PRIORITIZE_ACMG` (the ACMG part of `assign_acmg_criteria`, now `acmg_assign_table`). The script finds its helpers from `BASH_SOURCE`, since the step wrapper sources it. Family tasks no longer check the annotation-only folders (VEP cache, VEP plugin folders, CADD, LOFTEE); every other resource is still checked before and after each task.

`priva annotation-hub lookup|fill|add` (`rust/src/annotation_hub.rs`) keeps one step's results per variant in `HUB/SOURCE/LINEAGE/CHROM/part-SESSION.parquet`; validation records the five lineages in the contract (`rust/src/nextflow_task/hub_lineage.rs`).

| Job | Checked | Result |
| --- | --- | --- |
| 1539590 | Five prioritization tasks against the Snakefile, three families and no PED | Every table and ACMG matrix identical; resume 26/26 |
| 1539604 | Same, after leaving annotation-only folders out of family checks | Identical; three families 1,491 s (2,752 s before, 1,524 s with one task) |
| 1539608 | `--annotation_hub`: empty hub, then a new session with the filled hub, then `-resume` | Both runs identical to the Snakefile (annotated VCF, CADD table, every family table and ACMG matrix); second run 7,761 of 7,761 variants stored for gnomAD, ClinVar, VEP and CADD, no step but control ran; resume 32/32 |

## Annotated mode and the old hub options removed

The whole-script annotation task (`ANNOTATE_COHORT`), annotated mode (`--mode annotated --cadd`) and `--vcf_hub_mode`/`--cadd_hub_mode` are removed from the workflow; setting them stops the launch with a message that points to `--annotation_hub`. Every run annotates one task per step; `--annotation_hub` reuses stored annotation and `-resume` finishes an interrupted run. The `priva task` adapter still accepts the earlier mode values, which the workflow now always passes as `raw` and `off`. `--action hub_refresh` (`priva annotation-hub refresh-plan`) re-annotates the variants of a step's earlier hub folders into its new folder after a resource changes; trio job 1539801 (after hub lineages left paths out and hashed bundled files) filled a fresh hub, changed only the ClinVar resource, and the refresh re-annotated ClinVar alone (7,761 variants; gnomAD, VEP and CADD stayed current); the next cohort run found every variant stored and its outputs were identical to the Snakefile, and a second refresh found nothing to do. Job 1539721, before that change, had also refreshed VEP, because two bundled files changed path and time between checkout copies.

## Hub lineages from the install lock

`--resources_lock resources.lock.json` (from `--action install`) names each annotation step's hub folder from the lineages of the install steps that produced the resources it reads, one record per resource even for a folder such as the VEP cache; resources the lock does not list keep their own identity (content for bundled files; size with digest or modification time otherwise). Trio job 1540195, with the lock of install job 1539537: empty hub then filled hub, both identical to the Snakefile, second run 7,761 of 7,761 stored for gnomAD, ClinVar, VEP and CADD, resume 32/32; every lineage differed from the run without a lock, as expected.

## Installation as one task per step

`nextflow run main.nf --action install` runs the installer as Nextflow tasks (`workflows/install.nf`, generated from `conf/install_steps.json`; `tests/install_workflow_check.py` fails when the two differ). Nextflow 25.10's strict syntax refuses `-entry`, so the main workflow selects it by `--action`.

| Part | What it does |
| --- | --- |
| `INSTALL_RELEASE` | `scripts/install_task.sh release`, once per step with a changing source (registry `release`): the current git commit, or the server's ETag, Last-Modified and size (an answer without ETag and Last-Modified, such as an API, by the SHA-256 of its sorted lines: the ClinGen evidence repository returns the same rows in varying order). Pinned steps run it only when `--refresh` names them or their pin can expire (`max_age_months`) |
| `INSTALL_PLAN` | Copies your configuration without any key an install step writes, and records each step's code digest (the text of every installer function its calls reach, plus the script files they name), its settings and the files you supply |
| `INSTALL_STEP` (one alias per step) | `scripts/install_task.sh run`: a private configuration whose `cache_root` is `resource_root/STEP/LINEAGE/data`, the keys of the input receipts, then the unchanged installer functions with `PRIVA_CACHE_INSTALL_INTERNAL=1`; the output keys are checked and `receipt.json` is written; the folder is then made read-only unless the step is marked writable (CADD) |
| `INSTALL_PUBLISH` | One configuration from every receipt, `priva cache-install prepare` and `publish` for the native statistics files, and `resources.lock.json` |

A step's lineage hashes its code digest, settings, supplied files, the release of a changing source (only for steps that have one), the `priva` source identity (only for steps whose results `priva` computes) and the lineage of every step it needs. SpliceVault (B5), MaveDB (B6), UTRannotator (B1) and LOFTEE (B8) are pinned: their lineage uses the release recorded in `resource_root/pins.json` (`initial` before any refresh) until `--refresh B5,B8` checks them again, and publish then records the new release there with the day it was taken (`{"release": ..., "recorded": "YYYY-MM-DD"}`; a release already pinned keeps its day). The ClinVar download (D0) is pinned with `max_age_months: 3`, the original installer's refresh interval: its release is kept until three calendar months after the day it was taken (the day clamped to the month's length, as `add_calendar_months` does), then the current release is used. `storeDir` reuses a lineage folder that holds a receipt, across runs and after the work folder is deleted. A step that rewrites an input in place (`copies`) works on a copy in its own folder. `--adopt_config` with `--adopt all --run_steps ID,...` writes receipts for an existing installation's files instead of building them.

Most steps are one existing installer function. Two are split further:

| Steps | What they do |
| --- | --- |
| `C1_chr1` … `C1_chrY` | `gnomAD_download_chrom`: one gnomAD v4.1 joint sites chromosome each, fetched and indexed in its own lineage folder (1 CPU, 4 GB; `withName` selector). `--adopt` finds an older installation's downloads through the registry's `adopt_paths`, because its configuration names only `gnomad_vcf_chrX` |
| `C1` | `gnomAD_collect_chroms` links the 24 downloads into its folder, then the unchanged `gnomAD_install` finds them complete and records the key, or lifts them over for GRCh37. `gnomAD_install` keeps its own loop through the same per-chromosome function (`_gnomAD_fetch_chrom`) |
| `D0` | `ClinVar_download` (`priva clinvar-gene-refresh download`): detects the latest ClinVar release and acquires its VCF, index, variant_summary and MANE with the existing acquisition code; writes `clinvar_release_inputs` |
| `D1` | `ClinVar_install` as before, but the refresh plan judges D0's release instead of asking NCBI, the acquisition checks D0's files by MD5 instead of downloading, and the branch that keeps the bundled gene statistics deploys D0's VCF |

The runner turns `pipefail` off while an installer function runs, as the original installer does: checks written as `bcftools view FILE | head -200` otherwise fail on any VCF longer than 200 lines.

Wiring test with stub steps (job 1540683, 57 steps): 58 folders on the first run, none new on the second; a changed `mavedb_url` gave new folders for exactly B6, C10, C11, D1, D6, D8, D9, F9 and publish; `--refresh B5` the same set with B5 and recorded B5 in `pins.json`, and the next run without it nothing new; D0 pinned to another release taken today gave new folders for D0 and the steps it reaches, and the same pin dated 2020 had expired, so the current release (and the earlier folders) came back. Live release checks of all eleven sources: `validation/install_workflow_20261110/release_live2`, stable on repeat (`release_live4`, `release_live5`).

Real gnomAD test (job 1540694, `validation/install_workflow_20261110/gnomad_1540694`): with every other step adopted from the existing installation, `C1_chrY` downloaded chrY into its own folder with the same bytes as the existing chrY file, `C1` linked all 24 chromosomes (chrY to the new download, the others to the existing files) and ran `gnomAD_install`, and a second launch submitted no step task.

Real ClinVar test (job 1540983, `validation/install_workflow_20261110/clinvar_1540983`, 24 CPUs): with every other step adopted from the existing installation, `D0` found the release dated 2026-11-08 and acquired it, and `D1` annotated it (VEP 34 hours, gnomAD and the AA/splice statistics about 2.5 hours; 37 h 16 min in all). The annotated ClinVar VCF carries D0's fileDate; `D1` kept the bundled gene statistics, so it deployed D0's VCF without an acquisition record. `pins.json` recorded D0 with the day it was taken, and a second launch submitted no step task (128 seconds). Like the existing installation, `D1`'s folder keeps the 24 per-chromosome VEP files (about 5.5 GB). An earlier attempt with 8 CPUs (job 1540692) ran out of its 12-hour limit during VEP.

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
