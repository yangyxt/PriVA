# PriVA scripts

RefSeq candidate acquisition, 2026-10-19: `cadd_refseq_resources.py` exposes
`plan`, `acquire` and `check --config CONFIG` for the pinned GRCh38/merged-VEP116
RefSeq source under configured `cache_root/references/RefSeq/GRCh38/`.
Source identity, publisher MD5/length, local SHA256 and owned partial recovery
are covered by 23 synthetic tests per checkout. The closed log is
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/deployment/cadd-refseq-candidate-helper-tests-20261019-v1.log`;
the mirrored source/tests remain under
`/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/{PriVA,PriVA-nextflow}/{scripts,tests}`.
Acquisition is standalone; transcript/coordinate validation and profile activation
are separate. Job 1515736 completed 0:0 on 2026-10-19: all four original files
(242277346 B) passed publisher MD5/length and local SHA256 checks under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/caches/references/RefSeq/GRCh38/GCF_000001405.40-RS_2025_08/`.
DataLad 0d159924 records the closed bundle/logs while preserving physical identities.
The candidate remains unqualified and inactive; no model ran.

Current installer recovery checkpoint, 2026-10-18: InterPro consumers honor the
selected raw/gzip mapping; unexpected chromosome failures abort publication.
HPO/MONDO and DAS multi-output/intermediate caches now have checked recovery.
DAS CLI and PM1 regeneration require explicit completion after the original
writes; early exits cannot register stale results. Twelve labelled boundary
cases and 24 existing reader tests passed per checkout; full-model checks remain separate.
DAS metadata binds to the selected VCF, assembly and extraction source; only
managed derived metadata regenerates when that association is missing/stale.
Ordinary valid pickle values/order and bundled existence/nonzero policy remain.
Config-file arguments are quoted through installer dispatch. The root README
has exact focused-test receipts and separate pending full/scientific acceptance.
Source/tests live here; all logs and synthetic runtime fixtures are under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/{logs,deployment,validation,tmp}`;
maintained harnesses are under
`/stornext/snfs170/next-gen/scratch/yangyxt/projects/PriVA_upgrade/scripts/`.


Updated: 2026-09-14. This directory contains PriVA's installation, annotation,
cache-building, and interpretation scripts. Existing script headers and the
repository README document their interfaces. This inventory records the files
changed for installer startup; other upstream scripts retain their prior roles.

| File | Purpose and change provenance | Status |
| --- | --- | --- |
| install_utils.sh | Resource installer based on commit 810e60165a663bd2e1223809af9f32d664a9d543. Preserves the 2026-09-13 manager discovery and fixed priva_acmg bootstrap. Adds cache_root plus six optional folder names, validates/resolves paths after environment activation, updates saved descendant paths and creates selected directories. A --check entry point previews without writes. CADD follows the selected parent without its legacy storage prompt. | Updated 2026-09-14; 14 focused checks passed on macOS and RHEL 8.10; full resource deployment not performed. |
| spliceai2_converter.rs | Standalone installation-only, standard-library Rust projection extracted from the existing converter; `project` preserves native-row validation, gene/strand joins, source order and intact MNV alleles. Built directly with rustc, without production PriVA/Cargo dependencies. | Added 2026-10-17. |
| common_bash_utils.sh | Shared shell helpers; quote source-path discovery so loading the installer works from directories containing spaces. | Targeted compatibility change, 2026-09-13. |
| annotation_vcf.sh | Annotation workflow functions loaded by the installer; quote source-path discovery. | Targeted compatibility change, 2026-09-13. |
| README.md | Change inventory and pointers to existing script documentation. | Added 2026-09-13. |

Installer repair, 2026-10-16: `install_utils.sh` and the narrowly scoped
`{installer_prerequisites,cadd_cache_install,spliceai_cache_input}.py` helpers
validate original inputs and downloaded/indexed resources before use. The shell
fixes configured references, all-chromosome gnomAD routing, Conservation ordering,
cache-directory creation and download/parser/index failure propagation. The
portable ClinGen `awk` invocation retains its original expression. AlphaMissense
and PM1 receive the original mapping explicitly; no scientific Python is edited.
CADD retains checksum-verified archives/partials and records extraction identities;
existing mismatches require explicit recovery, never overwrite. SpliceAI2's existing
TSV adapter and annotation-only guard are shared by both installers, without new
YAML switches. `common_bash_utils.sh` accepts official `.vcf.bgz` inputs and stops
failed BGZF preparation before CrossMap. The main README records required manual
inputs and unresolved scientific gates. Existing uncommitted work was preserved.

The original installer also now uses the previously tested unattended VEP cache/plugin route (`--NO_UPDATE`, selected cache release), verifies required cache/plugin files after a zero exit, and retains the Conda-provided official VEP API.

BioMart mapping update, 2026-10-16: `biomart_nm_to_enst.py` retrieves/verifies the
explicit GRCh37.p13 mart116 crosswalk, preserving five-column bytes, ordering and
multiplicity; `installer_prerequisites.py` resolves the bundled or explicit input.
`install_utils.sh` passes that path to both unchanged motif/PM1 Python scripts and
binds derived-cache reuse to mapping/input provenance. The bundle and sidecar are
`data/hgnc/ensembl_GRCh37_mart116_nm_to_enst.tsv*`; exact size/SHA256, query/source
identity, acquisition jobs, partial HCSeeker coverage and scientific limits are in
[the current mapping record](../README.md#biomart-mapping). Earlier statements that
only the unavailable original table could unblock acquisition are superseded by
the user's explicit authorization of this new source. No original-profile identity
or held-job release is implied.

Reference acquisition, 2026-10-16: `reference_resources.py` and the small
`install_utils.sh` reference wrappers provide selected/both-assembly FASTA and
GoFCards forward-chain routes. Exact UCSC source pins, resumed checksum-checked
staging, samtools indexes and provenance protect automatic defaults. Explicit
original references are preserved; unresolved supplied paths fail without
substitution. Configuration updates quote paths and occur atomically only after
all requested resources validate. The main installer requests both assemblies
only for a fresh GoFCards build; the opposite merged VEP cache remains manual.
Source URLs and limits are in the root README. No scientific Python changed.

VEP release selection, 2026-10-16: the existing `vep_installed_version` scalar now
selects the numeric major for `main_install CONFIG` and `conda_install_vep CONFIG`.
`install_utils.sh` verifies executable/API/installer agreement, selects official
release plugins and assembly caches, and preserves all other installed Conda
package versions/builds during an in-place major change. `vep_plugins_release.py`
validates module names and hashes, records completed installation, and preserves
transition intent for safe retries after partial publication. Custom plugin bytes
are protected. Cache checks require assembly/species metadata, tabix conversion
and a nonempty CSI or TBI index. These are installation consistency checks, not
scientific acceptance. The live VEP113 environment and frozen Nextflow deployment
recipe were not changed; see the root README and tests README for validation.

Actual VEP 116 replacement, 2026-10-16: both installers omit the redundant
`--freeze-installed` selector while retaining exact constraints on every other
installed Conda package. Real solver checks and upgrade job 1514975 changed only
VEP 113.2 to 116.3; 527 other package identities remained unchanged. The current
cluster recipe/configs request 116.3/116; frozen historical snapshots remain
preserved. Official plugin job 1514978 loaded all 12 required modules with API 116
and preserved LoFtee. Core cache deployment and scientific acceptance are separate.

SpliceAI2 MNV recovery, 2026-10-16: `spliceai2_to_vep_vcf.py` now explicitly
qualifies official113/116 plugin identities; `test_spliceai2_plugin_identity.py`
passes 5 cases per checkout. The shared optimized Rust converter accepts intact
multi-base substitutions in publisher indel shards. Full deployment/test limits
and recovery job1515009 are recorded in the repository README. No scientific
Python, YAML resource keys or annotation sequence was changed by this correction.

ClinVar regeneration repair, 2026-10-17: `clinvar_stat_variants.py` now requires
explicit `--clinvar-summary` and `--mane-summary` inputs and passes empty unused
CADD/HPO maps to the table reader. `clinvar_annotation_contract.py` checks the
actual VEP/plugin and population fields required by each statistics stage. The
shared reader keeps its strict CLNCSQ default; only explicit ClinVar source mode
may omit that unused field. Nullable splice scores, LoF, domain and event scalars
remain missing instead of raising exceptions. Transcript mapping, review filters,
thresholds and category order are preserved, including the existing NMD label
substring behavior. The bundled `clinvar_2star_stats.pkl.gz` stays the baseline;
new pinned ClinVar/VEP116 outputs are separate candidates pending qualification.

Unified CADD runtime installation, 2026-10-17: `cadd_runtime{,_install}.py` and
`cadd_runtime_pins.json` share the existing image/interpreter checks between both
installers. The helper selects `annotate`, `calculate` or `both` (default), verifies
source/image/assets and preserves the standalone 4194be3 interface with original
Python through PHRED. Calculation-only requires neither annotation archives nor
feature-model environments. The Nextflow deployment module re-exports its prior
runtime APIs. `cadd_cache_install.py` reuses an existing `cache_root/.downloads/GRCh*_v1.7.tar.gz`
in place when no canonical archive exists; a present canonical mismatch still
fails. Existing sources, wrappers, original caches and historical receipts remain
unchanged. Operational path overrides use `--help` and `PRIVA_CADD_*`; cluster
receipts use the explicitly selected data-side state directory. The retained
wrapper contains the existing HGSC bind roots.

Real calculation-only installs passed for GRCh37 and GRCh38 in Slurm1514363,
including full image hashes, original calculation assets and the exact frozen
standalone executable copy. Receipts match `deployment/cadd-install.GRCh*.calculate.*.ready.json`
under `yangyxt_dt/projects/PriVA_upgrade`. Frozen continuation1515129 waits for
1515113, then verifies both complete annotation trees and selected runtimes.
No VEP or model inference was performed, and no scientific profile was activated.
Legacy110 resources and PriVA116 shared-feature qualification remain distinct.

ClinVar installer completion, 2026-10-17: `ClinVar_VCF_deploy` now uses the
complete original VEP operation before gnomAD annotation, resumes a completed
VEP stage and checks actual annotation fields before success. The corrected
quoted plugin-directory list loads the configured LoFtee source; the original
comma argument instead loaded different Conda modules. AA/splice and gene/exon
AF outputs publish only after complete pickle/schema validation. Existing
`ClinVar_Gene_stat` still reuses the bundled baseline. The explicit separate
`ClinVar_Gene_stat_regenerate CONFIG SUMMARY MANE OUTPUT.pkl.gz` installer action
passes pinned reference inputs to the repaired reader and rejects baseline output
paths. Its output is a candidate and does not activate a new gene-stat baseline.

CADD fresh source acquisition, 2026-10-17: a missing configured source directory
now uses `cadd_runtime_install.py --source-only` to acquire the reviewed exact
Git revision under the configured parent. The configured archive URL must match
`cadd_runtime_pins.json`; existing non-Git source trees are retained and rejected
by the runtime source check. Source acquisition does not download image or score
assets. No existing source tree is replaced or reset.

Installer shell follow-up, 2026-10-17: `install_utils.sh` declares the LoFtee cache
directory before use and resolves GRCh37/GRCh38 aliases locally to the existing
hg19/hg38 motif, domain-filter and PM1 bundles/arguments/destinations. Canonical
hg19/hg38 behavior and scientific operation order are retained. The bundled
ClinVar baseline diagnostic now describes the separate explicit-input candidate
regeneration route. Queued frozen job snapshots were unchanged.

Conservation integrity follow-up, 2026-10-17: `conservation_resource_check.py`
uses the already-declared `perl-bio-bigfile` reader and the actual configured
reference FASTA index. Reuse/download checks require the assembly's chr1 length
and every shared normalized chromosome length; extra track contigs are reported.
Existing custom files retain their own SHA256 identity without a release claim
from the download URL. Fresh official release114 downloads additionally require
9,611,134,022 bytes and publisher BSD checksum `30633 9385874` before atomic
publication and YAML registration. Invalid configured inputs stop; failed partials
are retained. The hg19 source and runtime `Conservation,...,MAX` remain unchanged.
No scientific qualification or old-current identity is established.

Installer edge-failure follow-up, 2026-10-17: LoFtee acquisition/Git/configuration failures now stop installation; reuse validates the ancestral FASTA gzip/index, SQLite integrity and same-assembly GRCh38 BigWig before success. Git operations preserve the caller directory and the original update condition; failed downloads remain for explicit recovery. ClinVar canonical complete-cache reuse registers its resolved path, while explicit custom inputs are checked without repair or substitution. InterPro retains the exact downloaded XML plus its observed SHA256 for reuse/rebuilds and rejects changed retained bytes. Scientific source URLs, plugin options and parsing/statistics order are preserved.

Configured HPO path follow-up, 2026-10-17: `prioritization_vcf_per_fam.sh` now passes the existing `hpo_assertions` selection to `combine_annotations.py`; an unset key retains the legacy source-relative file. Explicit missing files stop before conversion. The selected file and configuration timestamp participate in table freshness. Installer preflight permits a data-side builder output while retaining old-host/type/inspection guards. Counterpart shell: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts/prioritization_vcf_per_fam.sh; counterpart preflight: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts/installer_prerequisites.py. Nextflow's matching manifest selection is in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/bin/priva_task.py. HPO matching, pandas operations and stage order are unchanged.

Explicit generated-output routing, 2026-10-17: `InterPro_parsing` honors an absolute `interpro_mapping_pickle` exactly; retained XML and observed SHA256 cache metadata stay beside that pickle. Its format-comparison report is for manual curation and keeps the established source-side `data/InterPro/Interpro_format_report.txt` destination. `ClinGen_deploy` honors absolute `clingen_map` and `gene_dosage_sensitivity` independently; the eRepo source stays beside the selected map and dosage staging stays beside the selected dosage table. Bare-name legacy destinations, source URLs, parser/filter options and operation order are preserved. Existing complete bare InterPro caches register their resolved path without rebuilding. Counterpart installer: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts/install_utils.sh. No production paths or frozen snapshots were changed.

Managed-reference cache routing, 2026-10-17: installer wrappers explicitly opt into `reference_resources.py --use-cache-root`. With the existing `cache_root` key set, new managed FASTAs/indexes use `cache_root/reference/ucsc.hg{19,38}.fasta` and the forward chain uses `cache_root/liftOver/hg19ToHg38.over.chain`. Existing selected paths remain unchanged; old managed defaults keep their publisher/receipt checks. Only absent exact legacy defaults without retained index/archive/staging/receipt/lock state can redirect. Missing custom inputs stop. Preflight and installation share the same planner; configuration paths are saved only after every requested resource succeeds. Direct helper calls retain explicit `--data-dir` semantics unless opted in. Supplied external FASTAs still use adjacent indexes/repair and the existing adjacent-lock contract; no relocation feature was added. Counterpart helper: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts/reference_resources.py; counterpart wrappers/preflight: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts/install_utils.sh and /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts/installer_prerequisites.py. No source URLs, release pins or scientific order changed.

Final installer input/I/O repair, 2026-10-17: UTR acquisition retains the caller directory. Configured CDS FASTAs retain their exact path spelling, including aliases, after complete gzip/plain-FASTA and declared-assembly checks; missing/invalid custom inputs stop without a replacement download. New managed CDS downloads validate before publication and retain failed staging bytes. Alpha manual index retries and mechanism/pext/MaveDB configuration failures now propagate. GRCh38 and hg38 use the same existing MaveDB assay-metadata route.

Explicit ClinVar VCV selection uses the optional operational `--clinvar-vcv` builder argument in both source and full phases. Supplied XML is fully read/structurally checked and identified by local path, size, MD5 and SHA256, without publisher/release qualification. Custom GoFCards reuse/injection binds both XML spelling and SHA256; changed input fails before publication, and the full builder requires matching injected provenance. Literal managed weekly paths retain the original weekly/lag policy; differently named aliases are explicit inputs. No XML classification/review filters or source order changed.

CDS runtime cache naming now uses only the filename suffix: standard `.fa.gz` keeps its existing adjacent cache name; other filenames append `.cds_cache.pkl.gz`, preventing FASTA overwrite. Shell CDS arguments remain one argv item, and Nextflow inventories/copies the same safe cache path into private task inputs. This is an I/O-only change to the detector constructor; sequence parsing, alternative-start calculations and criteria are unchanged. Counterpart installer/reader: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts/install_utils.sh and /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts/alternative_start_codon.py. Nextflow adapter: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/bin/priva_task.py. Frozen job snapshots/configurations were not edited.

Portable installation configuration and bundled-data policy, 2026-10-17: `config.yaml` is now a host-independent first-install template. Set `cache_root` and `tmp_dir`, select assembly/VEP, and leave workflow inputs until a workflow run. Basic scalars are checked before environment changes; the full effective configuration is preflighted in selected temporary storage before updating the source config. The original six-folder managed-root remapping is unchanged, and explicit external paths are retained. `installer_config_paths.py resolve --config CONFIG --source-root CHECKOUT [--check]` fills unset storage paths and resolves bundled builders/curated tables relative to the checkout. Existing canonical InterPro, ClinVar gene statistics, GoFCards/mechanism and assembly-specific ClinGen bundles use existence and nonzero size only; LOEUF checks the selected bundled file before registration. Missing automatic bundles and forced unset outputs use managed storage. Explicit force against a canonical bundled output fails without overwriting it. External/downloaded resources retain their existing checks and freshness rules; paired non-LOF/GoFCards reuse skips only the unnecessary VCV source phase. Scientific Python readers and annotation/resource step order are unchanged. Counterpart helper: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts/installer_config_paths.py; counterpart installer/preflight: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts/install_utils.sh and /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts/installer_prerequisites.py. No live or frozen-job configuration was edited.

2026-10-17: `cache_resources.py` inspects the root lifecycle map without touching payloads. `clinvar_gene_refresh.py` reads its requested three-calendar-month policy, preserves bundled statistics and source scientific order, and selects a completed managed output by full path. Maintained source is `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts`; closed validation logs/receipts are at `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/{logs,deployment}` with `clinvar-gene-refresh*20261017*` and `installer-policy-*-20261017-resource-map-*` patterns. The data-side README points back to both checkouts. These checks do not establish scientific or family acceptance.

Standalone SpliceAI2 installer and MaveDB cleanup, 2026-10-17: `spliceai2_to_vep_vcf.py` exposes `prepare`, `build-converter --tools-dir DIR [--rustc PATH]` and `gene-context` procedures. Fresh conversion builds only `spliceai2_converter.rs` into the derived `SpliceAI2/.tools/spliceai2-converter/SOURCE_SHA256/` cache; build receipts bind its own source/binary/compiler identities. Rust 1.70+ is sufficient for the standalone source (edition2021 and `is_some_and`/`is_ok_and`); absent a compatible host compiler, the existing Conda/Mamba manager acquires private `rust=1.90.0` under that tool root, with a private package cache. Completed adapter receipts retain their original identity and reuse path without compiling. The unused `MaveDB_deployment` bulk-download/neutral-band helper and its active resource-registry entries are removed; production `MaveDB_install`, bundled curated assay bands and approved hg19 omission remain. No full score conversion or scientific activation is part of this change.

Closed standalone-converter check, 2026-10-18 UTC: Slurm 1514363 step 651 passed 128 focused cases per checkout and both registry/shell checks. Publisher prefixes (64 SNVs and 65 indel-shard rows, including an intact MNV) have byte-identical projected rows and diagnostics versus the retained production converter. Receipt and samples: `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/standalone-spliceai2-installer-20261017-v4/`; closed logs: the matching `standalone-spliceai2-installer-*-1514363-20261017-v4.log` files under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/logs/`. Source hashes remained unchanged during the run. The source/scripts stay in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts; generated binaries, samples and logs remain data-side.

ClinVar installation candidates, 2026-10-18: `annotation_vcf.sh anno_VEP_cache_candidate --config CONFIG --input_vcf MANAGED_VCF [--threads N]` prepares the existing VEP plugin sequence as an explicitly unqualified managed cache. `clinvar_gene_refresh.py candidate-route --source-root ROOT` checks source wiring before inputs exist; `candidate-begin`/`candidate-complete` bind actual inputs and publish `*.installation-candidate.{pending,ready}.json` for the VEP stage before gnomAD enrichment. Production entry points retain the original SpliceAI2 guard. Missing/empty default ClinVar gene statistics rebuild from the existing pinned same-assembly ClinVar/MANE route into managed storage; the bundled file is never recreated or overwritten, and missing explicit inputs still fail. Source remains in `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts`; bounded validation logs/receipts live under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/clinvar-cache-candidate-20261018-v*`.

ClinVar profile workspaces, 2026-10-18: managed regeneration now uses `gene_statistics/ASSEMBLY/FILEDATE/vepMAJOR-PROFILE_HASH/`. The hash derives from the existing scientific-profile identity; recorded SHA256 replaces mtime for hashed code/small inputs. Changed profiles build separately without altering prior cache bytes; unchanged profiles and matching legacy unsuffixed completions reuse their existing results. Input release selection, original scientific order and explicit-input conflict checks are unchanged.

CADD native installation, 2026-10-18: `cadd_runtime_install.py` now creates original Python environments at `cache_root/runtime/cadd-1.7.2` by default. `--runtime-dir` explicitly selecting an existing pinned image reuses the reference container; an explicit container failure does not select native. `cadd_native_exec.sh run --runtime-root ROOT --name NAME -- ARGS...` exposes the named original interpreters without Snakemake. The original four-stage Python calculation wrapper and source/model pins are unchanged; cadd-rs remains suspended, superseding the earlier standalone description above. Native package caches, pip caches, temporary requirements and managed environments remain under the configured cache root. The measured RegSeq distribution is seqiolib 0.2.4, although its module reports 0.2.0. Existing source/model/annotation assets can be reused with completed external member receipts without reacquisition.

Slurm 1515248 proved exact native/reference/frozen-original outputs for all four calculation stages on both 18-row assembly fixtures. Ten additional original PHRED CLI edge cases matched stdout, stderr and exits, including the original malformed NaN behavior; this does not accept malformed artifacts or clear the historical experimental-wrapper failure. The 39 focused installer cases passed per checkout. Proof and installation records match `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/cadd-native-runtime-20261018/*`; closed test logs match `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/logs/cadd-native-runtime-installer-*-v3-20261018.log`. Maintained helpers remain in `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts` and the reciprocal `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts`. Feature-environment deployment and shared VEP116 scientific qualification are separate; no production shared-VEP route is enabled by installation.

The native helper also discovers the existing initialized manager through `MAMBA_EXE`, `CONDA_EXE` or PATH. Conda/Mamba use their supported create/install interface, the same upstream channel and Conda specifications, and a separate original pip phase within the selected prefix. This compatibility branch passed focused command/ownership/cache-routing checks; the real native deployment used micromamba 2.9.0. The expanded suite passed 41 cases per checkout in `cadd-native-runtime-installer-*-v4-20261018.log` under the same data-side log folder. Completed native prefixes reuse without requiring a manager or solving packages again.

CADD source-cache handoff checks passed in the expanded 44-case suite per checkout (`cadd-native-runtime-installer-*-v5-20261018.log`). A byte-identical source at a new managed path reuses existing native package receipts unchanged. ESM checkpoint links can retarget atomically only when the old and selected model paths identify the same hardlinked file; all five links are checked before any change, and different objects remain rejected even if their bytes match. The native interpreter and original calculation wrapper bytes are unchanged. Nextflow's 15 existing deployment API checks also passed (`cadd-native-runtime-publication-v2-20261018.log`).

Physical prepared SpliceAI2 adoption, 2026-10-18: `spliceai_cache_input.py adopt-prepared --source-dir SOURCE --output-dir DEST --cache-root ROOT` copies the completed two-VCF/index/gene-context bundle into managed storage without conversion. The importer verifies original and copied SHA256 values, resumes only matching owned partial prefixes, and publishes the directory atomically. `adapter.original.ready.json` preserves the original receipt bytes; the new `adapter.ready.json` changes output/context locations and records `relocation.ready.json` provenance while preserving all scores, scientific metadata and the unqualified state. Explicit root directory aliases resolve to canonical receipts; linked child directories, shared staging files and missing original locks fail without changing the source. Maintained source is `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts`; fixture logs/receipts are under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/spliceai2-physical-adoption-20261018-v*`.

Selected plugin inputs and interrupted acquisition, 2026-10-18: `indexed_plugin_cache.py` exposes `ensure-index` and `format-primate`. Index repair preserves table bytes and readable publisher TBI coordinate flags, columns, comment and skipped-header settings, with bounded first/adjacent-coordinate query comparisons. Raw PrimateAI uses the original transformation/sort order and records its exact source SHA256 and declared assembly under the managed plugin cache. UTRAnnotator, SpliceVault, MaveDB, pext and LoFtee honor selected paths and reject missing custom inputs. New UTR/LoFtee clones and LoFtee resources stage before validation/publication; existing LoFtee revisions stay fixed during cache recovery. GRCh37 SpliceVault uses configured `tmp_dir`. Mirrored sources are in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts; closed fixture logs, commands and identity receipts are in /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/plugin-selected-inputs-20261018-v*.

Mapping-bound motif/PM1 upgrades, 2026-10-18: `biomart_nm_to_enst.py artifact-route` accepts named artifact, mapping, anchor-input, cache-root, output-directory, assembly and kind arguments. It verifies the old artifact against its receipt before comparing the existing mapping identity. Changed mapping/provenance/anchor hashes select a deterministic managed namespace; completed matching targets can be reused and registered without rebuilding. Missing/invalid receipts, altered outputs and unsafe linked/occupied destinations remain errors. Old artifacts and receipts are retained; scientific options, calculation order, freshness checks and receipt-write fields are unchanged. Mirrored source is in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts; closed focused logs, commands and identity receipts are in /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/mapping-artifact-routes-20261018-v{1,2}.

Original AlphaMissense producer handoff, 2026-10-24 CT: `alphamissense_original_cache.py` exposes `snapshot`, `import` and `check` with named options. It copies an externally attested completed original producer output without loading/recalculating it. The source content hash, selected VCF/index identities, native completion evidence and exact output size/SHA are retained in a separate `.original.complete.json`; retrospective source/input bindings must state their evidence limits. An owned pending journal supports interrupted import; nonowned outputs/receipts, input aliases and cache-root escapes are preserved/refused. `AlphaMissense_stat` checks this receipt before its existing chromosome builder route; invalid/pending original receipts stop. Identical producer content across checkouts reuses the same output; warm checks use unchanged output metadata rather than rehashing the full payload. This is execution/publication proof, not scientific qualification. Mirrored helper: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts/alphamissense_original_cache.py. Closed synthetic fixtures/logs/receipts are under /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/alphamissense-original-handoff-20261024-v{1,2,3}; the original scientific producer and chromosome helper are unchanged by this handoff work.

DAS receipt source portability, 24 October 2026 CT: `das_install_state.py` version3 adds named `--source NAME=PATH`, `--bundled NAME=PATH` and ordered `--parameter NAME=VALUE` roles to its existing `capture`, `check` and `record` procedures. Sources compare content SHA256 across checkouts; dynamic scientific inputs and outputs keep exact path/target/inode/size/mtime identities with node-local device excluded. Actual bundled snapshots use stable logical roles and existence/nonzero checks only; no bundled contents are opened. `AlphaMissense_pick_intolerant_domains` classifies a selected InterPro mapping as bundled only at its genuine canonical bundle origin; external copies and links to external targets remain strict inputs. Thread selection and coldspot/benign filter enablement are recorded in order. The scientific command/options/default thresholds and source calculations are unchanged. Incompatible version2 receipts fail verification and are not promoted. Mirrored helper: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts/das_install_state.py. Closed synthetic logs/results/source snapshots are under /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/das-state-portability-20261024-v1; this is reuse/orchestration evidence, not actual DAS or installation acceptance.
