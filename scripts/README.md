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

Cluster deployment, 2026-10-15: `install_cluster.sh` and
`check_cluster_environment.py` create/check an isolated `priva_acmg` environment,
export a resolved package lock, and deploy an explicitly committed scientific
source snapshot. `envs/priva-cluster.yml` preserves the key scientific versions
from the bundled environment export and includes missing cache-builder imports.
CADD retains its separate upstream Python environments. A successful software
check does not establish annotation-cache readiness or family-result equality.

`deploy_cluster_{resources,plugin_inputs,gnomad}.py` installs independent,
resumable resource groups from the reviewed `envs/cluster-*.json` manifests.
`prepare_cluster_config.py` rewrites old-host paths into a separate deployment
config with empty family inputs. `cadd_container_exec.sh` provides original
CADD Python environment entry points; its VEP110 entry is for reference
qualification only. Production annotation shares the official PriVA VEP113.
The resource stages keep SpliceAI/PrimateAI imports, unsupported-assembly policy,
shared CADD transcript qualification, and downstream cache builds explicit.
Successful downloads alone do not constitute a qualified complete pipeline.

Reference/VEP deployment, 2026-10-16: **1514579 COMPLETED (0:0)** at
08:21:30Z after 8 h 39 m 49 s, using 2 CPUs / 8 GB; batch MaxRSS was
28,420 KiB. Frozen source `fc2186901029acc88cd41bb004828778c1fd1cb8`
produced deployed `deployment/{reference,vep}.ready.json` receipts with manifest
SHA256 `b61cbeb328bceac5b72b8a7f55205976848de58cb988b78c816c38cdb3769b84`.
The 08:46:04Z audit passed full archive hashes, pinned publisher checksums and full FASTA SHA256,
93 FASTA headers / 279 indexed positions, all 6,278 expected canonical-cache
entries, 146 sampled transcript/regulatory chunk headers and 25 variation/CSI
magic pairs. The extracted cache payloads were not fully rehashed, and no real
VEP/plugin annotation was run. The existing isolated config already selects the
reference and VEP cache; no config edit was required. The closed-log command and
frozen code are verified, while the submitted wrapper was unavailable.
Scientific/shared-CADD, public-sample and family qualification remain pending.
AlphaMissense stays held for the original mapping as described below.
DataLad `80bb1f0` records the receipts, closed
`logs/cluster-reference-vep_1514579.log`,
`validation/vep-reference-receipt-verified-1514579.json` and data-project README:
five selected paths, four added/updated, with metadata recorded and selected
paths clean. The read-only audit used Slurm allocation 1514363.

CADD prescore recovery, 2026-10-16: source
`20774bd56a706d8b82a5640bf23b758216544136` binds the Washington GRCh37 v1.7
download root and `MD5SUMs` to the historical indel data MD5
`3d88f8fde9df07060b4f8fe07a0d5d5d` and matching index MD5
`c583ed9f852c9896081205c3523f9f23`. `deploy_cluster_resources.py` fully hashes
an existing partial before attempting HTTP resume; only a matching MD5 can
skip HTTP. Existing final files are fully rehashed, size gates remain active,
and publication stays atomic after verification. Job 1514744 COMPLETED (0:0)
at 2026-10-16T08:02:48Z after 13 m 15 s on analysis-cpa1, using the frozen
`priva-cluster-tools-20774bd` snapshot. It fully rehashed and reused the
602,910,809,417-byte SNV file, existing SNV index and 11,003,494,781-byte indel
partial, then atomically published the indel file. The only HTTP transfer was
the matching 2,644,685-byte Washington indel index, which passed its MD5.
`cadd-prescores.ready.json` records deployed status, job 1514744 and manifest
SHA256 `add789be1ce0180fd0d246175c4289f5f4f2989acd8364e8aa0f36b38c0ff59a`.
Bounded real indexed-pair validation passed in Slurm step 1514363.262 at
08:03:30Z: pysam 0.22.1 verified all 24 contigs in each SNV/indel pair, matching
135-column schemas and exact-coordinate refetches. It did not repeat full-body
MD5/CRC. DataLad `9e70725` records the closed receipt/log/report evidence;
real VEP/model, public-sample, shared-annotation and family acceptance remain pending.
The earlier `fc21869` snapshot and failed 1514581 evidence remain unchanged.
All 20 deployment and 31 licensed-import fixture tests pass; the missing frozen
`054b494` helper fixture was restored with its existing SHA256 unchanged.

`build_cluster_alpha.sbatch` requires five arguments: `LOG CONFIG SOURCE_PREFIX
APPROVED_NM_TO_ENST_TSV SHA256`, plus an explicit scheduler dependency on the
VEP cache job. Missing or changed mapping data fails before activation or
annotation. The wrapper verifies the frozen `09223ef` installer and motif
Python hashes, runs original `AlphaMissense_anno`, then preserves the motif
helper's option order/defaults with explicit `--nm_to_enst_tsv`. It retains the
original timestamp reuse predicates and reports reuse separately; reuse does
not establish the existing output's mapping provenance. Fresh calculations
recheck mapping hash and file identity before completion/config updates.
The mapping is neither selected nor downloaded here: its positional ENST/NM
columns and scientific provenance still need approval. Scientific Python is
unchanged, and held job `1514589` remains an older frozen submission: PENDING,
`JobHeldUser`, even after VEP job 1514579 completed. Updating this wrapper does
not update that job or authorize its release. The original mapping remains
required; a RefSeq/MANE substitute is not qualified.

`diagnose_cluster_install.sbatch`
runs the unchanged historical `main_install` in a separate writable prefix,
using real copies of the exact package lock and a bounded runtime. A timeout is
an operator checkpoint, not successful installation. The restorable package
lock is `environment.conda-explicit.txt`; the earlier `environment.explicit.txt`
is retained only as a historical package URL inventory.

The real historical `main_install` diagnostic exposed VEP's newer-release
prompt exiting zero on noninteractive EOF and PriVA continuing without a cache.
`install_utils.sh` now preserves the installed Conda API, requests cache/plugins
for that exact release with updates disabled, and validates required files
before recording success. The original diagnostic source remains unchanged.

Subsequent installer corrections pass the configured GRCh37 FASTA to gnomAD
liftover, schedule all chr1–22/X/Y, stop after transfer/liftover failures, and
install Conservation before AlphaMissense annotation. Conservation reuse no
longer fails or replaces an existing configured path. AlphaMissense motif
output uses the cluster configuration instead of the retired host prefix.
The scientific Python calculations and per-chromosome transforms are unchanged;
full derived-cache qualification still awaits upstream resource completion.

`build_cluster_gnomad.sbatch` waits for the verified source-download job, archives
an exact installer commit, checks the frozen config and reference identity, then
runs all 24 original per-chromosome liftover transformations. GRCh38 source
`.vcf.bgz{,.tbi}` files remain separate from new hg19 outputs and intermediates.
The completion receipt requires every final VCF/index, its expected indexed
chromosome, and a successful indexed query; this does not establish family or
ACMG parity. Initial job 1514601 used installer
`054b4944055054148005016192647a8ed0384eef` and failed `.vcf.bgz` extension
validation. Replacement 1514729 is running with the minimal validation fix; eight
focused checks passed. Scientific transformations are unchanged.

`build_cluster_splicevault.sbatch` runs the original GRCh38-to-GRCh37 conversion
against the release113 SpliceVault input after checking its publisher checksum
and byte size. Job 1514603 COMPLETED (0:0) at 2026-10-16T07:03:13Z after
6 h 24 m 15 s, using installer `054b494` and scientific baseline `09223ef`.
The `splicevault-hg19.ready.json` receipt records `derived_resource_built` for
GRCh37. Seven SHA256 checks passed with stable identities: source, output,
index, reference, chain, converter and job config. Publisher BSD sum
`[57027, 863191]`, full `bgzip -t` CRC and pysam 0.22.1 queries across all 24
indexed contigs passed. The 12-field source/output headers agree with the plugin
layout; no real VEP/plugin run was performed. The output has 144,086,759
rows from 144,094,769 source rows minus 8,010 original CrossMap misses;
the original mapping policy was preserved. The main config already points to
this output, so no config promotion was needed. Its source URL uses
`current_variation`; the actual job config pinned release113. Frozen
implementation and receipt were verified, but no submitted wrapper-hash archive
was found, so the full submitted wrapper was not reconstructed. Real VEP,
scientific and family qualification remain pending. The newer publisher hg19
table remains separate from this baseline conversion.

`import_cluster_licensed.sbatch` registers explicitly confirmed acquisitions under
`caches/{SpliceAI_v1,PrimateAI_v0_2}`. The optional ninth argument is the resource
selection: `all` (the original default), `spliceai`, or `primateai`. It still runs
the frozen `install_utils.sh` functions; selecting SpliceAI calls only
`SpliceAI_install` and requires both SNV and indel VCFs with their indexes.
Both original collections are acquired. Hg19 SpliceAI job 1514615 and PrimateAI
job 1514620 completed with valid import receipts; the separate three-key hg19
main-config promotion is also complete with a receipt. This script does not download or authenticate.

The confirmation JSON requires `status: downloads_complete`, `confirmed_by`,
matching `assembly`, and exactly the selected `files` roles: `spliceai_snv` plus
`spliceai_indel` for SpliceAI, `primateai` for PrimateAI, or all three for `all`.
Each file entry contains an absolute `path` and optional acquisition
`bytes`/`sha256`. Optional `resources` must match the selection (`["spliceai"]`
for the current partial import); completion refers only to those resources.
The job reports each hash, gzip and index-check stage and verifies source
size/timestamps/inode from before the first checksum through import completion.

Registration writes a separate config and verifies every unselected key remains
unchanged. SpliceAI stays in its shared location; selected raw PrimateAI TSVs
retain their original bytes while the original formatter creates sibling BGZF
and index files. Receipt names are `licensed-caches-ASSEMBLY.ready.json` for
`all` and `licensed-caches-ASSEMBLY-RESOURCE.ready.json` for a subset. Subset
receipts use `selected_licensed_inputs_imported`, list the registered resources,
and mark others `not_selected_unverified`; they do not establish complete
licensed-resource deployment. Main-config promotion and scientific VEP/family
qualification remain separate. A SpliceAI-only submission can override the
combined formatter allocation with `--cpus-per-task=2 --mem=8G --time=1-00:00:00`.

Raw PrimateAI imports additionally require exact external reconstruction before
publishing a receipt. The original `054b494` formatter is unchanged. The
validator preserves the first 11 records after the original `sed` operations,
adds AWK's final newline where needed, and runs the exact GNU sort arguments
with the inherited locale. It explicitly checks each child, then compares the
entire uncompressed derived stream byte-for-byte with `cmp`, including EOF.
This catches the demonstrated partial-sort false success that gzip/index
checks alone miss. One additional sort and a temporary uncompressed reference
under the job's `TMPDIR` are required; normal success/failure paths reap owned
children and remove temporary files. Slurm cancellation remains separate.

Receipts record matched bytes, validation command/tool/locale identities,
source/output stability and the SHA256 of the submitted batch script. A supplied
preformatted BGZF without its raw source is marked `raw_parity_verified: false`.
Raw originals and successfully formatted bytes are never rewritten by this
validator; failed validation leaves no completion receipt or main-config change.
The existing shared importer lock is retained.

## SpliceAI2 annotation-only preparation

Current installation direction: keep only the original spliceai_snv_prescore and spliceai_indel_prescore YAML keys. Choose assembly once at installation: fresh GRCh37 uses v1 and fresh GRCh38 uses v2. Inspect decompressed headers automatically: an existing compatible VCF passes through; native 84-column SpliceAI2 TSV invokes the narrow adapter. A prepared pair is registered under the two original keys and installation succeeds. Resource/model identity, hashes, schema and qualification remain in sidecars and cache headers, without extra selector or converter YAML switches. Existing compatible GRCh38 v1 VCF pairs pass through; fresh GRCh38 installations select v2. Wrong assembly, mixed input types and failed v2 preparation never trigger a fallback.

`spliceai2_to_vep_vcf.py` verifies a complete pinned source manifest, qualifies
scored-gene intervals against exact GENCODE48, and invokes the standalone
installation converter's `project` command through mature `gzip`/`bgzip`/`tabix` processes. It preserves originals and publishes
the receipt last; it performs no acquisition or scientific activation.
SpliceAI_install invokes scripts/spliceai_cache_input.py install automatically, retaining the existing prompts and two resource keys. The helper detects gzip by magic bytes, then classifies a bounded decompressed header as vcf_v1, vcf_spliceai2 or the exact native 84-column SpliceAI2 TSV schema; filenames, extensions and score magnitudes do not select the predictor. An unmarked VCF receives the legacy-compatible vcf_v1 routing label; that label alone does not prove predictor identity. Existing compatible VCF pairs pass through, including existing GRCh38 v1; fresh GRCh38 selects v2. Wrong assembly, mixed types or failing v2 preparation fail without fallback. SpliceAI2 VCF passthrough additionally requires a completed adjacent adapter.ready.json with exact SNV/indel/index membership and complete published-catalogue status. This checks completed-publication metadata and recorded membership, not a new full VCF/index-body checksum; it performs no conversion. Reject .partial or .building path components, including symlink targets, and refuse in-progress preparation or staged outputs. A final VCF renamed before receipt publication is not ready. V1 VCF passthrough requires no additional adapter receipt.

```bash
python scripts/spliceai_cache_input.py install --config CONFIG --cache-dir CACHE_DIR --snv SNV --indel INDEL
```

For native TSV, the helper discovers the fixed adjacent SpliceAI2_GRCh38.downloads-complete.json sidecar. The selected SNV/indel inputs must match their manifest roles and the full 47-shard catalogue. The existing low-level SpliceAI2_install preparation/reuse path performs full checksum verification and prepares or reuses a validated receipt, then updates only spliceai_snv_prescore and spliceai_indel_prescore. Fresh conversion builds/caches the bundled standalone installation converter; no YAML converter or manifest option is added. Provenance and scientific qualification remain in sidecars/cache headers. The current adjacent sidecar at yangyxt/caches/SpliceAI2_GRCh38.downloads-complete.json has SHA256 bf546817c03689359138d9230e3587fe225764a06546b3914c8819afbbd2bdde: all 47 original source identities remain unchanged; the qualified plugin pin is116. The historical manifest remains preserved. The historical complete preparation used yangyxt/softwares/bin/priva selecting priva-tools-spliceai2-mnv-20261016/bin/priva, SHA25670eefad75c6ad5a8001fcd87ca6970c3b0b1495f21b55f8c7775db67122c567a; that production executable is retained as a comparison reference. This changes no YAML keys beyond the original two paths.

Two-key interface qualification: code 35779c4578509ce01218a052a8e7e1d8a785a75b passed 239 Python tests, zero failures/errors/skips, in 190.848 seconds: 39 helper, 25 installer, 14 runtime-guard, 21 Nextflow-adapter and 140 existing deployment/driver/BGZ tests. The installer-shell cases mock the helper boundary; they are not real conversion evidence. Log: yangyxt_dt/projects/PriVA_upgrade/logs/spliceai-two-key-integration-20261016.log. This supersedes the earlier 197-test interface checkpoint without adding its count; unchanged Rust evidence remains 41 tests and the earlier real-projection evidence remains 512 unique rows.

Actual resource checks used final helper SHA256 f69d6ea8d6179d68ae44cad849b3c6c7912a879f05e305a36acffba09c81bb88: all 47 native TSV and four v1 VCF headers classified correctly, including assembly metadata. Real SpliceAI_install reused the existing v1 pair for hg19/GRCh37 and hg38/GRCh38 with configuration bytes, inode and mtime unchanged and resource stat identities unchanged. Receipt: yangyxt_dt/projects/PriVA_upgrade/validation/spliceai-two-key-real-install-20261016/receipt.json. These are bounded header/reuse checks, not new body checksums, full VEP execution or model/scientific qualification.

Nextflow stub verification passed all four scenarios: raw input, resume with all four tasks cached, preannotated input and stale-environment rejection. Log: yangyxt_dt/projects/PriVA_upgrade/logs/spliceai-two-key-nextflow-20261016.log; retained validation/log artifacts: validation/spliceai-two-key-nextflow-20261016/ under that data project. The run used code 35779c4 plus the documentation state at execution; its recorded source_digest identifies that exact state. Subsequent documentation edits are not attributed to the run.

Historical frozen-run boundary: job 1514732 used installer f8de4b3 and converter225aa4a, then failed at the publisher MNV on indel chr1 line11929. It produced no completed adapter receipt. Its explicit-manifest/six-key interface and197 Python tests remain historical provenance. Recovery retains all47 original download identities and the two-key public installer interface; only validated complete outputs can later be imported. Full scientific acceptance remains pending.

The [route reference](../docs/nextflow_migration.md#spliceai2-preparation) provides
commands, the exact JSON manifest contract, all 47 required shards, verified
GENCODE48 and plugin identities, output patterns and the explicit subset-test
boundary. Downloads are complete (1514724; 47 final shards, zero partials).
All 20,083 real regions passed the stable-ID/GENCODE48 join, retaining full
versions and `_PAR_Y`; ten same-strand symbol groups retain all IDs/intervals
and coalesce labels only per prediction. Opposite-strand symbol contexts fail.
The execution manifest is complete and saved in DataLad commit 5f5946b.
All 47 full source SHA256 checks passed in frozen job 1514732. At the bounded
2026-10-16 ~08:45Z scheduler/log check it was RUNNING on c81o-03 and its log
had reached SNV chr4; complete conversion and indexed-VEP qualification remain pending. Publisher expected SHA256 values and VERIFIED acquisition receipts for
all 47 shards are now supplied in [caches/README.md](../../../caches/README.md);
there is no missing-publisher-checksum blocker. Historical Rust 225aa4a passed 32
projection + 9 existing tests; Python f8de4b3 passed 197 tests with zero failures/skips.
A separate 512-row real prefix check (chr1 positive strand only) preserved every
projected record/order and passed bgzip/tabix queries; this is not full-cache or
full-VEP qualification. Optimized source 225aa4aa774b4be7f53576e1cde56af71e67bd2c
is installed at `yangyxt/softwares/priva-tools-225aa4a/bin/priva`, SHA256
`bbc3b7288ed2ea63f1f8455ed87b8a9d08c7bc1fe0a562c1d6134525e2a12dc6`.
Installed-release reruns passed the same 25 driver cases and same 512 source rows;
the frozen checkpoint had 197 Python/41 Rust tests and 512 rows; current interface coverage is 239 Python tests, with unchanged 41 Rust tests and 512 unique rows. A stale shared Cargo
release artifact was detected before deployment; forced rebuilding and installed
smoke checks corrected it, with the failure retained in the receipt. Historical preparation job 1514732 ran on c81o-03 (4 CPUs, 16 GB, 7 days) and failed on the MNV, using frozen
installer f8de4b3, this converter and the isolated GRCh38 config. The manifest/
handoff were saved in DataLad commit 5f5946b before launch; all 47 full source SHA256 checks passed. The bounded 2026-10-16 ~08:45Z check confirmed RUNNING and a log through SNV chr4; complete conversion remains pending.
Original licensed source and derived views remain local under their resource terms.

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

CADD runtime completion checks, 2026-10-16: `deploy_cluster_resources.py` now
verifies five unchanged wrapper routes, all four named Python environments with
network-disabled imports, packaged MMSplice models, calculation models/PHRED
assets and bounded per-assembly annotation-asset readiness. It installs a zero-copy
ESM checkpoint view without changing scientific profiles. An optional separate
`--runtime-receipt` preserves earlier deployment evidence. Frozen Slurm1514901
passed; GRCh38 annotation assets and all scientific qualification remain pending.
The original `cadd_container_exec.sh` and numerical/model pins are unchanged.

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

Configured HPO path follow-up, 2026-10-17: `prioritization_vcf_per_fam.sh` now passes the existing `hpo_assertions` selection to `combine_annotations.py`; an unset key retains the legacy source-relative file. Explicit missing files stop before conversion. The selected file and configuration timestamp participate in table freshness. Installer preflight permits a data-side builder output while retaining old-host/type/inspection guards. Counterpart shell: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts/prioritization_vcf_per_fam.sh; counterpart preflight: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts/installer_prerequisites.py. Nextflow's matching manifest selection is in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/bin/priva_task.py. HPO matching, pandas operations and stage order are unchanged.

Explicit generated-output routing, 2026-10-17: `InterPro_parsing` honors an absolute `interpro_mapping_pickle` exactly; retained XML and observed SHA256 cache metadata stay beside that pickle. Its format-comparison report is for manual curation and keeps the established source-side `data/InterPro/Interpro_format_report.txt` destination. `ClinGen_deploy` honors absolute `clingen_map` and `gene_dosage_sensitivity` independently; the eRepo source stays beside the selected map and dosage staging stays beside the selected dosage table. Bare-name legacy destinations, source URLs, parser/filter options and operation order are preserved. Existing complete bare InterPro caches register their resolved path without rebuilding. Counterpart installer: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts/install_utils.sh. No production paths or frozen snapshots were changed.

Managed-reference cache routing, 2026-10-17: installer wrappers explicitly opt into `reference_resources.py --use-cache-root`. With the existing `cache_root` key set, new managed FASTAs/indexes use `cache_root/reference/ucsc.hg{19,38}.fasta` and the forward chain uses `cache_root/liftOver/hg19ToHg38.over.chain`. Existing selected paths remain unchanged; old managed defaults keep their publisher/receipt checks. Only absent exact legacy defaults without retained index/archive/staging/receipt/lock state can redirect. Missing custom inputs stop. Preflight and installation share the same planner; configuration paths are saved only after every requested resource succeeds. Direct helper calls retain explicit `--data-dir` semantics unless opted in. Supplied external FASTAs still use adjacent indexes/repair and the existing adjacent-lock contract; no relocation feature was added. Counterpart helper: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts/reference_resources.py; counterpart wrappers/preflight: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts/install_utils.sh and /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts/installer_prerequisites.py. No source URLs, release pins or scientific order changed.

Final installer input/I/O repair, 2026-10-17: UTR acquisition retains the caller directory. Configured CDS FASTAs retain their exact path spelling, including aliases, after complete gzip/plain-FASTA and declared-assembly checks; missing/invalid custom inputs stop without a replacement download. New managed CDS downloads validate before publication and retain failed staging bytes. Alpha manual index retries and mechanism/pext/MaveDB configuration failures now propagate. GRCh38 and hg38 use the same existing MaveDB assay-metadata route.

Explicit ClinVar VCV selection uses the optional operational `--clinvar-vcv` builder argument in both source and full phases. Supplied XML is fully read/structurally checked and identified by local path, size, MD5 and SHA256, without publisher/release qualification. Custom GoFCards reuse/injection binds both XML spelling and SHA256; changed input fails before publication, and the full builder requires matching injected provenance. Literal managed weekly paths retain the original weekly/lag policy; differently named aliases are explicit inputs. No XML classification/review filters or source order changed.

CDS runtime cache naming now uses only the filename suffix: standard `.fa.gz` keeps its existing adjacent cache name; other filenames append `.cds_cache.pkl.gz`, preventing FASTA overwrite. Shell CDS arguments remain one argv item, and Nextflow inventories/copies the same safe cache path into private task inputs. This is an I/O-only change to the detector constructor; sequence parsing, alternative-start calculations and criteria are unchanged. Counterpart installer/reader: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts/install_utils.sh and /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts/alternative_start_codon.py. Nextflow adapter: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/bin/priva_task.py. Frozen job snapshots/configurations were not edited.

Portable installation configuration and bundled-data policy, 2026-10-17: `config.yaml` is now a host-independent first-install template. Set `cache_root` and `tmp_dir`, select assembly/VEP, and leave workflow inputs until a workflow run. Basic scalars are checked before environment changes; the full effective configuration is preflighted in selected temporary storage before updating the source config. The original six-folder managed-root remapping is unchanged, and explicit external paths are retained. `installer_config_paths.py resolve --config CONFIG --source-root CHECKOUT [--check]` fills unset storage paths and resolves bundled builders/curated tables relative to the checkout. Existing canonical InterPro, ClinVar gene statistics, GoFCards/mechanism and assembly-specific ClinGen bundles use existence and nonzero size only; LOEUF checks the selected bundled file before registration. Missing automatic bundles and forced unset outputs use managed storage. Explicit force against a canonical bundled output fails without overwriting it. External/downloaded resources retain their existing checks and freshness rules; paired non-LOF/GoFCards reuse skips only the unnecessary VCV source phase. Scientific Python readers and annotation/resource step order are unchanged. Counterpart helper: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts/installer_config_paths.py; counterpart installer/preflight: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts/install_utils.sh and /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts/installer_prerequisites.py. No live or frozen-job configuration was edited.

2026-10-17: `cache_resources.py` inspects the root lifecycle map without touching payloads. `clinvar_gene_refresh.py` reads its requested three-calendar-month policy, preserves bundled statistics and source scientific order, and selects a completed managed output by full path. Maintained source is `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts`; closed validation logs/receipts are at `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/{logs,deployment}` with `clinvar-gene-refresh*20261017*` and `installer-policy-*-20261017-resource-map-*` patterns. The data-side README points back to both checkouts. These checks do not establish scientific or family acceptance.

Standalone SpliceAI2 installer and MaveDB cleanup, 2026-10-17: `spliceai2_to_vep_vcf.py` exposes `prepare`, `build-converter --tools-dir DIR [--rustc PATH]` and `gene-context` procedures. Fresh conversion builds only `spliceai2_converter.rs` into the derived `SpliceAI2/.tools/spliceai2-converter/SOURCE_SHA256/` cache; build receipts bind its own source/binary/compiler identities. Rust 1.70+ is sufficient for the standalone source (edition2021 and `is_some_and`/`is_ok_and`); absent a compatible host compiler, the existing Conda/Mamba manager acquires private `rust=1.90.0` under that tool root, with a private package cache. Completed adapter receipts retain their original identity and reuse path without compiling. The unused `MaveDB_deployment` bulk-download/neutral-band helper and its active resource-registry entries are removed; production `MaveDB_install`, bundled curated assay bands and approved hg19 omission remain. No full score conversion or scientific activation is part of this change.

Closed standalone-converter check, 2026-10-18 UTC: Slurm 1514363 step 651 passed 128 focused cases per checkout and both registry/shell checks. Publisher prefixes (64 SNVs and 65 indel-shard rows, including an intact MNV) have byte-identical projected rows and diagnostics versus the retained production converter. Receipt and samples: `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/standalone-spliceai2-installer-20261017-v4/`; closed logs: the matching `standalone-spliceai2-installer-*-1514363-20261017-v4.log` files under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/logs/`. Source hashes remained unchanged during the run. The source/scripts stay in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts; generated binaries, samples and logs remain data-side.

ClinVar installation candidates, 2026-10-18: `annotation_vcf.sh anno_VEP_cache_candidate --config CONFIG --input_vcf MANAGED_VCF [--threads N]` prepares the existing VEP plugin sequence as an explicitly unqualified managed cache. `clinvar_gene_refresh.py candidate-route --source-root ROOT` checks source wiring before inputs exist; `candidate-begin`/`candidate-complete` bind actual inputs and publish `*.installation-candidate.{pending,ready}.json` for the VEP stage before gnomAD enrichment. Production entry points retain the original SpliceAI2 guard. Missing/empty default ClinVar gene statistics rebuild from the existing pinned same-assembly ClinVar/MANE route into managed storage; the bundled file is never recreated or overwritten, and missing explicit inputs still fail. Source remains in `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts`; bounded validation logs/receipts live under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/clinvar-cache-candidate-20261018-v*`.

ClinVar profile workspaces, 2026-10-18: managed regeneration now uses `gene_statistics/ASSEMBLY/FILEDATE/vepMAJOR-PROFILE_HASH/`. The hash derives from the existing scientific-profile identity; recorded SHA256 replaces mtime for hashed code/small inputs. Changed profiles build separately without altering prior cache bytes; unchanged profiles and matching legacy unsuffixed completions reuse their existing results. Input release selection, original scientific order and explicit-input conflict checks are unchanged.

CADD native installation, 2026-10-18: `cadd_runtime_install.py` now creates original Python environments at `cache_root/runtime/cadd-1.7.2` by default. `--runtime-dir` explicitly selecting an existing pinned image reuses the reference container; an explicit container failure does not select native. `cadd_native_exec.sh run --runtime-root ROOT --name NAME -- ARGS...` exposes the named original interpreters without Snakemake. The original four-stage Python calculation wrapper and source/model pins are unchanged; cadd-rs remains suspended, superseding the earlier standalone description above. Native package caches, pip caches, temporary requirements and managed environments remain under the configured cache root. The measured RegSeq distribution is seqiolib 0.2.4, although its module reports 0.2.0. Existing source/model/annotation assets can be reused with completed external member receipts without reacquisition.

Slurm 1515248 proved exact native/reference/frozen-original outputs for all four calculation stages on both 18-row assembly fixtures. Ten additional original PHRED CLI edge cases matched stdout, stderr and exits, including the original malformed NaN behavior; this does not accept malformed artifacts or clear the historical experimental-wrapper failure. The 39 focused installer cases passed per checkout. Proof and installation records match `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/cadd-native-runtime-20261018/*`; closed test logs match `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/logs/cadd-native-runtime-installer-*-v3-20261018.log`. Maintained helpers remain in `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts` and the reciprocal `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts`. Feature-environment deployment and shared VEP116 scientific qualification are separate; no production shared-VEP route is enabled by installation.

The native helper also discovers the existing initialized manager through `MAMBA_EXE`, `CONDA_EXE` or PATH. Conda/Mamba use their supported create/install interface, the same upstream channel and Conda specifications, and a separate original pip phase within the selected prefix. This compatibility branch passed focused command/ownership/cache-routing checks; the real native deployment used micromamba 2.9.0. The expanded suite passed 41 cases per checkout in `cadd-native-runtime-installer-*-v4-20261018.log` under the same data-side log folder. Completed native prefixes reuse without requiring a manager or solving packages again.

CADD source-cache handoff checks passed in the expanded 44-case suite per checkout (`cadd-native-runtime-installer-*-v5-20261018.log`). A byte-identical source at a new managed path reuses existing native package receipts unchanged. ESM checkpoint links can retarget atomically only when the old and selected model paths identify the same hardlinked file; all five links are checked before any change, and different objects remain rejected even if their bytes match. The native interpreter and original calculation wrapper bytes are unchanged. Nextflow's 15 existing deployment API checks also passed (`cadd-native-runtime-publication-v2-20261018.log`).

Physical prepared SpliceAI2 adoption, 2026-10-18: `spliceai_cache_input.py adopt-prepared --source-dir SOURCE --output-dir DEST --cache-root ROOT` copies the completed two-VCF/index/gene-context bundle into managed storage without conversion. The importer verifies original and copied SHA256 values, resumes only matching owned partial prefixes, and publishes the directory atomically. `adapter.original.ready.json` preserves the original receipt bytes; the new `adapter.ready.json` changes output/context locations and records `relocation.ready.json` provenance while preserving all scores, scientific metadata and the unqualified state. Explicit root directory aliases resolve to canonical receipts; linked child directories, shared staging files and missing original locks fail without changing the source. Maintained source is `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/scripts`; fixture logs/receipts are under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/spliceai2-physical-adoption-20261018-v*`.

Selected plugin inputs and interrupted acquisition, 2026-10-18: `indexed_plugin_cache.py` exposes `ensure-index` and `format-primate`. Index repair preserves table bytes and readable publisher TBI coordinate flags, columns, comment and skipped-header settings, with bounded first/adjacent-coordinate query comparisons. Raw PrimateAI uses the original transformation/sort order and records its exact source SHA256 and declared assembly under the managed plugin cache. UTRAnnotator, SpliceVault, MaveDB, pext and LoFtee honor selected paths and reject missing custom inputs. New UTR/LoFtee clones and LoFtee resources stage before validation/publication; existing LoFtee revisions stay fixed during cache recovery. GRCh37 SpliceVault uses configured `tmp_dir`. Mirrored sources are in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts; closed fixture logs, commands and identity receipts are in /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/plugin-selected-inputs-20261018-v*.

Mapping-bound motif/PM1 upgrades, 2026-10-18: `biomart_nm_to_enst.py artifact-route` accepts named artifact, mapping, anchor-input, cache-root, output-directory, assembly and kind arguments. It verifies the old artifact against its receipt before comparing the existing mapping identity. Changed mapping/provenance/anchor hashes select a deterministic managed namespace; completed matching targets can be reused and registered without rebuilding. Missing/invalid receipts, altered outputs and unsafe linked/occupied destinations remain errors. Old artifacts and receipts are retained; scientific options, calculation order, freshness checks and receipt-write fields are unchanged. Mirrored source is in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts; closed focused logs, commands and identity receipts are in /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/mapping-artifact-routes-20261018-v{1,2}.

<!-- nextflow-main-driver-20261019-start -->
Nextflow driver automation (2026-10-19): `main_install` calls
`Nextflow_runtime_install CONFIG`, which derives `cache_root/runtime/priva-nextflow-25.10.4`.
The helper exposes `install|check --prefix`; existing pins/launcher are unchanged.
152 affected/dispatch tests passed, followed by six dispatch tests with the added
snapshot-failure case;16 helper fixtures use explicitly fake download/JVM assets.
Real cold/warm/check passed in1514363.1291. These checks do not run a workflow or
qualify models. Maintained source is this checkout; closed logs/receipts are under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/nextflow-driver-{integration,real-install}-20261019-v1/`,
`logs/nextflow-bootstrap-unit-20261019-v{1,2,3}.log`, and
`deployment/nextflow-driver-main-integration-real-20261019-v1.json`.
The later root-alias correction passed nine focused cases plus real canonical
and alias warm reuse; it preserves lexical normalization and config bytes.
Closed evidence: `deployment/nextflow-driver-root-alias-20261019-v1.json` and
`validation/nextflow-driver-root-alias-20261019-v1/*.log` under the same data project.
<!-- nextflow-main-driver-20261019-end -->

Original AlphaMissense producer handoff, 2026-10-24 CT: `alphamissense_original_cache.py` exposes `snapshot`, `import` and `check` with named options. It copies an externally attested completed original producer output without loading/recalculating it. The source content hash, selected VCF/index identities, native completion evidence and exact output size/SHA are retained in a separate `.original.complete.json`; retrospective source/input bindings must state their evidence limits. An owned pending journal supports interrupted import; nonowned outputs/receipts, input aliases and cache-root escapes are preserved/refused. `AlphaMissense_stat` checks this receipt before its existing chromosome builder route; invalid/pending original receipts stop. Identical producer content across checkouts reuses the same output; warm checks use unchanged output metadata rather than rehashing the full payload. This is execution/publication proof, not scientific qualification. Mirrored helper: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts/alphamissense_original_cache.py. Closed synthetic fixtures/logs/receipts are under /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/alphamissense-original-handoff-20261024-v{1,2,3}; the original scientific producer and chromosome helper are unchanged by this handoff work.

DAS receipt source portability, 24 October 2026 CT: `das_install_state.py` version3 adds named `--source NAME=PATH`, `--bundled NAME=PATH` and ordered `--parameter NAME=VALUE` roles to its existing `capture`, `check` and `record` procedures. Sources compare content SHA256 across checkouts; dynamic scientific inputs and outputs keep exact path/target/inode/size/mtime identities with node-local device excluded. Actual bundled snapshots use stable logical roles and existence/nonzero checks only; no bundled contents are opened. `AlphaMissense_pick_intolerant_domains` classifies a selected InterPro mapping as bundled only at its genuine canonical bundle origin; external copies and links to external targets remain strict inputs. Thread selection and coldspot/benign filter enablement are recorded in order. The scientific command/options/default thresholds and source calculations are unchanged. Incompatible version2 receipts fail verification and are not promoted. Mirrored helper: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/scripts/das_install_state.py. Closed synthetic logs/results/source snapshots are under /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/das-state-portability-20261024-v1; this is reuse/orchestration evidence, not actual DAS or installation acceptance.
