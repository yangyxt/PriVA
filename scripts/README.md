# PriVA scripts

Updated: 2026-09-14. This directory contains PriVA's installation, annotation,
cache-building, and interpretation scripts. Existing script headers and the
repository README document their interfaces. This inventory records the files
changed for installer startup; other upstream scripts retain their prior roles.

| File | Purpose and change provenance | Status |
| --- | --- | --- |
| install_utils.sh | Resource installer based on commit 810e60165a663bd2e1223809af9f32d664a9d543. Preserves the 2026-09-13 manager discovery and fixed priva_acmg bootstrap. Adds cache_root plus six optional folder names, validates/resolves paths after environment activation, updates saved descendant paths and creates selected directories. A --check entry point previews without writes. CADD follows the selected parent without its legacy storage prompt. | Updated 2026-09-14; 14 focused checks passed on macOS and RHEL 8.10; full resource deployment not performed. |
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
unchanged, and held job `1514589` remains an older frozen submission; updating
this wrapper does not update that job or authorize its release.

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
scored-gene intervals against exact GENCODE48, and invokes `priva spliceai2-vcf`
through mature `gzip`/`bgzip`/`tabix` processes. It preserves originals and publishes
the receipt last; it performs no acquisition or scientific activation.
SpliceAI_install invokes scripts/spliceai_cache_input.py install automatically, retaining the existing prompts and two resource keys. The helper detects gzip by magic bytes, then classifies a bounded decompressed header as vcf_v1, vcf_spliceai2 or the exact native 84-column SpliceAI2 TSV schema; filenames, extensions and score magnitudes do not select the predictor. An unmarked VCF receives the legacy-compatible vcf_v1 routing label; that label alone does not prove predictor identity. Existing compatible VCF pairs pass through, including existing GRCh38 v1; fresh GRCh38 selects v2. Wrong assembly, mixed types or failing v2 preparation fail without fallback. SpliceAI2 VCF passthrough additionally requires a completed adjacent adapter.ready.json with exact SNV/indel/index membership and complete published-catalogue status. This checks completed-publication metadata and recorded membership, not a new full VCF/index-body checksum; it performs no conversion. Reject .partial or .building path components, including symlink targets, and refuse in-progress preparation or staged outputs. A final VCF renamed before receipt publication is not ready. V1 VCF passthrough requires no additional adapter receipt.

```bash
python scripts/spliceai_cache_input.py install --config CONFIG --cache-dir CACHE_DIR --snv SNV --indel INDEL
```

For native TSV, the helper discovers the fixed adjacent SpliceAI2_GRCh38.downloads-complete.json sidecar. The selected SNV/indel inputs must match their manifest roles and the full 47-shard catalogue. The existing low-level SpliceAI2_install preparation/reuse path performs full checksum verification and prepares or reuses a validated receipt, then updates only spliceai_snv_prescore and spliceai_indel_prescore. The pinned priva converter is resolved through PATH; no YAML converter or manifest option is added. Provenance and scientific qualification remain in sidecars/cache headers. The fixed adjacent sidecar is published at yangyxt/caches/SpliceAI2_GRCh38.downloads-complete.json, SHA256 fbc6c9fcd1859d6f1605f9343d5659b8420d7397c4152da19daab34c6194e49d. The PATH entry yangyxt/softwares/bin/priva resolves to the verified optimized 225aa4a converter. This changes no YAML keys beyond the original two paths.

Two-key interface qualification: code 35779c4578509ce01218a052a8e7e1d8a785a75b passed 239 Python tests, zero failures/errors/skips, in 190.848 seconds: 39 helper, 25 installer, 14 runtime-guard, 21 Nextflow-adapter and 140 existing deployment/driver/BGZ tests. The installer-shell cases mock the helper boundary; they are not real conversion evidence. Log: yangyxt_dt/projects/PriVA_upgrade/logs/spliceai-two-key-integration-20261016.log. This supersedes the earlier 197-test interface checkpoint without adding its count; unchanged Rust evidence remains 41 tests and the earlier real-projection evidence remains 512 unique rows.

Actual resource checks used final helper SHA256 f69d6ea8d6179d68ae44cad849b3c6c7912a879f05e305a36acffba09c81bb88: all 47 native TSV and four v1 VCF headers classified correctly, including assembly metadata. Real SpliceAI_install reused the existing v1 pair for hg19/GRCh37 and hg38/GRCh38 with configuration bytes, inode and mtime unchanged and resource stat identities unchanged. Receipt: yangyxt_dt/projects/PriVA_upgrade/validation/spliceai-two-key-real-install-20261016/receipt.json. These are bounded header/reuse checks, not new body checksums, full VEP execution or model/scientific qualification.

Nextflow stub verification passed all four scenarios: raw input, resume with all four tasks cached, preannotated input and stale-environment rejection. Log: yangyxt_dt/projects/PriVA_upgrade/logs/spliceai-two-key-nextflow-20261016.log; retained validation/log artifacts: validation/spliceai-two-key-nextflow-20261016/ under that data project. The run used code 35779c4 plus the documentation state at execution; its recorded source_digest identifies that exact state. Subsequent documentation edits are not attributed to the run.

Frozen-run boundary: job 1514732 continues unchanged with installer f8de4b3 and the previously verified converter. Its prior explicit-manifest/six-key interface is historical execution provenance, not the new installation contract. After successful completion, validated output artifacts may be imported through the two original keys; do not attribute those artifacts or the 197 Python tests to the new dispatcher. The new two-key interface is qualified separately at 35779c4; full scientific acceptance remains pending.

The [route reference](../docs/nextflow_migration.md#spliceai2-preparation) provides
commands, the exact JSON manifest contract, all 47 required shards, verified
GENCODE48 and plugin identities, output patterns and the explicit subset-test
boundary. Downloads are complete (1514724; 47 final shards, zero partials).
All 20,083 real regions passed the stable-ID/GENCODE48 join, retaining full
versions and `_PAR_Y`; ten same-strand symbol groups retain all IDs/intervals
and coalesce labels only per prediction. Opposite-strand symbol contexts fail.
The execution manifest is complete and saved in DataLad commit 5f5946b.
All 47 full source SHA256 checks passed in frozen job 1514732, which has begun
projecting SNV chr1; complete conversion and indexed-VEP qualification remain pending. Publisher expected SHA256 values and VERIFIED acquisition receipts for
all 47 shards are now supplied in [caches/README.md](../../../caches/README.md);
there is no missing-publisher-checksum blocker. Current Rust 225aa4a passed 32
projection + 9 existing tests; Python f8de4b3 passed 197 tests with zero failures/skips.
A separate 512-row real prefix check (chr1 positive strand only) preserved every
projected record/order and passed bgzip/tabix queries; this is not full-cache or
full-VEP qualification. Optimized source 225aa4aa774b4be7f53576e1cde56af71e67bd2c
is installed at `yangyxt/softwares/priva-tools-225aa4a/bin/priva`, SHA256
`bbc3b7288ed2ea63f1f8455ed87b8a9d08c7bc1fe0a562c1d6134525e2a12dc6`.
Installed-release reruns passed the same 25 driver cases and same 512 source rows;
the frozen checkpoint had 197 Python/41 Rust tests and 512 rows; current interface coverage is 239 Python tests, with unchanged 41 Rust tests and 512 unique rows. A stale shared Cargo
release artifact was detected before deployment; forced rebuilding and installed
smoke checks corrected it, with the failure retained in the receipt. Full preparation job 1514732 is running on c81o-03 (4 CPUs, 16 GB, 7 days) using frozen
installer f8de4b3, this converter and the isolated GRCh38 config. The manifest/
handoff were saved in DataLad commit 5f5946b before launch; all 47 full source SHA256 checks passed and SNV chr1 projection started; complete conversion remains pending.
Original licensed source and derived views remain local under their resource terms.
