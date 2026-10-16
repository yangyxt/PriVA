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
ACMG parity. The initial queued job is 1514601 after source job 1514584, using
installer `054b4944055054148005016192647a8ed0384eef`.

`build_cluster_splicevault.sbatch` runs the original GRCh38-to-GRCh37 conversion
against the release113 SpliceVault input after checking its publisher checksum
and byte size. Job 1514603 freezes the same installer revision, uses a separate
config, and records reference/chain/source/output identities on completion.
It reserves eight CPUs for the original sort command; the newer publisher hg19
table is kept separate from this baseline conversion.

`import_cluster_licensed.sbatch` registers explicitly confirmed acquisitions under
`caches/{SpliceAI_v1,PrimateAI_v0_2}`. The optional ninth argument is the resource
selection: `all` (the original default), `spliceai`, or `primateai`. It still runs
the frozen `install_utils.sh` functions; selecting SpliceAI calls only
`SpliceAI_install` and requires both SNV and indel VCFs with their indexes.
SpliceAI v1 downloads are now owner-confirmed complete; PrimateAI acquisition is
still pending. This script does not download or authenticate.

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
