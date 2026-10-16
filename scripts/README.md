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

`build_cluster_alpha.sbatch` runs the original AlphaMissense cache stages after
an explicit scheduler dependency on the VEP cache job. `diagnose_cluster_install.sbatch`
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
