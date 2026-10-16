# PriVA

A comprehensive pipeline for ACMG/AMP-based variant interpretation and variant prioritization in genetic analysis.

## Nextflow / Rust migration branch

This branch develops a Nextflow workflow and one opt-in `priva` Rust/Polars
binary, with faithful function order and exact comparison to source `09223ef`.
See [execution and deployment](docs/nextflow_migration.md) and
[Rust scope and exact comparison](docs/rust-migration.md). The
[frozen ranking ledger](docs/ranking-source-ledger-09223ef.md) maps remaining
control/ranking mutations, SciPy dependencies and parser/failure contracts.

```text
One official VEP run → PriVA annotation view → family interpretation / Rust ports
                    ↘ CADD feature enrichment → complete annotation artifact
                       → pinned standalone cadd-rs → original Python raw + PHRED
```

`main.nf` isolates existing PriVA tasks; `cadd.nf` develops independent annotation
and calculation stages. Shared-VEP transcript selection and complete feature
mapping must be qualified before production activation. CADD retains its original
Python calculation, scientific dependencies and models. Profile first; enable
Python optimizations only after exact score/row/order/missingness/failure checks
and measured timing/memory comparisons against the pinned unmodified reference.

Rust currently covers final ordering and a bounded criteria-summary/posterior
interface, not all criterion assignment. Cluster software installation has passed;
cache deployment and calibration are still in progress. Synthetic checks, real
CADD model tests, public benchmark acceptance and family acceptance are distinct.
The public benchmark and supplied-family exact ACMG/ranking gates remain pending.

The actual installer reproduced and fixed VEP's unattended false success; the
corrected run reached the pinned cache download. Full cache deployment continues.
Both 18-row CADD assembly fixtures now match original Python exactly through
Nextflow. One `NaN` edge still changes the wrapper exit status when original
Python emits a malformed table; this mismatch remains reported and unaccepted.

The separate experimental CADD batching branch has passed 44 exact original
transformation comparisons and 3 helper cancellation cases. Completed Slurm
benchmark 1514607 passed all 24 paired synthetic trials, including 96 exact
intermediate/final-table comparisons. At 100,008 rows, median calculation time
fell from 556–578 seconds to about 279 seconds, with similar sampled summed
process-tree RSS (this estimate can count shared pages repeatedly). These results
do not qualify shared annotation, failure behavior, public data or families.
The tested integration remains pinned to standalone `4194be3`.
The later wrapper audit failed PID-only SIGTERM cleanup for the experimental
shell, profiler and standalone CLI; process-group termination passed and the
audit reaped all remaining children. The candidate remains disabled.

SpliceAI v1 originals/indexes and both PrimateAI v0.2 originals are acquired.
The hg19 SpliceAI import 1514615 and PrimateAI import 1514620 completed with
exit 0 and valid receipts. The PrimateAI receipt includes exact external
source-to-formatted validation; cancelled job 1514618 never ran. Original hg38
files remain retained and unselected. Separate main-config promotion updated only
the three selected hg19 resource-path keys with a receipt; scientific qualification
remains pending. gnomAD job 1514601 failed its `.vcf.bgz` extension validation;
replacement 1514729 is running with the minimal fix (eight focused tests passed).
AlphaMissense 1514589
remains held for the original NM-to-ENST mapping. No evidence is silently omitted.
The [installation dependency ledger](docs/installation-dependency-ledger-09223ef.md)
records the ordered calibration recipe and separates bundled ClinVar baseline
statistics from its currently unqualified fresh-build path.

## Overview

This repository contains scripts and configuration files for an automated variant prioritization pipeline based on the ACMG/AMP guidelines. PriVA combines transcript-level annotation, population and functional data, gene-disease mechanism and inheritance context, and family information to evaluate ACMG/AMP criteria and produce ranked variant lists for each family. Annotation resources are deployed once and can be refreshed; each cohort run then performs annotation, family-level filtering, criteria matching, and prioritization.

![PriVA workflow (Figure 1)](PriVA_figure1.png)

## Requirements

- Micromamba, Mamba, or Conda for environment management
- Sufficient disk space for annotation resources (>1TB)
- Reference genome file

## Pedigree File Requirements

If a valid pedigree table is provided, please note:
1. The family sample ID should be consistent between the pedigree table and the VCF file.
2. The proband of the family should stay as the first record of the family in the pedigree table
3. Every family should have a unique identifier in the pedigree table

## SpliceAI2 GRCh38 preparation

Current installation direction: keep only the original spliceai_snv_prescore and spliceai_indel_prescore YAML keys. Choose assembly once at installation: fresh GRCh37 uses v1 and fresh GRCh38 uses v2. Inspect decompressed headers automatically: an existing compatible VCF passes through; native 84-column SpliceAI2 TSV invokes the narrow adapter. A prepared pair is registered under the two original keys and installation succeeds. Resource/model identity, hashes, schema and qualification remain in sidecars and cache headers, without extra selector or converter YAML switches. Existing compatible GRCh38 v1 VCF pairs pass through; fresh GRCh38 installations select v2. Wrong assembly, mixed input types and failed v2 preparation never trigger a fallback.

The [preparation route](docs/nextflow_migration.md#spliceai2-preparation) keeps the
original 84-column tables and creates a separate top0 VCF for the official plugin,
using exact scored-region/GENCODE48 gene context. The Rust converter streams rows;
the Python deployment driver verifies the 47-shard manifest and controls
`gzip`/`bgzip`/`tabix`. It downloads nothing. GENCODE48 has passed its resource
checks. Download job 1514724 completed (47 final shards, zero partials), and the
complete 20,083-region stable-ID join passed with exact bounds/strand/biotype.
Full GTF IDs/versions and `_PAR_Y` remain in context; ten same-strand symbol
groups coalesce labels only per prediction, while opposite-strand symbol contexts
fail. The execution manifest is complete and saved in DataLad commit 5f5946b.
All 47 full source SHA256 checks passed in frozen job 1514732, which has started
projecting SNV chr1; complete conversion and indexed-VEP checks remain pending.

```text
Assembly → decompressed headers → VCF passthrough / native84 adapter
                                → indexed VCF pair / two original YAML keys
                                → sidecars + cache-header scientific gate
```

The simplified installer succeeds after preparing/registering the indexed pair
under the two original resource keys; provenance stays in sidecars.

SpliceAI_install invokes scripts/spliceai_cache_input.py install automatically, retaining the existing prompts and two resource keys. The helper detects gzip by magic bytes, then classifies a bounded decompressed header as vcf_v1, vcf_spliceai2 or the exact native 84-column SpliceAI2 TSV schema; filenames, extensions and score magnitudes do not select the predictor. An unmarked VCF receives the legacy-compatible vcf_v1 routing label; that label alone does not prove predictor identity. Existing compatible VCF pairs pass through, including existing GRCh38 v1; fresh GRCh38 selects v2. Wrong assembly, mixed types or failing v2 preparation fail without fallback. SpliceAI2 VCF passthrough additionally requires a completed adjacent adapter.ready.json with exact SNV/indel/index membership and complete published-catalogue status; it performs no conversion. Reject .partial or .building path components, including symlink targets, and refuse in-progress preparation or staged outputs. A final VCF renamed before receipt publication is not ready. V1 VCF passthrough requires no additional adapter receipt.

```bash
python scripts/spliceai_cache_input.py install --config CONFIG --cache-dir CACHE_DIR --snv SNV --indel INDEL
```

For native TSV, the helper discovers the fixed adjacent SpliceAI2_GRCh38.downloads-complete.json sidecar. The selected SNV/indel inputs must match their manifest roles and the full 47-shard catalogue. It calls the existing low-level SpliceAI2_install to prepare or reuse a fully validated receipt, then updates only spliceai_snv_prescore and spliceai_indel_prescore. The pinned priva converter is resolved through PATH; no YAML converter or manifest option is added. Provenance and scientific qualification remain in sidecars/cache headers. The fixed adjacent sidecar is published at yangyxt/caches/SpliceAI2_GRCh38.downloads-complete.json, SHA256 fbc6c9fcd1859d6f1605f9343d5659b8420d7397c4152da19daab34c6194e49d. The PATH entry yangyxt/softwares/bin/priva resolves to the verified optimized 225aa4a converter. This changes no YAML keys beyond the original two paths.

Two-key interface qualification: code 35779c4578509ce01218a052a8e7e1d8a785a75b passed 239 Python tests, zero failures/errors/skips, in 190.848 seconds: 39 helper, 25 installer, 14 runtime-guard, 21 Nextflow-adapter and 140 existing deployment/driver/BGZ tests. The installer-shell cases mock the helper boundary; they are not real conversion evidence. Log: yangyxt_dt/projects/PriVA_upgrade/logs/spliceai-two-key-integration-20261016.log. This supersedes the earlier 197-test interface checkpoint without adding its count; unchanged Rust evidence remains 41 tests and the earlier real-projection evidence remains 512 unique rows.

Earlier remote CI at eb16493 passed [Actions 38852306750](https://github.com/yangyxt/PriVA/actions/runs/38852306750). The new 35779c4 source has passed local checks; its push/remote CI are not yet recorded at this checkpoint.

Actual resource checks used final helper SHA256 f69d6ea8d6179d68ae44cad849b3c6c7912a879f05e305a36acffba09c81bb88: all 47 native TSV and four v1 VCF headers classified correctly, including assembly metadata. Real SpliceAI_install reused the existing v1 pair for hg19/GRCh37 and hg38/GRCh38 with configuration bytes, inode and mtime unchanged and resource stat identities unchanged. Receipt: yangyxt_dt/projects/PriVA_upgrade/validation/spliceai-two-key-real-install-20261016/receipt.json. These are bounded header/reuse checks, not new body checksums, full VEP execution or model/scientific qualification.

Nextflow stub verification passed all four scenarios: raw input, resume with all four tasks cached, preannotated input and stale-environment rejection. Log: yangyxt_dt/projects/PriVA_upgrade/logs/spliceai-two-key-nextflow-20261016.log; retained validation/log artifacts: validation/spliceai-two-key-nextflow-20261016/ under that data project. The run used code 35779c4 plus the documentation state at execution; its recorded source_digest identifies that exact state. Subsequent documentation edits are not attributed to the run.

Frozen-run boundary: job 1514732 continues unchanged with installer f8de4b3 and the previously verified converter. Its prior explicit-manifest/six-key interface is historical execution provenance, not the new installation contract. After successful completion, validated output artifacts may be imported through the two original keys; do not attribute those artifacts or the 197 Python tests to the new dispatcher. The new two-key interface is qualified separately at 35779c4; full scientific acceptance remains pending.
GRCh37 v1 and legacy CADD masked-v1 inputs remain unchanged. The real-plugin
boundary harness passed 11 subtests / 38 assertions with synthetic tabix/VEP
objects. Historical checks before the stable-ID/same-strand-symbol update: Rust
passed 30 new plus 9 existing tests; the integrated installer/driver run passed 165 tests with zero failures/skips, including real Rust/bgzip/tabix
on synthetic inputs. These historical totals remain tied to their producing source.
Current Rust 225aa4a passed 32 projection + 9 existing tests; Python f8de4b3 passed
197 tests with zero failures/skips (169 installer/deployment + 20 adapter + 8 BGZ).
A bounded real-cache check passed 256 SNV + 256 indel chr1 positive-strand rows
with exact full-record/order and bgzip/tabix checks. It does not establish model
inference, whole-cache, full VEP, public-sample or family acceptance.

All 47 publisher SHA256 values and VERIFIED acquisition receipts are now supplied
in [caches/README.md](../../caches/README.md). The complete execution manifest and handoff were saved in DataLad commit 5f5946b before launch. Full preparation job
1514732 is running on c81o-03 (4 CPUs, 16 GB, 7 days), using frozen installer f8de4b3 and
the isolated GRCh38 preparation config; all 47 full source SHA256 checks passed and SNV chr1 projection started; complete conversion is pending. Optimized
source 225aa4aa774b4be7f53576e1cde56af71e67bd2c is installed at
`yangyxt/softwares/priva-tools-225aa4a/bin/priva`, SHA256
`bbc3b7288ed2ea63f1f8455ed87b8a9d08c7bc1fe0a562c1d6134525e2a12dc6`.
The installed release passed the same 25 driver cases and same 512 real-prefix
rows again with zero mismatches; the frozen checkpoint had 197 Python, 41 Rust and 512 rows; the new interface
has 239 Python tests, with unchanged 41 Rust tests and 512 unique rows. A stale shared Cargo artifact was caught before deployment, rebuilt and
checked at the installed path; the receipt retains that failed attempt. Full
conversion and real VEP/scientific acceptance remain pending. See the linked route
for the exact manifest, GTF/plugin identities, commands and pending gates.

Scientific boundary for the simplified installer: inspect the actual selected cache headers and sidecar provenance to prevent unqualified v2 scientific use, including preannotated import and resume; no new enable/disable YAML toggle is introduced. Successful installation is separate from v2 ACMG/ClinVar, full VEP, public-sample and supplied-family acceptance. Unmarked external annotations still require pinned scientific provenance because predictor version cannot be inferred from numeric values. The runtime guards now inspect actual cache headers at the two original keys, including preannotated launch/task paths; the current 21 Nextflow-adapter and 14 runtime-guard tests are included once in the 239-test total. Earlier f8de4b3 checks remain historical; scientific acceptance is still pending.

## Installation

The github repository includes an one-stop installation BASH script (`scripts/install_utils.sh`) that helps deploy all required dependencies on your host system, with the help from a configuration YAML file (e.g., `config.yaml`). This YAML file is specifying the paths to external annotation dependencies.

**Carefully read comments in `config.yaml` to prepare for execution of install_utils.sh and PriVA**

`main_install` detects an available Micromamba, Mamba, or Conda executable,
then checks for the fixed environment `priva_acmg`. If absent, it creates that
environment from the repository's `acmg_conda_rough.yml`; if present, it reuses
it. The installer activates and checks the environment before reading the
project configuration with its `yq` and continuing resource installation.
The bundled environment includes `yq`, so it is not needed beforehand.

```bash
bash scripts/install_utils.sh main_install path/to/config.yaml
```

The environment file and name are built into the installer; the only argument
is the project configuration, which specifies resource paths and settings.
The legacy `conda_env_yaml` configuration value is ignored. Existing environments
are checked for required tools and VEP packages; an incomplete environment stops
installation with an error. Miniconda installations are detected through their
`conda` executable.

Executable discovery checks `MAMBA_EXE`, then `micromamba`/`mamba` on `PATH`,
then `CONDA_EXE`/`conda`. It retains the selected executable or wrapper path.
Conda-backed Mamba installations use their sibling Conda for activation and saved
environment variables. Keep the manager's root, package cache, and environment
directories configured to the intended storage location. The installer does not
initialize shell profiles or choose a new manager root.

To prepare only the environment, the existing entry point remains available:

```bash
bash scripts/install_utils.sh conda_install_vep
```

Activation inside a child installer does not activate the caller's shell. Enter
the environment yourself before invoking individual resource functions directly.
The first resource installation remains assembly-specific; configure and prepare
its reference and other settings before running the full installer. With
`cache_root`, the six selected resource directories are created automatically.

Startup regression tests use isolated fake managers and do not install packages
or download resources:

```bash
python -m unittest discover -s tests -p test_install_utils_bootstrap.py -v
```

Installer startup updated 2026-09-13. Implementation and helper metadata are in
`scripts/README.md`; regression-test metadata are in `tests/README.md`.

`main_install` also deploys the gene-mechanism resources used by ACMG step 3 in
dependency order: pinned HPO/MONDO data, the broader LOF/non-LOF mechanism
cache, exact normalized GoFCards variants, the canonical non-LOF/ClinVar VCV
cache, and the integrated HPO-framed condition cache. All builders, parsers,
schemas, and outputs are inside the PriVA tree; `llm_gene_reranker` is not an
installation dependency.

The complete canonical JSON and ClinVar VCV/RCV/SCV field contract is documented
in `GENE_NONLOF_MECHANISM_JSON_STRUCTURE.md` and enforced by the bundled JSON
Schema.

The supported full manual refresh is the same entry point as installation:

```bash
bash scripts/install_utils.sh main_install path/to/config.yaml
```

Users who want automatic refreshes should schedule that command with their own
cron or cluster scheduler. PriVA does not install a separate weekly job. Each
run checks source-specific refresh intervals and rebuilds downstream caches only
when an input checksum or implementation file changes. A forced mechanism
refresh can be requested with `PRIVA_FORCE_ALL_CACHES=1`.

The mechanism stages can also be run independently for debugging:

```bash
bash scripts/install_utils.sh gene_pathogenic_mechanism_cache_install path/to/config.yaml
bash scripts/install_utils.sh gofcards_exact_gof_cache_install path/to/config.yaml
bash scripts/install_utils.sh gene_nonlof_mechanism_cache_install path/to/config.yaml sources-only
bash scripts/install_utils.sh gofcards_clinvar_injection_install path/to/config.yaml
bash scripts/install_utils.sh gene_nonlof_mechanism_cache_install path/to/config.yaml
bash scripts/install_utils.sh mondo_hpo_scope_install path/to/config.yaml
bash scripts/install_utils.sh hpo_condition_mechanism_cache_install path/to/config.yaml
bash scripts/install_utils.sh mechanism_resource_install path/to/config.yaml
```

The GoFCards steps must run in the order shown. `gofcards_exact_gof_cache_install`
writes only a build intermediate; `gofcards_clinvar_injection_install` turns that
into the one deployed cache every consumer reads. The `sources-only` call in
between downloads the weekly ClinVar VCV XML that injection needs, because the
component that owns that download is the non-LOF builder, and the non-LOF
builder's ordinary run cannot start until injection has finished. See
`docs/GOFCARDS_CACHE_SCHEMA.md` for the full ordering.

The integrated JSON is published by writing and validating a complete temporary
snapshot before atomically replacing the live cache. It nests external
mechanism histories and variants under an HPO gene/condition only through exact
OMIM, ORPHA, or MONDO identifiers; unresolved records remain under the gene's
`unmapped_evidence` section.

The MONDO/HPO installer downloads pinned releases, builds one disease-scope
registry, and adds the resulting scope and provenance fields to each HPO
assertion. Confirmed non-monogenic diseases are excluded from PriVA's flattened
inheritance/penetrance context. Neoplastic diseases without positive hereditary
evidence remain review-required rather than being assumed somatic-only. Manual
decisions are stored in `data/mondo/disease_scope_overrides.tsv` and take
precedence when the registry is rebuilt.

## Smoke test

The smoke test is a small end-to-end demo for checking that PriVA can run on the current machine. It uses curated demo assets in `example/smoke_data/` and does not run the full production annotation-resource installation.

From a fresh clone, first create and activate the PriVA conda environment:

```bash
git clone <PriVA_repo_url>
cd PriVA
mamba env create -f acmg_conda.yml
mamba activate priva_acmg
```

If `mamba` is unavailable, use `conda env create -f acmg_conda.yml` instead.

Then run the smoke test for one or both supported assemblies:

```bash
./example/run_smoke_test.sh hg19
./example/run_smoke_test.sh hg38
```

The wrapper checks required command-line tools and bundled smoke resources, downloads the configured UCSC reference FASTA if needed, creates the FASTA index with `samtools faidx`, and runs the PriVA Snakemake workflow on the 20-variant demo set. Expected outputs are written under `example/smoke_runs/`.

The smoke test does not require Git LFS on current `main`; all current smoke assets are stored as normal Git files.

## Configuration File
The `config.yaml` file in the PriVA repo contains all the parameters you need. Most of them can be automatically handled by the BASH script `install_utils.sh`. Here we list the items need to be specified by users prior to executing `install_utils.sh`

- TEMPORARY DIRECTORY: `tmp_dir`, used to store the temp folder and file.
- REF GENOME ASSEMBLY: `assembly` (e.g., hg19, hg38), `ref_genome`, path to the corresponding FASTA file.
- ANNOTATION CACHE ROOT: `cache_root`, an absolute storage path for the six resource groups below.

```yaml
cache_root: /stornext/snfs170/next-gen/scratch/yangyxt/caches

# Optional; omit the mapping or individual entries to use these defaults.
cache_folders:
  vep_cache: VEP_cache
  vep_plugin: VEP_plugin
  vep_plugin_cache: VEP_plugin_cache
  gnomad: gnomAD
  clinvar: ClinVar
  cadd: CADD
```

Each directory is `cache_root` plus its folder name. Use single folder names
containing letters, digits, underscores, dots or hyphens, starting with a letter,
digit or underscore. The root must be an absolute path without spaces or shell
characters, as required by the existing resource commands. The folder mapping
must be indented as shown; unknown entries and nested/absolute folder paths fail.

After environment activation, the installer resolves the six directories, creates
them as needed and saves `vep_cache_dir`, `vep_plugins_dir`,
`vep_plugins_cachedir`, `gnomad_vcf_dir`, `clinvar_vcf_dir` and `cadd_parent_dir`
for the existing resource and annotation functions. Treat those full paths as
installer-maintained values when `cache_root` is set. Saved resource file paths
beneath the previous directories are updated to the selected root/folder names;
explicit resource paths outside those directories remain unchanged. Files are
not moved or copied. Rerun installation after changing the root or folder names
to prepare/reuse resources at the new destinations before running annotation.

CADD uses the selected `CADD` folder without asking to store data in the package
manager's environment. Its required `CADD-scripts-<version>/data/` structure is
preserved beneath that folder. SpliceAI/PrimateAI manual-download prompts remain.
Reference FASTA, temporary storage and resource groups outside these six keep
their existing settings. Configurations without `cache_root` retain legacy full
paths and the previous CADD storage prompt.

After activating the environment, preview resolution without changing files or
creating folders, or prepare only the folder configuration before invoking an
individual resource function:

```bash
bash scripts/install_utils.sh priva_configure_cache_paths config.yaml --check
bash scripts/install_utils.sh priva_configure_cache_paths config.yaml
```

Updated 2026-09-14: `config.yaml` documents the cache root and optional folder
names and contains the corresponding HGSC resource paths. `scripts/install_utils.sh`
resolves those paths and makes CADD follow the selected directory; the environment
startup work from 2026-09-13 is preserved. No annotation resources were installed
as part of this code change. Additional metadata is in `scripts/README.md` and
`tests/README.md`.

**Special annotation resources need manual download**
There are 2 annotation resource (SpliceAI, PrimateAI) need manual download from Illumina
`primateai_prescore`: prompt will present during the execution of `install_utils.sh`, please just follow the prompt to download the file
`spliceai_snv_prescore` & `spliceai_indel_prescore`: each one takes space up to several hundreds GB, so please follow the prompt to download the file


There are also items in `config.yaml` that needs to be specified for each PriVA exectuion:

- INPUT_VCF: `input_vcf`, which supposed to be a cohort VCF file
- PEDIGREE_FILE: `ped_file`, which is the pedigree file corresponding to the cohort VCF file. 
- OUTPUT_DIR: `output_dir`, which is the directory containing all the output annotation files
- HUB_CACHE_FILEs: `hub_file_dir`, parent folder storing the hub cache files; `hub_vcf_file`, hub annotated VCF file for annotation caching; `hub_cadd_file`, hub annotated CADD TSV for CADD annotation caching
- ALLELE_FREQUENCY_THRESHOLDS: `extreme_rare_PAF`, the extreme rare population allele frequency threshold for PM2; `af_cutoff`, population allele frequency threshold distinguishing common and rare variants
- DISEASE_INCIDENCE: `exp_disease_incidence`, the expected incidence of the disease
- ALT_DISEASE_PAT_VARs: `alt_disease_vcf`, Abs path to the VCF file containing the variants from alternative disease patients, used for BP5
- CONTROL_VARs: `control_vcf`, abs path to the VCF containing variants from control samples
- PP1_RESOURCES: `pp1_vcf`, abs path to the VCF containing patients (as well as family members) from the same disease, for family consegregation; `pp1_ped`, pedigree file corresponding to `pp1_vcf`
- COMPUTATION_RESOURCES: `threads`, total CPUs assigned to PriVA execution; `threads_per_fam`, per fam processing (step2 and step3) are executed in parallel across all families, this specify how many CPUs to use per family.


## Running PriVA

The recommended way to run the pipeline is using Snakemake, which manages dependencies and execution order.

**Prerequisites:**

1. Activate Conda Environment: Ensure the `priva_acmg` conda environment is activated:
   ```bash
   conda activate priva_acmg
   ```
2. Configuration File: Have your `config.yaml` file ready with all necessary paths and parameters correctly specified.

**Running the Pipeline:**

Navigate to your working directory (or specify absolute paths in the command).
The config file can be duplicated and customized for batch-specific runs.

**Basic Command:**

```bash
snakemake --snakefile /path/to/your/PriVA/Snakefile \
          --configfile /path/to/your/config.yaml \
          --cores <number_of_cores>
```

* `--snakefile`: Path to the pipeline's `Snakefile`.
* `--configfile`: Path to your specific configuration file.
* `--cores` (or `-j`): Maximum number of CPU cores Snakemake can use.

**Bash Function Wrapper Example:**

You can wrap the Snakemake command in a bash function for easier execution. Add this function to your `.bashrc` or `.bash_profile`, or save it in a script file and source it.

```bash
run_priva_pipeline() {
    # --- Configuration ---
    local snakefile_path="/path/to/your/PriVA/Snakefile" # CHANGE THIS to the actual path
    local default_config="/path/to/your/default_config.yaml" # Optional: Set a default config
    local default_cores=8
    local default_log_dir="snakemake_logs"

    # --- Argument Parsing ---
    local config_file="${1:-$default_config}"
    local num_cores="${2:-$default_cores}"
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local log_dir="${3:-$default_log_dir}"
    local log_file="${log_dir}/priva_run_${timestamp}.log"

    # --- Input Validation ---
    if [[ ! -f "$snakefile_path" ]]; then
        echo "ERROR: Snakefile not found at $snakefile_path"
        return 1
    fi
    if [[ -z "$config_file" ]] || [[ ! -f "$config_file" ]]; then
        echo "ERROR: Config file not found or not specified: '$config_file'"
        echo "Usage: run_priva_pipeline <config_file_path> [num_cores] [log_directory]"
        return 1
    fi
    if ! [[ "$num_cores" =~ ^[0-9]+$ ]] || [[ "$num_cores" -lt 1 ]]; then
        echo "ERROR: Number of cores must be a positive integer: '$num_cores'"
        return 1
    fi

    # --- Environment Check ---
    if [[ -z "$CONDA_DEFAULT_ENV" ]] || [[ "$CONDA_DEFAULT_ENV" != "priva_acmg" ]]; then
        echo "WARNING: Conda environment 'priva_acmg' might not be active. Attempting to activate..."
        conda activate priva_acmg
        if [[ $? -ne 0 ]]; then
            echo "ERROR: Failed to activate conda environment 'priva_acmg'. Please activate it manually."
            return 1
        fi
        echo "Activated conda environment 'priva_acmg'."
    fi

    # --- Execution ---
    echo "Starting PriVA pipeline..."
    echo "Snakefile: $snakefile_path"
    echo "Config File: $config_file"
    echo "Cores: $num_cores"
    mkdir -p "$log_dir"
    echo "Log File: $log_file"
    echo "-----------------------------------------"

    nohup snakemake --snakefile "$snakefile_path" \
                    --configfile "$config_file" \
                    --cores "$num_cores" \
                    --printshellcmds \
                    --rerun-incomplete \
                    --verbose \
                    > "$log_file" 2>&1 &

    local job_pid=$!
    echo "Pipeline submitted in the background. PID: $job_pid"
    echo "Monitor progress with: tail -f $log_file"
    echo "-----------------------------------------"
}

# --- How to use the function ---
# 1. Source the script or add the function to your .bashrc/.bash_profile
# 2. Run it:
#    run_priva_pipeline /path/to/my_experiment_config.yaml 16 /path/to/my/logs
#    # Or using defaults:
#    run_priva_pipeline /path/to/my_experiment_config.yaml
```


**Explanation of the Function:**

1. **Configuration:** Set the path to your `Snakefile` and optionally a default config file path.
2. **Arguments:** Takes the config file path (required), number of cores (optional, defaults to 8), and log directory (optional, defaults to `snakemake_logs`) as arguments.
3. **Validation:** Checks if the Snakefile and config file exist and if the core count is valid.
4. **Environment Check:** Verifies if the `priva_acmg` conda environment is active and tries to activate it if not.
5. **Execution:** Creates a log directory, constructs the `snakemake` command with useful options (`--printshellcmds`, `--verbose`, `--rerun-incomplete`), and runs it in the background using `nohup`. It prints the PID and log file location for monitoring.

**Comprehensive Example Command (from original README):**

This shows a full command often used for running on a cluster or server:

```bash
# Activate environment first!
# conda activate priva_acmg

nohup snakemake --snakefile /path/to/PriVA/Snakefile \
                --cores 50 \
                --configfile /path/to/config.yaml \
                --printshellcmds \
                --verbose \
                > /path/to/snakemake_pipeline.log 2>&1 &
```

* `nohup ... &`: Runs the command in the background, detached from the terminal, logging output to the specified file.
* `--printshellcmds`: Shows the actual shell commands being executed by Snakemake.
* `--verbose`: Provides more detailed output from Snakemake.

For family-specific analysis, ensure your pedigree file (`ped_file` in `config.yaml`) is properly formatted and specified in the config file.


## Annotation Resources

The pipeline integrates multiple annotation resources to comprehensively evaluate variant pathogenicity:

1. **Splicing disruption effect**: 
   - SpliceAI (VEP_plugin)
   - SpliceVault (VEP_Plugin)

2. **5-UTR uORF disruption effect**: 
   - UTRannotator (VEP_plugin)

3. **General Deleterious effect prediction**: 
   - CADD (standalone client)
   - PrimateAI (VEP_plugin)
   - AlphaMissense (VEP_plugin)
   - LOFTEE (VEP_plugin)

4. **Haplo-insufficiency**: 
   - LOEUF (VEP_plugin)
   - AlphaMissense_mean_score_per_gene

5. **Transcript disruption effect prediction**: 
   - VEP

6. **Clinical variants**: 
   - ClinVar VCF (bcftools annotate)

7. **Population-wise allele frequency + number of homozygous carriers**: 
   - gnomADv4 VCFs (bcftools annotate)

8. **Conservation**: 
   - Conservation (VEP_plugin)


## Pipeline Structure

The pipeline consists of three main steps:

1. **Annotation**: Annotates variants with multiple resources (Rule `annotate_variants`)
   - VEP annotation with plugins
   - CADD scoring
   - gnomAD frequency annotation
   - ClinVar annotation

2. **Filtration**: Filters variants based on family information (if available) (Rules `filter_variants_per_family` / `filter_variants_no_family`)
   - Family-specific filtering when pedigree is provided
   - General filtering without family information

3. **Prioritization**: Prioritizes variants according to ACMG guidelines (Rules `prioritize_variants_per_family` / `prioritize_variants_no_family`)
   - Generates TSV output with prioritized variants
   - Produces ACMG-specific reports


## Output Files

The pipeline produces several output files in the directory specified by `output_dir` in your `config.yaml`:

1. Annotated VCF files (`<input_base>.anno.vcf.gz`)
2. CADD score files (`<input_base>.anno.cadd.tsv`)
3. Filtered VCF files (`<input_base>.anno.filtered.vcf.gz` or `<input_base>.anno.<family>.filtered.vcf.gz`)
4. Prioritized variant tables (`<input_base>.anno.filtered.tsv` or `<input_base>.anno.<family>.filtered.tsv`)
5. ACMG-classified variant tables (`<input_base>.anno.filtered.acmg.tsv` or `<input_base>.anno.<family>.filtered.acmg.tsv`)

Where `<input_base>` is derived from your `input_vcf` filename, and `<family>` corresponds to family IDs found in your `ped_file`.

## Workflow Modes

The pipeline automatically adapts based on the presence of the `ped_file` specified in the `config.yaml`:

- **With family information**: If `ped_file` exists, uses pedigree data to perform family-specific filtering and prioritization for each family listed.
- **Without family information**: If `ped_file` does not exist or is not specified, performs general filtering and prioritization on the entire VCF.

## License

PriVA is licensed under the Apache License 2.0 with Commons Clause. This means:
- You can use, modify, and share it freely for academic and clinical research.
- Commercial use is prohibited without my permission, including selling PriVA as a product or involving PriVA as a part of charging service.

## Contact

Maintainer: yangyxt@hku.hk, yangyxt@gmail.com
