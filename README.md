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
transformation comparisons and 3 helper cancellation cases. Full scoring and
timing/memory benchmarks continue under Slurm; production remains pinned to
standalone `4194be3`. These helper checks do not qualify complete scoring,
outer-wrapper cancellation, shared annotation or family interpretation.

SpliceAI v1 originals and indexes for both assemblies are confirmed downloaded.
The current hg19 pair is undergoing full verification and original-installer
registration in Slurm 1514615. Both PrimateAI v0.2 originals are now acquired
under `yangyxt/caches/PrimateAI_v0_2/` and passed acquisition gzip/header checks.
The hg19 import is pending: job 1514618 was held before execution after a fixture
showed the original formatter could register a truncated result following a
child-sort failure. An independent exact source-to-formatted comparison is being
added before any import receipt or production registration. The hg38 original
remains available without changing the selected scientific profile.
The AlphaMissense build is held for a separate missing original NM-to-ENST
mapping table. No annotation or evidence branch has been silently disabled.
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
