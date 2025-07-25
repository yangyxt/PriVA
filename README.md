# PriVA

A comprehensive pipeline for ACMG-based variant prioritization in genetic analysis.

## Overview

This repository contains scripts and configuration files for running an automated variant prioritization pipeline based on ACMG guidelines. The pipeline integrates multiple annotation resources to evaluate variant pathogenicity and prioritize variants for clinical interpretation.

![PriVA Workflow](PriVA_figure1.png)

## Requirements

- Snakemake workflow management system
- Conda/Mamba for environment management
- Sufficient disk space for annotation resources (~1TB)
- Reference genome files

## Pedigree File Requirements

If a valid pedigree table is provided, please note:
1. The family sample ID should be consistent between the pedigree table and the VCF file.
2. The proband of the family should stay as the first record of the family in the pedigree table
3. Every family should have a unique identifier in the pedigree table

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

## Installation

The pipeline includes an installation utility (`scripts/install_utils.sh`) that helps set up all required components:

1. Create and configure the conda environment (expected name `priva_acmg`):
   ```bash
   bash scripts/install_utils.sh conda_install_vep path/to/env.yaml
   # Make sure the environment name in env.yaml is 'priva_acmg' or activate it manually later
   ```

2. Install VEP and its plugins:
   ```bash
   bash scripts/install_utils.sh vep_cache_api_install --VEP_CACHEDIR /path/to/cache --VEP_ASSEMBLY GRCh37/GRCh38
   ```

3. Install annotation resources:
   ```bash
   bash scripts/install_utils.sh main_install path/to/config.yaml
   ```

## Configuration

The pipeline is configured through a YAML file (e.g., `config.yaml`). This file is crucial and specifies:

- Input Files: `input_vcf`, `ped_file` (optional), `ref_genome`.
- Genome Assembly: `assembly` (e.g., hg19, hg38).
- Output & Directories: `output_dir`, `base_dir`.
- Resources: `threads`.
- Annotation Paths: Paths to VEP cache (`vep_cache_dir`), plugins (`vep_plugins_dir`, `vep_plugins_cachedir`), gnomAD (`gnomad_vcf_chrX`), ClinVar (`clinvar_vcf`), CADD (`cadd_base_dir`), and various plugin-specific data files (e.g., `alphamissense_prescore`, `spliceai_snv_prescore`). See the example `config.yaml` for the full list.
- Filtering Thresholds: `af_cutoff`.

**Ensure all paths in your `config.yaml` are correct and accessible from where you run Snakemake.**

## Usage

The recommended way to run the pipeline is using Snakemake, which manages dependencies and execution order.

**Prerequisites:**

1. Activate Conda Environment: Ensure the `priva_acmg` conda environment (or the one created during installation) is activated:
   ```bash
   conda activate priva_acmg
   ```
2. Configuration File: Have your `config.yaml` file ready with all necessary paths and parameters correctly specified.
3. Snakefile Location: Know the path to the `Snakefile`.

**Running the Pipeline:**

Navigate to your working directory (or specify absolute paths in the command).
The config file can be duplicated and customized for batch-specific runs.

**Basic Command:**

```bash
snakemake --snakefile /path/to/your/Snakefile \
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
