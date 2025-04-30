# PriVA

A comprehensive pipeline for ACMG-based variant prioritization in genetic analysis.

## Overview

This repository contains scripts and configuration files for running an automated variant prioritization pipeline based on ACMG guidelines. The pipeline integrates multiple annotation resources to evaluate variant pathogenicity and prioritize variants for clinical interpretation.

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

4. **Haplo-insufficiency**: 
   - LOEUF (VEP_plugin)
   - AlphaMissense_mean_score_per_gene

5. **Transcript disruption effect prediction**: 
   - VEP

6. **Clinical variants**: 
   - ClinVar (bcftools annotate)

7. **Population-wise allele frequency + number of homozygous carriers**: 
   - gnomADv4 (bcftools annotate)

8. **Conservation**: 
   - Conservation (VEP_plugin)

9. **Codon based evaluation**: 
   - SameCodon (VEP_plugin, only used with internet connection in database mode, needs separate running)

## Pipeline Structure

The pipeline consists of three main steps:

1. **Annotation**: Annotates variants with multiple resources
   - VEP annotation with plugins
   - CADD scoring
   - gnomAD frequency annotation
   - ClinVar annotation

2. **Filtration**: Filters variants based on family information (if available)
   - Family-specific filtering when pedigree is provided
   - General filtering without family information

3. **Prioritization**: Prioritizes variants according to ACMG guidelines
   - Generates TSV output with prioritized variants
   - Produces ACMG-specific reports

## Installation

The pipeline includes an installation utility (`scripts/install_utils.sh`) that helps set up all required components:

1. Create and configure the conda environment:
   ```bash
   bash scripts/install_utils.sh conda_install_vep path/to/env.yaml
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

The pipeline is configured through a YAML file (`config.yaml`) that specifies:

- Input/output paths
- Reference genome and assembly version
- Annotation resource locations
- Computational parameters (threads, etc.)
- Filtering thresholds

## Usage

To run the pipeline:

```bash
snakemake --configfile config.yaml -j <threads>
```

A more comprehensive example:

```bash
nohup snakemake --snakefile /paedyl01/disk1/yangyxt/PriVA/Snakefile --cores 50 --configfile /paedyl01/disk1/yangyxt/test_acmg_auto/test_clinvar/config_clinvar.hg38.yaml --printshellcmds --verbose > /paedyl01/disk1/yangyxt/test_acmg_auto/test_clinvar/clinvar_2star_snakemake_pipeline.hg38.log 2>&1 &
```

For family-specific analysis, ensure your pedigree file is properly formatted and specified in the config file.

## Output Files

The pipeline produces several output files:

1. Annotated VCF files (`.anno.vcf.gz`)
2. CADD score files (`.anno.cadd.tsv`)
3. Filtered VCF files (`.anno.filtered.vcf.gz` or `.anno.<family>.filtered.vcf.gz`)
4. Prioritized variant tables (`.anno.filtered.tsv` or `.anno.<family>.filtered.tsv`)
5. ACMG-classified variant tables (`.anno.filtered.acmg.tsv` or `.anno.<family>.filtered.acmg.tsv`)

## Workflow Modes

The pipeline can run in two modes:

- **With family information**: Uses pedigree data to perform family-specific filtering and prioritization
- **Without family information**: Performs general filtering and prioritization

## License

[Specify your license information here]

## Contact

Maintainer: yangyxt@hku.hk, yangyxt@gmail.com
