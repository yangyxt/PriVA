#!/usr/bin/env snakemake

import os
from pathlib import Path
import math

configfile: "config.yaml"

# Define variables from config with defaults
def get_config(key, default=None):
    return config.get(key, default)

# Required inputs
INPUT_VCF = get_config("input_vcf")
PED_FILE = get_config("ped_file")
ASSEMBLY = get_config("assembly")
REF_GENOME = get_config("ref_genome")

# Get input basename for output naming
INPUT_BASE = os.path.splitext(os.path.basename(INPUT_VCF))[0]
if INPUT_BASE.endswith('.vcf'):
    INPUT_BASE = os.path.splitext(INPUT_BASE)[0]

# Optional configs with defaults
OUTPUT_DIR = get_config("output_dir", "results")
THREADS = get_config("threads", 1)

# Get list of families from PED file
def get_families():
    with open(PED_FILE) as f:
        families = {line.split()[0] for line in f if not line.startswith('#')}
    return sorted(list(families))

FAMILIES = get_families()
N_FAMILIES = len(FAMILIES)

# Calculate threads per family (minimum 1 thread)
def get_threads_per_family():
    return max(1, math.floor(THREADS / N_FAMILIES))

THREADS_PER_FAMILY = get_threads_per_family()

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define the final outputs for each family
def get_family_outputs(wildcards):
    return [
        os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{fam}.filtered.vcf.gz") for fam in FAMILIES
    ] + [
        os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{fam}.filtered.tsv") for fam in FAMILIES
    ] + [
        os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{fam}.filtered.acmg.tsv") for fam in FAMILIES
    ]

rule all:
    input:
        get_family_outputs

# Step 1: Annotation (shared step for all families)
rule annotate_variants:
    input:
        vcf=INPUT_VCF,
        ped=PED_FILE,
        ref=REF_GENOME
    output:
        vcf=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.vcf.gz")
    params:
        script=os.path.join("scripts", "annotation_vcf.sh"),
        config="config.yaml"
    threads: THREADS  # Uses full thread allocation as it's a shared step
    shell:
        """
        bash {params.script} main_workflow \
            -i {input.vcf} \
            -p {input.ped} \
            --config {params.config}
        """

# Checkpoint to validate the pedigree file
checkpoint validate_pedigree:
    input:
        ped=PED_FILE
    output:
        valid=os.path.join(OUTPUT_DIR, "pedigree_valid.txt")
    run:
        if os.path.isfile(input.ped) and os.path.getsize(input.ped) > 0:
            with open(output.valid, 'w') as f:
                f.write("valid")
        else:
            with open(output.valid, 'w') as f:
                f.write("invalid")

# Function to determine the next steps based on pedigree validation
def determine_next_steps(wildcards):
    valid_file = checkpoints.validate_pedigree.output.valid
    if os.path.isfile(valid_file):
        with open(valid_file) as f:
            status = f.read().strip()
        if status == "valid":
            return expand(os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{{family}}.filtered.vcf.gz"), family=FAMILIES)
    return [os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.vcf.gz")]

# Step 2: Filter variants
rule filter_variants:
    input:
        vcf=rules.annotate_variants.output.vcf,
        ped=PED_FILE,
        valid=checkpoints.validate_pedigree.output.valid
    output:
        vcf=determine_next_steps
    params:
        script=os.path.join("scripts", "filtration_vcf_per_fam.sh"),
        config="config.yaml"
    threads: THREADS_PER_FAMILY
    shell:
        """
        if [[ $(cat {input.valid}) == "valid" ]]; then
            bash {params.script} \
                -v {input.vcf} \
                -f {wildcards.family} \
                -c {params.config}
        else
            bash {params.script} \
                -v {input.vcf} \
                -c {params.config}
        fi
        """

# Function to determine output paths for prioritization
def determine_prioritization_outputs(wildcards):
    valid_file = checkpoints.validate_pedigree.output.valid
    if os.path.isfile(valid_file):
        with open(valid_file) as f:
            status = f.read().strip()
        if status == "valid":
            return {
                'tsv': os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{wildcards.family}.filtered.tsv"),
                'mat': os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{wildcards.family}.filtered.acmg.tsv")
            }
    return {
        'tsv': os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.tsv"),
        'mat': os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.acmg.tsv")
    }

# Step 3: Prioritize variants
rule prioritize_variants:
    input:
        vcf=rules.filter_variants.output.vcf,
        ped=PED_FILE,
        valid=checkpoints.validate_pedigree.output.valid
    output:
        unpack(determine_prioritization_outputs)
    params:
        script=os.path.join("scripts", "prioritization_vcf_per_fam.sh"),
        config="config.yaml"
    threads: THREADS_PER_FAMILY
    shell:
        """
        if [[ $(cat {input.valid}) == "valid" ]]; then
            bash {params.script} \
                {input.vcf} \
                {params.config} \
                {wildcards.family}
        else
            bash {params.script} \
                {input.vcf} \
                {params.config}
        fi
        """
