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

# Step 2: Filter variants per family
rule filter_variants:
    input:
        vcf=rules.annotate_variants.output.vcf,
        ped=PED_FILE
    output:
        vcf=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{{family}}.filtered.vcf.gz")
    params:
        script=os.path.join("scripts", "filtration_vcf_per_fam.sh"),
        config="config.yaml"
    threads: THREADS_PER_FAMILY  # Dynamically allocated threads per family
    shell:
        """
        bash {params.script} \
            -v {input.vcf} \
            -f {wildcards.family} \
            -c {params.config}
        """

# Step 3: Prioritize variants per family
rule prioritize_variants:
    input:
        vcf=rules.filter_variants.output.vcf,
        ped=PED_FILE
    output:
        tsv=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{{family}}.filtered.tsv"),
        mat=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{{family}}.filtered.acmg.tsv")
    params:
        script=os.path.join("scripts", "prioritization_vcf_per_fam.sh"),
        config="config.yaml"
    threads: THREADS_PER_FAMILY  # Dynamically allocated threads per family
    shell:
        """
        bash {params.script} \
            {input.vcf} \
            {params.config} \
            {wildcards.family}
        """
