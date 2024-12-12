#!/usr/bin/env snakemake

import os
import sys
from pathlib import Path
import math

config: "config.yaml"

# Define variables from config with defaults
def get_config(key, default=None):
    return config.get(key, default)

# Required inputs
INPUT_VCF = get_config("input_vcf")
PED_FILE = get_config("ped_file")
ASSEMBLY = get_config("assembly")
REF_GENOME = get_config("ref_genome")
OUTPUT_DIR = get_config("output_dir", "results")
THREADS = get_config("threads", 1)
BASE_DIR = get_config("base_dir", ".")
SCRIPT_DIR = os.path.join(BASE_DIR, "scripts")


# Get input basename for output naming
INPUT_BASE = os.path.splitext(os.path.basename(INPUT_VCF))[0]
if INPUT_BASE.endswith('.vcf'):
    INPUT_BASE = os.path.splitext(INPUT_BASE)[0]

print(f"INPUT_BASE: {INPUT_BASE}", file=sys.stderr)

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)


def get_final_outputs():
    if not os.path.exists(PED_FILE):
        outputs = [
            os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.vcf.gz"),
            os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.tsv"),
            os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.acmg.tsv")
        ]
        print(f"Expected outputs: {outputs}", file=sys.stderr)
        return outputs

    with open(PED_FILE, 'r') as f:
        families = set(line.split()[0] for line in f if not line.startswith('#'))

    outputs = []
    for family in families:
        outputs.append(os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{family}.filtered.vcf.gz"))
        outputs.append(os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{family}.filtered.tsv"))
        outputs.append(os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{family}.filtered.acmg.tsv"))

    print(f"Expected outputs: {outputs}", file=sys.stderr)
    return outputs


rule all:
    input:
        get_final_outputs()


# Step 1: Annotation (shared step)
rule annotate_variants:
    input:
        vcf=INPUT_VCF,
        ped=PED_FILE,
        ref=REF_GENOME
    output:
        vcf=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.vcf.gz")
    params:
        script=os.path.join(SCRIPT_DIR, "annotation_vcf.sh"),
        config=os.path.join(BASE_DIR, "config.yaml")
    threads: THREADS
    shell:
        """
        bash {params.script} main_workflow \
            -i {input.vcf} \
            -p {input.ped} \
            --config {params.config}
        """

# Step 2: Filter variants with family information
rule filter_variants_per_family:
    input:
        vcf=rules.annotate_variants.output.vcf,
        ped=PED_FILE
    output:
        vcf=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{{family}}.filtered.vcf.gz")
    params:
        script=os.path.join(SCRIPT_DIR, "filtration_vcf_per_fam.sh"),
        config=os.path.join(BASE_DIR, "config.yaml")
    threads: THREADS
    shell:
        """
        bash {params.script} \
            -v {input.vcf} \
            -f {wildcards.family} \
            -c {params.config}
        """

# Step 2 alternative: Filter variants without family information
rule filter_variants_no_family:
    input:
        vcf=rules.annotate_variants.output.vcf
    output:
        vcf=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.vcf.gz")
    params:
        script=os.path.join(SCRIPT_DIR, "filtration_vcf_per_fam.sh"),
        config=os.path.join(BASE_DIR, "config.yaml")
    threads: THREADS
    shell:
        """
        bash {params.script} \
            -v {input.vcf} \
            -c {params.config}
        """

# Step 3: Prioritize variants with family information
rule prioritize_variants_per_family:
    input:
        vcf=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{{family}}.filtered.vcf.gz"),
        ped=PED_FILE
    output:
        tsv=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{{family}}.filtered.tsv"),
        acmg=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{{family}}.filtered.acmg.tsv")
    params:
        script=os.path.join(SCRIPT_DIR, "prioritization_vcf_per_fam.sh"),
        config=os.path.join(BASE_DIR, "config.yaml")
    threads: THREADS
    shell:
        """
        bash {params.script} \
            {input.vcf} \
            {params.config} \
            {wildcards.family}
        """

# Step 3 alternative: Prioritize variants without family information
rule prioritize_variants_no_family:
    input:
        vcf=rules.filter_variants_no_family.output.vcf
    output:
        tsv=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.tsv"),
        acmg=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.acmg.tsv")
    params:
        script=os.path.join(SCRIPT_DIR, "prioritization_vcf_per_fam.sh"),
        config=os.path.join(BASE_DIR, "config.yaml")
    threads: THREADS
    shell:
        """
        bash {params.script} \
            {input.vcf} \
            {params.config}
        """
