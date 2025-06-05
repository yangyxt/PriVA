#!/usr/bin/env snakemake

import os
import sys
from pathlib import Path
import math

configfile: "config.yaml"


# Required inputs
INPUT_VCF = config.get("input_vcf")
PED_FILE = config.get("ped_file")
ASSEMBLY = config.get("assembly")
REF_GENOME = config.get("ref_genome")
OUTPUT_DIR = config.get("output_dir", "results")
THREADS = config.get("threads", 8) # Total threads for annotation/non-family steps (default 8)
THREADS_PER_FAM = config.get("threads_per_fam", 1) # Threads per family job (default 1)
BASE_DIR = config.get("base_dir", ".")
SCRIPT_DIR = os.path.join(BASE_DIR, "scripts")


# Get input basename for output naming
INPUT_BASE = os.path.splitext(os.path.basename(INPUT_VCF))[0]
if INPUT_BASE.endswith('.vcf'):
    INPUT_BASE = os.path.splitext(INPUT_BASE)[0]

# Use f-string for cleaner output formatting
print(f"INPUT_BASE: {INPUT_BASE}, SCRIPT_DIR: {SCRIPT_DIR}, OUTPUT_DIR: {OUTPUT_DIR}", file=sys.stderr)
print(f"THREADS - Total/Annotation: {THREADS}, Per Family: {THREADS_PER_FAM}", file=sys.stderr) # Updated print statement


# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Ensure logs directory exists
os.makedirs(os.path.join(OUTPUT_DIR, "logs"), exist_ok=True)


def get_final_outputs():
    # Check if PED_FILE exists and is not empty or just whitespace
    ped_exists = PED_FILE and os.path.exists(PED_FILE) and os.path.getsize(PED_FILE) > 0
    families = set()

    if ped_exists:
        try:
            with open(PED_FILE, 'r') as f:
                # Ensure lines are split correctly and handle potential empty lines/comments robustly
                families = set(line.split()[0] for line in f if line.strip() and not line.startswith('#') and line.split())
        except Exception as e:
            print(f"Warning: Could not read pedigree file {PED_FILE}: {e}", file=sys.stderr)
            ped_exists = False # Treat as if no ped file if reading fails

    if not ped_exists or not families:
        print("No valid pedigree file found or families detected. Running in single VCF mode.", file=sys.stderr)
        outputs = [
            os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.vcf.gz"),
            os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.tsv"),
            os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.acmg.tsv")
        ]
    else:
        print(f"Found families: {families}. Running in multi-family mode.", file=sys.stderr)
        outputs = []
        for family in families:
            outputs.append(os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{family}.filtered.vcf.gz"))
            outputs.append(os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{family}.filtered.tsv"))
            outputs.append(os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{family}.filtered.acmg.tsv"))

    print(f"Expected final outputs: {outputs}", file=sys.stderr)
    return outputs


rule all:
    input:
        get_final_outputs()


# Step 1: Annotation (shared step)
rule annotate_variants:
    input:
        vcf=INPUT_VCF,
        # ped=PED_FILE, # Ped file not needed as direct input if script reads from config
        ref=REF_GENOME
    output:
        vcf=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.vcf.gz"),
        cadd_tsv=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.cadd.tsv")
    params:
        script=os.path.join(SCRIPT_DIR, "annotation_vcf.sh"),
        config=workflow.configfiles[0]
        # Pedigree file path is expected to be read from the config file by the script if needed
    threads: THREADS # Use total threads for annotation
    log:
        os.path.join(OUTPUT_DIR, "logs", f"{INPUT_BASE}.annotation.log")
    shell:
        """
        echo "[$(date)] Starting annotation step. Detailed log: {log}" >&2
        bash {params.script} main_workflow \
            -i {input.vcf} \
            --config {params.config} > {log} 2>&1
        """

# Step 2: Filter variants with family information
rule filter_variants_per_family:
    input:
        vcf=rules.annotate_variants.output.vcf,
        ped=PED_FILE # Ped file is needed here to know which families exist
    output:
        vcf=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{{family}}.filtered.vcf.gz")
    params:
        script=os.path.join(SCRIPT_DIR, "filtration_vcf_per_fam.sh"),
        config=workflow.configfiles[0]
    threads: THREADS_PER_FAM # Use per-family threads
    log:
        os.path.join(OUTPUT_DIR, "logs", f"{INPUT_BASE}.{{family}}.filtration.log")
    shell:
        """
        echo "[$(date)] Starting filtration for family {wildcards.family}. Detailed log: {log}" >&2
        bash {params.script} \
            -v {input.vcf} \
            -f {wildcards.family} \
            -c {params.config} > {log} 2>&1
        """

# Step 2 alternative: Filter variants without family information
rule filter_variants_no_family:
    input:
        vcf=rules.annotate_variants.output.vcf
    output:
        vcf=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.vcf.gz")
    params:
        script=os.path.join(SCRIPT_DIR, "filtration_vcf_per_fam.sh"),
        config=workflow.configfiles[0]
    threads: THREADS # Use total threads for non-family filtering
    log:
        os.path.join(OUTPUT_DIR, "logs", f"{INPUT_BASE}.no_family.filtration.log")
    shell:
        """
        echo "[$(date)] Starting filtration (no family mode). Detailed log: {log}" >&2
        bash {params.script} \
            -v {input.vcf} \
            -c {params.config} > {log} 2>&1
        """

# Step 3: Prioritize variants with family information
rule prioritize_variants_per_family:
    input:
        vcf=rules.filter_variants_per_family.output.vcf,
        cadd_tsv=rules.annotate_variants.output.cadd_tsv,
        ped=PED_FILE # Ped file might be needed by script for context
    output:
        tsv=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{{family}}.filtered.tsv"),
        acmg=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{{family}}.filtered.acmg.tsv")
    params:
        script=os.path.join(SCRIPT_DIR, "prioritization_vcf_per_fam.sh"),
        config=workflow.configfiles[0]
    threads: THREADS_PER_FAM # Use per-family threads
    log:
        os.path.join(OUTPUT_DIR, "logs", f"{INPUT_BASE}.{{family}}.prioritization.log")
    shell:
        """
        echo "[$(date)] Starting prioritization for family {wildcards.family}. Detailed log: {log}" >&2
        bash {params.script} \
            -i {input.vcf} \
            -c {params.config} \
            -f {wildcards.family} \
            -t {input.cadd_tsv} > {log} 2>&1
        """

# Step 3 alternative: Prioritize variants without family information
rule prioritize_variants_no_family:
    input:
        vcf=rules.filter_variants_no_family.output.vcf,
        cadd_tsv=rules.annotate_variants.output.cadd_tsv
    output:
        tsv=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.tsv"),
        acmg=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.filtered.acmg.tsv")
    params:
        script=os.path.join(SCRIPT_DIR, "prioritization_vcf_per_fam.sh"),
        config=workflow.configfiles[0]
    threads: THREADS # Use total threads for non-family prioritization
    log:
        os.path.join(OUTPUT_DIR, "logs", f"{INPUT_BASE}.no_family.prioritization.log")
    shell:
        """
        echo "[$(date)] Starting prioritization (no family mode). Detailed log: {log}" >&2
        bash {params.script} \
            -i {input.vcf} \
            -c {params.config} \
            -t {input.cadd_tsv} > {log} 2>&1
        """

