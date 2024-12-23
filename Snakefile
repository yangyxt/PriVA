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
THREADS = config.get("threads", 1)
BASE_DIR = config.get("base_dir", ".")
SCRIPT_DIR = os.path.join(BASE_DIR, "scripts")


# Get input basename for output naming
INPUT_BASE = os.path.splitext(os.path.basename(INPUT_VCF))[0]
if INPUT_BASE.endswith('.vcf'):
    INPUT_BASE = os.path.splitext(INPUT_BASE)[0]

print(f"INPUT_BASE: {INPUT_BASE}, SCRIPT_DIR: {SCRIPT_DIR}, OUTPUT_DIR: {OUTPUT_DIR}, THREADS: {THREADS}, BASE_DIR: {BASE_DIR}", file=sys.stderr)

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
        vcf=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.vcf.gz"),
        cadd_tsv=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.cadd.tsv")
    params:
        script=os.path.join(SCRIPT_DIR, "annotation_vcf.sh"),
        config=workflow.configfiles[0]
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
        config=workflow.configfiles[0]
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
        config=workflow.configfiles[0]
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
        vcf=rules.filter_variants_per_family.output.vcf,
        cadd_tsv=rules.annotate_variants.output.cadd_tsv,
        ped=PED_FILE
    output:
        tsv=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{{family}}.filtered.tsv"),
        acmg=os.path.join(OUTPUT_DIR, f"{INPUT_BASE}.anno.{{family}}.filtered.acmg.tsv")
    params:
        script=os.path.join(SCRIPT_DIR, "prioritization_vcf_per_fam.sh"),
        config=workflow.configfiles[0]
    threads: THREADS
    shell:
        """
        bash {params.script} \
            -i {input.vcf} \
            -c {params.config} \
            -f {wildcards.family} \
            -t {input.cadd_tsv}
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
    threads: THREADS
    shell:
        """
        bash {params.script} \
            -i {input.vcf} \
            -c {params.config} \
            -t {input.cadd_tsv}
        """

