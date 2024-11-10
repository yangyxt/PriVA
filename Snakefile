#!/usr/bin/env snakemake

import os
from pathlib import Path

configfile: "config.yaml"

# Define variables from config with defaults
def get_config(key, default=None):
    return config.get(key, default)

INPUT_VCF = get_config("input_vcf")
PED_FILE = get_config("ped_file")
FAM_NAME = get_config("fam_name")
ASSEMBLY = get_config("assembly")
REF_GENOME = get_config("ref_genome")
OUTPUT_DIR = get_config("output_dir", "results")
THREADS = get_config("threads", 1)
AF_CUTOFF = get_config("af_cutoff", 0.01)
GNOMAD_VCF_DIR = str(Path(get_config("gnomad_vcf_chrX")).parent)
CLINVAR_VCF = get_config("clinvar_vcf")
VEP_CACHE_DIR = get_config("vep_cache_dir")
VEP_PLUGINS_DIR = get_config("vep_plugins_dir")
VEP_PLUGINS_CACHEDIR = get_config("vep_plugins_cachedir")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define the final outputs for each major step
def get_output_path(filename):
    return os.path.join(OUTPUT_DIR, f"{FAM_NAME}.{filename}.vcf.gz")

rule all:
    input:
        get_output_path("annotated"),  # Final annotated VCF
        get_output_path("filtered"),   # After variant filtration
        get_output_path("prioritized") # After variant prioritization

# Step 1: Annotation
rule annotate_variants:
    input:
        vcf=INPUT_VCF,
        ped=PED_FILE,
        ref=REF_GENOME
    output:
        vcf=get_output_path("annotated")
    params:
        script=os.path.join("scripts", "annotation_vcf.sh"),
        fam_name=FAM_NAME,
        gnomad_dir=GNOMAD_VCF_DIR,
        clinvar=CLINVAR_VCF,
        vep_cache=VEP_CACHE_DIR,
        vep_plugins=VEP_PLUGINS_DIR,
        vep_cache_plugins=VEP_PLUGINS_CACHEDIR
    threads: THREADS
    shell:
        """
        bash {params.script} main_workflow \
            --input_vcf {input.vcf} \
            --ped_file {input.ped} \
            --fam_name {params.fam_name} \
            --assembly {ASSEMBLY} \
            --ref_genome {input.ref} \
            --output_dir {OUTPUT_DIR} \
            --threads {threads} \
            --af_cutoff {AF_CUTOFF} \
            --gnomad_vcf_chrX {params.gnomad_dir}/gnomad.genomes.v3.1.2.sites.chr1.vcf.gz \
            --clinvar_vcf {params.clinvar} \
            --vep_cache_dir {params.vep_cache} \
            --vep_plugins_dir {params.vep_plugins} \
            --vep_plugins_cachedir {params.vep_cache_plugins}
        """

# Step 2: Variant Filtration
rule filter_variants:
    input:
        vcf=rules.annotate_variants.output.vcf,
        ped=PED_FILE
    output:
        vcf=get_output_path("filtered")
    params:
        script=os.path.join("scripts", "filter_variants.sh"),
        fam_name=FAM_NAME
    threads: THREADS
    shell:
        "bash {params.script} -i {input.vcf} -p {input.ped} -f {params.fam_name} -o {output.vcf}"

# Step 3: Variant Prioritization
rule prioritize_variants:
    input:
        vcf=rules.filter_variants.output.vcf,
        ped=PED_FILE
    output:
        vcf=get_output_path("prioritized")
    params:
        script=os.path.join("scripts", "prioritize_variants.sh"),
        fam_name=FAM_NAME
    threads: THREADS
    shell:
        "bash {params.script} -i {input.vcf} -p {input.ped} -f {params.fam_name} -o {output.vcf}"
