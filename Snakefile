#!/usr/bin/env snakemake
# Snakefile

import os

configfile: "config.yaml"

# Define variables from the config file for convenience
input_vcf = config["input_vcf"]
ped_file = config["ped_file"]
fam_name = config["fam_name"]
assembly = config["assembly"]
ref_genome = config["ref_genome"]
output_dir = config["output_dir"]
threads = config["threads"]
af_cutoff = config["af_cutoff"]
gnomad_vcf_dir = config["gnomad_vcf_dir"]
clinvar_vcf = config["clinvar_vcf"]
vep_cache_dir = config["vep_cache_dir"]
script_dir = config["script_dir"]

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Define the final output
final_vcf = os.path.join(output_dir, f"{fam_name}.annotated.vcf.gz")

rule all:
    input:
        final_vcf

rule preprocess_vcf:
    input:
        vcf=input_vcf,
        ped=ped_file,
        ref=ref_genome
    output:
        vcf=os.path.join(output_dir, f"{fam_name}.preprocessed.vcf.gz")
    params:
        script=os.path.join(script_dir, "preprocess_vcf.sh"),
        fam_name=fam_name
    threads: threads
    shell:
        """
        bash {params.script} \
            -i {input.vcf} \
            -p {input.ped} \
            -f {params.fam_name} \
            -r {input.ref} \
            -o {output.vcf}
        """

rule filter_pedigree:
    input:
        vcf=rules.preprocess_vcf.output.vcf,
        ped=ped_file
    output:
        vcf=os.path.join(output_dir, f"{fam_name}.filtered.vcf.gz")
    params:
        script=os.path.join(script_dir, "filter_allele_based_on_pedigree_with_py.sh"),
        fam_name=fam_name
    threads: threads
    shell:
        """
        bash {params.script} \
            -i {input.vcf} \
            -p {input.ped} \
            -f {params.fam_name} \
            -o {output.vcf}
        """

rule annotate_gnomad:
    input:
        vcf=rules.filter_pedigree.output.vcf,
        gnomad_dir=gnomad_vcf_dir
    output:
        vcf=os.path.join(output_dir, f"{fam_name}.gnomad_annotated.vcf.gz")
    params:
        script=os.path.join(script_dir, "anno_agg_gnomAD_data.sh"),
        assembly=assembly,
        af_cutoff=af_cutoff
    threads: threads
    shell:
        """
        bash {params.script} \
            -i {input.vcf} \
            -g {input.gnomad_dir} \
            -a {params.assembly} \
            -t {threads} \
            -o {output.vcf} \
            --af_cutoff {params.af_cutoff}
        """

rule annotate_clinvar:
    input:
        vcf=rules.annotate_gnomad.output.vcf,
        clinvar_vcf=clinvar_vcf
    output:
        vcf=final_vcf
    params:
        script=os.path.join(script_dir, "anno_clinvar_data.sh")
    threads: threads
    shell:
        """
        bash {params.script} \
            -i {input.vcf} \
            -c {input.clinvar_vcf} \
            -o {output.vcf}
        """
