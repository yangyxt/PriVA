#!/usr/bin/env bash
# The annotation step wrapper the Nextflow tasks ran up to cd91ec2: one
# function of annotation_vcf.sh with the arguments main_workflow gave it.
# Kept so tests/annotation_step_compare.py can run the shell functions the
# Rust steps (rust/src/pipeline_steps) replaced.
# Arguments: SOURCE_ROOT CONFIG VCF STEP [OUTPUT].
source "$1/scripts/annotation_vcf.sh" || exit 1
config=$2 vcf=$3 step=$4 output=${5:-}
value() { read_yaml "$config" "$1"; }
case $step in
    preprocess)
        preprocess_vcf -i "$vcf" -o "$output" -r "$(value ref_genome)" -t "$(value threads)" ;;
    gnomad)
        anno_agg_gnomAD_data "$vcf" "$(value threads)" "$(value assembly)" "$(value gnomad_vcf_chrX)" ;;
    clinvar)
        anno_clinvar_data "$vcf" "$(value clinvar_vcf)" "$(value threads)" ;;
    control)
        optional=$(value control_vcf_optional)
        anno_control_vcf_allele "$vcf" "$(value control_vcf)" "$(value ped_file)" "$(value threads)" "${optional:-false}" ;;
    vep)
        anno_VEP_data --input_vcf "$vcf" --config "$config" --assembly "$(value assembly)"             --ref_genome "$(value ref_genome)" --vep_cache_dir "$(value vep_cache_dir)"             --vep_plugins_dir "$(value vep_plugins_dir)"             --vep_plugins_cachedir "$(value vep_plugins_cachedir)" --threads "$(value threads)" ;;
    cadd)
        Calculate_CADD "$vcf" "$config" "$output" ;;
    *)
        echo "Unknown annotation step: $step" >&2; exit 2 ;;
esac
