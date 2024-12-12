#!/usr/bin/env bash

# This script is also not part of the pipeline
# It is only used to prepare the clinvar benchmark variant callset

SELF_SCRIPT="$(realpath ${BASH_SOURCE[0]})"
SCRIPT_DIR="$(dirname "${SELF_SCRIPT}")"
BASE_DIR="$(dirname ${SCRIPT_DIR})"
if [[ ${BASE_DIR} == "/" ]]; then BASE_DIR=""; fi
DATA_DIR="${BASE_DIR}/data"

if [[ -z $TMPDIR ]]; then TMPDIR=/tmp; fi

# Source the other script
source "${SCRIPT_DIR}/common_bash_utils.sh"


function add_pseudo_gt_to_clinvar() {
    local input_vcf=$1
    local output_vcf=$2
    local sample_name=${3:-"CLINVAR_TEST"}  # Default sample name is "CLINVAR_TEST"

    # Add GT format field and sample with heterozygous genotypes
    bcftools view -Ov ${input_vcf} | \
    bcftools annotate --header-lines <(echo -e '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">') | \
    bcftools reheader -h <(bcftools view -h ${input_vcf} | \
        sed '$d'; \
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample_name}") | \
    awk -v OFS="\t" '{if($0 !~ /^#/) {print $0,"GT","0/1"} else {print $0}}' | \
    bcftools view -Oz -o "${output_vcf}" && \
    tabix -f -p vcf "${output_vcf}"
}

if [[ ${BASH_SOURCE[0]} == ${0} ]]; then
    "$@"
fi