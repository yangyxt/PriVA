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

    # First, create a complete header with FORMAT fields
    (bcftools view -h ${input_vcf} | grep -v "#CHROM"; \
     echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'; \
     echo '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'; \
     echo '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">'; \
     bcftools view -h ${input_vcf} | grep "#CHROM" | \
        sed "s/\$/\tFORMAT\t${sample_name}/") > header.tmp

    # Then use this header to create the new VCF
    (cat header.tmp; \
     bcftools view -H ${input_vcf} | \
     awk -v OFS="\t" '{print $0,"GT:AD:DP","0/1:15,15:30"}') | \
    bcftools view -Oz -o "${output_vcf}"

    # Clean up and index
    rm header.tmp
    tabix -f -p vcf "${output_vcf}"
}



function filter_clinvar_records() {
    local input_vcf=$1
    local output_vcf=$2
    
    # Filter:
    # 1. Remove empty ALT with -e 'ALT!="."'
    # 2. Keep records with review status matching exactly the high confidence statuses
    bcftools filter -e 'ALT=="."' -Ou "${input_vcf}" | \
	bcftools filter -i 'INFO/CLNREVSTAT ~ "practice_guideline" || INFO/CLNREVSTAT ~ "reviewed_by_expert_panel" || INFO/CLNREVSTAT ~ "_no_conflicts" || INFO/CLNREVSTAT ~ "_single_submitter"' -Ou | \
	bcftools sort -Oz -o "${output_vcf}" && \
    tabix -f -p vcf "${output_vcf}"

    # Display summary
    log "Filtered VCF file created: ${output_vcf}"
    log "Original record count: $(bcftools view -H ${input_vcf} | wc -l)"
    log "Filtered record count: $(bcftools view -H ${output_vcf} | wc -l)"
}


function run_InterVar() {
	local input_vcf=${1}
	local output_dir=${2}

	local intervar_py="/paedyl01/disk1/yangyxt/Tools/InterVar-2.2.1/Intervar.py"
	local intervar_config="/paedyl01/disk1/yangyxt/Tools/InterVar-2.2.1/config.ini"

	# You need to download the mim2gene.txt from OMIM and place it in the corresponding position in the intervar_directory
	# You also need to put the input and output directory values in the config.ini file
	# Here is the link to the github page: https://github.com/WGLab/InterVar

	# Run InterVar and there is no thread settings so I guess this is single threaded.
	python ${intervar_py} -c ${intervar_config}
}



function run_TAPES() {
	local input_vcf=${1}
	local output_dir=${2}

	local tapes_py="/paedyl01/disk1/yangyxt/Tools/tapes/tapes.py"
	local tapes_config="/paedyl01/disk1/yangyxt/Tools/tapes/config.ini"

	# You kind of need to first tell tapes.py where you store the annovar and then let it scan and download necessary annotation resources.
	# The link of the TAPES github page is: https://github.com/a-xavier/tapes

	# Run TAPES
	python3 ${tapes_py} annotate -i ${input_vcf} -o ${output_dir}/tapes_annotated_clinvar.vcf --acmg && \
	log "TAPES annotation completed. Output file: ${output_dir}/tapes_annotated_clinvar.vcf" && \
	python3 ${tapes_py} sort -i ${output_dir}/tapes_annotated_clinvar.hg19_multianno.vcf -o ${output_dir}/tapes_annotated_clinvar.hg19_multianno.sorted.tsv --tab && \
	log "TAPES sorting completed. Output file: ${output_dir}/tapes_annotated_clinvar.hg19_multianno.sorted.tsv" && \
	display_table ${output_dir}/tapes_annotated_clinvar.hg19_multianno.sorted.tsv
}



function run_GeneBe() {
	local input_vcf=${1}
	local output_dir=${2}

	local genebe_py="/paedyl01/disk1/yangyxt/Tools/GeneBe/GeneBe.py"
	local genebe_config="/paedyl01/disk1/yangyxt/Tools/GeneBe/config.ini"

	# Run GeneBe
	python ${genebe_py} -c ${genebe_config}
}


if [[ ${BASH_SOURCE[0]} == ${0} ]]; then
    "$@"
fi