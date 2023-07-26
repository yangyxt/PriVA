#!/usr/bin/env bash
# This script is used to define all the annotation process functions and run the main pipeline
# For now, the script is only used to process short variants (Indels and SNVs). Large structural variants and CNVs are not included.

# The script expect a Family VCF input along with a pedigree table contains the Family pedigree information.
# Notice that the family sample ID should be consistent between the pedigree table and the VCF file.

# For the same set of configurations (arguments), the pipeline should start with the position it ends last time, unless user specifically asked to forcefully rerun the pipeline

# Maintainer: yangyxt@gmail.com


FORCE_RERUN="False"
SELF_SCRIPT="${BASH_SOURCE[0]}"
SCRIPT_DIR="$(dirname "${SELF_SCRIPT}")"

# Source the other script
source "$SCRIPT_DIR/common_bash_utils.sh"


function log() {
    local msg="$1"
    local script_name="${BASH_SOURCE[1]##*/}"
    local func_name="${FUNCNAME[1]}"
    local line_num="${BASH_LINENO[0]}"
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    >&2 echo "[$timestamp] [Script $script_name: Line $line_num] [Func: $func_name] $msg"
}


function main_function_help_message() {

}


function preprocess_vcf() {
	local input_vcf=${1}
	local output_vcf=${input_vcf/.vcf/.preanno.vcf}
	local ref_genome=${2}


	# Test if output_vcf is already valid
	if [[ ${output_vcf} -nt ${input_vcf} ]] && \
		check_vcf_validity ${output_vcf} && \
		check_vcf_multiallelics ${output_vcf} && \
		[[ $(bcftools view -Ov -H ${output_vcf} | cut -f 5 | grep -c '\.') -eq 0 ]] && \
		whether_same_varset ${input_vcf} ${output_vcf} && \
		contain_only_variants ${output_vcf} && \
		[[ ${FORCE_RERUN} ! =~ [Tt][Rr][Uu][Ee]$ ]]; then
		log "The ${output_vcf} is valid and udpated. Skip this function"
		return 0;
    fi

	# First sort the input_vcf,
	# Then normalize the input_vcf with bcftools 
	normalize_vcf ${input_vcf} ${input_vcf/.vcf/.norm.vcf} ${ref_genome} $TMPDIR

	# Then filter the variants where DP value less than or equalt to 5
	bcftools filter -e 'FORMAT/DP[0] <= 5' -Oz -o ${input_vcf/.vcf*/.filtered.vcf.gz} ${input_vcf/.vcf/.norm.vcf} && \
	announce_remove_tmps ${input_vcf/.vcf/.norm.vcf}

	# Then filter out the variants where no patients in this family has alternative alleles

}


function filter_records_with_missingGT_on_pat_vcf {
	local input_vcf=${1}
	local output_vcf=${3}
	local tmp_dir=${4}

	if [[ -z ${tmp_dir} ]]; then local tmp_dir=$TMPDIR; fi

	log "Input vcf is ${input_vcf}"
	local -a patient_IDs=($(echo ${2}))

	if [[ ${patient_IDs} =~ ^[0-9]+$ ]]; then
		local -a patient_cols=($(echo "${patient_IDs[*]}"))
	elif [[ ${input_vcf} =~ \.vcf ]]; then
		local -a patient_cols=($(bcftools view -h ${input_vcf} | \
								 tail -1 | \
								 mawk -F '\t' -v ids="${patient_IDs[*]}" 'BEGIN{split(ids,id_arr," ");} \
								 										  $1 ~ /^#CHROM/{for(i in id_arr) \
																		  					{for(h=1;h<=NF;h++) \
																								{if($h == id_arr[i]) printf "%s ", h;}}}'))
	fi

	log "patient IDs are ${patient_IDs[*]}"
	log "patient IDs corresponding column indices are ${patient_cols[*]}"
	log "Before filtering the records where patients GT are all ref or null, ${input_vcf} has $(wc -l ${input_vcf} | awk '{printf $1}') rows."
	
    vcf_filter_by_GT \
	${input_vcf} \
	"${patient_IDs[*]}" \
	${output_vcf} && \
	display_vcf ${output_vcf} && \
	log "Finish filtering the records where patients GT are all ref or null. ${input_vcf} has $(wc -l ${input_vcf} | awk '{printf $1}') rows." $'\n\n'
}