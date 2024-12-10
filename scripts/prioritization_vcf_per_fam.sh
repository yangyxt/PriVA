#!/usr/bin/env bash

# This script is used to prioritize the variants in the VCF file for each family
# The script expect a Family VCF input along with a pedigree table contains the Family pedigree information.
# Notice that the family sample ID should be consistent between the pedigree table and the VCF file.
# Notice that the proband of the family should stay as the first record of the family in the pedigree table

# Maintainer: yangyxt@gmail.com, yangyxt@hku.hk

SCRIPT_DIR=$(dirname $(readlink -f $0))
source ${SCRIPT_DIR}/common_bash_utils.sh

# In this script, we have several steps:
# 1. Convert the annotated VCF to a Table with transcript specific annotations as rows
# 2. Assign the ACMG criterias to each variant-transcript annotation record
# 3. Prioritize the variants based on the ACMG criterias


function prepare_combined_tab () {
    local input_vcf=${1}  # This input VCF file should be annotated with many things
    local CADD_anno_file=${2}
    local threads=${3}

    local output_tab=${input_vcf/.vcf*/.tsv}

    [[ -f ${output_tab} ]] && \
    [[ ${output_tab} -nt ${input_vcf} ]] && \
    log "The combined annotation table ${output_tab} is up to date, skip the conversion" && \
    return 0

	local hpo_tab=${BASE_DIR}/data/hpo/genes_to_phenotype.collapse.tsv.gz

    python ${SCRIPT_DIR}/combine_annotations.py \
    -i ${input_vcf} \
    -c ${CADD_anno_file} \
    -o ${output_tab} \
	-p ${hpo_tab} \
    -t ${threads} && \
    display_table ${output_tab}
}


function interpret_splicing_annotations () {
    local input_tab=${1}
    local config_file=${2}

    local splicing_py=${SCRIPT_DIR}/splicing_var_analysis.py

    local alphamissense_tranx_domain_map=$(read_yaml ${config_file} "alphamissense_tranx_domain_map")
    local intolerant_domains=$(read_yaml ${config_file} "all_intolerant_domains")
    local threads=$(read_yaml ${config_file} "threads")

    python ${splicing_py} \
    --anno_table ${input_tab} \
    --transcript_domain_map ${alphamissense_tranx_domain_map} \
    --intolerant_domains ${intolerant_domains} \
    --threads ${threads} && \
    display_table ${input_tab} && \
    log "The splicing interpretations are saved to ${input_tab}, added with two columns: splicing_lof and splicing_len_changing"
}



function interpret_utr_annotations () {
    local input_tab=${1}
    local config_file=${2}

    local utr_py=${SCRIPT_DIR}/utr_anno_interpret.py

    local alphamissense_tranx_domain_map=$(read_yaml ${config_file} "alphamissense_tranx_domain_map")
    local intolerant_domains=$(read_yaml ${config_file} "all_intolerant_domains")
    local threads=$(read_yaml ${config_file} "threads")

    python ${utr_py} \
    --variants_table ${input_tab} \
    --domain_map ${alphamissense_tranx_domain_map} \
    --intolerant_domains ${intolerant_domains} \
    --threads ${threads} && \
    display_table ${input_tab} && \
    log "The UTR interpretations are saved to ${input_tab}, added with two columns: 5UTR_lof and 5UTR_length_changing"
}



function assign_acmg_criteria () {
    local input_tab=${1}
    local config_file=${2}
    local fam_name=${3}
    
    # Preprocess step, interpret the splicing annotations
    interpret_splicing_annotations ${input_tab} ${config_file} || \
    { log "Failed to interpret splicing annotations for ${input_tab}"; return 1; }

    # Preprocess step, interpret the UTR annotations
    interpret_utr_annotations ${input_tab} ${config_file} || \
    { log "Failed to interpret UTR annotations for ${input_tab}"; return 1; }
    
    local mean_am_score_table=${BASE_DIR}/data/alphamissense/alphamissense_mean_score.tsv
    local ped_table=$(read_yaml ${config_file} "ped_file")
    local clinvar_aa_dict_pkl=$(read_yaml ${config_file} "clinvar_aa_stat")
    local intolerant_domains_pkl=$(read_yaml ${config_file} "all_intolerant_domains")
    local domain_mechanism_tsv=$(read_yaml ${config_file} "clinvar_intolerance_mechanisms")
    local alt_disease_vcf=$(read_yaml ${config_file} "alt_disease_vcf")
    local gnomAD_extreme_rare_threshold=$(read_yaml ${config_file} "extreme_rare_PAF")
    local expected_incidence=$(read_yaml ${config_file} "exp_disease_incidence")

    # Test whether the function can be skipped
    [[ -f ${input_tab} ]] && \
    [[ ${input_tab} -nt ${mean_am_score_table} ]] && \
    [[ ${input_tab} -nt ${ped_table} ]] && \
    [[ ${input_tab} -nt ${clinvar_aa_dict_pkl} ]] && \
    [[ ${input_tab} -nt ${intolerant_domains_pkl} ]] && \
    [[ ${input_tab} -nt ${domain_mechanism_tsv} ]] && \
    [[ ${input_tab} -nt ${alt_disease_vcf} ]] && \
    check_table_column ${input_tab} "ACMG_quant_score" && \
    check_table_column ${input_tab} "ACMG_class" && \
    check_table_column ${input_tab} "ACMG_criteria" && \
    log "The ACMG criterias are already assigned for ${input_tab}, skip the assignment" && \
    return 0
    
    local acmg_py=${SCRIPT_DIR}/acmg_criteria_assign.py

    log "Running the following command to assign the ACMG criterias: python ${acmg_py} --anno_table ${input_tab} --am_score_table ${mean_am_score_table} --ped_table ${ped_table} --fam_name ${fam_name} --clinvar_aa_dict_pkl ${clinvar_aa_dict_pkl} --intolerant_domains_pkl ${intolerant_domains_pkl} --domain_mechanism_tsv ${domain_mechanism_tsv} --alt_disease_vcf ${alt_disease_vcf} --gnomAD_extreme_rare_threshold ${gnomAD_extreme_rare_threshold} --expected_incidence ${expected_incidence} --threads ${threads}"
    python ${acmg_py} \
    --anno_table ${input_tab} \
    --am_score_table ${mean_am_score_table} \
    --ped_table ${ped_table} \
    --fam_name ${fam_name} \
    --clinvar_aa_dict_pkl ${clinvar_aa_dict_pkl} \
    --intolerant_domains_pkl ${intolerant_domains_pkl} \
    --domain_mechanism_tsv ${domain_mechanism_tsv} \
    --alt_disease_vcf ${alt_disease_vcf} \
    --gnomAD_extreme_rare_threshold ${gnomAD_extreme_rare_threshold} \
    --expected_incidence ${expected_incidence} \
    --threads ${threads} && \
    display_table ${input_tab} && \
    local output_acmg_mat=${input_tab::-4}.acmg.tsv && \
    log "The ACMG criterias are assigned for ${input_tab}, added with three columns: ACMG_quant_score, ACMG_class, ACMG_criteria, and the output matrix is saved to ${output_acmg_mat}" && \
    display_table ${output_acmg_mat} || \
    { log "Failed to assign the ACMG criterias for ${input_tab}"; return 1; }
}



function main_prioritization () {
    local input_vcf=${1}
    local config_file=${2}
    local fam_name=${3}
    
	# Read the expected number of threads
	local threads=$(read_yaml ${config_file} "threads")

    # Preprocess step, convert the annotated VCF to a Table with transcript specific annotations as rows
    local CADD_anno_file=$(read_yaml ${config_file} "cadd_output_file")
	prepare_combined_tab ${input_vcf} ${CADD_anno_file} ${threads} && \
    local input_tab=${input_vcf/.vcf*/.tsv} || \
    { log "Failed to prepare the combined annotation table"; return 1; }

    # Assign the ACMG criterias to each variant-transcript annotation record
    assign_acmg_criteria ${input_tab} ${config_file} ${fam_name} || \
    { log "Failed to assign the ACMG criterias"; return 1; }

}




if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
    declare -a func_names=($(typeset -f | awk '!/^main[ (]/ && /^[^ {}]+ *\(\)/ { gsub(/[()]/, "", $1); printf $1" ";}'))
    declare -a input_func_names=($(return_array_intersection "${func_names[*]}" "$*"))
    declare -a arg_indices=($(get_array_index "${input_func_names[*]}" "$*"))
    log "All function names are: ${func_names[*]}"
    if [[ ${#input_func_names[@]} -gt 0 ]]; then
        log "Seems like the command is trying to directly run a function in this script ${BASH_SOURCE[0]}."
        log "The identified function names are: ${input_func_names[*]}"
        
        first_func_ind=${arg_indices}
        log "The identified first func name is at the ${first_func_ind}th input argument, while the total input arguments are: $*"
        following_arg_ind=$((first_func_ind + 1))
        log "Executing: ${*:${following_arg_ind}}"
        "${@:${following_arg_ind}}"
    else
		log "Directly run main_prioritization with input args: $*"
		main_prioritization "$@"
	fi
fi


