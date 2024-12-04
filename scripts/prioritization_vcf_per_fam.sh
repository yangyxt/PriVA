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
    local config_file=${3}

    local threads=$(read_yaml ${config_file} "threads")
    local output_tab=${input_vcf/.vcf*/.tsv}

    [[ -f ${output_tab} ]] && \
    [[ ${output_tab} -nt ${input_vcf} ]] && \
    log "The combined annotation table ${output_tab} is up to date, skip the conversion" && \
    return 0

    python ${SCRIPT_DIR}/combine_annotations.py \
    -i ${input_vcf} \
    -c ${CADD_anno_file} \
    -o ${output_tab} \
    -t ${threads} && \
    display_table ${output_tab}
}


function interpret_splicing_annotations () {
    local input_tab=${1}
    local config_file=${2}

    local splicing_py=${SCRIPT_DIR}/splicing_var_analysis.py

    local alphamissense_tranx_domain_map=$(read_yaml ${config_file} "alphamissense_tranx_domain_map")
    local alphamissense_intolerant_domains=$(read_yaml ${config_file} "alphamissense_intolerant_domains")
    local threads=$(read_yaml ${config_file} "threads")

    python ${splicing_py} \
    --anno_table ${input_tab} \
    --transcript_domain_map ${alphamissense_tranx_domain_map} \
    --intolerant_domains ${alphamissense_intolerant_domains} \
    --threads ${threads} && \
    display_table ${input_tab} && \
    log "The splicing interpretations are saved to ${input_tab}, added with two columns: splicing_lof and splicing_len_changing"
}



function interpret_utr_annotations () {
    local input_tab=${1}
    local config_file=${2}

    local utr_py=${SCRIPT_DIR}/utr_anno_interpret.py

    local alphamissense_tranx_domain_map=$(read_yaml ${config_file} "alphamissense_tranx_domain_map")
    local alphamissense_intolerant_domains=$(read_yaml ${config_file} "alphamissense_intolerant_domains")
    local threads=$(read_yaml ${config_file} "threads")

    python ${utr_py} \
    --variants_table ${input_tab} \
    --domain_map ${alphamissense_tranx_domain_map} \
    --intolerant_domains ${alphamissense_intolerant_domains} \
    --threads ${threads} && \
    display_table ${input_tab} && \
    log "The UTR interpretations are saved to ${input_tab}, added with two columns: 5UTR_lof and 5UTR_length_changing"
}



function assign_acmg_criteria () {
    local input_tab=${1}
    local config_file=${2}

    
    # Preprocess step, interpret the splicing annotations
    interpret_splicing_annotations ${input_tab} ${config_file} || \
    { log "Failed to interpret splicing annotations"; return 1; }

    # Preprocess step, interpret the UTR annotations
    interpret_utr_annotations ${input_tab} ${config_file} || \
    { log "Failed to interpret UTR annotations"; return 1; }
    
    
    local clinvar_pd_stat=$(read_yaml ${config_file} "clinvar_pd_stat")
    local clinvar_aa_stat=$(read_yaml ${config_file} "clinvar_aa_stat")
    local am_pd_stat=$(read_yaml ${config_file} "alphamissense_pd_stat")

    
    
    local acmg_py=${SCRIPT_DIR}/acmg_criteria_assign.py






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


