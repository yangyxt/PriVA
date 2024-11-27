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

    local threads=$(read_config ${config_file} "threads")
    local output_tab=${input_vcf/.vcf*/.tsv}

    [[ -f ${output_tab} ]] && \
    [[ ${output_tab} -nt ${input_vcf} ]] && \
    log "The combined annotation table ${output_tab} is up to date, skip the conversion" && \
    return 0

    ${SCRIPT_DIR}/combine_annotations.py \
    -i ${input_vcf} \
    -c ${CADD_anno_file} \
    -o ${output_tab} \
    -t ${threads} && \
    display_table ${output_tab}
}




if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
    declare -a func_names=($(typeset -f | awk '!/^main[ (]/ && /^[^ {}]+ *\(\)/ { gsub(/[()]/, "", $1); printf $1" ";}'))
    declare -a input_func_names=($(return_array_intersection "${func_names[*]}" "$*"))
    declare -a arg_indices=($(get_array_index "${input_func_names[*]}" "$*"))
    if [[ ${#input_func_names[@]} -gt 0 ]]; then
        log "Seems like the command is trying to directly run a function in this script ${BASH_SOURCE[0]}."
        log "The identified function names are: ${input_func_names[*]}"
        log "All function names are: ${func_names[*]}"
        first_func_ind=${arg_indices}
        log "The identified first func name is at the ${first_func_ind}th input argument, while the total input arguments are: $*"
        following_arg_ind=$((first_func_ind + 1))
        log "Executing: ${*:${following_arg_ind}}"
        "${@:${following_arg_ind}}"
    else
		log "Directly run main_filtration with input args: $*"
		main_filtration "$@"
	fi
fi


