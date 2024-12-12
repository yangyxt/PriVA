#! /usr/bin/env bash

# This script is used to filter the variants in the VCF file for each family
# The script expect a Family VCF input along with a pedigree table contains the Family pedigree information.
# Notice that the family sample ID should be consistent between the pedigree table and the VCF file.
# Notice that the proband of the family should stay as the first record of the family in the pedigree table

# Maintainer: yangyxt@gmail.com, yangyxt@hku.hk

SCRIPT_DIR=$(dirname $(readlink -f $0))
source ${SCRIPT_DIR}/common_bash_utils.sh


# Just two very simple functions to filter out the variants based on the pedigree information
# The first function is to filter out the variants where patients GT are the same with controls GT info
# The second function is to filter out the variants based on the allele frequency


function filter_allele_based_on_pedigree_with_py {
    local OPTIND i o p f
    while getopts i:o::p::f:: args
    do
        case ${args} in
            i) local input_vcf=$OPTARG ;;
            p) local pedfile=$OPTARG ;;
            o) local output_vcf=$OPTARG ;;
            f) local family_name=$OPTARG ;;
            *) echo "No argument passed. At least pass an argument specifying the family ID"
        esac
    done

    if [[ -z ${output_vcf} ]]; then
        local output_vcf=${input_vcf/.vcf/.filtered.vcf}
    fi

    if [[ -z ${family_name} ]]; then
        local family_name=$(basename ${input_vcf} | awk -F '.' '{printf "%s", $1;}')
    fi

    if [[ ${output_vcf} -nt ${input_vcf} ]] && \
	   [[ $(check_vcf_validity ${output_vcf}) ]] && \
       [[ $(check_vcf_lineno ${output_vcf}) -le $(check_vcf_lineno ${input_vcf}) ]]; then
       log "The output vcf ${output_vcf} is valid and updated"
       return 0;
    fi

    log "Start filtering the records where patients GT are the same with controls GT info"
    log "Running commands: python3 ${SCRIPT_DIR}/pedigree_filtration_per_fam.py \
    -v ${input_vcf} \
    -p ${pedfile} \
    -f ${family_name} \
    -o ${output_vcf}"

    python3 ${SCRIPT_DIR}/pedigree_filtration_per_fam.py \
    -v ${input_vcf} \
    -p ${pedfile} \
    -f ${family_name} \
    -o ${output_vcf} && \
    tabix -p vcf -f ${output_vcf} && \
	display_vcf ${output_vcf} || { \
    log "Failed to filter the records where patients GT are the same with controls GT info"; \
    return 1; }
    
    log "Finish filtering the records where control sample has homozygous or hemizygous GTs or the records no patients carrying the variant allele."$'\n\n'
}



function filter_af () {
	local config_file
	local input_vcf
    local threads

    # Use getopts to parse the input arguments
    while getopts i:o::c:t:: args
    do
        case ${args} in
            i) local input_vcf=$OPTARG ;;
            o) local output_vcf=$OPTARG ;;
            c) local config_file=$OPTARG ;;
            t) local threads=$OPTARG ;;
            *) echo "No argument passed. At least pass an argument specifying the family ID"
        esac
    done

    [[ -z ${threads} ]] && threads=1

    local af_cutoff=$(read_yaml ${config_file} "af_cutoff")

    [[ -z ${output_vcf} ]] && \
    local tmp_tag=$(randomID) && \
    local output_vcf=${input_vcf/.vcf*/.${tmp_tag}.vcf.gz} && \
    local replace="TRUE" || \
    { log "Failed to generate a temporary output VCF file name. Quit with error"; return 1; }

    check_vcf_validity ${input_vcf} && \
    [[ $(bcftools view -h ${input_vcf} | grep -c "AF_joint > ${af_cutoff}") -gt 0 ]] && \
    log "The input VCF ${input_vcf} already been filtered on allele frequency at the cutoff ${af_cutoff}" && \
    return 0

    bcftools view --threads ${threads} -e "AF_joint > ${af_cutoff}" -Ou ${input_vcf} | \
    bcftools sort -Oz -o ${output_vcf} || \
    { log "Failed to filter the variants based on the allele frequency"; return 1; }
    
    if [[ ${replace} == "TRUE" ]]; then
        mv ${output_vcf} ${input_vcf} && \
        tabix -p vcf -f ${input_vcf} && \
        display_vcf ${input_vcf}
    else
        tabix -p vcf -f ${output_vcf} && \
        display_vcf ${output_vcf}
    fi
}


function extract_fam_vcf () {
    local input_vcf=${1}
    local ped_file=${2}
    local family_name=${3}
    local threads=${4}

    [[ -z ${threads} ]] && threads=1

    local fam_members=$(awk -v family_name=${family_name} '$1 == family_name {print $2}' ${ped_file} | tr '\n' ',')
    local tmp_tag=$(randomID)
    local output_vcf=${input_vcf/.vcf*/.${tmp_tag}.vcf.gz}

    check_vcf_validity ${input_vcf} || { \
    log "The input VCF ${input_vcf} is not valid"; \
    return 1; }

    check_vcf_validity ${input_vcf/.vcf*/.${family_name}.vcf.gz} 1 ${fam_members::-1} && \
    [[ ${input_vcf/.vcf*/.${family_name}.vcf.gz} -nt ${input_vcf} ]] && \
    { log "The output VCF ${input_vcf/.vcf*/.${family_name}.vcf.gz} is valid and updated"; \
    return 0; }

    log "Running command: bcftools view --threads ${threads} --force-samples -s ${fam_members::-1} -Oz -o ${output_vcf} ${input_vcf}"
    bcftools view --threads ${threads} --force-samples -s ${fam_members::-1} -Oz -o ${output_vcf} ${input_vcf} && \
    mv ${output_vcf} ${input_vcf/.vcf*/.${family_name}.vcf.gz} && \
    tabix -p vcf -f ${input_vcf/.vcf*/.${family_name}.vcf.gz} && \
    display_vcf ${input_vcf/.vcf*/.${family_name}.vcf.gz}
}


function main_filtration () {
    local OPTIND v f c
    while getopts v:f::c: args
    do
        case ${args} in
            v) local input_vcf=$OPTARG ;;
            f) local fam_name=$OPTARG ;;
            c) local config=$OPTARG ;;
            *) echo "No argument passed. At least pass an argument specifying the family ID"
        esac
    done

    local ped_file=$(read_yaml ${config} "ped_file")
    local threads=$SNAKEMAKE_THREADS
    [[ -z ${threads} ]] && threads=$(read_yaml ${config} "threads")

    # Extract the family samples from the ped vcf file
    if [[ -f ${ped_file} ]]; then
        # Optional filtration based on family and pedigree information
        extract_fam_vcf ${input_vcf} ${ped_file} ${fam_name} ${threads} && \
        local fam_vcf=${input_vcf/.vcf*/.${fam_name}.vcf.gz} || \
        { log "Failed to extract the family VCF"; return 1; }
    fi

    # Filter the variants based on the pedigree information\
    if [[ -f ${ped_file} ]]; then
        # Optional operation based on pedigree information
        filter_allele_based_on_pedigree_with_py \
        -i ${fam_vcf} \
        -p ${ped_file} \
        -f ${fam_name} \
        -o ${fam_vcf/.vcf*/.filtered.vcf.gz} && \
        local filtered_vcf=${fam_vcf/.vcf*/.filtered.vcf.gz} || \
        { log "Failed to filter the variants based on the pedigree information"; return 1; }
    else
        log "No pedigree table and family name provided, directly add filtered tag to the input VCF file as the output VCF file"
        local filtered_vcf=${input_vcf/.vcf*/.filtered.vcf.gz}
    fi

    # Filter the variants based on the allele frequency
    if [[ -f ${ped_file} ]]; then
        filter_af \
        -c ${config} \
        -i ${filtered_vcf} \
        -t ${threads} || \
        { log "Failed to filter the variants based on the allele frequency"; return 1; }
    else
        log "No pedigree table and family name provided, directly filter the variants based on the allele frequency to get the filtered VCF file"
        filter_af \
        -c ${config} \
        -i ${input_vcf} \
        -o ${filtered_vcf} \
        -t ${threads} || \
        { log "Failed to filter the variants based on the allele frequency"; return 1; }
    fi
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

