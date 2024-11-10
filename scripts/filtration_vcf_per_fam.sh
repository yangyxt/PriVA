#! /usr/bin/env bash

# This script is used to filter the variants in the VCF file for each family
# The script expect a Family VCF input along with a pedigree table contains the Family pedigree information.
# Notice that the family sample ID should be consistent between the pedigree table and the VCF file.
# Notice that the proband of the family should stay as the first record of the family in the pedigree table

# Maintainer: yangyxt@gmail.com, yangyxt@hku.hk



function filter_allele_based_on_pedigree_with_py {
    local OPTIND i o p f
    while getopts i:o::p:f:: args
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
        local output_vcf=${input_vcf/.vcf/.rmfit.vcf}
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
    -o ${output_vcf}
    check_return_code

	display_vcf ${output_vcf}
    log "Finish filtering the records where control sample has homozygous or hemizygous GTs or the records no patients carrying the variant allele."$'\n\n'
}

# Specifically, we only filter out variants carried by healthy parents in homozygous form
local ped_filter_vcf=${pre_anno_vcf/.vcf/.ped.vcf}
filter_allele_based_on_pedigree_with_py \
-i ${pre_anno_vcf} \
-o ${ped_filter_vcf} \
-p ${ped_file} \
-f ${fam_name} && \
log "Successfully filter on pedigree information on ${pre_anno_vcf}. The result is ${ped_filter_vcf}" && \
display_vcf ${ped_filter_vcf} || { \
log "Failed to filter on pedigree information on ${pre_anno_vcf}. Quit now"
return 1; }


