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
    [[ ${output_tab} -nt ${SCRIPT_DIR}/combine_annotations.py ]] && \
    log "The combined annotation table ${output_tab} is up to date, skip the conversion" && \
    return 0 || \
    log "The combined annotation table ${output_tab} is outdated, start the conversion"

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
    local threads=${3}
    local splicing_py=${SCRIPT_DIR}/splicing_var_analysis.py

    local alphamissense_tranx_domain_map=$(read_yaml ${config_file} "alphamissense_tranx_domain_map")
    local intolerant_domains=$(read_yaml ${config_file} "all_intolerant_domains")
    local interpro_entry_map_pkl=$(read_yaml ${config_file} "interpro_mapping_pickle")
    
    [[ -f ${input_tab} ]] && \
    [[ ${input_tab} -nt ${interpro_entry_map_pkl} ]] && \
    [[ ${input_tab} -nt ${alphamissense_tranx_domain_map} ]] && \
    [[ ${input_tab} -nt ${intolerant_domains} ]] && \
    [[ ${input_tab} -nt ${splicing_py} ]] && \
    check_table_column ${input_tab} "splicing_lof" && \
    check_table_column ${input_tab} "splicing_len_changing" && \
    log "The input table ${input_tab} is up to date, skip the splicing interpretation" && \
    return 0 || \
    log "The input table ${input_tab} is outdated, start the splicing interpretation"

    python ${splicing_py} \
    --anno_table ${input_tab} \
    --transcript_domain_map ${alphamissense_tranx_domain_map} \
    --intolerant_domains ${intolerant_domains} \
    --interpro_entry_map_pkl ${interpro_entry_map_pkl} \
    --threads ${threads} && \
    display_table ${input_tab} && \
    log "The splicing interpretations are saved to ${input_tab}, added with two columns: splicing_lof and splicing_len_changing"
}



function interpret_utr_annotations () {
    local input_tab=${1}
    local config_file=${2}
    local threads=${3}

    local utr_py=${SCRIPT_DIR}/utr_anno_interpret.py

    local alphamissense_tranx_domain_map=$(read_yaml ${config_file} "alphamissense_tranx_domain_map")
    local intolerant_domains=$(read_yaml ${config_file} "all_intolerant_domains")

    [[ -f ${input_tab} ]] && \
    [[ ${input_tab} -nt ${alphamissense_tranx_domain_map} ]] && \
    [[ ${input_tab} -nt ${intolerant_domains} ]] && \
    [[ ${input_tab} -nt ${utr_py} ]] && \
    check_table_column ${input_tab} "5UTR_lof" && \
    check_table_column ${input_tab} "5UTR_length_changing" && \
    log "The input table ${input_tab} is up to date, skip the UTR interpretation" && \
    return 0 || \
    log "The input table ${input_tab} is outdated, start the UTR interpretation"

    python ${utr_py} \
    --variants_table ${input_tab} \
    --domain_map ${alphamissense_tranx_domain_map} \
    --intolerant_domains ${intolerant_domains} && \
    display_table ${input_tab} && \
    log "The UTR interpretations are saved to ${input_tab}, added with two columns: 5UTR_lof and 5UTR_length_changing"
}



function assign_acmg_criteria () {
    local input_tab
    local config_file
    local fam_name
    local threads

    # Use getopts to parse the input arguments based on the defined local variable above
    local OPTIND i c f t
    while getopts i:c:f::t: args
    do
        case ${args} in
            i) input_tab=$OPTARG ;;
            c) config_file=$OPTARG ;;
            f) fam_name=$OPTARG ;;
            t) threads=$OPTARG ;;
            *) echo "No argument passed. At least pass an argument specifying the input table" ;;
        esac
    done

    # Preprocess step, interpret the splicing annotations
    interpret_splicing_annotations ${input_tab} ${config_file} ${threads} || \
    { log "Failed to interpret splicing annotations for ${input_tab}"; return 1; }

    # Preprocess step, interpret the UTR annotations
    interpret_utr_annotations ${input_tab} ${config_file} ${threads} || \
    { log "Failed to interpret UTR annotations for ${input_tab}"; return 1; }
    
    local mean_am_score_table=${BASE_DIR}/data/alphamissense/alphamissense_mean_score.tsv
    local ped_table=$(read_yaml ${config_file} "ped_file")
    local clinvar_aa_dict_pkl=$(read_yaml ${config_file} "clinvar_aa_stat")
    local clinvar_splice_dict_pkl=$(read_yaml ${config_file} "clinvar_splice_stat")
    local intolerant_domains_pkl=$(read_yaml ${config_file} "all_intolerant_domains")
    local domain_mechanism_tsv=$(read_yaml ${config_file} "clinvar_intolerance_mechanisms")
    local tranx_exon_domain_map_pkl=$(read_yaml ${config_file} "alphamissense_tranx_domain_map")
    local intolerant_motifs_pkl=$(read_yaml ${config_file} "alphamissense_intolerant_motifs")
    local alt_disease_vcf=$(read_yaml ${config_file} "alt_disease_vcf")
    local gnomAD_extreme_rare_threshold=$(read_yaml ${config_file} "extreme_rare_PAF")
    local expected_incidence=$(read_yaml ${config_file} "exp_disease_incidence")
    local am_score_vcf=$(read_yaml ${config_file} "alphamissense_vcf")
    local clinvar_patho_af_stat=$(read_yaml ${config_file} "clinvar_patho_af_stat")
    local clinvar_patho_exon_af_stat=$(read_yaml ${config_file} "clinvar_patho_exon_af_stat")
    local interpro_entry_map_pkl=$(read_yaml ${config_file} "interpro_mapping_pickle")
    local output_acmg_mat=${input_tab::-4}.acmg.tsv
    local repeat_region_file_name=$(read_yaml ${config_file} "repeat_region_file_name")
    local repeat_region_file=${BASE_DIR}/data/repeats/${repeat_region_file_name}

    local has_error=0
    check_path ${clinvar_aa_dict_pkl} "file" "clinvar_aa_stat" || has_error=1
    check_path ${clinvar_splice_dict_pkl} "file" "clinvar_splice_stat" || has_error=1
    check_path ${intolerant_domains_pkl} "file" "all_intolerant_domains" || has_error=1
    check_path ${domain_mechanism_tsv} "file" "clinvar_intolerance_mechanisms" || has_error=1
    check_path ${intolerant_motifs_pkl} "file" "alphamissense_intolerant_motifs" || has_error=1
    check_path ${am_score_vcf} "file" "alphamissense_vcf" || has_error=1
    check_path ${tranx_exon_domain_map_pkl} "file" "alphamissense_tranx_domain_map" || has_error=1
    check_path ${clinvar_patho_af_stat} "file" "clinvar_patho_af_stat" || has_error=1
    check_path ${clinvar_patho_exon_af_stat} "file" "clinvar_patho_exon_af_stat" || has_error=1
    check_path ${interpro_entry_map_pkl} "file" "interpro_mapping_pickle" || has_error=1
    check_path ${repeat_region_file} "file" "repeat_region_file" || has_error=1

    local assembly=$(read_yaml ${config_file} "assembly")
    if [[ ${assembly} == "hg38" ]]; then
        local mavedb_metadata_tsv=${BASE_DIR}/data/MaveDB/mavedb_assay_summary.tsv
        check_path ${mavedb_metadata_tsv} "file" "mavedb_assay_summary" || has_error=1
        local mavedb_arg="--mavedb_metadata_tsv ${mavedb_metadata_tsv}"
    else
        local mavedb_arg=""
    fi

    [[ ${has_error} -eq 1 ]] && \
    { log "Failed to offer the valid required files for the ACMG criteria assignment"; return 1; }
    
    # Test whether the function can be skipped
    [[ -f ${input_tab} ]] && \
    [[ ${input_tab} -nt ${mean_am_score_table} ]] && \
    [[ ${input_tab} -nt ${clinvar_aa_dict_pkl} ]] && \
    [[ ${input_tab} -nt ${clinvar_splice_dict_pkl} ]] && \
    [[ ${input_tab} -nt ${intolerant_domains_pkl} ]] && \
    [[ ${input_tab} -nt ${domain_mechanism_tsv} ]] && \
    [[ ${input_tab} -nt ${alt_disease_vcf} ]] && \
    [[ ${input_tab} -nt ${am_score_vcf} ]] && \
    [[ ${input_tab} -nt ${clinvar_patho_af_stat} ]] && \
    [[ ${input_tab} -nt ${clinvar_patho_exon_af_stat} ]] && \
    [[ ${input_tab} -nt ${interpro_entry_map_pkl} ]] && \
    [[ ${input_tab} -nt ${tranx_exon_domain_map_pkl} ]] && \
    [[ ${input_tab} -ot ${output_acmg_mat} ]] && \
    [[ ${input_tab} -nt ${repeat_region_file} ]] && \
    check_table_column ${input_tab} "ACMG_quant_score" && \
    check_table_column ${input_tab} "ACMG_class" && \
    check_table_column ${input_tab} "ACMG_criteria" && \
    log "The ACMG criterias are already assigned for ${input_tab}, skip the assignment" && \
    return 0
    
    local acmg_py=${SCRIPT_DIR}/acmg_criteria_assign.py
    
    [[ -z ${ped_table} ]] && local ped_arg="" || local ped_arg="--ped_table ${ped_table}"
    [[ -z ${fam_name} ]] && local fam_arg="" || local fam_arg="--fam_name ${fam_name}"
    [[ -z ${alt_disease_vcf} ]] && local alt_disease_arg="" || local alt_disease_arg="--alt_disease_vcf ${alt_disease_vcf}"

    log "Running the following command to assign the ACMG criterias: python ${acmg_py} --anno_table ${input_tab} --am_score_table ${mean_am_score_table} --clinvar_aa_dict_pkl ${clinvar_aa_dict_pkl} --intolerant_domains_pkl ${intolerant_domains_pkl} --intolerant_motifs_pkl ${intolerant_motifs_pkl} --domain_mechanism_tsv ${domain_mechanism_tsv} --gnomAD_extreme_rare_threshold ${gnomAD_extreme_rare_threshold} --expected_incidence ${expected_incidence} --am_score_vcf ${am_score_vcf} --threads ${threads} --tranx_exon_domain_map_pkl ${tranx_exon_domain_map_pkl} ${ped_arg} ${fam_arg} ${alt_disease_arg} ${mavedb_arg}"
    python ${acmg_py} \
    --anno_table ${input_tab} \
    --am_score_table ${mean_am_score_table} \
    --clinvar_patho_af_stat ${clinvar_patho_af_stat} \
    --clinvar_patho_exon_af_stat ${clinvar_patho_exon_af_stat} \
    --clinvar_aa_dict_pkl ${clinvar_aa_dict_pkl} \
    --clinvar_splice_dict_pkl ${clinvar_splice_dict_pkl} \
    --interpro_entry_map_pkl ${interpro_entry_map_pkl} \
    --intolerant_domains_pkl ${intolerant_domains_pkl} \
    --intolerant_motifs_pkl ${intolerant_motifs_pkl} \
    --repeat_region_file ${repeat_region_file} \
    --domain_mechanism_tsv ${domain_mechanism_tsv} \
    --gnomAD_extreme_rare_threshold ${gnomAD_extreme_rare_threshold} \
    --expected_incidence ${expected_incidence} \
    --am_score_vcf ${am_score_vcf} \
    --threads ${threads} \
    --tranx_exon_domain_map_pkl ${tranx_exon_domain_map_pkl} ${ped_arg} ${fam_arg} ${alt_disease_arg} ${mavedb_arg} && \
    display_table ${input_tab} && \
    log "The ACMG criterias are assigned for ${input_tab}, added with three columns: ACMG_quant_score, ACMG_class, ACMG_criteria, and the output matrix is saved to ${output_acmg_mat}" && \
    display_table ${output_acmg_mat} || \
    { log "Failed to assign the ACMG criterias for ${input_tab}"; return 1; }
}



function main_prioritization () {
    local input_vcf
    local config_file
    local fam_name
    local cadd_tsv

    # Use getopts to parse the input arguments based on the defined local variable above
    local OPTIND i c f t
    while getopts i:c:f::t:: args
    do
        case ${args} in
            i) input_vcf=$OPTARG ;;
            c) config_file=$OPTARG ;;
            f) fam_name=$OPTARG ;;
            t) cadd_tsv=$OPTARG ;;
            *) echo "No argument passed. At least pass an argument specifying the input table" ;;
        esac
    done


    # Read the expected number of threads
    local threads=$SNAKEMAKE_THREADS
    [[ -z ${threads} ]] && threads=$(read_yaml ${config_file} "threads")
    log "The number of threads is set to ${threads} for the current family ${fam_name} and input VCF file ${input_vcf}"

    # Preprocess step, convert the annotated VCF to a Table with transcript specific annotations as rows
    if [[ -z ${cadd_tsv} ]]; then
        cadd_tsv=$(read_yaml ${config_file} "cadd_output_file")
    fi
    prepare_combined_tab \
    ${input_vcf} \
    ${cadd_tsv} \
    ${threads} && \
    local input_tab=${input_vcf/.vcf*/.tsv} && \
    [[ -f ${input_tab} ]] || \
    { log "Failed to prepare the combined annotation table"; return 1; }

    # Assign the ACMG criterias to each variant-transcript annotation record
    [[ -z ${fam_name} ]] && local fam_arg="" || local fam_arg="-f ${fam_name}"
    assign_acmg_criteria \
    -i ${input_tab} \
    -c ${config_file} \
    -t ${threads} ${fam_arg} || \
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


