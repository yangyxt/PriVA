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


function filter_allele_based_on_pedigree_data {
    local OPTIND i o p f
    while getopts i:o::p::f:: args
    do
        case ${args} in
            i) local input_vcf=$OPTARG ;;
            p) local pedfile=$OPTARG ;;
            o) local output_vcf=$OPTARG ;;
            f) local family_name=$OPTARG ;;
            *) echo "No argument passed. At least pass an argument specifying the family ID" ;;
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
       [[ $(count_vcf_records ${output_vcf}) -le $(count_vcf_records ${input_vcf}) ]]; then
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
    local output_vcf
    local threads

    # Use getopts to parse the input arguments
    local OPTIND i o c t
    while getopts i:o::c:t:: args
    do
        case ${args} in
            i) input_vcf=$OPTARG ;;
            o) output_vcf=$OPTARG ;;
            c) config_file=$OPTARG ;;
            t) threads=$OPTARG ;;
            *) echo "No argument passed. At least pass an argument specifying the family ID" ;;
        esac
    done

    [[ -z ${threads} ]] && threads=1 || log "Using ${threads} threads to filter the variants based on the allele frequency"
    [[ ! -f ${config_file} ]] && { log "The config file ${config_file} does not exist"; return 1; }

    local af_cutoff=$(read_yaml ${config_file} "af_cutoff")

    [[ -z ${output_vcf} ]] && \
    local tmp_tag=$(randomID) && \
    local output_vcf=${input_vcf/.vcf*/.${tmp_tag}.vcf.gz} && \
    local replace="TRUE" || \
    { log "Failed to generate a temporary output VCF file name. Quit with error"; return 1; }

    check_vcf_validity ${input_vcf} && \
    [[ $(bcftools view -h ${input_vcf} | grep -c "AF_grpmax_joint > ${af_cutoff}") -gt 0 ]] && \
    log "The input VCF ${input_vcf} already been filtered on allele frequency at the cutoff ${af_cutoff}" && \
    return 0

    bcftools view --threads ${threads} -e "AF_grpmax_joint > ${af_cutoff}" -Ou ${input_vcf} | \
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
            c) local config_file=$OPTARG ;;
            *) echo "No argument passed. At least pass an argument specifying the family ID" ;;
        esac
    done

    log "The input config file is ${config_file}"

    local ped_file=$(read_yaml ${config_file} "ped_file")
    local threads=$(read_yaml ${config_file} "threads_per_fam")

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
        filter_allele_based_on_pedigree_data \
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
        log "Filter the variants based on the allele frequency, the config file is ${config_file}"
        filter_af \
        -c ${config_file} \
        -i ${filtered_vcf} \
        -t ${threads} || \
        { log "Failed to filter the variants based on the allele frequency"; return 1; }
    else
        log "No pedigree table and family name provided, directly filter the variants based on the allele frequency to get the filtered VCF file"
        filter_af \
        -c ${config_file} \
        -i ${input_vcf} \
        -o ${filtered_vcf} \
        -t ${threads} || \
        { log "Failed to filter the variants based on the allele frequency"; return 1; }
    fi
    
    # Step 2: Filter CADD table based on filtered VCF to reduce memory usage in prioritization
    # Generate family-specific CADD file that follows same naming as VCF
    if [[ -f ${filtered_vcf} ]]; then
        # Derive CADD file names based on VCF naming pattern
        local base_name=$(basename ${filtered_vcf} .filtered.vcf.gz)  # e.g., INPUT_BASE.anno.family1
        local family_specific_cadd=$(dirname ${filtered_vcf})/"${base_name}.filtered.cadd.tsv"
        
        # Get the original CADD file from config for filtering source
        local original_cadd=$(read_yaml ${config_file} "cadd_output_file")
        
        if [[ -f ${original_cadd} ]]; then
            log "Step 2: Filtering CADD table based on filtered variants to optimize memory usage"
            log "Creating family-specific CADD file: ${family_specific_cadd}"
            
            filter_cadd_based_on_vcf \
            -i ${filtered_vcf} \
            -c ${original_cadd} \
            -o ${family_specific_cadd} \
            -t ${threads} || \
            { log "Failed to filter CADD table based on VCF variants"; return 1; }
            
            log "CADD filtering complete. Family-specific CADD file: ${family_specific_cadd}"
        else
            log "Skipping CADD filtering: Original CADD file (${original_cadd}) not found"
        fi
    else
        log "Skipping CADD filtering: Filtered VCF (${filtered_vcf}) not found"
    fi
}


function filter_cadd_based_on_vcf () {
    local OPTIND i c o t
    while getopts i:c:o::t:: args
    do
        case ${args} in
            i) local input_vcf=$OPTARG ;;
            c) local cadd_table=$OPTARG ;;
            o) local output_cadd=$OPTARG ;;
            t) local threads=$OPTARG ;;
            *) echo "No argument passed. At least pass VCF input and CADD table arguments" ;;
        esac
    done

    [[ -z ${threads} ]] && threads=1
    [[ -z ${output_cadd} ]] && output_cadd=${cadd_table/.tsv/.filtered.tsv}

    # Check if output is already up to date
    [[ -f ${output_cadd} ]] && \
    [[ ${output_cadd} -nt ${input_vcf} ]] && \
    [[ ${output_cadd} -nt ${cadd_table} ]] && \
    log "The filtered CADD table ${output_cadd} is up to date, skip filtering" && \
    return 0

    log "Filtering CADD table ${cadd_table} based on variants in ${input_vcf}"
    
    # Create temporary files
    local tmp_dir=$(mktemp -d)
    local vcf_variants="${tmp_dir}/vcf_variants.txt"
    local filtered_cadd="${tmp_dir}/filtered_cadd.tsv"
    
    # Extract variant positions from VCF (remove chr prefix to match CADD format)
    log "Extracting variant positions from VCF..."
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${input_vcf} | \
    sed 's/^chr//' > ${vcf_variants}
    
    local variant_count=$(wc -l < ${vcf_variants})
    log "Found ${variant_count} variants in VCF file"
    
    # Get header from CADD table
    head -n 1 ${cadd_table} > ${filtered_cadd}
    
    # Filter CADD table using mawk for efficiency
    log "Filtering CADD table (this may take a few minutes)..."
    
    # Create an associative array lookup in mawk
    mawk -v vcf_file="${vcf_variants}" '
    BEGIN {
        # Load VCF variants into associative array for fast lookup
        while ((getline line < vcf_file) > 0) {
            split(line, fields, "\t")
            chrom = fields[1]
            pos = fields[2] 
            ref = fields[3]
            alt = fields[4]
            key = chrom "\t" pos "\t" ref "\t" alt
            vcf_variants[key] = 1
        }
        close(vcf_file)
        print "Loaded", length(vcf_variants), "variants from VCF" > "/dev/stderr"
    }
    NR > 1 {  # Skip header (already written)
        # Create key from CADD table row
        key = $1 "\t" $2 "\t" $3 "\t" $4
        if (key in vcf_variants) {
            print $0
        }
    }' ${cadd_table} >> ${filtered_cadd}
    
    # Move filtered result to final location
    mv ${filtered_cadd} ${output_cadd}
    
    # Clean up temporary directory
    rm -rf ${tmp_dir}
    
    # Report results
    local original_lines=$(wc -l < ${cadd_table})
    local filtered_lines=$(wc -l < ${output_cadd})
    local reduction_percent=$(echo "scale=1; (${original_lines} - ${filtered_lines}) * 100 / ${original_lines}" | bc -l)
    
    log "CADD table filtering complete:"
    log "  Original lines: ${original_lines}"
    log "  Filtered lines: ${filtered_lines}" 
    log "  Reduction: ${reduction_percent}%"
    log "  Output: ${output_cadd}"
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

