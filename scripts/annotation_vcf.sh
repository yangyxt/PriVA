#!/usr/bin/env bash
# This script is used to define all the annotation process functions and run the main pipeline
# For now, the script is only used to process short variants (Indels and SNVs). Large structural variants and CNVs are not included.
# Since there will be unexpected indels in the input VCF file. So prescores are not enough to cover them. We need to use some tools to address this issue.
# In this workflow, we use VEP, CADD, SpliceAI, VEP-plugin(UTRannotator). We abandon the ANNOVAR for its related releasing and distribution restrictions.

# For the same set of configurations (arguments), the pipeline should start with the position it ends last time, unless user specifically asked to forcefully rerun the pipeline

# Maintainer: yangyxt@gmail.com, yangyxt@hku.hk

SELF_SCRIPT="$(realpath ${BASH_SOURCE[0]})"
SCRIPT_DIR="$(dirname "${SELF_SCRIPT}")"
BASE_DIR="$(dirname ${SCRIPT_DIR})"
if [[ ${BASE_DIR} == "/" ]]; then BASE_DIR=""; fi
DATA_DIR="${BASE_DIR}/data"

if [[ -z $TMPDIR ]]; then TMPDIR=/tmp; fi

# Source the other script
source "${SCRIPT_DIR}/common_bash_utils.sh"
# conda activate priva_acmg

log "The folder storing scripts is ${SCRIPT_DIR}, the base folder for used scripts and data is ${BASE_DIR}"



function main_workflow() {
    # Define local variables first
    local input_vcf \
          config \
          assembly \
          ref_genome \
          output_dir \
          threads \
          af_cutoff \
          gnomad_vcf_chrX \
          clinvar_vcf \
          vep_cache_dir \
          vep_plugins_dir \
          vep_plugins_cachedir \
          hub_vcf_file \
          hub_cadd_file

    # Source the argparse.bash script
    source ${SCRIPT_DIR}/argparse.bash || { log "Failed to source argparse.bash"; return 1; }

    argparse "$@" < ${SCRIPT_DIR}/anno_main_args || { log "Failed to parse arguments"; return 1; }

    # Process config file
    local config_file="${config:-${BASE_DIR}/config.yaml}"
    local -A config_args

    if [[ -f "$config_file" ]]; then
        log "Using config file: $config_file"
        local -a config_keys=(input_vcf assembly vep_cache_dir output_dir ref_genome threads af_cutoff gnomad_vcf_chrX clinvar_vcf vep_plugins_dir vep_plugins_cachedir hub_vcf_file hub_cadd_file control_vcf)
        for key in "${config_keys[@]}"; do
            # Need to remove the quotes around the value
            config_args[$key]="$(read_yaml ${config_file} ${key})"
        done
    else
        log "No config file found at $config_file, will use command line arguments only"
    fi

    # Set variables with command line arguments taking precedence over config file
    local input_vcf="${input_vcf:-${config_args[input_vcf]}}"

    # Set the optional arguments
    local assembly="${assembly:-${config_args[assembly]}}"
    local vep_cache_dir="${vep_cache_dir:-${config_args[vep_cache_dir]}}"
    local output_dir="${output_dir:-${config_args[output_dir]}}"
    local ref_genome="${ref_genome:-${config_args[ref_genome]}}"
    local threads="${threads:-${config_args[threads]}}"
    local af_cutoff="${af_cutoff:-${config_args[af_cutoff]:-0.05}}"
    local gnomad_vcf_chrX="${gnomad_vcf_chrX:-${config_args[gnomad_vcf_chrX]}}"
    local clinvar_vcf="${clinvar_vcf:-${config_args[clinvar_vcf]}}"
    local control_vcf="${control_vcf:-${config_args[control_vcf]}}"
    local vep_plugins_dir="${vep_plugins_dir:-${config_args[vep_plugins_dir]}}"
    local vep_plugins_cachedir="${vep_plugins_cachedir:-${config_args[vep_plugins_cachedir]}}"
    
    # Simplify hub caching settings, remove hub_enabled
    local hub_vcf_file="${hub_vcf_file:-${config_args[hub_vcf_file]}}"
    local hub_cadd_file="${hub_cadd_file:-${config_args[hub_cadd_file]}}"

    # Create hub directories if needed, simplified
    if [[ -n "${hub_vcf_file}" ]]; then
        local hub_dir=$(dirname "${hub_vcf_file}")
        if [[ ! -d "${hub_dir}" ]]; then
            mkdir -p "${hub_dir}" || {
                log "Error: Failed to create hub directory: ${hub_dir}"
                has_error=1
            }
        fi
    fi

    if [[ -n "${hub_cadd_file}" ]]; then
        local cadd_dir=$(dirname "${hub_cadd_file}")
        if [[ ! -d "${cadd_dir}" ]]; then
            mkdir -p "${cadd_dir}" || {
                log "Error: Failed to create CADD hub directory: ${cadd_dir}"
                has_error=1
            }
        fi
    fi

    # If assembly not specified, try to extract it from the input VCF
    [[ -z ${assembly} ]] && assembly=$(check_vcf_assembly_version ${input_vcf})
    [[ -z ${assembly} ]] && assembly=$(extract_assembly_from_fasta ${ref_genome})

    # Exit if any errors were found
    if [[ $has_error -eq 1 ]]; then
        return 1
    fi

    # Log the parsed arguments
    log "Using the following parameters:"
    log "  input_vcf: $input_vcf"
    log "  assembly: $assembly"
    log "  vep_cache_dir: $vep_cache_dir"
    log "  output_dir: $output_dir"
    log "  ref_genome: $ref_genome"
    log "  threads: $threads"
    log "  af_cutoff: $af_cutoff"
    log "  gnomad_vcf_chrX: $gnomad_vcf_chrX"
    log "  clinvar_vcf: $clinvar_vcf"
    log "  vep_plugins_dir: $vep_plugins_dir"
    log "  vep_plugins_cachedir: $vep_plugins_cachedir"
    log "  hub_vcf_file: $hub_vcf_file"
    log "  hub_cadd_file: $hub_cadd_file"
    
    # Prepare a comma separated list of samples in the input vcf
    local samples=$(bcftools query -l "${input_vcf}" | tr '\n' ',')
    log "The samples in the input vcf are ${samples}"

    # Preprocess the input vcf to:
    # Remove the variants not located in primary chromsomes
    # Convert the contig names to UCSC style. Meaning mapping VCF to hg19 or hg38
    # Sort and normalize (including indel left alignment)
    # Remove variants where proband DP < 5
    # Remove VCF records where all patients in the family has either missing or homo ref genotype.

    local final_anno_vcf=${output_dir}/$(basename ${input_vcf/.vcf*/.anno.vcf.gz}) && \
	local anno_vcf=${input_vcf/.vcf*/.anno.vcf.gz} && \
	ls -lhtr ${anno_vcf} || log "Expect the intermediate vcf file at ${anno_vcf} and the final output vcf file at ${final_anno_vcf}"
    preprocess_vcf \
    -i ${input_vcf} \
    -o ${anno_vcf} \
    -r ${ref_genome} \
    -t ${threads} && \
    log "Successfully preprocess the input vcf ${input_vcf} for annotation. The result is ${anno_vcf}" && \
    display_vcf ${anno_vcf} || { \
    log "Failed to preprocess the input vcf ${input_vcf} for annotation. Quit with error."; \
    return 1; }

    # If hub VCF file exists, check for variants
    if [[ -f "${hub_vcf_file}" ]] && check_vcf_validity "${hub_vcf_file}"; then
        local covered_vcf="${anno_vcf/.vcf*/.covered.vcf.gz}"
        local uncovered_vcf="${anno_vcf/.vcf*/.uncovered.vcf.gz}"
        
        log "Checking for cached annotations in hub VCF"
        local total_count=$(count_vcf_records "${anno_vcf}")
        log "The total number of variants in the input vcf is ${total_count}"
        local tmp_dir=$(read_yaml ${config_file} "tmp_dir")
		log "Using command to find cached annotation records: use_hub_vcf_annotations -i ${anno_vcf} -h ${hub_vcf_file} -o ${covered_vcf} -u ${uncovered_vcf} -t ${threads} -c ${total_count} -d ${tmp_dir}"
        if use_hub_vcf_annotations -i "${anno_vcf}" -h "${hub_vcf_file}" -o "${covered_vcf}" -u "${uncovered_vcf}" -t "${threads}" -c "${total_count}" -d "${tmp_dir}"; then
            log "All variants found in hub VCF, skipping annotation pipeline"

            if check_vcf_validity "${covered_vcf}" && \
               clean_vcf_multiallelics "${covered_vcf}" "${ref_genome}" "${threads}" && \
               [[ $(count_vcf_records "${covered_vcf}") -ge ${total_count} ]] && \
               [[ ${covered_vcf} -nt ${anno_vcf} ]] && \
               [[ "$(bcftools query -l ${covered_vcf})" == "$(bcftools query -l ${anno_vcf})" ]]; then
                log "The covered vcf file ${covered_vcf} contains all variants in the input vcf ${input_vcf}."
            
                # Use covered variants as final result
                cp "${covered_vcf}" "${anno_vcf}"
                cp "${covered_vcf}.tbi" "${anno_vcf}.tbi"
                
                # Skip to CADD annotation
                local skip_annotation=true
            fi
        else
            if check_vcf_validity "${covered_vcf}" && \
               clean_vcf_multiallelics "${covered_vcf}" "${ref_genome}" "${threads}" && \
               [[ ${covered_vcf} -nt ${anno_vcf} ]] && \
               [[ "$(bcftools query -l ${covered_vcf})" == "$(bcftools query -l ${anno_vcf})" ]] && \
               [[ "$(bcftools query -l ${uncovered_vcf})" == "$(bcftools query -l ${anno_vcf})" ]] && \
               [[ ${uncovered_vcf} -nt ${anno_vcf} ]] && \
               check_vcf_validity "${uncovered_vcf}" && \
               clean_vcf_multiallelics "${uncovered_vcf}" "${ref_genome}" "${threads}"; then
                log "Found cached variants and uncovered variants, will only annotate uncovered variants"
                
                # Save the original anno_vcf path and work with uncovered variants
                local tmp_anno_vcf="${anno_vcf}"
                local anno_vcf="${uncovered_vcf}"
            else
                log "No cached variants found, will annotate all variants"
            fi
        fi
    fi

    # Skip annotation if all variants were in hub VCF
    if [[ "${skip_annotation}" != "true" ]]; then
        # First annotate gnomAD aggregated frequency and number of homozygous ALT allele carriers
        # Tag the variants with gnomAD_common and gnomAD_BA FILTER tags
        anno_agg_gnomAD_data \
        ${anno_vcf} \
        ${threads} \
        ${assembly} \
        ${gnomad_vcf_chrX} && \
        log "Successfully add aggregated gnomAD annotation on ${anno_vcf}. The result is ${anno_vcf}" || { \
        log "Failed to add aggregated gnomAD annotation on ${anno_vcf}. Quit now"
        return 1; }

        # Now we annotate ClinVar variants and return result as VCF file
        # Annotate the variants with
        anno_clinvar_data \
        ${anno_vcf} \
        ${clinvar_vcf} \
        ${threads} || { \
        log "Failed to add ClinVar annotation on ${anno_vcf}. Quit now"
        return 1; }

        anno_control_vcf_allele \
        ${anno_vcf} \
        ${control_vcf} || { \
        log "Failed to add control VCF allele annotation on ${anno_vcf}. Quit now"
        return 1; }

        # Now we annotate the variants with VEP
        anno_VEP_data \
        --input_vcf ${anno_vcf} \
        --config ${config_file} \
        --assembly ${assembly} \
        --ref_genome ${ref_genome} \
        --vep_cache_dir ${vep_cache_dir} \
        --vep_plugins_dir ${vep_plugins_dir} \
        --vep_plugins_cachedir ${vep_plugins_cachedir} \
        --threads ${threads} || { \
        log "Failed to add VEP annotation on ${anno_vcf}. Quit now"
        return 1; }
        
        # If using hub caching, merge covered and newly annotated variants
        if [[ -n "${tmp_anno_vcf}" ]]; then
            local newly_annotated_vcf="${anno_vcf}"
            # Update hub VCF with newly annotated variants
            update_hub_vcf "${tmp_anno_vcf}" "${hub_vcf_file}" "${threads}"

            merge_annotated_vcfs "${covered_vcf}" "${newly_annotated_vcf}" "${tmp_anno_vcf}"
            
            # Set anno_vcf back to final result for downstream steps
            local anno_vcf="${tmp_anno_vcf}"
        fi
    fi

    # Now we annotate the VCF with CADD
    if [[ -n "${hub_cadd_file}" ]]; then
        local cadd_covered_tsv="${anno_vcf/.vcf*/.cadd.covered.tsv}"
        local cadd_uncovered_vcf="${anno_vcf/.vcf*/.cadd.uncovered.vcf.gz}"
        
        log "Checking for cached CADD scores in ${hub_cadd_file} with command: find_cached_cadd_variants ${anno_vcf} ${hub_cadd_file} ${cadd_covered_tsv} ${cadd_uncovered_vcf}"
        if find_cached_cadd_variants "${anno_vcf}" "${hub_cadd_file}" "${cadd_covered_tsv}" "${cadd_uncovered_vcf}" "${threads}"; then
            # All variants are covered - copy to output
            cp "${cadd_covered_tsv}" "${anno_vcf/.vcf*/.cadd.tsv}"
            update_yaml "${config_file}" "cadd_output_file" "${anno_vcf/.vcf*/.cadd.tsv}"
            log "All variants have cached CADD scores, skipping CADD calculation"
        else
            # Some variants need CADD calculation
            if check_vcf_validity "${cadd_uncovered_vcf}" 1 ${samples}; then
                # Calculate CADD for uncovered variants
                Calculate_CADD "${cadd_uncovered_vcf}" "${config_file}" "${anno_vcf/.vcf*/.cadd.new.tsv}" || {
                    log "Failed to calculate CADD scores. Quit now"
                    return 1
                }
                
                # Merge covered and newly calculated scores
                merge_cadd_results "${cadd_covered_tsv}" "${anno_vcf/.vcf*/.cadd.new.tsv}" "${anno_vcf/.vcf*/.cadd.tsv}" && \
                update_hub_cadd "${anno_vcf/.vcf*/.cadd.new.tsv}" "${hub_cadd_file}" && \
                update_yaml "${config_file}" "cadd_output_file" "${anno_vcf/.vcf*/.cadd.tsv}" && \
				mv ${anno_vcf} ${final_anno_vcf} && \
				check_vcf_validity ${final_anno_vcf} || \
				{
					log "Failed to update the output vcf file with CADD scores. Quit now"
					return 1
				}
                
                # Clean up temporary files
                rm -f "${anno_vcf/.vcf*/.cadd.new.tsv}"
            else
                # Only covered variants (should not happen but just in case)
                log "There should be uncovered variants existed but the VCF file ${cadd_uncovered_vcf} is invalid. Quit with error"
                return 1
            fi
            
            # Clean up temporary files
            rm -f "${cadd_covered_tsv}" "${cadd_uncovered_vcf}" "${cadd_uncovered_vcf}.tbi"
        fi
    else
        # Original CADD calculation without caching
        Calculate_CADD "${anno_vcf}" "${config_file}" && \
		mv ${anno_vcf} ${final_anno_vcf} && \
		check_vcf_validity ${final_anno_vcf} || \
		{
            log "Failed to add CADD annotation on ${anno_vcf}. Quit now"
            return 1
        }
    fi
}


function preprocess_vcf() {
    local input_vcf
    local ref_genome
    local output_vcf
    local threads

    # Use getopts to parse the arguments
    local OPTIND=1
    while getopts "i:o:r:t::" opt; do
        case ${opt} in
            i) input_vcf=${OPTARG};;
            o) output_vcf=${OPTARG};;
            r) ref_genome=${OPTARG};;
            t) threads=${OPTARG};;
        esac
    done

    # Make sure several things
    # 1. The VCF file is having a chromosome notation with chr prefix (UCSC style instead of the NCBI style)
    # 2. The VCF file is using the same assembly as the ref genome fasta file
    # 3. The VCF file is sorted and normalized, no multi-allelics allowed which should be splitted into multiple records
    # 4. Temporarily do not deal with variant records located in alternative contigs

    check_path ${input_vcf} "file" "input_vcf" || return 1
    check_path ${ref_genome} "file" "ref_genome" || return 1

    local vcf_assembly=$(check_vcf_assembly_version ${input_vcf})

    if [[ -z ${threads} ]]; then
        threads=1
    fi

    if [[ -z ${ref_genome} ]]; then
        log "User does not specify the ref genome fasta file. Cant proceed now. Quit with Error!"
        return 1;
    elif [[ ! -z ${vcf_assembly} ]] && [[ ! ${ref_genome} =~ ${vcf_assembly} ]]; then
        log "User specified ref_genome fasta file ${ref_genome} seems not match with the VCF used assembly ${vcf_assembly}. Quit with Error"
        return 1;
    fi

	log "Expect the output vcf file at ${output_vcf}"

    # Test if output_vcf is already valid
    if [[ ${output_vcf} -nt ${input_vcf} ]] && \
        check_vcf_validity ${output_vcf} && \
        clean_vcf_multiallelics ${output_vcf} "${ref_genome}" "${threads}" && \
        [[ "$(bcftools query -l ${output_vcf})" == "$(bcftools query -l ${input_vcf})" ]]; then
        log "The ${output_vcf} is valid and udpated. Check the sizes of the output vcf and the input vcf ${input_vcf}"
        local output_vcf_size=$(ls -l ${output_vcf} | awk '{print $5}')
        local input_vcf_size=$(ls -l ${input_vcf} | awk '{print $5}')
        # If output is smaller than 80% (use awk to calculate the percentage) of the input, then we need to continue running the subsequent steps
        if [[ $(awk -v a=${output_vcf_size} -v b=${input_vcf_size} 'BEGIN{printf "%d", (a/b)*100}') -lt 80 ]]; then
            log "The output vcf ${output_vcf} $(ls -lh ${output_vcf} | cut -d ' ' -f 5) is smaller than 80% of the input vcf ${input_vcf} $(ls -lh ${input_vcf} | cut -d ' ' -f 5). Continue running the subsequent steps"
        else
            log "The output vcf ${output_vcf} $(ls -lh ${output_vcf} | cut -d ' ' -f 5) is larger than 80% of the input vcf ${input_vcf} $(ls -lh ${input_vcf} | cut -d ' ' -f 5). Skip the subsequent steps"
            return 0;
        fi        
    fi

    log "Start to preprocess the input vcf ${input_vcf} to make it: 1. Biallelic, 2. Use UCSC style chromosome notation, 3. Only contain variants in primary chromosomes, 5. Sorted and normalized (left align indels)"
    display_vcf ${input_vcf}

    # First filter out variants that not in primary chromosomes (chrM to chrY)(MT to Y)
    bcftools view --threads ${threads} -r "$(cat ${BASE_DIR}/data/liftover/ucsc_GRC.primary.contigs.tsv | tr '\n' ',')" -Ou ${input_vcf}| \
    bcftools sort -Oz -o ${input_vcf/.vcf*/.primary.vcf.gz} && \
    tabix -f -p vcf ${input_vcf/.vcf*/.primary.vcf.gz} && \
    display_vcf ${input_vcf/.vcf*/.primary.vcf.gz} || \
    { log "Failed to generate a VCF file that only contains records in primary chromosomes"; \
      return 1; }

    # First check whether the VCF is using NCBI or UCSC assembly
    local assembly_version=$(check_vcf_contig_version ${input_vcf/.vcf*/.primary.vcf.gz})
    if [[ ${assembly_version} =~ "ncbi" ]]; then
        log "The input vcf is detected to map variants to GRC assemblies instead of UCSC assemblies" && \
        liftover_from_GRCh_to_ucsc \
        ${input_vcf/.vcf*/.primary.vcf.gz} \
        ${BASE_DIR}/data/liftover/ucsc_to_GRC.contig.map.tsv \
        ${input_vcf/.vcf*/.ucsc.vcf.gz}
    else
        cp -f ${input_vcf/.vcf*/.primary.vcf.gz} ${input_vcf/.vcf*/.ucsc.vcf.gz}
    fi

    # First sort the input_vcf,
    # Then normalize the input_vcf with bcftools
    normalize_vcf ${input_vcf/.vcf*/.ucsc.vcf.gz} ${input_vcf/.vcf*/.norm.vcf.gz} ${ref_genome} ${threads}
    announce_remove_tmps ${input_vcf/.vcf*/}*tmp*vcf*

    # Then we add uniq IDs to these variants
    prepare_vcf_add_varID \
    ${input_vcf/.vcf/.norm.vcf} \
    ${output_vcf} && \
    tabix -f -p vcf ${output_vcf} && \
    announce_remove_tmps ${input_vcf/.vcf/.norm.vcf} && \
    announce_remove_tmps ${input_vcf/.vcf*/.ucsc.vcf.gz} && \
    announce_remove_tmps ${input_vcf/.vcf*/.primary.vcf.gz} && \
    display_vcf ${output_vcf}
}




function anno_agg_gnomAD_data () {
    local input_vcf=${1}
    local threads=${2}
    local assembly=${3}
    local gnomad_vcf_chrX=${4}
    local tmp_tag=$(randomID)
    local output_vcf=${input_vcf/.vcf/.${tmp_tag}.vcf}

    # Check the compression format of the input vcf file
    [[ ! ${input_vcf} =~ \.vcf\.gz$ ]] && {
        log "The input vcf ${input_vcf} is not a gzipped vcf file. Quit with Error!"
        return 1;
    }

    # Core frequency fields
    local core_fields="CHROM,POS,REF,ALT,.INFO/AC_joint,.INFO/AN_joint,.INFO/AF_joint,.INFO/nhomalt_joint"

    # Maximum values across populations (no sex-specific versions available)
    local max_fields=",.INFO/AC_grpmax_joint,.INFO/AF_grpmax_joint,.INFO/AN_grpmax_joint,.INFO/nhomalt_grpmax_joint"

    # Sex-specific fields (overall)
    local sex_fields=",.INFO/AF_joint_XX,.INFO/AF_joint_XY,.INFO/nhomalt_joint_XX,.INFO/nhomalt_joint_XY"

    # Major population frequencies
    local -a pop_codes=("nfe" "eas" "afr" "amr" "asj" "fin" "sas" "mid" "remaining")
    local pop_af_fields=""
    local pop_nhomalt_fields=""
    for pop in ${pop_codes[@]}; do
        pop_af_fields+=",.INFO/AF_joint_${pop}_XX,.INFO/AF_joint_${pop}_XY"
        pop_nhomalt_fields+=",.INFO/nhomalt_joint_${pop}_XX"
    done

	check_vcf_validity ${input_vcf} || \
	{
		log "The input vcf ${input_vcf} is not valid. Quit now"
		return 1
	}

    check_vcf_infotags ${input_vcf} "${core_fields}${max_fields}${sex_fields}${pop_af_fields}${pop_nhomalt_fields}" && \
    log "The input vcf ${input_vcf} already contains the INFO tags ${core_fields}${max_fields}${sex_fields}${pop_af_fields}${pop_nhomalt_fields}. We do not need to add them again" && \
    return 0 || \
	log "The input vcf ${input_vcf} does not contain the INFO tags ${core_fields}${max_fields}${sex_fields}${pop_af_fields}${pop_nhomalt_fields}. We need to add them"

    # We already make sure the input VCF is sorted and normalized
    # We also make sure no variants in alternative contigs are included in the input VCF
    # The used gnomAD vcf files should be the ones with joint AF information from both exome and genome datasets

    local -a chr_chroms
    # Step 1: Extract all the chromosome names from the input vcf
    mapfile -t chr_chroms < <(bcftools query -f '%CHROM\n' ${input_vcf} | uniq -)
    log "The chromosome names extracted from the input vcf are: ${chr_chroms[*]}"

    # Step 2: Further split the main-contig-variants file to multiple single-contig variant files
    local -a valid_chroms
    for chr in "${chr_chroms[@]}"; do
        if [[ ${chr} =~ chr[0-9XY]+ ]]; then
            if [[ ! -f ${input_vcf/.vcf*/.${chr}.vcf.gz} ]] || [[ ${input_vcf/.vcf*/.${chr}.vcf.gz} -ot ${input_vcf} ]]; then
                bcftools view -r ${chr} -Ou ${input_vcf} | \
                bcftools sort -Oz -o ${input_vcf/.vcf*/.${chr}.vcf.gz}
            fi
            if [[ ! -f ${input_vcf/.vcf*/.${chr}.vcf.gz.tbi} ]] || [[ ${input_vcf/.vcf*/.${chr}.vcf.gz.tbi} -ot ${input_vcf/.vcf*/.${chr}.vcf.gz} ]]; then
                tabix -f -p vcf ${input_vcf/.vcf*/.${chr}.vcf.gz}
            fi
            check_vcf_validity ${input_vcf/.vcf*/.${chr}.vcf.gz} && \
            valid_chroms+=( ${chr} )
        else
            log "WARNING!!! The chromosome ${chr} is not supported by gnomAD now. Currently only support chr1 to chrY. Skip this chromosome"
        fi
    done
    log "The valid chromosomes are: ${valid_chroms[*]}"

    # Step 3: Perform annotation with bcftools annotate to add the INFO fields from gnomAD vcfs to the input splitted vcfs
    export gnomad_vcf_chrX
    local gnomad_vcf_suffix=${gnomad_vcf_chrX/*.chrX.vcf/}
    parallel -j ${threads} --dry-run "if [[ ! -f ${input_vcf/.vcf*/}.{}.gnomAD.vcf.gz ]] || [[ ${input_vcf/.vcf*/}.{}.gnomAD.vcf.gz -ot ${input_vcf/.vcf*/}.{}.vcf.gz ]] || [[ ${input_vcf/.vcf*/}.{}.gnomAD.vcf.gz -ot ${gnomad_vcf_chrX/.chrX.vcf*gz/}.{}.vcf${gnomad_vcf_suffix} ]]; then bcftools annotate -a ${gnomad_vcf_chrX/.chrX.vcf*gz/}.{}.vcf${gnomad_vcf_suffix} -c ${core_fields}${max_fields}${sex_fields}${pop_af_fields}${pop_nhomalt_fields} -Oz -o ${input_vcf/.vcf*/}.{}.gnomAD.vcf.gz ${input_vcf/.vcf*/}.{}.vcf.gz; fi" ::: "${valid_chroms[@]}" && \
    parallel -j ${threads} --joblog ${input_vcf/.vcf*/.gnomAD.log} "if [[ ! -f ${input_vcf/.vcf*/}.{}.gnomAD.vcf.gz ]] || [[ ${input_vcf/.vcf*/}.{}.gnomAD.vcf.gz -ot ${input_vcf/.vcf*/}.{}.vcf.gz ]] || [[ ${input_vcf/.vcf*/}.{}.gnomAD.vcf.gz -ot ${gnomad_vcf_chrX/.chrX.vcf*gz/}.{}.vcf${gnomad_vcf_suffix} ]]; then bcftools annotate -a ${gnomad_vcf_chrX/.chrX.vcf*gz/}.{}.vcf${gnomad_vcf_suffix} -c ${core_fields}${max_fields}${sex_fields}${pop_af_fields}${pop_nhomalt_fields} -Oz -o ${input_vcf/.vcf*/}.{}.gnomAD.vcf.gz ${input_vcf/.vcf*/}.{}.vcf.gz; fi" ::: "${valid_chroms[@]}"; \
    check_parallel_joblog ${input_vcf/.vcf*/.gnomAD.log} || { \
    log "Failed to add aggregated gnomAD annotation on ${input_vcf}. Quit now"
    return 1; }

    # Merge the annotated vcfs together (including the alt-contig-variants file) to form into the output vcf
    # We can filter on AF here
    bcftools concat -Ou ${input_vcf/.vcf*/.chr*.gnomAD.vcf.gz} | \
    bcftools sort -Oz -o ${output_vcf} && \
    tabix -f -p vcf ${output_vcf} && \
    mv ${output_vcf} ${input_vcf} && \
    mv ${output_vcf}.tbi ${input_vcf}.tbi && \
    rm -f ${input_vcf/.vcf*/.chr*.gnomAD.vcf.gz} && \
    rm -f ${input_vcf/.vcf*/.chr*.gnomAD.vcf.gz.tbi}

}



function anno_clinvar_data () {
    local input_vcf=${1}
    local clinvar_vcf=${2}
    local threads=${3}

    if [[ -z ${threads} ]]; then
        threads=1
    fi

    local tmp_tag=$(randomID)
    local output_vcf=${input_vcf/.vcf/.${tmp_tag}.vcf}

    [[ ! ${input_vcf} =~ \.vcf\.gz$ ]] && {
        log "The input vcf ${input_vcf} is not a gzipped vcf file. Quit with Error!"
        return 1;
    }

    check_vcf_infotags ${input_vcf} "CLNDN,CLNHGVS,CLNREVSTAT,CLNSIG,GENEINFO,CLNCSQ" && \
    log "The input vcf ${input_vcf} already contains the INFO tags CLNDN,CLNHGVS,CLNREVSTAT,CLNSIG,GENEINFO,CLNCSQ. We do not need to add them again" && \
    return 0 || \
    log "The input vcf ${input_vcf} does not contain the INFO tags CLNDN,CLNHGVS,CLNREVSTAT,CLNSIG,GENEINFO,CLNCSQ. We need to add them"

    bcftools annotate --threads ${threads} -a ${clinvar_vcf} -c CHROM,POS,REF,ALT,.INFO/CLNDN,.INFO/CLNHGVS,.INFO/CLNREVSTAT,.INFO/CLNSIG,.INFO/GENEINFO,.INFO/CLNCSQ:=INFO/CSQ -Ou ${input_vcf} | \
    bcftools sort -Oz -o ${output_vcf} && \
    tabix -f -p vcf ${output_vcf} && \
    mv ${output_vcf} ${input_vcf} && \
    mv ${output_vcf}.tbi ${input_vcf}.tbi && \
    check_vcf_infotags ${input_vcf} "CLNDN,CLNHGVS,CLNREVSTAT,CLNSIG,GENEINFO,CLNCSQ" && \
    display_vcf ${input_vcf} || \
    { log "Failed to add ClinVar annotation on ${input_vcf}. Quit now"; return 1; }
}



function anno_control_vcf_allele() {
    local input_vcf=$1
    local control_vcf=$2
    local tmp_tag=$(randomID)
    # Output will initially be temporary, then moved to replace input_vcf
    local output_vcf=${input_vcf/.vcf.gz/.${tmp_tag}.vcf.gz}
    local processed_control_vcf="" # Will hold path to control vcf with required tags

    # --- Validate Inputs ---
    check_path "${input_vcf}" "file" "input_vcf" || { log "ERROR: input_vcf is not a valid file"; return 1; }
    check_path "${control_vcf}" "file" "control_vcf" || { log "WARNING: control_vcf is not a valid file but this function is optional"; return 0; }
    [[ "${input_vcf}" =~ \.vcf\.gz$ ]] || { log "ERROR: input_vcf must be bgzipped (.vcf.gz)"; return 1; }
    [[ "${control_vcf}" =~ \.vcf\.gz$ ]] || { log "ERROR: control_vcf must be bgzipped (.vcf.gz)"; return 1; }
    check_vcf_validity "${input_vcf}" || return 1
    check_vcf_validity "${control_vcf}" || return 1

    # --- Check if annotation already done ---
    local control_tags="control_AC,control_AN,control_AF,control_nhomalt"
    if check_vcf_infotags "${input_vcf}" "${control_tags}"; then
        log "Input VCF ${input_vcf} already contains control allele info tags: ${control_tags}. Skipping."
        return 0
    fi

    # --- Prepare Control VCF ---
    # Check for AC, AN, AF, and AC_Hom (assuming AC_Hom includes relevant hemizygous cases formatted as 1/1)
    local source_tags="AC,AN,AF,AC_Hom"
    if check_vcf_infotags "${control_vcf}" "${source_tags}"; then
        log "Control VCF ${control_vcf} already has required tags (${source_tags})."
        processed_control_vcf="${control_vcf}"
    else
        log "Control VCF ${control_vcf} missing one or more required tags (${source_tags}). Calculating..."
        local temp_control_vcf=$(mktemp --tmpdir="$TMPDIR" control_processed.XXXXXX.vcf.gz)
        announce_remove_tmps "${temp_control_vcf}" "${temp_control_vcf}.tbi" # Schedule cleanup

        # Ensure control VCF is indexed if needed
        if [[ ! -f "${control_vcf}.tbi" ]]; then
           log "Indexing control VCF: ${control_vcf}"
           tabix -p vcf -f "${control_vcf}" || { log "ERROR: Failed to index control VCF ${control_vcf}"; return 1; }
        fi

        # Calculate the required tags
        bcftools +fill-tags "${control_vcf}" -Oz -o "${temp_control_vcf}" -- -t ${source_tags} && \
        tabix -p vcf -f "${temp_control_vcf}" || {
            log "ERROR: Failed to calculate tags for control VCF ${control_vcf}"
            return 1
        }
        log "Calculated tags written to temporary file: ${temp_control_vcf}"
        processed_control_vcf="${temp_control_vcf}"
    fi

    # --- Annotate Input VCF ---
    log "Annotating ${input_vcf} with control allele info from ${processed_control_vcf}"
    bcftools annotate \
        -a "${processed_control_vcf}" \
        -c "CHROM,POS,REF,ALT,INFO/control_AC:=INFO/AC,INFO/control_AN:=INFO/AN,INFO/control_AF:=INFO/AF,INFO/control_nhomalt:=INFO/AC_Hom" \
        -Oz -o "${output_vcf}" \
        "${input_vcf}" && \
    tabix -p vcf -f "${output_vcf}" || {
        log "ERROR: Failed to annotate ${input_vcf} with control allele info."
        # Explicitly remove temp file on error if announce_remove_tmps hasn't run yet
        [[ "${processed_control_vcf}" != "${control_vcf}" ]] && rm -f "${processed_control_vcf}" "${processed_control_vcf}.tbi"
        rm -f "${output_vcf}" "${output_vcf}.tbi" # Remove potentially failed output
        return 1
    }

    # --- Finalize ---
    log "Successfully annotated ${input_vcf}. Moving temporary output ${output_vcf} to replace input."
    mv "${output_vcf}" "${input_vcf}" && \
    mv "${output_vcf}.tbi" "${input_vcf}.tbi" || {
        log "ERROR: Failed to move temporary output ${output_vcf} to ${input_vcf}"
        # Explicitly remove temp file on error
        [[ "${processed_control_vcf}" != "${control_vcf}" ]] && rm -f "${processed_control_vcf}" "${processed_control_vcf}.tbi"
        rm -f "${output_vcf}" "${output_vcf}.tbi" # Clean up the temp output we couldn't move
        return 1
    }

    # Cleanup temp file if created
    [[ "${processed_control_vcf}" != "${control_vcf}" ]] && rm -f "${processed_control_vcf}" "${processed_control_vcf}.tbi"

    display_vcf "${input_vcf}"
    log "Annotation with control cohort allele frequencies complete for ${input_vcf}."
    return 0
}




function prepare_vcf_add_varID {
    local input_vcf=${1}
    local output_vcf=${2}

    if [[ -z ${output_vcf} ]]; then
        local output_vcf=${input_vcf/.vcf/.id.vcf}
    fi

    python3 ${SCRIPT_DIR}/generate_varID_for_vcf.py \
    -m encode \
    -i ${input_vcf} \
    -o ${output_vcf}
    check_return_code

    display_vcf ${output_vcf}
}



function anno_VEP_data() {
    # Declare local variables first
    local input_vcf \
          config \
          assembly \
          ref_genome \
          vep_cache_dir \
          vep_plugins_dir \
          vep_plugins_cachedir \
          threads \
          utr_annotator_file \
          loeuf_prescore \
          alphamissense_prescore \
          spliceai_snv_prescore \
          spliceai_indel_prescore \
          primateai_prescore \
          conservation_file \
          loftee_conservation_file \
          gerp_bigwig \
          human_ancestor_fasta \
          splicevault_prescore

    # Process arguments using anno_vep_args
    source ${SCRIPT_DIR}/argparse.bash || { log "Failed to source argparse.bash"; return 1; }
    argparse "$@" < ${SCRIPT_DIR}/anno_vep_args || { log "Failed to parse arguments"; return 1; }

    # Process config file
    local config_file="${config:-${BASE_DIR}/config.yaml}"
    local -A config_args

    if [[ -f "$config_file" ]]; then
        log "Using config file: $config_file to assign default values to the following variables"
        local -a config_keys=(assembly ref_genome vep_cache_dir vep_plugins_dir vep_plugins_cachedir threads
                              utr_annotator_file loeuf_prescore alphamissense_prescore
                              spliceai_snv_prescore spliceai_indel_prescore splicevault_prescore primateai_prescore 
                              loftee_repo human_ancestor_fasta loftee_conservation_file gerp_bigwig 
                              conservation_file)
        for key in "${config_keys[@]}"; do
            # Need to strip the double quotes enclosing the string value
            config_args[$key]="$(read_yaml ${config_file} ${key})"
        done
    else
        log "No config file found at $config_file, will use command line arguments only"
    fi
    log "The splicevault_prescore is ${config_args[splicevault_prescore]}"

    # Set variables with command line arguments taking precedence over config file
    local assembly="${assembly:-${config_args[assembly]}}"
    local ref_genome="${ref_genome:-${config_args[ref_genome]}}"
    local vep_cache_dir="${vep_cache_dir:-${config_args[vep_cache_dir]}}"
    local vep_plugins_dir="${vep_plugins_dir:-${config_args[vep_plugins_dir]}}"
    local vep_plugins_cachedir="${vep_plugins_cachedir:-${config_args[vep_plugins_cachedir]}}"
    local threads="${threads:-${config_args[threads]:-4}}"

    # Plugin cache files
    local utr_annotator_file="${utr_annotator_file:-${config_args[utr_annotator_file]}}"
    local loeuf_prescore="${loeuf_prescore:-${config_args[loeuf_prescore]}}"
    local alphamissense_prescore="${alphamissense_prescore:-${config_args[alphamissense_prescore]}}"
    local spliceai_snv_prescore="${spliceai_snv_prescore:-${config_args[spliceai_snv_prescore]}}"
    local spliceai_indel_prescore="${spliceai_indel_prescore:-${config_args[spliceai_indel_prescore]}}"
    local primateai_prescore="${primateai_prescore:-${config_args[primateai_prescore]}}"
    local conservation_file="${conservation_file:-${config_args[conservation_file]}}"
    local splicevault_prescore="${splicevault_prescore:-${config_args[splicevault_prescore]}}"
    local loftee_repo="${loftee_repo:-${config_args[loftee_repo]}}"
    local human_ancestor_fasta="${human_ancestor_fasta:-${config_args[human_ancestor_fasta]}}"
    local loftee_conservation_file="${loftee_conservation_file:-${config_args[loftee_conservation_file]}}"
    local gerp_bigwig="${gerp_bigwig:-${config_args[gerp_bigwig]}}"


    # Validate inputs
    local has_error=0
    check_path "$input_vcf" "file" "input_vcf" || has_error=1
    check_path "$ref_genome" "file" "ref_genome" || has_error=1
    check_path "$vep_cache_dir" "dir" "vep_cache_dir" || has_error=1
    check_path "$vep_plugins_dir" "dir" "vep_plugins_dir" || has_error=1
    check_path "$vep_plugins_cachedir" "dir" "vep_plugins_cachedir" || has_error=1

    # Check plugin cache files
    check_path "$utr_annotator_file" "file" "utr_annotator_file" || has_error=1
    check_path "$loeuf_prescore" "file" "loeuf_prescore" || has_error=1
    check_path "$alphamissense_prescore" "file" "alphamissense_prescore" || has_error=1
    check_path "$spliceai_snv_prescore" "file" "spliceai_snv_prescore" || has_error=1
    check_path "$spliceai_indel_prescore" "file" "spliceai_indel_prescore" || has_error=1
    check_path "$primateai_prescore" "file" "primateai_prescore" || has_error=1
    check_path "$conservation_file" "file" "conservation_file" || has_error=1
    check_path "$splicevault_prescore" "file" "splicevault_prescore" || has_error=1
    check_path "$loftee_repo" "dir" "loftee_repo" || has_error=1
    check_path "$human_ancestor_fasta" "file" "human_ancestor_fasta" || has_error=1
    check_path "$loftee_conservation_file" "file" "loftee_conservation_file" || has_error=1

    # Try to determine assembly if not specified
    [[ -z ${assembly} ]] && assembly=$(check_vcf_assembly_version ${input_vcf})
    [[ -z ${assembly} ]] && assembly=$(extract_assembly_from_fasta ${ref_genome})
    [[ ${assembly} == "hg19" ]] && assembly="GRCh37" # Convert to NCBI assembly name
    [[ ${assembly} == "hg38" ]] && assembly="GRCh38" # Convert to NCBI assembly name

    if [[ ${assembly} == "GRCh37" ]]; then
        # VEP does not provide gerp_bigwig for GRCh37
        local loftee_gerp_bw_argument=""
        local loftee_conservation_argument=",conservation_file:${loftee_conservation_file}"
        local mavedb_argument=""
    elif [[ ${assembly} == "GRCh38" ]]; then
        local loftee_gerp_bw_argument=",gerp_bigwig:${gerp_bigwig}"
        check_path "$gerp_bigwig" "file" "gerp_bigwig" || has_error=1
        # Conservation file (SQL) is broken in GRCh38
        local loftee_conservation_argument=""
        local mavedb_file=$(read_yaml ${config_file} "mavedb_file")
        check_path "$mavedb_file" "file" "mavedb_file" || has_error=1
        local mavedb_argument="-plugin MaveDB,file=${mavedb_file},cols=urn:score:nt:pro:se:sd:df:epsilon:vtype:high_conf:pvalue:ci95_lower:ci95_upper"
    fi

    # Exit if any errors were found
    if [[ $has_error -gt 0 ]]; then
        return 1
    fi

    # We can check that whether the input vcf is already annotated by VEP
    # The annotated INFO tags should include:
    # - CSQ, from VEP
    # - NMD, from VEP plugin NMD
    # - am_pathogenicity & am_class, from VEP plugin AlphaMissense
    # - SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL, SpliceAI_pred_DS_DG, SpliceAI_pred_DS_DL, SpliceAI_pred_DP_AG, SpliceAI_pred_DP_AL, SpliceAI_pred_DP_DG, SpliceAI_pred_DP_DL, from VEP plugin SpliceAI
    check_vcf_infotags ${input_vcf} "CSQ" && \
    [[ ${input_vcf} -nt ${SELF_SCRIPT} ]] && \
    log "The input vcf ${input_vcf} is already annotated by VEP. Skip this step" && \
    return 0

    local tmp_tag=$(randomID)
    local output_vcf=${input_vcf/.vcf*/.${tmp_tag}.vcf}

    log "Running VEP annotation with the command below:"
    log "vep -i ${input_vcf} --format vcf --verbose --vcf --species homo_sapiens --use_transcript_ref --assembly ${assembly} --cache --offline --merged --domains --hgvs --numbers --symbol --canonical --total_length --variant_class --gene_phenotype --stats_file ${input_vcf/.vcf*/.vep.stats.html} --fork ${threads} --buffer_size 10000 --fasta ${ref_genome} --dir_cache ${vep_cache_dir} --dir_plugins ${vep_plugins_dir} -plugin UTRAnnotator,file=${utr_annotator_file} -plugin LOEUF,file=${loeuf_prescore},match_by=transcript -plugin AlphaMissense,file=${alphamissense_prescore} -plugin LoF,loftee_path:${loftee_repo},human_ancestor_fa:${human_ancestor_fasta}${loftee_conservation_argument}${loftee_gerp_bw_argument} -plugin SpliceAI,snv=${spliceai_snv_prescore},indel=${spliceai_indel_prescore},cutoff=0.5 -plugin PrimateAI,${primateai_prescore} -plugin SpliceVault,file=${splicevault_prescore} -plugin Conservation,${conservation_file},MAX -plugin NMD ${mavedb_argument} --force_overwrite -o ${output_vcf}"

    # Run VEP annotation
    vep -i ${input_vcf} \
    --format vcf \
    --verbose \
    --vcf \
    --species homo_sapiens \
    --use_transcript_ref \
    --assembly ${assembly} \
    --cache \
    --offline \
    --merged \
    --domains \
    --hgvs \
    --numbers \
    --symbol \
    --canonical \
    --total_length \
    --variant_class \
    --gene_phenotype \
    --stats_file ${input_vcf/.vcf*/.vep.stats.html} \
    --fork ${threads} \
    --buffer_size 10000 \
    --fasta ${ref_genome} \
    --dir_cache ${vep_cache_dir} \
    --dir_plugins ${vep_plugins_dir},${loftee_repo} \
    -plugin UTRAnnotator,file=${utr_annotator_file} \
    -plugin LOEUF,file=${loeuf_prescore},match_by=transcript \
    -plugin AlphaMissense,file=${alphamissense_prescore} \
    -plugin LoF,loftee_path:${loftee_repo},human_ancestor_fa:${human_ancestor_fasta}${loftee_conservation_argument}${loftee_gerp_bw_argument} \
    -plugin SpliceAI,snv=${spliceai_snv_prescore},indel=${spliceai_indel_prescore},cutoff=0.5 \
    -plugin PrimateAI,${primateai_prescore} \
    -plugin SpliceVault,file=${splicevault_prescore} \
    -plugin Conservation,${conservation_file},MAX \
    -plugin NMD ${mavedb_argument} \
    --force_overwrite \
    -o ${output_vcf} && \
    bcftools sort -Oz -o ${input_vcf/.vcf*/.vcf.gz} ${output_vcf} && \
    tabix -f -p vcf ${input_vcf/.vcf*/.vcf.gz} && \
    announce_remove_tmps ${output_vcf} && \
    display_vcf ${input_vcf/.vcf*/.vcf.gz}
}



function Calculate_CADD {
    local input_vcf=${1}
    local config_file=${2}
    local output_file=${3}
    local threads=$(read_yaml ${config_file} "threads")
    local cadd_base_dir=$(read_yaml ${config_file} "cadd_base_dir")
    local assembly=$(read_yaml ${config_file} "assembly")
    local tmp_dir=$(read_yaml ${config_file} "tmp_dir")
    local cadd_script=${cadd_base_dir}/CADD.sh

    check_vcf_validity ${input_vcf} || { \
    log "The input vcf file ${input_vcf} is not valid. Please check the file"; return 1; }

    #Run CADD and CADD only output gzipped file.
    [[ -z ${output_file} ]] && output_file=${input_vcf/.vcf*/.cadd.tsv.gz}
    [[ ! ${output_file} =~ \.gz$ ]] && output_file=${output_file}.gz  # Ensure the output file is gzipped

    # Check if the output file is valid and updated. If so, skip this function
    [[ -f ${output_file/.gz/} ]] && [[ ${output_file/.gz/} -nt ${input_vcf} ]] && { \
    log "Output file ${output_file} is valid and updated. Skip this function for now."; \
    return 0; }

    # First convert UCSC-style chr names to GRCh-style
    local nochr_vcf=${input_vcf/.vcf.gz/.nochr.vcf.gz}
    liftover_from_ucsc_to_GRCh \
    ${input_vcf} \
    "${BASE_DIR}/data/liftover/ucsc_to_GRC.contig.map.tsv" \
    ${nochr_vcf} || { \
    log "Failed to convert chromosome names from UCSC to GRCh style for ${input_vcf}. Quit now"; 
    return 1; }

    # Determine the genome tag based on the assembly
    if [[ ${assembly} == "GRCh37" ]] || [[ ${assembly} == "hg19" ]]; then
        local genome_tag="GRCh37"
    elif [[ ${assembly} == "GRCh38" ]] || [[ ${assembly} == "hg38" ]]; then
        local genome_tag="GRCh38"
    else
        log "The assembly ${assembly} is not supported. Please check the config file ${config_file}"
        return 1;
    fi

    # Run CADD
    log "Running CADD script ${cadd_script} with the following command"
    log "conda run -n priva_acmg ${cadd_script} -c ${threads} -a -p -m -d -g ${genome_tag} -o ${output_file} ${nochr_vcf}"
    export TMPDIR=${tmp_dir}
    conda run -n priva_acmg ${cadd_script} \
    -c ${threads} \
    -a -p -m -d \
    -g ${genome_tag} \
    -o ${output_file} \
    ${nochr_vcf} && \
    gunzip -f ${output_file} && \
    mawk -F '\t' '$1 !~ /^##/{print;}' ${output_file/.gz/} > ${output_file/.gz/}.tmp && \
    mv ${output_file/.gz/}.tmp ${output_file/.gz/} && \
    display_table ${output_file/.gz/} && \
    update_yaml ${config_file} "cadd_output_file" ${output_file/.gz/} || { \
    log "Failed to run CADD on ${input_vcf}. Quit now"; return 1; }
}


function use_hub_vcf_annotations() {
    local input_vcf
    local hub_vcf
    local output_covered_vcf
    local output_uncovered_vcf
    local threads
    local total_count
    local tmp_dir

    # Parse arguments using getopts
    local OPTIND=1
    while getopts "i:h:o:u:t:c::d::" opt; do
        case ${opt} in
            i) input_vcf=${OPTARG} ;;
            h) hub_vcf=${OPTARG} ;;
            o) output_covered_vcf=${OPTARG} ;;
            u) output_uncovered_vcf=${OPTARG} ;;
            t) threads=${OPTARG} ;;
            c) total_count=${OPTARG} ;;
            d) tmp_dir=${OPTARG} ;;
            *) log "Invalid option: -${OPTARG}" && return 1 ;;
        esac
    done

    log "Attempting to use hub VCF: ${hub_vcf}"

    # Check if hub VCF exists and is valid (basic check)
    if [[ ! -f "${hub_vcf}" ]] || ! check_vcf_validity "${hub_vcf}"; then
        log "Hub VCF '${hub_vcf}' is missing or invalid. All variants marked as uncovered."
        return 1
    fi

    if [[ -z ${total_count} ]]; then
        total_count=$(count_vcf_records "${input_vcf}")
    fi

    if [[ -z ${tmp_dir} ]]; then
        tmp_dir=${TMPDIR}
    fi

    local tmp_isec_dir=$(mktemp -d --tmpdir="$tmp_dir" isec_hub.XXXXXX)
    log "Created temporary intersection directory: ${tmp_isec_dir}"

    # Use isec -p to efficiently find intersection and complements
    log "Running bcftools isec to find covered/uncovered sites using command:"
    >&2 echo "bcftools isec --threads ${threads} -Oz -p ${tmp_isec_dir} ${input_vcf} ${hub_vcf}"
    bcftools isec --threads ${threads} -Oz -p "${tmp_isec_dir}" "${input_vcf}" "${hub_vcf}" || \
    { log "ERROR: bcftools isec failed." && \
      rm -rf "${tmp_isec_dir}" && \
      return 1; }

    local uncovered_sites_vcf="${tmp_isec_dir}/0000.vcf.gz"
    local covered_sites_vcf="${tmp_isec_dir}/0002.vcf.gz"

    # --- Handle Uncovered Variants ---
    if check_vcf_validity "${uncovered_sites_vcf}" 1; then
        log "Extracting uncovered variants from ${input_vcf} using command:"
        bcftools query -f '%CHROM\t%POS\n' ${uncovered_sites_vcf} > ${uncovered_sites_vcf/.vcf*/.pos.tsv} && \
        local uncovered_sites_file=${uncovered_sites_vcf/.vcf*/.pos.tsv}
        >&2 echo "parallel_extract_variants ${uncovered_sites_file} ${input_vcf} ${threads} ${output_uncovered_vcf}"
		parallel_extract_variants \
		${uncovered_sites_file} \
		${input_vcf} \
		${threads} \
		${output_uncovered_vcf} && \
		display_vcf ${output_uncovered_vcf} || \
		{ log "Failed to extract uncovered variants from ${input_vcf} using ${uncovered_sites_file}. Return with error." && return 1; }
        tabix -p vcf "${output_uncovered_vcf}" && \
        display_vcf "${output_uncovered_vcf}" && \
        local uncovered_count=$(count_vcf_records "${output_uncovered_vcf}")
        log "Found ${uncovered_count} uncovered variants, saved to ${output_uncovered_vcf}"
    else
        log "No uncovered variants found in ${tmp_isec_dir}"
        >&2 echo "$(ls -lh ${tmp_isec_dir})"
        local uncovered_count=0
    fi

    # --- Handle Covered Variants ---
	local covered_count=$(count_vcf_records "${covered_sites_vcf}")
    if [[ ${covered_count} -gt 0 ]]; then
        log "Annotating extracted covered variants with INFO from ${hub_vcf}..."
        # Annotate using only INFO fields from the hub
        bcftools annotate --threads ${threads} -a "${hub_vcf}" -c INFO -Ou "${input_vcf}" | \
		bcftools filter --threads ${threads} -i 'INFO/CSQ != "."' -Oz -o "${output_covered_vcf}" && \
        tabix -p vcf "${output_covered_vcf}"
        log "Found ${covered_count} covered variants, annotated and saved to ${output_covered_vcf}"
    else
        local uncovered_count=${total_count}
        log "No covered variants found. Directly return error to proceed annotation across all variants in the input VCF file ${input_vcf}, which contains ${total_count} variants"
    fi

    # Clean up temporary directory
    log "Cleaning up temporary directory: ${tmp_isec_dir}"
    rm -rf "${tmp_isec_dir}"

    # Determine return code: 0 if all variants were covered, 1 if uncovered variants exist
    if [[ ${uncovered_count} -eq 0 ]]; then
        log "All variants were found in the hub VCF. Skipping annotation pipeline for uncovered variants."
        return 0
    else
        log "Proceeding to annotate ${uncovered_count} uncovered variants."
        return 1
    fi
}


function update_hub_vcf() {
    local newly_annotated_vcf=$1
    local hub_vcf=$2
    local threads=${3}

    if [[ -z ${threads} ]]; then
        threads=1
    fi
    
    # Skip if no new annotations
    local new_variant_count=$(count_vcf_records "${newly_annotated_vcf}")
    if [[ ! -f "${newly_annotated_vcf}" ]] || [[ ${new_variant_count} -eq 0 ]]; then
        log "No new variants to add to hub VCF"
        return 0
    fi

    # Create a temp newly_annotated_vcf with only INFO fields (remove all FORMAT fields information)
    local tmp_tag=$(randomID)
    local temp_newly_annotated_vcf=${newly_annotated_vcf/.vcf/.tmp.${tmp_tag}.vcf}
    remove_format_fields ${newly_annotated_vcf} ${temp_newly_annotated_vcf}
    
    # If hub VCF doesn't exist or is invalid, create it
    if [[ ! -f "${hub_vcf}" ]] || ! check_vcf_validity "${hub_vcf}"; then
        log "Hub VCF was missing or invalid, creating new hub with ${new_variant_count} variants"
        mv "${temp_newly_annotated_vcf}" "${hub_vcf}" && \
        mv "${temp_newly_annotated_vcf}.tbi" "${hub_vcf}.tbi" && \
        return 0 || \
        { log "Failed to create new hub VCF by directly copying the newly annotated VCF ${temp_newly_annotated_vcf}. Return with error." && return 1; }
    fi
    
    # Merge hub VCF with newly annotated VCF
    local temp_hub="${hub_vcf/.vcf*/.${tmp_tag}.vcf.gz}"
    remove_format_fields ${hub_vcf}
    
    # Concatenate files and remove duplicates (keeping the newer version)
    bcftools concat -a "${hub_vcf}" "${temp_newly_annotated_vcf}" -Ou | \
    bcftools sort -Oz -o "${temp_hub}" && tabix -p vcf -f "${temp_hub}" && \
    normalize_vcf ${temp_hub} ${hub_vcf} ${ref_genome} ${threads} && \
    announce_remove_tmps "${temp_hub}" "${temp_newly_annotated_vcf}" && \
    log "Updated hub VCF with $(count_vcf_records "${newly_annotated_vcf}") new variants" || \
    { log "Failed to update hub VCF with $(count_vcf_records "${newly_annotated_vcf}") new variants. Return with error." && return 1; }
}

function merge_annotated_vcfs() {
    local covered_vcf=$1
    local newly_annotated_vcf=$2
    local output_vcf=$3
    
    # Handle edge cases
    if [[ ! -s "${newly_annotated_vcf}" ]] || [[ $(count_vcf_records "${newly_annotated_vcf}") -eq 0 ]]; then
		log "No newly annotated variants in ${newly_annotated_vcf}"
        if [[ -s "${covered_vcf}" ]] && [[ $(count_vcf_records "${covered_vcf}") -gt 0 ]]; then
            mv "${covered_vcf}" "${output_vcf}"
            tabix -p vcf "${output_vcf}"
            announce_remove_tmps "${covered_vcf}.*"
            log "Using only covered variants, no newly annotated variants"
        else
            log "ERROR: No covered variants in ${covered_vcf} or newly annotated variants in ${newly_annotated_vcf} available"
            return 1
        fi
        return 0
    fi
    
    if [[ ! -f "${covered_vcf}" ]] || [[ $(count_vcf_records "${covered_vcf}") -eq 0 ]]; then
        [[ ${newly_annotated_vcf} != ${output_vcf} ]] && \
		mv "${newly_annotated_vcf}" "${output_vcf}" && \
        tabix -p vcf "${output_vcf}" || :
        log "Using only newly annotated variants, no covered variants"
        return 0
    fi
    
    # Merge the files
    bcftools concat -a "${covered_vcf}" "${newly_annotated_vcf}" -Ou | \
    bcftools sort -Oz -o "${output_vcf/.vcf*/.tmp.vcf.gz}" && \
	mv "${output_vcf/.vcf*/.tmp.vcf.gz}" "${output_vcf}" && \
    tabix -p vcf "${output_vcf}" && \
    announce_remove_tmps "${covered_vcf}*" && \
    log "Merged $(count_vcf_records "${covered_vcf}") covered variants with $(count_vcf_records "${newly_annotated_vcf}") newly annotated variants"
}



function find_cached_cadd_variants() {
    local input_vcf=$1
    local hub_cadd_file=$2
    local output_covered_tsv=$3
    local output_uncovered_vcf=$4
    local threads=${5}

	if [[ -z ${threads} ]]; then
		threads=1
	fi

    # Create a temporary "positions file" from input VCF for lookup
    local pos_file=$(mktemp --tmpdir="$TMPDIR" pos_file.XXXXXX)
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${input_vcf}" > "${pos_file}" && \
	display_table ${pos_file}
    
    # Create a temporary file for variants not found in hub
    local uncovered_pos=$(mktemp --tmpdir="$TMPDIR" uncovered_pos.XXXXXX.tsv)
    touch "${uncovered_pos}"
    
    # Extract matching variants from CADD TSV using Python script
    python ${SCRIPT_DIR}/cadd_hub_util.py match \
        --hub "${hub_cadd_file}" \
        --pos "${pos_file}" \
        --covered "${output_covered_tsv}" \
        --uncovered "${uncovered_pos}" && \
	ls -lhtr ${output_covered_tsv} ${uncovered_pos} || \
	{ log "Failed to match variants in CADD hub TSV" && return 1; }
	if [[ $(cat ${uncovered_pos} | wc -l) -eq 0 ]]; then
		local match_result=0
	else
		local match_result=1
	fi
    
    # Clean up temporary positions file
    rm -f "${pos_file}"
    
    # Check if we found any uncovered variants
    if [[ ${match_result} -eq 0 ]]; then
        log "All variants found in CADD hub TSV"
        rm -f "${uncovered_pos}"
        return 0
    fi

	local -a uncovered_chrs=($(cut -f 1 ${uncovered_pos} | sort -u))
    
    # Extract uncovered variants from VCF using exact position matching
	parallel_extract_variants \
	${uncovered_pos} \
	${input_vcf} \
	${threads} \
	${output_uncovered_vcf} && \
	display_vcf ${output_uncovered_vcf} || \
	{ log "Failed to extract uncovered variants from ${input_vcf} using ${uncovered_pos}. Return with error." && return 1; }
    
    # Clean up temporary positions file
    rm -f "${uncovered_pos}"
    
    # Check if any uncovered variants
    if [[ $(count_vcf_records "${output_uncovered_vcf}") -eq 0 ]]; then
        log "No uncovered variants"
        return 0
    else
        log "Will process $(count_vcf_records "${output_uncovered_vcf}") uncached variants"
        return 1
    fi
}


function update_hub_cadd() {
    local new_cadd_tsv=$1
    local hub_cadd_file=$2
    
    # Skip if no new annotations
    if [[ ! -f "${new_cadd_tsv}" ]] || [[ $(wc -l < "${new_cadd_tsv}") -le 1 ]]; then
        log "No new CADD scores to add to hub"
        return 0
    fi
    
    # Ensure hub directory exists
    local hub_dir=$(dirname "${hub_cadd_file}")
    [[ -d "${hub_dir}" ]] || {
        log "Error: Hub directory ${hub_dir} does not exist to store the hub cadd anno file ${hub_cadd_file}"
        return 1
    }
    
    # Update hub CADD TSV with new scores
    python ${SCRIPT_DIR}/cadd_hub_util.py update --new "${new_cadd_tsv}" --hub "${hub_cadd_file}"
    log "Updated CADD hub TSV ${hub_cadd_file}"
}


function merge_cadd_results() {
    local covered_tsv=$1
    local new_cadd_tsv=$2
    local output_tsv=$3
    
    # Merge covered and new CADD scores
    python ${SCRIPT_DIR}/cadd_hub_util.py merge \
        --covered "${covered_tsv}" \
        --new "${new_cadd_tsv}" \
        --output "${output_tsv}"
    
    if [[ $? -eq 0 ]]; then
        log "Successfully merged CADD results"
        return 0
    else
        log "Failed to merge CADD results"
        return 1
    fi
}

if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
    declare -a func_names=($(typeset -f | awk '!/^main[ (]/ && /^[^ {}]+ *\(\)/ { gsub(/[()]/, "", $1); printf $1" ";}'))
    declare -a input_func_names=($(return_array_intersection "${func_names[*]}" "$*"))
    declare -a arg_indices=($(get_array_index "${input_func_names[*]}" "$*"))
    if [[ ${#input_func_names[@]} -gt 0 ]]; then
        log "Seems like the command is trying to directly run a function in this script ${BASH_SOURCE[0]}."
        first_func_ind=${arg_indices}
        log "The identified first func name is at the ${first_func_ind}th input argument, while the total input arguments are: $*"
        following_arg_ind=$((first_func_ind + 1))
        log "Executing: ${*:${following_arg_ind}}"
        "${@:${following_arg_ind}}"
    else
        log "Directly run main_workflow with input args: $*"
        main_workflow "$@"
    fi
fi

