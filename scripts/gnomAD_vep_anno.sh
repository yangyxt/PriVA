#! /usr/bin/env bash

# This script is used to annotate the variants with gnomAD and VEP
# Note that this script is not part of the pipeline, it is intended to be used for once to gather statistics about homozygous LoF variants in gnomAD to evaluate the functional essentiality of the gene(transcripts)

SELF_SCRIPT="$(realpath ${BASH_SOURCE[0]})"
SCRIPT_DIR="$(dirname "${SELF_SCRIPT}")"
BASE_DIR="$(dirname ${SCRIPT_DIR})"
if [[ ${BASE_DIR} == "/" ]]; then BASE_DIR=""; fi
DATA_DIR="${BASE_DIR}/data"

if [[ -z $TMPDIR ]]; then TMPDIR=/tmp; fi

# Source the other script
source "${SCRIPT_DIR}/common_bash_utils.sh"


function annotate_vep_gnomAD_per_chromosome() {
	local vcf_file
	local config_file
	local output_vcf
	local output_dir
	local threads

	# Use getopts to parse the command line arguments
	while getopts "v:c:o::d::t::" opt; do
		case ${opt} in
			v) input_vcf=${OPTARG} ;;
			c) config_file=${OPTARG} ;;
			o) output_vcf=${OPTARG} ;;
			d) output_dir=${OPTARG} ;;
			t) threads=${OPTARG} ;;
		esac
	done


	if [[ -z ${output_vcf} ]] && [[ -z ${output_dir} ]]; then
		output_vcf="${input_vcf/.vcf*/.vep.vcf.gz}" && \
		log "No output VCF file is provided, use the default name: ${output_vcf}"
	elif [[ -z ${output_vcf} ]] && [[ -n ${output_dir} ]]; then
		output_vcf="${output_dir}/$(basename ${input_vcf/.vcf*/.vep.vcf.gz})" && \
		log "No output VCF file is provided, use the default name: ${output_vcf}"
	else
		log "Output VCF file is provided: ${output_vcf}"
	fi

    
	# Set variables with command line arguments taking precedence over config file
    local assembly=$(read_yaml ${config_file} "assembly")
    local ref_genome=$(read_yaml ${config_file} "ref_genome")
    local vep_cache_dir=$(read_yaml ${config_file} "vep_cache_dir")
    local vep_plugins_dir=$(read_yaml ${config_file} "vep_plugins_dir")
    local vep_plugins_cachedir=$(read_yaml ${config_file} "vep_plugins_cachedir")
    if [[ -z ${threads} ]]; then
        threads=$(read_yaml ${config_file} "threads")
    fi

    # Plugin cache files
    local utr_annotator_file=$(read_yaml ${config_file} "utr_annotator_file")
    local loeuf_prescore=$(read_yaml ${config_file} "loeuf_prescore")
    local alphamissense_prescore=$(read_yaml ${config_file} "alphamissense_prescore")
    local spliceai_snv_prescore=$(read_yaml ${config_file} "spliceai_snv_prescore")
    local spliceai_indel_prescore=$(read_yaml ${config_file} "spliceai_indel_prescore")
    local primateai_prescore=$(read_yaml ${config_file} "primateai_prescore")
    local conservation_file=$(read_yaml ${config_file} "conservation_file")
	local splicevault_prescore=$(read_yaml ${config_file} "splicevault_prescore")
	local loftee_repo=$(read_yaml ${config_file} "loftee_repo")
	local human_ancestor_fasta=$(read_yaml ${config_file} "human_ancestor_fasta")
	local loftee_conservation_file=$(read_yaml ${config_file} "loftee_conservation_file")
	local gerp_bigwig=$(read_yaml ${config_file} "gerp_bigwig")


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

	[[ ${assembly} == "GRCh37" ]] && local gerp_bw_argument=""
	if [[ ${assembly} == "GRCh38" ]]; then
		local gerp_bw_argument=",gerp_bigwig ${gerp_bigwig}" && \
		check_path "$gerp_bigwig" "file" "gerp_bigwig" || has_error=1
	fi

    # Exit if any errors were found
    if [[ $has_error -gt 0 ]]; then
        return 1
    fi

	check_vcf_infotags ${output_vcf} "CSQ" && \
	log "The output vcf ${output_vcf} is already annotated by VEP. Skip this step" && \
	return 0

	local tmp_tag=$(randomID)
	local tmp_output=${output_vcf/.vcf.gz/.${tmp_tag}.vcf.gz}

	log "Start to annotate the vcf ${input_vcf} with VEP and the command is: vep -i ${input_vcf} --format vcf --verbose --vcf --species homo_sapiens --use_transcript_ref --assembly ${assembly} --cache --offline --merged --domains --hgvs --numbers --symbol --canonical --total_length --variant_class --gene_phenotype --stats_file ${input_vcf/.vcf*/.vep.stats.html} --fork ${threads} --buffer_size 10000 --fasta ${ref_genome} --dir_cache ${vep_cache_dir} --dir_plugins ${vep_plugins_dir} -plugin UTRAnnotator,file=${utr_annotator_file} -plugin LOEUF,file=${loeuf_prescore},match_by=transcript -plugin AlphaMissense,file=${alphamissense_prescore} -plugin LoF,loftee_path:${loftee_repo},human_ancestor_fa:${human_ancestor_fasta},conservation_file:${loftee_conservation_file}${gerp_bw_argument} -plugin SpliceAI,snv=${spliceai_snv_prescore},indel=${spliceai_indel_prescore},cutoff=0.5 -plugin PrimateAI,${primateai_prescore} -plugin SpliceVault,file=${splicevault_prescore} -plugin Conservation,${conservation_file},MAX -plugin NMD --force_overwrite"

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
    --dir_plugins ${vep_plugins_dir} \
    -plugin UTRAnnotator,file=${utr_annotator_file} \
    -plugin LOEUF,file=${loeuf_prescore},match_by=transcript \
    -plugin AlphaMissense,file=${alphamissense_prescore} \
	-plugin LoF,loftee_path:${loftee_repo},human_ancestor_fa:${human_ancestor_fasta},conservation_file:${loftee_conservation_file}${gerp_bw_argument} \
    -plugin SpliceAI,snv=${spliceai_snv_prescore},indel=${spliceai_indel_prescore},cutoff=0.5 \
    -plugin PrimateAI,${primateai_prescore} \
	-plugin SpliceVault,file=${splicevault_prescore} \
    -plugin Conservation,${conservation_file},MAX \
	-plugin NMD \
    --force_overwrite \
	-o ${tmp_output} && \
	bcftools sort -Oz -o ${output_vcf} ${tmp_output} && \
	tabix -f -p vcf ${output_vcf} && \
	announce_remove_tmps ${tmp_output} && \
	display_vcf ${output_vcf}
}



function parallel_annotate_vep_gnomAD() {
	local config_file=${1}
	local output_dir=${2}

	local threads=$(read_yaml ${config_file} "threads")
	local tmp_dir=$(read_yaml ${config_file} "tmp_dir")

	local gnomAD_chrX_vcf=$(read_yaml ${config_file} "gnomad_vcf_chrX")
	[[ ! -f ${gnomAD_chrX_vcf} ]] && log "The gnomAD VCF file ${gnomAD_chrX_vcf} does not exist" && return 1 || :

	local -a chromosomes=( "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y" )
	log "The chromosomes are: ${chromosomes[*]}"
	if [[ ${#chromosomes[@]} -lt 24 ]]; then
		log "The number of chromosomes is less than 24, which means there are something wrong with the definitions of the BASH array"
		return 1
	fi

	local -a gnomAD_vcf_files
	for chrom in "${chromosomes[@]}"; do
		gnomAD_vcf_files+=( "${gnomAD_chrX_vcf/.chrX./.chr${chrom}.}" )
	done

	log "All the gnomAD VCF files are: ${gnomAD_vcf_files[*]}"
	# Print the biggest integer we got with awk (using total threads divided by the number of chromosomes)
	local job_thread_num=$(awk -v t="${threads}" 'BEGIN {printf "%i", t/8; }' )
	log "The number of jobs is: 8 and the number of threads per job is: ${job_thread_num}"

	parallel -j 8 --dry-run --link --halt soon,fail=1 \
	--joblog ${output_dir}/vep_anno_gnomAD.log \
	--tmpdir ${tmp_dir} \
	bash ${SCRIPT_DIR}/gnomAD_vep_anno.sh annotate_vep_gnomAD_per_chromosome -v {1} -c ${config_file} -d ${output_dir} -t ${job_thread_num} '>' ${output_dir}/vep_anno_gnomAD_{2}.log '2>&1' ::: ${gnomAD_vcf_files[@]} ::: ${chromosomes[@]} && \
	export TMPDIR=${tmp_dir}
	parallel -j 8 \
	--halt soon,fail=1 \
	--link \
	--joblog ${output_dir}/vep_anno_gnomAD.log \
	--tmpdir ${tmp_dir} \
	bash ${SCRIPT_DIR}/gnomAD_vep_anno.sh annotate_vep_gnomAD_per_chromosome -v {1} -c ${config_file} -d ${output_dir} -t ${job_thread_num} '>' ${output_dir}/vep_anno_gnomAD_{2}.log '2>&1' ::: ${gnomAD_vcf_files[@]} ::: ${chromosomes[@]}
	check_parallel_joblog ${output_dir}/vep_anno_gnomAD.log
}



if [[ ${BASH_SOURCE[0]} == ${0} ]]; then
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
		parallel_annotate_vep_gnomAD "$@"
	fi
fi

