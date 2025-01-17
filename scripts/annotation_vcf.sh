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
# conda activate acmg

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
		  vep_plugins_cachedir

    # Source the argparse.bash script
    source ${SCRIPT_DIR}/argparse.bash || { log "Failed to source argparse.bash"; return 1; }

	argparse "$@" < ${SCRIPT_DIR}/anno_main_args || { log "Failed to parse arguments"; return 1; }

    # Process config file
    local config_file="${config:-${BASE_DIR}/config.yaml}"
    local -A config_args

    if [[ -f "$config_file" ]]; then
        log "Using config file: $config_file"
        local -a config_keys=(input_vcf assembly vep_cache_dir output_dir ref_genome threads af_cutoff gnomad_vcf_chrX clinvar_vcf vep_plugins_dir vep_plugins_cachedir)
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
	local vep_plugins_dir="${vep_plugins_dir:-${config_args[vep_plugins_dir]}}"
	local vep_plugins_cachedir="${vep_plugins_cachedir:-${config_args[vep_plugins_cachedir]}}"

    # Set and check required paths
    local has_error=0
    check_path "$input_vcf" "file" "input_vcf" || has_error=1

    # Check paths that must exist
    check_path "$ref_genome" "file" "ref_genome" || has_error=1
    check_path "$gnomad_vcf_chrX" "file" "gnomad_vcf_chrX" || has_error=1
    check_path "$clinvar_vcf" "file" "clinvar_vcf" || has_error=1
    check_path "$vep_cache_dir" "dir" "vep_cache_dir" || has_error=1
    check_path "$vep_plugins_dir" "dir" "vep_plugins_dir" || has_error=1
    check_path "$vep_plugins_cachedir" "dir" "vep_plugins_cachedir" || has_error=1

    # Create output directory if it doesn't exist
	local output_dir
	[[ -z ${output_dir} ]] && output_dir=$(dirname ${input_vcf})
    if [[ ! -d "$output_dir" ]]; then
        mkdir -p "$output_dir" || {
            log "Error: Failed to create output_dir: $output_dir"
            has_error=1
        }
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

	# Preprocess the input vcf to:
	# Remove the variants not located in primary chromsomes
	# Convert the contig names to UCSC style. Meaning mapping VCF to hg19 or hg38
	# Sort and normalize (including indel left alignment)
	# Remove variants where proband DP < 5
	# Remove VCF records where all patients in the family has either missing or homo ref genotype.

	local anno_vcf=${input_vcf/.vcf*/.anno.vcf.gz}
	preprocess_vcf \
	-i ${input_vcf} \
	-o ${anno_vcf} \
	-r ${ref_genome} && \
	log "Successfully preprocess the input vcf ${input_vcf} for annotation. The result is ${anno_vcf}" && \
	display_vcf ${anno_vcf} || { \
	log "Failed to preprocess the input vcf ${input_vcf} for annotation. Quit with error."; \
	return 1; }


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
	${clinvar_vcf} || { \
	log "Failed to add ClinVar annotation on ${anno_vcf}. Quit now"
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


	# Now we annotate the VCF with CADD
	Calculate_CADD \
	${anno_vcf} \
	${config_file} || { \
	log "Failed to add CADD annotation on ${anno_vcf}. Quit now"
	return 1; }

}


function preprocess_vcf() {
	local input_vcf
	local ref_genome
	local output_vcf

	# Use getopts to parse the arguments
	local OPTIND=1
	while getopts "i:o:r:" opt; do
		case ${opt} in
			i) input_vcf=${OPTARG};;
			o) output_vcf=${OPTARG};;
			r) ref_genome=${OPTARG};;
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

	if [[ -z ${ref_genome} ]]; then
		log "User does not specify the ref genome fasta file. Cant proceed now. Quit with Error!"
		return 1;
	elif [[ ! -z ${vcf_assembly} ]] && [[ ! ${ref_genome} =~ ${vcf_assembly} ]]; then
		log "User specified ref_genome fasta file ${ref_genome} seems not match with the VCF used assembly ${vcf_assembly}. Quit with Error"
		return 1;
	fi

	# Test if output_vcf is already valid
	if [[ ${output_vcf} -nt ${input_vcf} ]] && \
		check_vcf_validity ${output_vcf} && \
		check_vcf_multiallelics ${output_vcf}; then
		log "The ${output_vcf} is valid and udpated. Skip this function"
		return 0;
    fi

	# First filter out variants that not in primary chromosomes (chrM to chrY)(MT to Y)
	bcftools view -r "$(cat ${BASE_DIR}/data/liftover/ucsc_GRC.primary.contigs.tsv | tr '\n' ',')" -Ou ${input_vcf}| \
	bcftools sort -Oz -o ${input_vcf/.vcf*/.primary.vcf.gz} && \
	tabix -f -p vcf ${input_vcf/.vcf*/.primary.vcf.gz} || \
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
	normalize_vcf ${input_vcf/.vcf*/.ucsc.vcf.gz} ${input_vcf/.vcf*/.norm.vcf.gz} ${ref_genome} $TMPDIR
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

	check_vcf_infotags ${input_vcf} "AC_joint,AN_joint,AF_joint,nhomalt_joint_XX,nhomalt_joint_XY" && \
	log "The input vcf ${input_vcf} already contains the INFO tags AC_joint,AN_joint,AF_joint,nhomalt_joint_XX,nhomalt_joint_XY. We do not need to add them again" && \
	return 0

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
			bcftools view -r ${chr} -Ou ${input_vcf} | \
			bcftools sort -Oz -o ${input_vcf/.vcf*/.${chr}.vcf.gz} && \
			tabix -f -p vcf ${input_vcf/.vcf*/.${chr}.vcf.gz}
			valid_chroms+=( ${chr} )
		else
			log "WARNING!!! The chromosome ${chr} is not supported by gnomAD now. Currently only support chr1 to chrY. Skip this chromosome"
		fi
	done
	log "The valid chromosomes are: ${valid_chroms[*]}"

	# Step 3: Perform annotation with bcftools annotate to add the INFO fields from gnomAD vcfs to the input splitted vcfs
	export gnomad_vcf_chrX
	local gnomad_vcf_suffix=${gnomad_vcf_chrX/*.chrX.vcf/}
	parallel -j ${threads} --dry-run "bcftools annotate -a ${gnomad_vcf_chrX/.chrX.vcf*gz/}.{}.vcf${gnomad_vcf_suffix} -c CHROM,POS,REF,ALT,.INFO/AC_joint,.INFO/AN_joint,.INFO/AF_joint,.INFO/AF_joint_XX,.INFO/AF_joint_XY,.INFO/nhomalt_joint_XX,.INFO/nhomalt_joint_XY -Oz -o ${input_vcf/.vcf*/}.{}.gnomAD.vcf.gz ${input_vcf/.vcf*/}.{}.vcf.gz" ::: "${valid_chroms[@]}" && \
	parallel -j ${threads} --joblog ${input_vcf/.vcf*/.gnomAD.log} "bcftools annotate -a ${gnomad_vcf_chrX/.chrX.vcf*gz/}.{}.vcf${gnomad_vcf_suffix} -c CHROM,POS,REF,ALT,.INFO/AC_joint,.INFO/AN_joint,.INFO/AF_joint,.INFO/AF_joint_XX,.INFO/AF_joint_XY,.INFO/nhomalt_joint_XX,.INFO/nhomalt_joint_XY -Oz -o ${input_vcf/.vcf*/}.{}.gnomAD.vcf.gz ${input_vcf/.vcf*/}.{}.vcf.gz" ::: "${valid_chroms[@]}"; \
	check_parallel_joblog ${input_vcf/.vcf*/.gnomAD.log} || { \
	log "Failed to add aggregated gnomAD annotation on ${input_vcf}. Quit now"
	return 1; }

	# Merge the annotated vcfs together (including the alt-contig-variants file) to form into the output vcf
	# We can filter on AF here
	bcftools concat -Ou ${input_vcf/.vcf*/.chr*.gnomAD.vcf.gz} | \
	bcftools sort -Oz -o ${output_vcf} && \
	tabix -f -p vcf ${output_vcf} && \
	mv ${output_vcf} ${input_vcf} && \
	mv ${output_vcf}.tbi ${input_vcf}.tbi

}



function anno_clinvar_data () {
	local input_vcf=${1}
	local clinvar_vcf=${2}
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

	bcftools annotate -a ${clinvar_vcf} -c CHROM,POS,REF,ALT,.INFO/CLNDN,.INFO/CLNHGVS,.INFO/CLNREVSTAT,.INFO/CLNSIG,.INFO/GENEINFO,.INFO/CLNCSQ:=INFO/CSQ -Ou ${input_vcf} | \
	bcftools sort -Oz -o ${output_vcf} && \
	tabix -f -p vcf ${output_vcf} && \
	mv ${output_vcf} ${input_vcf} && \
	mv ${output_vcf}.tbi ${input_vcf}.tbi && \
	check_vcf_infotags ${input_vcf} "CLNDN,CLNHGVS,CLNREVSTAT,CLNSIG,GENEINFO,CLNCSQ" && \
	display_vcf ${input_vcf} || \
	{ log "Failed to add ClinVar annotation on ${input_vcf}. Quit now"; return 1; }
}




function prepare_vcf_add_varID {
    local input_vcf=${1}
    local output_vcf=${2}

    if [[ -z ${output_vcf} ]]; then
        local output_vcf=${input_vcf/.vcf/.id.vcf}
    fi

    python3 ${SCRIPT_DIR}/generate_varID_for_vcf.py \
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

	[[ ${assembly} == "GRCh37" ]] && local gerp_bw_argument=""
	if [[ ${assembly} == "GRCh38" ]]; then
		local gerp_bw_argument=",gerp_bigwig:${gerp_bigwig}"
		check_path "$gerp_bigwig" "file" "gerp_bigwig" || has_error=1
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
	log "The input vcf ${input_vcf} is already annotated by VEP. Skip this step" && \
	return 0

	local tmp_tag=$(randomID)
	local output_vcf=${input_vcf/.vcf*/.${tmp_tag}.vcf}

	log "Running VEP annotation with the command below:"
	log "vep -i ${input_vcf} --format vcf --verbose --vcf --species homo_sapiens --use_transcript_ref --assembly ${assembly} --cache --offline --merged --domains --hgvs --numbers --symbol --canonical --total_length --variant_class --gene_phenotype --stats_file ${input_vcf/.vcf*/.vep.stats.html} --fork ${threads} --buffer_size 10000 --fasta ${ref_genome} --dir_cache ${vep_cache_dir} --dir_plugins ${vep_plugins_dir} -plugin UTRAnnotator,file=${utr_annotator_file} -plugin LOEUF,file=${loeuf_prescore},match_by=transcript -plugin AlphaMissense,file=${alphamissense_prescore} -plugin LoF,loftee_path:${loftee_repo},human_ancestor_fa:${human_ancestor_fasta},conservation_file:${loftee_conservation_file}${gerp_bw_argument} -plugin SpliceAI,snv=${spliceai_snv_prescore},indel=${spliceai_indel_prescore},cutoff=0.5 -plugin PrimateAI,${primateai_prescore} -plugin SpliceVault,file=${splicevault_prescore} -plugin Conservation,${conservation_file},MAX -plugin NMD --force_overwrite -o ${output_vcf}"

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
	-plugin LoF,loftee_path:${loftee_repo},human_ancestor_fa:${human_ancestor_fasta},conservation_file:${loftee_conservation_file}${gerp_bw_argument} \
    -plugin SpliceAI,snv=${spliceai_snv_prescore},indel=${spliceai_indel_prescore},cutoff=0.5 \
    -plugin PrimateAI,${primateai_prescore} \
	-plugin SpliceVault,file=${splicevault_prescore} \
    -plugin Conservation,${conservation_file},MAX \
	-plugin NMD \
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
	log "conda run -n acmg ${cadd_script} -c ${threads} -a -p -m -d -g ${genome_tag} -o ${output_file} ${nochr_vcf}"
	export TMPDIR=${tmp_dir}
	conda run -n acmg ${cadd_script} \
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

