#!/usr/bin/env bash
# This script is used to define all the annotation process functions and run the main pipeline
# For now, the script is only used to process short variants (Indels and SNVs). Large structural variants and CNVs are not included.
# Since there will be unexpected indels in the input VCF file. So prescores are not enough to cover them. We need to use some tools to address this issue.
# In this workflow, we use VEP, CADD, SpliceAI, VEP-plugin(UTRannotator). We abandon the ANNOVAR for its related releasing and distribution restrictions.
# The script expect a Family VCF input along with a pedigree table contains the Family pedigree information.
# Notice that the family sample ID should be consistent between the pedigree table and the VCF file.
# Notice that the proband of the family should stay as the first record of the family in the pedigree table

# For the same set of configurations (arguments), the pipeline should start with the position it ends last time, unless user specifically asked to forcefully rerun the pipeline

# Maintainer: yangyxt@gmail.com


FORCE_RERUN="False"
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


function log() {
    local msg="$1"
    local script_name="${BASH_SOURCE[1]##*/}"
    local func_name="${FUNCNAME[1]}"
    local line_num="${BASH_LINENO[0]}"
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    >&2 echo "[$timestamp] [Script $script_name: Line $line_num] [Func: $func_name] $msg"
}


function main_workflow() {
	local input_vcf
    local ped_file
	local fam_name
	local assembly
    local annovar_dir
	local output_dir
	local ref_genome
	local gnomAD_db_dir
	local threads
	local af_cutoff

    local TEMP
	log "Raw input arguments: $#"
    TEMP=$(getopt -o a:i:p:f:d:o:r:g:t: --long assembly:,input_vcf:,ped_file:,fam_name:,annovar_dir:,output_dir:,ref_genome:,gnomAD_db_dir:,threads:,af_cutoff: -- "$@")

	log "TEMP: $TEMP"
    # if getopt failed, return an error
    [[ $? != 0 ]] && return 1

    eval set -- "$TEMP"

    while true; do
        case "$1" in
            -a|--assembly)
                assembly="$2"
                shift 2
                ;;
            -p|--ped_file)
                ped_file="$2"
                shift 2
                ;;
			-f|--fam_name)
                fam_name="$2"
                shift 2
                ;;
			-i|--input_vcf)
                input_vcf="$2"
                shift 2
                ;;
			-d|--annovar_dir)
                annovar_dir="$2"
                shift 2
                ;;
			-o|--output_dir)
                output_dir="$2"
                shift 2
                ;;
			-r|--ref_genome)
                ref_genome="$2"
                shift 2
                ;;
			-g|--gnomAD_db_dir)
				gnomAD_db_dir="$2"
				shift 2
				;;
			-t|--threads)
				threads="$2"
				shift 2
				;;
			--af_cutoff)
				af_cutoff="$2"
				shift 2
				;;
            --)
                shift
                break
                ;;
            *)
                log "Invalid option"
                return 2
                ;;
        esac
    done

	if [[ -z ${af_cutoff} ]]; then
		local af_cutoff=0.05
	fi

	# Prepare proband_ID and family_name
	if [[ -z ${fam_name} ]]; then
		local proband_ID=$(bcftools query -l ${input_vcf} | head -1)
		local fam_name=$(awk -F '\t' '$2 == "'${proband_ID}'"{printf "%s", $1; exit 0;}' ${ped_file})
    	local -a patient_IDs=($(awk -F '\t' '$1 == "'${fam_name}'" && $6 == "2" {printf "%s ", $2}' ${ped_file}))
	else
		local -a patient_IDs=($(awk -F '\t' '$1 == "'${fam_name}'" && $6 == "2" {printf "%s ", $2}' ${ped_file}))
		local proband_ID=${patient_IDs[0]}
	fi
    
	if [[ -z ${output_dir} ]]; then local output_dir=$(dirname ${input_vcf}); fi
	if [[ -z ${assembly} ]]; then local assembly=$(check_vcf_assembly_version ${input_vcf}); fi
	if [[ -z ${assembly} ]]; then
		# Here expecting the specified reference genome to have a name like ucsc.hg19.fasta or ucsc.hg38.fasta
		local ref_gen_extraction=$(basename ${ref_genome} | awk -F '.' '{printf "%s", $2;}')
		if [[ ${ref_gen_extraction} =~ ^hg[13][98]$ ]]; then
			local assembly=${ref_gen_extraction}
		else
			log "Please specify the assembly as hg19 or hg38. Failed to extract the assembly version either from input VCF or input ref_genome fasta file. Quit now."
			return 1;
		fi
	fi

	# Preprocess the input vcf to:  
	# Remove the variants not located in primary chromsomes
	# Convert the contig names to UCSC style. Meaning mapping VCF to hg19 or hg38
	# Sort and normalize (including indel left alignment)
	# Remove variants where proband DP < 5
	# Remove VCF records where all patients in the family has either missing or homo ref genotype. 

	local pre_anno_vcf=${input_vcf/.vcf*/.preanno.vcf.gz}
	preprocess_vcf \
	${input_vcf} \
	${ped_file} \
	${fam_name} \
	${ref_genome} \
	${pre_anno_vcf} && \
	display_vcf ${pre_anno_vcf} || { \
	log "Failed to preprocess the input vcf ${input_vcf} for annotation. Quit with error."; \
	return 1; }


	log "The assembly using is ${assembly}."
	log "The pedigree file used is ${ped_file}."
	log "The Family Name of the inspecting family is ${fam_name}."
	log "The VCF preprocessed for annotation is ${pre_anno_vcf}."
	log "The output folder for annotation results is ${output_dir}"



	# First annotate gnomAD aggregated frequency and number of homozygous ALT allele carriers
	# Filter on the allele frequency but first assign inhouse gnomAD AN AC AF data to the variants
	# Since the archive of the aggregated gnomAD database is too big. The tar file will be uploaded to our google drive and open for download.
	local gnomAD_tab=${reformat_tab/.tsv/.gnomAD.tsv}
	add_agg_gnomAD_data \
	${reformat_tab} \
	${threads} \
	${assembly} \
	${gnomAD_db_dir} \
	${gnomAD_tab} || { \
	log "Failed to add aggregated gnomAD annotation on ${reformat_tab}. Quit now"
	return 1; }


	local lowfreq_tab=${gnomAD_tab/.tsv/.lowfreq.tsv}
	filter_on_AF_HOMOALT \
	-i ${gnomAD_tab} \
	-o ${lowfreq_tab} \
	--af_cutoff ${af_cutoff} || { \
	log "Failed to filter on gnomAD annotation on ${gnomAD_tab}. Quit now"
	return 1; }


	# Now after filtering on the AF and HOMOALT no of variants. We need to filter on the pedigree information. 
	# Specifically, we only filter out variants carried by healthy parents in homozygous form
	local pedigree_filtered_tab=${lowfreq_tab/.tsv/.rmfit.tsv}
	filter_allele_based_on_pedigree_with_py \
    -i ${lowfreq_tab} \
    -o ${pedigree_filtered_tab} \
    -p ${ped_file} \
    -f ${fam_name} || { \
	log "Failed to filter on pedigree information on ${lowfreq_tab}. Quit now"
	return 1; }


	# Now we've done the variants filtration
	local filtered_vcf=${annovar_vcf/.vcf/.filtered.vcf}
	prepare_vcf_add_varID \
	${pedigree_filtered_tab} \
	${assembly} \
	${filtered_vcf}


	# Next we need to perform VEP and SpliceAI and 
}


function preprocess_vcf() {
	local input_vcf=${1}
	local output_vcf=${5}
	local ped_file=${2}
	local fam_name=${3}
	local ref_genome=${4}

	local proband=($(awk -F '\t' '$1 == "'${fam_name}'" && $6 == "2" {printf "%s\n", $2}' ${ped_file} | head -1))
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
		check_vcf_multiallelics ${output_vcf} && \
		[[ $(bcftools view -Ov -H ${output_vcf} | cut -f 5 | grep -c '\.') -eq 0 ]] && \
		whether_same_varset ${input_vcf} ${output_vcf} && \
		contain_only_variants ${output_vcf} ${proband} && \
		[[ ! ${FORCE_RERUN} =~ [Tt][Rr][Uu][Ee]$ ]]; then
		log "The ${output_vcf} is valid and udpated. Skip this function"
		return 0;
    fi

	# First filter out variants that not in primary chromosomes (chrM to chrY)(MT to Y)
	bcftools sort -Oz -o ${input_vcf/.vcf*/.sorted.vcf.gz} ${input_vcf} && \
	tabix -f -p vcf ${input_vcf/.vcf*/.sorted.vcf.gz} && \
	bcftools view -r $(cat ${BASE_DIR}/data/liftover/ucsc_GRC.primary.contigs.tsv | tr '\n' ',') -Ou ${input_vcf/.vcf*/.sorted.vcf.gz} | \
	bcftools sort -Oz -o ${input_vcf/.vcf*/.primary.vcf.gz} && \
	tabix -f -p vcf ${input_vcf/.vcf*/.primary.vcf.gz} || \
	{ log "Failed to generate a VCF file that only contains records in primary chromosomes"; \
	  return 1; }
 
	# First check whether the VCF is using NCBI or UCSC assembly
	local assembly_version=$(check_vcf_contig_version ${input_vcf/.vcf*/.primary.vcf.gz})
	if [[ ${assembly_version} == "ncbi" ]]; then
		log "The input vcf is detected to map variants to GRC assemblies instead of UCSC assemblies" && \
		liftover_from_GRCh_to_hg \
		${input_vcf/.vcf*/.primary.vcf.gz} \
		${BASE_DIR}/data/liftover/ucsc_to_GRC.contig.map.tsv \
		${input_vcf/.vcf*/.ucsc.vcf.gz}
	else
		if [[ $(vcf_content_sha256 ${input_vcf/.vcf*/.primary.vcf.gz} ) != $(vcf_content_sha256 ${input_vcf/.vcf*/.ucsc.vcf.gz}) ]]; then
			cp -f ${input_vcf/.vcf*/.primary.vcf.gz} ${input_vcf/.vcf*/.ucsc.vcf.gz}
		fi
	fi
	
	# First sort the input_vcf,
	# Then normalize the input_vcf with bcftools 
	normalize_vcf ${input_vcf/.vcf*/.ucsc.vcf.gz} ${input_vcf/.vcf*/.norm.vcf.gz} ${ref_genome} $TMPDIR

	# Then filter the variants where DP value less than or equalt to 5
	bcftools filter -i 'FORMAT/DP[*] >= 5' -Oz -o ${input_vcf/.vcf*/.filtered.vcf.gz} ${input_vcf/.vcf/.norm.vcf} && \
	tabix -f -p vcf ${input_vcf/.vcf*/.filtered.vcf.gz} && \
	announce_remove_tmps ${input_vcf/.vcf/.norm.vcf} && \
	announce_remove_tmps ${input_vcf/.vcf*/}*tmp*vcf*

	# Then filter out the variants where no patients in this family has alternative alleles
	filter_records_with_missingGT_on_pat_vcf \
	${input_vcf/.vcf*/.filtered.vcf.gz} \
	${ped_file} \
	${output_vcf} \
	${fam_name} && \
	tabix -f -p vcf ${output_vcf}
}


function filter_records_with_missingGT_on_pat_vcf {
	local input_vcf=${1}
	local ped_file=${2}
	local output_vcf=${3}
	local fam_name=${4}

	log "Input vcf is ${input_vcf}, the pedigree file storing this family's pedigree info is ${ped_file}"
	local -a patient_IDs=($(awk -F '\t' '$1 == "'${fam_name}'" && $6 == "2" {printf "%s ", $2}' ${ped_file}))

	log "patient IDs are ${patient_IDs[*]}"
	log "Before filtering the records where patients GT are all ref or null, ${input_vcf} has $(wc -l ${input_vcf} | awk '{printf $1}') rows."
	
    filter_vcf_by_GT \
	${input_vcf} \
	"$(echo ${patient_IDs[*]} | awk 'BEGIN{FS=OFS="\t";} {gsub(" ", ","); printf "%s", $0}')" \
	${output_vcf} && \
	display_vcf ${output_vcf} && \
	log "Finish filtering the records where patients GT are all ref or null. ${input_vcf} has $(wc -l ${input_vcf} | awk '{printf $1}') rows." $'\n\n'
}


function remove_redundant_AF_cols () {
	# This function is specifically designed for gnomAD211 and gnomAD312 databases from ANNOVAR
	local annovar_table=${1}
	local output_table=${2}
	local -a del_cols=( "AF_popmax" "AF_male" "AF_female" "AF_raw" "AF_afr" "AF_sas" "AF_amr" "AF_eas" "AF_nfe" "AF_fin" "AF_asj" "AF_oth" "non_topmed_AF_popmax" "non_neuro_AF_popmax" "non_cancer_AF_popmax" "controls_AF_popmax" )
	local -a del_col_ids

	# Delete the second AFs and rename the AF to gnomAD_exome_ALL and gnomAD_genome_ALL separately 
	for del_col in "${del_cols[@]}"; do
		local indices=$(awk -F '\t' 'NR == 1{for (i=1;i<=NF;i++) {if ($i == "'${del_col}'") printf "%s,", i;}}' ${annovar_table})
		del_col_ids=${del_col_ids}${indices}
	done

	log "The indices of the tobe deleted columns are ${del_col_ids}" && \
	display_table ${annovar_table} && \
	cut -f ${del_col_ids::-1} --complement ${annovar_table} | \
	mawk 'BEGIN{FS=OFS="\t";}
		  NR == 1{count = 0; \
				for (i=1;i<=NF;i++) { \
					if ($i == "AF") { \
						count++; \
						if (count == 1) { \
							$i = "gnomAD_exome_ALL"; } \
						else if (count == 2) { \
							$i = "gnomAD_genome_ALL"; } \
					} \
				} \
				print; } \
		  NR > 1 {print;}' > ${output_table} && \
	display_table ${output_table}
}


function reformat_annovar_table () {
	local annovar_table=${1}
	local tmp_table=${annovar_table::-4}.tmp.tsv

	remove_redundant_AF_cols \
	${annovar_table} \
	${tmp_table}

	# Then using an in-house python script to reformat the table to achieve these goals:
	# 1. Use VCF coordinates and REF/ALT alleles
	# 2. Split the FORMAT fields to seperate columns for each sample in this family
	# 3. Update the HGNC gene symbol
	python3 ${SCRIPT_DIR}/reformat_annovar_table.py \
	-at ${tmp_table} \
	-av ${annovar_vcf} \
	-ped ${ped_file} \
	-ot ${annovar_table::-4}.reformat.tsv \
	-fam ${fam_name} && \
	update_HGNC_symbol \
    -i ${annovar_table::-4}.reformat.tsv && \
    display_table ${annovar_table::-4}.reformat.tsv
}


function add_agg_gnomAD_data () {
	local anno_table=${1}
	local threads=${2}
	local assembly=${3}
	local gnomAD_db_dir=${4}
	local output_table=${5}

	if [[ -z ${output_table} ]]; then
		local output_table=${anno_table::-4}.gnomAD.tsv
	fi

	python3 ${SCRIPT_DIR}/annotate_agg_gnomAD.py \
	-v ${anno_table} \
	-a "${assembly}" \
	-t ${threads} \
	-d ${gnomAD_db_dir} \
	-o ${output_table} && \
	display_table ${output_table}
}



function filter_on_AF_HOMOALT () {
	local input_table
	local output_table
	local af_cutoff
	local excl_tags   # Should be the filter tag, if two or more are included, then use comma to separate them
	local af_cols	# If two or more are included, then use comma to separate them
	
	local TEMP
	log "Raw input arguments: $#"
    TEMP=$(getopt -o i:o: --long input_table:,output_table:,af_cutoff:,excl_tags:,af_cols: -- "$@")

	log "TEMP: $TEMP"
    # if getopt failed, return an error
    [[ $? != 0 ]] && return 1

    eval set -- "$TEMP"

    while true; do
        case "$1" in
            -i|--input_table)
				input_table="$2"
                shift 2
                ;;
			-o|--output_table)
				output_table="$2"
				shift 2
				;;
			--af_cutoff)
				af_cutoff="$2"
				shift 2
				;;
			--excl_tags)
				excl_tags="$2"
				shift 2
				;;
			--af_cols)
				af_cols="$2"
				shift 2
				;;
			--)
                shift
                break
                ;;
            *)
                log "Invalid option"
                return 2
                ;;
        esac
    done

	if [[ -z ${output_table} ]]; then
		local output_table=${input_table::-4}.lowfreq.tsv
	fi

	if [[ -z ${af_cutoff} ]]; then
		local af_cutoff="0.05"
	fi

	if [[ -z ${af_cols} ]]; then
		local af_cols="gnomAD_exome_ALL,gnomAD_genome_ALL,AF_gnomAD_Controls"
	fi

	python3 ${SCRIPT_DIR}/filter_on_AF_and_nhomoalt.py \
	-i ${input_table} \
	-o ${output_table} \
	-c ${af_cutoff} \
	-l ${af_cols} && \
	display_table ${output_table}
}


function update_HGNC_symbol {
    local OPTIND i g d
    while getopts i:g::d:: args
    do
        case ${args} in
            i) local input_table=$OPTARG ;;
            g) local gene_symbol_col=$OPTARG ;;
            d) local delimiter=$OPTARG ;;
            *) echo "No argument passed. At least pass an argument specifying the family ID"
        esac
    done

    local update_HGNC_py=${SCRIPT_DIR}/HGNC_symbol_updating.py
    local HGNC_anno_table=${DATA_DIR}/hgnc/non_alt_loci_set.tsv

    if [[ -z ${gene_symbol_col} ]]; then
        local gene_symbol_col="Gene.refGene"
    fi

    if [[ -z ${delimiter} ]]; then
        local delimiter=';'
    fi

    time python3 ${update_HGNC_py} \
    -ap ${HGNC_anno_table} \
    -v ${input_table} \
    -ch "${gene_symbol_col}" \
    -d "${delimiter}"
}




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
		log "Directly run main_workflow with input args: $#"
		main_workflow "$@"
	fi
fi

