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
	local output_dir
	local vep_cache_dir
	local ref_genome
	local threads
	local af_cutoff
	local gnomAD_chr1_vcf
	local clinvar_vcf
    local TEMP
	log "Raw input arguments: $#"
    TEMP=$(getopt -o a:i:p:f:c:o:r:t: --long assembly:,input_vcf:,ped_file:,fam_name:,vep_cache_dir:,output_dir:,ref_genome:,threads:,af_cutoff:,gnomAD_chr1_vcf:,clinvar_vcf: -- "$@")

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
			-d|--vep_cache_dir)
                vep_cache_dir="$2"
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
			-t|--threads)
				threads="$2"
				shift 2
				;;
			--af_cutoff)
				af_cutoff="$2"
				shift 2
				;;
			--gnomAD_chr1_vcf)
				gnomAD_chr1_vcf="$2"
				shift 2
				;;
			--clinvar_vcf)
				clinvar_vcf="$2"
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

	# Specifically, we only filter out variants carried by healthy parents in homozygous form
	local ped_filter_vcf=${pre_anno_vcf/.vcf/.ped.vcf}
	filter_allele_based_on_pedigree_with_py \
    -i ${pre_anno_vcf} \
    -o ${ped_filter_vcf} \
    -p ${ped_file} \
    -f ${fam_name} || { \
	log "Failed to filter on pedigree information on ${pre_anno_vcf}. Quit now"
	return 1; }


	log "The assembly using is ${assembly}."
	log "The pedigree file used is ${ped_file}."
	log "The Family Name of the inspecting family is ${fam_name}."
	log "The VCF preprocessed for annotation is ${ped_filter_vcf}."
	log "The output folder for annotation results is ${output_dir}"

	# First annotate gnomAD aggregated frequency and number of homozygous ALT allele carriers
	# Filter on the allele frequency but first assign inhouse gnomAD AN AC AF data to the variants
	# Since the archive of the aggregated gnomAD database is too big. The tar file will be uploaded to our google drive and open for download.
	local gnomAD_anno_vcf=${ped_filter_vcf/.vcf/.gnomAD.vcf}
	anno_agg_gnomAD_data \
	${ped_filter_vcf} \
	${threads} \
	${assembly} \
	${gnomAD_chr1_vcf} \
	${gnomAD_anno_vcf} && \
	log "Successfully add aggregated gnomAD annotation on ${ped_filter_vcf} and filtered on AF >= 0.05. The result is ${gnomAD_anno_vcf}" || { \
	log "Failed to add aggregated gnomAD annotation on ${ped_filter_vcf}. Quit now"
	return 1; }


	# Now we annotate ClinVar variants
	local clinvar_anno_vcf=${gnomAD_anno_vcf/.vcf/.clinvar.vcf}
	anno_clinvar_data \
	${gnomAD_anno_vcf} \
	${clinvar_vcf} \
	${clinvar_anno_vcf} || { \
	log "Failed to add ClinVar annotation on ${gnomAD_anno_vcf}. Quit now"
	return 1; }
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
	if [[ ${assembly_version} =~ "GRCh" ]]; then
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

	# Then filter the variants where DP value less than or equal to 5
	bcftools filter -i 'FORMAT/DP[*] >= 5' -Oz -o ${input_vcf/.vcf*/.filtered.vcf.gz} ${input_vcf/.vcf/.norm.vcf} && \
	tabix -f -p vcf ${input_vcf/.vcf*/.filtered.vcf.gz} && \
	announce_remove_tmps ${input_vcf/.vcf/.norm.vcf} && \
	announce_remove_tmps ${input_vcf/.vcf*/}*tmp*vcf*

	# Then we add uniq IDs to these variants
	prepare_vcf_add_varID \
	${input_vcf/.vcf/.norm.vcf} \
	${output_vcf} && \
	tabix -f -p vcf ${output_vcf}
}




function anno_agg_gnomAD_data () {
	local input_vcf=${1}
	local threads=${2}
	local assembly=${3}
	local gnomAD_vcf_chr1=${4}
	local output_vcf=${5}

	if [[ -z ${output_table} ]]; then
		local output_table=${anno_table::-4}.gnomAD.tsv
	fi

	local -a main_chroms=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
	# Step 1: Split the input vcf to alt-contig-variants file and main-contig-variants file
	bcftools view -r $(echo ${main_chroms[@]} | tr ' ' ',') -Ou ${input_vcf} | \
	bcftools sort -Oz -o ${input_vcf/.vcf*/.primary.vcf.gz} && \
	tabix -f -p vcf ${input_vcf/.vcf*/.primary.vcf.gz} && \
	bcftools isec -C -Ou ${input_vcf} ${input_vcf/.vcf*/.primary.vcf.gz} | \
	bcftools sort -Oz -o ${input_vcf/.vcf*/.alt.vcf.gz} && \
	tabix -f -p vcf ${input_vcf/.vcf*/.alt.vcf.gz}

	# Step 2: Further split the main-contig-variants file to multiple single-contig variant files
	for chr in "${main_chroms[@]}"; do
		bcftools view -r chr${chr} ${input_vcf} -Oz -o ${input_vcf/.vcf*/.chr${chr}.vcf.gz}
	done

	# Step 3: Perform annotation with bcftools annotate to add the INFO fields from gnomAD vcfs to the input splitted vcfs
	export gnomAD_vcf_chr1
	parallel -j ${threads} "bcftools annotate -a ${gnomAD_vcf_chr1/.chr1.vcf.bgz/}.chr{}.vcf.bgz -c CHROM,POS,REF,ALT,.INFO/AC_joint,.INFO/AN_joint,.INFO/AF_joint,.INFO/nhomo_joint_XX,.INFO/nhomo_joint_XY -Oz -o ${input_vcf/.vcf*/}.chr{}.gnomAD.vcf.gz ${input_vcf/.vcf*/}.chr{}.vcf.gz" ::: "${main_chroms[@]}"

	# Merge the annotated vcfs together (including the alt-contig-variants file) to form into the output vcf
	# We can filter on AF here
	bcftools concat -Ou ${input_vcf/.vcf*/.alt.vcf.gz} ${input_vcf/.vcf*/.chr*.gnomAD.vcf.gz} | \
	bcftools sort -Ou - | \
	bcftools filter -e 'INFO/AF_joint >= 0.05' -Oz -o ${output_vcf} && \
	tabix -f -p vcf ${output_vcf}
}



function anno_clinvar_data () {
	local input_vcf=${1}
	local clinvar_vcf=${2}
	local output_vcf=${3}

	if [[ -z ${output_vcf} ]]; then
		local output_vcf=${input_vcf/.vcf/.clinvar.vcf}
	fi

	bcftools annotate -a ${clinvar_vcf} -c CHROM,POS,REF,ALT,.INFO/CLNDN,.INFO/CLNHGVS,.INFO/CLNREVSTAT,.INFO/CLNSIG,.INFO/GENEINFO -Ou ${input_vcf} | \
	bcftools sort -Oz -o ${output_vcf} && \
	tabix -f -p vcf ${output_vcf} && \
	display_vcf ${output_vcf}
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

