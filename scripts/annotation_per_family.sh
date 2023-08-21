#!/usr/bin/env bash
# This script is used to define all the annotation process functions and run the main pipeline
# For now, the script is only used to process short variants (Indels and SNVs). Large structural variants and CNVs are not included.

# The script expect a Family VCF input along with a pedigree table contains the Family pedigree information.
# Notice that the family sample ID should be consistent between the pedigree table and the VCF file.
# Notice that the proband of the family should stay as the first record of the family in the pedigree table

# For the same set of configurations (arguments), the pipeline should start with the position it ends last time, unless user specifically asked to forcefully rerun the pipeline

# Maintainer: yangyxt@gmail.com


FORCE_RERUN="False"
SELF_SCRIPT="$(realpath ${BASH_SOURCE[0]})"
SCRIPT_DIR="$(dirname "${SELF_SCRIPT}")"
BASE_DIR="$(dirname ${SCRIPT_DIR})"

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

    local TEMP
	log "Raw input arguments: $#"
    TEMP=$(getopt -o a:i:p:f:d:o:r: --long assembly:,input_vcf:,ped_file:,fam_name:,annovar_dir:,output_dir:,ref_genome: -- "$@")

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
		local ref_gen_extraction==$(basename ${ref_genome} | awk -F '.' '{printf "%s", $2;}')
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
	log "The VCF preprocessed for ANNOVAR annotation is ${pre_anno_vcf}."
	log "The output folder for annotation results is ${output_dir}"
	log "The ANNOVAR perl executable file is located in ${annovar_dir}"


	# Start to run ANNOVAR on the preprocessed VCF files. Predefine the output annovar results
	local annovar_tab=${output_dir}/${fam_name}.${assembly}_multianno.txt
	local annovar_vcf=${annovar_tab/.txt/.vcf}
	log "Expected ANNOVAR annotation table is ${annovar_tab} and the annotation VCF file is ${annovar_vcf}"

	run_ANNOVAR \
	-a ${assembly} \
	-p ${ped_file} \
	-f ${fam_name} \
	-i ${pre_anno_vcf} \
	-o ${output_dir} \
	-d ${annovar_dir} && \
	display_table ${annovar_tab} || { \
	log "Failed to perform ANNOVAR on ${pre_anno_vcf}. Quit now"
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
	bcftools view -r $(cat /data/liftover/ucsc_GRC.primary.contigs.tsv | tr '\n' ',') -Ou ${input_vcf/.vcf*/.sorted.vcf.gz} | \
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
		/data/liftover/ucsc_to_GRC.contig.map.tsv \
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
	bcftools filter -e 'FORMAT/DP[0] < 5' -Oz -o ${input_vcf/.vcf*/.filtered.vcf.gz} ${input_vcf/.vcf/.norm.vcf} && \
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


function run_ANNOVAR {
    local input_vcf=""
    local ped_file=""
	local fam_name=""
	local assembly="hg19"
    local annovar_dir=""
	local output_dir=""

    local TEMP
    TEMP=$(getopt -o a:i:p:f:d:o: --long assembly:,input_vcf:,ped_file:,fam_name:,annovar_dir:,output_dir: -- "$@")

    # if getopt failed, return an error
    [[ $? != 0 ]] && return 1

    eval set -- "$TEMP"
	log "TEMP: $TEMP"

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
    local output_table=${output_dir}/${fam_name}.${assembly}_multianno.txt
	local output_vcf=${output_dir}/${fam_name}.${assembly}_multianno.vcf
	local table_annovar=${annovar_dir}/table_annovar.pl

    # Check whether we can skip this function
    if [[ ${output_table} -nt ${input_vcf} ]] && \
       [[ $(check_table_lineno ${output_table}) -eq $(check_vcf_lineno ${input_vcf}) ]] && \
       [[ ${output_vcf} -nt ${input_vcf} ]]; then
        log "The output table is already existed and updated. Skip the following steps"
        return 0;
    fi

    # Check whether the input is valid
    if check_vcf_validity ${input_vcf}; then
        log "Input vcf ${input_vcf} is updated and valid"
    else
        log "Input vcf ${input_vcf} is not valid"
        return 1;
    fi

	local af_dbs
	if [[ ${assembly} == "hg19" ]]; then
		local af_dbs="gnomad211_exome,gnomad211_genome"
	elif [[ ${assembly} == "hg38" ]]; then
		local af_dbs="gnomad211_exome,gnomad312_genome"
	fi

	log "Start running ANNOVAR core program on the family VCF ${input_vcf}"
	time perl ${table_annovar} \
    ${input_vcf} \
    ${annovar_dir}/humandb/ \
    -buildver ${assembly} \
    -out ${output_dir}/${fam_name} \
    -remove \
    -protocol refGene,${af_dbs},dbnsfp42a,clinvar_20221231 \
    -operation g,f,f,f,f \
    -otherinfo \
    -nastring . \
    -vcfinput
	check_return_code

	log "Finish running ANNOVAR core program on the sample ${fam_name}:${proband_ID} "
	log "**********************This is the separate line of ANNOVAR main program and customed program***************************"
	log "proband_ID=${proband_ID}"

	if [[ $(check_table_lineno ${output_table}) -eq $(check_vcf_lineno ${output_vcf} 2> /dev/null) ]] && \
       [[ $(check_table_lineno ${output_table}) -eq $(check_vcf_lineno ${input_vcf}) ]]; then
        log "The annotation seems successfully finished"
	else
		log "The ANNOVAR annotation process seems problematic since the output table ${output_table} has different number of records with the input vcf file: ${input_vcf}"
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
		log "Directly run main_workflow with input args: $#"
		main_workflow "$@"
	fi
fi

