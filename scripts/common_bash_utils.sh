#!/usr/bin/env bash


# This script is just a function pool, used to store commonly shared functions
function log() {
    local msg="$1"
    local script_name="${BASH_SOURCE[1]##*/}"
    local func_name="${FUNCNAME[1]}"
    local line_num="${BASH_LINENO[0]}"
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    >&2 echo "[$timestamp] [Script $script_name: Line $line_num] [Func: $func_name] $msg"
}


function randomID {
    dd bs=24 count=1 status=none if=/dev/urandom | base64 | tr +/ _
}


function display_table {
    if [[ -z $2 ]]; then local rows=10; else local rows=${2}; fi
    if [[ -z $3 ]]; then local delimiter="\t"; else local delimiter=${3}; fi

    if [[ ${delimiter} == "\t" ]]; then
        local del_arg=""
    else
        local del_arg="-d ${delimiter}"
    fi

    if [[ ${1} =~ \.vcf$ ]] || \
       [[ ${1} =~ \.vcf\.gz$ ]] || \
       [[ ${1} =~ \.bcf$ ]]; then
        log "Input an VCF file ${1}. Using bcftools to extract the records and view the content of first ${rows} lines."
        display_vcf ${1} ${rows}
        return;
    fi

    local row_num=$(tail -n +2 ${1} | wc -l)
    local col_num=$(head -1 ${1} | awk '{print NF;}')

    if [[ ${1} =~ \.gz$ ]]; then
        local tmp_tag=$(randomID)
        zcat ${1} > ${1/.gz/}.${tmp_tag} && \
        local input=${1/.gz/}.${tmp_tag}
    else
        local input=${1}
    fi

    
    log "${1} has ${row_num} rows and ${col_num} columns. It looks like:"
    if [[ ${rows} -le 0 ]]; then
        tsv-pretty -u ${del_arg} -m 5000 -l 200 -a ${input} >&2 2> /dev/null 
    else
        tsv-pretty -u ${del_arg} -m 5000 -l 200 -a ${input} | \
        head -n ${rows} - >&2 2> /dev/null || >&2 echo ""
    fi

    if [[ ${input} != ${1} ]]; then
        silent_remove_tmps ${input}
    fi
}


function display_vcf {
    local input_vcf=${1}
    local head_lines=${2}
    local tmp_tab="$TMPDIR/$(randomID).tsv"

    if [[ -z ${head_lines} ]]; then
        local head_lines=10
    else
        local head_lines=$((head_lines + 1))
    fi

    bcftools view -H ${input_vcf} | tail -n +${head_lines} > ${tmp_tab} && \
    display_table ${tmp_tab} && \
    silent_remove_tmps ${tmp_tab}
}


function silent_remove_tmps {
    local -a targets=($(echo "$@"))

    for target in "${targets[@]}"; do
        if [[ -f ${target} ]]; then
            rm ${target} || \
			log "File ${target} failed to be deleted (either its not existed or occupied)"
        elif [[ -d ${target} ]]; then
            if [[ ${target} != $TMPDIR ]]; then
                rm -r ${target} || \
				log "Folder ${target} failed to be deleted (either its not existed or occupied)"
            fi
        fi
    done
}


function filter_vcf_by_GT() {
    # This function will filter a multi-sample VCF file to exclude variants
    # where a specified subset of samples contain all ref calls or missing calls.
    # Arguments:
    #   $1 - VCF file to filter
    #   $2 - An array containing sample names, delimited by comma

    local vcf_file="$1"
    local -a sample_array=($(echo ${2} | awk 'BEGIN {RS=ORS=",";} {printf "%s ", $1;}')) # Access the array passed as argument
    local output_vcf="$3"

	local tmp_sample_list=$TMPDIR/$(randomID).txt

	if [[ ${output_vcf} =~ \.vcf ]]; then
		local place_holder="vcf"
	elif [[ ${output_vcf} =~ \.bcf ]]; then
		local place_holder="bcf"
	else
		log "Cannot identify the specified output VCF format: ${output_vcf}" 
	fi

    printf "%s\n" "${sample_array[@]}" > ${tmp_sample_list}
    
	# First extract the subset samples to a vcf file and make sure all samples are not missing or homo ref
	bcftools view -S ${tmp_sample_list} -i 'GT="alt"' ${vcf_file} -Ou | \
	bcftools sort -Oz -o ${output_vcf/.${place_holder}*/.subset.vcf.gz} && \
	tabix -f -p vcf ${output_vcf/.${place_holder}*/.subset.vcf.gz}

	bcftools annotate -a ${output_vcf/.${place_holder}*/.subset.vcf.gz} -m VALID_REC -Ou ${vcf_file} | \
	bcftools filter -i 'INFO/VALID_REC=1' -Ou | \
	bcftools sort -Oz -o ${output_vcf/.${place_holder}*/.vcf.gz} && \
	consistent_vcf_format ${output_vcf} ${output_vcf/.${place_holder}*/.vcf.gz} && \
    silent_remove_tmps ${tmp_sample_list}
}


function announce_remove_tmps {
    local -a targets=($(echo "$@"))

    for target in "${targets[@]}"; do
        if [[ -f ${target} ]]; then
            rm ${target} && \
			log "${target} has been succesfully removed" || \
			log "File ${target} failed to be deleted (either its not existed or occupied)"
        elif [[ -d ${target} ]]; then
            if [[ ${target} != $TMPDIR ]]; then
                rm -r ${target} && \
				log "${target} has been succesfully removed" || \
				log "Folder ${target} failed to be deleted (either its not existed or occupied)"
            fi
        fi
		>&2 echo ""
    done
}



function consistent_vcf_format () {
	local input_vcf=${1}
	local output_vcf=${2}

	# Convert output_vcf format to make it the same format with the input vcf
	if [[ ${input_vcf} == ${output_vcf} ]]; then
		true
	fi

	if [[ ${input_vcf/*.vcf/vcf} == ${output_vcf/*.vcf/vcf} ]]; then
		:
	elif [[ ${input_vcf/*.bcf/bcf} == ${output_vcf/*.bcf/bcf} ]]; then
		:
	elif [[ ${input_vcf} =~ \.vcf\.gz$ ]]; then
		local op_format="z"
		local op_suffix=".vcf.gz"
	elif [[ ${input_vcf} =~ \.vcf$ ]]; then
		local op_format="v"
		local op_suffix=".vcf"
	elif [[ ${input_vcf} =~ \.bcf$ ]]; then
		local op_format="b"
		local op_suffix=".bcf"
	else
		log "The input VCF file ${input_vcf} does not have a recognizable format judging by its file name suffix"
		return 10
	fi

	if [[ -z ${op_format} ]]; then
		if [[ ${input_vcf} == ${output_vcf} ]]; then
			return
		else
			mv ${input_vcf} ${output_vcf}
			mv ${input_vcf}.tbi ${output_vcf}.tbi 2> /dev/null || :
		fi
	fi

	if [[ ${output_vcf} =~ \.bcf$ ]]; then
		local ori_output_format=bcf
	elif [[ ${output_vcf} =~ \.vcf ]]; then
		local ori_output_format=vcf
	else
		log "The input output VCF file ${output_vcf} does not have a recognizable format judging by its file name suffix"
		return 10
	fi

	local final_output_vcf=${output_vcf/.${ori_output_format}*/${op_suffix}}
	bcftools sort ${output_vcf} | \
	bcftools view -O${op_format} --no-version -o ${final_output_vcf}
	
	if [[ ${final_output_vcf} =~ \.vcf\.gz$ ]]; then
		tabix -f -p vcf ${final_output_vcf}
	fi

	echo ${final_output_vcf}
}


function check_vcf_validity {
    local input_vcf=${1}
    local expected_lines=${2}
    local expected_samples=${3}

    if [[ -z ${expected_lines} ]]; then expected_lines=1; fi
    if [[ ! -f ${input_vcf} ]]; then
        log "${input_vcf} does not even exist."
        return 10
    fi

    if [[ ${input_vcf} =~ \.gz$ ]]; then
        log "${input_vcf} should be bgzipped, check gz vcf validity"
        if check_gz_vcf ${input_vcf} ${expected_lines}; then
            if check_vcf_samples ${input_vcf} ${expected_samples}; then
                log "${input_vcf} has solid sample names as expected."
            else
                return 10
            fi
        elif check_plain_vcf ${input_vcf::-3} ${expected_lines} && [[ ${input_vcf::-3} -nt ${input_vcf} ]]; then
            if check_vcf_samples ${input_vcf} ${expected_samples}; then
                bgzip -f -c ${input_vcf::-3} > ${input_vcf} && \
                tabix -f -p vcf ${input_vcf} && return || \
                return 10
            else
                return 10
            fi
        else
            return 10
        fi
    elif [[ ${input_vcf} =~ \.vcf$ ]]; then
        log "${input_vcf} should be plain text, check plain vcf validity"
        if check_plain_vcf ${input_vcf} ${expected_lines}; then
            if check_vcf_samples ${input_vcf} ${expected_samples}; then
                log "${input_vcf} has solid sample names as expected."
            else
                return 10
            fi
        elif check_gz_vcf ${input_vcf}.gz ${expected_lines} && [[ ${input_vcf}.gz -nt ${input_vcf} ]]; then
            if check_vcf_samples ${input_vcf} ${expected_samples}; then
                gunzip -f -c ${input_vcf}.gz > ${input_vcf} && return || \
                return 10
            else
                return 10
            fi
        else
            return 10
        fi
    else
        return 20
    fi
}


function check_gz_vcf {
    local input_vcf=${1}
    local expected_lines=${2}
    if [[ -z ${expected_lines} ]]; then expected_lines=1; fi
    if [[ ! -f ${input_vcf} ]]; then return 10; fi

    if check_vcf_format ${input_vcf}; then  # Check input_vcf format
        log "${input_vcf} has solid vcf format. Check whether contains enough record number"
        if [ $(zcat ${input_vcf} | mawk '$0 !~ /^#/{print;}' | wc -l | awk '{print $1}') -ge ${expected_lines} ]; then  # Check input_vcf content
            log "${input_vcf} has enough records."
        else
            log "${input_vcf} has $(zcat ${input_vcf} | mawk '$0 !~ /^#/{print;}' | wc -l | awk '{print $1}') valid lines while expected to have ${expected_lines} lines."
            return 30
        fi
    else
        return 20
    fi
}


function check_plain_vcf {
    local input_vcf=${1}
    local expected_lines=${2}
    if [[ -z ${expected_lines} ]]; then expected_lines=1; fi
    if [[ ! -f ${input_vcf} ]]; then return 10; fi

    if check_vcf_format ${input_vcf}; then
        log "${input_vcf} has solid vcf format. Check whether contains enough record number"
        if [ $(mawk '$0 !~ /^#/{print;}' ${input_vcf} | wc -l | awk '{print $1}') -ge ${expected_lines} ]; then
            log "${input_vcf} has enough records."
        else
            log "${input_vcf} has $(mawk '$0 !~ /^#/{print;}' ${input_vcf} | wc -l | awk '{print $1}') valid lines while expected to have ${expected_lines} lines."
            return 30
        fi
    else
        return 20
    fi
}


function check_vcf_samples {
    local input_vcf=${1}
    local expected_samples=${2}
    # expected samples should be delimited by comma

    if [[ ${input_vcf} =~ \.vcf$ ]]; then
        bgzip -f -c ${input_vcf} > ${input_vcf/.vcf/.tmp.vcf.gz} && \
        tabix -f -p vcf ${input_vcf/.vcf/.tmp.vcf.gz} && \
        local gz_vcf=${input_vcf/.vcf/.tmp.vcf.gz} && \
        local asks=$(bcftools query -l ${gz_vcf} | sort - | awk '{printf "%s\n", $1;}') && \
        rm -f ${gz_vcf} || { log "Failed to compress ${input_vcf}" && return 10; }
        if [[ $? != "0" ]]; then return 10; fi
    else
        local gz_vcf=${input_vcf} && \
        local asks=$(bcftools query -l ${gz_vcf} | sort - | awk '{printf "%s\n", $1;}')
    fi

    if [[ -z ${expected_samples} ]]; then
        log "No expected samples specified. Quit checking samples"
    else
        local expected_samples=$(echo ${expected_samples} | awk 'BEGIN{RS=",";} {printf "%s\n", $1;}' | sort - | awk '{printf "%s\n", $1;}')
        if [[ ${expected_samples} == ${asks} ]]; then
            log "Expected samples are identical with actual samples in vcf file ${input_vcf}"
        else
            log "${input_vcf} has these samples:"
            log "${asks}"
            log "While ${input_vcf} is expected to have these samples:"
            log "${expected_samples}"
            return 10
        fi
    fi
}


function check_vcf_format {
    local input_vcf=${1}
    if [[ ! -f ${input_vcf} ]]; then
		log "Input VCF file ${input_vcf} seems not exist"
		return 10
	fi
    
	bcftools view ${input_vcf} > /dev/null || { \
	log "The input VCF file ${input_vcf} has some format issue detected by bcftools. Quit this function with error."$'\n'; \
	return 10; }

    # If VCF file does not have variant format issues. Check if it has meaningless variant records. 
    filter_out_noalleles ${input_vcf}
}


function filter_out_noalleles {
    local input_vcf=${1}
    local tmp_tag=$(randomID)
    local empty_var_lines=$(bcftools query -f '%ALT{0}\n' ${input_vcf} | awk '($1 == ".") || ($1 ~ /^[ ]*$/){print;}' - | wc -l)

    if [[ ${empty_var_lines} -eq 0 ]]; then
        log "Input vcf format is not recognized to have empty variant records: ${input_vcf}. Return with success"
        return 0;
    fi

    if [[ ${input_vcf} =~ \.vcf$ ]]; then
        bcftools filter -e 'ALT=="."' -Ov -o ${input_vcf/.vcf/}.${tmp_tag}.vcf ${input_vcf} && \
        if [[ $(sha256sum ${input_vcf} | awk '{printf "%s", $1;}') != $(sha256sum ${input_vcf/.vcf/}.${tmp_tag}.vcf | awk '{printf "%s", $1;}') ]]; then
            mv ${input_vcf/.vcf/}.${tmp_tag}.vcf ${input_vcf}
        else
            log "No empty variant records occurred in ${input_vcf}"
        fi
    elif [[ ${input_vcf} =~ \.vcf\.gz$ ]]; then
        bcftools filter -e 'ALT=="."' -Ov ${input_vcf} | \
        bcftools sort -Oz -o ${input_vcf/.vcf/}.${tmp_tag}.vcf.gz - && \
        if [[ $(sha256sum ${input_vcf} | awk '{printf "%s", $1;}') != $(sha256sum ${input_vcf/.vcf/}.${tmp_tag}.vcf.gz | awk '{printf "%s", $1;}') ]]; then
            mv ${input_vcf/.vcf/}.${tmp_tag}.vcf.gz ${input_vcf} && \
            tabix -f -p vcf ${input_vcf}
        else
            log "No empty variant records occurred in ${input_vcf}"
        fi
    elif [[ ${input_vcf} =~ \.bcf$ ]]; then
        bcftools filter -e 'ALT=="."' -Ov ${input_vcf} | \
        bcftools sort -Ob -o ${input_vcf/.vcf/}.${tmp_tag}.bcf - && \
        if [[ $(sha256sum -b ${input_vcf} | awk '{printf "%s", $1;}') != $(sha256sum -b ${input_vcf/.vcf/}.${tmp_tag}.bcf | awk '{printf "%s", $1;}') ]]; then
            mv ${input_vcf/.vcf/}.${tmp_tag}.bcf ${input_vcf}
        else
            log "No empty variant records occurred in ${input_vcf}"
        fi
    else
        log "Input vcf format is not recognized by its suffix: ${input_vcf}"
        return 10;
    fi
}


function contain_only_variants {
    local input_vcf=${1}
    local sample_ID=${2}

    # Test whether the input sample contain only variants (No ref calls or missing calls) in the given VCF file for a given sample.

    if [[ -z ${sample_ID} ]];then
        local no_var_lines=$(bcftools query -f '[%GT\t]\n' ${input_vcf} | uniq - | sort - | uniq - | awk -F '\t' '$1 ~ /^[\.0][\|\/][\.0]/{print;}' | wc -l)
        local sample_ID=$(bcftools query -l ${input_vcf} | tail -n +2)
    else
        local samp_ind=$(bcftools query -l ${input_vcf} | awk '$1 == "'${sample_ID}'"{printf "%s", NR;}')
        local no_var_lines=$(bcftools query -f "[%GT\t]\n" ${input_vcf} | uniq - | sort - | uniq - | awk -F '\t' '$'${samp_ind}' ~ /^[\.0][\|\/][\.0]/{print;}' | wc -l)
    fi

    if [[ ${no_var_lines} -gt 0 ]]; then
        log "VCF file ${input_vcf} has missing/ref GTs records for sample ${sample_ID}. Return false."
        false
    else
        log "VCF file ${input_vcf} has 0 missing/ref GTs records for sample ${sample_ID}. Return true."
        true
    fi
}


function check_vcf_multiallelics { 
    local input_vcf=${1}

    local multi_lines=$( \
    bcftools query -f '%ALT\n' ${input_vcf} | \
    mawk '$1 ~ /,/{print;}' | \
    wc -l ) 

    if [[ ${multi_lines} -gt 0 ]]; then
        log "${input_vcf} has ${multi_lines} lines of multi-allelic records, by default return false"
        false
    else
        log "${input_vcf} does not have multi-allelic records"
        true
    fi
}


function whether_same_varset {
	# Input can be TSV or VCF
	# Test whether two file contains the same set of variants (excluding the inspection of the GTs)
	local tab1=${1}
	local tab2=${2}

	local tab1_ids=$(grab_varids_sha256 ${tab1})
	local tab2_ids=$(grab_varids_sha256 ${tab2})

	if [[ ${tab1_ids} == ${tab2_ids} ]]; then
		log "${tab1} and ${tab2} share the same set of variants"
        return 0;
	elif [[ $(echo ${tab2_ids} | awk -F ',' '{printf "%s", $NF;}') -gt $(echo ${tab1_ids} | awk -F ',' '{printf "%s", $NF;}') ]]; then
		log "${tab1} has more variants than ${tab2}."
        return 20
	else
		log "${tab2} has different variants with ${tab1}."
        return 10
	fi
}


function grab_varids_sha256 {
    local file=${1}
    if [[ ! -f ${file} ]]; then
        log "${file} not even exists, so cannot fetch the variant ID sha256sum check code."
        echo "$(randomID),0"
    fi

    if [[ ${file} =~ \.vcf$ ]] || \
       [[ ${file} =~ \.vcf\.gz$ ]] || \
       [[ ${file} =~ \.bcf$ ]]; then
        local -a id_lst=($(bcftools query -f '%ID\n' ${file} | sort - | uniq -))
        if [[ $(echo ${id_lst} | mawk '{printf "%s", length($1);}') -gt 50 ]]; then
            local var_count=${#id_lst[@]}
            echo "${id_lst[@]}" | sha256sum - | mawk -v vc="${var_count}" '{printf "%s,%s", $1, vc;}'
        else
            local -a id_lst=($(bcftools query -f '%CHROM:%POS:%REF->%ALT\n' ${file} | sort - | uniq -))
            local var_count=${#id_lst[@]}
            echo "${id_lst[@]}" | sha256sum - | mawk -v vc="${var_count}" '{printf "%s,%s", $1, vc;}'
        fi
    elif check_text_or_binary ${file}; then
        local id_col=$(grab_column_label_index ${file} "uniq_ID")
        if [[ ! -z ${id_col} ]]; then
            local -a id_lst=($(cut -f ${id_col} ${file} | tail -n +2 - | sort - | uniq -))
            local var_count=${#id_lst[@]}
            echo "${id_lst[@]}" | sha256sum - | mawk -v vc="${var_count}" '{printf "%s,%s", $1, vc;}'
        else
            echo "$(randomID),0"
        fi
    else
        log "${file} seems to be a binary file, so cannot fetch the variant ID sha256sum check code."
        echo "$(randomID),0"
    fi
}


function normalize_vcf () {
	local input_vcf=${1}
	local output_vcf=${2}
	local ref_genome=${3}
	local tmp_dir=${4}

	log "Try to normalize the vcf file ${input_vcf}, and let the temporary folder set as ${tmp_dir}"

	bcftools norm -m -both -f ${ref_genome} --atom-overlaps "." --keep-sum AD -a ${input_vcf} | \
    bcftools norm -d exact - | \
    bcftools view -i 'ALT!="*"' - | \
	bcftools sort -Oz -o ${output_vcf/.vcf*/.vcf.gz} && \
	consistent_vcf_format ${output_vcf} ${output_vcf/.vcf*/.vcf.gz}
	
	display_vcf ${output_vcf}
	>&2 echo ""  # Append a new line to the end of the logs from this function
}


function install_annovar () {
	local assembly=""
    local parent_dir=""
	local download_url=""

    local TEMP
    TEMP=$(getopt -o a:p:u: --long assembly:,parent_dir:,download_url: -- "$@")

    # if getopt failed, return an error
    [[ $? != 0 ]] && return 1

    eval set -- "$TEMP"

    while true; do
        case "$1" in
            -a|--assembly)
                local assembly="$2"
                shift 2
                ;;
            -p|--parent_dir)
                local parent_dir="$2"
                shift 2
                ;;
			-u|--download_url)
                local download_url="$2"
                shift 2
                ;;	
            --)
                shift
                break
                ;;
            *)
                echo "Invalid option"
                return 2
                ;;
        esac
    done

	wget -O ${parent_dir}/annovar_latest.tar.gz ${download_url} && \
	tar -xvzf ${parent_dir}/annovar_latest.tar.gz && \
	local annovar_dir=${parent_dir}/annovar && \
	local table_annovar=${annovar_dir}/table_annovar.pl
	local main_annovar=${annovar_dir}/annotate_variation.pl

	# Start download databases, several databases are required for downstream analysis
	# Clinvar records (clinvar20221231)
	# refGene annotation ( refGene )
	# gnomAD_exome, gnomAD_genome (optional, can replace these annotations with annotations from an in-house prepared gnomAD dataset)
	# FATHMM, PROVEAN, MetaSVM, MetaLR, DANN ( combined in dbnsfp42a )

	if [[ -z ${assembly} ]]; then
		local -a assemblies=( "hg19" "hg38" )
		for assembly in "${assemblies[@]}"; do
			download_basic_annovar_resources ${main_annovar} ${assembly}
		done
	else
		download_basic_annovar_resources ${main_annovar} ${assembly}
	fi

	# These annotations are acquired elsewhere
	# Protein domain information, which will be added to the data directory in the acmg_auto github (from Prot2hg database, now the URL www.prot2hg.com not accessible)
	# CADD
	# gnomAD v2 and v3 combined
	# ClinVAR	
}


function download_basic_annovar_resources () {
	local annovar_pl=${1}
	local assembly=${2}
	local annovar_dir=$(dirname ${annovar_pl})

	# refGene
	# dbnsfp42a
	# clinvar20221231
	# gnomad_exome
	# gnomad_genome

	${annovar_pl} -downdb refGene ${annovar_dir}/humandb -buildver ${assembly} && \
	${annovar_pl} -downdb dbnsfp42a ${annovar_dir}/humandb -buildver ${assembly} -webfrom annovar && \
	${annovar_pl} -downdb gerp++gt2 ${annovar_dir}/humandb -buildver ${assembly} -webfrom annovar && \
	${annovar_pl} -downdb clinvar_20221231 ${annovar_dir}/humandb -buildver ${assembly}

	if [[ ${assembly} == "hg19" ]]; then
		${annovar_pl} -downdb gnomad211_exome ${annovar_dir}/humandb -buildver ${assembly} -webfrom annovar && \
		${annovar_pl} -downdb gnomad211_genome ${annovar_dir}/humandb -buildver ${assembly} -webfrom annovar
	elif [[ ${assembly} == "hg38" ]]; then
		${annovar_pl} -downdb gnomad211_exome ${annovar_dir}/humandb -buildver ${assembly} -webfrom annovar && \
		${annovar_pl} -downdb gnomad312_genome ${annovar_dir}/humandb -buildver ${assembly} -webfrom annovar
	else
		log "Failed to have a valid assembly version: ${assembly}"
	fi

	log "Succesfully download refGene, dbnsfp42a, clinvar20221231, gnomad_exome, gnomad_genome to the $(dirname ${annovar_pl})/humandb for assembly ${assembly}"
}


function check_table_lineno {
    local input_table=${1}

    awk -F '\t' 'NR > 1{print;}' ${input_table} | wc -l
}


function check_vcf_lineno {
    local input_vcf=${1}
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${input_vcf} | wc -l
}


function check_vcf_contig_version () {
	local input_vcf=${1}

	local example_contig_name=$(bcftools query -f '%CHROM\n' ${input_vcf} | head -1)

	if [[ ${example_contig_name} =~ ^chr ]]; then
		echo "ucsc"
	else
		echo "ncbi"
	fi
}


function check_vcf_contig_size { 
    local input_vcf=${1}
    local contig_name=${2}

    if [[ -z ${contig_name} ]]; then
        local contig_name="chrM"
    fi

    bcftools view -h ${input_vcf} | \
    grep "^##contig" | \
    grep "ID=${contig_name}" | \
    awk -F '[,=]' '{for(i=1;i<=NF;++i) if($i=="length") print $(i+1)}'
}


function check_vcf_assembly_version () {
	# INput VCF should contain UCSC contig names instead of GRC contig names
	local input_vcf=${1}
	if ! check_vcf_validity ${input_vcf}; then
		log "Cant locate the required contig name in the VCF header line to tell the assembly version of this VCF file ${input_vcf}"
		return 1;
	fi

	local chr1_size=$(check_vcf_contig_size ${input_vcf} "chr1")
	
	local random_contig_presence=$(bcftools view -h ${input_vcf} | \
					  			   awk -F '=' '$1 ~ /##contig/ && $3 ~ /chr7_gl000195_random/{print;}' | wc -l)

	local alt_contig_presence=$(bcftools view -h ${input_vcf} | \
					  			   awk -F '=' '$1 ~ /##contig/ && $3 ~ /chr6_GL000256v2_alt/{print;}' | wc -l)
	
	if [[ ${chr1_size} == "249250621" ]]; then
		echo "hg19"
	elif [[ ${chr1_size} == "248956422" ]]; then
		echo "hg38"
	else
		if [[ ${random_contig_presence} -ge 1 ]]; then
			echo "hg19"
		elif [[ ${alt_contig_presence} -ge 1 ]]; then
			echo "hg38"
		else
			log "Cant locate the required contig name in the VCF header line to tell the assembly version of this VCF file ${input_vcf}"
			echo ""
		fi
	fi
}


function liftover_from_GRCh_to_hg () {
	# Only designed for lifting over GRCh37 to hg19 or GRCh38 to hg38
    local input_vcf=${1}
    local contig_map=${2}
	local output_vcf=${3}
    local tmp_vcf=${input_vcf/.vcf*/.tmp.vcf}
    if [[ -z ${output_vcf} ]]; then local output_vcf=${input_vcf/.vcf/.addchr.vcf}; fi

	# Upon test, we do not need to escape < in awk regex
	bcftools view -Ov ${input_vcf} | \
	awk 'BEGIN{OFS=FS="\t";} \
		NR == FNR {arr[$2] = $1;} \
		NR > FNR && $0 ~ /^##contig/{old_contig = gensub(/##contig=<ID=([a-zA-Z_0-9\.]+),*(.+)>$/, "\\1", "g", $1); \
									if (old_contig !~ /^[0-9MTXT]$/) next; \
									else if (length(arr[old_contig]) == 0) new_contig = old_contig; \
									else new_contig = arr[old_contig]; \
									mod_line = gensub(/^(.+contig=<ID=)[a-zA-Z_0-9\.]+(,*.*>)/, "\\1"new_contig"\\2", "g", $0); \
									printf "%s\n", mod_line;} \
		NR > FNR && $0 !~ /^##contig/ && $0 ~ /^##/{print;} \
		NR > FNR && $0 ~ /^#CHROM/{print;} \
		NR > FNR && $0 !~ /^#/{ if (length(arr[$1]) == 0) new_contig = $1; \
								else new_contig = arr[$1]; \
								gsub(/.+/, new_contig, $1); \
								printf "%s\n", $0;}' ${contig_map} - > ${tmp_vcf}

	if [[ ${output_vcf} =~ \.[b]*gz$ ]]; then
		bgzip -f ${tmp_vcf} && \
		bcftools sort --temp-dir $TMPDIR -Oz -o ${tmp_vcf}.gz ${tmp_vcf}.gz && \
		tabix -f -p vcf ${tmp_vcf}.gz
		if check_vcf_validity ${tmp_vcf}.gz; then
			mv ${tmp_vcf}.gz ${output_vcf} && \
			mv ${tmp_vcf}.gz.tbi ${output_vcf}.tbi && \
			ls -lh ${output_vcf}*
		fi
	elif check_vcf_validity ${tmp_vcf}; then
		mv ${tmp_vcf} ${output_vcf} && \
		ls -lh ${output_vcf}*
	fi
}


function return_array_intersection {
    # Both arrays should not contain item values with space
    local -a array1=($(echo ${1}))
    local -a array2=($(echo ${2}))
    local -a special_char=("." "+" "?" "^" "$" "(" ")" "[" "]" "{" "}" "|" ":")
    local spe=0

    for item in "${array1[@]}"; do
        local new_item=${item}
        for char in "${special_char[@]}"; do
            if [[ ${item} == *"${char}"* ]]; then local spe=1 && local new_item=$(echo ${new_item} | awk '{gsub("\\'${char}'","\\'${char}'",$0); print;}'); fi
        done
        # if [[ ${spe} -gt 0 ]]; then echo "Line "${LINENO}": In function "${FUNCNAME}: After adding escape symbol to special characters, iterating item ${item} now looks like ${new_item}; fi
        if [[ ${array2[*]} =~ ${new_item} ]]; then
            result+=(${item})
        fi
    done
    echo "${result[*]}"
}


function get_array_index {
    # Only applied to non-associative array
    local -a values=($(echo ${1}))
    local -a array=($(echo ${2}))
    local -a indices

    for value in "${values[@]}"; do
        for i in "${!array[@]}"; do
            if [[ "${array[$i]}" == "${value}" ]]; then
                indices+=("${i}")
            fi
        done
    done

    echo "${indices[*]}"
}


function vcf_content_sha256 () {
	local input_vcf=${1}

	if check_vcf_validity ${input_vcf}; then
		bcftools view --no-version -Ov ${input_vcf} | \
		sha256sum | \
		awk '{printf "%s", $1;}'
	else
		echo $(randomID)
	fi
}



if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
    "$@"
fi
