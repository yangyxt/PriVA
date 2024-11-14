#!/usr/bin/env bash

SELF=$(basename ${BASH_SOURCE[0]})
SELF_DIR=$(dirname ${SELF})
BASE_DIR=$(dirname ${SELF_DIR})
ARGPARSE=${SELF_DIR}/argparse.bash
# This script is just a function pool, used to store commonly shared functions


function randomID {
    dd bs=24 count=1 status=none if=/dev/urandom | base64 | tr +/ _
}


function log() {
    local msg="$1"
    local script_name="${BASH_SOURCE[1]##*/}"
    local func_name="${FUNCNAME[1]}"
    local line_num="${BASH_LINENO[0]}"
    local timestamp
	timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    >&2 echo "[$timestamp] [Script $script_name: Line $line_num] [Func: $func_name] $msg"
}


# Function to read YAML configuration
function read_yaml() {
    local yaml_file=$1
    local key=$2
    # Remove the quotes enclosing the value
    yq ".$key" "$yaml_file" | sed 's/^"\(.*\)"$/\1/'
}


# Function to update YAML configuration
function update_yaml() {
    local yaml_file=$1
    local key=$2
    local value=$3

    # Replace the value while preserving comments using sed
    sed -i "s/^\($key *: *\).*/\1$value/" "${yaml_file}"
}


function mamba_activate {
    local env_name=${1}
    mamba activate ${env_name} 2> /dev/null || { \
    local conda_path=$(dirname $(dirname ${CONDA_EXE})) && \
    source ${conda_path}/etc/profile.d/conda.sh && source ${conda_path}/etc/profile.d/mamba.sh && \
    mamba activate ${env_name} || \
    return 1; }
}


function mamba_deactivate {
    mamba deactivate 2> /dev/null || { \
    local conda_path=$(dirname $(dirname ${CONDA_EXE})) && \
    source ${conda_path}/etc/profile.d/conda.sh && source ${conda_path}/etc/profile.d/mamba.sh && \
    mamba deactivate || \
    return 1; }
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
    if [[ -z ${TMPDIR} ]]; then
        TMPDIR="/tmp"
    fi

    local tmp_tab="$TMPDIR/$(randomID).tsv"
    log "Using ${tmp_tab} as the temporary file to store the VCF records from ${input_vcf}"

    if [[ -z ${head_lines} ]]; then
        local head_lines=10
    else
        local head_lines=$((head_lines + 1))
    fi

    bcftools view -H ${input_vcf} | head -${head_lines} > ${tmp_tab} && \
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
    tabix -f -p vcf ${output_vcf/.${place_holder}*/.vcf.gz} && \
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


function check_vcf_validity {
    local input_vcf=${1}
    local expected_lines=${2}
    local expected_samples=${3}

    if [[ -z ${expected_lines} ]]; then expected_lines=1; fi
    if [[ ! -f ${input_vcf} ]]; then
        log "${input_vcf} does not even exist."
        return 10
    fi

    if [[ ${input_vcf} =~ gz$ ]]; then
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
        if [ $(bcftools view -H ${input_vcf} | head -${expected_lines} | wc -l | awk '{print $1}') -ge ${expected_lines} ]; then  # Check input_vcf content
            log "${input_vcf} has enough records."
        else
            log "${input_vcf} has $(bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ${input_vcf} | wc -l | awk '{print $1}') valid lines while expected to have ${expected_lines} lines."
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
        if [ $(bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ${input_vcf} | wc -l | awk '{print $1}') -ge ${expected_lines} ]; then
            log "${input_vcf} has enough records."
        else
            log "${input_vcf} has $(bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ${input_vcf} | wc -l | awk '{print $1}') valid lines while expected to have ${expected_lines} lines."
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

    bcftools view ${input_vcf} | head -200 > /dev/null && \
    log "The input VCF file ${input_vcf} has no format issue detected by bcftools" && \
    return 0 || \
    { log "The input VCF file ${input_vcf} has some format issue detected by bcftools. Quit this function with error."$'\n'; return 10; }
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

	[[ ! ${output_vcf} =~ \.vcf\.gz$ ]] && \
	log "Output vcf file ${output_vcf} is not ending with .vcf.gz. Only support .vcf.gz format now." && \
	return 1;

    bcftools norm -m -both -f ${ref_genome} --atom-overlaps "." --keep-sum AD -a ${input_vcf} | \
    bcftools norm -d exact - | \
    bcftools view -i 'ALT!="*" && GT="alt"' - | \
    bcftools sort -Oz -o ${output_vcf} && \
    tabix -f -p vcf ${output_vcf}

    display_vcf ${output_vcf}
    >&2 echo ""  # Append a new line to the end of the logs from this function
}


function check_table_lineno {
    local input_table=${1}

    awk -F '\t' 'NR > 1{print;}' ${input_table} | wc -l
}


function check_table_column {
    local input_table=${1}
    local column=${2}

    if [[ ! -f ${input_table} ]]; then
        log "Input table ${input_table} does not even exist."
        false;
    fi

    local -a col_inds=($(head -1 ${input_table} | \
                         awk -F '\t' -v col="${column}" '{for(i=1;i<=NF;i++) if ($i == col) printf "%s ", i;}'))

    if [[ ${#col_inds[@]} -eq 1 ]]; then
        true;
    elif [[ ${#col_inds[@]} -eq 0 ]]; then
        false;
    elif [[ ${#col_inds[@]} -gt 1 ]]; then
        log "Multiple columns have the same column label ${column}. The column indices are: ${col_inds[*]}"
        echo "${col_inds[*]}"
        true;
    fi
}


function check_vcf_lineno {
    local input_vcf=${1}
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${input_vcf} | wc -l
}


function check_vcf_contig_version () {
    local input_vcf=${1}

    local example_contig_name
	example_contig_name=$(bcftools query -f '%CHROM\n' ${input_vcf} | head -1)

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


function liftover_from_GRCh_to_ucsc () {
    # Only designed for lifting over GRCh37 to hg19 or GRCh38 to hg38
    local input_vcf=${1}
    local contig_map=${2}
    local output_vcf=${3}
    [[ -z ${contig_map} ]] && local contig_map=${BASE_DIR}/data/liftover/GRC_to_ucsc.contig.map.tsv
    [[ -z ${output_vcf} ]] && local output_vcf=${input_vcf/.vcf/.addchr.vcf}
    [[ ! ${output_vcf} =~ \.vcf\.gz$ ]] && local output_vcf=${output_vcf/.vcf*/.vcf.gz}

    # Upon test, we do not need to escape < in awk regex
    bcftools annotate --rename-chrs ${contig_map} ${input_vcf} -Ou | \
    bcftools sort -Oz -o ${output_vcf} && \
    tabix -f -p vcf ${output_vcf} && \
    log "The output vcf file is ${output_vcf}" && \
    display_vcf ${output_vcf}
}


function liftover_from_ucsc_to_GRCh () {
    local input_vcf=${1}
    local contig_map=${2}
    local output_vcf=${3}
    [[ -z ${contig_map} ]] && local contig_map=$(dirname ${SELF_DIR})/data/liftover/ucsc_to_GRC.contig.map.tsv
    [[ -z ${output_vcf} ]] && local output_vcf=${input_vcf/.vcf/.rmchr.vcf}
    [[ ! ${output_vcf} =~ \.vcf\.gz$ ]] && local output_vcf=${output_vcf/.vcf*/.vcf.gz}

    bcftools annotate --rename-chrs ${contig_map} ${input_vcf} -Ou | \
    bcftools sort -Oz -o ${output_vcf} && \
    tabix -f -p vcf ${output_vcf} && \
    log "The output vcf file is ${output_vcf}" && \
    display_vcf ${output_vcf}
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


function check_return_code {
    local return_code=$(echo $?)
    log ": The last step's return code is ${return_code}";
    if [[ ${return_code} -gt 0 ]]; then
        log "The last step does not finish in normal way. Exiting the whole script now." && \
        exit 1
    else
        log ": The last step finished properly. Continue";
    fi
}



function crossmap_liftover_hg382hg19 {
    local chain_file=""
    local input_vcf=""
    local output_vcf=""
    local hg19_fasta=""

    local TEMP
    TEMP=$(getopt -o hc:i:o:f: --long help,chain_file:,input_vcf:,output_vcf:,hg19_fasta: -- "$@")

    if [[ $? != 0 ]]; then return 1; fi

    eval set -- "$TEMP"

    while true; do
        case "$1" in
            -h|--help)
                echo "Usage: crossmap_liftover_hg382hg19 [options]"
                echo "Options:"
                echo "  -c, --chain_file        Path to the chain file for liftover"
                echo "  -i, --input_vcf         Path to the input VCF file"
                echo "  -o, --output_vcf        Path to the output VCF file (default: <input_vcf>.hg19.vcf)"
                echo "  -f, --hg19_fasta        Path to the hg19 fasta file"
                return 0
                ;;
            -c|--chain_file)
                chain_file="$2"
                log "The chain file is ${chain_file}"
                shift 2
                ;;
            -i|--input_vcf)
                input_vcf="$2"
                log "The input vcf is ${input_vcf}"
                shift 2
                ;;
            -o|--output_vcf)
                output_vcf="$2"
                log "The output vcf is ${output_vcf}"
                shift 2
                ;;
            -f|--hg19_fasta)
                hg19_fasta="$2"
                log "The hg19 fasta is ${hg19_fasta}"
                shift 2
                ;;
            --)
                shift
                break
                ;;
            *)
                echo "Invalid option: $1" >&2
                return 1
                ;;
        esac
    done


    if [[ ${input_vcf} =~ \.bgz$ ]]; then
        bcftools view -Oz -o ${input_vcf/.bgz/.gz} ${input_vcf} && \
        display_vcf ${input_vcf/.bgz/.gz} && \
        local input_vcf=${input_vcf/.bgz/.gz} && \
        bcftools index -f -t ${input_vcf} && \
        ls -lh ${input_vcf}*
    fi

    if [[ -z ${output_vcf} ]]; then
        local output_vcf=${input_vcf/.vcf/.hg19.vcf}
    fi

    if [[ ! ${CONDA_PREFIX} =~ acmg ]]; then
        mamba_activate acmg
    fi


    CrossMap vcf \
    --chromid l \
    ${chain_file} \
    ${input_vcf} \
    ${hg19_fasta} \
    ${output_vcf/.gz/} && \
    bcftools sort -Oz -o ${output_vcf} ${output_vcf/.gz/} && \
    bcftools index -f -t ${output_vcf} && \
    display_vcf ${output_vcf} && \
    rm ${output_vcf/.gz/}
}



function bcftools_concatvcfs {
    local OPTIND v o e s a t m
    while getopts v:o::e::s::a::t::m:: args
    do
        case ${args} in
            v) local input_vcfs=$OPTARG ;;
            o) local merged_vcf=$OPTARG ;;
            e) local ignore_error=$OPTARG ;;
            a) local other_args=$OPTARG ;;
            t) local threads=$OPTARG ;;
            s) local samples=$OPTARG ;;  # Should be delimited by comma
            m) local tmp_dir=$OPTARG ;;
            *) echo "No argument passed, Pls at least pass sample path." ;;
        esac
    done

    if [[ -z ${tmp_dir} ]]; then
        local tmp_dir=${TMPDIR}
        [[ -z ${tmp_dir} ]] && local tmp_dir=/tmp
    fi
    log "The tmp dir is ${tmp_dir}"

    if [[ ${input_vcfs}  =~ \/ ]] && [[ ! ${input_vcfs} =~ \.vcf(\.[b]*gz)*$ ]] && [[ ! ${input_vcfs} =~ , ]]; then
        local -a vcfs=($(awk '{printf "%s ", $1;}' < ${input_vcfs}))
    else
        local -a vcfs=($(echo ${input_vcfs} | awk 'BEGIN{FS=",";} {for(i=1;i<=NF;i++) printf $i" ";}'))
    fi

    if [[ -z ${threads} ]]; then
        local threads=1
    fi

    log "the vcfs are ${vcfs[*]}"
    log "The merged vcf is ${merged_vcf}"
    # Check file existence.
    if [[ -z ${ignore_error} ]]; then
        local -a invalid_vcfs
        for vcf in "${vcfs[@]}"; do
            if check_vcf_validity ${vcf} 1 2> /dev/null; then
                log "To be merged ${vcf} is valid."
            else
                log "${vcf} not existed or corrupted. Run check_vcf_validity ${vcf} to see for yourself"
                invalid_vcfs+=( ${vcf} )
            fi
            if [[ ${vcf} =~ \.vcf$ ]]; then
                log "Since ${vcf} is plain text format and bcftools -f need to use bgzipped vcfs, we compress the vcf and index it with tabix."
                bcftools sort --temp-dir ${tmp_dir} -Oz -o ${vcf}.gz ${vcf} && tabix -f -p vcf ${vcf}.gz && \
                ls -lh ${vcf}.gz
            fi
            if [[ ! -z ${samples} ]]; then
                bcftools view -s "${samples}" -Oz -o ${vcf/.vcf*/.samp.vcf.gz} ${vcf/.vcf*/.vcf.gz} && \
                ls -lht ${vcf/.vcf*/.samp.vcf.gz} && \
                check_vcf_validity ${vcf/.vcf*/.samp.vcf.gz} 1 && \
                mv ${vcf/.vcf*/.samp.vcf.gz} ${vcf/.vcf*/.vcf.gz} && \
                tabix -f -p vcf ${vcf/.vcf*/.vcf.gz}
            fi
        done

        if [[ ${#invalid_vcfs[@]} -gt 0 ]]; then
            log "The following vcfs are invalid: ${invalid_vcfs[*]}"
            return 1
        fi
    fi

    local tmp_file_list=${tmp_dir}/$(randomID).lst
    echo "${vcfs[*]}" | \
    awk -F '\t' 'BEGIN{RS=" ";} length($1) > 0{gsub(/\n$/, ""); if ($1 ~ /\.gz$/) printf "%s\n", $1; else printf "%s.gz\n", $1;}' > ${tmp_file_list}
    log "Here is the temp list file storing the paths of to be concat vcfs:"
    ls -lh ${tmp_file_list}
    cat ${tmp_file_list}


    if [[ $(cat ${tmp_file_list} | wc -l) -eq 1 ]]; then
        cp -f $(cat ${tmp_file_list}) ${merged_vcf} && \
        bcftools index -f -t ${merged_vcf} && \
        ls -lh ${merged_vcf}*
    else
        # If using file list for input a list of vcfs, each one of them need to be bgzipped and tabix indexed
        bcftools concat -o ${merged_vcf} ${other_args} -a --threads ${threads} -d exact -Oz -f ${tmp_file_list} && \
        tabix -f -p vcf ${merged_vcf} && \
        bcftools sort --temp-dir ${tmp_dir} -o ${merged_vcf/.vcf.gz/.sorted.vcf.gz} -Oz ${merged_vcf} && \
        mv ${merged_vcf/.vcf.gz/.sorted.vcf.gz} ${merged_vcf} && \
        bcftools index -f -t ${merged_vcf} && \
        ls -lh ${merged_vcf}*
    fi

    display_table ${merged_vcf} 20
}


# Simplified helper function to check path existence
function check_path() {
    local path="$1"
    local type="$2"  # "file" or "dir"
    local var_name="$3"

    if [[ -z "$path" ]]; then
        log "Error: $var_name not specified (in either command line or config)"
        return 1
    fi

    if [[ "$type" == "file" && ! -f "$path" ]]; then
        log "Error: $var_name file not found: $path"
        return 1
    elif [[ "$type" == "dir" && ! -d "$path" ]]; then
        log "Error: $var_name directory not found: $path"
        return 1
    fi

    return 0
}


function parse_ped_file() {
	# Extract the proband ID and all the patient IDs from the pedigree file
	local ped_file=${1}
	local fam_vcf=${2}
	local fam_name=${3}

	# Extract the proband ID, usually it is the first sample in the VCF file
	local proband_ID
	proband_ID=$(bcftools query -l ${fam_vcf} | head -1) || {
			log "Failed to get proband ID from ${input_vcf}. Quit now"
			return 1
		}

	# Extract the family name from the pedigree file
	if [[ -z ${fam_name} ]]; then
		fam_name=$(awk -F '\t' '$2 == "'${proband_ID}'"{printf "%s", $1; exit 0;}' ${ped_file}) || {
			log "Failed to get family name from ${ped_file} using proband ID ${proband_ID}. Quit now"
			return 1
		}
	fi

	# Extract all the patient IDs from the pedigree file
	local -a patient_IDs
	mapfile -t patient_IDs < <(awk -F '\t' '$1 == "'${fam_name}'" && $6 == "2" {printf "%s\n", $2}' "${ped_file}") || {
		log "Failed to get patient IDs from ${ped_file} for family ${fam_name}. Quit now"
		return 1
	}

	echo "${fam_name}"
	echo "${proband_ID}"
	echo "${patient_IDs[@]}"
}


function check_parallel_joblog {
    local ret_code=$(echo $?)
    local job_log=${1}

	# compatible with the GNU parallel joblog format
    # print failed job commands to stdout

    local job_num=$(tail -n +2 ${job_log} | wc -l)
    local -a fail_job_ids=($(tail -n +2 ${job_log} | awk -F '\t' '$7 != "0"{print $1;}'))
    log ": Parallel job log is $(ls -lh ${job_log})"
    if [[ ${#fail_job_ids[@]} -gt 0 ]]; then
        >&2 cat ${job_log}
        log ": According to ${job_log}: totally ${job_num} jobs were running in parallel, ${#fail_job_ids[@]} of the parallel jobs failed, get a return code of ${ret_code}, the commands are:"
        for id in "${fail_job_ids[@]}"; do
            >&2 awk -F '\t' '$1 == "'${id}'"{printf "%s\n", $NF;}' ${job_log}
        done
    else
        log ": All parellel jobs are finished without an error."
        >&2 cat ${job_log}
    fi

	return ${#fail_job_ids[@]}
}


function extract_assembly_from_fasta() {
	local fasta_file=${1}
	# Here expecting the specified reference genome to have a name like ucsc.hg19.fasta or ucsc.hg38.fasta
	assembly=$(basename ${fasta_file} | awk -F '.' '{printf "%s", $2;}')
	[[ ! ${assembly} =~ ^hg[13][98]$ ]] && \
	log "Failed to extract the assembly version either from the fasta file name. Quit now." && \
	return 1
}


function check_vcf_infotags() {
    # Check if one or more INFO tags are defined in VCF header
    # Arguments:
    #   $1 - Input VCF file
    #   $2 - Comma-separated list of INFO tags to check
    # Returns:
    #   0 if all tags are found
    #   1 if any tag is missing or input is invalid

    local input_vcf="$1"
    local -a tags=($(echo "$2" | tr ',' ' '))

    # Validate inputs
    if [[ ! -f ${input_vcf} ]]; then
        log "Input VCF file ${input_vcf} does not exist"
        return 1
    fi

    if [[ ${#tags[@]} -eq 0 ]]; then
        log "No INFO tags specified to check"
        return 1
    fi

    # Extract INFO tags from header
    local header_tags=$(bcftools view -h "${input_vcf}" | grep -E "^##INFO=<ID=" | cut -d',' -f1 | cut -d'=' -f3)

    # Check each requested tag
    local missing_tags=()
    for tag in "${tags[@]}"; do
        if ! echo "${header_tags}" | grep -q "^${tag}$"; then
            missing_tags+=("${tag}")
        fi
    done

    if [[ ${#missing_tags[@]} -gt 0 ]]; then
        log "The following INFO tags are not defined in ${input_vcf} header: ${missing_tags[*]}"
        return 1
    else
        log "All requested INFO tags are present in ${input_vcf} header: ${tags[*]}"
        return 0
    fi
}



if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
    "$@"
fi
