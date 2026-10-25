#!/usr/bin/env bash
# Immediate pre-change E8 serial control; only this function is overridden.
# annotation_vcf.sh SHA256: 0029d817e6b43792eed944b763094d8ff4acff4bb9872407c904a072dbdc5743
# anno_agg_gnomAD_data SHA256: 2fd4717defe5ed472ad4d27b698c56c6ab721711bcca4cfaf1f08ac870c959cb
# Retain exact commands, selection, guards, reuse and concat ordering.
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

    local publication_pending="$input_vcf.gnomAD.publish.pending"
    [[ ! -L $publication_pending && ( ! -e $publication_pending || -f $publication_pending ) ]] || return 1
    if [[ -e $publication_pending ]]; then
        anno_agg_gnomAD_index --input-vcf "$input_vcf" --force || return $?
        check_vcf_validity "$input_vcf" || return $?
        rm -- "$publication_pending" || return $?
    else
        anno_agg_gnomAD_index --input-vcf "$input_vcf" || return $?
    fi

    # Core frequency fields
    local core_fields=".INFO/AC_joint,.INFO/AN_joint,.INFO/AF_joint,.INFO/nhomalt_joint"

    # Maximum values across populations (no sex-specific versions available)
    local max_fields=",.INFO/AC_grpmax_joint,.INFO/AF_grpmax_joint,.INFO/AN_grpmax_joint,.INFO/nhomalt_grpmax_joint"

    # Sex-specific fields (overall). AC/AN by sex are required by the sex-aware
    # ACMG evidence gates added in commit 111952f (combine_annotations.py reads
    # AC_joint_XX/AN_joint_XX/AC_joint_XY/AN_joint_XY); without them pysam raises
    # "Invalid header" during prioritization.
    local sex_fields=",.INFO/AC_joint_XX,.INFO/AN_joint_XX,.INFO/AC_joint_XY,.INFO/AN_joint_XY,.INFO/AF_joint_XX,.INFO/AF_joint_XY,.INFO/nhomalt_joint_XX,.INFO/nhomalt_joint_XY"

    # Major population frequencies
    local -a pop_codes=("nfe" "eas" "afr" "amr" "asj" "fin" "sas" "mid" "remaining")
    local pop_af_fields=""
    local pop_nhomalt_fields=""
    for pop in "${pop_codes[@]}"; do
        pop_af_fields+=",.INFO/AF_joint_${pop}_XX,.INFO/AF_joint_${pop}_XY"
        pop_nhomalt_fields+=",.INFO/nhomalt_joint_${pop}_XX"
    done

    check_vcf_validity "${input_vcf}" || \
    {
        log "The input vcf ${input_vcf} is not valid. Quit now"
        return 1
    }

    # Header IDs are bare names; .INFO/ belongs only to bcftools -c columns.
    local readiness_fields="${core_fields}${max_fields}${sex_fields}${pop_af_fields}${pop_nhomalt_fields}"
    readiness_fields=${readiness_fields//.INFO\//}
    check_vcf_infotags "$input_vcf" "$readiness_fields" && \
    log "The input vcf ${input_vcf} already contains the INFO tags ${core_fields}${max_fields}${sex_fields}${pop_af_fields}${pop_nhomalt_fields}. We do not need to add them again" && \
    return 0 || \
    log "The input vcf ${input_vcf} does not contain the INFO tags ${core_fields}${max_fields}${sex_fields}${pop_af_fields}${pop_nhomalt_fields}. We need to add them"

    # We already make sure the input VCF is sorted and normalized
    # We also make sure no variants in alternative contigs are included in the input VCF
    # The used gnomAD vcf files should be the ones with joint AF information from both exome and genome datasets

    local chromosome_list
    local -a chr_chroms=()
    # A process substitution hides the producer's exit status from mapfile.
    chromosome_list=$(set -o pipefail; bcftools query -f '%CHROM\n' "${input_vcf}" | uniq -) || return $?
    [[ -z "${chromosome_list}" ]] || mapfile -t chr_chroms <<< "${chromosome_list}"
    log "The chromosome names extracted from the input vcf are: ${chr_chroms[*]}"

    # Step 2: Further split the main-contig-variants file to multiple single-contig variant files
    local -a valid_chroms=()
    for chr in "${chr_chroms[@]}"; do
        if [[ ${chr} =~ chr[0-9XY]+ ]]; then
            [[ ! -d ${input_vcf/.vcf*/.${chr}.vcf.gz.tbi} ]] || return 1
            if [[ ! -f ${input_vcf/.vcf*/.${chr}.vcf.gz} ]] || [[ ${input_vcf/.vcf*/.${chr}.vcf.gz} -ot ${input_vcf} ]]; then
                local split_stage split_output
                split_output="${input_vcf/.vcf*/.${chr}.vcf.gz}"
                split_stage=$(mktemp -d "$(dirname "${input_vcf}")/.gnomad-split.XXXXXX") || return $?
                (set -o pipefail
                 bcftools view -r "${chr}" -Ou "${input_vcf}" | \
                 bcftools sort -Oz -o "${split_stage}/${split_output##*/}") || return $?
                tabix -f -p vcf "${split_stage}/${split_output##*/}" || return $?
                check_vcf_validity "${split_stage}/${split_output##*/}" || return $?
                mv "${split_stage}/${split_output##*/}" "${split_output}" && \
                mv -T "${split_stage}/${split_output##*/}.tbi" "${split_output}.tbi" && \
                rmdir "${split_stage}" || return $?
            fi
            if [[ ! -f ${input_vcf/.vcf*/.${chr}.vcf.gz.tbi} ]] || [[ ${input_vcf/.vcf*/.${chr}.vcf.gz.tbi} -ot ${input_vcf/.vcf*/.${chr}.vcf.gz} ]]; then
                tabix -f -p vcf "${input_vcf/.vcf*/.${chr}.vcf.gz}" || return $?
            fi
            check_vcf_validity "${input_vcf/.vcf*/.${chr}.vcf.gz}" || return $?
            valid_chroms+=( "${chr}" )
        else
            log "WARNING!!! The chromosome ${chr} is not supported by gnomAD now. Currently only support chr1 to chrY. Skip this chromosome"
        fi
    done
    log "The valid chromosomes are: ${valid_chroms[*]}"
    [[ ${#valid_chroms[@]} -gt 0 ]] || {
        log "No supported chromosomes are available for gnomAD annotation."
        return 1
    }

    # Step 3: Perform annotation with bcftools annotate to add the INFO fields from gnomAD vcfs to the input splitted vcfs
    export gnomad_vcf_chrX
    local gnomad_vcf_suffix=${gnomad_vcf_chrX/*.chrX.vcf/}
    local annotation_job
    printf -v annotation_job 'bash %q anno_agg_gnomAD_chromosome --input-vcf %q.{}.vcf.gz --gnomad-vcf %q.{}.vcf%q --output-vcf %q.{}.gnomAD.vcf.gz --columns %q' \
        "${SELF_SCRIPT}" "${input_vcf/.vcf*/}" "${gnomad_vcf_chrX/.chrX.vcf*gz/}" \
        "${gnomad_vcf_suffix}" "${input_vcf/.vcf*/}" \
        "${core_fields}${max_fields}${sex_fields}${pop_af_fields}${pop_nhomalt_fields}"
    parallel -j "${threads}" --dry-run "${annotation_job}" ::: "${valid_chroms[@]}" || return $?
    local job_log=${input_vcf/.vcf*/.gnomAD.log}
    : > "${job_log}" || return $?
    parallel -j "${threads}" --joblog "${job_log}" "${annotation_job}" ::: "${valid_chroms[@]}" || return $?
    check_parallel_joblog --job-log "${job_log}" \
        --expected-jobs "${#valid_chroms[@]}" --exit-status 0 || return $?

    # Keep the original shell-glob concat order, but accept only the shards
    # selected for this invocation. Unrelated leftovers must not be published.
    local shard
    local -A expected_shards=()
    local -a annotated_vcfs=()
    for chr in "${valid_chroms[@]}"; do
        shard="${input_vcf/.vcf*/}.${chr}.gnomAD.vcf.gz"
        check_vcf_validity "${shard}" || return $?
        expected_shards["${shard}"]=1
    done
    for shard in "${input_vcf/.vcf*/}".chr*.gnomAD.vcf.gz; do
        [[ ${expected_shards["${shard}"]:-} == 1 ]] && annotated_vcfs+=( "${shard}" )
    done
    [[ ${#annotated_vcfs[@]} -eq ${#expected_shards[@]} ]] || {
        log "Expected gnomAD annotation shards are missing from the concat inputs."
        return 1
    }

    (set -o pipefail
     bcftools concat -Ou "${annotated_vcfs[@]}" | \
     bcftools sort -Oz -o "${output_vcf}") || return $?
    tabix -f -p vcf "${output_vcf}" || return $?
    check_vcf_validity "${output_vcf}" || return $?
    printf '%s\n' "$input_vcf" > "$output_vcf.publish.pending" && \
    mv "$output_vcf.publish.pending" "$publication_pending" && \
    mv "${output_vcf}" "${input_vcf}" && \
    mv -T "${output_vcf}.tbi" "${input_vcf}.tbi" && \
    rm -- "$publication_pending" && \
    rm -f -- "${annotated_vcfs[@]}" "${annotated_vcfs[@]/%/.tbi}"

}
