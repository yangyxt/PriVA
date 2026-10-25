# Frozen v7 serial VEP body; exact source SHA256 5be1b9469c06d5ea519c9f66bed2ed286e092c52150c4dc442d979d9769c7f90.
# Function SHA256 f887aa2298360af1779eb9ad28a0fbd737093c66d73e10658af4950bd385415e.
# Source this after current annotation_vcf.sh only as the serial test oracle.
function _priva_anno_VEP_data() {
    # Private shared command body. Direct calls retain the runtime guard.
    local priva_vep_purpose=runtime
    [[ "${FUNCNAME[1]:-}" != anno_VEP_cache_candidate ]] || priva_vep_purpose=cache_candidate
    local cadd_profile=""
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
                              conservation_file
                              pext_enabled pext_mean_bw pext_tissue_bw_dir pext_tissues)
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
    local loeuf_prescore="${DATA_DIR}/loeuf/$(basename ${loeuf_prescore:-${config_args[loeuf_prescore]}})"
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

    # pext (proportion expressed across transcripts) configuration
    local pext_enabled="${pext_enabled:-${config_args[pext_enabled]:-false}}"
    local pext_mean_bw="${pext_mean_bw:-${config_args[pext_mean_bw]}}"
    local pext_tissue_bw_dir="${pext_tissue_bw_dir:-${config_args[pext_tissue_bw_dir]}}"
    local pext_tissues="${pext_tissues:-${config_args[pext_tissues]}}"

    if [[ "$priva_vep_purpose" == cache_candidate ]]; then
        python "${SCRIPT_DIR}/clinvar_gene_refresh.py" candidate-begin \
            --config "$config_file" --source-root "$BASE_DIR" --input-vcf "$input_vcf" || return 1
    else
        priva_require_legacy_spliceai "${config_file}" "${spliceai_snv_prescore}" "${spliceai_indel_prescore}" || return 1
    fi

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

    # Build pext custom annotation arguments (VEP --custom bigwig)
    # pext = proportion expressed across transcripts from gnomAD (GTEx expression data)
    local pext_arguments=""
    if [[ "$pext_enabled" == "true" ]]; then
        # Add mean pext annotation if file exists
        if [[ -n "$pext_mean_bw" ]] && [[ -f "$pext_mean_bw" ]]; then
            pext_arguments="--custom ${pext_mean_bw},PEXT_MEAN,bigwig,overlap,0"
            log "Adding mean pext annotation from: $pext_mean_bw"
        else
            log "WARNING: pext_enabled=true but pext_mean_bw not found: $pext_mean_bw"
        fi
        
        # Add tissue-specific pext annotations if configured
        # Note: Tissue files are stored in a "tissues" subdirectory under pext_tissue_bw_dir
        if [[ -n "$pext_tissues" ]] && [[ "$pext_tissues" != "null" ]] && [[ -d "$pext_tissue_bw_dir" ]]; then
            local tissue_files_dir="${pext_tissue_bw_dir}/tissues"
            if [[ -d "$tissue_files_dir" ]]; then
                local -a tissues
                mapfile -t tissues < <(get_pext_tissues_array "$pext_tissues")
                
                for tissue in "${tissues[@]}"; do
                    local tissue_bw="${tissue_files_dir}/${tissue}.bw"
                    if [[ -f "$tissue_bw" ]]; then
                        # Convert tissue name to VEP-safe INFO tag (replace special chars with _)
                        local tag="PEXT_$(echo "$tissue" | tr '-' '_')"
                        # Use ${var:+$var } pattern to avoid leading space when pext_arguments is empty
                        pext_arguments="${pext_arguments:+$pext_arguments }--custom ${tissue_bw},${tag},bigwig,overlap,0"
                        log "Adding tissue pext annotation: $tissue -> $tag"
                    else
                        log "WARNING: Tissue pext file not found: $tissue_bw"
                    fi
                done
            else
                log "WARNING: Tissue pext directory not found: $tissue_files_dir"
            fi
        fi
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
    if [[ "$priva_vep_purpose" != cache_candidate ]] && check_vcf_infotags "$input_vcf" CSQ; then
        if [[ -n $cadd_profile ]]; then
            # Old CSQs cannot acquire native ordinals after the fact. Reuse a
            # complete qualified source or fail; never launch a second VEP.
            python "${SCRIPT_DIR}/cadd_shared_workflow.py" verify-source --config "$config_file" \
                --profile "$cadd_profile" --input "$input_vcf" || return 1
            return 0
        fi
        if [[ ${input_vcf} -nt ${SELF_SCRIPT} ]]; then
            log "The input vcf ${input_vcf} is already annotated by VEP. Skip this step"
            return 0
        fi
    fi

    local -a cadd_vep_args=() cadd_vep_environment=()
    local cadd_plugin_dir="" cadd_settings
    if [[ -n $cadd_profile ]]; then
        [[ "$priva_vep_purpose" == runtime ]] || { log "CADD activation is not an installer-cache operation"; return 1; }
        cadd_settings=$(python "${SCRIPT_DIR}/cadd_shared_workflow.py" vep-options \
            --config "$config_file" --profile "$cadd_profile") || return 1
        local -a cadd_values
        mapfile -t cadd_values <<< "$cadd_settings"
        [[ ${#cadd_values[@]} -eq 4 ]] || return 1
        cadd_plugin_dir=" ${cadd_values[0]}"
        cadd_vep_args=(--regulatory --sift b --polyphen b --ccds --cache_version "${cadd_values[2]}"
            --plugin "CaddSelection,cadd_transcript_sources=${cadd_values[1]}")
        cadd_vep_environment=(env "PERL_HASH_SEED=${cadd_values[3]}" PERL_PERTURB_KEYS=0)
    fi

    local tmp_tag=$(randomID)
    local output_vcf=${input_vcf/.vcf*/.${tmp_tag}.vcf}

    log "Running VEP annotation with the command below:"
    log "vep -i ${input_vcf} --format vcf --verbose --vcf --species homo_sapiens --use_transcript_ref --assembly ${assembly} --cache --offline --merged --domains --hgvs --numbers --symbol --canonical --total_length --variant_class --gene_phenotype --stats_file ${input_vcf/.vcf*/.vep.stats.html} --fork ${threads} --buffer_size 10000 --fasta ${ref_genome} --dir_cache ${vep_cache_dir} --dir_plugins '${vep_plugins_dir} ${loftee_repo}' -plugin UTRAnnotator,file=${utr_annotator_file} -plugin LOEUF,file=${loeuf_prescore},match_by=transcript -plugin AlphaMissense,file=${alphamissense_prescore} -plugin LoF,loftee_path:${loftee_repo},human_ancestor_fa:${human_ancestor_fasta}${loftee_conservation_argument}${loftee_gerp_bw_argument} -plugin SpliceAI,snv=${spliceai_snv_prescore},indel=${spliceai_indel_prescore},cutoff=0.5 -plugin PrimateAI,${primateai_prescore} -plugin SpliceVault,file=${splicevault_prescore} -plugin Conservation,${conservation_file},MAX -plugin NMD ${mavedb_argument} ${pext_arguments} --force_overwrite -o ${output_vcf}"
    if [[ -n $cadd_profile ]]; then
        local cadd_logged_options
        printf -v cadd_logged_options '%q ' "${cadd_vep_environment[@]}" vep "${cadd_vep_args[@]}"
        log "Qualified CADD profile ${cadd_profile}; additional plugin directory:${cadd_plugin_dir}; exact added invocation settings: ${cadd_logged_options}"
    fi

    # Run VEP annotation
    "${cadd_vep_environment[@]}" vep -i ${input_vcf} \
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
    --dir_plugins "${vep_plugins_dir} ${loftee_repo}${cadd_plugin_dir}" \
    -plugin UTRAnnotator,file=${utr_annotator_file} \
    -plugin LOEUF,file=${loeuf_prescore},match_by=transcript \
    -plugin AlphaMissense,file=${alphamissense_prescore} \
    -plugin LoF,loftee_path:${loftee_repo},human_ancestor_fa:${human_ancestor_fasta}${loftee_conservation_argument}${loftee_gerp_bw_argument} \
    -plugin SpliceAI,snv=${spliceai_snv_prescore},indel=${spliceai_indel_prescore},cutoff=0.5 \
    -plugin PrimateAI,${primateai_prescore} \
    -plugin SpliceVault,file=${splicevault_prescore} \
    -plugin Conservation,${conservation_file},MAX \
    -plugin NMD ${mavedb_argument} \
    ${pext_arguments} \
    "${cadd_vep_args[@]}" \
    --force_overwrite \
    -o ${output_vcf} && \
    bcftools sort -Oz -o ${input_vcf/.vcf*/.vcf.gz} ${output_vcf} && \
    tabix -f -p vcf ${input_vcf/.vcf*/.vcf.gz} && \
    announce_remove_tmps ${output_vcf} && \
    display_vcf ${input_vcf/.vcf*/.vcf.gz} || return $?

    if [[ -n $cadd_profile ]]; then
        python "${SCRIPT_DIR}/cadd_shared_workflow.py" verify-source --config "$config_file" \
            --profile "$cadd_profile" --input "${input_vcf/.vcf*/.vcf.gz}" || return 1
    fi

    if [[ "$priva_vep_purpose" == cache_candidate ]]; then
        python "${SCRIPT_DIR}/clinvar_gene_refresh.py" candidate-complete \
            --config "$config_file" --source-root "$BASE_DIR" --input-vcf "$input_vcf" || return 1
    fi
}

