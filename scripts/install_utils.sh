#! /usr/bin/env bash
source $(dirname $(realpath ${BASH_SOURCE[0]}))/common_bash_utils.sh



function conda_install_vep() {
    local env_yaml=${1}

    # Extract the env name from the env.yaml file
    local env_name=$(head -1 ${env_yaml} | awk -F ': ' '{print $2;}')
    log "The environment name to-be setup according to the ${env_yaml} file is ${env_name}"
    [[ -z ${env_name} ]] && { log "Failed to extract the env name from the ${env_yaml} file"; return 1; }

    # Test if mamba is available
    if ! command -v mamba &> /dev/null; then
        log "mamba is not installed, please install mamba first"
        return 1
    fi

    # Before do anything, we can check whether the env_name is already installed and check several key packages in the env if already installed
    if [[ $(mamba env list | grep ${env_name}) =~ ${env_name} ]]; then
        log "The environment ${env_name} is already installed"
        local skip_env_creation=1
    else
        log "The environment ${env_name} is not installed"
        local skip_env_creation=0
    fi

    # Test if several key packages are installed
    if [[ ${skip_env_creation} -eq 1 ]]; then
        mamba_activate ${env_name}
        # The to be checked packages are: ensembl-vep, perl-bio-procedural, perl-bioperl
        if [[ $(mamba list | grep ensembl-vep) =~ ensembl-vep ]] && \
           [[ $(mamba list | grep perl-bio-procedural) =~ perl-bio-procedural ]] && \
           [[ $(mamba list | grep perl-bioperl) =~ perl-bioperl ]] && \
           [[ $(which vep) =~ ${CONDA_PREFIX} ]]; then
            log "The key packages are already installed in the environment ${env_name}"
        else
            log "The key packages are not installed in the environment ${env_name}"
            return 1
        fi
    fi

    # You need to separately install the ensembl-vep and perl-bio-procedural packages
    # because Bio::Perl module is migrated from perl-bioperl to perl-bio-procedural since perl-bioperl 1.7.3
    # Or you can directly install the env from our env.yaml file
    if [[ ${skip_env_creation} -eq 0 ]]; then
        mamba env create -f ${env_yaml} && \
        mamba_activate ${env_name} && \
        [[ $(which vep) =~ ${CONDA_PREFIX} ]] && log "The VEP binaries are successfully linked to the conda env ${env_name}" || \
        { log "The VEP binaries are not linked to the conda env ${env_name}, please run the function vep_bin_check manually"; return 1; }
    fi
}


function vep_install_wrapper() {
    local VEP_CACHEDIR=""
    local VEP_ASSEMBLY=""
    local VEP_PLUGINS=""
    local VEP_PLUGINSDIR=""
    local VEP_PLUGINSCACHEDIR=""

    local TEMP
    TEMP=$(getopt -o hc:y:r:p: --long help,VEP_CACHEDIR:,VEP_ASSEMBLY:,VEP_PLUGINSDIR:,VEP_PLUGINSCACHEDIR:,VEP_PLUGINS: -- "$@")

    if [[ $? != 0 ]]; then return 1; fi

    eval set -- "$TEMP"

    while true; do
        case "$1" in
            -h|--help)
                echo "Usage: vep_install [options]"
                echo "Options:"
                echo "  -c, --VEP_CACHEDIR          Set destination directory for cache files"
                echo "  -y, --VEP_ASSEMBLY          Assembly name to use if more than one"
                echo "  -g, --VEP_PLUGINS           Comma-separated list of plugins to install"
                echo "  -r, --VEP_PLUGINSDIR        Set destination directory for VEP plugins files"
                echo "  -p, --VEP_PLUGINSCACHEDIR  Set destination direcotry for VEP plugins caches"
                return 0
                ;;
            -c|--VEP_CACHEDIR)
                VEP_CACHEDIR="$2"
                shift 2
                ;;
            -y|--VEP_ASSEMBLY)
                VEP_ASSEMBLY="$2"
                shift 2
                ;;
            -g|--VEP_PLUGINS)
                VEP_PLUGINS="$2"
                shift 2
                ;;
            -r|--VEP_PLUGINSDIR)
                VEP_PLUGINSDIR="$2"
                shift 2
                ;;
            -p|--VEP_PLUGINSCACHEDIR)
                VEP_PLUGINSCACHEDIR="$2"
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

    # The VEP_DESTDIR is the directory used to access installed PERL dependencies. Usually you need to specify the current conda env perl installation location for site-packages
    [[ $CONDA_PREFIX =~ envs ]] && \
    log "currently in a conda env $(basename $CONDA_PREFIX)" || \
    { log "Not encouraged to install the VEP dependencies directly in the conda base env. So quit with error."; return 1; }

    # Test the PERL5LIB value and PATH value to see if they already include the VEP_DESTDIR and VEP_PLUGINSDIR
    # If yes for both, then we can directly skip the follow up installation of VEP API
    local VEP_DESTDIR=$(perl -e 'print join("\n", @INC);' | head -1)
	[[ $(echo $PERL5LIB) =~ "${VEP_DESTDIR}" ]] && \
    [[ $(echo $PATH) =~ "${VEP_DESTDIR}" ]] && \
    [[ $(echo $PERL5LIB}) =~ "${VEP_PLUGINSDIR}" ]] && \
    { log "The PERL5LIB value and PATH value already include ${VEP_DESTDIR} and ${VEP_PLUGINSDIR}, indicating the following installation process has already been succesfully performed. Skip the function for now."; return 0; }

    local conda_env_name=$(basename $CONDA_PREFIX)
    [[ ${VEP_DESTDIR} =~ ${conda_env_name} ]] && log "The dest dir is set to the directory (${VEP_DESTDIR}) where perl modules are installed by conda" || \
    { log "Since the function is designed to perform follow up installation of VEP upon the installation of VEP dependencies via conda, we only accept installing VEP API at the perl module installation location previously used by conda"; return 1; }

    # Construct the command
    [[ ${VEP_ASSEMBLY} =~ "hg19" ]] && VEP_ASSEMBLY="GRCh37"
    [[ ${VEP_ASSEMBLY} =~ "hg38" ]] && VEP_ASSEMBLY="GRCh38"
    local cmd="vep_install -d ${VEP_DESTDIR} --AUTO acp -s homo_sapiens_merged --NO_HTSLIB --NO_BIOPERL --CONVERT"
    [[ -n "$VEP_CACHEDIR" ]] && cmd+=" --CACHEDIR $VEP_CACHEDIR"
    [[ -n "$VEP_ASSEMBLY" ]] && cmd+=" --ASSEMBLY $VEP_ASSEMBLY"
    [[ -n "$VEP_PLUGINS" ]] && [[ $VEP_PLUGINS != "empty" ]] && cmd+=" --PLUGINS $VEP_PLUGINS,SpliceAI,UTRAnnotator,AlphaMissense,GO,Phenotypes,Conservation,LOEUF,SameCodon,NMD,PrimateAI" || cmd+=" --PLUGINS SpliceAI,UTRAnnotator,AlphaMissense,GO,Phenotypes,Conservation,LOEUF,SameCodon,NMD,PrimateAI"
    [[ -n "$VEP_PLUGINSDIR" ]] && cmd+=" --PLUGINSDIR $VEP_PLUGINSDIR"

    # Execute the command
    # After executing the command, we need to bind the new PERL5LIB value with the current conda env
    log "Now we start to install the VEP api and downloading caches (which might take a while to finish). So pls be patient, here is the command we are going to execute: ${cmd}"
    $cmd && \
    conda env config vars set PERL5LIB="$VEP_DESTDIR:$VEP_PLUGINSDIR" && \
    conda env config vars set PATH="$VEP_DESTDIR:$VEP_DESTDIR/htslib:$PATH" && \
    log "Now we have bound the new PERL5LIB value (adding ${VEP_DESTDIR} and ${VEP_PLUGINSDIR} to the PERL5LIB) with the current conda env ${CONDA_ENV_NAME}" && \
    log "Now we have bound the new PATH value (adding ${VEP_DESTDIR} to the PATH) with the current conda env ${CONDA_ENV_NAME}" && \
    mamba_deactivate && \
    mamba_activate ${conda_env_name} || \
    { log "Failed to install the VEP api and download caches"; return 1; }
}


function VEP_plugins_install() {
    local VEP_PLUGINSDIR=${1}
    local VEP_PLUGINSCACHEDIR=${2}
    local ASSEMBLY_VERSION=${3}
    local config_file=${4}
    local conda_env_name=${5}


    [[ $CONDA_PREFIX =~ "envs" ]] && \
    log "currently in a conda env $(basename $CONDA_PREFIX)" || \
    { [[ -z ${conda_env_name} ]] && log "Not encouraged to install the VEP plugins directly in the conda base env. So quit with error." && return 1; }

    if [[ -z ${conda_env_name} ]]; then
        local conda_env_name=$(basename $CONDA_PREFIX)
    fi

    # Install VEP plugins
    # First enter the VEP plugins directory
    cd $VEP_PLUGINSDIR

    # Then install UTRAnnotator
    local utr_annotator_files=$(UTRAnnotator_install ${VEP_PLUGINSCACHEDIR}) || \
    { log "Failed to install UTRAnnotator"; return 1; }

    if [[ ${ASSEMBLY_VERSION} =~ "GRCh37" ]] || [[ ${ASSEMBLY_VERSION} =~ "hg19" ]]; then
        update_yaml "$config_file" "utr_annotator_file" "$(echo ${utr_annotator_files} | head -1)"
    elif [[ ${ASSEMBLY_VERSION} =~ "GRCh38" ]] || [[ ${ASSEMBLY_VERSION} =~ "hg38" ]]; then
        update_yaml "$config_file" "utr_annotator_file" "$(echo ${utr_annotator_files} | tail -1)"
    else
        log "Not supported assembly version: ${ASSEMBLY_VERSION}"
        return 1
    fi


    # Then install LOEUF
    local loeuf_prescore=$(LOEUF_install ${VEP_PLUGINSCACHEDIR} ${ASSEMBLY_VERSION} ${VEP_PLUGINSDIR}) || \
    { log "Failed to install LOEUF"; return 1; }
    update_yaml "$config_file" "loeuf_prescore" "${loeuf_prescore}"


    # Then install AlphaMissense
    local alphamissense_prescore=$(AlphaMissense_install ${VEP_PLUGINSCACHEDIR} ${ASSEMBLY_VERSION}) || \
    { log "Failed to install AlphaMissense"; return 1; }
    update_yaml "$config_file" "alphamissense_prescore" "${alphamissense_prescore}"


    # Then install SpliceAI
    local spliceai_prescores=$(SpliceAI_install ${VEP_PLUGINSCACHEDIR} ${VEP_PLUGINSDIR}) || \
    { log "Failed to install SpliceAI"; return 1; }
    update_yaml "$config_file" "spliceai_snv_prescore" "$(echo ${spliceai_prescores} | head -1)" && \
    update_yaml "$config_file" "spliceai_indel_prescore" "$(echo ${spliceai_prescores} | tail -1)"


    # Last, install PrimateAI
    local primateai_prescore=$(PrimateAI_install ${VEP_PLUGINSCACHEDIR} ${VEP_PLUGINSDIR}) || \
    { log "Failed to install PrimateAI"; return 1; }
    update_yaml "$config_file" "primateai_prescore" "${primateai_prescore}"


    # Install Conservation file if assembly version is hg38
    if [[ ${ASSEMBLY_VERSION} =~ "GRCh38" ]] || [[ ${ASSEMBLY_VERSION} =~ "hg38" ]]; then
        log "Now we start to install the conservation file for hg38"
        local conservation_file=$(Conservation_install ${VEP_PLUGINSCACHEDIR} ${ASSEMBLY_VERSION}) || \
        { log "Failed to install Conservation file"; return 1; }
        update_yaml "$config_file" "conservation_file" "${conservation_file}"
    else
        log "Not supported assembly version: ${ASSEMBLY_VERSION}, so skip installing the conservation file"
    fi
}


function UTRAnnotator_install() {
    local PLUGIN_CACHEDIR=${1}

    if [[ ! -d ${PLUGIN_CACHEDIR}/UTRannotator ]] || \
       [[ ! -f ${PLUGIN_CACHEDIR}/UTRannotator/uORF_starts_ends_GRCh38_PUBLIC.txt ]] || \
       [[ ! -f ${PLUGIN_CACHEDIR}/UTRannotator/uORF_starts_ends_GRCh37_PUBLIC.txt ]]; then
        cd ${PLUGIN_CACHEDIR} && \
        git clone https://github.com/Ensembl/UTRannotator && \
        log "${PLUGIN_CACHEDIR}/UTRannotator/uORF_5UTR_GRCh37_PUBLIC.txt and ${PLUGIN_CACHEDIR}/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt are well-downloaded, remember to specify them when using UTRAnnotator"
        log "A typical command syntax would be: vep -i variations.vcf --plugin UTRAnnotator,file=${PLUGIN_CACHEDIR}/UTRannotator/uORF_starts_ends_GRCh38_PUBLIC.txt"
    else
        log "The UTRannotator plugin is already downloaded for both hg19 and hg38 assemblies"
    fi

    echo ${PLUGIN_CACHEDIR}/UTRannotator/uORF_5UTR_GRCh37_PUBLIC.txt && \
    echo ${PLUGIN_CACHEDIR}/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt
}


function SpliceAI_install() {
    local PLUGIN_CACHEDIR=${1}
    local PLUGIN_DIR=${2}

    log "You need to download the prescores yourself to ${PLUGIN_CACHEDIR}/SpliceAI by registering a free BaseSequence Account from Illumina. Note that the downloaded files are pretty large (around 200 GB) and might take hours to days to complete the downloading process."
    log "For details, read ${PLUGIN_DIR}/SpliceAI.pm"
    log "The following steps are necessary before running this plugin:

        The files with the annotations for all possible substitutions (snv), 1 base insertions
        and 1-4 base deletions (indel) within genes are available here:
        https://basespace.illumina.com/s/otSPW8hnhaZR

        GRCh37:
        tabix -p vcf spliceai_scores.raw.snv.hg37.vcf.gz
        tabix -p vcf spliceai_scores.raw.indel.hg37.vcf.gz

        GRCh38:
        tabix -p vcf spliceai_scores.raw.snv.hg38.vcf.gz
        tabix -p vcf spliceai_scores.raw.indel.hg38.vcf.gz

        The plugin can then be run:
        ./vep -i variations.vcf --plugin SpliceAI,snv=/path/to/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/path/to/spliceai_scores.raw.indel.hg38.vcf.gz
        ./vep -i variations.vcf --plugin SpliceAI,snv=/path/to/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/path/to/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 "
    # Prepare a command waiting for the user to download the prescores and put them into the corresponding directory
    # Waiting for user's confirmation (yes or no) to finish the installation properly, we dont need to offer the cmd directly because such downloading will need them to register account and login, which is not convenient to put in the script

    read -p "Have you finished downloading the prescores for the SpliceAI plugin? (yes or no)"
    if [[ ${REPLY} =~ "yes" ]] || [[ ${REPLY} =~ "y" ]] || [[ ${REPLY} =~ "Y" ]] || [[ ${REPLY} =~ "Yes" ]] || [[ ${REPLY} =~ "YES" ]]; then
        log "Now the installation of the SpliceAI plugin is completed. Please specify the paths to the prescore files (one for snv and one for indel)"
        read -p "Please specify the absolute path to the snv prescore file, remember the file should correspond to the assembly version of your VEP installation"
        local snv_pre_score_file=${REPLY}
        if [[ ! -f ${snv_pre_score_file} ]]; then
            log "The file ${snv_pre_score_file} does not exist"
            return 1
        fi
        read -p "Please specify the absolute path to the indel prescore file, remember the file should correspond to the assembly version of your VEP installation"
        local indel_pre_score_file=${REPLY}
        if [[ ! -f ${indel_pre_score_file} ]]; then
            log "The file ${indel_pre_score_file} does not exist"
            return 1
        fi
        log "Now we start to install the SpliceAI plugin"
        # We need to return the prescore file paths to the outside of the function
        echo ${snv_pre_score_file}
        echo ${indel_pre_score_file}
    else
        log "Please download the prescores for the SpliceAI plugin first"
        return 1
    fi
}



function convert_splicevault_hg38_to_hg19() {
    local input_tsv=$1
    local output_tsv=$2
    local chain_file=$3
    local reference_fasta=$4  # GRCh38 reference
    local temp_dir=$(mktemp -d)

    log "Starting conversion process..."

    # 1. Convert TSV to VCF with all columns as INFO fields
    log "Converting TSV to VCF..."
    bcftools convert --tsv2vcf "$input_tsv" \
        -f "$reference_fasta" \
        -s "SAMPLE" \
        --columns CHROM,POS,REF,ALT \
        --annot-fields TRANSCRIPT_ID,SITE_TYPE,SITE_LOC,SPLICEAI_DELTA,OUT_OF_FRAME,TOP4_EVENTS,SAMPLE_COUNT,MAX_DEPTH \
        | bcftools annotate \
            --set-id '%CHROM\_%POS\_%REF\_%ALT' \
            -Oz -o "${temp_dir}/grch38.vcf.gz" && \
    display_vcf "${temp_dir}/grch38.vcf.gz"

    # Add INFO field descriptions to header
    bcftools header -h <(cat << EOF
##INFO=<ID=TRANSCRIPT_ID,Number=1,Type=String,Description="Ensembl transcript ID">
##INFO=<ID=SITE_TYPE,Number=1,Type=String,Description="SpliceVault site type">
##INFO=<ID=SITE_LOC,Number=1,Type=String,Description="SpliceVault site location">
##INFO=<ID=SPLICEAI_DELTA,Number=1,Type=Float,Description="SpliceAI delta score">
##INFO=<ID=OUT_OF_FRAME,Number=1,Type=String,Description="Out of frame events">
##INFO=<ID=TOP4_EVENTS,Number=1,Type=String,Description="Top 4 splicing events">
##INFO=<ID=SAMPLE_COUNT,Number=1,Type=Integer,Description="Sample count">
##INFO=<ID=MAX_DEPTH,Number=1,Type=Integer,Description="Maximum depth">
EOF
    ) -o "${temp_dir}/grch38.header.vcf.gz"

    # 2. Liftover using CrossMap
    log "Lifting over coordinates..."
    CrossMap.py vcf "$chain_file" \
        "${temp_dir}/grch38.vcf.gz" \
        "$reference_fasta" \
        "${temp_dir}/hg19.vcf" && \
    display_vcf "${temp_dir}/hg19.vcf"

    # Compress and index the lifted VCF
    bgzip -f "${temp_dir}/hg19.vcf" && \
    bcftools index -t "${temp_dir}/hg19.vcf.gz" && \
    display_vcf "${temp_dir}/hg19.vcf.gz"

    # 3. Convert back to TSV format preserving all INFO fields
    log "Converting back to TSV format..."
    # First, save header
    head -n 1 "$input_tsv" > "$output_tsv"

    # Extract data from VCF including all INFO fields
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/TRANSCRIPT_ID\t%INFO/SITE_TYPE\t%INFO/SITE_LOC\t%INFO/SPLICEAI_DELTA\t%INFO/OUT_OF_FRAME\t%INFO/TOP4_EVENTS\t%INFO/SAMPLE_COUNT\t%INFO/MAX_DEPTH\n' \
        "${temp_dir}/hg19.vcf.gz" >> "$output_tsv"

    # Clean up
    rm -rf "$temp_dir"
    log "Conversion complete. Output written to $output_tsv"
}



function CADD_install() {
    local config_file=${1}
    local assembly=${2}

    [[ -z ${assembly} ]] && assembly=$(read_yaml "${config_file}" "assembly")
    [[ ${assembly} == "hg19" ]] && assembly="GRCh37"
    [[ ${assembly} == "hg38" ]] && assembly="GRCh38"
    local CADD_script

    # Running CADD requires 4 big parts of preparations:
    # 1. snakemake (solved by conda environment)
    # 2. CADD dependencies (solved by conda environment)
    # 3. CADD genome annotation file (need to download)
    # 4. CADD prescores (need to download)

    CADD_script=$(find ${CONDA_PREFIX}/ -type f -name "CADD.sh")
    if [[ -z ${CADD_script} ]] || [[ ! -f ${CADD_script} ]]; then
        log "Cannot find the CADD script in the conda environment"
        return 1
    fi

    local CADD_cache_dir=$(dirname ${CADD_script})/data
    # Note that by default, CADD will use the CADD cache file in the conda environment
    # Ask the user whether they want to use the default CADD cache directory to store the genome annotation file and the prescores.
    read -p "Do you want to use the default CADD cache directory ${CADD_cache_dir} to store the genome annotation file and the prescores which could be take up to 1TB space? (yes or no)"
    if [[ ${REPLY} =~ "yes" ]] || [[ ${REPLY} =~ "y" ]] || [[ ${REPLY} =~ "Y" ]] || [[ ${REPLY} =~ "Yes" ]] || [[ ${REPLY} =~ "YES" ]]; then
        local CADD_prescore_dir=${CADD_cache_dir}/prescored
        local CADD_anno_dir=${CADD_cache_dir}/annotations
    else
        read -p "In this case, you need to download the CADD repo as a zip file and unzip it to a local directory and all the cache files will be store in that directory. Please specify the absolute path to the directory: "
        local CADD_parent_dir=${REPLY}
        local CADD_zip_download_url=$(read_yaml "${config_file}" "cadd_zip_download_url")
        wget ${CADD_zip_download_url} -O ${CADD_parent_dir}/CADD-scripts.zip && \
        unzip ${CADD_parent_dir}/CADD-scripts.zip -d ${CADD_parent_dir}/ && \
        rm ${CADD_parent_dir}/CADD-scripts.zip && \
        # Get the base folder name from the unzipped directory
        local cadd_version=$(read_yaml "${config_file}" "cadd_version")
        local CADD_base_dir=$(find ${CADD_parent_dir}/ -maxdepth 1 -type d -name "CADD-scripts-*${cadd_version#v}" -print | head -n1) && \
        [[ -z ${CADD_base_dir} ]] && { log "Could not find CADD scripts directory"; return 1; }
        local CADD_script=${CADD_base_dir}/CADD.sh && \
        local CADD_cache_dir=${CADD_base_dir}/data && \
        local CADD_prescore_dir=${CADD_cache_dir}/prescored && \
        local CADD_anno_dir=${CADD_cache_dir}/annotations
    fi


    # Now we start downloading the CADD genome annotations corresponding to the assembly version
    if [[ ${assembly} == "GRCh37" ]] || [[ ${assembly} == "hg19" ]]; then
        wget -c https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh37/GRCh37_v1.7.tar.gz -O ${CADD_anno_dir}/GRCh37_v1.7.tar.gz && \
        md5sum ${CADD_anno_dir}/GRCh37_v1.7.tar.gz | grep -q $(read_yaml "${config_file}" "cadd_GRCh37_anno_md5") && \
        tar -xzvf ${CADD_anno_dir}/GRCh37_v1.7.tar.gz -C ${CADD_anno_dir}/ && \
        rm ${CADD_anno_dir}/GRCh37_v1.7.tar.gz
    elif [[ ${assembly} == "GRCh38" ]] || [[ ${assembly} == "hg38" ]]; then
        wget -c https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/GRCh38_v1.7.tar.gz -O ${CADD_anno_dir}/GRCh38_v1.7.tar.gz && \
        md5sum ${CADD_anno_dir}/GRCh38_v1.7.tar.gz | grep -q $(read_yaml "${config_file}" "cadd_GRCh38_anno_md5") && \
        tar -xzvf ${CADD_anno_dir}/GRCh38_v1.7.tar.gz -C ${CADD_anno_dir}/ && \
        rm ${CADD_anno_dir}/GRCh38_v1.7.tar.gz
    else
        log "Not supported assembly version: ${assembly}"
        return 1
    fi


    # Now we start downloading the CADD prescores corresponding to the assembly version
    local snv_prescore_url
    local indel_prescore_url
    local snv_prescore_md5
    local indel_prescore_md5
    if [[ ${assembly} == "GRCh37" ]] || [[ ${assembly} == "hg19" ]]; then
        snv_prescore_url=$(read_yaml "${config_file}" "cadd_GRCh37_snv_anno_url")
        indel_prescore_url=$(read_yaml "${config_file}" "cadd_GRCh37_indel_anno_url")
        snv_prescore_md5=$(read_yaml "${config_file}" "cadd_GRCh37_snv_anno_md5")
        indel_prescore_md5=$(read_yaml "${config_file}" "cadd_GRCh37_indel_anno_md5")
    elif [[ ${assembly} == "GRCh38" ]] || [[ ${assembly} == "hg38" ]]; then
        snv_prescore_url=$(read_yaml "${config_file}" "cadd_GRCh38_snv_anno_url")
        indel_prescore_url=$(read_yaml "${config_file}" "cadd_GRCh38_indel_anno_url")
        snv_prescore_md5=$(read_yaml "${config_file}" "cadd_GRCh38_snv_anno_md5")
        indel_prescore_md5=$(read_yaml "${config_file}" "cadd_GRCh38_indel_anno_md5")
    fi

    local snv_file_name=$(basename ${snv_prescore_url})
    local indel_file_name=$(basename ${indel_prescore_url})
    wget -c ${snv_prescore_url} -O ${CADD_prescore_dir}/${snv_file_name} && \
    md5sum ${CADD_prescore_dir}/${snv_file_name} | grep -q ${snv_prescore_md5} && \
    wget -c ${snv_prescore_url}.tbi -O ${CADD_prescore_dir}/${snv_file_name}.tbi && \
    wget -c ${indel_prescore_url} -O ${CADD_prescore_dir}/${indel_file_name} && \
    md5sum ${CADD_prescore_dir}/${indel_file_name} | grep -q ${indel_prescore_md5} && \
    wget -c ${indel_prescore_url}.tbi -O ${CADD_prescore_dir}/${indel_file_name}.tbi
}



function PrimateAI_install() {
    local PLUGIN_CACHEDIR=${1}
    local PLUGIN_DIR=${2}

    log "You need to download the prescores yourself to ${PLUGIN_CACHEDIR}/PrimateAI by registering a free BaseSequence Account from Illumina. Note that the downloaded files are a bit large (around 2 GB) and might take minutes to complete the downloading process."
    log "For details, read ${PLUGIN_DIR}/PrimateAI.pm"
    log "In brief, common missense mutations in non-human primate species are usually
        benign in humans. Thousands of common variants from six non-human primate
        species were used to train a deep neural network to identify pathogenic
        mutations in humans with a rare disease.

        This plugin uses files generated by the PrimateAI software, which is
        available from https://github.com/Illumina/PrimateAI. The files containing
        predicted pathogenicity scores can be downloaded from
        https://basespace.illumina.com/s/yYGFdGih1rXL (a free BaseSpace account may
        be required):
            PrimateAI_scores_v0.2.tsv.gz (for GRCh37/hg19)
            PrimateAI_scores_v0.2_hg38.tsv.gz (for GRCh38/hg38)

        Before running the plugin for the first time, the following steps must be
        taken to format the downloaded files:

        1.  Unzip the score files
        2.  Add '#' in front of the column description line
        3.  Remove any empty lines.
        4.  Sort the file by chromosome and position
        5.  Compress the file in .bgz format
        6.  Create tabix index (requires tabix to be installed).

        Command line examples for formatting input files:
            > gunzip -cf PrimateAI_scores_v0.2.tsv.gz | sed '12s/.*/#&/' | sed '/^$/d' | awk 'NR<12{print $0;next}{print $0 | "sort -k1,1 -k 2,2n -V"}' | bgzip > PrimateAI_scores_v0.2_GRCh37_sorted.tsv.bgz
            > tabix -s 1 -b 2 -e 2 PrimateAI_scores_v0.2_GRCh37_sorted.tsv.bgz

            > gunzip -cf PrimateAI_scores_v0.2_hg38.tsv.gz | sed '12s/.*/#&/' | sed '/^$/d' | awk 'NR<12{print $0;next}{print $0 | "sort -k1,1 -k 2,2n -V"}' | bgzip > PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz
            > tabix -s 1 -b 2 -e 2 PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz

        Example running command:
          ./vep -i variations.vcf --plugin PrimateAI,PrimateAI_scores_v0.2_GRCh37_sorted.tsv.bgz
          ./vep -i variations.vcf --plugin PrimateAI,PrimateAI_scores_v0.2_GRCh38_sorted.tsv.bgz"

    read -p "Have you finished downloading the files? (yes or no)"
    if [[ ${REPLY} =~ "yes" ]] || [[ ${REPLY} =~ "y" ]] || [[ ${REPLY} =~ "Y" ]] || [[ ${REPLY} =~ "Yes" ]] || [[ ${REPLY} =~ "YES" ]]; then
        read -p "Please specify the absolute path to the tsv prescore file, choose the right assembly version"
        local prescore_file=${REPLY}
        if [[ -f ${prescore_file} ]] && [[ ${prescore_file} =~ \.tsv\.gz$ ]]; then
            # We need to return the prescore file paths to the outside of the function
			gunzip -cf ${prescore_file} | \
			sed '12s/.*/#&/' | \
			sed '/^$/d' | \
			awk 'NR<12{print $0;next}{print $0 | "sort -k1,1 -k 2,2n -V"}' | \
			bgzip > ${prescore_file/.tsv.gz/.sorted.tsv.bgz} && \
			tabix -s 1 -b 2 -e 2 ${prescore_file/.tsv.gz/.sorted.tsv.bgz} && \
            echo ${prescore_file/.tsv.gz/.sorted.tsv.bgz}
		elif [[ -f ${prescore_file} ]] && [[ ${prescore_file} =~ \.tsv.bgz$ ]]; then
			echo ${prescore_file}
        else
            log "The file ${prescore_file} does not exist"
            return 1
        fi
    else
        log "Please finish formatting the downloaded files first"
        return 1
    fi
}


function AlphaMissense_install() {
    local PLUGIN_CACHEDIR=${1}
    local ASSEMBLY_VERSION=${2}

    log "For details, read ${PLUGIN_DIR}/AlphaMissense.pm"
    local hg19_url="https://console.cloud.google.com/storage/browser/_details/dm_alphamissense/AlphaMissense_hg19.tsv.gz"
    local hg38_url="https://storage.cloud.google.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"


    if [[ ! -d ${PLUGIN_CACHEDIR} ]]; then
        mkdir -p ${PLUGIN_CACHEDIR}
    fi

    if [[ ! -d ${PLUGIN_CACHEDIR}/AlphaMissense ]]; then
        mkdir -p ${PLUGIN_CACHEDIR}/AlphaMissense
    fi

    if [[ -f ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg19.tsv.gz ]] && [[ ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg19.tsv.gz -ot ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg19.tsv.gz.tbi ]]; then
        log "The AlphaMissense_hg19.tsv.gz is already downloaded to ${PLUGIN_CACHEDIR}"
    else
        wget ${hg19_url} -O ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg19.tsv.gz && \
        tabix -s 1 -b 2 -e 2 -f -S 1 ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg19.tsv.gz || \
        { log "Failed to index the AlphaMissense_hg19.tsv.gz file, this might be due to the downloading file format not correct. Now we need you to download the file manually by opening the url ${hg19_url} in your browser and save the file to ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg19.tsv.gz"; \
          read -p "Have you finished downloading the file? (yes or no)"; \
          if [[ ${REPLY} =~ "yes" ]] || [[ ${REPLY} =~ "y" ]] || [[ ${REPLY} =~ "Y" ]] || [[ ${REPLY} =~ "Yes" ]] || [[ ${REPLY} =~ "YES" ]]; then
            log "Now we start to index the file"
            tabix -s 1 -b 2 -e 2 -f -S 1 ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg19.tsv.gz
          else
            log "Please download the file first"
            return 1
          fi
        }
    fi

    if [[ -f ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg38.tsv.gz ]] && [[ ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg38.tsv.gz -ot ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg38.tsv.gz.tbi ]]; then
        log "The AlphaMissense_hg38.tsv.gz is already downloaded to ${PLUGIN_CACHEDIR}"
    else
        wget ${hg38_url} -O ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg38.tsv.gz && \
        tabix -s 1 -b 2 -e 2 -f -S 1 ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg38.tsv.gz || \
        { log "Failed to index the AlphaMissense_hg38.tsv.gz file, this might be due to the downloading file format not correct. Now we need you to download the file manually by opening the url ${hg38_url} in your browser and save the file to ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg38.tsv.gz"; \
          read -p "Have you finished downloading the file? (yes or no)"; \
          if [[ ${REPLY} =~ "yes" ]] || [[ ${REPLY} =~ "y" ]] || [[ ${REPLY} =~ "Y" ]] || [[ ${REPLY} =~ "Yes" ]] || [[ ${REPLY} =~ "YES" ]]; then
            log "Now we start to index the file"
            tabix -s 1 -b 2 -e 2 -f -S 1 ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg38.tsv.gz
          else
            log "Please download the file first"
            return 1
          fi
        }
    fi

    if [[ ${ASSEMBLY_VERSION} == "GRCh37" ]] || [[ ${ASSEMBLY_VERSION} == "hg19" ]]; then
        echo ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg19.tsv.gz
    elif [[ ${ASSEMBLY_VERSION} == "GRCh38" ]] || [[ ${ASSEMBLY_VERSION} == "hg38" ]]; then
        echo ${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg38.tsv.gz
    else
        log "Not supported assembly version: ${ASSEMBLY_VERSION}"
        return 1
    fi
}


function LOEUF_install() {
    local PLUGIN_CACHEDIR=${1}
    local assembly_version=${2}
	local PLUGIN_DIR=${3}

    if [[ -z ${assembly_version} ]]; then
        local assembly_version="hg19"
    fi

    if [[ ! -d ${PLUGIN_CACHEDIR} ]]; then
        mkdir -p ${PLUGIN_CACHEDIR}
    fi

    if [[ ! -d ${PLUGIN_CACHEDIR}/LOEUF ]]; then
        mkdir -p ${PLUGIN_CACHEDIR}/LOEUF
    fi

    log "Starting to install the annotation files for plugin LOEUF, for details, read ${PLUGIN_DIR}/LOEUF.pm"

    if [[ ! -d ${PLUGIN_CACHEDIR}/supplement ]]; then
        local total_package_url="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7334197/bin/41586_2020_2308_MOESM4_ESM.zip"
        wget ${total_package_url} -O ${PLUGIN_CACHEDIR}/41586_2020_2308_MOESM4_ESM.zip && \
        unzip ${PLUGIN_CACHEDIR}/41586_2020_2308_MOESM4_ESM.zip -d ${PLUGIN_CACHEDIR} || \
        { log "Failed to download and unzip the LOEUF files"; return 1; }
    else
        log "The LOEUF files are already downloaded to ${PLUGIN_CACHEDIR}/supplement"
    fi

    if [[ ${assembly_version} == "GRCh37" ]] || [[ ${assembly_version} == "hg19" ]]; then
        if [[ ! -f ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset.tsv.gz ]] || \
        [[ ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset.tsv.gz -ot ${PLUGIN_CACHEDIR}/supplement/supplementary_dataset_11_full_constraint_metrics.tsv.gz ]] || \
        [[ ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset.tsv.gz.tbi -ot ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset.tsv.gz ]]; then
            zcat ${PLUGIN_CACHEDIR}/supplement/supplementary_dataset_11_full_constraint_metrics.tsv.gz | (head -n 1 && tail -n +2  | sort -t$'\t' -k 76,76 -k 77,77n ) > ${PLUGIN_CACHEDIR}/supplement/loeuf_temp.tsv && \
            sed '1s/.*/#&/' ${PLUGIN_CACHEDIR}/supplement/loeuf_temp.tsv > ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset.tsv && \
            bgzip ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset.tsv && \
            tabix -f -s 76 -b 77 -e 78 ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset.tsv.gz || \
			{ log "Failed to index the LOEUF file"; return 1; }
        else
            log "The LOEUF files are already processed and indexed for hg19 assembly"
        fi
        echo ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset.tsv.gz
    fi

    if [[ ${assembly_version} == "GRCh38" ]] || [[ ${assembly_version} == "hg38" ]]; then
        log "Since the LOEUF scores are primarily based on the hg19 assembly, we need to re-process the LOEUF files for hg38 assembly"
        if [[ ! -f ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset_grch38.tsv.gz ]] || \
        [[ ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset_grch38.tsv.gz -nt ${PLUGIN_CACHEDIR}/supplement/supplementary_dataset_11_full_constraint_metrics_grch38.tsv ]] || \
        [[ ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset_grch38.tsv.gz.tbi -ot ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset_grch38.tsv.gz ]]; then
            # Extract the VEP version from VEP binary executable first
            local VEP_VERSION=$(vep_install --help | grep ensembl-vep | head -1 | awk '{print $NF;}' | awk -F '.' '{print $1;}')
            wget https://ftp.ensembl.org/pub/release-${VEP_VERSION}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${VEP_VERSION}.chr.gtf.gz -O ${PLUGIN_CACHEDIR}/supplement/Homo_sapiens.GRCh38.${VEP_VERSION}.chr.gtf.gz
            # Find the path of the current shell script, pwd is not working
            local SCRIPT_DIR=$(dirname $(readlink -f $0))
            python3 ${SCRIPT_DIR}/python_utils.py -f extract_transcript_coordinates -a "${PLUGIN_CACHEDIR}/supplement/Homo_sapiens.GRCh38.${VEP_VERSION}.chr.gtf.gz;${PLUGIN_CACHEDIR}/supplement/Homo_sapiens.GRCh38.${VEP_VERSION}.chr.tranx.span.tsv" && \
            python3 ${SCRIPT_DIR}/python_utils.py -f update_tranx_span -a "${PLUGIN_CACHEDIR}/supplement/supplementary_dataset_11_full_constraint_metrics.tsv.gz;${PLUGIN_CACHEDIR}/supplement/Homo_sapiens.GRCh38.${VEP_VERSION}.chr.tranx.span.tsv;${PLUGIN_CACHEDIR}/supplement/supplementary_dataset_11_full_constraint_metrics.grch38.tsv.gz" && \
            zcat ${PLUGIN_CACHEDIR}/supplement/supplementary_dataset_11_full_constraint_metrics_grch38.tsv | (head -n 1 && tail -n +2  | sort -t$'\t' -k 76,76 -k 77,77n ) > ${PLUGIN_CACHEDIR}/supplement/loeuf_grch38_temp.tsv && \
            sed '1s/.*/#&/' ${PLUGIN_CACHEDIR}/supplement/loeuf_grch38_temp.tsv > ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset_grch38.tsv && \
            bgzip ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset_grch38.tsv && \
            tabix -f -s 76 -b 77 -e 78 ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset_grch38.tsv.gz
        else
            log "The LOEUF files are already processed and indexed for hg38 assembly"
        fi
        echo ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset_grch38.tsv.gz
    fi
}


function gnomAD_install() {
    local CACHEDIR=${1}
    local assembly=${2}

    # Assembly is the fasta file path
    # If assembly is not provided, we will use GRCh38 assembly

    if [[ -z ${assembly} ]]; then
        local assembly="GRCh38"
    fi

    # We need to iterate over all chromosomes to check if the files are already downloaded and valid
    # local -a chromosomes=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
    local -a chromosomes=( "22" "X" "Y" )
    for chr in "${chromosomes[@]}"; do
        log "About to download the gnomAD v4.1 files for chr${chr} from the url https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz"
        if check_vcf_validity ${CACHEDIR}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz; then
            log "The gnomAD v4.1 file for chr${chr} is already downloaded to ${CACHEDIR} and updated"
            if [[ ${CACHEDIR}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz.tbi -ot ${CACHEDIR}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz ]]; then
                log "The index file is not updated, now we start to index the file"
                tabix -f -p vcf ${CACHEDIR}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz
            fi
        else
            log "The gnomAD v4.1 file for chr${chr} (${CACHEDIR}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz) is not downloaded or not updated, now we start to download the files"
            wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz -O ${CACHEDIR}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz && \
            wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz.tbi -O ${CACHEDIR}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz.tbi
        fi
    done

    log "The gnomAD v4.1 files are already downloaded to ${CACHEDIR}, remember that the VCF records are 1-based and currently mapped to hg38 assembly"

    if [[ ${assembly} =~ "GRCh37" ]] || [[ ${assembly} =~ "hg19" ]]; then
        local self_script=$(realpath ${BASH_SOURCE[0]})
        local self_dir=$(dirname ${self_script})
        local chain_file=$(dirname ${self_dir})/data/liftover/hg38ToHg19.over.chain.gz
        log "Since the gnomAD v4.1 files are currently mapped to hg38 assembly, we need to liftover the files to hg19 assembly, the chain file is ${chain_file}, self_dir is ${self_dir}, self_script is ${self_script}"

        gnomAD_liftover \
        ${CACHEDIR}/gnomad.joint.v4.1.sites.chr1.vcf.bgz \
        ${chain_file} \
        ${assembly} 4 && \
        log "The gnomAD v4.1 files are liftovered to hg19 assembly and saved to ${CACHEDIR}" && \
        echo ${CACHEDIR}/gnomad.joint.v4.1.sites.hg19.chrX.vcf.gz || \
        { log "Failed to liftover the gnomAD v4.1 files to hg19 assembly"; return 1; }
    else
        echo ${CACHEDIR}/gnomad.joint.v4.1.sites.chrX.vcf.bgz
    fi
}


function gnomAD_liftover_per_chromosome() {
    local hg38_vcf=${1}
    local output_dir=${2}
    local hg19_fasta=${3}
    local chain_file=${4}
    local hg19_vcf_name=$(basename ${hg38_vcf/.vcf*/.hg19.vcf.gz}) # make sure the suffix is .vcf.gz instead of .vcf or .vcf.bgz

    local chrom=$(basename ${hg38_vcf/.vcf*/} | awk -F '.' '{print $NF;}')
    log "About to liftover the gnomAD v4.1 VCF file ${hg38_vcf} on chromosome ${chrom} to hg19 assembly, the chain file is ${chain_file}"

    local hg19_vcf=${output_dir}/${hg19_vcf_name}

    if check_vcf_validity ${hg19_vcf}; then
        if [[ ${hg19_vcf}.tbi -nt ${hg19_vcf} ]]; then
            log "The hg19 VCF file is already indexed and the VCF file is valid and updated"
            return 0
        else
            tabix -f -p vcf ${hg19_vcf} && \
            log "The hg19 VCF file is already indexed and the VCF file is valid and updated" && \
            return 0 || \
            { log "Failed to index the hg19 VCF file, try to regenerate the VCF file"; }
        fi
    fi

    check_vcf_validity ${hg38_vcf} || \
    { log "The input hg38 VCF file is not valid, please check the file"; return 1; }

    log "Calling crossmap_liftover_hg382hg19 --chain_file ${chain_file} --input_vcf ${hg38_vcf} --output_vcf ${hg19_vcf/.${chrom}/.mixed} --hg19_fasta ${hg19_fasta}"

    crossmap_liftover_hg382hg19 \
    --chain_file ${chain_file} \
    --input_vcf ${hg38_vcf} \
    --output_vcf ${hg19_vcf/.${chrom}/.mixed} \
    --hg19_fasta ${hg19_fasta} || \
    { log "Failed to liftover the gnomAD v4.1 VCF file ${hg38_vcf} from hg38 to hg19 assembly"; return 1; }

    # Now we need to split the VCF file by chromosome as some variants from different chromosomes are mixed in the result hg19 VCF file
    local -a hg19_chrs=($(bcftools query -f '%CHROM\n' ${hg19_vcf/.${chrom}/.mixed} | sort - | uniq - ))
    log "There are ${#hg19_chrs[@]} chromosomes in the hg19 VCF file ${hg19_vcf/.${chrom}/.mixed}: ${hg19_chrs[*]}"

    for chr in "${hg19_chrs[@]}"; do
        if [[ ${chr} == "${chrom}" ]]; then
            bcftools view -Oz -o ${hg19_vcf} ${hg19_vcf/.${chrom}/.mixed} "${chr}" && \
            bcftools index -f -t ${hg19_vcf} && \
            display_vcf ${hg19_vcf}
        else
            bcftools view -Oz -o ${hg19_vcf/.${chrom}/.${chr}} ${hg19_vcf/.${chrom}/.mixed} "${chr}" && \
            bcftools index -f -t ${hg19_vcf/.${chrom}/.${chr}} && \
            display_vcf ${hg19_vcf/.${chrom}/.${chr}}
        fi
    done

    local -a generated_vcfs=($(ls ${output_dir}/*.vcf.gz))
    log "The following VCF files are generated: ${generated_vcfs[*]} after the liftover from hg38 to hg19 assembly on VCF file ${hg38_vcf}"
}



function gnomAD_liftover() {
    local hg38_vcf_chr1=${1}
    local chain_file=${2}
    local hg19_fasta=${3}
    local threads=${4}

    local self_script=$(realpath ${BASH_SOURCE[0]})

    if [[ -z ${threads} ]]; then
        local threads=4
    fi

    local basename_hg38_vcf=${hg38_vcf_chr1/.chr1*/}
    local suffix_hg38_vcf=$(basename ${hg38_vcf_chr1/*.chr1./})

    local hg19_dir=$(dirname ${hg38_vcf_chr1})/hg19
    local -a hg38_chrs=( "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY" )
    # local -a hg38_chrs=( "chr22" "chrX" "chrY" )

    for chr in "${hg38_chrs[@]}"; do
        [[ -d ${hg19_dir}/${chr} ]] || mkdir -p ${hg19_dir}/${chr}
    done

    parallel -j ${threads} --dry-run \
    bash ${self_script} \
    gnomAD_liftover_per_chromosome \
    ${basename_hg38_vcf}.{}.${suffix_hg38_vcf} \
    ${hg19_dir}/{} \
    ${hg19_fasta} \
    ${chain_file} ::: "${hg38_chrs[@]}" && \
    parallel -j ${threads} \
    bash ${self_script} \
    gnomAD_liftover_per_chromosome \
    ${basename_hg38_vcf}.{}.${suffix_hg38_vcf} \
    ${hg19_dir}/{} \
    ${hg19_fasta} \
    ${chain_file} ::: "${hg38_chrs[@]}"

    for chr in "${hg38_chrs[@]}"; do
        local -a tmp_vcfs=($(ls ${hg19_dir}/*/*${chr}.*vcf.gz))
        if [[ ${#tmp_vcfs[@]} -gt 0 ]]; then
            local hg38_vcf=${hg38_vcf_chr1/.chr1/.${chr}}
            local hg19_vcf=${hg38_vcf/.${chr}/.hg19.${chr}}
            local hg19_vcf=${hg19_vcf/.bgz/.gz}
            log "About to concat the VCF files ${tmp_vcfs[*]} to ${hg19_vcf}, the original hg38 VCF file is ${hg38_vcf}"
            [[ -d $(dirname ${hg19_vcf}) ]] && \
            bcftools_concatvcfs -v "${tmp_vcfs[*]}" -o ${hg19_vcf} || \
            { log "Failed to concat the VCF files ${tmp_vcfs[*]} to ${hg19_vcf}, please check the folder $(dirname ${hg19_vcf})"; return 1; }
        else
            log "No VCF files found for chromosome ${chr} at hg19 in folder ${hg19_dir}"
        fi
    done
}



function ClinVar_VCF_deploy() {
    local CACHEDIR=${1}
    local assembly_version=${2}
    local contig_map=${BASE_DIR}/data/liftover/GRC_to_ucsc.contig.map.tsv

    if [[ -z ${assembly_version} ]]; then
        local assembly_version="GRCh37"
    elif [[ ${assembly_version} == "hg19" ]] || [[ ${assembly_version} == "GRCh37" ]]; then
        local assembly_version="GRCh37"
    elif [[ ${assembly_version} == "hg38" ]] || [[ ${assembly_version} == "GRCh38" ]]; then
        local assembly_version="GRCh38"
    else
        log "Currently we only support hg19 and hg38 assemblies"
        return 1
    fi

    if [[ ${assembly_version} == "GRCh37" ]] || [[ ${assembly_version} == "hg19" ]]; then
        local ucsc_assembly_version="hg19"
    elif [[ ${assembly_version} == "GRCh38" ]] || [[ ${assembly_version} == "hg38" ]]; then
        local ucsc_assembly_version="hg38"
    else
        log "Currently we only support hg19 and hg38 assemblies"
        return 1
    fi

    if [[ ! -d ${CACHEDIR} ]]; then
        mkdir -p ${CACHEDIR}
    fi

    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${assembly_version}/clinvar.vcf.gz -O ${CACHEDIR}/clinvar.${assembly_version}.vcf.gz && \
    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${assembly_version}/clinvar.vcf.gz.tbi -O ${CACHEDIR}/clinvar.${assembly_version}.vcf.gz.tbi && \
    log "The ClinVar VCF file is downloaded to ${CACHEDIR}/clinvar.${assembly_version}.vcf.gz" && \
    liftover_from_ucsc_to_GRCh \
    ${CACHEDIR}/clinvar.${assembly_version}.vcf.gz \
    ${contig_map} \
    ${CACHEDIR}/clinvar.${ucsc_assembly_version}.vcf.gz && \
    echo ${CACHEDIR}/clinvar.${ucsc_assembly_version}.vcf.gz
}


function Conservation_install() {
    local CACHEDIR=${1}
    local assembly_version=${2}

    if [[ -z ${assembly_version} ]]; then
        local assembly_version="GRCh38"
    elif [[ ${assembly_version} == "hg38" ]] || [[ ${assembly_version} == "GRCh38" ]]; then
        local assembly_version="GRCh38"
	elif [[ ${assembly_version} == "hg19" ]] || [[ ${assembly_version} == "GRCh37" ]]; then
		local assembly_version="GRCh37"
	else
		log "Currently we only support GRCh37 and GRCh38 assembly versions"
		return 1
	fi

    if [[ ! -d ${CACHEDIR} ]]; then
        mkdir -p ${CACHEDIR}
    fi

    if [[ ! -d ${CACHEDIR}/Conservation ]]; then
        mkdir -p ${CACHEDIR}/Conservation
    fi

	if [[ ${assembly_version} == "GRCh38" ]]; then
		[[ ! -f ${CACHEDIR}/Conservation/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw ]] && \
		wget http://ftp.ensembl.org/pub/current_compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw -O ${CACHEDIR}/Conservation/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw && \
		log "The conservation scores for ${assembly_version} assembly version are downloaded to ${CACHEDIR}/Conservation" && \
		echo ${CACHEDIR}/Conservation/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw || \
		{ log "Failed to download the conservation scores for ${assembly_version} assembly version"; return 1; }
	elif [[ ${assembly_version} == "GRCh37" ]]; then
		wget http://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw -O ${CACHEDIR}/Conservation/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw
		log "The conservation scores for ${assembly_version} assembly version are downloaded to ${CACHEDIR}/Conservation"
		echo ${CACHEDIR}/Conservation/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw
	else
		log "Currently we only support GRCh37 and GRCh38 assembly versions"
		return 1
	fi
}



function LoFtee_install() {
	# Temporarily disabled due to inconsistent cache file support across different assemblies
	# Also deprecated due to the redundant functionality based on CADD and SpliceAI
    local PLUGIN_CACHEDIR=${1}
    local assembly_version=${2}

    log "You need to download the prescores yourself to ${PLUGIN_CACHEDIR}. And the detailed info can be found in the corresponding github repo https://github.com/konradjk/loftee"
    log "Example running command: perl variant_effect_predictor.pl [--other options to VEP] --plugin LoF,loftee_path:/path/to/loftee,human_ancestor_fa:/path/to/human_ancestor.fa.gz"
    if [[ ${assembly_version} == "GRCh37" ]] || [[ ${assembly_version} == "hg19" ]]; then
        local human_ancestor_fasta_url="https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz"
        local human_ancestor_fasta_fai_url="https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai"
        local conservation_file_url="https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/phylocsf_gerp.sql.gz"
    elif [[ ${assembly_version} == "GRCh38" ]] || [[ ${assembly_version} == "hg38" ]]; then
        local gerp_bigwig_url="https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
        local human_ancestor_fasta_url="https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz"
        local human_ancestor_fasta_fai_url="https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai"
        local conservation_file_url="https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/loftee.sql.gz"
    else
        log "Currently we only support GRCh37 and GRCh38 assembly versions"
        return 1
    fi

    mkdir -p ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}
    # Check if the files are already downloaded
    if [[ ! -f ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/human_ancestor.fa.gz ]]; then
        wget ${human_ancestor_fasta_url} -O ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/human_ancestor.fa.gz
    fi

    if [[ ! -f ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/human_ancestor.fa.gz.fai ]]; then
        wget ${human_ancestor_fasta_fai_url} -O ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/human_ancestor.fa.gz.fai
    fi

    if [[ ! -f ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/phylocsf_gerp.sql.gz ]]; then
        wget ${conservation_file_url} -O ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/phylocsf_gerp.sql.gz
    fi

    if [[ ${gerp_bigwig_url} =~ \.gz$ ]] && [[ ! -f ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/gerp_conservation_scores.homo_sapiens.bw ]]; then
        wget ${gerp_bigwig_url} -O ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/gerp_conservation_scores.homo_sapiens.bw
    elif [[ -z ${gerp_bigwig_url} ]]; then
        log "Specifying the hg19/GRCh37 assembly version. So the bigwig file ${gerp_bigwig_url} is not available"
    else
        log "The bigwig file ${gerp_bigwig_url} is already downloaded"
    fi

    # Now we need to return the paths of the downloaded files as a BASH variable which can be used outside this function
    echo "${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/human_ancestor.fa.gz"
    echo "${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/human_ancestor.fa.gz.fai"
    echo "${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/phylocsf_gerp.sql.gz"
    if [[ ${gerp_bigwig_url} =~ \.gz$ ]]; then
        echo "${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/gerp_conservation_scores.homo_sapiens.bw"
    else
        echo ""
    fi
}


# Main installation function
function main_install() {
    local config_file=${1}

    # Read configuration
    local conda_env_yaml=$(read_yaml "$config_file" "conda_env_yaml")
    local vep_cache_dir=$(read_yaml "$config_file" "vep_cache_dir")
    local vep_plugins_dir=$(read_yaml "$config_file" "vep_plugins_dir")
    local assembly=$(read_yaml "$config_file" "assembly")
    local ref_fasta=$(read_yaml "$config_file" "ref_genome")
    local vep_plugins_cachedir=$(read_yaml "$config_file" "vep_plugins_cachedir")
    # ... read other config values ...

    # Perform installation steps
    # 1. Install the conda env
    conda_install_vep "${conda_env_yaml}" || \
    { log "Failed to install the conda env"; return 1; }

    # Test whether currently the conda env is activated
    local conda_env_name=$(head -1 ${conda_env_yaml} | awk -F ': ' '{print $2;}')
    if [[ ${CONDA_PREFIX} =~ ${conda_env_name} ]] && [[ ! -z ${conda_env_name} ]]; then
        log "The conda env $conda_env_name is already activated"
    else
        conda activate $conda_env_name
    fi

    # 2. Install VEP
    # Update config with installation results
    local vep_version=$(vep --help | grep "ensembl-vep" | awk '{print $NF}' | awk -F '.' '{print $1}')
    # Update config with installation results
    update_yaml "$config_file" "vep_installed_version" "$vep_version"

    # Install VEP API and caches and plugins first
    vep_install_wrapper \
    --VEP_CACHEDIR "$vep_cache_dir" \
    --VEP_PLUGINSDIR "$vep_plugins_dir" \
    --VEP_ASSEMBLY "$assembly" \
    --VEP_PLUGINSCACHEDIR "$vep_plugins_cachedir" || \
    { log "Failed to install VEP API and caches"; return 1; }

    # Install VEP plugins caches
    # 3. Install VEP plugins
    VEP_plugins_install \
    ${vep_plugins_dir} \
    ${vep_plugins_cachedir} \
    ${assembly} \
    ${config_file} \
    ${conda_env_name} || \
    { log "Failed to install VEP plugins"; return 1; }


    # 4. Install gnomAD VCF (basically download bgzipped VCF files)
    local gnomad_vcf_dir=$(read_yaml "$config_file" "gnomad_vcf_dir")
    if [[ ${assembly} =~ "GRCh37" ]] || [[ ${assembly} =~ "hg19" ]]; then
        # Chain file is small enough to be included in the git repo
        local chain_file=$(read_yaml "$config_file" "chain_file")
    fi
    local gnomAD_chrX_vcf=$(gnomAD_install ${gnomad_vcf_dir} ${ref_fasta} | tail -1)
    [[ -f ${gnomAD_chrX_vcf} ]] && \
    update_yaml "$config_file" "gnomad_vcf_chrX" "${gnomAD_chrX_vcf}" || \
    { log "Failed to install gnomAD VCF"; return 1; }

    # 5. Install ClinVar VCF
    local clinvar_vcf_dir=$(read_yaml "$config_file" "clinvar_vcf_dir")
    local prepared_clinvar_vcf=$(ClinVar_VCF_deploy ${clinvar_vcf_dir} ${assembly} | tail -1)
    [[ -f ${prepared_clinvar_vcf} ]] && \
    update_yaml "$config_file" "clinvar_vcf" "${prepared_clinvar_vcf}" || \
    { log "Failed to install ClinVar VCF"; return 1; }

    # 6. Install CADD prescores
    CADD_install \
    ${config_file} \
    ${assembly} || \
    { log "Failed to install CADD prescores"; return 1; }
}



if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
    "$@"
fi
