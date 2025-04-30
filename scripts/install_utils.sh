#! /usr/bin/env bash
SCRIPT_DIR=$(dirname $(realpath ${BASH_SOURCE[0]}))
source ${SCRIPT_DIR}/common_bash_utils.sh
source ${SCRIPT_DIR}/annotation_vcf.sh
log "The base directory is ${BASE_DIR}, the scripts directory is ${SCRIPT_DIR}"


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


function vep_cache_api_install() {
    local VEP_CACHEDIR=""
    local VEP_ASSEMBLY=""
    local VEP_PLUGINS=""
    local VEP_PLUGINSDIR=""
    local VEP_VERSION=""
    local VEP_PLUGINSCACHEDIR=""

    local TEMP
    TEMP=$(getopt -o hc:y:r:p:v: --long help,VEP_CACHEDIR:,VEP_ASSEMBLY:,VEP_PLUGINSDIR:,VEP_PLUGINSCACHEDIR:,VEP_PLUGINS:,VEP_VERSION: -- "$@")

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
                echo "  -v, --VEP_VERSION         Set the version of VEP to install"
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
            -v|--VEP_VERSION)
                VEP_VERSION="$2"
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

    [[ ${VEP_ASSEMBLY} =~ "hg19" ]] && VEP_ASSEMBLY="GRCh37"
    [[ ${VEP_ASSEMBLY} =~ "hg38" ]] && VEP_ASSEMBLY="GRCh38"

    # Test the PERL5LIB value and PATH value to see if they already include the VEP_DESTDIR and VEP_PLUGINSDIR
    # If yes for both, then we can directly skip the follow up installation of VEP API
    local VEP_DESTDIR=$(perl -e 'print join("\n", @INC);' | grep site_perl | head -1)
    [[ $(echo $PERL5LIB) =~ "${VEP_DESTDIR}" ]] && \
    [[ $(echo $PATH) =~ "${VEP_DESTDIR}" ]] && \
    [[ $(echo $PERL5LIB}) =~ "${VEP_PLUGINSDIR}" ]] && \
    [[ -d ${VEP_CACHEDIR}/homo_sapiens_merged/${VEP_VERSION}_${VEP_ASSEMBLY} ]] && \
    [[ -d ${VEP_CACHEDIR}/homo_sapiens_merged/${VEP_VERSION}_${VEP_ASSEMBLY}/1 ]] && \
    [[ -f ${VEP_CACHEDIR}/homo_sapiens_merged/${VEP_VERSION}_${VEP_ASSEMBLY}/1/all_vars.gz ]] && \
    [[ -f ${VEP_CACHEDIR}/homo_sapiens_merged/${VEP_VERSION}_${VEP_ASSEMBLY}/chr_synonyms.txt ]] && \
    [[ -f ${VEP_PLUGINSDIR}/SpliceVault.pm ]] && \
    [[ -f ${VEP_PLUGINSDIR}/SpliceAI.pm ]] && \
    [[ -f ${VEP_PLUGINSDIR}/UTRAnnotator.pm ]] && \
    [[ -f ${VEP_PLUGINSDIR}/AlphaMissense.pm ]] && \
    [[ -f ${VEP_PLUGINSDIR}/GO.pm ]] && \
    [[ -f ${VEP_PLUGINSDIR}/Phenotypes.pm ]] && \
    [[ -f ${VEP_PLUGINSDIR}/Conservation.pm ]] && \
    [[ -f ${VEP_PLUGINSDIR}/LOEUF.pm ]] && \
    [[ -f ${VEP_PLUGINSDIR}/SameCodon.pm ]] && \
    [[ -f ${VEP_PLUGINSDIR}/NMD.pm ]] && \
    [[ -f ${VEP_PLUGINSDIR}/PrimateAI.pm ]] && \
    [[ -f ${VEP_PLUGINSDIR}/MaveDB.pm ]] && \
    { log "The PERL5LIB value and PATH value already include ${VEP_DESTDIR} and ${VEP_PLUGINSDIR}, indicating the following installation process has already been succesfully performed. Skip the function for now."; return 0; }

    local conda_env_name=$(basename $CONDA_PREFIX)
    [[ ${VEP_DESTDIR} =~ ${conda_env_name} ]] && log "The dest dir is set to the directory (${VEP_DESTDIR}) where perl modules are installed by conda" || \
    { log "Since the function is designed to perform follow up installation of VEP upon the installation of VEP dependencies via conda, we only accept installing VEP API at the perl module installation location previously used by conda"; return 1; }

    # Construct the command
    local cmd="vep_install -d ${VEP_DESTDIR} --AUTO acp -s homo_sapiens_merged --NO_HTSLIB --NO_BIOPERL --CONVERT"
    [[ -n "$VEP_CACHEDIR" ]] && cmd+=" --CACHEDIR $VEP_CACHEDIR"
    [[ -n "$VEP_ASSEMBLY" ]] && cmd+=" --ASSEMBLY $VEP_ASSEMBLY"
    [[ -n "$VEP_PLUGINS" ]] && [[ $VEP_PLUGINS != "empty" ]] && cmd+=" --PLUGINS $VEP_PLUGINS,SpliceAI,SpliceVault,UTRAnnotator,AlphaMissense,GO,Phenotypes,Conservation,LOEUF,SameCodon,NMD,PrimateAI,MaveDB" || cmd+=" --PLUGINS SpliceAI,SpliceVault,UTRAnnotator,AlphaMissense,GO,Phenotypes,Conservation,LOEUF,SameCodon,NMD,PrimateAI,MaveDB"
    [[ -n "$VEP_PLUGINSDIR" ]] && cmd+=" --PLUGINSDIR $VEP_PLUGINSDIR"

    # Execute the command
    # After executing the command, we need to bind the new PERL5LIB value with the current conda env
    log "Now we start to install the VEP api and downloading caches (which might take a while to finish). So pls be patient, here is the command we are going to execute: ${cmd}"
    $cmd && \
    [[ ! ${PERL5LIB} == *${VEP_DESTDIR}* ]] && conda env config vars set PERL5LIB="$VEP_DESTDIR:$VEP_PLUGINSDIR" || log "The PERL5LIB value is already bound with the current conda env ${CONDA_ENV_NAME}"
    [[ ! ${PATH} == *${VEP_DESTDIR}* ]] && conda env config vars set PATH="$VEP_DESTDIR:$VEP_DESTDIR/htslib:$PATH" || log "The PATH value is already bound with the current conda env ${CONDA_ENV_NAME}"
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
    # First install UTRAnnotator
    log "The Cache directory to store the requried annotation files for VEP Plugins is ${VEP_PLUGINSCACHEDIR}"
    UTRAnnotator_install ${config_file} ${VEP_PLUGINSCACHEDIR} || \
    { log "Failed to install UTRAnnotator"; return 1; }

    # Then install LOEUF
    LOEUF_install ${config_file} ${ASSEMBLY_VERSION} ${BASE_DIR} || \
    { log "Failed to install LOEUF"; return 1; }

    # Then install AlphaMissense
    AlphaMissense_install ${config_file} ${VEP_PLUGINSCACHEDIR} ${ASSEMBLY_VERSION} ${VEP_PLUGINSDIR} || \
    { log "Failed to install AlphaMissense"; return 1; }
    # Then annotate the AlphaMissense prescore file to vcf file
    AlphaMissense_anno ${config_file} || \
    { log "Failed to convert the AlphaMissense prescore file to vcf file and annotate protein domains to it"; return 1; }
    AlphaMissense_pick_intolerant_motifs ${config_file} || \
    { log "Failed to pick intolerant motifs from the AlphaMissense prescore file"; return 1; }

    # Then install SpliceAI
    SpliceAI_install ${config_file} ${VEP_PLUGINSCACHEDIR} ${VEP_PLUGINSDIR} || \
    { log "Failed to install SpliceAI"; return 1; }

    # Followed by SpliceVault
    SpliceVault_install ${config_file} ${VEP_PLUGINSCACHEDIR} ${SCRIPT_DIR} || \
    { log "Failed to install SpliceVault"; return 1; }

    # Followed by MaveDB
    MaveDB_install ${config_file} ${VEP_PLUGINSCACHEDIR} || \
    { log "Failed to install MaveDB"; return 1; }

    # Last, install PrimateAI
    PrimateAI_install ${config_file} ${VEP_PLUGINSCACHEDIR} ${VEP_PLUGINSDIR} || \
    { log "Failed to install PrimateAI"; return 1; }


    # Install Conservation file if assembly version is hg38
    log "Now we start to install the conservation file"
    local conservation_file=$(Conservation_install ${config_file} ${VEP_PLUGINSCACHEDIR} ${ASSEMBLY_VERSION}) || \
    { log "Failed to install Conservation file"; return 1; }
}


function UTRAnnotator_install() {
    local config_file=${1}
    local PLUGIN_CACHEDIR=${2}
    local assembly_version=$(read_yaml ${config_file} "assembly")

    local utr_annotator_file=$(read_yaml ${config_file} "utr_annotator_file")
    [[ -f ${utr_annotator_file} ]] && \
    log "The UTRAnnotator plugin is already downloaded for both hg19 and hg38 assemblies" && return 0

    [[ ! -d ${PLUGIN_CACHEDIR} ]] && { log "The cache directory ${PLUGIN_CACHEDIR} is not found, please check the file"; return 1; }

    if [[ ! -d ${PLUGIN_CACHEDIR}/UTRannotator ]] || \
       [[ ! -f ${PLUGIN_CACHEDIR}/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt ]] || \
       [[ ! -f ${PLUGIN_CACHEDIR}/UTRannotator/uORF_5UTR_GRCh37_PUBLIC.txt ]]; then
        cd ${PLUGIN_CACHEDIR} && \
        git clone https://github.com/Ensembl/UTRannotator && \
        log "${PLUGIN_CACHEDIR}/UTRannotator/uORF_5UTR_GRCh37_PUBLIC.txt and ${PLUGIN_CACHEDIR}/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt are well-downloaded, remember to specify them when using UTRAnnotator"
        log "A typical command syntax would be: vep -i variations.vcf --plugin UTRAnnotator,file=${PLUGIN_CACHEDIR}/UTRannotator/uORF_starts_ends_GRCh38_PUBLIC.txt"
    else
        log "The UTRannotator plugin is already downloaded for both hg19 and hg38 assemblies"
    fi

    if [[ ${assembly_version} =~ "GRCh37" ]] || [[ ${assembly_version} =~ "hg19" ]]; then
        update_yaml ${config_file} "utr_annotator_file" "${PLUGIN_CACHEDIR}/UTRannotator/uORF_5UTR_GRCh37_PUBLIC.txt"
    elif [[ ${assembly_version} =~ "GRCh38" ]] || [[ ${assembly_version} =~ "hg38" ]]; then
        update_yaml ${config_file} "utr_annotator_file" "${PLUGIN_CACHEDIR}/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt"
    else
        log "Not supported assembly version: ${assembly_version}"
        return 1
    fi
}


function SpliceAI_install() {
    local config_file=${1}
    local PLUGIN_CACHEDIR=${2}
    local PLUGIN_DIR=${3}

    local snv_prescore=$(read_yaml ${config_file} "spliceai_snv_prescore")
    local indel_prescore=$(read_yaml ${config_file} "spliceai_indel_prescore")

    [[ -f ${snv_prescore} ]] && [[ -f ${indel_prescore} ]] && \
    log "The prescores for the SpliceAI plugin are already downloaded, skip the installation process" && return 0

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
    [[ ! -d ${PLUGIN_CACHEDIR} ]] && { log "The cache directory ${PLUGIN_CACHEDIR}/SpliceAI is not found, please check the file"; return 1; }
    [[ ! -d ${PLUGIN_DIR} ]] && { log "The plugins directory ${PLUGIN_DIR} is not found, please check the file"; return 1; }

    read -p "Have you finished downloading the prescores for the SpliceAI plugin? (yes or no)"
    if [[ ${REPLY} =~ "yes" ]] || [[ ${REPLY} =~ "y" ]] || [[ ${REPLY} =~ "Y" ]] || [[ ${REPLY} =~ "Yes" ]] || [[ ${REPLY} =~ "YES" ]]; then
        log "Now the installation of the SpliceAI plugin is completed. Please specify the paths to the prescore files (one for snv and one for indel)"
        read -p "Please specify the absolute path to the snv prescore file, remember the file should correspond to the assembly version of your VEP installation"
        local snv_pre_score_file=${REPLY}
        if [[ ! -f ${snv_pre_score_file} ]]; then
            log "The file ${snv_pre_score_file} does not exist"
            return 1
        fi
        update_yaml ${config_file} "spliceai_snv_prescore" "${snv_pre_score_file}"
        read -p "Please specify the absolute path to the indel prescore file, remember the file should correspond to the assembly version of your VEP installation"
        local indel_pre_score_file=${REPLY}
        if [[ ! -f ${indel_pre_score_file} ]]; then
            log "The file ${indel_pre_score_file} does not exist"
            return 1
        fi
        log "Now we start to install the SpliceAI plugin"
        # We need to return the prescore file paths to the outside of the function
        update_yaml ${config_file} "spliceai_indel_prescore" "${indel_pre_score_file}"
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
    local cadd_version=$(read_yaml "${config_file}" "cadd_version")
    # Note that by default, CADD will use the CADD cache file in the conda environment
    # Ask the user whether they want to use the default CADD cache directory to store the genome annotation file and the prescores.
    read -p "Do you want to use the default CADD cache directory ${CADD_cache_dir} to store the genome annotation file and the prescores which could be take up to 1TB space? (yes or no)"
    if [[ ${REPLY} =~ "yes" ]] || [[ ${REPLY} =~ "y" ]] || [[ ${REPLY} =~ "Y" ]] || [[ ${REPLY} =~ "Yes" ]] || [[ ${REPLY} =~ "YES" ]]; then
        local CADD_prescore_dir=${CADD_cache_dir}/prescored
        local CADD_anno_dir=${CADD_cache_dir}/annotations
        update_yaml ${config_file} "cadd_base_dir" "$(dirname ${CADD_script})"
    else
        read -p "In this case, you need to download the CADD repo as a zip file and unzip it to a local directory and all the cache files will be store in that directory. Please specify the absolute path to the directory: "
        local CADD_parent_dir=${REPLY}

        # Check if the CADD scripts directory exists, if not, download the zip file and unzip it
        local CADD_base_dir=$(find ${CADD_parent_dir}/ -maxdepth 1 -type d -name "CADD-scripts-*${cadd_version#v}*" -print | head -n1)
        [[ ! -d ${CADD_base_dir} ]] && \
        local CADD_zip_download_url=$(read_yaml "${config_file}" "cadd_zip_download_url") && \
        wget ${CADD_zip_download_url} -O ${CADD_parent_dir}/CADD-scripts.zip && \
        unzip ${CADD_parent_dir}/CADD-scripts.zip -d ${CADD_parent_dir}/ && \
        rm ${CADD_parent_dir}/CADD-scripts.zip

        # Get the base folder name from the unzipped directory
        [[ -z ${CADD_base_dir} ]] && { log "Could not find CADD scripts directory"; return 1; }
        update_yaml ${config_file} "cadd_base_dir" "${CADD_base_dir}" && \
        log "Now the CADD base directory is ${CADD_base_dir} and it's updated to the config file ${config_file}"
        local CADD_script="${CADD_base_dir}/CADD.sh" && \
        local CADD_cache_dir="${CADD_base_dir}/data" && \
        local CADD_prescore_dir="${CADD_cache_dir}/prescored" && \
        local CADD_anno_dir="${CADD_cache_dir}/annotations" && \
        log "Now the CADD script is ${CADD_script}, the CADD cache directory is ${CADD_cache_dir}, the CADD prescore directory is ${CADD_prescore_dir}, the CADD annotation directory is ${CADD_anno_dir}"
    fi


    # Now we start downloading the CADD genome annotations corresponding to the assembly version
    if [[ ${assembly} == "GRCh37" ]] || [[ ${assembly} == "hg19" ]]; then
        if [[ ! -d ${CADD_anno_dir}/GRCh37_${cadd_version} ]]; then
            wget -c https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh37/GRCh37_v1.7.tar.gz -O ${CADD_anno_dir}/GRCh37_v1.7.tar.gz && \
            md5sum ${CADD_anno_dir}/GRCh37_v1.7.tar.gz | grep -q $(read_yaml "${config_file}" "cadd_GRCh37_anno_md5") && \
            tar -xzvf ${CADD_anno_dir}/GRCh37_v1.7.tar.gz -C ${CADD_anno_dir}/ && \
            rm ${CADD_anno_dir}/GRCh37_v1.7.tar.gz
        else
            log "The CADD genome annotation for GRCh37_${cadd_version} is already downloaded, skip the downloading process"
        fi
    elif [[ ${assembly} == "GRCh38" ]] || [[ ${assembly} == "hg38" ]]; then
        if [[ ! -d ${CADD_anno_dir}/GRCh38_${cadd_version} ]]; then
            wget -c https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/GRCh38_v1.7.tar.gz -O ${CADD_anno_dir}/GRCh38_v1.7.tar.gz && \
            md5sum ${CADD_anno_dir}/GRCh38_v1.7.tar.gz | grep -q $(read_yaml "${config_file}" "cadd_GRCh38_anno_md5") && \
            tar -xzvf ${CADD_anno_dir}/GRCh38_v1.7.tar.gz -C ${CADD_anno_dir}/ && \
            rm ${CADD_anno_dir}/GRCh38_v1.7.tar.gz
        else
            log "The CADD genome annotation for GRCh38_${cadd_version} is already downloaded, skip the downloading process"
        fi
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
    log "Now the CADD prescore URL for ${assembly} is ${snv_prescore_url} and ${indel_prescore_url}"

    local snv_file_name=$(basename ${snv_prescore_url})
    local indel_file_name=$(basename ${indel_prescore_url})

    if [[ ! -f ${CADD_prescore_dir}/${assembly}_${cadd_version}/incl_anno/${snv_file_name} ]]; then
        wget -c ${snv_prescore_url} -O ${CADD_prescore_dir}/${assembly}_${cadd_version}/incl_anno/${snv_file_name} && \
        md5sum ${CADD_prescore_dir}/${assembly}_${cadd_version}/incl_anno/${snv_file_name} | grep -q ${snv_prescore_md5} && \
        wget -c ${snv_prescore_url}.tbi -O ${CADD_prescore_dir}/${assembly}_${cadd_version}/incl_anno/${snv_file_name}.tbi
    else
        log "The CADD prescore for ${assembly} is already downloaded at ${CADD_prescore_dir}/${assembly}_${cadd_version}/incl_anno/${snv_file_name}, skip the downloading process"
    fi

    if [[ ! -f ${CADD_prescore_dir}/${assembly}_${cadd_version}/incl_anno/${indel_file_name} ]]; then
        wget -c ${indel_prescore_url} -O ${CADD_prescore_dir}/${assembly}_${cadd_version}/incl_anno/${indel_file_name} && \
        md5sum ${CADD_prescore_dir}/${assembly}_${cadd_version}/incl_anno/${indel_file_name} | grep -q ${indel_prescore_md5} && \
        wget -c ${indel_prescore_url}.tbi -O ${CADD_prescore_dir}/${assembly}_${cadd_version}/incl_anno/${indel_file_name}.tbi
    else
        log "The CADD prescore for ${assembly} is already downloaded at ${CADD_prescore_dir}/${assembly}_${cadd_version}/incl_anno/${indel_file_name}, skip the downloading process"
    fi
}



function PrimateAI_install() {
    local config_file=${1}
    local PLUGIN_CACHEDIR=${2}
    local PLUGIN_DIR=${3}

    local prescore_file=$(read_yaml ${config_file} "primateai_prescore")
    [[ -f ${prescore_file} ]] && \
    log "The prescore file for the PrimateAI plugin is already downloaded, skip the installation process" && return 0

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
            update_yaml ${config_file} "primateai_prescore" "${prescore_file/.tsv.gz/.sorted.tsv.bgz}"
        elif [[ -f ${prescore_file} ]] && [[ ${prescore_file} =~ \.tsv.bgz$ ]]; then
            update_yaml ${config_file} "primateai_prescore" "${prescore_file}"
        else
            log "The file ${prescore_file} does not exist"
            return 1
        fi
    else
        log "Please finish formatting the downloaded files first"
        return 1
    fi
}


function MaveDB_install() {
    local config_file=${1}
    local PLUGIN_CACHEDIR=${2}

    if [[ -z ${PLUGIN_CACHEDIR} ]]; then
        local PLUGIN_CACHEDIR=$(read_yaml ${config_file} "vep_plugins_dir")
    fi

    local assembly=$(read_yaml ${config_file} "assembly")
    local mavedb_url=$(read_yaml ${config_file} "mavedb_url")
    local mavedb_file=$(read_yaml ${config_file} "mavedb_file")

    if [[ ! -d ${PLUGIN_CACHEDIR}/MaveDB ]]; then
        mkdir -p ${PLUGIN_CACHEDIR}/MaveDB
    fi

    if [[ ${assembly} == "GRCh37" ]] || [[ ${assembly} == "hg19" ]]; then
        log "The MaveDB file is not available for ${assembly}, skip the downloading process"
        update_yaml ${config_file} "mavedb_file" ""
        return 0
    fi

    if [[ -f ${mavedb_file} ]] && \
       [[ -f ${mavedb_file}.tbi ]] && \
       [[ ${mavedb_file} -ot ${mavedb_file}.tbi ]]; then
        log "The MaveDB file is already downloaded, skip the downloading process"
    else
        wget -c ${mavedb_url} -O ${PLUGIN_CACHEDIR}/MaveDB/MaveDB_variants.${assembly}.tsv.gz && \
        wget -c ${mavedb_url}.tbi -O ${PLUGIN_CACHEDIR}/MaveDB/MaveDB_variants.${assembly}.tsv.gz.tbi && \
        update_yaml ${config_file} "mavedb_file" "${PLUGIN_CACHEDIR}/MaveDB/MaveDB_variants.${assembly}.tsv.gz"
    fi

}



function AlphaMissense_install() {
    local config_file=${1}
    local PLUGIN_CACHEDIR=${2}
    local ASSEMBLY_VERSION=${3}
    local PLUGIN_DIR=${4}

    log "For details, read ${PLUGIN_DIR}/AlphaMissense.pm"
    local hg19_url="https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz"
    local hg38_url="https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"

    local alphamissense_prescore=$(read_yaml ${config_file} "alphamissense_prescore")
    [[ -f ${alphamissense_prescore} ]] && \
    log "The prescore file for the AlphaMissense plugin is already downloaded, skip the installation process" && return 0

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
        update_yaml ${config_file} "alphamissense_prescore" "${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg19.tsv.gz"
    elif [[ ${ASSEMBLY_VERSION} == "GRCh38" ]] || [[ ${ASSEMBLY_VERSION} == "hg38" ]]; then
        update_yaml ${config_file} "alphamissense_prescore" "${PLUGIN_CACHEDIR}/AlphaMissense/AlphaMissense_hg38.tsv.gz"
    else
        log "Not supported assembly version: ${ASSEMBLY_VERSION}"
        return 1
    fi
}



function AlphaMissense_anno() {
    local config_file=${1}

    local alphamissense_prescore=$(read_yaml ${config_file} "alphamissense_prescore")
    local alphamissense_vcf=$(read_yaml ${config_file} "alphamissense_vcf")
    local alphamissense_vep_vcf=$(read_yaml ${config_file} "alphamissense_vep_vcf")
    local conservation_file=$(read_yaml ${config_file} "conservation_file")


    check_vcf_validity ${alphamissense_vcf} && \
    [[ ${alphamissense_vcf} -nt ${alphamissense_prescore} ]] && \
    check_vcf_validity ${alphamissense_vep_vcf} && \
    [[ ${alphamissense_vep_vcf} -nt ${alphamissense_vcf} ]] && \
    check_vcf_infotags ${alphamissense_vep_vcf} "CSQ" && \
    [[ ${alphamissense_vep_vcf} -nt ${SCRIPT_DIR}/install_utils.sh ]] && \
    log "The AlphaMissense VCF file ${alphamissense_vcf} is already annotated by VEP to ${alphamissense_vep_vcf}. Skip this step" && return 0
    

    local convert_py=${SCRIPT_DIR}/alphmis_tsv2vcf.py
    python ${convert_py} \
    ${alphamissense_prescore} \
    ${alphamissense_prescore/.tsv/.vcf} && \
    local vep_cache_dir=$(read_yaml ${config_file} "vep_cache_dir") && \
    local threads=$(read_yaml ${config_file} "threads") && \
    local ref_genome=$(read_yaml ${config_file} "ref_genome") && \
    local assembly=$(read_yaml ${config_file} "assembly") && \
    local vep_plugins_dir=$(read_yaml ${config_file} "vep_plugins_dir") && \
    basic_vep_annotation \
    -i ${alphamissense_prescore/.tsv/.vcf} \
    -a ${assembly} \
    -r ${ref_genome} \
    -c ${vep_cache_dir} \
    -t ${threads} \
    -p ${vep_plugins_dir} \
    -v ${conservation_file} && \
    display_vcf ${alphamissense_prescore/.tsv/.vep.vcf} && \
    update_yaml ${config_file} "alphamissense_vep_vcf" "${alphamissense_prescore/.tsv/.vep.vcf}" && \
    log "The AlphaMissense VCF file ${alphamissense_prescore/.tsv/.vep.vcf} is annotated by VEP and saved to ${alphamissense_vcf}"
}


function AlphaMissense_pick_intolerant_motifs() {
    local config_file=${1}
    local intolerant_motif_prefix=${2}
    local pick_intolerant_motifs_py=${SCRIPT_DIR}/am_pick_intolerant_motifs.py

    if [[ -z ${intolerant_motif_prefix} ]]; then
        log "The intolerant motif prefix is not provided, we will use the default prefix: /paedyl01/disk1/yangyxt/public_data/VEP_plugins_caches/AlphaMissense/AlphaMissense_hg19.kde"
        local assembly=$(read_yaml ${config_file} "assembly")
        local intolerant_motif_prefix="/paedyl01/disk1/yangyxt/public_data/VEP_plugins_caches/AlphaMissense/AlphaMissense_${assembly}.kde"
    fi

    local alphamissense_vep_vcf=$(read_yaml ${config_file} "alphamissense_vep_vcf")
    local alphamissense_intolerant_motifs=$(read_yaml ${config_file} "alphamissense_intolerant_motifs")
    local threads=$(read_yaml ${config_file} "threads")
    [[ -f ${alphamissense_intolerant_motifs} ]] && \
    [[ ${alphamissense_intolerant_motifs} -nt ${alphamissense_vep_vcf} ]] && \
    [[ ${alphamissense_intolerant_motifs} -nt ${pick_intolerant_motifs_py} ]] && \
    log "The AlphaMissense intolerant motifs file ${alphamissense_intolerant_motifs} is already generated, skip this step" && \
    return 0 || \
    log "The AlphaMissense intolerant motifs file ${alphamissense_intolerant_motifs} is not generated, we will generate it"

    [[ -d $(dirname ${alphamissense_intolerant_motifs}) ]] || \
    log "The directory $(dirname ${alphamissense_intolerant_motifs}) does not exist, we will use the prefix ${intolerant_motif_prefix}" && \
    local alphamissense_intolerant_motifs="${intolerant_motif_prefix}.pkl"

    python ${pick_intolerant_motifs_py} \
    --vcf_path ${alphamissense_vep_vcf} \
    --output ${alphamissense_intolerant_motifs} \
    --processes ${threads} && \
    update_yaml ${config_file} "alphamissense_intolerant_motifs" "${alphamissense_intolerant_motifs}" && \
    log "The AlphaMissense intolerant motifs file is generated and saved to ${alphamissense_intolerant_motifs}"
}


function AlphaMissense_stat() {
    local config_file=${1}
    local alphamissense_stat=$(read_yaml ${config_file} "alphamissense_pd_stat")
    local alphamissense_vcf=$(read_yaml ${config_file} "alphamissense_vep_vcf")
    local stat_py=${SCRIPT_DIR}/stat_protein_domain_amscores.py

    [[ -f ${alphamissense_stat} ]] && \
    [[ ${alphamissense_stat} -nt ${alphamissense_vcf} ]] && \
    [[ ${alphamissense_stat} -nt ${stat_py} ]] && \
    log "The AlphaMissense statistics JSON file ${alphamissense_stat} is already generated, skip this step" && \
    return 0 || \
    log "The AlphaMissense statistics JSON file ${alphamissense_stat} is not generated or updated, we will generate it"

    python ${stat_py} \
    ${alphamissense_vcf} \
    ${alphamissense_vcf/.vcf*/.prot.domain.stats.pkl} && \
    update_yaml ${config_file} "alphamissense_pd_stat" "${alphamissense_vcf/.vcf*/.prot.domain.stats.pkl}" && \
    log "The AlphaMissense statistics pickle file ${alphamissense_vcf/.vcf*/.prot.domain.stats.pkl} is generated and saved to ${alphamissense_stat}"
}


function AlphaMissense_pick_intolerant_domains() {
    local config_file=${1}
    
    local assembly=$(read_yaml ${config_file} "assembly")
    local alphamissense_pd_stat=$(read_yaml ${config_file} "alphamissense_pd_stat")
    local alphamissense_tranx_domain_map=$(read_yaml ${config_file} "alphamissense_tranx_domain_map")
    local alphamissense_intolerant_domains=$(read_yaml ${config_file} "alphamissense_intolerant_domains")
    local pick_intolerant_domains_py=${SCRIPT_DIR}/am_pick_intolerant_domains.py

    [[ -f ${alphamissense_intolerant_domains} ]] && \
    [[ ${alphamissense_intolerant_domains} -nt ${alphamissense_pd_stat} ]] && \
    [[ -f ${alphamissense_tranx_domain_map} ]] && \
    [[ ${alphamissense_tranx_domain_map} -nt ${alphamissense_pd_stat} ]] && \
    [[ ${alphamissense_intolerant_domains} -nt ${pick_intolerant_domains_py} ]] && \
    log "The AlphaMissense intolerant domains file ${alphamissense_intolerant_domains} is already generated, skip this step" && return 0
    
    local threads=$(read_yaml ${config_file} "threads")
    local alphamissense_dir=$(dirname ${alphamissense_pd_stat})
    
    log "Running the AlphaMissense intolerant domains analysis with this command: python ${pick_intolerant_domains_py} --pickle_file ${alphamissense_pd_stat} --threads ${threads} --output_dir ${alphamissense_dir}"
    python ${pick_intolerant_domains_py} \
    --assembly ${assembly} \
    --pickle_file ${alphamissense_pd_stat} \
    --threads ${threads} \
    --output_dir ${alphamissense_dir} && \
    update_yaml ${config_file} "alphamissense_intolerant_domains" "${alphamissense_dir}/domain_tolerance_analysis.${assembly}.tsv" && \
    update_yaml ${config_file} "alphamissense_tranx_domain_map" "${alphamissense_dir}/transcript_exon_domain_mapping.${assembly}.pkl" && \
    update_yaml ${config_file} "all_intolerant_domains" "${alphamissense_dir}/all_intolerant_domains.${assembly}.pkl" && \
    log "The AlphaMissense intolerant domains file ${alphamissense_intolerant_domains} and transcript exon domain mapping file ${alphamissense_tranx_domain_map} are generated and saved to ${alphamissense_dir}"
}



function LOEUF_install() {
    # The LOEUF plugin cache files are stored in the github repo because the original url is not accessible anymore
    local config_file=${1}
    local assembly_version=${2}
    local git_repo_dir=${3}

    local loeuf_hg19_file=${git_repo_dir}/data/loeuf/loeuf_dataset.tsv.gz
    local loeuf_hg38_file=${git_repo_dir}/data/loeuf/loeuf_dataset_hg38.tsv.gz

    if [[ ${assembly_version} =~ "GRCh37" ]] || [[ ${assembly_version} =~ "hg19" ]]; then
        update_yaml ${config_file} "loeuf_prescore" "${loeuf_hg19_file}"
    elif [[ ${assembly_version} =~ "GRCh38" ]] || [[ ${assembly_version} =~ "hg38" ]]; then
        update_yaml ${config_file} "loeuf_prescore" "${loeuf_hg38_file}"
    else
        log "Not supported assembly version: ${assembly_version}"
        return 1
    fi
}


function gnomAD_install() {
    local config_file=${1}
    local CACHEDIR=${2}
    local assembly=${3}
    # Assembly is the fasta file path
    # If assembly is not provided, we will use GRCh38 assembly

    if [[ -z ${assembly} ]]; then
        local assembly="GRCh38"
    fi

    # We need to iterate over all chromosomes to check if the files are already downloaded and valid
    # local -a chromosomes=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
    local -a chromosomes=( "22" "X" "Y" )
    local gnomAD_chrX_vcf=$(read_yaml ${config_file} "gnomad_vcf_chrX")

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
        parallel -j4 bash ${SCRIPT_DIR}/common_bash_utils.sh check_vcf_validity ${CACHEDIR}/gnomad.joint.v4.1.sites.hg19.chr{}.vcf.gz ::: "${chromosomes[@]}" && \
        update_yaml ${config_file} "gnomad_vcf_chrX" "${CACHEDIR}/gnomad.joint.v4.1.sites.hg19.chrX.vcf.gz" || \
        { log "Failed to liftover the gnomAD v4.1 files to hg19 assembly"; return 1; }
    else
        update_yaml ${config_file} "gnomad_vcf_chrX" "${CACHEDIR}/gnomad.joint.v4.1.sites.chrX.vcf.bgz"
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
            log "The hg19 VCF file ${hg19_vcf} is already indexed and the VCF file is valid and updated"
            return 0
        else
            tabix -f -p vcf ${hg19_vcf} && \
            log "The hg19 VCF file ${hg19_vcf} is already indexed and the VCF file is valid and updated" && \
            return 0 || \
            { log "Failed to index the hg19 VCF file ${hg19_vcf}, try to regenerate the VCF file"; }
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
    # local -a hg38_chrs=( "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY" )
    local -a hg38_chrs=( "chr22" "chrX" "chrY" )

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
            check_vcf_validity ${hg19_vcf} || { \
            [[ -d $(dirname ${hg19_vcf}) ]] && \
            bcftools_concatvcfs -v "${tmp_vcfs[*]}" -o ${hg19_vcf} || \
            { log "Failed to concat the VCF files ${tmp_vcfs[*]} to ${hg19_vcf}, please check the folder $(dirname ${hg19_vcf})"; return 1; } }
        else
            log "No VCF files found for chromosome ${chr} at hg19 in folder ${hg19_dir}"
        fi
    done
}



function basic_vep_annotation() {
    # This function can be used to annotate ClinVar VCF file
    # This function can be used to annotate AlphaMissense VCF file (converted from TSV file)

    local input_vcf
    local assembly
    local ref_genome
    local vep_cache_dir
    local threads
    local vep_plugins_dir
    local spliceai_snv_prescore
    local spliceai_indel_prescore
    local config_file

    # Parse arguments using getopts
    local OPTIND=1
    while getopts "i:a::r::c::t::p::s::d::f::v::" opt; do
        case ${opt} in
            i) input_vcf=${OPTARG} ;;
            a) assembly=${OPTARG} ;;
            r) ref_genome=${OPTARG} ;;
            c) vep_cache_dir=${OPTARG} ;;
            t) threads=${OPTARG} ;;
            p) vep_plugins_dir=${OPTARG} ;;
            s) spliceai_snv_prescore=${OPTARG} ;;
            d) spliceai_indel_prescore=${OPTARG} ;;
            v) conservation_file=${OPTARG} ;;
            f) config_file=${OPTARG} ;;
            *) log "Invalid option: -${OPTARG}" && return 1 ;;
        esac
    done

    # If config file is provided, read values from it if not already set
    if [[ -n "${config_file}" ]] && [[ -f "${config_file}" ]]; then
        [[ -z "${assembly}" ]] && assembly=$(read_yaml ${config_file} "assembly")
        [[ -z "${ref_genome}" ]] && ref_genome=$(read_yaml ${config_file} "ref_genome")
        [[ -z "${vep_cache_dir}" ]] && vep_cache_dir=$(read_yaml ${config_file} "vep_cache_dir")
        [[ -z "${threads}" ]] && threads=$(read_yaml ${config_file} "threads")
        [[ -z "${vep_plugins_dir}" ]] && vep_plugins_dir=$(read_yaml ${config_file} "vep_plugins_dir")
    fi

    # Validate required parameters
    if [[ -z "${input_vcf}" ]]; then
        log "Error: Required parameter -i (input_vcf) is missing"
        return 1
    fi

    # Convert hg19/hg38 to GRCh37/GRCh38
    [[ ${assembly} == "hg19" ]] && assembly="GRCh37"
    [[ ${assembly} == "hg38" ]] && assembly="GRCh38"

    local tmp_tag=$(randomID)
    local tmp_output=${input_vcf/.vcf*/.${tmp_tag}.vep.vcf}
    local output_vcf=${input_vcf/.vcf*/.vep.vcf.gz}

    # Check if VCF is already annotated
    check_vcf_validity ${output_vcf} && \
    check_vcf_infotags ${output_vcf} "CSQ" && \
    [[ ${output_vcf} -nt ${input_vcf} ]] && \
    log "The output vcf ${output_vcf} is already annotated by VEP. Skip this step" && \
    return 0

    if [[ -f ${spliceai_snv_prescore} ]] && [[ -f ${spliceai_indel_prescore} ]]; then
        local spliceai_arg="-plugin SpliceAI,snv=${spliceai_snv_prescore},indel=${spliceai_indel_prescore},cutoff=0.5"
    else
        log "The spliceai prescore files are not found, skip the SpliceAI plugin annotation"
        local spliceai_arg=""
    fi

    if [[ -n ${conservation_file} ]]; then
        local conservation_arg="-plugin Conservation,${conservation_file},MAX"
    else
        local conservation_arg=""
    fi

    log "Running basic VEP annotation with the command below:"
    log "vep -i ${input_vcf} --format vcf --vcf --species homo_sapiens --assembly ${assembly} --cache --offline --merged --hgvs --symbol --canonical --numbers --stats_file ${input_vcf/.vcf*/.vep.stats.html} --fork ${threads} --buffer_size 10000 --fasta ${ref_genome} --dir_cache ${vep_cache_dir} ${spliceai_arg} ${conservation_arg} --force_overwrite -o ${output_vcf}"

    # Run VEP annotation with basic options
    vep -i ${input_vcf} \
    --format vcf \
    --vcf \
    --species homo_sapiens \
    --assembly ${assembly} \
    --cache \
    --offline \
    --merged \
    --use_transcript_ref \
    --hgvs \
    --symbol \
    --canonical \
    --total_length \
    --domains \
    --numbers \
    --stats_file ${input_vcf/.vcf*/.vep.stats.html} \
    --fork ${threads} \
    --buffer_size 10000 \
    --fasta ${ref_genome} \
    --dir_cache ${vep_cache_dir} \
    --dir_plugins ${vep_plugins_dir} ${spliceai_arg} ${conservation_arg} \
    --force_overwrite \
    -o ${tmp_output} && \
    bcftools sort -Oz -o ${output_vcf} ${tmp_output} && \
    tabix -f -p vcf ${output_vcf} && \
    announce_remove_tmps ${tmp_output} && \
    display_vcf ${output_vcf}
}



function ClinVar_VCF_deploy() {
    local config_file=${1}
    local CACHEDIR=${2}
    local assembly_version=${3}
    local contig_map=${BASE_DIR}/data/liftover/GRC_to_ucsc.contig.map.tsv

    [[ ! -f ${contig_map} ]] && { log "The contig map file ${contig_map} is not found, please check the file"; return 1; }

    if [[ -z ${assembly_version} ]]; then
        local assembly_version=$(read_yaml ${config_file} "assembly")
    fi

    if [[ ${assembly_version} == "hg19" ]] || [[ ${assembly_version} == "GRCh37" ]]; then
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

    local clinvar_vcf=$(read_yaml ${config_file} "clinvar_vcf")
    local vep_cache_dir=$(read_yaml ${config_file} "vep_cache_dir")
    local vep_plugins_dir=$(read_yaml ${config_file} "vep_plugins_dir")
    local threads=$(read_yaml ${config_file} "threads")
    local ref_genome=$(read_yaml ${config_file} "ref_genome")
    local assembly=$(read_yaml ${config_file} "assembly")
    local spliceai_snv_prescore=$(read_yaml ${config_file} "spliceai_snv_prescore")
    local spliceai_indel_prescore=$(read_yaml ${config_file} "spliceai_indel_prescore")

    check_vcf_validity ${clinvar_vcf} && \
    check_vcf_infotags ${clinvar_vcf} "CSQ" && \
    check_vcf_infotags ${clinvar_vcf} "AF_grpmax_joint" && \
    log "The ClinVar VCF file ${clinvar_vcf} is already downloaded and annotated by VEP. Skip the following steps" && return 0

    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${assembly_version}/clinvar.vcf.gz -O ${CACHEDIR}/clinvar.${assembly_version}.vcf.gz && \
    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${assembly_version}/clinvar.vcf.gz.tbi -O ${CACHEDIR}/clinvar.${assembly_version}.vcf.gz.tbi && \
    log "The ClinVar VCF file is downloaded to ${CACHEDIR}/clinvar.${assembly_version}.vcf.gz" && \
    liftover_from_ucsc_to_GRCh \
    ${CACHEDIR}/clinvar.${assembly_version}.vcf.gz \
    ${contig_map} \
    ${CACHEDIR}/clinvar.${ucsc_assembly_version}.vcf.gz && \
    check_vcf_validity ${CACHEDIR}/clinvar.${ucsc_assembly_version}.vcf.gz && \
    basic_vep_annotation \
    -i ${CACHEDIR}/clinvar.${ucsc_assembly_version}.vcf.gz \
    -a ${assembly} \
    -r ${ref_genome} \
    -c ${vep_cache_dir} \
    -t ${threads} \
    -p ${vep_plugins_dir} \
    -s ${spliceai_snv_prescore} \
    -d ${spliceai_indel_prescore} && \
    check_vcf_validity ${CACHEDIR}/clinvar.${ucsc_assembly_version}.vep.vcf.gz && \
    local gnomad_vcf_chrX=$(read_yaml ${config_file} "gnomad_vcf_chrX") && \
    bash ${SCRIPT_DIR}/annotation_vcf.sh \
    anno_agg_gnomAD_data \
    ${CACHEDIR}/clinvar.${ucsc_assembly_version}.vep.vcf.gz \
    ${threads} \
    ${assembly} \
    ${gnomad_vcf_chrX} && \
    update_yaml ${config_file} "clinvar_vcf" "${CACHEDIR}/clinvar.${ucsc_assembly_version}.vep.vcf.gz"
}


function ClinVar_Gene_stat () {
    local config_file=${1}

    local clinvar_stat=$(read_yaml ${config_file} "clinvar_gene_stat")
    local clinvar_vcf=$(read_yaml ${config_file} "clinvar_vcf")
    local threads=$(read_yaml ${config_file} "threads")

    [[ -f ${clinvar_stat} ]] && \
    log "The clinvar stat file ${clinvar_stat} is already generated. Skip this step" && return 0

    local clinvar_stat=$(dirname ${SCRIPT_DIR})/data/ClinVar/clinvar_2star_stats.pkl.gz
    local clinvar_stat_py=${SCRIPT_DIR}/clinvar_stat_variants.py

    check_path ${clinvar_vcf} || \
    { log "The clinvar vcf file ${clinvar_vcf} is not found, please check the file"; return 1; }

    log "Running command: python ${clinvar_stat_py} -v ${clinvar_vcf} -o ${clinvar_stat} -t ${threads}"
    python ${clinvar_stat_py} \
    -v ${clinvar_vcf} \
    -o ${clinvar_stat} \
    -t ${threads} && \
    update_yaml ${config_file} "clinvar_gene_stat" "${clinvar_stat}"
}


function ClinVar_AA_stat () {
    local config_file=${1}

    local clinvar_aa_stat=$(read_yaml ${config_file} "clinvar_aa_stat")
    local clinvar_splice_stat=$(read_yaml ${config_file} "clinvar_splice_stat")
    local clinvar_vcf=$(read_yaml ${config_file} "clinvar_vcf")
    local threads=$(read_yaml ${config_file} "threads")
    
    [[ -f ${clinvar_aa_stat} ]] && \
    [[ ${clinvar_aa_stat} -nt ${clinvar_vcf} ]] && \
    [[ -f ${clinvar_splice_stat} ]] && \
    [[ ${clinvar_splice_stat} -nt ${clinvar_vcf} ]] && \
    log "The clinvar AA stat file ${clinvar_aa_stat} and splice stat file ${clinvar_splice_stat} are already generated. Skip this step" && \
    return 0 || \
    log "The clinvar AA stat file ${clinvar_aa_stat} and/or splice stat file ${clinvar_splice_stat} are not ready, start generating them"

    local clinvar_aa_stat_py=${SCRIPT_DIR}/stat_aachange_clinvar.py

    log "Running command: python ${clinvar_aa_stat_py} ${clinvar_vcf} ${clinvar_vcf/.vcf*/.aa_change.stats.pkl}"
    python ${clinvar_aa_stat_py} \
    ${clinvar_vcf} \
    ${clinvar_vcf/.vcf*/.aa_change.stats.pkl} \
    ${clinvar_vcf/.vcf*/.splice_change.stats.pkl} \
    ${threads} && \
    update_yaml ${config_file} "clinvar_aa_stat" "${clinvar_vcf/.vcf*/.aa_change.stats.pkl}" && \
    update_yaml ${config_file} "clinvar_splice_stat" "${clinvar_vcf/.vcf*/.splice_change.stats.pkl}"
}



function ClinVar_patho_AF_stat () {
    local config_file=${1}
    local clinvar_vcf=$(read_yaml ${config_file} "clinvar_vcf")
    local threads=$(read_yaml ${config_file} "threads")
    local assembly=$(read_yaml ${config_file} "assembly")
    local gene_pickle=$(read_yaml ${config_file} "clinvar_patho_af_stat")
    local exon_pickle=$(read_yaml ${config_file} "clinvar_patho_exon_af_stat")

    python ${SCRIPT_DIR}/stat_gene_patho_afs.py \
    ${clinvar_vcf} \
    ${clinvar_vcf/.vcf*/.patho_af_stat.pkl} \
    ${threads} && \
    update_yaml ${config_file} "clinvar_patho_af_stat" "${clinvar_vcf/.vcf*/.patho_af_stat.pkl}" && \
    python ${SCRIPT_DIR}/stat_exon_patho_afs.py \
    --vcf_path ${clinvar_vcf} \
    --output ${clinvar_vcf/.vcf*/.patho_exon_af_stat.pkl} \
    --threads ${threads} && \
    update_yaml ${config_file} "clinvar_patho_exon_af_stat" "${clinvar_vcf/.vcf*/.patho_exon_af_stat.pkl}"
}



function Conservation_install() {
    local config=${1}
    local CACHEDIR=${2}
    local assembly_version=${3}

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

    local existed_file=$(read_yaml "${config}" "conservation_file")

    if [[ ${assembly_version} == "GRCh38" ]]; then
        [[ ! -f ${existed_file} ]] && \
        wget http://ftp.ensembl.org/pub/current_compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw -O ${CACHEDIR}/Conservation/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw && \
        log "The conservation scores for ${assembly_version} assembly version are downloaded to ${CACHEDIR}/Conservation" && \
        update_yaml "${config}" "conservation_file" "${CACHEDIR}/Conservation/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw" || \
        { log "Failed to download the conservation scores for ${assembly_version} assembly version"; return 1; }
    elif [[ ${assembly_version} == "GRCh37" ]]; then
        [[ ! -f ${existed_file} ]] && \
        wget http://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw -O ${CACHEDIR}/Conservation/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw
        log "The conservation scores for ${assembly_version} assembly version are downloaded to ${CACHEDIR}/Conservation" && \
        update_yaml "${config}" "conservation_file" "${CACHEDIR}/Conservation/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw" || \
        { log "Failed to download the conservation scores for ${assembly_version} assembly version"; return 1; }
    else
        log "Currently we only support GRCh37 and GRCh38 assembly versions"
        return 1
    fi
}



function InterPro_parsing () {
    local config=${1}
    # Assembly independent

    local interpro_dir=$(read_yaml "${config}" "interpro_dir")
    local interpro_mapping_pickle=$(read_yaml "${config}" "interpro_mapping_pickle")
    local summary_mapping_py=${SCRIPT_DIR}/protein_domain_mapping.py
    local interpro_xml_gz_url=$(read_yaml "${config}" "interpro_xml_gz_url")
    local am_vep_vcf=$(read_yaml "${config}" "alphamissense_vep_vcf")

    [[ ! -d ${interpro_dir} ]] && { log "The InterPro directory ${interpro_dir} is not found, please check the file"; return 1; }
    [[ -f ${interpro_mapping_pickle} ]] && \
    [[ -f ${interpro_dir}/$(basename ${interpro_xml_gz_url}) ]] && \
    [[ ${interpro_mapping_pickle} -nt ${interpro_dir}/$(basename ${interpro_xml_gz_url}) ]] && \
    [[ ${interpro_mapping_pickle} -nt ${summary_mapping_py} ]] && \
    log "The InterPro mapping pickle file ${interpro_mapping_pickle} is already generated. Skip this step" && return 0 || \
    log "The InterPro mapping pickle file ${interpro_mapping_pickle} is not ready or updated, start generating it"


    log "Running command: wget -c ${interpro_xml_gz_url} -O ${interpro_dir}/$(basename ${interpro_xml_gz_url})"
    wget -c ${interpro_xml_gz_url} -O ${interpro_dir}/$(basename ${interpro_xml_gz_url}) && \
    log "Running command: python ${summary_mapping_py} \
    --xml ${interpro_dir}/$(basename ${interpro_xml_gz_url}) \
    --mapping_output ${interpro_dir}/Interpro_entry_mapping.pkl.gz \
    --vcf ${am_vep_vcf} \
    --report ${interpro_dir}/Interpro_format_report.txt"
    python ${summary_mapping_py} \
    --xml ${interpro_dir}/$(basename ${interpro_xml_gz_url}) \
    --mapping_output ${interpro_dir}/Interpro_entry_mapping.pkl.gz \
    --vcf ${am_vep_vcf} \
    --report ${interpro_dir}/Interpro_format_report.txt && \
    update_yaml "${config}" "interpro_mapping_pickle" "${interpro_dir}/Interpro_entry_mapping.pkl.gz"
}


function SpliceVault_install () {
    local config=${1}
    local PLUGIN_CACHEDIR=${2}
    local git_scripts_dir=${3}


    [[ -z ${git_scripts_dir} ]] && git_scripts_dir=$(dirname $(dirname $(dirname $(read_yaml "${config}" "hg38_hg19_chain"))))/scripts
    [[ ! -d ${PLUGIN_CACHEDIR} ]] && { log "The cache directory ${PLUGIN_CACHEDIR} is not found, please check the file"; return 1; }
    local splicevault_prescore=$(read_yaml "${config}" "splicevault_prescore")
    [[ -f ${splicevault_prescore} ]] && { log "The SpliceVault plugin is already installed at ${splicevault_prescore}"; return 0; } || \
    log "The SpliceVault plugin is not installed by checking the existence of specified file path ${splicevault_prescore} in the config file ${config}, now start installing it"

    local splicevault_url=$(read_yaml "${config}" "splicevault_url")
    local ASSEMBLY_VERSION=$(read_yaml "${config}" "assembly")
    if [[ ${ASSEMBLY_VERSION} =~ "GRCh37" ]] || [[ ${ASSEMBLY_VERSION} =~ "hg19" ]]; then
        local chain_file=$(read_yaml "${config}" "hg38_hg19_chain")
        local ref_genome=$(read_yaml "${config}" "ref_genome")
        local tmp_dir=$(read_yaml "${config}" "tmp_dir")
        [[ ! -f ${ref_genome} ]] && { log "The reference genome file ${ref_genome} is not found, please check the file"; return 1; }
        [[ ! -f ${chain_file} ]] && { log "The chain file ${chain_file} is not found, please check the file"; return 1; }
        [[ ! -d ${tmp_dir} ]] && { log "The temporary directory ${tmp_dir} is not found, please check the file"; return 1; }
        # Currently, we only have the SpliceVault data for GRCh38, so we need to convert it to GRCh37
        wget -c ${splicevault_url} -O ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh38.tsv.gz && \
        gunzip -c ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh38.tsv.gz > ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh38.tsv && \
        python ${git_scripts_dir}/SpliceVault_tsv_vcf_conversion.py \
        --input_tsv ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh38.tsv \
        --output_tsv ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv \
        --chain_file ${chain_file} \
        --reference_fasta ${ref_genome} && \
        head -1 ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv > ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv.header && \
        tail -n +2 ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv > ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv.body && \
        LC_ALL=C sort -k1,1V -k2,2n --parallel=8 --buffer-size=8G -T ${tmp_dir} ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv.body > \
        ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv && \
        cat ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv.header ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv | \
        bgzip > ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv.gz && \
        tabix -s 1 -b 2 -e 2 -c "#" ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv.gz && \
        update_yaml "${config}" "splicevault_prescore" "${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv.gz" && \
        announce_remove_tmps ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv.header ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh37.tsv.body
    elif [[ ${ASSEMBLY_VERSION} =~ "GRCh38" ]] || [[ ${ASSEMBLY_VERSION} =~ "hg38" ]]; then
        wget -c ${splicevault_url} -O ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh38.tsv.gz && \
        wget -c ${splicevault_url}.tbi -O ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh38.tsv.gz.tbi && \
        update_yaml "${config}" "splicevault_prescore" "${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh38.tsv.gz"
    else
        log "Not supported assembly version: ${ASSEMBLY_VERSION}, so skip installing the SpliceVault plugin"
    fi
}



function LoFtee_install() {
    # Temporarily disabled due to inconsistent cache file support across different assemblies
    # Also deprecated due to the redundant functionality based on CADD and SpliceAI
    local PLUGIN_CACHEDIR=${1}
    local config=${2}

    local assembly_version=$(read_yaml "${config}" "assembly")
    local loftee_parent_dir=$(read_yaml "${config}" "loftee_parent_dir")

    [[ ! -d ${loftee_parent_dir} ]] && { log "The parent directory for LOFTEE ${loftee_parent_dir} is not found, please check the directory before proceeding"; return 1; }

    local loftee_repo=$(read_yaml "${config}" "loftee_repo")
    local human_ancestor_fasta=$(read_yaml "${config}" "human_ancestor_fasta")
    local loftee_conservation_file=$(read_yaml "${config}" "loftee_conservation_file")
    local gerp_bigwig=$(read_yaml "${config}" "gerp_bigwig")

    if [[ ${assembly_version} == "GRCh37" ]] || [[ ${assembly_version} == "hg19" ]]; then
        [[ -d ${loftee_repo} ]] && [[ -f ${loftee_repo}/LoF.pm ]] && \
        [[ -f ${human_ancestor_fasta} ]] && \
        [[ -f ${loftee_conservation_file} ]] && \
        log "The LOFTEE repository for hg19/GRCh37 is already installed at ${loftee_repo}" && \
        return 0
    fi

    if [[ ${assembly_version} == "GRCh38" ]] || [[ ${assembly_version} == "hg38" ]]; then
        [[ -d ${loftee_repo} ]] && [[ -f ${loftee_repo}/LoF.pm ]] && \
        [[ -f ${human_ancestor_fasta} ]] && \
        [[ -f ${loftee_conservation_file} ]] && \
        [[ -f ${gerp_bigwig} ]] && \
        log "The LOFTEE repository for hg38/GRCh38 is already installed at ${loftee_repo}" && \
        return 0
    fi


    log "You need to download the prescores yourself to ${PLUGIN_CACHEDIR}. And the detailed info can be found in the corresponding github repo https://github.com/konradjk/loftee"
    log "Example running command: perl variant_effect_predictor.pl [--other options to VEP] --plugin LoF,loftee_path:/path/to/loftee,human_ancestor_fa:/path/to/human_ancestor.fa.gz"
    if [[ ${assembly_version} == "GRCh37" ]] || [[ ${assembly_version} == "hg19" ]]; then
        if [[ ! -d ${loftee_parent_dir}/loftee-hg19 ]]; then
            git clone https://github.com/konradjk/loftee.git ${loftee_parent_dir}/loftee-hg19 && \
            update_yaml "${config}" "loftee_repo" "${loftee_parent_dir}/loftee-hg19" && \
            conda env config vars set PERL5LIB="${loftee_parent_dir}/loftee-hg19/:$PERL5LIB" && \
            log "The LOFTEE repository for hg19/GRCh37 is installed at ${loftee_parent_dir}/loftee-hg19 and the PERL5LIB is set to ${PERL5LIB}"
        elif [[ -d ${loftee_parent_dir}/loftee-hg19 ]] && [[ -f ${loftee_parent_dir}/loftee-hg19/.git ]]; then
            cd ${loftee_parent_dir}/loftee-hg19 && \
            git pull && \
            log "The LOFTEE repository for hg19/GRCh37 is updated at ${loftee_parent_dir}/loftee-hg19 and the PERL5LIB is set to ${PERL5LIB}"
        elif [[ -d ${loftee_parent_dir}/loftee-hg19 ]] && [[ -f ${loftee_parent_dir}/loftee-hg19/LoF.pm ]]; then
            log "It seems the LOFTEE repository for hg19/GRCh37 is already installed at ${loftee_parent_dir}/loftee-hg19 but we cant update it by git pull since it is not a git repo"
        else
            log "Failed to install the LOFTEE repository for hg19/GRCh37"
            return 1
        fi
        local human_ancestor_fasta_url="https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz"
        local human_ancestor_fasta_fai_url="https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai"
        local conservation_file_url="https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/phylocsf_gerp.sql.gz"
    elif [[ ${assembly_version} == "GRCh38" ]] || [[ ${assembly_version} == "hg38" ]]; then
        if [[ ! -d ${loftee_parent_dir}/loftee-hg38 ]]; then
            git clone https://github.com/konradjk/loftee.git ${loftee_parent_dir}/loftee-hg38 && \
            cd ${loftee_parent_dir}/loftee-hg38 && \
            git checkout grch38 && \
            git pull && \
            update_yaml "${config}" "loftee_repo" "${loftee_parent_dir}/loftee-hg38" && \
            conda env config vars set PERL5LIB="${loftee_parent_dir}/loftee-hg38/:$PERL5LIB" && \
            log "The LOFTEE repository for hg38/GRCh38 is installed at ${loftee_parent_dir}/loftee-hg38 and the PERL5LIB is set to ${PERL5LIB}"
        elif [[ -d ${loftee_parent_dir}/loftee-hg38 ]] && [[ -f ${loftee_parent_dir}/loftee-hg38/.git ]]; then
            cd ${loftee_parent_dir}/loftee-hg38 && \
            git pull && \
            log "The LOFTEE repository for hg38/GRCh38 is updated at ${loftee_parent_dir}/loftee-hg38 and the PERL5LIB is set to ${PERL5LIB}"
        elif [[ -d ${loftee_parent_dir}/loftee-hg38 ]] && [[ -f ${loftee_parent_dir}/loftee-hg38/LoF.pm ]]; then
            log "It seems the LOFTEE repository for hg38/GRCh38 is already installed at ${loftee_parent_dir}/loftee-hg38 but we cant update it by git pull since it is not a git repo"
        else
            log "Failed to install the LOFTEE repository for hg38/GRCh38"
            return 1
        fi
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
        wget ${human_ancestor_fasta_url} -O ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/human_ancestor.fa.gz && \
        update_yaml "${config}" "human_ancestor_fasta" "${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/human_ancestor.fa.gz"
    fi

    if [[ ! -f ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/human_ancestor.fa.gz.fai ]]; then
        wget ${human_ancestor_fasta_fai_url} -O ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/human_ancestor.fa.gz.fai
    fi

    if [[ ! -f ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/phylocsf_gerp.sql.gz ]]; then
        wget ${conservation_file_url} -O ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/phylocsf_gerp.sql.gz && \
        gunzip -f ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/phylocsf_gerp.sql.gz && \
        update_yaml "${config}" "loftee_conservation_file" "${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/phylocsf_gerp.sql"
    fi

    if [[ ${gerp_bigwig_url} =~ \.bw$ ]] && [[ ! -f ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/gerp_conservation_scores.homo_sapiens.bw ]]; then
        wget ${gerp_bigwig_url} -O ${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/gerp_conservation_scores.homo_sapiens.bw && \
        update_yaml "${config}" "gerp_bigwig" "${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}/gerp_conservation_scores.homo_sapiens.bw"
    elif [[ -z ${gerp_bigwig_url} ]]; then
        log "Specifying the hg19/GRCh37 assembly version. So the bigwig file ${gerp_bigwig_url} is not available" && \
        update_yaml "${config}" "gerp_bigwig" ""
    else
        log "The bigwig file ${gerp_bigwig_url} is already downloaded"
    fi
}


function AlphaMissense_anno_gnomAD() {
    local config_file=${1}
    local alphamissense_vep_vcf=$(read_yaml "${config_file}" "alphamissense_vep_vcf")
    local gnomAD_vcf_chrX=$(read_yaml "${config_file}" "gnomad_vcf_chrX")
    local threads=$(read_yaml "${config_file}" "threads")
    local assembly=$(read_yaml "${config_file}" "assembly")

    [[ -f ${alphamissense_vep_vcf} ]] && [[ -f ${gnomAD_vcf_chrX} ]] || { log "The AlphaMissense VCF file ${alphamissense_vep_vcf} or the gnomAD VCF file ${gnomAD_vcf_chrX} is not found, please check the files"; return 1; }
    
    log "Running anno_agg_gnomAD_data ${alphamissense_vep_vcf} ${threads} ${assembly} ${gnomAD_vcf_chrX} to annotate gnomAD data to AlphaMissense VCF"
    anno_agg_gnomAD_data \
    ${alphamissense_vep_vcf} \
    ${threads} \
    ${assembly} \
    ${gnomAD_vcf_chrX} && \
    update_yaml "${config_file}" "alphamissense_vep_vcf" "${alphamissense_vep_vcf}"
}



# Main installation function
function main_install() {
    local config_file=${1}

    local has_error=0
    # Read configuration
    local conda_env_yaml=$(read_yaml "$config_file" "conda_env_yaml")
    local vep_cache_dir=$(read_yaml "$config_file" "vep_cache_dir")
    local vep_plugins_dir=$(read_yaml "$config_file" "vep_plugins_dir")
    local assembly=$(read_yaml "$config_file" "assembly")
    local ref_fasta=$(read_yaml "$config_file" "ref_genome")
    local vep_plugins_cachedir=$(read_yaml "$config_file" "vep_plugins_cachedir")
    [[ -f ${conda_env_yaml} ]] || { log "The conda env yaml file ${conda_env_yaml} is not found, please check the file"; has_error=1; }
    [[ -d ${vep_cache_dir} ]] || { log "The cache directory for VEP ${vep_cache_dir} is not found, please check the file"; has_error=1; }
    [[ -d ${vep_plugins_dir} ]] || { log "The plugins directory for VEP ${vep_plugins_dir} is not found, please check the file"; has_error=1; }
    [[ -d ${vep_plugins_cachedir} ]] || { log "The plugins cache directory for VEP ${vep_plugins_cachedir} is not found, please check the file"; has_error=1; }
    [[ -f ${ref_fasta} ]] || { log "The reference genome fasta file ${ref_fasta} is not found, please check the file"; has_error=1; }
    [[ ${has_error} -eq 1 ]] && { log "Please check the values in the configuration file ${config_file}"; return 1; }
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
    vep_cache_api_install \
    --VEP_VERSION "$vep_version" \
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

    # Install LOFTEE separately despite it is also a VEP plugin
    LoFtee_install ${vep_plugins_cachedir} ${config_file} || \
    { log "Failed to install LOFTEE"; return 1; }

    # 4. Install gnomAD VCF (basically download bgzipped VCF files)
    local gnomad_vcf_chrX=$(read_yaml "$config_file" "gnomad_vcf_chrX")
    local gnomad_vcf_dir=$(dirname ${gnomad_vcf_chrX})
    if [[ ${assembly} =~ "GRCh37" ]] || [[ ${assembly} =~ "hg19" ]]; then
        # Chain file is small enough to be included in the git repo
        local chain_file=$(read_yaml "$config_file" "chain_file")
    fi
    gnomAD_install ${config_file} ${gnomad_vcf_dir} ${ref_fasta} || \
    { log "Failed to install gnomAD VCF"; return 1; }
    AlphaMissense_anno_gnomAD ${config_file} || \
    { log "Failed to annotate gnomAD data to AlphaMissense VCF"; return 1; }

    # 5. Install InterPro mapping pickle
    InterPro_parsing ${config_file} || \
    { log "Failed to install InterPro mapping pickle"; return 1; }

    # 6. Install ClinVar VCF
    local clinvar_vcf_dir=$(read_yaml "$config_file" "clinvar_vcf_dir")
    ClinVar_VCF_deploy ${config_file} ${clinvar_vcf_dir} ${assembly} || \
    { log "Failed to install ClinVar VCF"; return 1; }
    ClinVar_AA_stat ${config_file} || \
    { log "Failed to stat ClinVar per AA change"; return 1; }
    ClinVar_Gene_stat ${config_file} || \
    { log "Failed to stat ClinVar records per gene"; return 1; }
    AlphaMissense_stat ${config_file} || \
    { log "Failed to generate the AlphaMissense statistics JSON file"; return 1; }
    AlphaMissense_pick_intolerant_domains ${config_file} || \
    { log "Failed to pick intolerant domains from AlphaMissense"; return 1; }
    ClinVar_patho_AF_stat ${config_file} || \
    { log "Failed to stat ClinVar pathogenic variants per allele frequency"; return 1; }


    # 7. Install CADD prescores
    CADD_install \
    ${config_file} \
    ${assembly} || \
    { log "Failed to install CADD prescores"; return 1; }

    log "Congratulations! The installation is completed!"
}



if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
    "$@"
fi
