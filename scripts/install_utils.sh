#!/usr/bin/env bash
SELF_SCRIPT=$(readlink -f "$0")
SELF_DIR=$(dirname ${SELF_SCRIPT})
source ${SELF_DIR}/common_bash_utils.sh

function conda_install_vep() {
    local env_yaml=${1}

	# Extract the env name from the env.yaml file
	local env_name=$(head -1 ${env_yaml} | awk -F ': ' '{print $2;}')
	log "The environment name to-be setup according to the env.yaml ${env_yaml} file is ${env_name}"

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
		mamba activate ${env_name}
		# The to be checked packages are: ensembl-vep, perl-bio-procedural, perl-bioperl
		if [[ $(mamba list | grep ensembl-vep) =~ ensembl-vep ]] && \
		   [[ $(mamba list | grep perl-bio-procedural) =~ perl-bio-procedural ]] && \
		   [[ $(mamba list | grep perl-bioperl) =~ perl-bioperl ]]; then
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
		mamba env create -p ${CONDA_PREFIX} -f ${env_yaml}
	fi
}


function vep_install_wrapper() {
    local VEP_CACHEDIR=""
    local VEP_ASSEMBLY=""
    local VEP_PLUGINS=""
    local VEP_PLUGINSDIR=""
    local VEP_PLUGINSCACHEDIR=""

    local TEMP
    TEMP=$(getopt -o hc:y:r:p: --long help,CACHEDIR:,ASSEMBLY:,PLUGINSDIR:,PLUGINS_CACHEDIR:,PLUGINS:: -- "$@")

    if [[ $? != 0 ]]; then return 1; fi

    eval set -- "$TEMP"

    while true; do
        case "$1" in
            -h|--help)
                echo "Usage: vep_install [options]"
                echo "Options:"
                echo "  -c, --CACHEDIR          Set destination directory for cache files"
                echo "  -y, --ASSEMBLY          Assembly name to use if more than one"
                echo "  -g, --PLUGINS           Comma-separated list of plugins to install"
                echo "  -r, --PLUGINSDIR        Set destination directory for VEP plugins files"
                echo "  -p, --PLUGINS_CACHEDIR  Set destination direcotry for VEP plugins caches"
                return 0
                ;;
            -c|--CACHEDIR)
                VEP_CACHEDIR="$2"
                shift 2
                ;;
            -y|--ASSEMBLY)
                VEP_ASSEMBLY="$2"
                shift 2
                ;;
            -g|--PLUGINS)
                VEP_PLUGINS="$2"
                shift 2
                ;;
            -r|--PLUGINSDIR)
                VEP_PLUGINSDIR="$2"
                shift 2
                ;;
            -p|--PLUGINSCACHEDIR)
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
    [[ $(echo $PERL5LIB) =~ "${VEP_DESTDIR}" ]] && \
    [[ $(echo $PATH) =~ "${VEP_DESTDIR}" ]] && \
    [[ $(echo $PERL5LIB}) =~ "${VEP_PLUGINSDIR}" ]] && \
    { log "The PERL5LIB value and PATH value already include ${VEP_DESTDIR} and ${VEP_PLUGINSDIR}, indicating the following installation process has already been succesfully performed. Skip the function for now."; return 0; }

    local CONDA_ENV_NAME=$(basename $CONDA_PREFIX)
    local VEP_DESTDIR=$(perl -e 'print join("\n", @INC);' | head -1)
    [[ ${VEP_DESTDIR} =~ ${CONDA_ENV_NAME} ]] && log "The dest dir is set to the directory (${VEP_DESTDIR})where perl modules are installed by conda" || \
    { log "Since the function is designed to perform follow up installation of VEP upon the installation of VEP dependencies via conda, we only accept installing VEP API at the perl module installation location previously used by conda"; return 1; }

    # Construct the command
    local cmd="vep_install -d ${VEP_DESTDIR} --AUTO acp -s homo_sapiens --NO_HTSLIB --NO_BIOPERL --CONVERT"
    [[ -n "$VEP_CACHEDIR" ]] && cmd+=" --CACHEDIR $VEP_CACHEDIR"
    [[ -n "$VEP_ASSEMBLY" ]] && cmd+=" --ASSEMBLY $VEP_ASSEMBLY"
    [[ -n "$VEP_PLUGINS" ]] && [[ $VEP_PLUGINS != "empty" ]] && cmd+=" --PLUGINS $VEP_PLUGINS,SpliceAI,UTRAnnotator,AlphaMissense,GO,Phenotypes,Conservation,LOEUF,SameCodon,NMD,PrimateAI" || cmd+=" --PLUGINS SpliceAI,UTRAnnotator,AlphaMissense,GO,Phenotypes,Conservation,LOEUF,SameCodon,NMD,PrimateAI"
    [[ -n "$VEP_PLUGINSDIR" ]] && cmd+=" --PLUGINSDIR $VEP_PLUGINSDIR"

    # Execute the command
    # After executing the command, we need to bind the new PERL5LIB value with the current conda env
    log "Now we start to install the VEP api and downloading caches (which might take a while to finish). So pls be patient"
    $cmd && \
    conda env config vars set PERL5LIB="$VEP_DESTDIR:$VEP_PLUGINSDIR" && \
    conda env config vars set PATH="$VEP_DESTDIR:$PATH" && \
    log "Now we have bound the new PERL5LIB value (adding ${VEP_DESTDIR} and ${VEP_PLUGINSDIR} to the PERL5LIB) with the current conda env ${CONDA_ENV_NAME}" && \
    log "Now we have bound the new PATH value (adding ${VEP_DESTDIR} to the PATH) with the current conda env ${CONDA_ENV_NAME}"
}


function VEP_plugins_install() {
    local VEP_PLUGINSDIR=${1}
    local VEP_PLUGINSCACHEDIR=${2}
    local ASSEMBLY_VERSION=${3}
    local config_file=${4}
    local conda_env_name=${5}


    [[ $CONDA_PREFIX =~ envs ]] && \
    log "currently in a conda env $(basename $CONDA_PREFIX)" || \
    { [[ -z ${conda_env_name} ]] && log "Not encouraged to install the VEP plugins directly in the conda base env. So quit with error." && return 1 || \
      mamba activate ${conda_env_name}; }

    if [[ -z ${conda_env_name} ]]; then
        local conda_env_name=$(basename $CONDA_PREFIX)
    fi


    # Install VEP plugins
    # First install LoFtee
    cd $VEP_PLUGINSDIR && \
    git clone https://github.com/konradjk/loftee.git && \
    conda env config vars set PERL5LIB="$VEP_PLUGINSDIR/loftee/:$PERL5LIB" && \
    local loftee_cache_paths=$(Loftee_install ${VEP_PLUGINSCACHEDIR} ${ASSEMBLY_VERSION}) || \
    { log "Failed to install LoFtee"; return 1; }

    if [[ ${ASSEMBLY_VERSION} =~ "GRCh37" ]] || [[ ${ASSEMBLY_VERSION} =~ "hg19" ]]; then
        # update the human_ancestor_fasta_url and conservation_file_url in the config.yaml file
        update_yaml "$config_file" "loftee_human_ancestor_fasta" "$(echo ${loftee_cache_paths} | head -1)" && \
        update_yaml "$config_file" "loftee_conservation_file" "$(echo ${loftee_cache_paths} | tail -1)"
    elif [[ ${ASSEMBLY_VERSION} =~ "GRCh38" ]] || [[ ${ASSEMBLY_VERSION} =~ "hg38" ]]; then
        # update the gerp_bigwig_url in the config.yaml file
        update_yaml "$config_file" "loftee_human_ancestor_fasta" "$(echo ${loftee_cache_paths} | head -1)" && \
        update_yaml "$config_file" "loftee_conservation_file" "$(echo ${loftee_cache_paths} | tail -2 | head -1)" && \
        update_yaml "$config_file" "loftee_gerp_bigwig" "$(echo ${loftee_cache_paths} | tail -1)"
    else
        log "Not supported assembly version: ${ASSEMBLY_VERSION}"
        return 1
    fi


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
    local loeuf_prescore=$(LOEUF_install ${VEP_PLUGINSCACHEDIR} ${ASSEMBLY_VERSION}) || \
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
        update_yaml "$config_file" "conservation_file_hg38" "${conservation_file}"
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
        read -p "Please specify the absolute path to the indel prescore file, remember the file should correspond to the assembly version of your VEP installation"
        local indel_pre_score_file=${REPLY}
        log "Now we start to install the SpliceAI plugin"
        # We need to return the prescore file paths to the outside of the function
        echo ${snv_pre_score_file}
        echo ${indel_pre_score_file}
    else
        log "Please download the prescores for the SpliceAI plugin first"
        return 1
    fi
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

    read -p "Have you finished formatting the downloaded files? (yes or no)"
    if [[ ${REPLY} =~ "yes" ]] || [[ ${REPLY} =~ "y" ]] || [[ ${REPLY} =~ "Y" ]] || [[ ${REPLY} =~ "Yes" ]] || [[ ${REPLY} =~ "YES" ]]; then
        read -p "Please specify the absolute path to the tsv prescore file, choose the right assembly version"
        local prescore_file=${REPLY}
        # We need to return the prescore file paths to the outside of the function
        echo ${prescore_file}
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

    if [[ -z ${assembly_version} ]]; then
        local assembly_version="hg19"
    fi

    if [[ ! -d ${PLUGIN_CACHEDIR} ]]; then
        mkdir -p ${PLUGIN_CACHEDIR}
    fi

    if [[ ! -d ${PLUGIN_CACHEDIR}/LOEUF ]]; then
        mkdir -p ${PLUGIN_CACHEDIR}/LOEUF
    fi

    log "For details, read ${PLUGIN_DIR}/LOEUF.pm"

    if [[ ! -d ${PLUGIN_CACHEDIR}/supplement ]]; then
        local total_package_url="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7334197/bin/41586_2020_2308_MOESM4_ESM.zip"
        wget ${total_package_url} -O ${PLUGIN_CACHEDIR}/41586_2020_2308_MOESM4_ESM.zip && \
        unzip ${PLUGIN_CACHEDIR}/41586_2020_2308_MOESM4_ESM.zip -d ${PLUGIN_CACHEDIR}
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
            tabix -f -s 76 -b 77 -e 78 ${PLUGIN_CACHEDIR}/LOEUF/loeuf_dataset.tsv.gz
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
    local assembly_version=${2}

    # We need to iterate over all chromosomes to check if the files are already downloaded and valid
    local -a chromosomes=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
    for chr in "${chromosomes[@]}"; do
        log "About to download the gnomAD v4.1 files for chr${chr} from the url https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz"
        if [[ ! -f ${CACHEDIR}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz ]] && \
           check_vcf_validity ${CACHEDIR}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz && \
           [[ ${CACHEDIR}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz.tbi -nt ${CACHEDIR}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz ]]; then
            log "The gnomAD v4.1 files for chr${chr} are already downloaded to ${CACHEDIR} and updated"
        else
            wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz -O ${CACHEDIR}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz && \
            wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz.tbi -O ${CACHEDIR}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz.tbi
        fi
    done

    log "The gnomAD v4.1 files are already downloaded to ${CACHEDIR}, remember that the VCF records are 1-based and currently mapped to hg38 assembly"
}


function gnomAD_liftover_per_chromosome() {
    local hg38_vcf=${1}
    local output_dir=${2}
    local chain_file=${3}
    local hg19_fasta=${4}
    local hg19_vcf_name=$(basename ${hg38_vcf/.vcf*/.hg19.vcf.gz}) # make sure the suffix is .vcf.gz instead of .vcf or .vcf.bgz

    local chrom=$(basename ${hg38_vcf/.vcf*/} | awk -F '.' '{print $NF;}')

    local hg19_vcf=${output_dir}/${hg19_vcf_name}

    crossmap_liftover_hg382hg19 \
    --chain_file ${chain_file} \
    --input_vcf ${hg38_vcf} \
    --output_vcf ${hg19_vcf/.${chrom}/.mixed} \
    --hg19_fasta ${hg19_fasta} && \
    display_vcf ${hg19_vcf} || \
    { log "Failed to liftover the gnomAD v4.1 VCF file from hg38 to hg19 assembly"; return 1; }

    # Now we need to split the VCF file by chromosome as some variants from different chromosomes are mixed in the result hg19 VCF file
    local -a hg19_chrs=($(bcftools query -f '%CHROM\n' ${hg19_vcf} | sort - | uniq - ))

    for chr in "${hg19_chrs[@]}"; do
        if [[ ${chr} == "${chrom}" ]]; then
            bcftools view -Oz -o ${hg19_vcf} ${hg19_vcf/.${chrom}/.mixed} --include 'CHROM="${chr}"'
        else
            bcftools view -Oz -o ${hg19_vcf/.${chrom}/.${chr}} ${hg19_vcf/.${chrom}/.mixed} --include 'CHROM="${chr}"'
        fi
    done
}



function gnomAD_liftover() {
    local hg38_vcf_chr1=${1}
    local chain_file=${2}
    local hg19_fasta=${3}
    local threads=${4}

    if [[ -z ${threads} ]]; then
        local threads=4
    fi

    local basename_hg38_vcf=${hg38_vcf_chr1/.chr1*/}
    local suffix_hg38_vcf=$(basename ${hg38_vcf_chr1/*.chr1./})

    local hg19_dir=$(dirname ${hg38_vcf_chr1})/hg19
    local -a hg38_chrs=( "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY" )

    for chr in "${hg38_chrs[@]}"; do
        mkdir -p ${hg19_dir}/${chr}
    done

    parallel -j ${threads} --dry-run \
    gnomAD_liftover_per_chromosome \
    ${basename_hg38_vcf}.chr{}.${suffix_hg38_vcf} \
    ${hg19_dir}/{} \
    ${hg19_fasta} \
    ${chain_file} ::: "${hg38_chrs[@]}" && \
    parallel -j ${threads} \
    gnomAD_liftover_per_chromosome \
    ${basename_hg38_vcf}.chr{}.${suffix_hg38_vcf} \
    ${hg19_dir}/{} \
    ${hg19_fasta} \
    ${chain_file} ::: "${hg38_chrs[@]}"

    for chr in "${hg38_chrs[@]}"; do
        local tmp_vcfs=($(ls ${hg19_dir}/*/*${chr}*.vcf.gz))
        local hg38_vcf=${hg38_vcf_chr1/.chr1/.${chr}}
        local hg19_vcf=${hg38_vcf/.${chr}/.hg19.${chr}}
        bcftools_concatvcfs -v "${tmp_vcfs[@]}" -o ${hg19_vcf}
    done
}



function ClinVar_VCF_deploy() {
    local CACHEDIR=${1}
    local assembly_version=${2}

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

    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${assembly_version}/clinvar.vcf.gz -O ${CACHEDIR}/${ucsc_assembly_version}/clinvar.vcf.gz && \
    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${assembly_version}/clinvar.vcf.gz.tbi -O ${CACHEDIR}/${ucsc_assembly_version}/clinvar.vcf.gz.tbi
}


function Conservation_install() {
    local CACHEDIR=${1}
    local assembly_version=${2}

    if [[ -z ${assembly_version} ]]; then
        local assembly_version="GRCh38"
    elif [[ ${assembly_version} == "hg38" ]] || [[ ${assembly_version} == "GRCh38" ]]; then
        local assembly_version="GRCh38"
    else
        log "Currently we only support hg38 assemblies, hg19 does not have corresponding GERP conservation bigwig file, you have to use database option to acquire corresponding GERP conservation scores over the internet"
        return 1
    fi

    if [[ ! -d ${CACHEDIR} ]]; then
        mkdir -p ${CACHEDIR}
    fi

    if [[ ! -d ${CACHEDIR}/Conservation ]]; then
        mkdir -p ${CACHEDIR}/Conservation
    fi

    wget http://ftp.ensembl.org/pub/current_compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw -O ${CACHEDIR}/Conservation/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw
    log "The conservation scores for ${assembly_version} assembly version are downloaded to ${CACHEDIR}/Conservation"
    echo ${CACHEDIR}/Conservation/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw
}



function LoFtee_install() {
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
    local vep_plugins_cachedir=$(read_yaml "$config_file" "vep_plugins_cachedir")
    # ... read other config values ...

    # Perform installation steps
    # 1. Install the conda env
    conda_install_vep "$config_file"
    # Test whether currently the conda env is activated
    local conda_env_name=$(conda env list | grep "$conda_env_yaml" | awk '{print $1;}')
    if [[ ${CONDA_PREFIX} =~ ${conda_env_name} ]]; then
        log "The conda env $conda_env_name is already activated"
    else
        conda activate $conda_env_name
    fi

    # 2. Install VEP
    # Update config with installation results
    local vep_version=$(vep --version)
    # Update config with installation results
    update_yaml "$config_file" "vep_installed_version" "$vep_version"

    # Install VEP API and caches and plugins first
    vep_install_wrapper \
    --CACHEDIR "$vep_cache_dir" \
    --PLUGINSDIR "$vep_plugins_dir" \
    --ASSEMBLY "$assembly" \
    --PLUGINSCACHEDIR "$vep_plugins_cachedir" || \
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
        local chain_file=$(read_yaml "$config_file" "chain_file")
    fi

    gnomAD_install ${gnomad_vcf_dir} ${assembly} || \
    { log "Failed to install gnomAD VCF"; return 1; }

    # 5. Install ClinVar VCF
    local clinvar_vcf_dir=$(read_yaml "$config_file" "clinvar_vcf_dir")
    ClinVar_VCF_deploy ${clinvar_vcf_dir} ${assembly} || \
    { log "Failed to install ClinVar VCF"; return 1; }

    # ... other installation steps ...
}



if [[ "${#BASH_SOURCE[@]}" -eq 1 ]]; then
    "$@"
fi
