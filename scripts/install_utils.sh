#! /usr/bin/env bash
SCRIPT_DIR=$(dirname "$(realpath "${BASH_SOURCE[0]}")")
BASE_DIR=$(dirname "${SCRIPT_DIR}")
DATA_DIR=${BASE_DIR}/data
source "${SCRIPT_DIR}/common_bash_utils.sh"
source "${SCRIPT_DIR}/annotation_vcf.sh"
log "The base directory is ${BASE_DIR}, the scripts directory is ${SCRIPT_DIR}"


function priva_find_env_manager() {
    local candidate
    # MAMBA_EXE also finds an initialized Micromamba wrapper outside PATH.
    # Keep its path: a wrapper may intentionally control caches and storage.
    for candidate in "${MAMBA_EXE:-}" "$(type -P micromamba)" \
                     "$(type -P mamba)" "${CONDA_EXE:-}" "$(type -P conda)"; do
        [[ -n ${candidate} && -x ${candidate} ]] || continue
        case "${candidate##*/}" in
            micromamba*) PRIVA_ENV_MANAGER_KIND=micromamba ;;
            mamba*) PRIVA_ENV_MANAGER_KIND=mamba ;;
            conda*) PRIVA_ENV_MANAGER_KIND=conda ;;
            *) continue ;;
        esac
        PRIVA_ENV_MANAGER=${candidate}
        PRIVA_ENV_ACTIVATOR=${candidate}
        PRIVA_ENV_ACTIVATOR_KIND=${PRIVA_ENV_MANAGER_KIND}
        # Older Mamba installations delegate activation/configuration to the
        # Conda executable in the same installation, not an unrelated PATH entry.
        if [[ ${PRIVA_ENV_MANAGER_KIND} == mamba && -x ${candidate%/*}/conda ]]; then
            PRIVA_ENV_ACTIVATOR=${candidate%/*}/conda
            PRIVA_ENV_ACTIVATOR_KIND=conda
        fi
        log "Using ${PRIVA_ENV_MANAGER_KIND}: ${PRIVA_ENV_MANAGER}"
        return 0
    done
    log "ERROR: Install micromamba, mamba, or conda and expose its executable through PATH, MAMBA_EXE, or CONDA_EXE."
    return 1
}


function priva_activate_env() {
    local env_name=${1}
    local activation
    [[ -n ${PRIVA_ENV_MANAGER:-} ]] || priva_find_env_manager || return 1

    # Evaluate activation code in this shell; no interactive alias or profile
    # initialization is required, and the selected wrapper remains in use.
    if [[ ${PRIVA_ENV_ACTIVATOR_KIND} == conda ]]; then
        activation=$("${PRIVA_ENV_ACTIVATOR}" shell.bash activate "${env_name}") || return 1
    else
        activation=$("${PRIVA_ENV_ACTIVATOR}" shell activate --shell bash --name "${env_name}") || return 1
    fi
    eval "${activation}" || return 1
    [[ -n ${CONDA_PREFIX:-} && ${CONDA_PREFIX##*/} == "${env_name}" ]] || {
        log "ERROR: Activation did not select the requested environment ${env_name}."
        return 1
    }
}


function priva_env_vars_set() {
    local assignment
    [[ -n ${CONDA_PREFIX:-} ]] || { log "ERROR: Activate a PriVA environment first."; return 1; }
    [[ -n ${PRIVA_ENV_MANAGER:-} ]] || priva_find_env_manager || return 1
    "${PRIVA_ENV_ACTIVATOR}" env config vars set --prefix "${CONDA_PREFIX}" "$@" || return 1
    # Make the saved settings effective immediately for the remaining steps.
    for assignment in "$@"; do
        export "${assignment}" || return 1
    done
}


function conda_install_vep() {
    # Retain this public function name for existing callers of the installer.
    local -r env_yaml="${BASE_DIR}/acmg_conda_rough.yml"
    local -r env_name=priva_acmg
    local environments packages tool executable
    [[ $# -eq 0 ]] || {
        log "ERROR: conda_install_vep takes no arguments; it uses PriVA's bundled environment."
        return 2
    }
    [[ -r ${env_yaml} ]] || { log "ERROR: Environment file not readable: ${env_yaml}"; return 1; }

    # The repository file and environment name are fixed. No YAML parser is
    # needed until this environment has supplied yq.
    priva_find_env_manager || return 1
    environments=$("${PRIVA_ENV_MANAGER}" env list) || return 1
    if ! awk -v name="${env_name}" '$1 == name { found = 1 } END { exit !found }' <<< "${environments}"; then
        log "Creating environment ${env_name} from ${env_yaml}"
        "${PRIVA_ENV_MANAGER}" env create --name "${env_name}" --file "${env_yaml}" --yes || return 1
    else
        log "Reusing existing environment ${env_name}"
    fi
    priva_activate_env "${env_name}" || return 1

    # yq and samtools must be usable before main_install reads configuration
    # or prepares the reference; do not accidentally use the parent environment.
    for tool in yq samtools vep; do
        executable=$(type -P "${tool}")
        if [[ -z ${executable} || ${executable} != "${CONDA_PREFIX}/"* ]]; then
            log "ERROR: ${env_name} is missing its own ${tool}; update it from ${env_yaml}."
            return 1
        fi
    done
    packages=$("${PRIVA_ENV_MANAGER}" list --name "${env_name}") || return 1
    for tool in ensembl-vep perl-bio-procedural perl-bioperl; do
        if ! awk -v name="${tool}" '$1 == name { found = 1 } END { exit !found }' <<< "${packages}"; then
            log "ERROR: ${env_name} is missing package ${tool}; update it from ${env_yaml}."
            return 1
        fi
    done
    PRIVA_INSTALL_ENV_NAME=${env_name}
    log "Environment ${env_name} is ready; yq is available before configuration parsing."
}


function priva_configure_cache_paths() {
    # Run after environment activation: yq supplies JSON; Python handles paths
    # without evaluating YAML values as shell code. --check makes no changes.
    local config_file=${1}
    local mode=${2:-write}
    local cache_root config_json
    cache_root=$(yq '.cache_root' "${config_file}") || return 1
    [[ -z ${cache_root} || ${cache_root} == null || ${cache_root} == '""' ]] && return 0
    config_json=$(yq '.' "${config_file}") || return 1
    python - "${config_file}" "${config_json}" "${mode}" <<'PRIVA_CACHE_PY'
import json
from pathlib import Path
import posixpath
import re
import sys


DEFAULT_CACHE_FOLDERS = {
    "vep_cache": ("VEP_cache", "vep_cache_dir"),
    "vep_plugin": ("VEP_plugin", "vep_plugins_dir"),
    "vep_plugin_cache": ("VEP_plugin_cache", "vep_plugins_cachedir"),
    "gnomad": ("gnomAD", "gnomad_vcf_dir"),
    "clinvar": ("ClinVar", "clinvar_vcf_dir"),
    "cadd": ("CADD", "cadd_parent_dir"),
}

CACHE_FILE_KEYS = {
    "gnomad_vcf_chrX", "clinvar_vcf", "clinvar_aa_stat",
    "clinvar_splice_stat", "clinvar_patho_af_stat", "clinvar_patho_exon_af_stat",
    "clinvar_gene_stat", "utr_annotator_file", "loftee_parent_dir", "loftee_repo",
    "human_ancestor_fasta", "loftee_conservation_file", "gerp_bigwig",
    "alphamissense_prescore", "alphamissense_vcf", "alphamissense_vep_vcf",
    "alphamissense_pd_stat", "alphamissense_tranx_domain_map",
    "alphamissense_intolerant_domains", "alphamissense_intolerant_motifs",
    "all_intolerant_domains", "pm1_regions_pkl", "conservation_file",
    "splicevault_prescore", "primateai_prescore", "spliceai_snv_prescore",
    "spliceai_indel_prescore", "cds_fasta_file", "mavedb_file",
    "pext_mean_bw", "pext_tissue_bw_dir", "cadd_base_dir",
}


def resolve_cache_paths(config):
    """Resolve six folders and relocate saved paths beneath their old roots."""
    if not isinstance(config, dict):
        raise ValueError("The configuration must be a YAML mapping")
    root = config.get("cache_root")
    if root is None or root == "":
        return {}, {}  # Preserve legacy configurations with full directory paths.
    if not isinstance(root, str) or not re.fullmatch(r"/[A-Za-z0-9_./-]+", root):
        raise ValueError("cache_root must be an absolute path without spaces or shell characters")
    root = posixpath.normpath(root)
    if not root.strip("/"):
        raise ValueError("cache_root must name a dedicated cache directory")

    folders = config.get("cache_folders")
    if folders is None:
        folders = {}
    if not isinstance(folders, dict):
        raise ValueError("cache_folders must be a YAML mapping of resource names to folder names")
    unknown = set(folders) - set(DEFAULT_CACHE_FOLDERS)
    if unknown:
        raise ValueError("Unknown cache_folders entries: " + ", ".join(map(str, unknown)))
    if set(config) & set(DEFAULT_CACHE_FOLDERS):
        raise ValueError("Indent resource folder names beneath cache_folders")

    directories = {}
    previous_roots = []
    for resource, (default_name, key) in DEFAULT_CACHE_FOLDERS.items():
        name = folders.get(resource, default_name)
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", name):
            raise ValueError(f"cache_folders.{resource} must be one folder name, not a path")
        directories[key] = posixpath.join(root, name)
        old = config.get(key)
        if isinstance(old, str) and old.startswith("/"):
            old = posixpath.normpath(old)
            if old != "/":
                previous_roots.append((old, directories[key]))

    # Check the longest old root first, in case plugin folders were nested.
    previous_roots.sort(key=lambda pair: len(pair[0]), reverse=True)
    updates = dict(directories)
    for key in CACHE_FILE_KEYS:
        value = config.get(key)
        if not isinstance(value, str) or not value.startswith("/"):
            continue
        for old, new in previous_roots:
            if value == old or value.startswith(old + "/"):
                updates[key] = new + value[len(old):]
                break
    return directories, updates


def render_cache_config(text, updates):
    """Replace top-level path scalars, retaining unrelated settings/comments."""
    for key, value in sorted(updates.items()):
        pattern = re.compile(r"^(" + re.escape(key) + r"[ \t]*:)[^\n]*$", re.MULTILINE)
        matches = list(pattern.finditer(text))
        if len(matches) > 1:
            raise ValueError(f"Duplicate configuration key: {key}")
        comment = ""
        if matches:
            old_value = matches[0].group(0).split(":", 1)[1].strip()
            scalar = re.fullmatch(r'''(?:"(?:\\.|[^"\\])*"|'(?:''|[^'])*'|[^#]*?)(\s+#.*)?''', old_value)
            if scalar and scalar.group(1):
                comment = scalar.group(1)
        replacement = key + ": " + json.dumps(value) + comment
        if matches:
            text = pattern.sub(lambda match: replacement, text)
        else:
            text = text.rstrip("\n") + "\n" + replacement + "\n"
    return text


def configure_cache_paths(config_file, config, mode):
    if mode not in ("write", "--check"):
        raise ValueError("Use priva_configure_cache_paths CONFIG [--check]")
    directories, updates = resolve_cache_paths(config)
    if mode == "--check":
        print(json.dumps(updates, indent=2, sort_keys=True))
        return
    if not updates:
        return
    path = Path(config_file)
    original = path.read_text()
    updated = render_cache_config(original, updates)
    for directory in directories.values():
        Path(directory).mkdir(parents=True, exist_ok=True)
    if updated != original:
        path.write_text(updated)
    for key, directory in directories.items():
        print(f"{key}: {directory}")


if __name__ == "__main__":
    try:
        configure_cache_paths(sys.argv[1], json.loads(sys.argv[2]), sys.argv[3])
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
PRIVA_CACHE_PY
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
    log "currently in a conda env $(basename "$CONDA_PREFIX")" || \
    { log "Not encouraged to install the VEP dependencies directly in the conda base env. So quit with error."; return 1; }

    [[ ${VEP_ASSEMBLY} =~ "hg19" ]] && VEP_ASSEMBLY="GRCh37"
    [[ ${VEP_ASSEMBLY} =~ "hg38" ]] && VEP_ASSEMBLY="GRCh38"

    # Test the PERL5LIB value and PATH value to see if they already include the VEP_DESTDIR and VEP_PLUGINSDIR
    # If yes for both, then we can directly skip the follow up installation of VEP API
    local VEP_DESTDIR=$(perl -e 'print join("\n", @INC);' | grep site_perl | head -1)
    [[ $(echo $PERL5LIB) =~ "${VEP_DESTDIR}" ]] && \
    [[ $(echo $PATH) =~ "${VEP_DESTDIR}" ]] && \
    [[ $(echo $PERL5LIB) =~ "${VEP_PLUGINSDIR}" ]] && \
    [[ -d ${VEP_CACHEDIR}/homo_sapiens_merged/${VEP_VERSION}_${VEP_ASSEMBLY} ]] && \
    [[ -d ${VEP_CACHEDIR}/homo_sapiens_merged/${VEP_VERSION}_${VEP_ASSEMBLY}/1 ]] && \
    [[ -s ${VEP_CACHEDIR}/homo_sapiens_merged/${VEP_VERSION}_${VEP_ASSEMBLY}/info.txt ]] && \
    [[ -s ${VEP_CACHEDIR}/homo_sapiens_merged/${VEP_VERSION}_${VEP_ASSEMBLY}/1/all_vars.gz ]] && \
    [[ -s ${VEP_CACHEDIR}/homo_sapiens_merged/${VEP_VERSION}_${VEP_ASSEMBLY}/chr_synonyms.txt ]] && \
    [[ -s ${VEP_PLUGINSDIR}/SpliceVault.pm ]] && \
    [[ -s ${VEP_PLUGINSDIR}/SpliceAI.pm ]] && \
    [[ -s ${VEP_PLUGINSDIR}/UTRAnnotator.pm ]] && \
    [[ -s ${VEP_PLUGINSDIR}/AlphaMissense.pm ]] && \
    [[ -s ${VEP_PLUGINSDIR}/GO.pm ]] && \
    [[ -s ${VEP_PLUGINSDIR}/Phenotypes.pm ]] && \
    [[ -s ${VEP_PLUGINSDIR}/Conservation.pm ]] && \
    [[ -s ${VEP_PLUGINSDIR}/LOEUF.pm ]] && \
    [[ -s ${VEP_PLUGINSDIR}/SameCodon.pm ]] && \
    [[ -s ${VEP_PLUGINSDIR}/NMD.pm ]] && \
    [[ -s ${VEP_PLUGINSDIR}/PrimateAI.pm ]] && \
    [[ -s ${VEP_PLUGINSDIR}/MaveDB.pm ]] && \
    { log "The PERL5LIB value and PATH value already include ${VEP_DESTDIR} and ${VEP_PLUGINSDIR}, indicating the following installation process has already been succesfully performed. Skip the function for now."; return 0; }

    local conda_env_name=$(basename "$CONDA_PREFIX")
    [[ ${VEP_DESTDIR} =~ ${conda_env_name} ]] && log "The dest dir is set to the directory (${VEP_DESTDIR}) where perl modules are installed by conda" || \
    { log "Since the function is designed to perform follow up installation of VEP upon the installation of VEP dependencies via conda, we only accept installing VEP API at the perl module installation location previously used by conda"; return 1; }

    # Conda already supplied VEP and its API. Install only cache/plugins, with
    # no newer-release prompt: VEP's update prompt exits zero on batch EOF.
    # Reinstalling the API would also overwrite Conda-managed Perl packages.
    local required_plugins=SpliceAI,SpliceVault,UTRAnnotator,AlphaMissense,GO,Phenotypes,Conservation,LOEUF,SameCodon,NMD,PrimateAI,MaveDB
    local -a cmd=(vep_install -d "$VEP_DESTDIR" --AUTO cp --NO_UPDATE
                  --CACHE_VERSION "$VEP_VERSION" -s homo_sapiens_merged
                  --NO_HTSLIB --NO_BIOPERL --CONVERT)
    [[ -n "$VEP_CACHEDIR" ]] && cmd+=(--CACHEDIR "$VEP_CACHEDIR")
    [[ -n "$VEP_ASSEMBLY" ]] && cmd+=(--ASSEMBLY "$VEP_ASSEMBLY")
    [[ -n "$VEP_PLUGINS" && "$VEP_PLUGINS" != empty ]] && required_plugins="$VEP_PLUGINS,$required_plugins"
    cmd+=(--PLUGINS "$required_plugins")
    [[ -n "$VEP_PLUGINSDIR" ]] && cmd+=(--PLUGINSDIR "$VEP_PLUGINSDIR")

    # Execute the command
    # After executing the command, we need to bind the new PERL5LIB value with the current conda env
    log "Installing cache/plugins for the existing VEP ${VEP_VERSION}: ${cmd[*]}"
    "${cmd[@]}" || { log "Failed to install the VEP API and download caches"; return 1; }
    local expected cache_directory="$VEP_CACHEDIR/homo_sapiens_merged/${VEP_VERSION}_${VEP_ASSEMBLY}"
    for expected in "$cache_directory/info.txt" "$cache_directory/chr_synonyms.txt" \
                    "$cache_directory/1/all_vars.gz"; do
        [[ -s "$expected" ]] || { log "ERROR: VEP installer exited without required cache file: $expected"; return 1; }
    done
    local -a plugin_names
    IFS=, read -r -a plugin_names <<< "$required_plugins"
    for expected in "${plugin_names[@]}"; do
        [[ -s "$VEP_PLUGINSDIR/$expected.pm" ]] || { log "ERROR: VEP installer exited without required plugin: $VEP_PLUGINSDIR/$expected.pm"; return 1; }
    done
    if [[ ":${PERL5LIB:-}:" != *":${VEP_DESTDIR}:"* || ":${PERL5LIB:-}:" != *":${VEP_PLUGINSDIR}:"* ]]; then
        priva_env_vars_set PERL5LIB="$VEP_DESTDIR:$VEP_PLUGINSDIR${PERL5LIB:+:$PERL5LIB}" || return 1
    fi
    if [[ ":${PATH}:" != *":${VEP_DESTDIR}:"* ]]; then
        priva_env_vars_set PATH="$VEP_DESTDIR:$VEP_DESTDIR/htslib:$PATH" || return 1
    fi
    log "VEP environment settings saved and applied to ${conda_env_name}."
}


function VEP_plugins_install() {
    local VEP_PLUGINSDIR=${1}
    local VEP_PLUGINSCACHEDIR=${2}
    local ASSEMBLY_VERSION=${3}
    local config_file=${4}
    local conda_env_name=${5}


    [[ $CONDA_PREFIX =~ "envs" ]] && \
    log "currently in a conda env $(basename "$CONDA_PREFIX")" || \
    { [[ -z ${conda_env_name} ]] && log "Not encouraged to install the VEP plugins directly in the conda base env. So quit with error." && return 1; }

    if [[ -z ${conda_env_name} ]]; then
        local conda_env_name=$(basename "$CONDA_PREFIX")
    fi

    # Install VEP plugins
    # First install UTRAnnotator
    log "The Cache directory to store the requried annotation files for VEP Plugins is ${VEP_PLUGINSCACHEDIR}"
    UTRAnnotator_install ${config_file} ${VEP_PLUGINSCACHEDIR} || \
    { log "Failed to install UTRAnnotator"; return 1; }

    # Then install LOEUF
    LOEUF_install ${config_file} ${ASSEMBLY_VERSION} || \
    { log "Failed to install LOEUF"; return 1; }

    # AlphaMissense_anno consumes Conservation; deploy it before that stage.
    Conservation_install "${config_file}" "${VEP_PLUGINSCACHEDIR}" "${ASSEMBLY_VERSION}" || \
    { log "Failed to install Conservation file"; return 1; }

    # Then install AlphaMissense
    AlphaMissense_install ${config_file} ${VEP_PLUGINSCACHEDIR} ${ASSEMBLY_VERSION} ${VEP_PLUGINSDIR} || \
    { log "Failed to install AlphaMissense"; return 1; }
    # Then annotate the AlphaMissense prescore file to vcf file
    AlphaMissense_anno ${config_file} || \
    { log "Failed to convert the AlphaMissense prescore file to vcf file and annotate protein domains to it"; return 1; }
    AlphaMissense_pick_missense_constrained_segments ${config_file} || \
    { log "Failed to pick missense constrained segments from the AlphaMissense prescore file"; return 1; }

    # Then install SpliceAI
    SpliceAI_install "${config_file}" "${VEP_PLUGINSCACHEDIR}" "${VEP_PLUGINSDIR}" || \
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

    local assembly_version=$(read_yaml "${config_file}" assembly)
    local snv_prescore=$(read_yaml "${config_file}" "spliceai_snv_prescore")
    local indel_prescore=$(read_yaml "${config_file}" "spliceai_indel_prescore")
    case "${assembly_version}" in
        hg19|GRCh37|hg38|GRCh38) ;;
        *) log "Not supported assembly version: ${assembly_version}"; return 1 ;;
    esac
    [[ -f ${SCRIPT_DIR}/spliceai_cache_input.py ]] || {
        log "The SpliceAI cache format helper is missing from this source snapshot."
        return 1
    }

    if [[ -f "${snv_prescore}" && -f "${indel_prescore}" ]]; then
        python "${SCRIPT_DIR}/spliceai_cache_input.py" install --config "${config_file}" \
            --cache-dir "${PLUGIN_CACHEDIR}" --snv "${snv_prescore}" --indel "${indel_prescore}"
        return $?
    fi

    if [[ ${assembly_version} == hg19 || ${assembly_version} == GRCh37 ]]; then
        log "GRCh37 installation uses SpliceAI v1. Download its assembly-matched SNV and indel VCF caches through your authorized Illumina BaseSpace account before continuing."
        log "Publisher data: https://basespace.illumina.com/s/otSPW8hnhaZR"
        log "Index each completed compressed VCF with tabix -p vcf before supplying its absolute path."
        log "For the official VEP plugin input fields, read ${PLUGIN_DIR}/SpliceAI.pm"
    else
        log "GRCh38 installation uses SpliceAI2 v2.0. Complete the publisher-authorized downloads of all native score shards and gene-context resources before continuing."
        log "Publisher data: https://huggingface.co/datasets/illumina-ai/SpliceAI2-data/tree/cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9"
        log "Provide an SNV and an indel input from the complete cache. Decompressed headers determine VCF pass-through or native TSV preparation; the download manifest remains a sidecar."
        log "For the official VEP plugin input fields, read ${PLUGIN_DIR}/SpliceAI.pm"
    fi
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
        read -p "Please specify the absolute path to the indel prescore file, remember the file should correspond to the assembly version of your VEP installation"
        local indel_pre_score_file=${REPLY}
        if [[ ! -f ${indel_pre_score_file} ]]; then
            log "The file ${indel_pre_score_file} does not exist"
            return 1
        fi
        log "Now we start to install the SpliceAI plugin"
        # Validate both staged inputs before atomically updating the original two keys.
        python "${SCRIPT_DIR}/spliceai_cache_input.py" install --config "${config_file}" \
            --cache-dir "${PLUGIN_CACHEDIR}" --snv "${snv_pre_score_file}" --indel "${indel_pre_score_file}" || return 1
    else
        log "Please download the prescores for the SpliceAI plugin first"
        return 1
    fi
}



function SpliceAI2_install() {
    # Preparation only. The official VEP plugin stays in place; scientific
    # criteria/CADD activation is a separate qualification gate.
    [[ $# -eq 3 || $# -eq 4 ]] || {
        log "ERROR: SpliceAI2_install CONFIG MANIFEST OUTPUT_DIR [PRIVA_BIN]"
        return 2
    }
    local config_file=$1 manifest=$2 output_dir=$3 converter=${4:-priva}
    local receipt="${output_dir}/adapter.ready.json" identities config_sha manifest_sha
    identities=$(python - "${config_file}" "${manifest}" "${output_dir}" <<'PY'
import hashlib
import json
import os
from pathlib import Path
import sys
import yaml

config, manifest, output = map(Path, sys.argv[1:])
content = config.read_bytes()
settings = yaml.safe_load(content)
if not isinstance(settings, dict) or settings.get('assembly') not in ('hg38', 'GRCh38'):
    raise SystemExit('SpliceAI2 preparation requires hg38/GRCh38; no fallback is performed.')
if not manifest.is_absolute() or not output.is_absolute():
    raise SystemExit('SpliceAI2 manifest and output directory must be absolute paths.')

def partial(path):
    return any('.partial' in part.lower() for part in Path(path).parts)

if (partial(manifest) or partial(output) or manifest.resolve() != manifest
        or output.resolve() != output or Path(os.path.abspath(output)) != output or not manifest.is_file()):
    raise SystemExit('SpliceAI2 requires a complete regular manifest and a non-partial output path.')
raw = manifest.read_bytes()
data = json.loads(raw)
if (not isinstance(data, dict) or data.get('status') != 'downloads_complete'
        or not isinstance(data.get('confirmed_by'), str) or not data['confirmed_by'].strip()
        or data.get('assembly') != 'GRCh38' or data.get('version') != '2.0'
        or data.get('source_revision') != 'cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9'):
    raise SystemExit('SpliceAI2 requires an owner-confirmed complete GRCh38 v2.0 source manifest at the pinned revision.')

def inspect_paths(value):
    if isinstance(value, dict):
        for key, item in value.items():
            if key == 'path' and (not isinstance(item, str) or partial(item)):
                raise SystemExit('SpliceAI2 does not consume .partial source paths.')
            inspect_paths(item)
    elif isinstance(value, list):
        for item in value:
            inspect_paths(item)

inspect_paths(data)
if output.exists():
    if not output.is_dir() or any(path.name.startswith('.building-') for path in output.iterdir()):
        raise SystemExit('SpliceAI2 output is incomplete or still being prepared; refusing reuse.')
    if not (output / 'adapter.ready.json').exists() and any(output.iterdir()):
        raise SystemExit('SpliceAI2 output has no completed receipt; refusing duplicate conversion.')
print(hashlib.sha256(content).hexdigest(), hashlib.sha256(raw).hexdigest())
PY
    ) || return 1
    read -r config_sha manifest_sha <<< "${identities}"
    if [[ -e ${receipt} || -L ${receipt} ]]; then
        log "Reusing the completed SpliceAI2 receipt after verification; no conversion is started."
    else
        [[ -f ${SCRIPT_DIR}/spliceai2_to_vep_vcf.py ]] || {
            log "ERROR: The SpliceAI2 deployment driver is absent from this source snapshot."
            return 1
        }
        log "Preparing SpliceAI2 GRCh38 from the complete source manifest; no download or v1 fallback is performed."
        python "${SCRIPT_DIR}/spliceai2_to_vep_vcf.py" \
            --manifest "${manifest}" --assembly GRCh38 --output-dir "${output_dir}" \
            --receipt "${receipt}" --priva-bin "${converter}" || return 1
    fi
    python - "${config_file}" "${manifest}" "${output_dir}" "${receipt}" "${config_sha}" "${manifest_sha}" <<'PY'
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import sys
import tempfile
import yaml

config, manifest, output, receipt = map(Path, sys.argv[1:5])
config_sha, manifest_sha = sys.argv[5:7]

def identity(path):
    value = path.stat()
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)

def sha256(path):
    print('SpliceAI2: verify SHA256 ' + str(path), flush=True)
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def no_links_or_partial(path):
    if not path.is_absolute() or Path(os.path.abspath(path)) != path:
        raise SystemExit('SpliceAI2 receipt paths must be absolute and canonical.')
    for part in (path, *path.parents):
        if '.partial' in part.name.lower() or part.is_symlink():
            raise SystemExit('SpliceAI2 refuses symlinks or .partial paths: ' + str(path))

no_links_or_partial(output)
no_links_or_partial(receipt)
if receipt.parent != output or not receipt.is_file() or receipt.stat().st_size == 0:
    raise SystemExit('SpliceAI2 requires a nonempty regular adapter receipt under OUTPUT_DIR.')
lock_path = output / '.install.lock'
no_links_or_partial(lock_path)
lock = lock_path.open('a')
try:
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
except BlockingIOError:
    raise SystemExit('SpliceAI2 preparation is active; refusing receipt reuse.')
if any(path.name.startswith('.building-') for path in output.iterdir()):
    raise SystemExit('SpliceAI2 output is incomplete or still being prepared; refusing reuse.')
tracked = {receipt: identity(receipt), manifest: identity(manifest), config: identity(config)}
data = json.loads(receipt.read_bytes())
if not isinstance(data, dict):
    raise SystemExit('SpliceAI2 adapter receipt must be an object.')
expected = dict(schema_version=1, resource='SpliceAI2', version='2.0', assembly='GRCh38',
                status='prepared_annotation_only', scientific_qualified=False,
                coverage='complete_published',
                source_revision='cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9',
                source_manifest_sha256=manifest_sha)
for key, value in expected.items():
    if data.get(key) != value or type(data.get(key)) is not type(value):
        raise SystemExit('SpliceAI2 receipt does not qualify full cache preparation: ' + key)
if sha256(manifest) != manifest_sha:
    raise SystemExit('SpliceAI2 source manifest changed during preparation.')
original = config.read_bytes()
if hashlib.sha256(original).hexdigest() != config_sha:
    raise SystemExit('PriVA config changed during SpliceAI2 preparation; refusing to overwrite it.')
outputs = data.get('outputs')
if not isinstance(outputs, dict) or set(outputs) != {'snv', 'indel'}:
    raise SystemExit('SpliceAI2 receipt requires both SNV and indel outputs.')
paths = {}
seen = set()
for role in ('snv', 'indel'):
    item = outputs[role]
    if not isinstance(item, dict):
        raise SystemExit('SpliceAI2 output entry must be an object: ' + role)
    for key, checksum_key in (('path', 'sha256'), ('index_path', 'index_sha256')):
        value, checksum = item.get(key), item.get(checksum_key)
        if not isinstance(value, str) or not isinstance(checksum, str) or not re.fullmatch('[0-9a-f]{64}', checksum):
            raise SystemExit('SpliceAI2 output path/checksum is missing or invalid: ' + role)
        path = Path(value)
        no_links_or_partial(path)
        if output not in path.parents or not path.is_file() or path.stat().st_size == 0 or path in seen:
            raise SystemExit('SpliceAI2 output must be a distinct nonempty regular file under OUTPUT_DIR: ' + str(path))
        seen.add(path)
        tracked[path] = identity(path)
        if key == 'path':
            if type(item.get('bytes')) is not int or item['bytes'] != path.stat().st_size:
                raise SystemExit('SpliceAI2 output size mismatch: ' + role)
            if not path.name.endswith(('.vcf.gz', '.vcf.bgz')):
                raise SystemExit('SpliceAI2 output must have a BGZF VCF filename: ' + role)
            paths[role] = str(path)
        elif value != item['path'] + '.tbi':
            raise SystemExit('SpliceAI2 tabix index must match its VCF path: ' + role)
        if sha256(path) != checksum:
            raise SystemExit('SpliceAI2 output SHA256 mismatch: ' + str(path))

updates = {
    'spliceai_snv_prescore': paths['snv'],
    'spliceai_indel_prescore': paths['indel'],
}
settings = yaml.safe_load(original)
text = original.decode()
for key, value in updates.items():
    pattern = re.compile(r'^' + re.escape(key) + r':[^\r\n]*', re.MULTILINE)
    if len(pattern.findall(text)) > 1:
        raise SystemExit('Duplicate SpliceAI2 config key: ' + key)
    replacement = key + ': ' + json.dumps(value)
    if pattern.search(text):
        text = pattern.sub(lambda match: replacement, text)
    else:
        text = text + ('' if not text or text.endswith('\n') else '\n') + replacement + '\n'
modified = yaml.safe_load(text)
if modified != dict(settings, **updates):
    raise SystemExit('SpliceAI2 config update would alter unrelated settings; refusing it.')
for path, before in tracked.items():
    if identity(path) != before:
        raise SystemExit('SpliceAI2 verification input changed: ' + str(path))
# Preserve all existing keys/comments except the original two cache paths.
# Preparation provenance and scientific status remain in the sidecar receipt.
name = None
try:
    with tempfile.NamedTemporaryFile(mode='w', dir=config.parent, prefix=config.name + '.', delete=False) as stream:
        name = stream.name
        stream.write(text)
        stream.flush()
        os.fsync(stream.fileno())
    os.chmod(name, stat.S_IMODE(config.stat().st_mode))
    os.replace(name, config)
finally:
    if name and os.path.exists(name):
        os.unlink(name)
print('SpliceAI2 preparation recorded; scientific activation remains pending.', flush=True)
PY
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
    local cache_root cadd_version reply
    local CADD_script CADD_base_dir CADD_parent_dir=""

    priva_configure_cache_paths "${config_file}" || return 1
    [[ -z ${assembly} ]] && assembly=$(read_yaml "${config_file}" "assembly")
    [[ ${assembly} == "hg19" ]] && assembly="GRCh37"
    [[ ${assembly} == "hg38" ]] && assembly="GRCh38"
    cadd_version=$(read_yaml "${config_file}" "cadd_version")
    cache_root=$(read_yaml "${config_file}" "cache_root")

    if [[ -n ${cache_root} && ${cache_root} != null ]]; then
        # A selected cache root is authoritative; do not prompt to store the
        # large CADD data inside the package-manager environment instead.
        CADD_parent_dir=$(read_yaml "${config_file}" "cadd_parent_dir")
    else
        # Preserve the previous interactive behavior for legacy configurations.
        CADD_script=$(find "${CONDA_PREFIX}/" -type f -name CADD.sh -print -quit)
        [[ -f ${CADD_script} ]] || {
            log "Cannot find CADD.sh in the active environment"
            return 1
        }
        CADD_base_dir=$(dirname "${CADD_script}")
        read -r -p "Store CADD data in ${CADD_base_dir}/data? (yes or no) " reply
        case "${reply}" in
            y|Y|yes|Yes|YES) ;;
            *)
                read -r -p "CADD parent directory (empty uses cadd_parent_dir): " CADD_parent_dir
                [[ -d ${CADD_parent_dir} ]] || CADD_parent_dir=$(read_yaml "${config_file}" "cadd_parent_dir")
                [[ -d ${CADD_parent_dir} ]] || {
                    log "The CADD parent directory does not exist: ${CADD_parent_dir}"
                    return 1
                }
                ;;
        esac
    fi

    if [[ -n ${CADD_parent_dir} ]]; then
        CADD_base_dir=$(find "${CADD_parent_dir}/" -maxdepth 1 -type d -name "CADD-scripts-*${cadd_version#v}*" -print -quit)
        if [[ -z ${CADD_base_dir} ]]; then
            local archive="${CADD_parent_dir}/CADD-scripts.zip"
            local url
            url=$(read_yaml "${config_file}" "cadd_zip_download_url")
            wget "${url}" -O "${archive}" &&
                unzip "${archive}" -d "${CADD_parent_dir}/" &&
                rm "${archive}" || return 1
            CADD_base_dir=$(find "${CADD_parent_dir}/" -maxdepth 1 -type d -name "CADD-scripts-*${cadd_version#v}*" -print -quit)
        fi
    fi

    CADD_script="${CADD_base_dir}/CADD.sh"
    [[ -f ${CADD_script} ]] || {
        log "CADD.sh was not found in the selected CADD installation: ${CADD_base_dir}"
        return 1
    }
    update_or_append_yaml "${config_file}" "cadd_base_dir" "${CADD_base_dir}" || return 1
    local CADD_cache_dir="${CADD_base_dir}/data"
    local CADD_prescore_dir="${CADD_cache_dir}/prescored"
    local CADD_anno_dir="${CADD_cache_dir}/annotations"
    mkdir -p "${CADD_anno_dir}" "${CADD_prescore_dir}/${assembly}_${cadd_version}/incl_anno" || return 1
    log "CADD data directory: ${CADD_cache_dir}"

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
        taken by executing the command offered below to format the downloaded files:

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
        PLUGIN_CACHEDIR="$(read_yaml "${config_file}" "vep_plugins_dir")"
    fi

    local assembly
    assembly="$(read_yaml "${config_file}" "assembly")"
    local mavedb_url
    mavedb_url="$(read_yaml "${config_file}" "mavedb_url")"
    local mavedb_file
    mavedb_file="$(read_yaml "${config_file}" "mavedb_file")"

    local MAVEDB_CACHE_DIR="${PLUGIN_CACHEDIR%/}/MaveDB"
    mkdir -p "${MAVEDB_CACHE_DIR}"

    # Existing plugin TSV download logic
    if [[ ${assembly} == "GRCh37" ]] || [[ ${assembly} == "hg19" ]]; then
        log "The MaveDB file is not available for ${assembly}, skip the downloading process"
        update_yaml "${config_file}" "mavedb_file" ""
        return 0
    fi

    if [[ -f ${mavedb_file} && -f ${mavedb_file}.tbi && ${mavedb_file} -ot ${mavedb_file}.tbi ]]; then
        log "The MaveDB file is already downloaded, skip the downloading process"
    else
        wget -c "${mavedb_url}" -O "${MAVEDB_CACHE_DIR}/MaveDB_variants.${assembly}.tsv.gz" && \
        wget -c "${mavedb_url}.tbi" -O "${MAVEDB_CACHE_DIR}/MaveDB_variants.${assembly}.tsv.gz.tbi" && \
        update_yaml "${config_file}" "mavedb_file" "${MAVEDB_CACHE_DIR}/MaveDB_variants.${assembly}.tsv.gz"
    fi
}




function MaveDB_deployment() {
    # This function is not used for production. Only a record of the commands used to prepare MaveDB metadata
    # Bulk ZIP / csv handling
    # - Expect CSVs under ${MAVEDB_CACHE_DIR}/csv
    # - If count of *scores.csv < 2680, (re)download the zenodo zip and extract into MAVEDB_CACHE_DIR
    local CSV_DIR="${MAVEDB_CACHE_DIR%/}/csv"
    local EXPECTED_COUNT=2680

    local n_scores
    n_scores="$(find ${CSV_DIR} -type f -name '*scores.csv' | wc -l)"
    if [[ ${n_scores} -lt ${EXPECTED_COUNT} ]]; then
        # try to get zenodo url from config
        local ZENODO_URL
        ZENODO_URL="$(read_yaml "${config_file}" "mavedb_zenodo_url" 2>/dev/null || true)"
        if [[ -n "${ZENODO_URL}" ]]; then
            local ZIPPATH="${MAVEDB_CACHE_DIR%/}/mavedb_all.zip"
            log "[mavedb] (re)downloading Zenodo bulk archive to ${ZIPPATH} ..."
            wget -O "${ZIPPATH}" "${ZENODO_URL}"
            log "[mavedb] extracting ${ZIPPATH} into ${MAVEDB_CACHE_DIR} ..."
            unzip -q -o "${ZIPPATH}" -d "${MAVEDB_CACHE_DIR}"
            rm -f "${ZIPPATH}"
        else
            log "[mavedb] mavedb_zenodo_url not found in config; skipping bulk download."
            [[ n_scores -eq 0 ]] && log "[mavedb] No score CSVs found; cannot compute neutral bands." && return 0 || log "[mavedb] proceeding with existing ${n_scores} score CSVs."
        fi
    fi

    log "[mavedb] detected ${n_scores} '*scores.csv' files under ${CSV_DIR}"

    # Run neutral-band python script if CSVs exist
    if [[ ${n_scores} -gt 0 ]]; then
        local SCRIPT_PATH="${SCRIPT_DIR}/compute_neutral_bands.py"
        local BANDS_OUT="${MAVEDB_CACHE_DIR%/}/neutral_bands.tsv"
        log "[mavedb] Running compute_neutral_bands.py --scores-dir \"${CSV_DIR}\" --output \"${BANDS_OUT}\" ..."
        python3 "${SCRIPT_PATH}" --scores-dir "${CSV_DIR}" --output "${BANDS_OUT}"
        log "[mavedb] neutral bands written to ${BANDS_OUT}"
    else
        log "[mavedb] No score CSVs found; skipping neutral-band computation."
    fi

    return 0
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


function AlphaMissense_pick_missense_constrained_segments() {
    local config_file=${1}
    local intolerant_motif_prefix=${2:-}
    local pick_intolerant_motifs_py=${SCRIPT_DIR}/am_pick_missense_constrained_segments.py
    local assembly=$(read_yaml ${config_file} "assembly")
    local alphamissense_intolerant_motifs=$(read_yaml ${config_file} "alphamissense_intolerant_motifs")

    if [[ -z ${intolerant_motif_prefix} ]]; then
        if [[ -n $alphamissense_intolerant_motifs && $alphamissense_intolerant_motifs != null ]]; then
            intolerant_motif_prefix=${alphamissense_intolerant_motifs%.pkl}
        else
            intolerant_motif_prefix="$(read_yaml "$config_file" vep_plugins_cachedir)/AlphaMissense/AlphaMissense_${assembly}.kde"
        fi
    fi

    local alphamissense_vep_vcf=$(read_yaml ${config_file} "alphamissense_vep_vcf")
    local threads=$(read_yaml ${config_file} "threads")
    [[ -f ${alphamissense_intolerant_motifs} ]] && \
    [[ ${alphamissense_intolerant_motifs} -nt ${alphamissense_vep_vcf} ]] && \
    [[ ${alphamissense_intolerant_motifs} -nt ${pick_intolerant_motifs_py} ]] && \
    log "The AlphaMissense intolerant motifs file ${alphamissense_intolerant_motifs} is already generated, skip this step" && \
    return 0 || \
    log "The AlphaMissense intolerant motifs file ${alphamissense_intolerant_motifs} is not generated, we will generate it"

    if [[ -z $alphamissense_intolerant_motifs || $alphamissense_intolerant_motifs == null || ! -d $(dirname "$alphamissense_intolerant_motifs") ]]; then
        alphamissense_intolerant_motifs="${intolerant_motif_prefix}.pkl"
    fi
    mkdir -p "$(dirname "$alphamissense_intolerant_motifs")" || return 1

    python ${pick_intolerant_motifs_py} \
    --vcf_path ${alphamissense_vep_vcf} \
    --rmc_tsv ${SCRIPT_DIR}/../data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.${assembly}.tsv \
    --hcseeker_tsv ${SCRIPT_DIR}/../data/HCSeeker/HC_spots.${assembly}.tsv \
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
    local threads=$(read_yaml ${config_file} "threads")

    [[ -f ${alphamissense_stat} ]] && \
    [[ ${alphamissense_stat} -nt ${alphamissense_vcf} ]] && \
    [[ ${alphamissense_stat} -nt ${stat_py} ]] && \
    log "The AlphaMissense statistics JSON file ${alphamissense_stat} is already generated, skip this step" && \
    return 0 || \
    log "The AlphaMissense statistics JSON file ${alphamissense_stat} is not generated or updated, we will generate it"

    python ${stat_py} \
    ${alphamissense_vcf} \
    ${alphamissense_vcf/.vcf*/.prot.domain.stats.pkl} \
    ${threads:-1} && \
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

    # Post-KS filter inputs
    local alphamissense_vep_vcf=$(read_yaml ${config_file} "alphamissense_vep_vcf")
    local hcseeker_spots_tsv=${DATA_DIR}/HCSeeker/HC_spots.${assembly}.tsv
    local clinvar_vcf=$(read_yaml ${config_file} "clinvar_vcf")

    # Build optional filter flags
    local filter_flags=""
    if [[ -f ${alphamissense_vep_vcf} ]] && [[ -f ${hcseeker_spots_tsv} ]]; then
        filter_flags+=" --am_vcf ${alphamissense_vep_vcf}"
        filter_flags+=" --hcseeker_spots_tsv ${hcseeker_spots_tsv}"
        log "Filter 1 (coldspot overlap) enabled: am_vcf=${alphamissense_vep_vcf}, hcseeker=${hcseeker_spots_tsv}"
    else
        log "Filter 1 (coldspot overlap) skipped: am_vcf or hcseeker_spots_tsv not found"
    fi
    if [[ -f ${clinvar_vcf} ]]; then
        filter_flags+=" --clinvar_vcf ${clinvar_vcf}"
        log "Filter 2 (benign ClinVar) enabled: clinvar_vcf=${clinvar_vcf}"
    else
        log "Filter 2 (benign ClinVar) skipped: clinvar_vcf not found"
    fi

    log "Running DAS intolerant domain analysis: python ${pick_intolerant_domains_py} --pickle_file ${alphamissense_pd_stat} --threads ${threads} --output_dir ${alphamissense_dir}${filter_flags}"
    python ${pick_intolerant_domains_py} \
    --assembly ${assembly} \
    --pickle_file ${alphamissense_pd_stat} \
    --threads ${threads} \
    --output_dir ${alphamissense_dir} \
    ${filter_flags} && \
    update_yaml ${config_file} "alphamissense_intolerant_domains" "${alphamissense_dir}/domain_tolerance_analysis.${assembly}.tsv" && \
    update_yaml ${config_file} "alphamissense_tranx_domain_map" "${alphamissense_dir}/transcript_pp_domain_mapping.${assembly}.pkl" && \
    update_yaml ${config_file} "all_intolerant_domains" "${alphamissense_dir}/all_intolerant_domains.${assembly}.pkl" && \
    log "DAS intolerant domains analysis complete. Outputs saved to ${alphamissense_dir}"
}


function prepare_pm1_regions() {
    local config_file=${1}
    local prepare_pm1_py=${SCRIPT_DIR}/prepare_pm1_regions.py
    local assembly=$(read_yaml ${config_file} "assembly")
    local threads=$(read_yaml ${config_file} "threads")
    local alphamissense_dir=$(dirname $(read_yaml ${config_file} "alphamissense_pd_stat"))
    local pm1_regions_pkl="${alphamissense_dir}/pm1_regions.${assembly}.pkl"
    local all_intolerant_pkl="${alphamissense_dir}/all_intolerant_domains.${assembly}.pkl"

    [[ -f ${pm1_regions_pkl} ]] && \
    [[ ${pm1_regions_pkl} -nt ${prepare_pm1_py} ]] && \
    [[ ${pm1_regions_pkl} -nt ${all_intolerant_pkl} ]] && \
    log "The PM1 regions file ${pm1_regions_pkl} is already generated, skip this step" && \
    return 0 || \
    log "Preparing PM1 regions (DAS + HCSeeker + RMC) …"

    python ${prepare_pm1_py} \
    --config ${config_file} \
    --output ${pm1_regions_pkl} \
    --assembly ${assembly} \
    --threads ${threads} \
    --skip_das && \
    update_yaml ${config_file} "pm1_regions_pkl" "${pm1_regions_pkl}" && \
    log "PM1 regions saved to ${pm1_regions_pkl}"
}



function LOEUF_install() {
    # The LOEUF plugin cache files are stored in the github repo because the original url is not accessible anymore
    local config_file=${1}
    local assembly_version=${2}

    local loeuf_hg19_file=${DATA_DIR}/loeuf/loeuf_dataset.tsv.gz
    local loeuf_hg38_file=${DATA_DIR}/loeuf/loeuf_dataset_hg38.tsv.gz

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
    local hg19_fasta=${2:-$(read_yaml "${config_file}" "ref_genome")}
    local assembly=$(read_yaml "${config_file}" "assembly")
    # Historical callers pass either a reference path or an assembly label.
    case "$hg19_fasta" in
        hg19|GRCh37|hg38|GRCh38)
            assembly=$hg19_fasta
            hg19_fasta=$(read_yaml "${config_file}" "ref_genome") ;;
    esac
    case "$assembly" in
        hg19|GRCh37)
            [[ -s "$hg19_fasta" ]] || { log "ERROR: GRCh37 reference is missing: $hg19_fasta"; return 1; } ;;
        hg38|GRCh38) ;;
        *) log "ERROR: Unsupported gnomAD target assembly: $assembly"; return 1 ;;
    esac

    # We need to iterate over all chromosomes to check if the files are already downloaded and valid
    local -a chromosomes=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
    local gnomAD_vcf_dir=$(read_yaml ${config_file} "gnomad_vcf_dir")

    for chr in "${chromosomes[@]}"; do
        log "About to download the gnomAD v4.1 files for chr${chr} from the url https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz"
        if check_vcf_validity ${gnomAD_vcf_dir}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz; then
            log "The gnomAD v4.1 file for chr${chr} is already downloaded to ${gnomAD_vcf_dir} and updated"
            if [[ ${gnomAD_vcf_dir}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz.tbi -ot ${gnomAD_vcf_dir}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz ]]; then
                log "The index file is not updated, now we start to index the file"
                tabix -f -p vcf ${gnomAD_vcf_dir}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz
            fi
        else
            log "The gnomAD v4.1 file for chr${chr} (${gnomAD_vcf_dir}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz) is not downloaded or not updated, now we start to download the files"
            wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz -O ${gnomAD_vcf_dir}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz && \
            wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz.tbi -O ${gnomAD_vcf_dir}/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz.tbi || \
            { log "Failed to download gnomAD chromosome ${chr}"; return 1; }
        fi
    done

    log "The gnomAD v4.1 files are already downloaded to ${gnomAD_vcf_dir}, remember that the VCF records are 1-based and currently mapped to hg38 assembly"

    if [[ ${assembly} =~ "GRCh37" ]] || [[ ${assembly} =~ "hg19" ]]; then
        local chain_file=${BASE_DIR}/data/liftover/hg38ToHg19.over.chain.gz
        log "Since the gnomAD v4.1 files are currently mapped to hg38 assembly, we need to liftover the files to hg19 assembly, the chain file is ${chain_file}"

        gnomAD_liftover \
        ${gnomAD_vcf_dir}/gnomad.joint.v4.1.sites.chr1.vcf.bgz \
        ${chain_file} \
        "${hg19_fasta}" 4 && \
        log "The gnomAD v4.1 files are liftovered to hg19 assembly and saved to ${gnomAD_vcf_dir}" && \
        parallel -j4 bash ${SCRIPT_DIR}/common_bash_utils.sh check_vcf_validity ${gnomAD_vcf_dir}/gnomad.joint.v4.1.sites.hg19.chr{}.vcf.gz ::: "${chromosomes[@]}" && \
        update_yaml ${config_file} "gnomad_vcf_chrX" "${gnomAD_vcf_dir}/gnomad.joint.v4.1.sites.hg19.chrX.vcf.gz" || \
        { log "Failed to liftover the gnomAD v4.1 files to hg19 assembly"; return 1; }
    else
        update_yaml ${config_file} "gnomad_vcf_chrX" "${gnomAD_vcf_dir}/gnomad.joint.v4.1.sites.chrX.vcf.bgz"
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
    local -a hg38_chrs=( chr{1..22} chrX chrY )

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
    ${chain_file} ::: "${hg38_chrs[@]}" || \
    { log "Failed to liftover one or more gnomAD chromosomes; refusing partial concatenation"; return 1; }

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
    local contig_map=${DATA_DIR}/liftover/GRC_to_ucsc.contig.map.txt

    [[ ! -f ${contig_map} ]] && { log "The contig map file ${contig_map} is not found, please check the file"; return 1; }

    local CACHEDIR=$(read_yaml ${config_file} "clinvar_vcf_dir")
    local assembly_version=$(read_yaml ${config_file} "assembly")

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

    local clinvar_vcf=${CACHEDIR}/$(basename $(read_yaml ${config_file} "clinvar_vcf"))
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

    local clinvar_stat=${DATA_DIR}/ClinVar/$(basename $(read_yaml ${config_file} "clinvar_gene_stat"))
    local clinvar_vcf=$(read_yaml ${config_file} "clinvar_vcf")
    local threads=$(read_yaml ${config_file} "threads")

    [[ -f ${clinvar_stat} ]] && \
    log "The clinvar stat file ${clinvar_stat} is already generated. Skip this step" && return 0

    local clinvar_stat=${DATA_DIR}/ClinVar/clinvar_2star_stats.pkl.gz
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

    local existed_file=$(read_yaml "${config}" "conservation_file")
    if [[ -s "$existed_file" ]]; then
        log "Reusing the configured Conservation file: $existed_file"
        return 0
    fi
    local target="${CACHEDIR}/Conservation/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw"
    local url
    if [[ ${assembly_version} == "GRCh38" ]]; then
        url="https://ftp.ensembl.org/pub/current_compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.${assembly_version}.bw"
    else
        url="https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw"
    fi
    mkdir -p "${CACHEDIR}/Conservation" || return 1
    wget -c "$url" -O "${target}.partial" || \
        { log "Failed to download Conservation for ${assembly_version}"; return 1; }
    [[ -s "${target}.partial" ]] || { log "ERROR: Empty Conservation download"; return 1; }
    mv "${target}.partial" "$target" || return 1
    update_yaml "$config" conservation_file "$target" || return 1
    log "Conservation for ${assembly_version} is deployed to $target"
}



function InterPro_parsing () {
    local config=${1}
    # Assembly independent

    local interpro_dir=${DATA_DIR}/InterPro
    local interpro_mapping_pickle=${interpro_dir}/$(basename $(read_yaml "${config}" "interpro_mapping_pickle"))
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
    update_yaml "${config}" "interpro_mapping_pickle" "${interpro_dir}/Interpro_entry_mapping.pkl.gz" && \
    rm ${interpro_dir}/$(basename ${interpro_xml_gz_url})
}


function SpliceVault_install () {
    local config=${1}
    local PLUGIN_CACHEDIR=${2}
    
    [[ ! -d ${PLUGIN_CACHEDIR} ]] && { log "The cache directory ${PLUGIN_CACHEDIR} is not found, please check the file"; return 1; }
    local splicevault_prescore=$(read_yaml "${config}" "splicevault_prescore")
    [[ -f ${splicevault_prescore} ]] && { log "The SpliceVault plugin is already installed at ${splicevault_prescore}"; return 0; } || \
    log "The SpliceVault plugin is not installed by checking the existence of specified file path ${splicevault_prescore} in the config file ${config}, now start installing it"

    local splicevault_url=$(read_yaml "${config}" "splicevault_url")
    local ASSEMBLY_VERSION=$(read_yaml "${config}" "assembly")
    if [[ ${ASSEMBLY_VERSION} =~ "GRCh37" ]] || [[ ${ASSEMBLY_VERSION} =~ "hg19" ]]; then
        local chain_file=${DATA_DIR}/liftover/$(basename $(read_yaml "${config}" "hg38_hg19_chain"))
        local ref_genome=$(read_yaml "${config}" "ref_genome")
        local tmp_dir=$(read_yaml "${config}" "tmp_dir")
        [[ ! -f ${ref_genome} ]] && { log "The reference genome file ${ref_genome} is not found, please check the file"; return 1; }
        [[ ! -f ${chain_file} ]] && { log "The chain file ${chain_file} is not found, please check the file"; return 1; }
        [[ ! -d ${tmp_dir} ]] && { log "The temporary directory ${tmp_dir} is not found, please check the file"; return 1; }
        # Currently, we only have the SpliceVault data for GRCh38, so we need to convert it to GRCh37
        wget -c ${splicevault_url} -O ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh38.tsv.gz && \
        gunzip -c ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh38.tsv.gz > ${PLUGIN_CACHEDIR}/SpliceVault/SpliceVault_data_GRCh38.tsv && \
        python ${SCRIPT_DIR}/SpliceVault_tsv_vcf_conversion.py \
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
    local vep_plugin_dir=$(read_yaml "${config}" "vep_plugins_dir")

    [[ ! -d ${PLUGIN_CACHEDIR} ]] && { log "The VEP plugin cache directory ${PLUGIN_CACHEDIR} is not found, please check the file"; return 1; }
    [[ ! -d ${vep_plugin_dir} ]] && { log "The VEP plugin directory ${vep_plugin_dir} is not found, please check the file"; return 1; }

    local loftee_parent_dir=${vep_plugin_dir}/LoFtee_repos

    [[ ! -d ${loftee_parent_dir} ]] && mkdir -p ${loftee_parent_dir} || log "The LOFTEE parent directory ${loftee_parent_dir} is prepared"

    local loftee_repo=${loftee_parent_dir}/"loftee-${assembly_version}"

    [[ ! -d ${loftee_cache_dir} ]] && mkdir -p ${loftee_cache_dir} || log "The LOFTEE cache directory ${loftee_cache_dir} is prepared"
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
        if [[ ! -d ${loftee_repo} ]]; then
            git clone https://github.com/konradjk/loftee.git ${loftee_repo} && \
            update_yaml "${config}" "loftee_repo" "${loftee_repo}" && \
            priva_env_vars_set PERL5LIB="${loftee_repo}/:${PERL5LIB:-}" && \
            log "The LOFTEE repository for hg19/GRCh37 is installed at ${loftee_repo} and the PERL5LIB is set to ${PERL5LIB}"
        elif [[ -d ${loftee_repo} ]] && [[ -f ${loftee_repo}/.git ]]; then
            cd ${loftee_repo} && \
            git pull && \
            log "The LOFTEE repository for hg19/GRCh37 is updated at ${loftee_repo} and the PERL5LIB is set to ${PERL5LIB}"
        elif [[ -d ${loftee_repo} ]] && [[ -f ${loftee_repo}/LoF.pm ]]; then
            log "It seems the LOFTEE repository for hg19/GRCh37 is already installed at ${loftee_repo} but we cant update it by git pull since it is not a git repo"
        else
            log "Failed to install the LOFTEE repository for hg19/GRCh37"
            return 1
        fi
        local human_ancestor_fasta_url="https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz"
        local human_ancestor_fasta_fai_url="https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai"
        local conservation_file_url="https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/phylocsf_gerp.sql.gz"
    elif [[ ${assembly_version} == "GRCh38" ]] || [[ ${assembly_version} == "hg38" ]]; then
        if [[ ! -d ${loftee_repo} ]]; then
            git clone https://github.com/konradjk/loftee.git ${loftee_repo} && \
            cd ${loftee_repo} && \
            git checkout grch38 && \
            git pull && \
            update_yaml "${config}" "loftee_repo" "${loftee_repo}" && \
            priva_env_vars_set PERL5LIB="${loftee_repo}/:${PERL5LIB:-}" && \
            log "The LOFTEE repository for hg38/GRCh38 is installed at ${loftee_repo} and the PERL5LIB is set to ${PERL5LIB}"
        elif [[ -d ${loftee_repo} ]] && [[ -f ${loftee_repo}/.git ]]; then
            cd ${loftee_repo} && \
            git pull && \
            log "The LOFTEE repository for hg38/GRCh38 is updated at ${loftee_repo} and the PERL5LIB is set to ${PERL5LIB}"
        elif [[ -d ${loftee_repo} ]] && [[ -f ${loftee_repo}/LoF.pm ]]; then
            log "It seems the LOFTEE repository for hg38/GRCh38 is already installed at ${loftee_repo} but we cant update it by git pull since it is not a git repo"
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

    local loftee_cache_dir=${PLUGIN_CACHEDIR}/LoFtee/${assembly_version}
    mkdir -p ${loftee_cache_dir}

    # Check if the files are already downloaded
    if [[ ! -f ${human_ancestor_fasta} ]]; then
        wget ${human_ancestor_fasta_url} -O ${loftee_cache_dir}/human_ancestor.fa.gz && \
        update_yaml "${config}" "human_ancestor_fasta" "${loftee_cache_dir}/human_ancestor.fa.gz"
    fi

    if [[ ! -f ${human_ancestor_fasta}.fai ]]; then
        wget ${human_ancestor_fasta_fai_url} -O ${loftee_cache_dir}/human_ancestor.fa.gz.fai
    fi

    if [[ ! -f ${loftee_conservation_file} ]]; then
        wget ${conservation_file_url} -O ${loftee_cache_dir}/phylocsf_gerp.sql.gz && \
        gunzip -f ${loftee_cache_dir}/phylocsf_gerp.sql.gz && \
        update_yaml "${config}" "loftee_conservation_file" "${loftee_cache_dir}/phylocsf_gerp.sql"
    fi

    if [[ ${gerp_bigwig_url} =~ \.bw$ ]] && [[ ! -f ${gerp_bigwig} ]]; then
        wget ${gerp_bigwig_url} -O ${loftee_cache_dir}/gerp_conservation_scores.homo_sapiens.bw && \
        update_yaml "${config}" "gerp_bigwig" "${loftee_cache_dir}/gerp_conservation_scores.homo_sapiens.bw"
    elif [[ -z ${gerp_bigwig_url} ]]; then
        log "Specified the hg19/GRCh37 assembly version. So the bigwig file ${gerp_bigwig_url} is not available" && \
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


function CDS_FASTA_install() {
    # Download Ensembl CDS FASTA file for alternative start codon detection
    # Per ClinGen SVI guidelines (PMC6185798), for start_lost variants we need to check
    # if there's a downstream in-frame ATG that could serve as alternative translation start
    local config_file=${1}
    local assembly=$(read_yaml "${config_file}" "assembly")
    local vep_plugins_cachedir=$(read_yaml "${config_file}" "vep_plugins_cachedir")
    
    local cds_fasta_dir="${vep_plugins_cachedir}/Ensembl_CDS"
    mkdir -p "${cds_fasta_dir}"
    
    local cds_fasta_file
    local cds_fasta_url
    
    # Set URLs based on assembly version
    # Using "current" URLs that automatically point to the latest Ensembl release
    # Note: Using HTTPS instead of FTP for better reliability
    if [[ ${assembly} == "GRCh38" ]] || [[ ${assembly} == "hg38" ]]; then
        # GRCh38 uses current_fasta directory (always latest release)
        cds_fasta_url="https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"
        cds_fasta_file="${cds_fasta_dir}/Homo_sapiens.GRCh38.cds.all.fa.gz"
    elif [[ ${assembly} == "GRCh37" ]] || [[ ${assembly} == "hg19" ]]; then
        # GRCh37 uses grch37/current directory (frozen annotations at release-75, but current structure)
        cds_fasta_url="https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.cds.all.fa.gz"
        cds_fasta_file="${cds_fasta_dir}/Homo_sapiens.GRCh37.cds.all.fa.gz"
    else
        log "Unsupported assembly version: ${assembly}. Only GRCh37/hg19 and GRCh38/hg38 are supported."
        return 1
    fi
    
    # Check if file already exists and is valid
    if [[ -f "${cds_fasta_file}" ]] && [[ -s "${cds_fasta_file}" ]]; then
        # Verify it's a valid gzipped file
        if gzip -t "${cds_fasta_file}" 2>/dev/null; then
            log "CDS FASTA file already exists and is valid: ${cds_fasta_file}"
            update_yaml "${config_file}" "cds_fasta_file" "${cds_fasta_file}"
            return 0
        else
            log "Existing CDS FASTA file is corrupted, re-downloading..."
            rm -f "${cds_fasta_file}"
        fi
    fi
    
    # Download the file
    log "Downloading Ensembl CDS FASTA for ${assembly} from: ${cds_fasta_url}"
    log "This file is ~50MB and will be used for alternative start codon detection in start_lost variants"
    
    wget -c "${cds_fasta_url}" -O "${cds_fasta_file}" || {
        log "Failed to download CDS FASTA file from ${cds_fasta_url}"
        log "You can manually download it and place it at: ${cds_fasta_file}"
        return 1
    }
    
    # Verify download
    if gzip -t "${cds_fasta_file}" 2>/dev/null; then
        log "CDS FASTA file downloaded and verified: ${cds_fasta_file}"
        update_yaml "${config_file}" "cds_fasta_file" "${cds_fasta_file}"
        
        # Pre-generate cache by running the Python script
        log "Pre-generating CDS sequence cache (this may take a minute)..."
        python "${SCRIPT_DIR}/alternative_start_codon.py" \
            --cds_fasta "${cds_fasta_file}" \
            --transcript "ENST00000000233" 2>/dev/null || \
            log "Cache generation will happen on first use"
        
        return 0
    else
        log "Downloaded file is corrupted: ${cds_fasta_file}"
        rm -f "${cds_fasta_file}"
        return 1
    fi
}


function ClinGen_deploy() {
    local config_file=${1}
    local assembly=$(read_yaml "${config_file}" "assembly")
    local clinvar_vcf=$(read_yaml "${config_file}" "clinvar_vcf")

    # if the local file is at least 1 month old, then download the latest version, otherwise use the local file
    local clingen_url="https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download"
    local clingen_dir=${DATA_DIR}/clingen
    if [[ ! -f ${clingen_dir}/clingen_erepo.tsv ]] || [[ ! $(find ${clingen_dir}/clingen_erepo.tsv -mtime -30) ]]; then
        local download=1
    else
        local download=0
    fi

    if [[ ${download} -eq 1 ]]; then
        wget -c ${clingen_url} -O ${clingen_dir}/clingen_erepo.tsv
    fi

    python ${SCRIPT_DIR}/prepare_clingen_map.py \
    -v ${clinvar_vcf} \
    -t ${clingen_dir}/clingen_erepo.tsv \
    -o ${clingen_dir}/clingen_map.${assembly}.pkl.gz && \
    update_yaml "${config_file}" "clingen_map" "${clingen_dir}/clingen_map.${assembly}.pkl.gz"

    # Download the dosage sensitivity file
    [[ ${assembly} =~ "GRCh38" ]] || [[ ${assembly} =~ "hg38" ]] && local url_assembly=GRCh38 || :
    [[ ${assembly} =~ "GRCh37" ]] || [[ ${assembly} =~ "hg19" ]] && local url_assembly=GRCh37 || :

    if [[ ! -f ${clingen_dir}/gene_dosage_sensitivity.${assembly}.tsv ]] || [[ ! $(find ${clingen_dir}/gene_dosage_sensitivity.${assembly}.tsv -mtime -30) ]]; then
        wget -c http://ftp.clinicalgenome.org/ClinGen_gene_curation_list_${url_assembly}.tsv -O ${clingen_dir}/gene_dosage_sensitivity.${assembly}.tsv
    else
        log "The dosage sensitivity file ${clingen_dir}/gene_dosage_sensitivity.${assembly}.tsv is already downloaded"
    fi

    mawk -F '\t' '$1 !~ /^#/ {print;} $1 ~ /^#Gene/{print;}' ${clingen_dir}/gene_dosage_sensitivity.${assembly}.tsv > ${clingen_dir}/gene_dosage_sensitivity.${assembly}.tsv.clean && \
    mv ${clingen_dir}/gene_dosage_sensitivity.${assembly}.tsv.clean ${clingen_dir}/gene_dosage_sensitivity.${assembly}.tsv && \
    update_yaml "${config_file}" "gene_dosage_sensitivity" "${clingen_dir}/gene_dosage_sensitivity.${assembly}.tsv"
}



# Main installation function
function pext_install() {
    # Downloads pext (proportion expressed across transcripts) BigWig files from UCSC
    # Mean pext file is always downloaded if pext_enabled=true
    # Tissue-specific files are downloaded if pext_tissues is specified
    # Files are stored in assembly-specific subdirectories (hg19/ or hg38/)
    local config_file="$1"
    
    # Read config values
    local pext_enabled=$(read_yaml "$config_file" "pext_enabled")
    local pext_tissues=$(read_yaml "$config_file" "pext_tissues")
    local vep_plugins_cachedir=$(read_yaml "$config_file" "vep_plugins_cachedir")
    local assembly=$(read_yaml "$config_file" "assembly")
    
    # Skip if pext not enabled
    if [[ "$pext_enabled" != "true" ]]; then
        log "pext annotation disabled (pext_enabled != true), skipping download"
        return 0
    fi
    
    # Validate and normalize assembly name
    local ucsc_assembly
    if [[ "$assembly" == "hg19" ]] || [[ "$assembly" == "GRCh37" ]]; then
        ucsc_assembly="hg19"
    elif [[ "$assembly" == "hg38" ]] || [[ "$assembly" == "GRCh38" ]]; then
        ucsc_assembly="hg38"
    else
        log "ERROR: Unsupported assembly for pext: $assembly"
        return 1
    fi
    
    # Get base URL and mean filename for this assembly
    local base_url=$(get_pext_base_url "$ucsc_assembly")
    local mean_filename=$(get_pext_mean_filename "$ucsc_assembly")
    
    # Create assembly-specific pext directory under VEP plugins cache
    # Structure: pext/{hg19,hg38}/mean.bw and pext/{hg19,hg38}/tissues/*.bw
    local pext_base_dir="${vep_plugins_cachedir}/pext"
    local pext_dir="${pext_base_dir}/${ucsc_assembly}"
    mkdir -p "$pext_dir" || {
        log "ERROR: Failed to create pext directory: $pext_dir"
        return 1
    }
    log "pext files for ${ucsc_assembly} will be stored in: $pext_dir"
    
    # Always download mean_proportion/exp_prop_mean BigWig
    local mean_bw="${pext_dir}/${mean_filename}"
    if [[ ! -f "$mean_bw" ]]; then
        log "Downloading mean pext file: ${mean_filename} for ${ucsc_assembly}..."
        wget -q --show-progress -O "$mean_bw" "${base_url}/${mean_filename}" || {
            log "ERROR: Failed to download ${mean_filename} from ${base_url}"
            rm -f "$mean_bw"
            return 1
        }
        log "Downloaded: $mean_bw ($(ls -lh "$mean_bw" | awk '{print $5}'))"
    else
        log "Mean pext file already exists: $mean_bw"
    fi
    
    # Update config with mean_bw path and tissue directory
    update_yaml "$config_file" "pext_mean_bw" "$mean_bw"
    update_yaml "$config_file" "pext_tissue_bw_dir" "$pext_dir"
    
    # Download tissue-specific files if specified
    if [[ -n "$pext_tissues" ]] && [[ "$pext_tissues" != "null" ]]; then
        # Validate tissue names first
        validate_pext_tissues "$pext_tissues" "$ucsc_assembly" || {
            log "ERROR: Invalid tissue names in pext_tissues config, aborting pext installation"
            return 1
        }
        
        # Create tissues subdirectory under assembly-specific folder
        local tissue_bw_dir="${pext_dir}/tissues"
        mkdir -p "$tissue_bw_dir"
        
        # Download each tissue file
        local -a tissues
        mapfile -t tissues < <(get_pext_tissues_array "$pext_tissues")
        
        for tissue in "${tissues[@]}"; do
            local tissue_bw="${tissue_bw_dir}/${tissue}.bw"
            
            if [[ ! -f "$tissue_bw" ]]; then
                log "Downloading tissue pext file: ${tissue}.bw for ${ucsc_assembly}..."
                wget -q --show-progress -O "$tissue_bw" "${base_url}/${tissue}.bw" || {
                    log "ERROR: Failed to download ${tissue}.bw from ${base_url}"
                    rm -f "$tissue_bw"
                    return 1
                }
                log "Downloaded: $tissue_bw ($(ls -lh "$tissue_bw" | awk '{print $5}'))"
            else
                log "Tissue pext file already exists: $tissue_bw"
            fi
        done
        
        log "Downloaded ${#tissues[@]} tissue-specific pext files to $tissue_bw_dir"
    fi
    
    log "pext data installation complete for ${ucsc_assembly}"
    return 0
}


function reference_genome_install() {
    local config_file=${1}
    local assembly=$(read_yaml "${config_file}" "assembly")
    local ref_fasta=$(read_yaml "${config_file}" "ref_genome")
    local ucsc_assembly=""
    local ref_hint=$(echo "${ref_fasta}" | tr '[:upper:]' '[:lower:]')

    if [[ ${assembly} == "hg19" ]] || [[ ${assembly} == "GRCh37" ]]; then
        ucsc_assembly="hg19"
    elif [[ ${assembly} == "hg38" ]] || [[ ${assembly} == "GRCh38" ]]; then
        ucsc_assembly="hg38"
    elif [[ ${ref_hint} =~ hg19 ]] || [[ ${ref_hint} =~ grch37 ]]; then
        ucsc_assembly="hg19"
    elif [[ ${ref_hint} =~ hg38 ]] || [[ ${ref_hint} =~ grch38 ]]; then
        ucsc_assembly="hg38"
    else
        log "ERROR: Unsupported assembly/ref label. assembly=${assembly}, ref_genome=${ref_fasta}. Supported labels: hg19, GRCh37, hg38, GRCh38"
        return 1
    fi

    local ref_fasta_url=""
    if [[ ${ucsc_assembly} == "hg19" ]]; then
        ref_fasta_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
    else
        ref_fasta_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    fi

    if [[ -z ${ref_fasta} ]] || [[ ${ref_fasta} == "null" ]]; then
        ref_fasta="${DATA_DIR}/reference/ucsc.${ucsc_assembly}.fasta"
        update_yaml "${config_file}" "ref_genome" "${ref_fasta}"
        log "ref_genome was empty; set to default path: ${ref_fasta}"
    fi

    local ref_dir=$(dirname "${ref_fasta}")
    mkdir -p "${ref_dir}" || { log "ERROR: Failed to create reference genome directory ${ref_dir}"; return 1; }

    if ! command -v samtools >/dev/null 2>&1; then
        log "ERROR: samtools is required but not found in PATH"
        return 1
    fi

    if [[ ! -f ${ref_fasta} ]]; then
        local downloader=""
        if command -v wget >/dev/null 2>&1; then
            downloader="wget"
        elif command -v curl >/dev/null 2>&1; then
            downloader="curl"
        else
            log "ERROR: Either wget or curl is required to download ${ucsc_assembly} reference FASTA"
            return 1
        fi

        if ! command -v gzip >/dev/null 2>&1; then
            log "ERROR: gzip is required but not found in PATH"
            return 1
        fi
        if ! command -v gunzip >/dev/null 2>&1; then
            log "ERROR: gunzip is required but not found in PATH"
            return 1
        fi

        local tmp_gz="${ref_fasta}.download.$$.gz"
        local tmp_fasta="${ref_fasta}.download.$$"
        rm -f "${tmp_gz}" "${tmp_fasta}"

        log "Downloading UCSC ${ucsc_assembly} reference FASTA from ${ref_fasta_url}"
        if [[ ${downloader} == "wget" ]]; then
            wget -c "${ref_fasta_url}" -O "${tmp_gz}" || {
                rm -f "${tmp_gz}" "${tmp_fasta}"
                log "ERROR: Failed to download reference FASTA from ${ref_fasta_url}"
                return 1
            }
        else
            curl -L "${ref_fasta_url}" -o "${tmp_gz}" || {
                rm -f "${tmp_gz}" "${tmp_fasta}"
                log "ERROR: Failed to download reference FASTA from ${ref_fasta_url}"
                return 1
            }
        fi

        gzip -t "${tmp_gz}" || {
            rm -f "${tmp_gz}" "${tmp_fasta}"
            log "ERROR: Downloaded reference FASTA archive is corrupted: ${tmp_gz}"
            return 1
        }

        gunzip -c "${tmp_gz}" > "${tmp_fasta}" || {
            rm -f "${tmp_gz}" "${tmp_fasta}"
            log "ERROR: Failed to decompress downloaded archive ${tmp_gz}"
            return 1
        }
        [[ -s ${tmp_fasta} ]] || {
            rm -f "${tmp_gz}" "${tmp_fasta}"
            log "ERROR: Decompressed FASTA is empty: ${tmp_fasta}"
            return 1
        }

        mv "${tmp_fasta}" "${ref_fasta}" || {
            rm -f "${tmp_gz}" "${tmp_fasta}"
            log "ERROR: Failed to move FASTA into place at ${ref_fasta}"
            return 1
        }
        rm -f "${tmp_gz}"
        log "Reference FASTA deployed to ${ref_fasta}"
    else
        log "Reference FASTA already exists: ${ref_fasta}"
    fi

    local ref_fai="${ref_fasta}.fai"
    if [[ ! -f ${ref_fai} ]] || [[ ${ref_fasta} -nt ${ref_fai} ]]; then
        log "Building FASTA index with samtools faidx: ${ref_fasta}"
        samtools faidx "${ref_fasta}" || {
            log "ERROR: Failed to create FASTA index ${ref_fai}"
            return 1
        }
    else
        log "Reference FASTA index already up-to-date: ${ref_fai}"
    fi
}


function yaml_value_or_default() {
    local config_file=${1}
    local key=${2}
    local default_value=${3}
    local value
    value=$(read_yaml "${config_file}" "${key}")
    if [[ -z ${value} ]] || [[ ${value} == "null" ]]; then
        printf "%s\n" "${default_value}"
    else
        printf "%s\n" "${value}"
    fi
}


function update_or_append_yaml() {
    local config_file=${1}
    local key=${2}
    local value=${3}
    local current_value

    current_value=$(read_yaml "${config_file}" "${key}")
    if [[ ${current_value} == "${value}" ]]; then
        return 0
    fi

    if grep -Eq "^${key}[[:space:]]*:" "${config_file}"; then
        update_yaml "${config_file}" "${key}" "${value}"
    else
        printf "\n%s: %s\n" "${key}" "${value}" >> "${config_file}"
    fi
}


function validate_gene_pathogenic_mechanism_evidence() {
    local evidence_tsv=${1}

    [[ -s ${evidence_tsv} ]] || { log "ERROR: Missing/empty DDG2P/G2P mechanism evidence TSV: ${evidence_tsv}"; return 1; }

    python - "${evidence_tsv}" <<'PY'
import csv
import sys

evidence_tsv = sys.argv[1]

required = {
    "gene_symbol",
    "source",
    "source_record_id",
    "source_condition_id",
    "mondo_id",
    "disease_scope",
    "priva_scope",
    "scope_review_status",
    "disease_label",
    "inheritance",
    "penetrance_raw",
    "penetrance_hpo_ids",
    "normalized_penetrance",
    "patho_mode_raw",
    "normalized_mechanisms",
    "mechanism_confidence",
    "disease_confidence",
    "pmids",
}
row_count = 0
g2p_rows = 0
g2p_lof_rows = 0
g2p_penetrance_rows = 0
orphadata_rows = 0
orphadata_mondo_rows = 0
unexpected_sources = set()
with open(evidence_tsv, encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    missing = sorted(required - set(reader.fieldnames or []))
    if missing:
        raise SystemExit(f"{evidence_tsv} missing required columns: {', '.join(missing)}")
    for row in reader:
        row_count += 1
        source = row.get("source", "")
        if source not in {"G2P_DDG2P", "Orphadata"}:
            unexpected_sources.add(source or "<blank>")
        if source == "G2P_DDG2P":
            g2p_rows += 1
            normalized = row.get("normalized_mechanisms", "").upper()
            raw = row.get("patho_mode_raw", "").lower()
            if "LOF" in normalized or "loss of function" in raw:
                g2p_lof_rows += 1
            penetrance = row.get("normalized_penetrance", "").strip().lower()
            if penetrance:
                if penetrance not in {"complete", "incomplete"}:
                    raise SystemExit(
                        f"{evidence_tsv} has unsupported penetrance {penetrance!r}"
                    )
                if not (row.get("source_condition_id") or row.get("mondo_id")):
                    raise SystemExit(
                        f"{evidence_tsv} has penetrance without a condition identifier"
                    )
                g2p_penetrance_rows += 1
        if source == "Orphadata":
            orphadata_rows += 1
            if row.get("mondo_id"):
                orphadata_mondo_rows += 1
if row_count == 0:
    raise SystemExit(f"{evidence_tsv} contains no evidence rows")
if unexpected_sources:
    raise SystemExit(
        f"{evidence_tsv} contains non-condition sources: "
        + ", ".join(sorted(unexpected_sources))
    )
if g2p_rows == 0:
    raise SystemExit(f"{evidence_tsv} contains no G2P_DDG2P rows")
if g2p_lof_rows == 0:
    raise SystemExit(f"{evidence_tsv} contains no DDG2P/G2P LoF rows")
if g2p_penetrance_rows == 0:
    raise SystemExit(f"{evidence_tsv} contains no condition-linked G2P penetrance rows")
if orphadata_rows == 0:
    raise SystemExit(f"{evidence_tsv} contains no explicit Orphadata mechanism rows")
if orphadata_mondo_rows == 0:
    raise SystemExit(f"{evidence_tsv} contains no MONDO-mapped Orphadata rows")
print(
    f"evidence_rows={row_count}; g2p_rows={g2p_rows}; "
    f"g2p_lof_rows={g2p_lof_rows}; "
    f"g2p_penetrance_rows={g2p_penetrance_rows}; "
    f"orphadata_rows={orphadata_rows}; "
    f"orphadata_mondo_rows={orphadata_mondo_rows}"
)
PY
}


function gene_pathogenic_mechanism_cache_install() {
    local config_file=${1}
    local cache_dir
    local raw_dir
    local builder_script
    local evidence_tsv
    local disease_scope_registry

    cache_dir=$(yaml_value_or_default "${config_file}" "gene_mechanism_cache_dir" "${DATA_DIR}/patho_mechanism")
    raw_dir=$(yaml_value_or_default "${config_file}" "gene_mechanism_raw_dir" "${cache_dir}/raw")
    builder_script=$(yaml_value_or_default "${config_file}" "gene_mechanism_builder_script" "${SCRIPT_DIR}/build_gene_pathogenic_mechanism_cache.py")
    evidence_tsv=$(yaml_value_or_default "${config_file}" "ddg2p_mechanism_evidence" "${cache_dir}/gene_pathogenic_mechanism_evidence.tsv")
    disease_scope_registry=$(yaml_value_or_default "${config_file}" "mondo_disease_scope_registry" "${DATA_DIR}/mondo/disease_scope.tsv.gz")

    mkdir -p "${cache_dir}" "${raw_dir}" "$(dirname "${evidence_tsv}")" || {
        log "ERROR: Failed to create gene mechanism cache directories"
        return 1
    }
    [[ -f ${builder_script} ]] || { log "ERROR: Gene mechanism builder not found: ${builder_script}"; return 1; }

    local -a builder_args=(
        "${builder_script}"
        --cache-dir "${cache_dir}"
        --output-tsv "${evidence_tsv}"
        --shared-raw-dir "${raw_dir}"
        --disease-scope-registry "${disease_scope_registry}"
        --timeout "${GENE_MECHANISM_TIMEOUT:-120}"
        --retries "${GENE_MECHANISM_RETRIES:-3}"
    )
    if [[ "${PRIVA_FORCE_GENE_MECHANISM_CACHE:-0}" == "1" ]] || \
       [[ "${PRIVA_FORCE_ALL_CACHES:-0}" == "1" ]]; then
        builder_args+=(--force)
    fi
    if [[ -n "${PROXY_URL:-}" ]]; then
        builder_args+=(--proxy-url "${PROXY_URL}" --download-tool auto)
    fi

    log "Checking and rebuilding PriVA gene mechanism/DDG2P cache as required"
    python "${builder_args[@]}" || {
        log "ERROR: Failed to build gene mechanism/DDG2P cache"
        return 1
    }

    validate_gene_pathogenic_mechanism_evidence "${evidence_tsv}" || {
        log "ERROR: Gene mechanism/DDG2P cache validation failed"
        return 1
    }
    update_or_append_yaml "${config_file}" "ddg2p_mechanism_evidence" "${evidence_tsv}"
    log "Gene condition-mechanism evidence deployed: ${evidence_tsv}"
}


function validate_gofcards_exact_gof_cache() {
    # Any argument after the path is a "<metadata key>=<expected value>" pair,
    # checked against the build provenance the cache records about itself. That
    # lets a caller ask the question a timestamp cannot answer: was this cache
    # built with the same inputs we are about to build with? A pair whose
    # expected value is empty is skipped, because an unset configuration value
    # says nothing about whether the cache is still fit to use.
    local gofcards_json=${1}
    shift
    [[ -s ${gofcards_json} ]] || { log "ERROR: Missing/empty GoFCards cache: ${gofcards_json}"; return 1; }

    python - "${gofcards_json}" "$@" <<'GOFPY'
import gzip
import json
import re
import sys

path = sys.argv[1]
expected_provenance = sys.argv[2:]
opener = gzip.open if path.endswith(".gz") else open
with opener(path, "rt", encoding="utf-8") as handle:
    cache = json.load(handle)

for key in ("metadata", "genes"):
    if key not in cache:
        raise SystemExit(f"{path} has no top-level {key!r} block")
meta = cache["metadata"]
if meta.get("source") != "GoFCards":
    raise SystemExit(f"{path} metadata does not declare GoFCards as the source")
if meta.get("mechanism") != "GOF":
    raise SystemExit(f"{path} metadata does not declare GOF as the mechanism")

for pair in expected_provenance:
    key, _, want = pair.partition("=")
    if not want:
        continue
    # A dotted key walks into a nested block, so a caller can ask about
    # "clinvar.injector_sha256" as easily as about a value recorded at the top.
    have = meta
    for part in key.split("."):
        have = have.get(part) if isinstance(have, dict) else None
    if str(have) != want:
        raise SystemExit(
            f"{path} was built with {key}={have!r}, but {want!r} is configured now"
        )

# loc_<chrom>:<start>:<ref>-><alt>_grch37
VARIANT_ID = re.compile(r"^loc_[0-9XYMT]+:\d+:[ACGTN-]+->[ACGTN-]+_grch37$")
ELIGIBLE = "eligible"
ALLOWED = {
    ELIGIBLE,
    "quarantined_gene_discordance",
    "quarantined_allele_gene_discordance",
    "quarantined_missing_source_gene",
    "quarantined_reviewed_lof",
    "quarantined_reviewed_mixed",
    "quarantined_reviewed_dominant_negative",
    "quarantined_reviewed_uncertain",
    "quarantined_reviewed_exclude",
    "quarantined_mechanism_review_required",
}
genes = variants = eligible = quarantined = with_protein = with_hg38 = evidence = 0
hgvsc_keys = hgvsp_keys = 0
for symbol, gene in cache["genes"].items():
    genes += 1
    if not symbol.strip():
        raise SystemExit(f"{path} has an empty gene symbol key")
    for vid, variant in (gene.get("variants") or {}).items():
        variants += 1
        if not VARIANT_ID.match(vid):
            raise SystemExit(f"{path}: {symbol}/{vid} is not a well-formed variant identifier")
        record = variant.get("record") or {}
        state = (record.get("eligibility") or {}).get("status", "")
        if state not in ALLOWED:
            raise SystemExit(f"{path}: {symbol}/{vid} has unknown eligibility {state!r}")
        assemblies = variant.get("assemblies") or {}
        if not assemblies:
            raise SystemExit(f"{path}: {symbol}/{vid} is present on no assembly")
        protein = coding = False
        for assembly, block in assemblies.items():
            if assembly not in ("hg19", "hg38"):
                raise SystemExit(f"{path}: {symbol}/{vid} has unknown assembly {assembly!r}")
            coords = block.get("genomic") or {}
            if not coords.get("chrom") or not isinstance(coords.get("pos"), int):
                raise SystemExit(f"{path}: {symbol}/{vid}/{assembly} has no usable coordinates")
            if not coords.get("ref") or not coords.get("alt"):
                raise SystemExit(f"{path}: {symbol}/{vid}/{assembly} has no reference or alternate")
            for transcript, view in (block.get("transcripts") or {}).items():
                by_c = view.get("by_hgvsc") or {}
                by_p = view.get("by_hgvsp") or {}
                hgvsc_keys += len(by_c)
                hgvsp_keys += len(by_p)
                if by_c: coding = True
                if by_p: protein = True
                for hgvsp, coding_keys in by_p.items():
                    if not isinstance(coding_keys, list) or not coding_keys:
                        raise SystemExit(f"{path}: {symbol}/{vid} protein key {hgvsp!r} maps to no coding change")
                    for ck in coding_keys:
                        if ck not in by_c:
                            raise SystemExit(f"{path}: {symbol}/{vid} protein key {hgvsp!r} points at unknown {ck!r}")
        if "hg38" in assemblies: with_hg38 += 1
        if protein: with_protein += 1
        if state == ELIGIBLE:
            eligible += 1
            if not (protein or coding or assemblies):
                raise SystemExit(f"{path}: {symbol}/{vid} is eligible but offers no match route")
        else:
            quarantined += 1
            if not (record.get("eligibility") or {}).get("reason"):
                raise SystemExit(f"{path}: {symbol}/{vid} is quarantined without a reason")
        evidence += len(record.get("evidence") or [])

if variants == 0:
    raise SystemExit(f"{path} contains no variants")
if eligible == 0:
    raise SystemExit(f"{path} contains no runtime-eligible variants")
if with_protein == 0:
    raise SystemExit(f"{path} contains no variant with a protein-change key")
print(
    f"gofcards_genes={genes}; variants={variants}; eligible={eligible}; "
    f"quarantined={quarantined}; with_protein_key={with_protein}; on_hg38={with_hg38}; "
    f"hgvsc_keys={hgvsc_keys}; hgvsp_keys={hgvsp_keys}; evidence_entries={evidence}"
)
GOFPY
}


function gofcards_clinvar_injection_install() {
    # Attach ClinVar VCV condition evidence to each normalized GoFCards variant.
    # This is its own step so a ClinVar refresh does not force the VEP pipeline
    # to rerun. See docs/GOFCARDS_CACHE_SCHEMA.md for the injected shape.
    local config_file=${1}
    [[ -f ${config_file} ]] || { log "ERROR: Config file not found: ${config_file}"; return 1; }

    local injector
    local source_cache
    local target_json
    local clinvar_xml
    local min_stars
    local stats_json
    local audit_tsv
    local workdir
    local lag_days

    injector=$(yaml_value_or_default "${config_file}" "gofcards_clinvar_injector_script" "${SCRIPT_DIR}/inject_clinvar_into_gofcards.py")
    workdir=$(yaml_value_or_default "${config_file}" "gofcards_workdir" "${DATA_DIR}/gofcards/build_work")
    # Reads the build intermediate and writes the one deployed GoFCards cache.
    # Everything downstream reads that single file, so there is never a question
    # of which of two caches a consumer meant.
    source_cache="${workdir}/gofcards_exact_gof.normalized.json.gz"
    target_json=$(yaml_value_or_default "${config_file}" "gofcards_exact_gof_cache" "${DATA_DIR}/gofcards/gofcards_exact_gof.json.gz")
    clinvar_xml=$(yaml_value_or_default "${config_file}" "clinvar_vcv_xml" "${DATA_DIR}/patho_mechanism/raw/clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz")
    min_stars=$(yaml_value_or_default "${config_file}" "clinvar_min_review_stars" "2")
    lag_days=$(yaml_value_or_default "${config_file}" "gofcards_clinvar_reinjection_lag_days" "90")
    stats_json="${workdir}/clinvar_injection_stats.json"
    audit_tsv=$(yaml_value_or_default "${config_file}" "clinvar_vcv_match_audit" "${DATA_DIR}/patho_mechanism/clinvar_vcv_gofcards_matches.tsv")

    [[ -f ${injector} ]] || { log "ERROR: ClinVar injector not found: ${injector}"; return 1; }
    [[ -s ${source_cache} ]] || {
        log "ERROR: Normalized GoFCards cache not found (run gofcards_exact_gof_cache_install first): ${source_cache}"
        return 1
    }
    [[ -s ${clinvar_xml} ]] || { log "ERROR: ClinVar VCV XML not found: ${clinvar_xml}"; return 1; }

    mkdir -p "$(dirname "${target_json}")" "${workdir}" || return 1

    # Reading the weekly VCV XML end to end is the most expensive step in the
    # mechanism chain and the installer runs this step on every pass, so the
    # question here is only ever "would rereading it change anything".
    #
    # It is answered from what the cache records about how it was made, not from
    # file timestamps. A timestamp says nothing useful: checking out a branch,
    # copying a tree or touching a file moves it without changing a byte, and
    # each of those would otherwise have cost hours of rereading. A recorded
    # content hash that differs means the injector or its input really did
    # change, which is worth acting on at once -- a fix to the injector should
    # not have to wait out a refresh interval before it reaches the cache.
    #
    # The XML itself is the exception. A new VCV release lands every week and
    # almost never alters the conditions attached here, so the cache is allowed
    # to trail it and is refreshed once the gap reaches
    # gofcards_clinvar_reinjection_lag_days, which is months rather than days.
    local reinject=0
    if [[ "${PRIVA_FORCE_GOFCARDS_CLINVAR_INJECTION:-0}" == "1" ]] || \
       [[ "${PRIVA_FORCE_ALL_CACHES:-0}" == "1" ]] || [[ ! -s ${target_json} ]]; then
        reinject=1
    else
        local injector_sha256
        local source_cache_sha256
        injector_sha256=$(sha256sum "${injector}" | awk '{print $1}')
        source_cache_sha256=$(sha256sum "${source_cache}" | awk '{print $1}')
        if ! validate_gofcards_exact_gof_cache "${target_json}" \
                 "clinvar.injector_sha256=${injector_sha256}" \
                 "clinvar.source_cache_sha256=${source_cache_sha256}" \
                 "clinvar.min_review_stars=${min_stars}" >/dev/null 2>&1; then
            log "Deployed GoFCards cache was not produced by the current injector, input or review-star setting"
            reinject=1
        fi
        if [[ ${reinject} -eq 0 ]]; then
            local lag_seconds=$(( $(stat -c %Y "${clinvar_xml}") - $(stat -c %Y "${target_json}") ))
            if (( lag_seconds > lag_days * 86400 )); then
                log "Deployed GoFCards cache trails the ClinVar VCV XML by $(( lag_seconds / 86400 )) days, past the ${lag_days}-day limit"
                reinject=1
            fi
        fi
    fi

    if [[ ${reinject} -eq 1 ]]; then
        log "Injecting ClinVar VCV conditions into GoFCards variants -> ${target_json}"
        python "${injector}" \
            --gofcards-cache "${source_cache}" \
            --clinvar-vcv "${clinvar_xml}" \
            --out-json "${target_json}" \
            --min-review-stars "${min_stars}" \
            --audit-tsv "${audit_tsv}" \
            --stats-json "${stats_json}" 2>&1 | tee "${workdir}/clinvar_injection.log"

        local status=${PIPESTATUS[0]}
        [[ ${status} -eq 0 ]] || {
            log "ERROR: ClinVar injection failed; see ${workdir}/clinvar_injection.log"
            return 1
        }
    else
        log "GoFCards cache already carries current ClinVar conditions: ${target_json}"
    fi

    validate_gofcards_exact_gof_cache "${target_json}" || {
        log "ERROR: GoFCards cache failed validation after ClinVar injection"
        return 1
    }
    update_or_append_yaml "${config_file}" "gofcards_exact_gof_cache" "${target_json}"
    update_or_append_yaml "${config_file}" "gofcards_clinvar_injector_script" "${injector}"
    log "GoFCards exact GOF cache deployed (with ClinVar conditions): ${target_json}"
}


function gofcards_exact_gof_cache_install() {
    # Build the GoFCards exact GOF evidence cache from the raw GoFCards sources.
    # The whole pipeline is one script; see docs/GOFCARDS_NORMALIZATION_METHODS.md
    # for why the reference check and the liftover are separate gates.
    local config_file=${1}
    [[ -f ${config_file} ]] || { log "ERROR: Config file not found: ${config_file}"; return 1; }

    local builder_script
    local target_tsv
    local workdir
    local hg19_fasta
    local hg38_fasta
    local chain_file
    local vep_cache_dir
    local vep_version
    local hgnc_table
    local review_tsv
    local stats_json
    local refresh_days

    builder_script=$(yaml_value_or_default "${config_file}" "gofcards_builder_script" "${SCRIPT_DIR}/build_gofcards_exact_gof_cache.py")
    workdir=$(yaml_value_or_default "${config_file}" "gofcards_workdir" "${DATA_DIR}/gofcards/build_work")
    refresh_days=$(yaml_value_or_default "${config_file}" "gofcards_refresh_days" "180")
    # A build intermediate, not the deployed cache. The single deployed file is
    # written by gofcards_clinvar_injection_install, which reads this one and
    # adds the ClinVar conditions. Keeping the intermediate in the build
    # directory is what makes the deployed artifact unambiguous: there is one
    # GoFCards cache, and it is the finished one.
    target_tsv="${workdir}/gofcards_exact_gof.normalized.json.gz"

    hg19_fasta=$(read_yaml "${config_file}" "ref_genome")
    hg38_fasta=$(yaml_value_or_default "${config_file}" "gofcards_hg38_fasta" "")
    chain_file=$(yaml_value_or_default "${config_file}" "gofcards_hg19_to_hg38_chain" "${DATA_DIR}/liftOver/hg19ToHg38.over.chain")
    vep_cache_dir=$(read_yaml "${config_file}" "vep_cache_dir")
    vep_version=$(yaml_value_or_default "${config_file}" "vep_installed_version" "")
    hgnc_table=$(yaml_value_or_default "${config_file}" "hgnc_table" "${DATA_DIR}/hgnc/non_alt_loci_set.tsv")
    review_tsv=$(yaml_value_or_default "${config_file}" "gofcards_mechanism_review_tsv" "${DATA_DIR}/gofcards/gofcards_mechanism_reviews.tsv")
    stats_json="${workdir}/gofcards_build_stats.json"

    # Age alone is not a fitness test. The refresh interval exists to pick up new
    # GoFCards releases, and it cannot see a change on this side of the build:
    # edited normalization code, a newly reviewed allele, a new HGNC table, a new
    # chain or FASTA, or a different VEP annotation cache. Polarity review is
    # applied only while building, so without these checks an allele newly marked
    # LOF or mixed-mechanism in the review table keeps its GOF eligibility until
    # the interval lapses or someone forces a rebuild.
    local rebuild=0
    if [[ "${PRIVA_FORCE_GOFCARDS_CACHE:-0}" == "1" ]] || \
       [[ "${PRIVA_FORCE_ALL_CACHES:-0}" == "1" ]] || [[ ! -s ${target_tsv} ]]; then
        rebuild=1
    else
        local source_file
        for source_file in "${builder_script}" "${review_tsv}" "${hgnc_table}" \
                           "${chain_file}" "${hg19_fasta}" "${hg38_fasta}"; do
            if [[ -e ${source_file} ]] && [[ ${source_file} -nt ${target_tsv} ]]; then
                log "GoFCards normalization input is newer than its cache: ${source_file}"
                rebuild=1
                break
            fi
        done
        # The cache records which VEP cache, version and chain produced it, which
        # is the only way to notice an annotation-cache change: those inputs are
        # named by configuration, not by a file whose timestamp moves.
        if [[ ${rebuild} -eq 0 ]] && \
           ! validate_gofcards_exact_gof_cache "${target_tsv}" \
                 "vep_cache_version=${vep_version}" \
                 "vep_cache_dir=${vep_cache_dir}" \
                 "liftover_chain=${chain_file}" >/dev/null 2>&1; then
            log "GoFCards normalized cache failed validation or was built against different inputs"
            rebuild=1
        fi
        if [[ ${rebuild} -eq 0 ]] && \
           [[ -z $(find "${target_tsv}" -mtime "-${refresh_days}" -print -quit 2>/dev/null) ]]; then
            log "GoFCards normalized cache is past its ${refresh_days}-day refresh interval"
            rebuild=1
        fi
    fi

    if [[ ${rebuild} -eq 0 ]]; then
        log "GoFCards normalized cache is current: every input unchanged and within its ${refresh_days}-day refresh interval: ${target_tsv}"
        return 0
    fi

    [[ -f ${builder_script} ]] || { log "ERROR: GoFCards builder not found: ${builder_script}"; return 1; }
    [[ -s ${hg19_fasta} ]] || { log "ERROR: GRCh37 FASTA not found for GoFCards build: ${hg19_fasta}"; return 1; }
    [[ -s ${hg38_fasta} ]] || { log "ERROR: GRCh38 FASTA not found for GoFCards build; set gofcards_hg38_fasta"; return 1; }
    [[ -s ${chain_file} ]] || { log "ERROR: hg19->hg38 chain file not found: ${chain_file}"; return 1; }
    [[ -s ${hgnc_table} ]] || { log "ERROR: HGNC complete-set table not found: ${hgnc_table}"; return 1; }
    command -v vep >/dev/null 2>&1 || { log "ERROR: vep is not on PATH"; return 1; }
    command -v bcftools >/dev/null 2>&1 || { log "ERROR: bcftools is not on PATH"; return 1; }
    command -v CrossMap >/dev/null 2>&1 || { log "ERROR: CrossMap is not on PATH; install it into the PriVA environment"; return 1; }

    mkdir -p "$(dirname "${target_tsv}")" "${workdir}" || return 1

    log "Building GoFCards exact GOF cache -> ${target_tsv}"
    python "${builder_script}" \
        --workdir "${workdir}" \
        --out-json "${target_tsv}" \
        --hg19-fasta "${hg19_fasta}" \
        --hg38-fasta "${hg38_fasta}" \
        --chain "${chain_file}" \
        --vep-cache-dir "${vep_cache_dir}" \
        --vep-cache-version "${vep_version}" \
        --hgnc-table "${hgnc_table}" \
        --mechanism-review-tsv "${review_tsv}" \
        --stats-json "${stats_json}" 2>&1 | tee "${workdir}/gofcards_build.log"

    local build_status=${PIPESTATUS[0]}
    [[ ${build_status} -eq 0 ]] || {
        log "ERROR: GoFCards cache build failed; see ${workdir}/gofcards_build.log"
        return 1
    }

    validate_gofcards_exact_gof_cache "${target_tsv}" || {
        log "ERROR: GoFCards normalized cache validation failed"
        return 1
    }
    update_or_append_yaml "${config_file}" "gofcards_builder_script" "${builder_script}"
    update_or_append_yaml "${config_file}" "gofcards_mechanism_review_tsv" "${review_tsv}"
    log "GoFCards normalized cache built: ${target_tsv}"
}


function validate_gene_nonlof_mechanism_cache() {
    local builder_script=${1}
    local cache_json=${2}
    local schema_json=${3}

    [[ -f ${builder_script} ]] || { log "ERROR: Non-LOF builder not found: ${builder_script}"; return 1; }
    [[ -s ${cache_json} ]] || { log "ERROR: Missing/empty non-LOF cache: ${cache_json}"; return 1; }
    [[ -s ${schema_json} ]] || { log "ERROR: Missing/empty non-LOF schema: ${schema_json}"; return 1; }
    python "${builder_script}" \
        --output-schema "${schema_json}" \
        --validate-only "${cache_json}"
}


function gene_nonlof_mechanism_cache_install() {
    # The second argument selects the phase. "build" (the default) produces the
    # cache. "sources-only" stops after the downloads, which is how the
    # installer puts the weekly ClinVar VCV XML in place before ClinVar
    # injection needs it: injection writes the GoFCards cache this builder
    # reads, so the full build cannot also be what fetches that XML.
    local config_file=${1}
    local mode=${2:-build}
    [[ -f ${config_file} ]] || { log "ERROR: Config file not found: ${config_file}"; return 1; }
    [[ ${mode} == "build" || ${mode} == "sources-only" ]] || {
        log "ERROR: Unknown gene_nonlof_mechanism_cache_install mode: ${mode} (expected build or sources-only)"
        return 1
    }

    local cache_dir
    local raw_dir
    local builder_script
    local cache_json
    local schema_json
    local gofcards_variants
    local hgnc_table

    cache_dir=$(yaml_value_or_default "${config_file}" "gene_mechanism_cache_dir" "${DATA_DIR}/patho_mechanism")
    raw_dir=$(yaml_value_or_default "${config_file}" "gene_mechanism_raw_dir" "${cache_dir}/raw")
    builder_script=$(yaml_value_or_default "${config_file}" "gene_nonlof_mechanism_builder_script" "${SCRIPT_DIR}/build_gene_nonlof_mechanism_cache.py")
    cache_json=$(yaml_value_or_default "${config_file}" "gene_nonlof_mechanism_json" "${cache_dir}/gene_nonlof_mechanism_curated_assertions.json.gz")
    schema_json=$(yaml_value_or_default "${config_file}" "gene_nonlof_mechanism_schema" "${cache_dir}/gene_nonlof_mechanism_curated_assertions.schema.json")
    # The one deployed GoFCards cache, which already carries its ClinVar blocks.
    gofcards_variants=$(yaml_value_or_default "${config_file}" "gofcards_exact_gof_cache" "${DATA_DIR}/gofcards/gofcards_exact_gof.json.gz")
    hgnc_table=$(yaml_value_or_default "${config_file}" "hgnc_table" "${DATA_DIR}/hgnc/non_alt_loci_set.tsv")

    [[ -f ${builder_script} ]] || { log "ERROR: Non-LOF builder not found: ${builder_script}"; return 1; }
    [[ -s ${schema_json} ]] || { log "ERROR: Non-LOF schema not found: ${schema_json}"; return 1; }
    [[ -s ${hgnc_table} ]] || { log "ERROR: HGNC table not found: ${hgnc_table}"; return 1; }
    # Only the build phase reads the GoFCards cache. The sources-only phase runs
    # before injection has written that cache, which is why the phase exists.
    if [[ ${mode} == "build" ]]; then
        [[ -s ${gofcards_variants} ]] || { log "ERROR: GoFCards exact GOF cache not found (run gofcards_clinvar_injection_install first): ${gofcards_variants}"; return 1; }
    fi
    mkdir -p "${cache_dir}" "${raw_dir}" "$(dirname "${cache_json}")" || return 1

    local -a builder_args=(
        "${builder_script}"
        --cache-dir "${cache_dir}"
        --shared-raw-dir "${raw_dir}"
        --gofcards-exact-variants "${gofcards_variants}"
        --hgnc-table "${hgnc_table}"
        --output-json "${cache_json}"
        --output-schema "${schema_json}"
        --clinvar-min-review-stars "${CLINVAR_GOF_MIN_REVIEW_STARS:-2}"
        --clinvar-max-download-seconds "${CLINVAR_VCV_MAX_DOWNLOAD_SECONDS:-86400}"
        --timeout "${GENE_MECHANISM_TIMEOUT:-120}"
        --retries "${GENE_MECHANISM_RETRIES:-3}"
    )
    if [[ "${PRIVA_FORCE_NONLOF_MECHANISM_CACHE:-0}" == "1" ]] || \
       [[ "${PRIVA_FORCE_ALL_CACHES:-0}" == "1" ]]; then
        builder_args+=(--force)
    fi
    if [[ -n "${PROXY_URL:-}" ]]; then
        builder_args+=(--proxy-url "${PROXY_URL}" --download-tool auto)
    fi

    if [[ ${mode} == "sources-only" ]]; then
        builder_args+=(--fetch-sources-only)
        log "Checking PriVA non-LOF/ClinVar mechanism sources only; no cache build in this phase"
    else
        log "Checking PriVA non-LOF/ClinVar mechanism sources and rebuilding only when inputs changed"
    fi
    python "${builder_args[@]}" || {
        log "ERROR: PriVA non-LOF mechanism builder failed (mode=${mode})"
        return 1
    }

    if [[ ${mode} == "sources-only" ]]; then
        # The builder downloads into <raw_dir>/<source filename>, so that is
        # where the XML actually is. Record it under the key the injection step
        # reads, and only when the configured value names a different file, so a
        # deliberate equivalent path is left untouched.
        local fetched_clinvar_xml="${raw_dir}/clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz"
        [[ -s ${fetched_clinvar_xml} ]] || {
            log "ERROR: ClinVar VCV XML missing after the source fetch: ${fetched_clinvar_xml}"
            return 1
        }
        local configured_clinvar_xml
        configured_clinvar_xml=$(yaml_value_or_default "${config_file}" "clinvar_vcv_xml" "")
        if [[ $(readlink -f "${configured_clinvar_xml}") != $(readlink -f "${fetched_clinvar_xml}") ]]; then
            update_or_append_yaml "${config_file}" "clinvar_vcv_xml" "${fetched_clinvar_xml}"
        fi
        log "PriVA non-LOF mechanism sources staged; ClinVar VCV XML ready: ${fetched_clinvar_xml}"
        return 0
    fi

    validate_gene_nonlof_mechanism_cache \
        "${builder_script}" "${cache_json}" "${schema_json}" || {
        log "ERROR: PriVA non-LOF mechanism cache validation failed"
        return 1
    }

    update_or_append_yaml "${config_file}" "gene_nonlof_mechanism_builder_script" "${builder_script}"
    update_or_append_yaml "${config_file}" "gene_nonlof_mechanism_schema" "${schema_json}"
    update_or_append_yaml "${config_file}" "gene_nonlof_mechanism_json" "${cache_json}"
    log "PriVA non-LOF mechanism cache deployed: ${cache_json}"
}


function validate_hpo_condition_mechanism_cache() {
    local builder_script=${1}
    local cache_json=${2}
    shift 2

    [[ -f ${builder_script} ]] || { log "ERROR: HPO condition mechanism builder not found: ${builder_script}"; return 1; }
    [[ -s ${cache_json} ]] || { log "ERROR: Missing/empty HPO condition mechanism cache: ${cache_json}"; return 1; }
    python "${builder_script}" --validate-only "${cache_json}" "$@"
}


function hpo_condition_mechanism_cache_install() {
    local config_file=${1}
    [[ -f ${config_file} ]] || { log "ERROR: Config file not found: ${config_file}"; return 1; }

    local builder_script
    local cache_json
    local hpo_assertions
    local mechanism_evidence
    local clingen_dosage
    local gofcards_variants
    local hpo_release
    local mondo_release

    builder_script=$(yaml_value_or_default "${config_file}" "hpo_condition_mechanism_builder_script" "${SCRIPT_DIR}/build_hpo_condition_mechanism_cache.py")
    cache_json=$(yaml_value_or_default "${config_file}" "hpo_condition_mechanism_json" "${DATA_DIR}/patho_mechanism/hpo_condition_mechanism_cache.json.gz")
    hpo_assertions=$(yaml_value_or_default "${config_file}" "hpo_assertions" "${DATA_DIR}/hpo/genes_to_phenotype.assertions.tsv.gz")
    mechanism_evidence=$(yaml_value_or_default "${config_file}" "ddg2p_mechanism_evidence" "${DATA_DIR}/patho_mechanism/gene_pathogenic_mechanism_evidence.tsv")
    clingen_dosage=$(yaml_value_or_default "${config_file}" "gene_dosage_sensitivity" "${DATA_DIR}/clingen/gene_dosage_sensitivity.hg38.tsv")
    # The one deployed GoFCards cache. Every condition link this builder emits
    # comes from the ClinVar blocks nested inside it.
    gofcards_variants=$(yaml_value_or_default "${config_file}" "gofcards_exact_gof_cache" "${DATA_DIR}/gofcards/gofcards_exact_gof.json.gz")
    hpo_release=$(yaml_value_or_default "${config_file}" "hpo_release" "")
    mondo_release=$(yaml_value_or_default "${config_file}" "mondo_release" "")

    [[ -f ${builder_script} ]] || { log "ERROR: HPO condition mechanism builder not found: ${builder_script}"; return 1; }
    [[ -s ${hpo_assertions} ]] || { log "ERROR: Scoped HPO assertions not found: ${hpo_assertions}"; return 1; }
    [[ -s ${mechanism_evidence} ]] || { log "ERROR: Condition mechanism evidence not found: ${mechanism_evidence}"; return 1; }
    [[ -s ${clingen_dosage} ]] || { log "ERROR: ClinGen dosage evidence not found: ${clingen_dosage}"; return 1; }
    [[ -s ${gofcards_variants} ]] || { log "ERROR: GoFCards exact GOF cache not found (run gofcards_clinvar_injection_install first): ${gofcards_variants}"; return 1; }
    mkdir -p "$(dirname "${cache_json}")" || return 1

    local -a provenance_args=(
        --hpo-assertions "${hpo_assertions}"
        --mechanism-evidence "${mechanism_evidence}"
        --clingen-dosage "${clingen_dosage}"
        --gofcards-variants "${gofcards_variants}"
        --hpo-release "${hpo_release}"
        --mondo-release "${mondo_release}"
    )

    local rebuild=0
    if [[ "${PRIVA_FORCE_HPO_MECHANISM_JSON:-0}" == "1" ]] || \
       [[ "${PRIVA_FORCE_ALL_CACHES:-0}" == "1" ]] || [[ ! -s ${cache_json} ]]; then
        rebuild=1
    elif ! validate_hpo_condition_mechanism_cache \
        "${builder_script}" "${cache_json}" "${provenance_args[@]}" \
        >/dev/null 2>&1; then
        rebuild=1
    fi

    if [[ ${rebuild} -eq 1 ]]; then
        log "Building atomic HPO-framed condition mechanism JSON: ${cache_json}"
        python "${builder_script}" \
            --hpo-assertions "${hpo_assertions}" \
            --mechanism-evidence "${mechanism_evidence}" \
            --clingen-dosage "${clingen_dosage}" \
            --gofcards-variants "${gofcards_variants}" \
            --hpo-release "${hpo_release}" \
            --mondo-release "${mondo_release}" \
            --output "${cache_json}" || {
                log "ERROR: Failed to build HPO condition mechanism cache"
                return 1
            }
    else
        log "HPO condition mechanism cache is already current: ${cache_json}"
    fi

    validate_hpo_condition_mechanism_cache \
        "${builder_script}" "${cache_json}" "${provenance_args[@]}" || {
        log "ERROR: HPO condition mechanism cache validation failed"
        return 1
    }
    update_or_append_yaml "${config_file}" "hpo_condition_mechanism_builder_script" "${builder_script}"
    update_or_append_yaml "${config_file}" "hpo_condition_mechanism_json" "${cache_json}"
    log "HPO condition mechanism cache deployed: ${cache_json}"
}


function download_resource_atomic() {
    local url=${1}
    local target=${2}
    local temporary="${target}.download.$$"

    mkdir -p "$(dirname "${target}")" || return 1
    rm -f "${temporary}"
    if command -v wget >/dev/null 2>&1; then
        wget -c "${url}" -O "${temporary}" || { rm -f "${temporary}"; return 1; }
    elif command -v curl >/dev/null 2>&1; then
        curl -fL "${url}" -o "${temporary}" || { rm -f "${temporary}"; return 1; }
    else
        log "ERROR: wget or curl is required to download ${url}"
        return 1
    fi
    [[ -s ${temporary} ]] || { log "ERROR: Empty download from ${url}"; rm -f "${temporary}"; return 1; }
    mv -f "${temporary}" "${target}"
}


function validate_mondo_hpo_scope_cache() {
    local mondo_obo=${1}
    local registry=${2}
    local hpo_assertions=${3}
    [[ -s ${mondo_obo} ]] || { log "ERROR: Missing MONDO OBO: ${mondo_obo}"; return 1; }
    [[ -s ${registry} ]] || { log "ERROR: Missing disease-scope registry: ${registry}"; return 1; }
    [[ -s ${hpo_assertions} ]] || { log "ERROR: Missing scoped HPO assertions: ${hpo_assertions}"; return 1; }
    grep -q '^data-version: mondo/' "${mondo_obo}" || {
        log "ERROR: MONDO OBO has no valid data-version header: ${mondo_obo}"
        return 1
    }
    grep -q '^id: MONDO:0003847$' "${mondo_obo}" || {
        log "ERROR: MONDO hereditary disease root is missing: ${mondo_obo}"
        return 1
    }

    python - "${registry}" "${hpo_assertions}" <<'PY'
import csv
import gzip
import sys

registry_path, assertions_path = sys.argv[1:]

def inspect(path, required):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise SystemExit(f"{path} missing required columns: {', '.join(missing)}")
        count = sum(1 for _ in reader)
    if count == 0:
        raise SystemExit(f"{path} contains no data rows")
    return count

registry_count = inspect(registry_path, {
    "disease_id", "mondo_id", "disease_scope", "priva_scope",
    "scope_evidence", "scope_review_status",
})
assertion_count = inspect(assertions_path, {
    "gene_symbol", "disease_id", "hpo_id", "frequency", "evidence",
    "reference", "mondo_id", "disease_scope", "priva_scope",
    "scope_evidence", "scope_review_status",
})
print(f"disease_scope_rows={registry_count}; scoped_hpo_assertions={assertion_count}")
PY
}


function mondo_hpo_scope_install() {
    local config_file=${1}
    [[ -f ${config_file} ]] || { log "ERROR: Config file not found: ${config_file}"; return 1; }

    local hpo_release
    local hpo_raw_dir
    local hpo_genes_url
    local hpo_hpoa_url
    local hpo_genes
    local hpo_hpoa
    local hpo_base_assertions
    local hpo_assertions
    local mondo_release
    local mondo_url
    local mondo_sha256
    local mondo_obo
    local scope_registry
    local scope_overrides
    local scope_builder
    local hpo_builder

    hpo_release=$(yaml_value_or_default "${config_file}" "hpo_release" "v2026-06-23")
    hpo_raw_dir=$(yaml_value_or_default "${config_file}" "hpo_raw_dir" "${DATA_DIR}/hpo/raw/${hpo_release}")
    hpo_genes_url=$(yaml_value_or_default "${config_file}" "hpo_genes_to_phenotype_url" "https://github.com/obophenotype/human-phenotype-ontology/releases/download/${hpo_release}/genes_to_phenotype.txt")
    hpo_hpoa_url=$(yaml_value_or_default "${config_file}" "hpo_phenotype_hpoa_url" "https://github.com/obophenotype/human-phenotype-ontology/releases/download/${hpo_release}/phenotype.hpoa")
    hpo_genes="${hpo_raw_dir}/genes_to_phenotype.txt"
    hpo_hpoa="${hpo_raw_dir}/phenotype.hpoa"
    hpo_base_assertions="${hpo_raw_dir}/genes_to_phenotype.assertions.base.tsv.gz"
    hpo_assertions=$(yaml_value_or_default "${config_file}" "hpo_assertions" "${DATA_DIR}/hpo/genes_to_phenotype.assertions.tsv.gz")
    hpo_builder=$(yaml_value_or_default "${config_file}" "hpo_assertion_builder_script" "${SCRIPT_DIR}/collapse_HPO_anno.py")

    mondo_release=$(yaml_value_or_default "${config_file}" "mondo_release" "v2026-07-06")
    mondo_url=$(yaml_value_or_default "${config_file}" "mondo_url" "https://github.com/monarch-initiative/mondo/releases/download/${mondo_release}/mondo-simple.obo")
    mondo_sha256=$(yaml_value_or_default "${config_file}" "mondo_sha256" "")
    mondo_obo=$(yaml_value_or_default "${config_file}" "mondo_obo" "${DATA_DIR}/mondo/mondo-simple.${mondo_release}.obo")
    scope_registry=$(yaml_value_or_default "${config_file}" "mondo_disease_scope_registry" "${DATA_DIR}/mondo/disease_scope.tsv.gz")
    scope_overrides=$(yaml_value_or_default "${config_file}" "mondo_disease_scope_overrides" "${DATA_DIR}/mondo/disease_scope_overrides.tsv")
    scope_builder=$(yaml_value_or_default "${config_file}" "mondo_scope_builder_script" "${SCRIPT_DIR}/build_mondo_disease_scope.py")

    [[ -f ${hpo_builder} ]] || { log "ERROR: HPO assertion builder not found: ${hpo_builder}"; return 1; }
    [[ -f ${scope_builder} ]] || { log "ERROR: MONDO scope builder not found: ${scope_builder}"; return 1; }
    mkdir -p "${hpo_raw_dir}" "$(dirname "${hpo_assertions}")" "$(dirname "${mondo_obo}")" || return 1

    if [[ "${PRIVA_FORCE_HPO_ASSERTIONS:-0}" == "1" ]] || \
       [[ "${PRIVA_FORCE_ALL_CACHES:-0}" == "1" ]] || [[ ! -s ${hpo_genes} ]]; then
        log "Downloading HPO ${hpo_release} genes_to_phenotype.txt"
        download_resource_atomic "${hpo_genes_url}" "${hpo_genes}" || return 1
    fi
    if [[ "${PRIVA_FORCE_HPO_ASSERTIONS:-0}" == "1" ]] || \
       [[ "${PRIVA_FORCE_ALL_CACHES:-0}" == "1" ]] || [[ ! -s ${hpo_hpoa} ]]; then
        log "Downloading HPO ${hpo_release} phenotype.hpoa"
        download_resource_atomic "${hpo_hpoa_url}" "${hpo_hpoa}" || return 1
    fi

    local mondo_valid=0
    if [[ -s ${mondo_obo} ]] && grep -q '^data-version: mondo/' "${mondo_obo}"; then
        mondo_valid=1
        if [[ -n ${mondo_sha256} ]]; then
            local actual_sha256
            actual_sha256=$(sha256sum "${mondo_obo}" | awk '{print $1}')
            [[ ${actual_sha256} == "${mondo_sha256}" ]] || mondo_valid=0
        fi
    fi
    if [[ "${PRIVA_FORCE_MONDO_CACHE:-0}" == "1" ]] || \
       [[ "${PRIVA_FORCE_ALL_CACHES:-0}" == "1" ]] || [[ ${mondo_valid} -eq 0 ]]; then
        log "Downloading MONDO ${mondo_release} ontology"
        download_resource_atomic "${mondo_url}" "${mondo_obo}" || return 1
        if [[ -n ${mondo_sha256} ]]; then
            local downloaded_sha256
            downloaded_sha256=$(sha256sum "${mondo_obo}" | awk '{print $1}')
            [[ ${downloaded_sha256} == "${mondo_sha256}" ]] || {
                log "ERROR: MONDO SHA256 mismatch: expected ${mondo_sha256}, found ${downloaded_sha256}"
                return 1
            }
        fi
    fi

    local rebuild_hpo=0
    if [[ "${PRIVA_FORCE_HPO_ASSERTIONS:-0}" == "1" ]] || \
       [[ "${PRIVA_FORCE_ALL_CACHES:-0}" == "1" ]] || [[ ! -s ${hpo_base_assertions} ]] || \
       [[ ${hpo_genes} -nt ${hpo_base_assertions} ]] || [[ ${hpo_hpoa} -nt ${hpo_base_assertions} ]] || \
       [[ ${hpo_builder} -nt ${hpo_base_assertions} ]]; then
        rebuild_hpo=1
        log "Building release-matched base HPO assertions"
        python "${hpo_builder}" "${hpo_genes}" "${hpo_hpoa}" "${hpo_base_assertions}" || return 1
    fi

    local rebuild_scope=0
    if [[ ${rebuild_hpo} -eq 1 ]] || [[ "${PRIVA_FORCE_MONDO_CACHE:-0}" == "1" ]] || \
       [[ "${PRIVA_FORCE_ALL_CACHES:-0}" == "1" ]] || \
       [[ ! -s ${scope_registry} ]] || [[ ! -s ${hpo_assertions} ]] || \
       [[ ${mondo_obo} -nt ${scope_registry} ]] || [[ ${hpo_base_assertions} -nt ${scope_registry} ]] || \
       [[ ${scope_builder} -nt ${scope_registry} ]] || [[ ${scope_overrides} -nt ${scope_registry} ]]; then
        rebuild_scope=1
    fi
    if [[ ${rebuild_scope} -eq 1 ]]; then
        log "Building MONDO/HPO disease-scope registry and scoped HPO assertions"
        local -a scope_args=(
            "${scope_builder}"
            --mondo-obo "${mondo_obo}"
            --hpo-assertions "${hpo_base_assertions}"
            --registry-output "${scope_registry}"
            --annotated-assertions-output "${hpo_assertions}"
        )
        [[ -s ${scope_overrides} ]] && scope_args+=(--manual-overrides "${scope_overrides}")
        python "${scope_args[@]}" || return 1
    fi

    validate_mondo_hpo_scope_cache "${mondo_obo}" "${scope_registry}" "${hpo_assertions}" || return 1
    update_yaml "${config_file}" "hpo_assertions" "${hpo_assertions}"
    update_yaml "${config_file}" "mondo_obo" "${mondo_obo}"
    update_yaml "${config_file}" "mondo_disease_scope_registry" "${scope_registry}"
    log "MONDO/HPO disease-scope resources deployed: ${scope_registry}; ${hpo_assertions}"
}


function mechanism_resource_install() {
    local config_file=${1}

    # Mechanism sources use OMIM and ORPHA identifiers. Build the pinned MONDO
    # registry first so the mechanism cache can resolve both identifier systems
    # to one condition without disease-name matching. This ordering is part of
    # the data contract: the mechanism builder must never fall back to a stale
    # or absent disease-scope registry and silently emit unscoped assertions.
    mondo_hpo_scope_install "${config_file}" || return 1
    gene_pathogenic_mechanism_cache_install "${config_file}" || return 1

    # The GoFCards cache is built in two steps, and the ClinVar VCV XML the
    # second step needs is downloaded by the non-LOF builder. That would be a
    # standstill -- injection waiting on the XML, the non-LOF build waiting on
    # injection's output -- so the downloads are pulled forward into their own
    # phase here. Order: normalize, fetch sources, inject, then build.
    gofcards_exact_gof_cache_install "${config_file}" || return 1
    gene_nonlof_mechanism_cache_install "${config_file}" sources-only || return 1
    gofcards_clinvar_injection_install "${config_file}" || return 1
    gene_nonlof_mechanism_cache_install "${config_file}" || return 1
    hpo_condition_mechanism_cache_install "${config_file}" || return 1
}


function main_install() {
    if [[ $# -eq 1 && ( ${1:-} == -h || ${1:-} == --help ) ]]; then
        printf 'Usage: bash scripts/install_utils.sh main_install CONFIG\n'
        printf 'Creates or reuses priva_acmg from %s/acmg_conda_rough.yml.\n' "${BASE_DIR}"
        printf 'The environment is activated before yq reads CONFIG.\n'
        return 0
    fi
    [[ $# -eq 1 && ${1:-} != -* && -r ${1:-} ]] || {
        log "ERROR: Supply exactly one readable CONFIG file; see main_install --help."
        return 2
    }
    local config_file
    config_file=$(realpath "$1") || return 1

    # Prepare the bundled environment before parsing the project configuration.
    conda_install_vep || return 1
    local conda_env_name=${PRIVA_INSTALL_ENV_NAME}

    # Resolve/create the selected folders before checking VEP destinations.
    # Existing resource functions and annotation commands read the saved keys.
    priva_configure_cache_paths "${config_file}" || return 1

    local has_error=0
    local skip_vep_install=0
    [[ "${PRIVA_SMOKE_SKIP_VEP:-0}" == "1" ]] && skip_vep_install=1
    # Read configuration
    local vep_cache_dir=$(read_yaml "$config_file" "vep_cache_dir")
    local vep_plugins_dir=$(read_yaml "$config_file" "vep_plugins_dir")
    local assembly=$(read_yaml "$config_file" "assembly")
    local ref_fasta=$(read_yaml "$config_file" "ref_genome")
    local vep_plugins_cachedir=$(read_yaml "$config_file" "vep_plugins_cachedir")
    if [[ ${skip_vep_install} -eq 0 ]]; then
        [[ -d ${vep_cache_dir} ]] || { log "The cache directory for VEP ${vep_cache_dir} is not found, please check the file"; has_error=1; }
        [[ -d ${vep_plugins_dir} ]] || { log "The plugins directory for VEP ${vep_plugins_dir} is not found, please check the file"; has_error=1; }
        [[ -d ${vep_plugins_cachedir} ]] || { log "The plugins cache directory for VEP ${vep_plugins_cachedir} is not found, please check the file"; has_error=1; }
    fi
    [[ ${has_error} -eq 1 ]] && { log "Please check the values in the configuration file ${config_file}"; return 1; }

    reference_genome_install "${config_file}" || \
    { log "Failed to deploy/index configured reference genome"; return 1; }
    ref_fasta=$(read_yaml "$config_file" "ref_genome")
    [[ -f ${ref_fasta} ]] || { log "The reference genome fasta file ${ref_fasta} is not found, please check the file"; return 1; }

    if [[ ${skip_vep_install} -eq 1 ]]; then
        log "PRIVA_SMOKE_SKIP_VEP=1: skipping VEP API/cache, VEP plugins, and LoFtee installation; smoke inputs must already contain complete CSQ annotations."
    else
        # 2. Install VEP
        # Update config with installation results
        local vep_help vep_version
        vep_help=$(vep --help 2>&1) || { log "ERROR: Installed Conda VEP/API is not runnable: ${vep_help}"; return 1; }
        vep_version=$(grep "ensembl-vep" <<< "$vep_help" | awk '{print $NF}' | awk -F '.' '{print $1}')
        [[ $vep_version =~ ^[0-9]+$ ]] || { log "ERROR: Could not identify the installed VEP release"; return 1; }
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
    fi

    # 4. Install gnomAD VCF (basically download bgzipped VCF files)
    gnomAD_install ${config_file} ${ref_fasta} || \
    { log "Failed to install gnomAD VCF"; return 1; }
    AlphaMissense_anno_gnomAD ${config_file} || \
    { log "Failed to annotate gnomAD data to AlphaMissense VCF"; return 1; }

    # 5. Install InterPro mapping pickle
    InterPro_parsing ${config_file} || \
    { log "Failed to install InterPro mapping pickle"; return 1; }

    # 6. Install ClinVar VCF
    ClinVar_VCF_deploy ${config_file} || \
    { log "Failed to install ClinVar VCF"; return 1; }
    ClinVar_AA_stat ${config_file} || \
    { log "Failed to stat ClinVar per AA change"; return 1; }
    ClinVar_Gene_stat ${config_file} || \
    { log "Failed to stat ClinVar records per gene"; return 1; }
    AlphaMissense_stat ${config_file} || \
    { log "Failed to generate the AlphaMissense statistics JSON file"; return 1; }
    AlphaMissense_pick_intolerant_domains ${config_file} || \
    { log "Failed to pick intolerant domains from AlphaMissense"; return 1; }
    prepare_pm1_regions ${config_file} || \
    { log "Failed to prepare PM1 regions pickle"; return 1; }
    ClinVar_patho_AF_stat ${config_file} || \
    { log "Failed to stat ClinVar pathogenic variants per allele frequency"; return 1; }


    # 7. Install CADD prescores
    CADD_install \
    ${config_file} \
    ${assembly} || \
    { log "Failed to install CADD prescores"; return 1; }

    # 8. ClinGen install
    ClinGen_deploy \
    ${config_file} || \
    { log "Failed to deploy ClinGen data"; return 1; }

    # 9. Install Ensembl CDS FASTA for alternative start codon detection
    CDS_FASTA_install \
    ${config_file} || \
    { log "Failed to install CDS FASTA for alternative start codon detection"; return 1; }

    # 10. Install pext (proportion expressed across transcripts) BigWig files
    pext_install \
    ${config_file} || \
    { log "Failed to install pext BigWig files"; return 1; }

    # 11. Install pinned HPO/MONDO, condition mechanisms, and exact GoFCards
    mechanism_resource_install \
        ${config_file} || \
    { log "Failed to install HPO/MONDO or mechanism resources"; return 1; }
    
    log "Congratulations! The installation is completed!"
}



if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    "$@"
fi
