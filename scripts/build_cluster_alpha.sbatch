#!/usr/bin/env bash
#SBATCH --job-name=PriVA_AlphaMissense_Cache
#SBATCH --partition=copy,long
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --no-requeue
# Supply an explicitly approved mapping; this wrapper never selects/downloads one.

alpha_validate_mapping() {
    local mapping=$1 expected_sha=$2 actual_sha before after
    [[ $mapping == /* && -f $mapping && -r $mapping && -s $mapping ]] || {
        echo 'An explicit absolute, readable, nonempty NM-to-ENST mapping file is required' >&2; return 67;
    }
    [[ $expected_sha =~ ^[0-9a-f]{64}$ ]] || {
        echo 'An explicit approved mapping SHA256 is required' >&2; return 67;
    }
    before=$(stat -Lc '%d:%i:%s:%y:%z' -- "$mapping") || return
    actual_sha=$(sha256sum -- "$mapping") || return
    after=$(stat -Lc '%d:%i:%s:%y:%z' -- "$mapping") || return
    [[ $before == "$after" ]] || {
        echo 'NM-to-ENST mapping changed during checksum verification' >&2; return 67;
    }
    [[ ${actual_sha%% *} == "$expected_sha" ]] || {
        echo 'Approved NM-to-ENST mapping checksum mismatch' >&2; return 67;
    }
}

alpha_validate_source() {
    local reference=$1 actual_sha name expected_sha
    [[ $(cat "$reference/REFERENCE_COMMIT") == 09223ef5a8b0ceb4c6d61055c6eefe62b51ba903 ]] || {
        echo 'Unexpected frozen scientific source revision' >&2; return 65;
    }
    for name in install_utils.sh am_pick_missense_constrained_segments.py; do
        case $name in
            install_utils.sh) expected_sha=7e4c5e6fa3617d0b74fd5428e7cf0b8e731bd9018de89bdeb6242b7f6972cc7f ;;
            am_pick_missense_constrained_segments.py) expected_sha=279b376f2a2b3e669852d8bc72b9934efc919b9b940e1371e76cc136060e5e33 ;;
        esac
        actual_sha=$(sha256sum -- "$reference/scripts/$name") || return
        [[ ${actual_sha%% *} == "$expected_sha" ]] || {
            echo "Frozen scientific source checksum mismatch: $name" >&2; return 65;
        }
    done
}

alpha_run_stages() {
    local config=$1 reference=$2 mapping=$3 mapping_sha=$4
    local before_mapping after_mapping prefix assembly vcf output threads motif_py
    # Preserve the original wrapper's explicit prefix, computed before annotation.
    alpha_validate_mapping "$mapping" "$mapping_sha" || return
    before_mapping=$(stat -Lc '%d:%i:%s:%y:%z' -- "$mapping") || return
    prefix=$(yq -r .alphamissense_intolerant_motifs "$config") || return
    prefix=${prefix%.pkl}
    printf 'Approved NM-to-ENST mapping: %s; SHA256 %s\n' "$mapping" "$mapping_sha"
    bash "$reference/scripts/install_utils.sh" AlphaMissense_anno "$config" || return

    # Same reads, reuse predicates and Python option order as the original
    # 09223ef helper. Its command block SHA256 is
    # 569619479a4267693bc2f9caca56cc12761d8abb87ee3a232138b0624ae57b99.
    # Only the existing mapping option is added; scientific defaults are retained.
    assembly=$(yq -r .assembly "$config") || return
    vcf=$(yq -r .alphamissense_vep_vcf "$config") || return
    output=$(yq -r .alphamissense_intolerant_motifs "$config") || return
    threads=$(yq -r .threads "$config") || return
    motif_py="$reference/scripts/am_pick_missense_constrained_segments.py"
    alpha_validate_mapping "$mapping" "$mapping_sha" || return
    after_mapping=$(stat -Lc '%d:%i:%s:%y:%z' -- "$mapping") || return
    [[ $after_mapping == "$before_mapping" ]] || {
        echo 'NM-to-ENST mapping changed during annotation' >&2; return 67;
    }
    if [[ -f $output && $output -nt $vcf && $output -nt $motif_py ]]; then
        echo "The AlphaMissense intolerant motifs file $output is already generated, skip this step" >&2
        echo 'Reused existing output under original timestamp predicates; its mapping provenance is unverified.'
        return 0
    fi
    echo "The AlphaMissense intolerant motifs file $output is not generated, we will generate it" >&2
    # The original helper receives this explicit prefix; its final destination
    # is prefix.pkl. Preserve that behavior and do not create/redirect outputs.
    [[ -d $(dirname "$output") ]] || echo "The directory $(dirname "$output") does not exist, we will use the prefix $prefix" >&2
    output="${prefix}.pkl"
    python "$motif_py" \
        --vcf_path "$vcf" \
        --rmc_tsv "$reference/scripts/../data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.${assembly}.tsv" \
        --hcseeker_tsv "$reference/scripts/../data/HCSeeker/HC_spots.${assembly}.tsv" \
        --output "$output" \
        --processes "$threads" \
        --nm_to_enst_tsv "$mapping" || return
    alpha_validate_mapping "$mapping" "$mapping_sha" || return
    after_mapping=$(stat -Lc '%d:%i:%s:%y:%z' -- "$mapping") || return
    [[ $after_mapping == "$before_mapping" ]] || {
        echo 'NM-to-ENST mapping changed during motif generation; output is unqualified' >&2; return 67;
    }
    bash "$reference/scripts/common_bash_utils.sh" update_yaml "$config" alphamissense_intolerant_motifs "$output" || return
    echo "The AlphaMissense intolerant motifs file is generated and saved to $output" >&2
    echo 'Fresh motif generation completed with the explicit mapping; scientific qualification remains separate.'
}

main() {
    local config=$1 reference=$2 mapping=$3 mapping_sha=$4 manager activation cache
    trap 'result=$?; printf "[%s] Exit status: %s\n" "$(date --iso-8601=seconds)" "$result"' EXIT
    printf '[%s] Job %s on %s\n' "$(date --iso-8601=seconds)" "$SLURM_JOB_ID" "$(hostname)"
    # Missing/wrong mapping fails before environment activation or annotation.
    alpha_validate_mapping "$mapping" "$mapping_sha" || return
    alpha_validate_source "$reference" || return
    manager=/stornext/snfs170/next-gen/scratch/yangyxt/mamba/bin/micromamba
    activation=$("$manager" shell activate --shell bash --name priva_acmg) || return
    eval "$activation"
    export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
    cache=$(yq -r .vep_cache_dir "$config") || return
    [[ -s "$cache/homo_sapiens_merged/113_GRCh37/info.txt" ]] || {
        echo 'VEP113 GRCh37 cache is not installed; refusing partial annotation' >&2; return 66;
    }
    alpha_run_stages "$config" "$reference" "$mapping" "$mapping_sha"
}

if [[ ${BASH_SOURCE[0]} == "$0" ]]; then
    set -euo pipefail
    [[ $# == 5 ]] || { echo 'Usage: sbatch --dependency=afterok:JOB build_cluster_alpha.sbatch LOG CONFIG SOURCE_PREFIX APPROVED_NM_TO_ENST_TSV SHA256' >&2; exit 64; }
    log_pattern=$1
    log_file=${log_pattern//%j/${SLURM_JOB_ID:?}}
    ( main "$2" "$3" "$4" "$5" ) > "$log_file" 2>&1
fi
