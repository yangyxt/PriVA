#!/usr/bin/env bash
#SBATCH --job-name=PriVA_AlphaMissense_Cache
#SBATCH --partition=copy,long
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --no-requeue
# Submit with --dependency=afterok:VEP_CACHE_JOB, then LOG CONFIG SOURCE_PREFIX.
set -euo pipefail
[[ $# == 3 ]] || { echo 'Usage: sbatch --dependency=afterok:JOB build_cluster_alpha.sbatch LOG CONFIG SOURCE_PREFIX' >&2; exit 64; }
log_pattern=$1
config=$2
reference=$3
log_file=${log_pattern//%j/${SLURM_JOB_ID:?}}
main() {
    trap 'result=$?; printf "[%s] Exit status: %s\n" "$(date --iso-8601=seconds)" "$result"' EXIT
    printf '[%s] Job %s on %s\n' "$(date --iso-8601=seconds)" "$SLURM_JOB_ID" "$(hostname)"
    [[ $(cat "$reference/REFERENCE_COMMIT") == 09223ef5a8b0ceb4c6d61055c6eefe62b51ba903 ]] || return 65
    manager=/stornext/snfs170/next-gen/scratch/yangyxt/mamba/bin/micromamba
    activation=$("$manager" shell activate --shell bash --name priva_acmg)
    eval "$activation"
    export PYTHONUNBUFFERED=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
    cache=$(yq -r .vep_cache_dir "$config")
    [[ -s "$cache/homo_sapiens_merged/113_GRCh37/info.txt" ]] || {
        echo 'VEP113 GRCh37 cache is not installed; refusing partial annotation' >&2; return 66;
    }
    prefix=$(yq -r .alphamissense_intolerant_motifs "$config")
    prefix=${prefix%.pkl}
    bash "$reference/scripts/install_utils.sh" AlphaMissense_anno "$config"
    bash "$reference/scripts/install_utils.sh" AlphaMissense_pick_missense_constrained_segments "$config" "$prefix"
}
( main ) > "$log_file" 2>&1
