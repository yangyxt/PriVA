#!/usr/bin/env bash
#SBATCH --job-name=PriVA_Installer_Rerun
#SBATCH --partition=copy,long
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --no-requeue
# Reuse a copied diagnostic environment only after exact-lock and smoke checks.
set -euo pipefail
[[ $# == 5 ]] || { echo 'Usage: sbatch rerun_cluster_install.sbatch LOG REPO REVISION PREVIOUS_DIAGNOSTIC NEW_DIAGNOSTIC' >&2; exit 64; }
log_pattern=$1; repo=$2; revision=$3; previous=$4; diagnostic=$5
log_file=${log_pattern//%j/${SLURM_JOB_ID:?}}
main() {
    printf '[%s] Job %s on %s\n' "$(date -Is)" "$SLURM_JOB_ID" "$(hostname)"
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" PYTHONUNBUFFERED=1
    [[ ! -e $diagnostic ]] || { echo 'Refusing existing diagnostic destination'; return 65; }
    mkdir -p "$diagnostic/source"
    git -C "$repo" archive "$revision" | tar -xf - -C "$diagnostic/source"
    printf '%s\n' "$revision" > "$diagnostic/source/REFERENCE_COMMIT"
    export MAMBA_EXE="$previous/bin/micromamba"
    "$MAMBA_EXE" env export --name priva_acmg --explicit > "$diagnostic/environment.conda-explicit.txt"
    # An explicit export describes package identities independently of prefix.
    cmp "$previous/environment.conda-explicit.txt" "$diagnostic/environment.conda-explicit.txt"
    activation=$("$MAMBA_EXE" shell activate --shell bash --name priva_acmg)
    eval "$activation"
    [[ $CONDA_PREFIX == "$previous/manager/envs/priva_acmg" ]] || return 66
    export MAMBA_EXE="$previous/bin/micromamba"
    python "$diagnostic/source/scripts/check_cluster_environment.py" --spec "$diagnostic/source/envs/priva-cluster.yml" --output "$diagnostic/environment.ready.json"
    python "$diagnostic/source/scripts/prepare_cluster_config.py" \
        --template "$diagnostic/source/config.yaml" --source-root "$diagnostic/source" \
        --cache-root "$diagnostic/caches" \
        --data-root /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade \
        --threads 2 --output "$diagnostic/config.yaml"
    python - "$diagnostic/config.yaml" <<'PY'
import pathlib,sys,yaml
p=pathlib.Path(sys.argv[1]); data=yaml.safe_load(p.read_text())
data['ref_genome']='/stornext/snfs170/next-gen/scratch/yangyxt/caches/references/ucsc.hg19.fasta'
p.write_text(yaml.safe_dump(data,sort_keys=False))
PY
    sha256sum "$diagnostic/source/scripts/install_utils.sh" "$diagnostic/environment.conda-explicit.txt"
    printf 'Actual corrected command: timeout --signal=TERM --kill-after=30s 120s bash %q main_install %q\n' "$diagnostic/source/scripts/install_utils.sh" "$diagnostic/config.yaml"
    unset PRIVA_SMOKE_SKIP_VEP
    timeout --signal=TERM --kill-after=30s 120s bash "$diagnostic/source/scripts/install_utils.sh" main_install "$diagnostic/config.yaml"
}
( trap 'result=$?; printf "[%s] Exit status: %s\n" "$(date -Is)" "$result"; exit "$result"' EXIT; main ) > "$log_file" 2>&1
