#!/usr/bin/env bash
#SBATCH --job-name=PriVA_Installer_Diagnostic
#SBATCH --partition=copy,long
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --no-requeue
# A real, bounded main_install run. Its API/cache writes are fully isolated.
set -euo pipefail
[[ $# == 5 ]] || { echo 'Usage: sbatch diagnose_cluster_install.sbatch LOG REPO REVISION DIAGNOSTIC_ROOT STATE_DIR' >&2; exit 64; }
log_pattern=$1
repo=$2
revision=$3
diagnostic=$4
state=$5
log_file=${log_pattern//%j/${SLURM_JOB_ID:?}}
main() {
    printf '[%s] Job %s on %s\n' "$(date -Is)" "$SLURM_JOB_ID" "$(hostname)"
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" PYTHONUNBUFFERED=1
    [[ ! -e $diagnostic ]] || { echo 'Refusing an existing diagnostic directory'; return 65; }
    mkdir -p "$diagnostic/source" "$diagnostic/bin" "$diagnostic/manager" "$diagnostic/tmp"
    git -C "$repo" archive "$revision" | tar -xf - -C "$diagnostic/source"
    printf '%s\n' "$revision" > "$diagnostic/source/REFERENCE_COMMIT"
    sha256sum "$diagnostic/source/scripts/install_utils.sh" "$state/environment.ready.json"
    /stornext/snfs170/next-gen/scratch/yangyxt/mamba/bin/micromamba env export --name priva_acmg --explicit > "$diagnostic/environment.conda-explicit.txt"
    grep -qx '@EXPLICIT' "$diagnostic/environment.conda-explicit.txt"
    sha256sum "$diagnostic/environment.conda-explicit.txt"
    cat > "$diagnostic/bin/micromamba" <<'MANAGER'
#!/usr/bin/env bash
set -euo pipefail
diagnostic=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
exec env MAMBA_ROOT_PREFIX="$diagnostic/manager" CONDA_ENVS_PATH="$diagnostic/manager/envs" CONDA_PKGS_DIRS=/stornext/snfs170/next-gen/scratch/yangyxt/mamba/pkgs XDG_CONFIG_HOME="$diagnostic/manager" XDG_CACHE_HOME="$diagnostic/manager" TMPDIR="$diagnostic/tmp" /stornext/snfs170/next-gen/scratch/yangyxt/mamba/libexec/micromamba --no-rc "$@"
MANAGER
    chmod +x "$diagnostic/bin/micromamba"
    export MAMBA_EXE="$diagnostic/bin/micromamba"
    "$MAMBA_EXE" create --name priva_acmg --file "$diagnostic/environment.conda-explicit.txt" --always-copy --offline --yes
    activation=$("$MAMBA_EXE" shell activate --shell bash --name priva_acmg)
    eval "$activation"
    [[ $CONDA_PREFIX == "$diagnostic/manager/envs/priva_acmg" ]] || return 66
    export MAMBA_EXE="$diagnostic/bin/micromamba"
    python - "$diagnostic" "$state" <<'PY'
import pathlib
import sys
import yaml
root = pathlib.Path(sys.argv[1])
source = root / "source"
project = pathlib.Path(sys.argv[2]).parent
template = yaml.safe_load((source / "config.yaml").read_text())
old = "/paedyl01/disk1/yangyxt/PriVA/"
for key, value in list(template.items()):
    if isinstance(value, str) and value.startswith(old):
        template[key] = str(source / value[len(old):])
template.update(cache_root=str(root / "caches"),
                ref_genome="/stornext/snfs170/next-gen/scratch/yangyxt/caches/references/ucsc.hg19.fasta",
                tmp_dir=str(project / "tmp/installer-diagnostic"), input_vcf="", ped_file="", threads=2,
                output_dir=str(project / "results/installer-diagnostic"), hub_file_dir=str(project / "annotation_hub/installer-diagnostic"))
(root / "config.yaml").write_text(yaml.safe_dump(template, sort_keys=False))
PY
    printf 'Reference source: %s\nDiagnostic environment: %s\n' "$revision" "$CONDA_PREFIX"
    printf 'Actual command: timeout --signal=TERM --kill-after=30s 300s bash %q main_install %q\n' "$diagnostic/source/scripts/install_utils.sh" "$diagnostic/config.yaml"
    # No skip flag, function replacement or fake cache marker: the installer
    # runs normally. A timeout is an operator checkpoint, never a passing test.
    unset PRIVA_SMOKE_SKIP_VEP
    timeout --signal=TERM --kill-after=30s 300s bash "$diagnostic/source/scripts/install_utils.sh" main_install "$diagnostic/config.yaml"
}
( trap 'result=$?; printf "[%s] Exit status: %s\n" "$(date -Is)" "$result"; exit "$result"' EXIT; main ) > "$log_file" 2>&1
