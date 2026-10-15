#!/usr/bin/env bash
# Per-environment entry points for the unmodified upstream CADD 1.7 container.
# Install as bin/cadd-container-exec and symlink the names handled below.
set -euo pipefail
runtime_root=$(cd -- "$(dirname -- "$(realpath -- "${BASH_SOURCE[0]}")")/.." && pwd)
image="$runtime_root/cadd-upstream-0.1.1.sif"
case ${0##*/} in
    cadd-python2) prefix=/conda-envs/fe7049a03d989a2faa6d7688f5987d09; executable=python ;;
    cadd-vep) prefix=/conda-envs/0b37582dec236e35970e555dc2c0f0b5; executable=vep ;;
    cadd-esm-python) prefix=/conda-envs/10a710eb7e199e21064fc58a16f9ea4f; executable=python ;;
    cadd-regseq-python) prefix=/conda-envs/f675d936396bfaffdd96e0282b787268; executable=python ;;
    cadd-mmsplice-python) prefix=/conda-envs/27e6f085c62573d05046d5330f07505d; executable=python ;;
    *) echo 'Use a cadd-python2/vep/esm-python/regseq-python/mmsplice-python entry point' >&2; exit 64 ;;
esac
[[ -s $image ]] || { echo "CADD runtime image has not finished installation: $image" >&2; exit 69; }
export APPTAINERENV_PREPEND_PATH="$prefix/bin"
export APPTAINERENV_CONDA_PREFIX="$prefix"
export APPTAINERENV_OMP_NUM_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
export APPTAINERENV_OPENBLAS_NUM_THREADS=1
export APPTAINERENV_PYTHONDONTWRITEBYTECODE=1
export APPTAINERENV_PYTHONNOUSERSITE=1
if [[ -n ${CADD:-} ]]; then export APPTAINERENV_CADD=$CADD; fi
if [[ -n ${TMPDIR:-} ]]; then export APPTAINERENV_TMPDIR=$TMPDIR; fi
# Source/cache files are shared read-only. Analysis staging/results are under
# the user's separately authorized data root; no copy of the caches is made.
exec apptainer exec --cleanenv --no-home \
    --bind /stornext/snfs170/next-gen/scratch/yangyxt:/stornext/snfs170/next-gen/scratch/yangyxt:ro \
    --bind /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade:/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade:rw \
    "$image" "$prefix/bin/$executable" "$@"
