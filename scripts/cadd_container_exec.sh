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
export APPTAINERENV_MKL_NUM_THREADS=1
export APPTAINERENV_NUMEXPR_NUM_THREADS=1
export APPTAINERENV_PYTHONDONTWRITEBYTECODE=1
export APPTAINERENV_PYTHONNOUSERSITE=1
if [[ -n ${CADD:-} ]]; then export APPTAINERENV_CADD=$CADD; fi
if [[ -n ${TMPDIR:-} ]]; then export APPTAINERENV_TMPDIR=$TMPDIR; fi
# Bind selected resource paths, including symlink destinations, instead of a
# site-specific home/data root. Analysis staging remains writable.
readonly_paths=("$runtime_root")
if [[ -n ${CADD:-} && -d $CADD ]]; then
    readonly_paths+=("$CADD" "$(dirname -- "$(dirname -- "$CADD")")")
fi
for name in ${!CADD_@}; do
    value=${!name}
    if [[ $value == /* && -e $value ]]; then
        resolved=$(realpath -e -- "$value")
        [[ -f $resolved ]] && resolved=$(dirname -- "$resolved")
        readonly_paths+=("$resolved" "$(dirname -- "$resolved")")
    fi
done
for value in "$@"; do
    if [[ $value == /* && -e $value ]]; then
        resolved=$(realpath -e -- "$value")
        [[ -f $resolved ]] && resolved=$(dirname -- "$resolved")
        readonly_paths+=("$resolved")
    fi
done
bind_args=()
declare -A bound=()
for value in "${readonly_paths[@]}"; do
    [[ -n ${bound[$value]:-} ]] && continue
    [[ $value != *:* && $value != *,* ]] || { echo 'Apptainer resource paths cannot contain colon or comma' >&2; exit 64; }
    bind_args+=(--bind "$value:$value:ro")
    bound[$value]=1
done
for value in "$PWD" "${TMPDIR:-$PWD}"; do
    [[ -d $value && $value != *:* && $value != *,* ]] || { echo 'Invalid CADD staging directory' >&2; exit 64; }
    bind_args+=(--bind "$value:$value:rw")
done
exec apptainer exec --cleanenv --no-home \
    "${bind_args[@]}" \
    "$image" "$prefix/bin/$executable" "$@"
