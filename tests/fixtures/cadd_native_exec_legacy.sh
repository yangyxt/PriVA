#!/usr/bin/env bash
# Named interpreter routes for unchanged CADD Python; scientific code stays upstream.
# Direct use: bash cadd_native_exec.sh run --runtime-root ROOT --name cadd-python2 -- PYTHON_ARGS...
set -euo pipefail

run() {
    local runtime_root= name= option
    if [[ ${0##*/} == cadd-native-exec || ${0##*/} == cadd_native_exec.sh ]]; then
        local parsed
        parsed=$(getopt -o '' --long runtime-root:,name: -- "$@") || return 64
        eval set -- "$parsed"
        while [[ $# -gt 0 ]]; do
            option=$1; shift
            case $option in
                --runtime-root) runtime_root=$1; shift ;;
                --name) name=$1; shift ;;
                --) break ;;
            esac
        done
    else
        runtime_root=$(cd -- "$(dirname -- "$(realpath -- "${BASH_SOURCE[0]}")")/.." && pwd)
        name=${0##*/}
    fi
    case $name in cadd-python2|cadd-esm-python|cadd-regseq-python|cadd-mmsplice-python) ;;
        *) printf 'Select a named CADD Python environment\n' >&2; return 64 ;;
    esac
    [[ -n $runtime_root ]] || { printf 'Supply --runtime-root\n' >&2; return 64; }
    local prefix="$runtime_root/envs/$name"
    [[ -x $prefix/bin/python ]] || { printf 'CADD environment is incomplete: %s\n' "$prefix" >&2; return 69; }
    unset PYTHONHOME PYTHONPATH
    export CONDA_PREFIX="$prefix" PATH="$prefix/bin:$PATH"
    export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
    export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
    export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    exec "$prefix/bin/python" "$@"
}

if [[ ${0##*/} == cadd-native-exec || ${0##*/} == cadd_native_exec.sh ]]; then
    [[ ${1:-} == run ]] || { printf 'Usage: %s run --runtime-root ROOT --name NAME -- PYTHON_ARGS...\n' "$0" >&2; exit 64; }
    shift
fi
run "$@"
