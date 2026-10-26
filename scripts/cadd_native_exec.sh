#!/usr/bin/env bash
# Named interpreter routes for unchanged CADD Python; scientific code stays upstream.
# Direct use: bash cadd_native_exec.sh run --runtime-root ROOT --name cadd-python2 -- PYTHON_ARGS...
set -euo pipefail

run() {
    local runtime_root= name= option
    if [[ ${0##*/} == cadd-native-exec || ${0##*/} == cadd_native_exec.sh ]]; then
        while [[ $# -gt 0 ]]; do
            option=$1; shift
            case $option in
                --runtime-root|--name)
                    [[ $# -gt 0 && -n $1 ]] || { printf 'Missing value for %s\n' "$option" >&2; return 64; }
                    if [[ $option == --runtime-root ]]; then runtime_root=$1; else name=$1; fi
                    shift ;;
                --) break ;;
                *) printf 'Unknown CADD route option: %s\n' "$option" >&2; return 64 ;;
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
    if [[ -e $runtime_root/native-prefixes || -L $runtime_root/native-prefixes ]]; then
        local metadata="$runtime_root/native-prefixes/$name.prefix" extra=
        [[ -d $runtime_root/native-prefixes && ! -L $runtime_root/native-prefixes &&
           -f $metadata && ! -L $metadata ]] || {
            printf 'CADD prefix metadata is incomplete or linked: %s\n' "$metadata" >&2; return 69;
        }
        # Plain, immutable path data: never source/eval metadata or activate in
        # the parent shell. The ready receipt independently pins this prefix.
        {
            IFS= read -r prefix && ! IFS= read -r extra && [[ -z $extra ]]
        } < "$metadata" || { printf 'Invalid CADD prefix metadata: %s\n' "$metadata" >&2; return 69; }
        [[ $prefix == /* && $prefix == */envs/"$name" && $prefix != *$'\r'* &&
           /$prefix/ != */../* && /$prefix/ != */./* ]] || {
            printf 'Invalid CADD prefix path: %s\n' "$metadata" >&2; return 69;
        }
    fi
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
