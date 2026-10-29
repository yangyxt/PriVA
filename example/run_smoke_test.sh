#!/usr/bin/env bash
# Compatibility wrapper; named run/prepare/validate entry points live in smoke_test.py.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

case "${1:-}" in
    run|prepare|validate) exec python "${SCRIPT_DIR}/smoke_test.py" "$@" ;;
    -h|--help) exec python "${SCRIPT_DIR}/smoke_test.py" --help ;;
    --*) exec python "${SCRIPT_DIR}/smoke_test.py" run "$@" ;;
    hg19|hg38|GRCh37|GRCh38) assembly=$1; shift ;;
    *) echo "Usage: $0 <hg19|hg38|GRCh37|GRCh38> [config] [cores] [--output-dir DIR]" >&2; exit 2 ;;
esac
args=(run --assembly "$assembly")
if [[ $# -gt 0 && $1 != --* ]]; then
    args+=(--config "$1")
    shift
fi
cores=${PRIVA_SMOKE_CORES:-2}
if [[ $# -gt 0 && $1 != --* ]]; then
    cores=$1
    shift
fi
exec python "${SCRIPT_DIR}/smoke_test.py" "${args[@]}" --cores "$cores" "$@"
