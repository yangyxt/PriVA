#!/usr/bin/env bash
# Synthetic wiring test of cadd.nf: -stub-run writes a header-only CADD table;
# no CADD data, environment or model is used.
# Usage: bash tests/nextflow/run_cadd_stub.sh /absolute/path/to/nextflow /absolute/run-root
set -euo pipefail
nextflow_bin=${1:?Pass the pinned Nextflow executable}
run_root=${2:?Pass a new durable test run root}
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
mkdir -p "${run_root}/logs" "${run_root}/work/driver-stub" "${run_root}/input"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' | gzip -c > "${run_root}/input/sample.vcf.gz"
cd "${run_root}/work/driver-stub"
if ! "${nextflow_bin}" -log "${run_root}/logs/cadd-stub.log" run "${repo_root}/cadd.nf" \
        -c "${repo_root}/conf/cadd.config" -profile test -stub-run -ansi-log false \
        -work-dir "${run_root}/work/cadd-stub" --input "${run_root}/input/sample.vcf.gz" \
        --assembly GRCh38 --cadd_base_dir "${run_root}/absent-cadd" --outdir "${run_root}/cadd" \
        > "${run_root}/logs/cadd-stub.stdout.log" 2>&1; then
    tail -80 "${run_root}/logs/cadd-stub.stdout.log"
    exit 1
fi
python3 - "${repo_root}" "${run_root}/cadd/sample.cadd.tsv" <<'PY'
import sys
sys.path.insert(0, sys.argv[1] + "/bin")
import priva_task
priva_task.validate_cadd(sys.argv[2])
print("cadd.nf stub run wrote a CADD table that main.nf accepts")
PY
