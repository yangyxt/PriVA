#!/usr/bin/env bash
# Synthetic wiring tests only: no mature variant tools or real caches invoked.
# Usage: bash tests/nextflow/run_stub.sh /absolute/path/to/nextflow /absolute/run-root
set -euo pipefail
nextflow_bin=${1:?Pass the pinned Nextflow executable}
run_root=${2:?Pass a new durable test run root}
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
mkdir -p "${run_root}/validation" "${run_root}/logs" "${run_root}/work/driver-stub"
fixture="${run_root}/validation/nextflow-fixture"
case_root="${run_root}/validation/nextflow-stub"
python3 "${repo_root}/tests/nextflow/make_stub_fixture.py" --outdir "${fixture}"
cd "${run_root}/work/driver-stub"
common=(run "${repo_root}" -profile test -stub-run -ansi-log false
    -work-dir "${run_root}/work/nextflow-stub"
    --input "${fixture}/cohort.vcf.gz" --ped "${fixture}/family.ped"
    --family FAM --assembly hg19 --legacy_config "${fixture}/legacy.yaml"
    --resource_manifest "${fixture}/resources.json" --software_id stub)
run_case() {
    local name=$1
    shift
    if ! "${nextflow_bin}" -log "${run_root}/logs/${name}.log" "${common[@]}" "$@" > "${run_root}/logs/${name}.stdout.log" 2>&1; then
        tail -80 "${run_root}/logs/${name}.stdout.log"
        return 1
    fi
}
run_case nextflow-stub-raw --outdir "${case_root}/raw"
run_case nextflow-stub-resume --outdir "${case_root}/raw" -resume
python3 - "${case_root}/raw" <<'PY'
import csv
from pathlib import Path
import sys
root = Path(sys.argv[1])
traces = sorted((root / 'pipeline_info').glob('trace-*.tsv'))
assert len(traces) >= 2, 'Fresh and resume traces are required'
rows = list(csv.DictReader(traces[-1].open(), delimiter='\t'))
assert len(rows) == 4 and all(row['status'] == 'CACHED' for row in rows), rows
assert len(list(root.rglob('*.filtered.acmg.tsv'))) == 1
assert len(list(root.rglob('*.filtered.vcf.gz'))) == 1
assert len(list(root.rglob('*.provenance.json'))) == 3
print('Fresh raw stub and all-four-task cached resume passed')
PY
run_case nextflow-stub-annotated --mode annotated --cadd "${fixture}/cohort.anno.cadd.tsv" --outdir "${case_root}/annotated"
# Same output directory with a different environment identity must fail before
# cached tasks are considered; this guards against silently retaining old tables.
if "${nextflow_bin}" -log "${run_root}/logs/nextflow-stub-reject-output.log" "${common[@]}" \
    --software_id different-environment --outdir "${case_root}/raw" -resume \
    > "${run_root}/logs/nextflow-stub-reject-output.stdout.log" 2>&1; then
    echo 'ERROR: mismatched output identity was accepted' >&2
    exit 1
fi
python3 - "${run_root}" <<'PY'
from pathlib import Path
import json
import sys
root = Path(sys.argv[1])
failure = (root / 'logs/nextflow-stub-reject-output.stdout.log').read_text()
assert 'different inputs/source/resources/settings' in failure, failure
contract = json.loads(next((root / 'validation/nextflow-stub/raw').rglob('contract.json')).read_text())
assert contract['samples'] == ['child', 'dad', 'mom', 'other']
assert contract['family_samples'] == ['child', 'dad', 'mom']
assert contract['stub'] is True
assert len(list((root / 'validation/nextflow-stub/annotated').rglob('*.filtered.acmg.tsv'))) == 1
print('Annotated stub, preserved cohort context, and stale-output rejection passed')
PY
