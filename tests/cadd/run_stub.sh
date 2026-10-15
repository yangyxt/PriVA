#!/usr/bin/env bash
# Usage: run_stub.sh PINNED_NEXTFLOW NEW_DURABLE_RUN_ROOT
set -euo pipefail
nextflow_bin=${1:?Pass pinned Nextflow executable}
run_root=${2:?Pass new durable run root}
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
mkdir -p "$run_root/logs" "$run_root/work/driver"
python3 "$repo_root/tests/cadd/make_stub_fixture.py" --output "$run_root/fixture"
cd "$run_root/work/driver"
common=(run "$repo_root/cadd.nf" -c "$repo_root/conf/cadd.config" -profile test -stub-run -ansi-log false
  -work-dir "$run_root/work/tasks" --cadd_annotation_memory '512 MB' --cadd_calculation_memory '512 MB')
run_case() {
    local name=$1; shift
    local input_mode=reference-vcf
    [[ "$name" != shared-* ]] || input_mode=shared-vep
    if ! "$nextflow_bin" -log "$run_root/logs/$name.log" "${common[@]}" --cadd_annotation_source "$input_mode" "$@" > "$run_root/logs/$name.stdout.log" 2>&1; then
        tail -n 65 "$run_root/logs/$name.stdout.log"
        return 1
    fi
}
run_case full-GRCh37 --cadd_mode full --cadd_profile "$run_root/fixture/GRCh37.json" \
    --assembly GRCh37 --input "$run_root/fixture/input.vcf.gz" --outdir "$run_root/full-GRCh37"
run_case resume --cadd_mode full --cadd_profile "$run_root/fixture/GRCh37.json" \
    --assembly GRCh37 --input "$run_root/fixture/input.vcf.gz" --outdir "$run_root/full-GRCh37" -resume
run_case full-GRCh38 --cadd_mode full --cadd_profile "$run_root/fixture/GRCh38.json" \
    --assembly GRCh38 --input "$run_root/fixture/input.vcf.gz" --outdir "$run_root/full-GRCh38"
run_case annotate-only --cadd_mode annotate --cadd_profile "$run_root/fixture/GRCh37.json" \
    --assembly GRCh37 --input "$run_root/fixture/input.vcf.gz" --outdir "$run_root/annotate-only"
run_case calculate-only --cadd_mode calculate --cadd_profile "$run_root/fixture/GRCh37.json" \
    --assembly GRCh37 --annotation_bundle "$run_root/annotate-only/annotations/annotation_bundle" --outdir "$run_root/calculate-only"
run_case shared-GRCh37 --cadd_mode full --cadd_profile "$run_root/fixture/GRCh37.json" \
    --assembly GRCh37 --shared_vep_bundle "$run_root/fixture/shared-GRCh37" --outdir "$run_root/shared-GRCh37"
run_case shared-GRCh38 --cadd_mode full --cadd_profile "$run_root/fixture/GRCh38.json" \
    --assembly GRCh38 --shared_vep_bundle "$run_root/fixture/shared-GRCh38" --outdir "$run_root/shared-GRCh38"
python3 - "$run_root" <<'PY'
import csv,json,sys
from pathlib import Path
root=Path(sys.argv[1])
traces=sorted((root/'full-GRCh37/pipeline_info').glob('trace-*.tsv'))
rows=list(csv.DictReader(traces[-1].open(),delimiter='\t'))
assert len(rows)==2 and all(r['status']=='CACHED' for r in rows),rows
for name in ['full-GRCh37','full-GRCh38','calculate-only','shared-GRCh37','shared-GRCh38']:
    receipt=json.loads((root/name/'calculation/calculation_bundle/manifest.json').read_text())
    assert receipt['stub'] is True and receipt['row_count']==1
assert (root/'annotate-only/annotations/annotation_bundle/manifest.json').is_file()
assert not (root/'annotate-only/calculation').exists()
print('Both assembly shared-VEP/reference routes, two-task cached resume, annotation-only and calculation-only passed')
PY
if "$nextflow_bin" -log "$run_root/logs/reject-assembly.log" "${common[@]}" --cadd_annotation_source reference-vcf --cadd_mode full \
    --cadd_profile "$run_root/fixture/GRCh37.json" --assembly GRCh38 \
    --input "$run_root/fixture/input.vcf.gz" --outdir "$run_root/reject-assembly" \
    > "$run_root/logs/reject-assembly.stdout.log" 2>&1; then
    echo 'ERROR: assembly mismatch accepted' >&2; exit 1
fi
if "$nextflow_bin" -log "$run_root/logs/reject-stale-output.log" "${common[@]}" --cadd_annotation_source reference-vcf --cadd_mode annotate \
    --cadd_profile "$run_root/fixture/GRCh37.json" --assembly GRCh37 \
    --input "$run_root/fixture/input.vcf.gz" --outdir "$run_root/full-GRCh37" -resume \
    > "$run_root/logs/reject-stale-output.stdout.log" 2>&1; then
    echo 'ERROR: stale output identity accepted' >&2; exit 1
fi
python3 - "$run_root" <<'PY'
from pathlib import Path
import sys
root=Path(sys.argv[1])
assert 'assembly differs' in (root/'logs/reject-assembly.stdout.log').read_text()
assert 'different inputs/profile/source' in (root/'logs/reject-stale-output.stdout.log').read_text()
print('Launch-time assembly and stale-output rejection passed')
PY
