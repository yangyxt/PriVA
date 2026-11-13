#!/usr/bin/env bash
# Synthetic wiring tests only: no mature variant tools or real caches invoked.
# Usage: bash tests/nextflow/run_stub.sh /absolute/path/to/nextflow /absolute/run-root [PRIVA]
set -euo pipefail
nextflow_bin=${1:?Pass the pinned Nextflow executable}
run_root=${2:?Pass a new durable test run root}
repo_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)
priva_bin=${3:-${repo_root}/rust/target/release/priva}
mkdir -p "${run_root}/validation" "${run_root}/logs" "${run_root}/work/driver-stub"
fixture="${run_root}/validation/nextflow-fixture"
case_root="${run_root}/validation/nextflow-stub"
python3 "${repo_root}/tests/nextflow/make_stub_fixture.py" --outdir "${fixture}"
cd "${run_root}/work/driver-stub"
# Without --ped or --family; the cases add them. The fixture PED has two
# families, FAM (a trio) and OTHER (a singleton).
shared=(run "${repo_root}" -profile test -stub-run -ansi-log false
    -work-dir "${run_root}/work/nextflow-stub"
    --input "${fixture}/cohort.vcf.gz" --assembly hg19 --legacy_config "${fixture}/legacy.yaml"
    --resource_manifest "${fixture}/resources.json" --software_id stub
    --priva_binary "${priva_bin}")
common=("${shared[@]}" --ped "${fixture}/family.ped" --family FAM)
run_with() {
    local name=$1
    shift
    if ! "${nextflow_bin}" -log "${run_root}/logs/${name}.log" "$@" > "${run_root}/logs/${name}.stdout.log" 2>&1; then
        tail -80 "${run_root}/logs/${name}.stdout.log"
        return 1
    fi
}
run_case() {
    local name=$1
    shift
    run_with "${name}" "${common[@]}" "$@"
}
# One annotation route: preprocessing; gnomAD, ClinVar, control, VEP and CADD
# side by side; a merge of their INFO fields. A resume reuses every task.
run_case nextflow-stub-raw --outdir "${case_root}/raw"
run_case nextflow-stub-resume --outdir "${case_root}/raw" -resume
python3 - "${case_root}/raw" <<'PY'
import csv
from pathlib import Path
import sys
root = Path(sys.argv[1])
traces = sorted((root / 'pipeline_info').glob('trace-*.tsv'))
assert len(traces) >= 2, 'Fresh and resume traces are required'
fresh = sorted(row['name'].split(' ')[0] for row in csv.DictReader(traces[-2].open(), delimiter='\t'))
PRIORITIZE = ['PRIORITIZE_ACMG', 'PRIORITIZE_COMBINE', 'PRIORITIZE_JOIN', 'PRIORITIZE_SPLICING', 'PRIORITIZE_UTR']
assert fresh == sorted(['ANNOTATE_CADD', 'ANNOTATE_MERGE', 'ANNOTATE_PREPROCESS'] + ['ANNOTATE_STEP'] * 4 + [
    'FILTER_FAMILY', 'VALIDATE_INPUTS'] + PRIORITIZE), fresh
resumed = list(csv.DictReader(traces[-1].open(), delimiter='\t'))
assert len(resumed) == 14 and all(row['status'] == 'CACHED' for row in resumed), resumed
assert sorted(p.name for p in (root / 'annotation').iterdir()) == [
    'cohort.anno.cadd.tsv', 'cohort.anno.vcf.gz', 'cohort.anno.vcf.gz.tbi']
names = sorted(p.name for p in (root / 'pipeline_info').glob('annotation-*.provenance.json'))
assert names == [f'annotation-{step}.provenance.json' for step in
                 ('cadd', 'clinvar', 'control', 'gnomad', 'merge', 'preprocess', 'vep')], names
assert len(list(root.rglob('*.filtered.acmg.tsv'))) == 1
assert len(list(root.rglob('*.filtered.vcf.gz'))) == 1
print('Annotation as one task per step, and its cached resume, passed')
PY
# Annotation hub: the first run stores every step's results; a second run, in
# a new session, finds every variant stored and runs no annotation step except
# the control cohort step, which has no hub folder without a control VCF.
hub_args=(--annotation_hub "${case_root}/hub")
run_case nextflow-stub-hub-first "${hub_args[@]}" --outdir "${case_root}/hub-first"
run_case nextflow-stub-hub-second "${hub_args[@]}" --outdir "${case_root}/hub-second"
python3 - "${case_root}" <<'PY'
import csv
import gzip
from pathlib import Path
import sys
root = Path(sys.argv[1])
def tasks(case):
    trace = sorted((root / case / 'pipeline_info').glob('trace-*.tsv'))[-1]
    names = []
    for row in csv.DictReader(trace.open(), delimiter='\t'):
        name = row['name'].split(' ')[0]
        if name == 'ANNOTATE_STEP':
            name += ':' + row['name'].split(':')[-1].rstrip(')')
        names.append(name)
    return sorted(names)
first, second = tasks('hub-first'), tasks('hub-second')
annotate = lambda names: [n for n in names if n.startswith('ANNOTATE_STEP') or n == 'ANNOTATE_CADD']
assert annotate(first) == ['ANNOTATE_CADD', 'ANNOTATE_STEP:clinvar', 'ANNOTATE_STEP:control',
                           'ANNOTATE_STEP:gnomad', 'ANNOTATE_STEP:vep'], first
assert annotate(second) == ['ANNOTATE_STEP:control'], second
for names in (first, second):
    assert names.count('HUB_LOOKUP') == 5 and names.count('HUB_FILL') == 5, names
assert first.count('HUB_ADD') == 4 and second.count('HUB_ADD') == 0, (first, second)
read = lambda case: gzip.open(root / case / 'annotation' / 'cohort.anno.vcf.gz', 'rt').read()
assert read('hub-first') == read('hub-second') == read('raw'), 'annotated VCFs differ'
stored = sorted({p.relative_to(root / 'hub').parts[0] for p in (root / 'hub').rglob('part-*.parquet')})
assert stored == ['cadd', 'clinvar', 'gnomad', 'vep'], stored
print('Annotation hub: first run stores, second run annotates nothing new, outputs identical')
PY
# Whole cohort: one validation and one annotation, then each family's
# filtration and prioritization tasks.
run_with nextflow-stub-cohort "${shared[@]}" --ped "${fixture}/family.ped" --outdir "${case_root}/cohort"
# No PED: one filtration and one prioritization without a family.
run_with nextflow-stub-no-pedigree "${shared[@]}" --outdir "${case_root}/no-pedigree"
python3 - "${case_root}" <<'PY'
import csv
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
def names(case, pattern):
    return sorted(p.name for p in (root / case).rglob(pattern))
ANNOTATE = ['ANNOTATE_CADD', 'ANNOTATE_MERGE', 'ANNOTATE_PREPROCESS'] + ['ANNOTATE_STEP'] * 4
PRIORITIZE = ['PRIORITIZE_ACMG', 'PRIORITIZE_COMBINE', 'PRIORITIZE_JOIN', 'PRIORITIZE_SPLICING', 'PRIORITIZE_UTR']
STEPS = ['acmg', 'combine', 'join', 'splicing', 'utr']
ANNOTATION_RECORDS = [f'annotation-{s}.provenance.json' for s in ('cadd', 'clinvar', 'control', 'gnomad', 'merge', 'preprocess', 'vep')]
def tasks(case):
    trace = sorted((root / case / 'pipeline_info').glob('trace-*.tsv'))[-1]
    return sorted(row['name'].split(' ')[0] for row in csv.DictReader(trace.open(), delimiter='\t'))
assert tasks('cohort') == sorted(ANNOTATE + ['FILTER_FAMILY', 'FILTER_FAMILY', 'VALIDATE_INPUTS'] + PRIORITIZE * 2), tasks('cohort')
assert names('cohort', '*.filtered.acmg.tsv') == ['cohort.anno.FAM.filtered.acmg.tsv', 'cohort.anno.OTHER.filtered.acmg.tsv'], names('cohort', '*.acmg.tsv')
assert names('cohort', '*.provenance.json') == sorted(ANNOTATION_RECORDS + ['filtration.FAM.provenance.json', 'filtration.OTHER.provenance.json'] + [f'prioritization-{s}.{f}.provenance.json' for s in STEPS for f in ('FAM', 'OTHER')]), names('cohort', '*.provenance.json')
contract = json.loads(next((root / 'cohort').rglob('contract.json')).read_text())
assert contract['families'] == ['OTHER', 'FAM'], contract['families']
assert contract['family_samples'] == {'OTHER': ['other'], 'FAM': ['child', 'dad', 'mom']}, contract['family_samples']
assert tasks('no-pedigree') == sorted(ANNOTATE + ['FILTER_FAMILY', 'VALIDATE_INPUTS'] + PRIORITIZE), tasks('no-pedigree')
assert names('no-pedigree', '*.filtered.acmg.tsv') == ['cohort.anno.filtered.acmg.tsv'], names('no-pedigree', '*.acmg.tsv')
assert names('no-pedigree', '*.filtered.vcf.gz') == ['cohort.anno.filtered.vcf.gz']
contract = json.loads(next((root / 'no-pedigree').rglob('contract.json')).read_text())
assert contract['families'] == [] and contract['ped_sha256'] is None, contract
print('Cohort stub with two families and no-pedigree stub passed')
PY
# The removed annotated mode and old hub options are refused with a message
# that points to the annotation hub.
if "${nextflow_bin}" -log "${run_root}/logs/nextflow-stub-removed-mode.log" "${common[@]}" \
    --mode annotated --cadd "${fixture}/cohort.anno.cadd.tsv" --outdir "${case_root}/removed-mode" \
    > "${run_root}/logs/nextflow-stub-removed-mode.stdout.log" 2>&1; then
    echo 'ERROR: the removed annotated mode was accepted' >&2
    exit 1
fi
# Same output directory with a different environment identity must fail before
# cached tasks are considered; this guards against silently retaining old tables.
if "${nextflow_bin}" -log "${run_root}/logs/nextflow-stub-reject-output.log" "${common[@]}" \
    --software_id different-environment --outdir "${case_root}/raw" -resume \
    > "${run_root}/logs/nextflow-stub-reject-output.stdout.log" 2>&1; then
    echo 'ERROR: mismatched output identity was accepted' >&2
    exit 1
fi
# An option value written after a space is read by Nextflow as a separate
# flag; the launch must refuse it before any task instead of passing 'true'.
if "${nextflow_bin}" -log "${run_root}/logs/nextflow-stub-flag-value.log" "${common[@]}" \
    --outdir "${case_root}/flag-value" --slurm_cluster_options --constraint=NAME \
    > "${run_root}/logs/nextflow-stub-flag-value.stdout.log" 2>&1; then
    echo 'ERROR: an option value read as a flag was accepted' >&2
    exit 1
fi
python3 - "${run_root}" <<'PY'
from pathlib import Path
import json
import sys
root = Path(sys.argv[1])
removed = (root / 'logs/nextflow-stub-removed-mode.stdout.log').read_text()
assert 'annotated mode was removed' in removed, removed
failure = (root / 'logs/nextflow-stub-reject-output.stdout.log').read_text()
assert 'different inputs/source/resources/settings' in failure, failure
flag = (root / 'logs/nextflow-stub-flag-value.stdout.log').read_text()
assert '--slurm_cluster_options needs a value' in flag, flag
contract = json.loads(next((root / 'validation/nextflow-stub/raw').rglob('contract.json')).read_text())
assert contract['samples'] == ['child', 'dad', 'mom', 'other']
assert contract['families'] == ['FAM']
assert contract['family_samples'] == {'FAM': ['child', 'dad', 'mom']}
assert contract['stub'] is True
print('Removed-option refusal, preserved cohort context, stale-output and flag-value rejection passed')
PY
