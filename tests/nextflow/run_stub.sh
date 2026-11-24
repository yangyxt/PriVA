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
        # Tasks of the ANNOTATE_NEW sub-workflow carry its name as a prefix.
        name = row['name'].split(' ')[0].split(':')[-1]
        if name == 'ANNOTATE_STEP':
            # "ANNOTATE_STEP (BASE:STEP)" or, for a VEP piece, "(BASE:vep:PIECE)".
            name += ':' + row['name'].split('(', 1)[1].rstrip(')').split(':')[1]
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
# VEP and CADD run on per-chromosome pieces (one piece in the stub), gathered back.
assert first.count('SPLIT_CHROMOSOMES') == 2 and first.count('GATHER_PIECES') == 2, first
assert 'SPLIT_CHROMOSOMES' not in second and 'GATHER_PIECES' not in second, second
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
# Phenotype gene bonus. It is on for each family whose proband has HPO terms:
# with terms for both, FAM and OTHER each get their genes ranked (a header-only
# stub table here) and their ACMG task the ranking. Terms of a family that is not
# in the PED file are reported and left out. With terms for FAM only, FAM gets
# its ranking and OTHER ranks as before, with a warning naming it.
# The ranking program is priva's own (priva hpo-gene-ranking); its reference
# folder is the installed hpo_gene_ranking_reference_dir, which main.nf takes
# from the resource manifest (here a stub manifest that lists only it).
mkdir -p "${fixture}/hpo_reference/raw" "${fixture}/hpo_reference/prepared" "${fixture}/hpo_reference/hpo_graph"
printf 'sample_id\thpo_term\nFAM\tHP:0001250\nFAM\tGlobal developmental delay\nOTHER\tHP:0001263\nABSENT\tHP:0000252\n' > "${fixture}/hpo_terms.tsv"
printf 'sample_id\thpo_term\nFAM\tHP:0001250\n' > "${fixture}/hpo_terms_partial.tsv"
python3 - "${fixture}" <<'PY'
import json
from pathlib import Path
import sys
fixture = Path(sys.argv[1])
reference = (fixture / 'hpo_reference').resolve()
record = {'path': str(reference), 'resolved_path': str(reference), 'type': 'directory', 'size': 4096, 'mtime_ns': 0}
(fixture / 'resources_phenotype.json').write_text(json.dumps(
    {'bundle_id': 'stub-only', 'stub': True, 'resources': {'hpo_gene_ranking_reference_dir': record}}) + '\n')
PY
phenotype_shared=("${shared[@]/#${fixture}\/resources.json/${fixture}/resources_phenotype.json}")
phenotype=(--ped "${fixture}/family.ped")
run_with nextflow-stub-phenotype "${phenotype_shared[@]}" "${phenotype[@]}" \
    --phenotype_hpo_terms "${fixture}/hpo_terms.tsv" --outdir "${case_root}/phenotype"
run_with nextflow-stub-phenotype-partial "${phenotype_shared[@]}" "${phenotype[@]}" \
    --phenotype_hpo_terms "${fixture}/hpo_terms_partial.tsv" --outdir "${case_root}/phenotype-partial"
# Without the installed reference in the manifest, and with the removed
# program and reference options, the run is refused at launch.
if "${nextflow_bin}" -log "${run_root}/logs/nextflow-stub-phenotype-no-reference.log" "${shared[@]}" "${phenotype[@]}" \
    --phenotype_hpo_terms "${fixture}/hpo_terms.tsv" --outdir "${case_root}/phenotype-no-reference" \
    > "${run_root}/logs/nextflow-stub-phenotype-no-reference.stdout.log" 2>&1; then
    echo 'ERROR: --phenotype_hpo_terms was accepted without an installed reference' >&2
    exit 1
fi
grep -q 'install step F10 HPO_GENE_RANKING' "${run_root}/logs/nextflow-stub-phenotype-no-reference.stdout.log"
if "${nextflow_bin}" -log "${run_root}/logs/nextflow-stub-phenotype-removed.log" "${phenotype_shared[@]}" "${phenotype[@]}" \
    --phenotype_hpo_terms "${fixture}/hpo_terms.tsv" --hpo_gene_ranking_binary "$(type -P true)" \
    --outdir "${case_root}/phenotype-removed" \
    > "${run_root}/logs/nextflow-stub-phenotype-removed.stdout.log" 2>&1; then
    echo 'ERROR: the removed --hpo_gene_ranking_binary was accepted' >&2
    exit 1
fi
grep -q 'hpo_gene_ranking_reference were removed' "${run_root}/logs/nextflow-stub-phenotype-removed.stdout.log"
python3 - "${case_root}" "${run_root}/logs" <<'PY'
import csv
import json
from pathlib import Path
import sys
cases, logs = Path(sys.argv[1]), Path(sys.argv[2])
def tasks(root):
    trace = sorted((root / 'pipeline_info').glob('trace-*.tsv'))[-1]
    return sorted(row['name'] for row in csv.DictReader(trace.open(), delimiter='\t'))
def provenance(root, family):
    return json.loads((root / 'pipeline_info' / f'prioritization-acmg.{family}.provenance.json').read_text())
# Every family has terms: the bonus is on for both.
root, log = cases / 'phenotype', (logs / 'nextflow-stub-phenotype.log').read_text()
names = tasks(root)
assert [t for t in names if t.startswith('RANK_PHENOTYPE_GENES')] == ['RANK_PHENOTYPE_GENES (FAM)', 'RANK_PHENOTYPE_GENES (OTHER)'], names
assert [t for t in names if t.startswith('PRIORITIZE_ACMG')] == ['PRIORITIZE_ACMG (FAM)', 'PRIORITIZE_ACMG (OTHER)'], names
for family in ('FAM', 'OTHER'):
    record = provenance(root, family)
    assert record['phenotype_gene_ranking']['name'] == f'{family}.gene_ranked.tsv', record
    assert sorted(p.name for p in (root / 'phenotype_gene_ranking' / family).iterdir()) == sorted(
        [f'{family}.{t}.tsv' for t in ('frontier_gene_evidence', 'frontier_p_terms', 'gene_ranked', 'gene_score_contributions', 'gene_summary', 'p_pairs')] + ['run_summary.tsv'])
assert 'left out: ABSENT' in log, 'the family missing from the PED file was not reported'
assert 'Phenotype gene bonus is on for 2 of 2 families' in log
# Terms for FAM only: the bonus is on for FAM and off for OTHER.
root, log = cases / 'phenotype-partial', (logs / 'nextflow-stub-phenotype-partial.log').read_text()
names = tasks(root)
assert [t for t in names if t.startswith('RANK_PHENOTYPE_GENES')] == ['RANK_PHENOTYPE_GENES (FAM)'], names
assert [t for t in names if t.startswith('PRIORITIZE_ACMG')] == ['PRIORITIZE_ACMG (FAM)', 'PRIORITIZE_ACMG (OTHER)'], names
assert provenance(root, 'FAM')['phenotype_gene_ranking']['name'] == 'FAM.gene_ranked.tsv'
assert 'phenotype_gene_ranking' not in provenance(root, 'OTHER')
assert sorted(p.name for p in (root / 'phenotype_gene_ranking').iterdir()) == ['FAM']
assert 'Phenotype gene bonus is on for 1 of 2 families' in log
assert 'Phenotype gene bonus is off for families without HPO terms: OTHER' in log, 'the family without terms was not reported'
print('Phenotype gene bonus wiring passed (on for each family with terms, off for one without, refused without the installed reference or with the removed options)')
PY
# SpliceAI2 scoring. With the model's weights in the configuration, each family
# gets a SCORE_SPLICEAI2 task between the combined table and the steps that read
# it; without them (every case above) there is none. The annotation hub's folder
# goes to the task.
cp "${fixture}/legacy.yaml" "${fixture}/legacy_spliceai2.yaml"
printf 'spliceai2_model_dir: "%s"  # weights\nspliceai2_python: /usr/bin/python3\n' "${fixture}/spliceai2_models" >> "${fixture}/legacy_spliceai2.yaml"
spliceai2_shared=("${shared[@]/#${fixture}\/legacy.yaml/${fixture}/legacy_spliceai2.yaml}")
run_with nextflow-stub-spliceai2 "${spliceai2_shared[@]}" --ped "${fixture}/family.ped" \
    --annotation_hub "${case_root}/hub-spliceai2" --outdir "${case_root}/spliceai2"
python3 - "${case_root}" <<'PY'
import csv
import json
from pathlib import Path
import sys
root = Path(sys.argv[1]) / 'spliceai2'
trace = sorted((root / 'pipeline_info').glob('trace-*.tsv'))[-1]
names = sorted(row['name'].split(':')[-1] for row in csv.DictReader(trace.open(), delimiter='\t'))
assert [n for n in names if n.startswith('SCORE_SPLICEAI2')] == ['SCORE_SPLICEAI2 (FAM)', 'SCORE_SPLICEAI2 (OTHER)'], names
for step in ('COMBINE', 'SPLICING', 'UTR', 'JOIN', 'ACMG'):
    assert len([n for n in names if n.startswith(f'PRIORITIZE_{step} ')]) == 2, names
for family in ('FAM', 'OTHER'):
    record = json.loads((root / 'pipeline_info' / f'prioritization-spliceai2.{family}.provenance.json').read_text())
    assert record['phase'] == 'prioritization-spliceai2' and record['family'] == family, record
assert len(list(root.rglob('*.filtered.acmg.tsv'))) == 2
print('SpliceAI2 scoring wiring passed (one task per family with the weights configured)')
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
