#!/usr/bin/env bash
#SBATCH --job-name=PriVA_gnomAD_hg19
#SBATCH --partition=copy,long
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --no-requeue
# Source downloads must finish first: submit with --dependency=afterok:SOURCE_JOB.
set -euo pipefail
[[ $# == 7 ]] || { echo 'Usage: sbatch --dependency=afterok:JOB build_cluster_gnomad.sbatch LOG REPO REVISION SOURCE_DEST CONFIG CONFIG_SHA256 STATE_DIR' >&2; exit 64; }
log_pattern=$1; repo=$2; revision=$3; source_dest=$4
config=$5; config_sha=$6; state_dir=$7
log_file=${log_pattern//%j/${SLURM_JOB_ID:?}}
main() {
    printf '[%s] Job %s on %s; installer revision %s\n' "$(date -Is)" "$SLURM_JOB_ID" "$(hostname)" "$revision"
    [[ $revision =~ ^[0-9a-f]{40}$ ]] || return 65
    [[ ! -e $source_dest ]] || { echo 'Refusing an existing source snapshot destination'; return 66; }
    mkdir -p "$source_dest"
    git -C "$repo" archive "$revision" | tar -xf - -C "$source_dest"
    printf '%s\n' "$revision" > "$source_dest/REFERENCE_COMMIT"
    manager=/stornext/snfs170/next-gen/scratch/yangyxt/mamba/bin/micromamba
    activation=$("$manager" shell activate --shell bash --name priva_acmg)
    eval "$activation"
    export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    export TMPDIR=/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/tmp/gnomad-liftover
    mkdir -p "$TMPDIR"
    # This is also the source downloader's lock: no simultaneous replacement.
    exec 9> "$state_dir/gnomad-source.install.lock"
    flock -n 9
    python - "$source_dest" "$revision" "$config" "$config_sha" "$state_dir" "$SLURM_CPUS_PER_TASK" <<'PY'
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import yaml

source, revision, config_path, config_sha, state, threads = sys.argv[1:]
source, config_path, state = Path(source), Path(config_path), Path(state)

def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

require(sha256(config_path) == config_sha, 'Frozen config changed after submission')
config = yaml.safe_load(config_path.read_text())
require(config['assembly'] in ('hg19', 'GRCh37'), 'This derived resource job is for GRCh37 only')
destination = Path(config['gnomad_vcf_dir']).resolve()
manifest_path = source / 'envs/cluster-gnomad.v4.1.json'
manifest = json.loads(manifest_path.read_text())
input_receipt = json.loads((state / 'gnomad-source.ready.json').read_text())
require(input_receipt['status'] == 'source_inputs_deployed', 'Source download did not complete')
require(input_receipt['assembly'] == 'GRCh38' and input_receipt['release'] == '4.1', 'Wrong source assembly/release')
require(input_receipt['manifest_sha256'] == sha256(manifest_path), 'Source manifest changed')
require(input_receipt['files'] == manifest['files'] and len(manifest['files']) == 48, 'Incomplete source manifest')
source_stats = {}
for entry in manifest['files']:
    path = destination / entry['path']
    require(path.is_file() and path.stat().st_size == entry['bytes'], f'Missing/truncated source: {path}')
    source_stats[str(path)] = (path.stat().st_size, path.stat().st_mtime_ns)

reference = Path(config['ref_genome']).resolve()
reference_receipt = json.loads((state / 'reference.ready.json').read_text())
require(reference_receipt['status'] == 'deployed', 'Reference deployment incomplete')
require(Path(reference_receipt['reference']).resolve() == reference, 'Reference receipt describes a different path')
require(sha256(reference) == reference_receipt['sha256'], 'Reference changed after verification')
require(Path(str(reference) + '.fai').is_file(), 'Missing reference index')
chain = source / 'data/liftover/hg38ToHg19.over.chain.gz'
require(chain.is_file() and chain.stat().st_size > 1024, 'Missing/LFS-pointer liftover chain')
subprocess.run(['gzip', '-t', str(chain)], check=True)
final_receipt = state / 'gnomad-hg19.ready.json'
require(not final_receipt.exists(), 'Derived receipt already exists; review before deliberate rerun')
chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
outputs = [destination / f'gnomad.joint.v4.1.sites.hg19.chr{c}.vcf.gz' for c in chromosomes]
require(not any(str(p) in source_stats for p in outputs), 'Input/output paths overlap')
command = ['bash', str(source / 'scripts/install_utils.sh'), 'gnomAD_liftover',
           str(destination / 'gnomad.joint.v4.1.sites.chr1.vcf.bgz'),
           str(chain), str(reference), threads]
print('Actual cache build command:', json.dumps(command), flush=True)
subprocess.run(command, check=True)
output_metadata = []
for chromosome, path in zip(chromosomes, outputs):
    require(path.is_file() and path.stat().st_size > 0, f'Missing derived VCF: {path}')
    index = Path(str(path) + '.tbi')
    require(index.is_file() and index.stat().st_size > 0, f'Missing derived index: {index}')
    contigs = subprocess.check_output(['tabix', '-l', str(path)], text=True).splitlines()
    require(contigs == [f'chr{chromosome}'], f'Unexpected indexed contigs in {path}: {contigs}')
    subprocess.run(['bcftools', 'view', '-h', str(path)], stdout=subprocess.DEVNULL, check=True)
    subprocess.run(['tabix', str(path), f'chr{chromosome}:1-1'], stdout=subprocess.DEVNULL, check=True)
    output_metadata.append({'path': str(path), 'bytes': path.stat().st_size,
                            'index_bytes': index.stat().st_size, 'indexed_contigs': contigs})
for name, before in source_stats.items():
    path = Path(name)
    require((path.stat().st_size, path.stat().st_mtime_ns) == before, f'Source file changed during liftover: {path}')
require(sha256(config_path) == config_sha, 'Frozen config changed during execution')
receipt = {
    'status': 'derived_resource_built', 'assembly': 'GRCh37', 'release': '4.1',
    'installer_revision': revision, 'scientific_baseline_revision': '09223ef5a8b0ceb4c6d61055c6eefe62b51ba903',
    'installer_sha256': sha256(source / 'scripts/install_utils.sh'),
    'chain_sha256': sha256(chain), 'reference_sha256': reference_receipt['sha256'],
    'config': str(config_path), 'config_sha256': config_sha,
    'input_manifest_sha256': input_receipt['manifest_sha256'],
    'source_download_job': input_receipt['slurm_job_id'], 'slurm_job_id': os.environ['SLURM_JOB_ID'],
    'command': command, 'outputs': output_metadata,
    'scope': 'All 24 hg19 liftover VCFs/indexes built with unchanged baseline per-chromosome transform; headers, indexed contigs and interval lookup checked. Family prioritization and calibration parity remain unqualified.',
    'source_integrity': 'Publisher MD5 validated by completed source downloader; exact manifest and file sizes checked before this job, size/mtime unchanged afterward. Source files were not rehashed during liftover.'}
temporary = final_receipt.with_name(final_receipt.name + '.partial')
temporary.write_text(json.dumps(receipt, indent=2) + '\n')
temporary.replace(final_receipt)
print('Published', final_receipt, flush=True)
PY
}
( trap 'result=$?; printf "[%s] Exit status: %s\n" "$(date -Is)" "$result"; exit "$result"' EXIT; main ) > "$log_file" 2>&1
