#!/usr/bin/env bash
#SBATCH --job-name=PriVA_SpliceVault_hg19
#SBATCH --partition=copy,long
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=7-00:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --no-requeue
set -euo pipefail
[[ $# == 7 ]] || { echo 'Usage: sbatch build_cluster_splicevault.sbatch LOG REPO REVISION SOURCE_DEST CONFIG CONFIG_SHA256 STATE_DIR' >&2; exit 64; }
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
    exec 9> "$state_dir/splicevault.install.lock"
    flock -n 9
    python - "$source_dest" "$revision" "$config" "$config_sha" "$state_dir" <<'PY'
import json
import os
from pathlib import Path
import subprocess
import sys
import yaml
source, revision, config_path, config_sha, state = sys.argv[1:]
source, config_path, state = Path(source), Path(config_path), Path(state)
sys.path.insert(0, str(source / 'scripts'))
from deploy_cluster_resources import checksum, download, write_json

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

require(checksum(config_path) == config_sha, 'Frozen config changed after submission')
config = yaml.safe_load(config_path.read_text())
url = 'https://ftp.ensembl.org/pub/release-113/variation/SpliceVault/SpliceVault_data_GRCh38.tsv.gz'
require(config['splicevault_url'] == url, 'Expected pinned release113 source URL')
require(config['assembly'] in ('hg19', 'GRCh37'), 'Expected GRCh37 output profile')
plugin_cache = Path(config['vep_plugins_cachedir']).resolve()
folder = plugin_cache / 'SpliceVault'
folder.mkdir(parents=True, exist_ok=True)
final = folder / 'SpliceVault_data_GRCh37.tsv.gz'
require(Path(config['splicevault_prescore']).resolve() == final, 'Unexpected output destination')
require(not final.exists(), 'Unqualified existing output; review before deliberate rerun')
receipt_path = state / 'splicevault-hg19.ready.json'
require(not receipt_path.exists(), 'A completion receipt already exists')
reference = Path(config['ref_genome']).resolve()
reference_receipt = json.loads((state / 'reference.ready.json').read_text())
require(reference_receipt['status'] == 'deployed', 'Reference is not deployed')
require(Path(reference_receipt['reference']).resolve() == reference, 'Wrong reference receipt')
require(checksum(reference) == reference_receipt['sha256'], 'Reference identity changed')
chain = source / 'data/liftover/hg38ToHg19.over.chain.gz'
require(chain.is_file() and chain.stat().st_size > 1024, 'Missing/LFS-pointer chain')
subprocess.run(['gzip', '-t', str(chain)], check=True)
input_path = download(url, folder / 'SpliceVault_data_GRCh38.tsv.gz', expected_size=883907365)
published = [57027, 863191]
actual = list(map(int, subprocess.check_output(['sum', '-r', str(input_path)], text=True).split()[:2]))
require(actual == published, 'SpliceVault input differs from Ensembl CHECKSUMS')
input_sha = checksum(input_path)
temporary = Path(config['tmp_dir'])
temporary.mkdir(parents=True, exist_ok=True)
command = ['bash', str(source / 'scripts/install_utils.sh'), 'SpliceVault_install', str(config_path), str(plugin_cache)]
print('Actual original conversion command:', json.dumps(command), flush=True)
subprocess.run(command, check=True, cwd=temporary)
require(checksum(input_path) == input_sha, 'Input changed during installation')
require(final.is_file() and final.stat().st_size > 0, 'No completed hg19 table')
require(Path(str(final) + '.tbi').is_file(), 'No completed hg19 index')
contigs = subprocess.check_output(['tabix', '-l', str(final)], text=True).splitlines()
require('chr1' in contigs or '1' in contigs, 'Expected chromosome1 in derived index')
subprocess.run(['tabix', str(final), contigs[0] + ':1-1'], stdout=subprocess.DEVNULL, check=True)
write_json(receipt_path, {
    'status': 'derived_resource_built', 'assembly': 'GRCh37', 'installer_revision': revision,
    'scientific_baseline_revision': '09223ef5a8b0ceb4c6d61055c6eefe62b51ba903',
    'source_url': url, 'source_bytes': input_path.stat().st_size, 'source_sha256': input_sha,
    'publisher_bsd_checksum': published,
    'checksum_source': 'https://ftp.ensembl.org/pub/release-113/variation/SpliceVault/CHECKSUMS',
    'reference_sha256': reference_receipt['sha256'], 'chain_sha256': checksum(chain),
    'converter_sha256': checksum(source / 'scripts/SpliceVault_tsv_vcf_conversion.py'),
    'output': str(final), 'output_sha256': checksum(final), 'index_sha256': checksum(Path(str(final) + '.tbi')),
    'config_initial_sha256': config_sha, 'config_final_sha256': checksum(config_path),
    'command': command, 'indexed_contigs': contigs, 'slurm_job_id': os.environ['SLURM_JOB_ID'],
    'scope': 'Original PriVA GRCh38-to-GRCh37 conversion retained. This does not substitute the newer direct hg19 distribution or establish family/ACMG parity.'})
print('Published', receipt_path, flush=True)
PY
}
( trap 'result=$?; printf "[%s] Exit status: %s\n" "$(date -Is)" "$result"; exit "$result"' EXIT; main ) > "$log_file" 2>&1
