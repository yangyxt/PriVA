#!/usr/bin/env bash
#SBATCH --job-name=PriVA_Licensed_Cache_Import
#SBATCH --partition=copy,long
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --no-requeue
# Do not submit until the acquisition owner explicitly confirms final paths.
# Twelve CPUs cover the original sort's default eight workers plus its pipes.
set -euo pipefail
[[ $# == 8 ]] || { echo 'Usage: sbatch import_cluster_licensed.sbatch LOG REPO REVISION SOURCE_DEST CONFIG_TEMPLATE SEPARATE_CONFIG CONFIRMED_READY_JSON STATE_DIR' >&2; exit 64; }
log_pattern=$1; repo=$2; revision=$3; source_dest=$4
template=$5; config=$6; confirmation=$7; state_dir=$8
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
    export TMPDIR=/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/tmp/licensed-cache-import
    mkdir -p "$TMPDIR"
    exec 9> "$state_dir/licensed-cache-import.install.lock"
    flock -n 9
    python - "$source_dest" "$revision" "$template" "$config" "$confirmation" "$state_dir" <<'PY'
import json
import os
from pathlib import Path
import subprocess
import sys
import yaml
source, revision, template, config_path, confirmation_path, state = map(Path, sys.argv[1:])
revision = str(revision)
sys.path.insert(0, str(source / 'scripts'))
from deploy_cluster_resources import checksum, write_json

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def assembly_name(value):
    require(value in ('hg19', 'GRCh37', 'hg38', 'GRCh38'), 'Unsupported assembly')
    return {'GRCh37': 'hg19', 'GRCh38': 'hg38'}.get(value, value)

confirmation = json.loads(confirmation_path.read_text())
require(confirmation.get('status') == 'downloads_complete' and confirmation.get('confirmed_by'),
        'Explicit acquisition-owner completion confirmation is required before accessing cache files')
settings = yaml.safe_load(template.read_text())
assembly = assembly_name(settings['assembly'])
require(assembly_name(confirmation['assembly']) == assembly, 'Confirmed files belong to another assembly')
require(not config_path.exists(), 'Refusing to replace an existing import config')
receipt_path = state / f'licensed-caches-{assembly}.ready.json'
require(not receipt_path.exists(), 'An import receipt already exists; review before deliberate rerun')
require(set(confirmation.get('files', {})) == {'spliceai_snv', 'spliceai_indel', 'primateai'},
        'All three matching-assembly resources are required')
cache_root = Path('/stornext/snfs170/next-gen/scratch/yangyxt/caches')
paths, identities, before = {}, {}, {}
for role, item in confirmation['files'].items():
    require(role in ('spliceai_snv', 'spliceai_indel', 'primateai'), 'Unknown resource role')
    root = cache_root / ('PrimateAI_v0_2' if role == 'primateai' else 'SpliceAI_v1')
    path = Path(item['path']).resolve()
    require(path.is_relative_to(root.resolve()), f'Unexpected shared-cache root: {path}')
    require(path.is_file() and path.stat().st_size > 0, f'Missing/empty confirmed input: {path}')
    require('.partial' not in path.name, 'An incomplete file cannot be imported')
    if role.startswith('spliceai_'):
        variant = role.removeprefix('spliceai_')
        require(path.name == f'spliceai_scores.raw.{variant}.{assembly}.vcf.gz', 'Unexpected SpliceAI assembly/score type')
    else:
        raw = 'PrimateAI_scores_v0.2' + ('_hg38' if assembly == 'hg38' else '')
        grch = 'GRCh37' if assembly == 'hg19' else 'GRCh38'
        allowed = (raw + '.tsv.gz', raw + '.sorted.tsv.bgz', f'PrimateAI_scores_v0.2_{grch}_sorted.tsv.bgz')
        require(path.name in allowed, 'Unexpected PrimateAI version/assembly filename')
    if 'bytes' in item:
        require(path.stat().st_size == item['bytes'], f'Acquisition byte-size mismatch: {path}')
    identity = checksum(path)
    if item.get('sha256'):
        require(identity == item['sha256'], f'Acquisition checksum mismatch: {path}')
    subprocess.run(['gzip', '-t', str(path)], check=True)
    paths[role] = path
    identities[role] = {'path': str(path), 'bytes': path.stat().st_size, 'sha256': identity}
    before[path] = (path.stat().st_size, path.stat().st_mtime_ns)
    if role.startswith('spliceai_') or path.name.endswith('.bgz'):
        index = Path(str(path) + '.tbi')
        require(index.is_file() and index.stat().st_size > 0, f'Missing matching index: {index}')
        identities[role]['index_sha256'] = checksum(index)
        before[index] = (index.stat().st_size, index.stat().st_mtime_ns)
        contigs = subprocess.check_output(['tabix', '-l', str(path)], text=True).splitlines()
        require(contigs, f'No indexed chromosomes: {path}')
        subprocess.run(['tabix', str(path), contigs[0] + ':1-1'], stdout=subprocess.DEVNULL, check=True)
primate = paths['primateai']
if primate.name.endswith('.tsv.gz'):
    formatted = primate.with_name(primate.name.replace('.tsv.gz', '.sorted.tsv.bgz'))
    require(not formatted.exists() and not Path(str(formatted) + '.tbi').exists(),
            'Formatting would overwrite an existing shared cache; supply that verified formatted file instead')
else:
    formatted = primate
for key in ('spliceai_snv_prescore', 'spliceai_indel_prescore', 'primateai_prescore'):
    settings[key] = ''  # Force original prompts in this separate registration config.
with config_path.open('x') as stream:
    yaml.safe_dump(settings, stream, sort_keys=False)
commands = []
for function, answers in (
    ('SpliceAI_install', ['yes', str(paths['spliceai_snv']), str(paths['spliceai_indel'])]),
    ('PrimateAI_install', ['yes', str(primate)]),
):
    command = ['bash', '-o', 'pipefail', str(source / 'scripts/install_utils.sh'), function,
               str(config_path), settings['vep_plugins_cachedir'], settings['vep_plugins_dir']]
    print('Actual installer import:', json.dumps(command), flush=True)
    subprocess.run(command, input='\n'.join(answers) + '\n', text=True, check=True)
    commands.append(command)
result = yaml.safe_load(config_path.read_text())
expected = {'spliceai_snv_prescore': paths['spliceai_snv'], 'spliceai_indel_prescore': paths['spliceai_indel'],
            'primateai_prescore': formatted}
for key, path in expected.items():
    require(Path(result[key]).resolve() == path, f'Installer registered unexpected {key}')
    require(path.is_file() and Path(str(path) + '.tbi').is_file(), f'Missing imported data/index: {path}')
    contigs = subprocess.check_output(['tabix', '-l', str(path)], text=True).splitlines()
    require(contigs, f'No indexed chromosomes: {path}')
    subprocess.run(['tabix', str(path), contigs[0] + ':1-1'], stdout=subprocess.DEVNULL, check=True)
for path, previous in before.items():
    require((path.stat().st_size, path.stat().st_mtime_ns) == previous, f'Original shared input changed: {path}')
write_json(receipt_path, {
    'status': 'licensed_inputs_imported', 'assembly': assembly, 'installer_revision': revision,
    'confirmation_sha256': checksum(confirmation_path), 'confirmed_by': confirmation['confirmed_by'],
    'source_files': identities, 'formatted_primateai': str(formatted),
    'formatted_primateai_sha256': checksum(formatted), 'formatted_primateai_index_sha256': checksum(Path(str(formatted) + '.tbi')),
    'config': str(config_path), 'config_sha256': checksum(config_path), 'registered_paths': {key: str(value) for key, value in expected.items()},
    'commands': commands, 'slurm_job_id': os.environ['SLURM_JOB_ID'],
    'scope': 'Original installer prompts register shared SpliceAI paths and retain original PrimateAI formatting order. Licensed originals preserved. Main config has not been changed; VEP/family scientific qualification remains separate.'})
print('Published', receipt_path, flush=True)
PY
}
( trap 'result=$?; printf "[%s] Exit status: %s\n" "$(date -Is)" "$result"; exit "$result"' EXIT; main ) > "$log_file" 2>&1
