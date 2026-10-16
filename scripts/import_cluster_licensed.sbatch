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
[[ $# == 8 || $# == 9 ]] || { echo 'Usage: sbatch import_cluster_licensed.sbatch LOG REPO REVISION SOURCE_DEST CONFIG_TEMPLATE SEPARATE_CONFIG CONFIRMED_READY_JSON STATE_DIR [all|spliceai|primateai]' >&2; exit 64; }
log_pattern=$1; repo=$2; revision=$3; source_dest=$4
template=$5; config=$6; confirmation=$7; state_dir=$8; resource_set=${9:-all}
case "$resource_set" in all|spliceai|primateai) ;; *) echo "Unknown resource set: $resource_set" >&2; exit 64 ;; esac
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
    python - "$source_dest" "$revision" "$template" "$config" "$confirmation" "$state_dir" "$resource_set" <<'PY'
import json
import os
from pathlib import Path
import subprocess
import sys
import yaml
CACHE_ROOT = Path('/stornext/snfs170/next-gen/scratch/yangyxt/caches')
RESOURCE_ROLES = {'spliceai': ('spliceai_snv', 'spliceai_indel'), 'primateai': ('primateai',)}
ROLE_KEYS = {'spliceai_snv': 'spliceai_snv_prescore', 'spliceai_indel': 'spliceai_indel_prescore',
             'primateai': 'primateai_prescore'}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def assembly_name(value):
    require(value in ('hg19', 'GRCh37', 'hg38', 'GRCh38'), 'Unsupported assembly')
    return {'GRCh37': 'hg19', 'GRCh38': 'hg38'}.get(value, value)


def progress(message):
    print(message, flush=True)


def file_identity(path):
    stat = path.stat()
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def require_stable_inputs(identities):
    for path, previous in identities.items():
        require(file_identity(path) == previous, f'Original shared input changed: {path}')


def import_resources(source, revision, template, config_path, confirmation_path, state,
                     resource_set='all', *, cache_root=CACHE_ROOT):
    # cache_root is injectable only for isolated fixture tests, never a CLI override.
    sys.path.insert(0, str(source / 'scripts'))
    from deploy_cluster_resources import checksum, write_json
    require(resource_set in ('all', 'spliceai', 'primateai'), 'Unknown resource set')
    selected = list(RESOURCE_ROLES) if resource_set == 'all' else [resource_set]
    required_roles = {role for resource in selected for role in RESOURCE_ROLES[resource]}
    selected_keys = {ROLE_KEYS[role] for role in required_roles}
    confirmation = json.loads(confirmation_path.read_text())
    require(confirmation.get('status') == 'downloads_complete' and confirmation.get('confirmed_by'),
            'Explicit acquisition-owner completion confirmation is required before accessing cache files')
    if 'resources' in confirmation:
        require(confirmation['resources'] == selected,
                'Confirmed resource selection does not match the requested resource set')
    settings = yaml.safe_load(template.read_text())
    original_settings = settings.copy()
    assembly = assembly_name(settings['assembly'])
    require(assembly_name(confirmation['assembly']) == assembly, 'Confirmed files belong to another assembly')
    require(not config_path.exists(), 'Refusing to replace an existing import config')
    suffix = '' if resource_set == 'all' else '-' + resource_set
    receipt_path = state / f'licensed-caches-{assembly}{suffix}.ready.json'
    require(not receipt_path.exists(), 'An import receipt already exists; review before deliberate rerun')
    role_message = ('All three matching-assembly resources are required' if resource_set == 'all' else
                    'Exactly the selected matching-assembly resource roles are required')
    require(set(confirmation.get('files', {})) == required_roles, role_message)
    progress(f'Selected resources: {selected}; assembly: {assembly}; only these resources will be imported')
    paths, identities, before = {}, {}, {}
    for role in sorted(required_roles):
        item = confirmation['files'][role]
        root = cache_root / ('PrimateAI_v0_2' if role == 'primateai' else 'SpliceAI_v1')
        path = Path(item['path'])
        require(path.is_absolute(), 'Confirmed cache paths must be absolute')
        path = path.resolve()
        require(path.is_relative_to(root.resolve()), f'Unexpected shared-cache root: {path}')
        require('.partial' not in path.name, 'An incomplete file cannot be imported')
        if role.startswith('spliceai_'):
            variant = role.removeprefix('spliceai_')
            require(path.name == f'spliceai_scores.raw.{variant}.{assembly}.vcf.gz',
                    'Unexpected SpliceAI assembly/score type')
        else:
            raw = 'PrimateAI_scores_v0.2' + ('_hg38' if assembly == 'hg38' else '')
            grch = 'GRCh37' if assembly == 'hg19' else 'GRCh38'
            allowed = (raw + '.tsv.gz', raw + '.sorted.tsv.bgz', f'PrimateAI_scores_v0.2_{grch}_sorted.tsv.bgz')
            require(path.name in allowed, 'Unexpected PrimateAI version/assembly filename')
        require(path.is_file() and path.stat().st_size > 0, f'Missing/empty confirmed input: {path}')
        before[path] = file_identity(path)  # Capture before any long checksum/gzip operation.
        if 'bytes' in item:
            require(before[path][2] == item['bytes'], f'Acquisition byte-size mismatch: {path}')
        progress(f'{role}: SHA256 start: {path}')
        identity = checksum(path)
        if item.get('sha256'):
            require(identity == item['sha256'], f'Acquisition checksum mismatch: {path}')
        progress(f'{role}: SHA256 complete; gzip integrity check start')
        subprocess.run(['gzip', '-t', str(path)], check=True)
        progress(f'{role}: gzip integrity check complete')
        paths[role] = path
        identities[role] = {'path': str(path), 'bytes': before[path][2], 'sha256': identity}
        if role.startswith('spliceai_') or path.name.endswith('.bgz'):
            index = Path(str(path) + '.tbi')
            require(index.resolve().is_relative_to(root.resolve()), f'Unexpected shared-cache index root: {index}')
            require(index.is_file() and index.stat().st_size > 0, f'Missing matching index: {index}')
            before[index] = file_identity(index)
            progress(f'{role}: index SHA256 start: {index}')
            identities[role]['index_sha256'] = checksum(index)
            progress(f'{role}: indexed chromosome/query checks start')
            contigs = subprocess.check_output(['tabix', '-l', str(path)], text=True).splitlines()
            require(contigs, f'No indexed chromosomes: {path}')
            subprocess.run(['tabix', str(path), contigs[0] + ':1-1'], stdout=subprocess.DEVNULL, check=True)
            progress(f'{role}: indexed checks complete')
    expected = {ROLE_KEYS[role]: path for role, path in paths.items()}
    invocations = []
    if 'spliceai' in selected:
        invocations.append(('SpliceAI_install', ['yes', str(paths['spliceai_snv']), str(paths['spliceai_indel'])]))
    if 'primateai' in selected:
        primate = paths['primateai']
        if primate.name.endswith('.tsv.gz'):
            formatted = primate.with_name(primate.name.replace('.tsv.gz', '.sorted.tsv.bgz'))
            require(not formatted.exists() and not Path(str(formatted) + '.tbi').exists(),
                    'Formatting would overwrite an existing shared cache; supply that verified formatted file instead')
        else:
            formatted = primate
        expected['primateai_prescore'] = formatted
        invocations.append(('PrimateAI_install', ['yes', str(primate)]))
    require_stable_inputs(before)  # Reject verification-time changes before writing the config.
    for key in selected_keys:
        settings[key] = ''  # Force only the selected original prompts in a separate config.
    with config_path.open('x') as stream:
        yaml.safe_dump(settings, stream, sort_keys=False)
    commands = []
    for function, answers in invocations:
        command = ['bash', '-o', 'pipefail', str(source / 'scripts/install_utils.sh'), function,
                   str(config_path), settings['vep_plugins_cachedir'], settings['vep_plugins_dir']]
        progress('Actual installer import: ' + json.dumps(command))
        subprocess.run(command, input='\n'.join(answers) + '\n', text=True, check=True)
        commands.append(command)
    result = yaml.safe_load(config_path.read_text())
    require({k: v for k, v in result.items() if k not in selected_keys} ==
            {k: v for k, v in original_settings.items() if k not in selected_keys},
            'Installer modified unselected configuration keys')
    for key, path in expected.items():
        require(Path(result[key]).resolve() == path, f'Installer registered unexpected {key}')
        require(path.is_file() and Path(str(path) + '.tbi').is_file(), f'Missing imported data/index: {path}')
        progress(f'{key}: post-import indexed checks start')
        contigs = subprocess.check_output(['tabix', '-l', str(path)], text=True).splitlines()
        require(contigs, f'No indexed chromosomes: {path}')
        subprocess.run(['tabix', str(path), contigs[0] + ':1-1'], stdout=subprocess.DEVNULL, check=True)
        progress(f'{key}: post-import indexed checks complete')
    receipt = {
        'status': 'licensed_inputs_imported' if resource_set == 'all' else 'selected_licensed_inputs_imported',
        'assembly': assembly, 'installer_revision': revision,
        'selected_resources': selected,
        'resource_status': {name: 'imported' if name in selected else 'not_selected_unverified'
                            for name in RESOURCE_ROLES},
        'confirmation_sha256': checksum(confirmation_path), 'confirmed_by': confirmation['confirmed_by'],
        'source_files': identities,
        'config': str(config_path), 'config_sha256': checksum(config_path),
        'registered_paths': {key: str(value) for key, value in expected.items()},
        'commands': commands, 'slurm_job_id': os.environ['SLURM_JOB_ID'],
        'scope': 'Original installer prompts import only the selected resources. Unselected resources are unverified. Licensed originals preserved. Main config has not been changed; VEP/family scientific qualification remains separate.'}
    if 'primateai' in selected:
        progress('primateai: formatted output and index SHA256 start')
        receipt.update(formatted_primateai=str(formatted), formatted_primateai_sha256=checksum(formatted),
                       formatted_primateai_index_sha256=checksum(Path(str(formatted) + '.tbi')))
    require_stable_inputs(before)  # Also reject changes during installer/output verification.
    write_json(receipt_path, receipt)
    progress('Published ' + str(receipt_path))
    return receipt_path


if __name__ == '__main__':
    source, revision, template, config_path, confirmation_path, state = sys.argv[1:7]
    import_resources(Path(source), revision, Path(template), Path(config_path), Path(confirmation_path),
                     Path(state), sys.argv[7] if len(sys.argv) > 7 else 'all')

PY
}
( trap 'result=$?; printf "[%s] Exit status: %s\n" "$(date -Is)" "$result"; exit "$result"' EXIT; main ) > "$log_file" 2>&1
