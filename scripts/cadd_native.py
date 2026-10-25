#!/usr/bin/env python3
"""Deploy and call the original-Python CADD stages without a second VEP.

install-profile records installed resources, never a transcript-pool decision or
scientific qualification. preflight/run require an explicit qualified profile.
run consumes only the qualified true-novel bundle prepared by the caller after
the original ordered publisher-prescore split; it does not replace that split.
Complete feature/score bundles are retained independently of legacy hub views.
"""
import argparse
import copy
import hashlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
import cadd_task as T
import cadd_shared_vep as S
_runtime_spec = importlib.util.spec_from_file_location('priva_native_runtime_integrity', ROOT / 'scripts/cadd_runtime.py')
R = importlib.util.module_from_spec(_runtime_spec); _runtime_spec.loader.exec_module(R)

FLAGS = dict(scientific_qualification=False, production_activation=False,
             full_public_benchmark=False, family_qualification=False)
PRESCORE_POLICY = 'upstream-v1.7.2-inclAnno-v1'
PRESCORE_CHECKS = dict(prepared_allele_mapping='exact', publisher_import='exact',
    novel_feature_joins='exact', source_order='preserved', missingness='preserved')


def require(ok, message):
    if not ok:
        raise T.ContractError(message)


def encoded(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()


def publish(path, value, verify_only=False):
    """Immutable, atomic publication; retain unrelated or interrupted bytes."""
    path = Path(path)
    require(not path.is_symlink(), 'Refusing linked native CADD metadata: ' + str(path))
    if path.exists():
        require(path.read_bytes() == encoded(value), 'Existing native CADD metadata differs; retained: ' + str(path))
        return
    require(not verify_only, 'Native CADD deployment is incomplete: ' + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    require(not path.parent.is_symlink(), 'Refusing linked native CADD metadata directory')
    with tempfile.NamedTemporaryFile(prefix='.' + path.name + '.', dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name); stream.write(encoded(value))
    try:
        try: os.link(temporary, path)
        except FileExistsError:
            require(not path.is_symlink() and path.read_bytes() == encoded(value),
                    'Concurrent native CADD metadata differs; retained: ' + str(path))
    finally:
        temporary.unlink(missing_ok=True)


def config_values(config):
    require(isinstance(config, dict), 'CADD configuration must be a mapping')
    assembly = T.ALIASES.get(config.get('assembly'))
    release = str(config.get('vep_installed_version', ''))
    require(assembly and re.fullmatch(r'[1-9][0-9]*', release), 'Select assembly and numeric VEP release')
    root = config.get('cache_root')
    require(isinstance(root, str) and Path(root).is_absolute() and Path(root) != Path('/'),
            'Set a dedicated absolute cache_root')
    return assembly, release, Path(root).absolute()


def file_pin(path):
    return {'path': str(path), 'sha256': T.digest(path)}


def verify_pin(pin):
    path = Path(pin['path'])
    require(path.is_absolute() and re.fullmatch('[0-9a-f]{64}', pin['sha256'])
            and T.digest(path) == pin['sha256'], 'Native CADD input proof artifact changed: ' + str(path))
    return path


def json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def prepared_keys(path):
    """Original VCF2vepVCF/extract_scored artifacts contain only five columns."""
    with S.text_file(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            require(len(fields) >= 5, 'Original prepared/novel VCF record is incomplete')
            yield fields[0], fields[1], fields[3], fields[4]


def verify_prescoring_profile(profile):
    block = profile['prescoring']; qualification = block.get('qualification')
    require(qualification, 'Publisher-prescore/novel profile policy is unqualified')
    evidence = T.read_json(verify_pin(qualification))
    expected = dict(schema_version=1, status='qualified', policy_id=PRESCORE_POLICY,
        profile_id=profile['profile_id'], resource_manifest_sha256=profile['contract']['resource_manifest_sha256'],
        stub=False, checks=PRESCORE_CHECKS)
    require(all(evidence.get(key) == value for key, value in expected.items()),
            'Publisher-prescore policy qualification does not match the profile')
    verify_pin(evidence['evidence'])


def verify_novel_input(profile, profile_path, bundle):
    """New native endpoint requires the caller's closed, ordered prescore split."""
    verify_prescoring_profile(profile)
    pipeline = T.read_json(Path(bundle) / 'manifest.json').get('pipeline_provenance', {})
    stage = 'upstream-prescore-novel-subset-v1'
    expected = dict(stage=stage, profile_sha256=T.digest(profile_path),
                    prescore_inputs_sha256=json_hash(profile['prescoring']))
    require(all(pipeline.get(key) == value for key, value in expected.items()),
            'Native CADD requires a profile-bound original publisher-prescore split proof')
    proof = T.read_json(verify_pin(pipeline['prescore_split']))
    require(all(proof.get(key) == value for key, value in
        dict(schema_version=1, status='completed', policy_id=PRESCORE_POLICY, **expected).items()),
        'Original publisher-prescore split is incomplete or belongs to another profile')
    for name in ('full_source', 'mapping'):
        require(proof[name] == pipeline[name], 'Publisher split/source mapping differs')
        verify_pin(proof[name])
    verify_pin(proof['prepared']); previous = proof['prepared']
    steps = proof['steps']; files = profile['prescoring']['files']
    require(len(steps) == len(files), 'Publisher split omitted an ordered cache file')
    for step, item in zip(steps, files):
        require(step['publisher_path'] == item['path'] and step['input'] == previous,
                'Publisher split cache order or remaining-input chain differs')
        for name in ('input', 'found', 'remaining'): verify_pin(step[name])
        previous = step['remaining']
    require(proof['novel'] == previous, 'Publisher split final novel artifact differs')
    novel = verify_pin(proof['novel'])
    keys = proof['novel_keys']
    require(isinstance(keys, dict) and type(keys.get('count')) is int and keys['count'] > 0,
            'A pinned streamed novel-key artifact and positive count are required')
    key_path = verify_pin(keys)
    count = 0
    with key_path.open() as handle:
        def recorded_keys():
            for line in handle:
                fields = line.rstrip('\r\n').split('\t')
                require(len(fields) == 4, 'Novel-key artifact must have four columns')
                yield tuple(fields)
        for view_key, novel_key, expected_key in itertools.zip_longest(
                S.variant_keys(Path(bundle) / 'vep.vcf.gz'), prepared_keys(novel), recorded_keys()):
            require(view_key is not None and view_key == novel_key == expected_key,
                    'Native CADD view is not the exact true-novel allele sequence')
            count += 1
    require(count == keys['count'] and pipeline.get('novel_keys_sha256') == keys['sha256'],
            'Native CADD true-novel key count or checksum differs')


def prescoring(source, assembly):
    """Inventory the literal upstream inclAnno branch, without reading tables."""
    import yaml
    configuration = source / ('config/config_' + assembly + '_v1.7.yml')
    settings = yaml.safe_load(configuration.read_text())
    folder = source / settings['PrescoredFolder']
    expected = source / ('data/prescored/' + assembly + '_v1.7/incl_anno')
    require(folder == expected, 'Original inclAnno prescore folder differs')
    tables = sorted(folder.glob('*.tsv.gz'), key=lambda p: os.fsencode(p.name))
    require(tables, 'Original publisher prescores are not installed')
    files = []
    for table in tables:
        index = Path(str(table) + '.tbi')
        require(table.is_file() and index.is_file(), 'Publisher prescore table/index is missing')
        files.append({'path': str(table), 'data': T.record(table), 'index': T.record(index)})
    return {'policy_id': PRESCORE_POLICY, 'folder': str(folder),
            'configuration': file_pin(configuration),
            'prepare': file_pin(source / 'src/scripts/VCF2vepVCF.py'),
            'extract': file_pin(source / 'src/scripts/extract_scored.py'),
            'snakefile': file_pin(source / 'Snakefile'), 'files': files, 'qualification': None}


def install_profile(config, receipt_path, verify_only=False):
    """Bind completed native components; perform no download, probe or inference."""
    assembly, release, cache_root = config_values(config)
    receipt_path = Path(receipt_path).absolute(); record = T.read_json(receipt_path)
    require(record.get('backend') == 'native' and record.get('selection') == 'both'
            and record.get('assembly') == assembly and record.get('source_commit') == T.PIN
            and record.get('status') == 'installed_not_scientifically_qualified',
            'A completed original-Python native both-stage installation is required')
    source = Path(record['source'])
    if config.get('cadd_base_dir'):
        require(source == Path(config['cadd_base_dir']).absolute(), 'Configured CADD source differs')
    candidate = record['selected_transcript_candidate']
    require(candidate['assembly'] == assembly and candidate['vep_release'] == release
            and T.digest(candidate['path']) == candidate['sha256'], 'Selected transcript candidate differs')
    candidate_record = T.read_json(candidate['path'])
    require(candidate_record.get('stage') == 'cadd-selected-transcript-candidate'
            and candidate_record.get('scientific_qualified') is False
            and candidate_record.get('cadd_enabled') is False, 'Transcript candidate is not installation-only')
    wrappers = record['runtime']['wrappers']
    routes = wrappers['entry_points']; runtime_root = Path(wrappers['path']).parent.parent
    calculation = record['calculation_runtime']
    require(calculation['implementation'] == 'all_four_calculation_stages_upstream_python'
            and calculation['sha256'] == T.PYTHON_CALCULATION_SHA256
            and T.digest(calculation['program']) == T.PYTHON_CALCULATION_SHA256,
            'Installed calculator must be the pinned original Python program')
    template = T.read_json(ROOT / 'conf/cadd' / (assembly + '.reference.json'))
    profile = copy.deepcopy(template)
    profile.update(profile_kind='adapted', software_id='cadd-native-' + record['fingerprint'],
                   deployment_status='installed_not_scientifically_qualified', **FLAGS)
    profile['profile_id'] = 'cadd-native-' + assembly.lower() + '-vep' + release + '-' + record['fingerprint'][:16]
    annotation = source / ('data/annotations/' + assembly + '_v1.7')
    env = dict(CADD_SOURCE_ROOT=str(source), CADD_PYTHON2=routes['cadd-python2'],
        CADD_ESM_PYTHON=routes['cadd-esm-python'], CADD_ESM_MODEL_DIR=record['assets']['esm_hub']['path'],
        CADD_ESM_MODELS=','.join('esm1v_t33_650M_UR90S_' + str(i) for i in range(1, 6)), CADD_ESM_BATCH_SIZE='1',
        CADD_ESM_TRANSCRIPTS=candidate['environment_candidates']['CADD_ESM_TRANSCRIPTS'],
        CADD_REGSEQ_PYTHON=routes['cadd-regseq-python'],
        CADD_REGSEQ_MODEL=str(annotation / 'regseq/Hyperopt400InclNegatives.json'),
        CADD_REGSEQ_WEIGHTS=str(annotation / 'regseq/Hyperopt400InclNegatives.h5'),
        CADD_REGSEQ_REFERENCE=str(annotation / 'regseq/reference.fa'),
        CADD_REGSEQ_GENOME=str(annotation / 'regseq/reference.fa.genome'),
        CADD_REFERENCE_CONFIG=str(source / ('config/references_' + assembly + '_v1.7.cfg')))
    # Use the already installed shared HTS executables, never a new CADD VEP env.
    import shutil
    for key, name in (('CADD_BGZIP', 'bgzip'), ('CADD_TABIX', 'tabix')):
        program = shutil.which(name)
        require(program, 'Install ' + name + ' before native CADD profile deployment')
        env[key] = str(Path(program).absolute())
    if assembly == 'GRCh38':
        env.update(CADD_MMSPLICE_PYTHON=routes['cadd-mmsplice-python'],
            CADD_MMSPLICE_GTF=candidate['environment_candidates']['CADD_MMSPLICE_GTF'],
            CADD_REFERENCE_FASTA=str(annotation / 'reference/reference.fa'))
    profile['annotation']['env'] = env
    profile['annotation']['shared_vep'] = dict(vep_release=release, transcript_release=release,
        esm_transcript_release=release, mmsplice_transcript_release=release if assembly == 'GRCh38' else None,
        adapter_id=S.NATIVE_ADAPTER, transcript_policy=None, selection_context=None, qualification=None)
    annotation_roots = [str(annotation), env['CADD_REFERENCE_CONFIG'], env['CADD_ESM_MODEL_DIR'], env['CADD_ESM_TRANSCRIPTS']]
    if assembly == 'GRCh38': annotation_roots.append(env['CADD_MMSPLICE_GTF'])
    profile['prescoring'] = prescoring(source, assembly)
    annotation_roots += [profile['prescoring']['folder']]
    if config.get('ref_genome'):
        annotation_roots += [str(Path(config['ref_genome']).absolute()), str(Path(str(config['ref_genome']) + '.fai').absolute())]
    profile['annotation']['resource_roots'] = annotation_roots
    profile['calculation']['env'] = calculation['environment']
    profile['calculation']['program'] = calculation['program']
    profile['calculation']['program_sha256'] = calculation['sha256']
    profile['calculation']['verify_files'] = [{'path': item['path'], 'sha256': item['sha256']}
        for item in record['calculation_assets'][:3]]
    profile['calculation']['resource_roots'] = [item['path'] for item in record['calculation_assets'][:3]]
    profile['installation'] = dict(receipt=file_pin(receipt_path), cache_root=str(cache_root),
        wrappers=wrappers, transcript_candidate=candidate, source_identity=T.source_identity('native'),
        esm_hub=record['assets']['esm_hub'],
        vep_fasta=T.record(config['ref_genome']) if config.get('ref_genome') else None,
        vep_fasta_index=T.record(str(config['ref_genome']) + '.fai') if config.get('ref_genome') else None)
    manifest = {'schema_version': 1, 'assembly': assembly,
                'resources': {stage: T.inventory(profile[stage]['resource_roots']) for stage in ('annotation', 'calculation')}}
    directory = cache_root / 'runtime/cadd-1.7.2/profiles'
    for parent in (directory, directory.parent, directory.parent.parent):
        require(not parent.is_symlink(), 'Refusing linked native CADD deployment directory: ' + str(parent))
    manifest_sha = hashlib.sha256(encoded(manifest)).hexdigest()
    manifest_path = directory / ('resources.' + assembly + '.' + manifest_sha[:16] + '.json')
    profile['resource_manifest_path'] = str(manifest_path)
    profile['contract']['resource_manifest_sha256'] = manifest_sha
    checksum = hashlib.sha256(encoded(profile)).hexdigest()
    path = directory / (assembly + '.vep' + release + '.' + checksum[:16] + '.installed.json')
    publish(manifest_path, manifest, verify_only); publish(path, profile, verify_only)
    return {'path': str(path), 'sha256': checksum, 'status': profile['deployment_status'], **FLAGS}


def check(config, profile_path):
    """Read-only component integrity; success explicitly does not permit execution."""
    assembly, release, cache_root = config_values(config)
    profile = T.load_profile(profile_path, assembly)
    binding = profile['installation']; pin = binding['receipt']
    require(binding['cache_root'] == str(cache_root) and T.digest(pin['path']) == pin['sha256'],
            'Installed native CADD/configuration binding differs')
    record = T.read_json(pin['path'])
    require(record['backend'] == 'native' and record['selection'] == 'both'
            and profile['annotation']['shared_vep']['vep_release'] == release,
            'Installed native CADD assembly/release/backend differs')
    require(binding['source_identity'] == T.source_identity('native'), 'Native CADD adapter source changed since installation')
    candidate = binding['transcript_candidate']
    require(candidate == record['selected_transcript_candidate'] and T.digest(candidate['path']) == candidate['sha256'],
            'Native CADD transcript candidate changed')
    if config.get('cadd_base_dir'):
        require(profile['annotation']['env']['CADD_SOURCE_ROOT'] == str(Path(config['cadd_base_dir']).absolute()),
                'Configured CADD source differs')
    wrappers = binding['wrappers']; wrapper = Path(wrappers['path'])
    require(not wrapper.is_symlink() and T.digest(wrapper) == wrappers['sha256'], 'Native interpreter wrapper changed')
    for name, value in wrappers['entry_points'].items():
        path = Path(value)
        require(path.is_symlink() and os.readlink(path) == wrapper.name and path.resolve() == wrapper.resolve(),
                'Native interpreter alias changed: ' + name)
    environments = record['runtime']['environments']
    require(set(environments) == set(wrappers['entry_points']), 'Native environment/entry-point membership differs')
    for name, value in environments.items():
        require(Path(value['prefix']) == wrapper.parent.parent / 'envs' / name,
                'Native environment prefix differs from the verified interpreter route')
        require(T.digest(value['receipt']) == value['sha256'], 'Native environment receipt changed: ' + name)
        R.verify_native_environment(value['prefix'], T.read_json(value['receipt']))
    for checkpoint in binding['esm_hub']['checkpoints']:
        link, target = Path(checkpoint['link']), Path(checkpoint['target'])
        require(link.is_symlink() and link.resolve() == target.resolve(), 'Native ESM checkpoint route changed')
    calculation = record['calculation_runtime']
    require(profile['calculation']['program'] == calculation['program']
            and profile['calculation']['program_sha256'] == calculation['sha256'] == T.PYTHON_CALCULATION_SHA256
            and T.digest(calculation['program']) == T.PYTHON_CALCULATION_SHA256
            and profile['calculation']['env'] == calculation['environment'],
            'Native original-Python calculation binding changed')
    require(profile['calculation']['verify_files'] == [{'path': item['path'], 'sha256': item['sha256']}
        for item in record['calculation_assets'][:3]], 'Native calculation asset pins changed')
    env = profile['annotation']['env']; annotation = Path(record['source']) / ('data/annotations/' + assembly + '_v1.7')
    for key, name in (('CADD_PYTHON2', 'cadd-python2'), ('CADD_ESM_PYTHON', 'cadd-esm-python'),
                      ('CADD_REGSEQ_PYTHON', 'cadd-regseq-python'), ('CADD_MMSPLICE_PYTHON', 'cadd-mmsplice-python')):
        if name in wrappers['entry_points']:
            require(env.get(key) == wrappers['entry_points'][name], 'Native annotation interpreter binding changed: ' + key)
    require(env['CADD_ESM_TRANSCRIPTS'] == candidate['environment_candidates']['CADD_ESM_TRANSCRIPTS']
            and env['CADD_ESM_MODEL_DIR'] == binding['esm_hub']['path'], 'Native ESM resource binding changed')
    if assembly == 'GRCh38':
        require(env['CADD_MMSPLICE_GTF'] == candidate['environment_candidates']['CADD_MMSPLICE_GTF']
                and env['CADD_REFERENCE_FASTA'] == str(annotation / 'reference/reference.fa'),
                'Native MMSplice resource binding changed')
    if config.get('ref_genome'):
        require(T.record(config['ref_genome']) == binding['vep_fasta']
                and T.record(str(config['ref_genome']) + '.fai') == binding['vep_fasta_index'],
                'Configured single-VEP reference/index differs')
    for name in ('configuration', 'prepare', 'extract', 'snakefile'):
        pin = profile['prescoring'][name]
        require(T.digest(pin['path']) == pin['sha256'], 'Original prescore source changed: ' + name)
    current = prescoring(Path(record['source']), assembly)
    require(all(profile['prescoring'].get(key) == value for key, value in current.items() if key != 'qualification'),
            'Original publisher prescore order or identities changed')
    for stage in ('annotation', 'calculation'): T.verify_resources(profile, stage)
    return profile


def preflight(config, profile_path):
    """Read-only execution gate; neither model execution nor output creation."""
    profile = T.load_profile(profile_path, T.ALIASES.get(config.get('assembly')))
    S.verify_profile(profile)
    require(profile['annotation']['shared_vep']['adapter_id'] == S.NATIVE_ADAPTER,
            'Native CADD requires qualified official native selection')
    verify_prescoring_profile(profile)
    profile = check(config, profile_path)
    for stage in ('annotation', 'calculation'): T.stage_preflight(profile, stage)
    return profile


def run(config, profile_path, shared_vep_bundle, output):
    """Run only caller-proven true-novel rows; retain all original scientific output."""
    profile_path = Path(profile_path).absolute()
    profile = preflight(config, profile_path)
    S.verify(shared_vep_bundle, profile)
    verify_novel_input(profile, profile_path, shared_vep_bundle)
    root = Path(output).absolute()
    require(not os.path.lexists(root), 'Fresh native CADD output directory required')
    temporary = Path(config.get('tmp_dir', ''))
    require(temporary.is_absolute() and temporary.is_dir(), 'Configured tmp_dir must be an existing absolute directory')
    root.mkdir(parents=True)
    work = Path(tempfile.mkdtemp(prefix='priva-cadd-native-', dir=temporary))
    result = dict(execution_status='failed', profile_sha256=T.digest(profile_path), work_directory=str(work), **FLAGS)
    try:
        for command, source, folder in (('annotate', shared_vep_bundle, 'annotation_bundle'),
                                       ('calculate', root / 'annotation_bundle', 'calculation_bundle')):
            args = argparse.Namespace(command=command, profile=str(profile_path), assembly=profile['contract']['assembly'],
                input=str(source), output=str(root / folder), annotation_source='shared-vep', stub=False)
            T.execute(args, source_scope='native', work_dir=work / folder)
        require(T.digest(profile_path) == result['profile_sha256'], 'Native CADD profile changed during execution')
        verify_novel_input(profile, profile_path, shared_vep_bundle)
        result.update(execution_status='completed', process_exit_code=0, annotation_manifest=file_pin(root / 'annotation_bundle/manifest.json'),
                      calculation_manifest=file_pin(root / 'calculation_bundle/manifest.json'))
    except BaseException as error:
        result['error'] = repr(error)
        if isinstance(error, subprocess.CalledProcessError): result['process_exit_code'] = error.returncode
        raise
    finally:
        T.write_json(root / 'run.json', result)
    return result


def arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('install-profile', 'check', 'preflight', 'run'):
        child = commands.add_parser(name); child.add_argument('--config', type=Path, required=True)
        if name == 'install-profile':
            child.add_argument('--runtime-receipt', type=Path, required=True)
            child.add_argument('--verify-only', action='store_true')
        else: child.add_argument('--profile', type=Path, required=True)
        if name == 'run':
            child.add_argument('--shared-vep-bundle', type=Path, required=True)
            child.add_argument('--output', type=Path, required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = arguments(argv)
    try:
        import yaml
        config = yaml.safe_load(args.config.read_text())
        if args.command == 'install-profile': result = install_profile(config, args.runtime_receipt, args.verify_only)
        elif args.command == 'run': result = run(config, args.profile, args.shared_vep_bundle, args.output)
        else:
            (check if args.command == 'check' else preflight)(config, args.profile)
            result = {'status': 'components_verified' if args.command == 'check' else 'execution_preflight_passed', **FLAGS}
        print(json.dumps(result, sort_keys=True)); return 0
    except subprocess.CalledProcessError as error:
        print('Native CADD process failed: ' + str(error), file=sys.stderr)
        return error.returncode if error.returncode > 0 else 128 - error.returncode
    except (T.ContractError, S.SharedVepError, OSError, KeyError, TypeError, ValueError, RuntimeError) as error:
        print('Native CADD contract: ' + str(error), file=sys.stderr); return 2


if __name__ == '__main__': raise SystemExit(main())
