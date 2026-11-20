#!/usr/bin/env python3
"""Freeze and validate the optional coerce-pext boundary in real Slurm steps.

prepare --repo REPO --oracle ROOT --oracle-manifest-sha256 SHA
        --work-root NEW_DATA_VALIDATION --build-root NEW_DATA_BUILD
execute --work-root DATA_VALIDATION --pins-sha256 SHA
The actual Python oracle must already exist, captured before Rust implementation.
No production backend, scientific Python, installer cache or full model is changed.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

Y = Path('/stornext/snfs170/next-gen/scratch/yangyxt')
D = Path('/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade')
TOOLCHAIN = Path('/users/u266734/.rustup/toolchains/1.90.0-x86_64-unknown-linux-gnu/bin')
PYTHON = Y / 'mamba/envs/priva_acmg/bin/python'
REGISTRY = Y / 'softwares/rust/cargo/registry'
SOURCE_SHA = 'a9efd0964f25e84dc01ad276d4f0905a10492a5d1f94c2b277573ac94244d354'
COMMIT = '09223ef5a8b0ceb4c6d61055c6eefe62b51ba903'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(path):
    stat = path.stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def native():
    assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_STEP_ID'), 'real Slurm step required'
    return dict(job_id=os.environ['SLURM_JOB_ID'], step_id=os.environ['SLURM_STEP_ID'])


def oracle_verify(root, expected):
    assert root.resolve().is_relative_to(D)
    assert digest(root / 'manifest.json') == expected, 'preimplementation oracle pin changed'
    manifest = json.loads((root / 'manifest.json').read_text())
    assert manifest['schema'] == 'priva.coerce-pext-python-oracle.v1'
    assert manifest['source_sha256'] == SOURCE_SHA
    assert manifest['source_commit'] == COMMIT
    assert manifest['runtime'] == dict(python='3.12.3', pandas='2.2.3', numpy='1.26.4', unicode='15.0.0')
    for name, wanted in manifest['files'].items():
        assert Path(name).name == name and not (root / name).is_symlink()
        assert digest(root / name) == wanted, name
    return manifest


def verify(root, expected):
    assert root.resolve().is_relative_to(D)
    assert digest(root / 'source-pins.json') == expected, 'source manifest changed'
    pins = json.loads((root / 'source-pins.json').read_text())
    for name, wanted in pins['source_files'].items():
        assert digest(root / 'source' / name) == wanted, name
    for path, wanted in pins['toolchain'].items():
        assert digest(Path(path)) == wanted, path
    oracle_verify(Path(pins['oracle']), pins['oracle_manifest_sha256'])
    return pins


def invoke(command, log, env=None):
    with log.open('x') as stream:
        process = subprocess.run(list(map(str, command)), stdout=stream, stderr=subprocess.STDOUT, env=env)
    assert process.returncode == 0, f'command failed ({process.returncode}); retained {log}'


def prepare(args):
    producing = native()
    assert args.repo.resolve().is_relative_to(Y)
    oracle = oracle_verify(args.oracle, args.oracle_manifest_sha256)
    baseline = subprocess.run(['git', '-C', str(args.repo), 'show', COMMIT + ':scripts/acmg_consequence.py'],
                              check=True, capture_output=True).stdout
    assert hashlib.sha256(baseline).hexdigest() == SOURCE_SHA
    for path in [args.work_root, args.build_root]:
        assert path.resolve().is_relative_to(D) and not path.exists(), 'fresh contained paths required'
    args.work_root.mkdir(parents=True)
    args.build_root.mkdir(parents=True)
    source = args.work_root / 'source'
    files = [args.repo / 'rust' / name for name in ['Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml', 'source_functions.tsv']]
    files += sorted((args.repo / 'rust/src').rglob('*.rs'))
    files += sorted((args.repo / 'rust/tests').glob('*.rs'))
    files += sorted((args.repo / 'rust/tests').glob('*.py'))
    files += [args.repo / 'rust/src/acmg_variant_mechanism/PANDAS_LICENSE']
    files += [args.repo / 'scripts' / name for name in ['acmg_consequence.py', 'acmg_variant_mechanism.py', 'acmg_scoring.py']]
    pins = dict(created_utc=utc(), prepare_native=producing, repo_path=str(args.repo),
                build_root=str(args.build_root), oracle=str(args.oracle),
                oracle_manifest_sha256=args.oracle_manifest_sha256, source_files={}, toolchain={},
                oracle_native=oracle['native'], baseline_commit=COMMIT,
                scope='independent opt-in Rust procedure; main1524420 uses separate frozen 6:39 PM CT sources and does not build/select this path',
                production_backend_selected=False, scientific_acceptance=False)
    for path in files:
        assert not path.is_symlink(), f'source alias: {path}'
        before = identity(path)
        wanted = digest(path)
        name = str(path.relative_to(args.repo))
        target = source / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        assert digest(target) == wanted and identity(path) == before and digest(path) == wanted, name
        pins['source_files'][name] = wanted
    assert pins['source_files']['scripts/acmg_consequence.py'] == SOURCE_SHA
    assert pins['source_files']['rust/tests/coerce_pext_oracle.py'] == oracle['harness_sha256'], 'actual oracle harness changed'
    for name in ['cargo', 'rustc', 'rustdoc', 'rustfmt']:
        path = TOOLCHAIN / name
        pins['toolchain'][str(path)] = digest(path)
    with (args.work_root / 'source-pins.json').open('x') as stream:
        json.dump(pins, stream, indent=2)
        stream.write('\n')
    expected = digest(args.work_root / 'source-pins.json')
    verify(args.work_root, expected)
    (args.work_root / 'logs').mkdir()
    print(json.dumps(dict(pins_sha256=expected, source_files=len(pins['source_files']),
                          oracle_manifest_sha256=args.oracle_manifest_sha256, native=producing)))


def execute(args):
    producing = native()
    pins = verify(args.work_root, args.pins_sha256)
    build, source, logs = Path(pins['build_root']), args.work_root / 'source', args.work_root / 'logs'
    cargo_home, temporary = build / 'cargo-home', build / 'tmp'
    cargo_home.mkdir()
    (cargo_home / 'registry').symlink_to(REGISTRY)
    temporary.mkdir()
    env = dict(os.environ, CARGO_HOME=str(cargo_home), CARGO_TARGET_DIR=str(build / 'target'),
               CARGO_BUILD_JOBS='2', RUSTC=str(TOOLCHAIN / 'rustc'), RUSTDOC=str(TOOLCHAIN / 'rustdoc'),
               PATH=str(TOOLCHAIN) + ':' + str(PYTHON.parent) + ':/usr/bin:/bin',
               POLARS_MAX_THREADS='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               PYTHONDONTWRITEBYTECODE='1', TMPDIR=str(temporary), RUSTFLAGS='-Dwarnings')
    report = dict(started_utc=utc(), native=producing, source_pins_sha256=args.pins_sha256,
                  oracle_manifest_sha256=pins['oracle_manifest_sha256'], completed=False,
                  production_backend_selected=False, scientific_acceptance=False, scope=pins['scope'])
    try:
        invoke([TOOLCHAIN / 'rustc', '--version', '--verbose'], logs / 'toolchain.log', env)
        invoke([TOOLCHAIN / 'cargo', 'fmt', '--manifest-path', source / 'rust/Cargo.toml', '--check'], logs / 'rustfmt.log', env)
        invoke([TOOLCHAIN / 'cargo', 'test', '--locked', '--offline', '--manifest-path',
                source / 'rust/Cargo.toml', '-j', '2'], logs / 'cargo-tests.log', env)
        binary = args.work_root / 'priva'
        shutil.copy2(build / 'target/debug/priva', binary)
        report['binary_sha256'] = digest(binary)
        invoke([PYTHON, '-B', source / 'rust/tests/coerce_pext_compare.py', 'compare', '--binary', binary,
                '--binary-sha256', report['binary_sha256'], '--oracle', pins['oracle'],
                '--oracle-manifest-sha256', pins['oracle_manifest_sha256'], '--output', args.work_root / 'comparison'],
               logs / 'coerce-pext-comparison.log', env)
        condition = args.work_root / 'condition-oracle'
        invoke([PYTHON, '-B', source / 'rust/tests/condition_masks_oracle.py', 'capture', '--source-file',
                source / 'scripts/acmg_variant_mechanism.py', '--output', condition], logs / 'condition-oracle.log', env)
        invoke([PYTHON, '-B', source / 'rust/tests/condition_masks_compare.py', 'compare', '--binary', binary,
                '--binary-sha256', report['binary_sha256'], '--oracle', condition,
                '--oracle-manifest-sha256', digest(condition / 'manifest.json'), '--output', args.work_root / 'condition-comparison'],
               logs / 'condition-comparison.log', env)
        mechanism = args.work_root / 'mechanism-oracle'
        invoke([PYTHON, '-B', source / 'rust/tests/mechanism_masks_oracle.py', 'capture', '--output', mechanism],
               logs / 'mechanism-oracle.log', env)
        mechanism_manifest = args.work_root / 'mechanism-oracle-manifest.json'
        mechanism_manifest.write_text(json.dumps({p.name: digest(p) for p in sorted(mechanism.iterdir()) if p.is_file()}, indent=2) + '\n')
        invoke([PYTHON, '-B', source / 'rust/tests/mechanism_masks_compare.py', 'compare', '--binary', binary,
                '--oracle', mechanism, '--oracle-manifest', mechanism_manifest, '--output', args.work_root / 'mechanism-comparison'],
               logs / 'mechanism-comparison.log', env)
        invoke([PYTHON, '-B', source / 'rust/tests/differential_scoring.py', binary], logs / 'scoring-comparison.log', env)
        verify(args.work_root, args.pins_sha256)
        assert digest(binary) == report['binary_sha256']
        report['completed'] = True
    except Exception as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['finished_utc'] = utc()
        with (args.work_root / 'validation.json').open('x') as stream:
            json.dump(report, stream, indent=2)
            stream.write('\n')
        print(json.dumps(report), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest='command', required=True)
    prepare_parser = command.add_parser('prepare')
    for name in ['repo', 'oracle', 'work-root', 'build-root']:
        prepare_parser.add_argument('--' + name, type=Path, required=True)
    prepare_parser.add_argument('--oracle-manifest-sha256', required=True)
    execute_parser = command.add_parser('execute')
    execute_parser.add_argument('--work-root', type=Path, required=True)
    execute_parser.add_argument('--pins-sha256', required=True)
    args = parser.parse_args()
    (prepare if args.command == 'prepare' else execute)(args)


if __name__ == '__main__':
    main()
