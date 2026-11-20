#!/usr/bin/env python3
"""Freeze/build and exactly validate only the new opt-in R5->R6 block.

prepare --repo REPO --oracle ROOT --oracle-manifest-sha256 SHA
        --work-root NEW_DATA_VALIDATION --build-root NEW_DATA_BUILD
execute --work-root DATA_VALIDATION --pins-sha256 SHA
        [--reuse-built-work-root CLOSED_ROOT --reuse-built-pins-sha256 SHA
         --binary-sha256 SHA]
Fresh Cargo tests include existing modules. Prior Python differentials are bound
to byte-identical old bodies/tests and their closed R4 result, not rerun here.
format --repo REPO --output NEW_DATA_ROOT freezes only the three owned Rust
files before formatting and records the actual native result/source hashes.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

from terminal_order_validate import D, Y, PYTHON, TOOLCHAIN, REGISTRY, COMMIT, SOURCE_SHA
from terminal_order_validate import digest, identity, utc, native, invoke

ORACLE = 'a585166765c2d6fe5859ffe380f2040f30d253bf8525be9be179a49f32312947'
PREVIOUS = D / 'validation/rust-ranking-compatibility-20261024-v4/validation-v2'
PREVIOUS_PINS = '684246afd4b69eb78d5457358e1d5f755e0afd9ce14919db08ae58d543e07b3e'
ADDITIONS = {
    'rust/src/acmg_scoring/mod.rs': b'pub mod ranking_gene_weights;\npub(crate) mod ranking_gene_weights_transport;\n',
    'rust/src/cli.rs': b'    if args.get(1).map(String::as_str) == Some("ranking-gene-weights") {\n        return crate::acmg_scoring::ranking_gene_weights_transport::dispatch(&args[2..]);\n    }\n',
}


def format_sources(args):
    producing = native()
    assert args.repo.resolve().is_relative_to(Y) and args.output.resolve().is_relative_to(D)
    args.output.mkdir(parents=True, exist_ok=False)
    names = ['rust/src/acmg_scoring/ranking_gene_weights.rs',
             'rust/src/acmg_scoring/ranking_gene_weights_transport.rs', 'rust/tests/ranking_gene_weights.rs']
    before = {}
    for name in names:
        path = args.repo / name
        before[name] = digest(path)
        target = args.output / 'source' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())
    process = subprocess.run([TOOLCHAIN / 'rustfmt', '--edition', '2024'] + [args.repo / name for name in names],
                             capture_output=True, text=True)
    (args.output / 'rustfmt.log').write_text(process.stdout + process.stderr)
    result = dict(native=producing, source_before=before, source_after={name: digest(args.repo / name) for name in names},
                  returncode=process.returncode, source_at_run_pre_hash_recorded=True)
    (args.output / 'format.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)
    process.check_returncode()


def prior_sources(repo):
    assert digest(PREVIOUS / 'source-pins.json') == PREVIOUS_PINS
    previous = json.loads((PREVIOUS / 'source-pins.json').read_text())
    report = json.loads((PREVIOUS / 'validation.json').read_text())
    assert report['completed'] and report['source_pins_sha256'] == PREVIOUS_PINS
    result = json.loads((PREVIOUS / 'comparison/comparison.json').read_text())
    assert result['successful'] and (result['actual_outputs'], result['original_failures'], result['deferred_to_python']) == (103, 9, 35)
    for name, wanted in previous['source_files'].items():
        old = PREVIOUS / 'source' / name
        assert digest(old) == wanted, name
        current = repo / name
        if name in ADDITIONS:
            addition = ADDITIONS[name]
            assert current.read_bytes().count(addition) == 1 and current.read_bytes().replace(addition, b'') == old.read_bytes(), name
        elif name == 'rust/source_functions.tsv':
            assert current.read_bytes().startswith(old.read_bytes())
        else:
            assert digest(current) == wanted, name
    return dict(source_pins_sha256=PREVIOUS_PINS, work_root=str(PREVIOUS),
                validation_sha256=digest(PREVIOUS / 'validation.json'), native=report['native'],
                all_prior_bodies_and_tests_unchanged=True, source_only_additions=list(ADDITIONS) + ['rust/source_functions.tsv'],
                earlier_python_differentials_not_repeated=True)


def oracle_verify(root, wanted):
    assert root.resolve().is_relative_to(D) and digest(root / 'manifest.json') == wanted == ORACLE
    for name, expected in json.loads((root / 'manifest.json').read_text()).items():
        path = root / name
        assert path.resolve().is_relative_to(root) and not path.is_symlink() and digest(path) == expected, name
    report = json.loads((root / 'oracle.json').read_text())
    assert report['source_commit'] == COMMIT and len(report['cases']) == 203
    assert sum(case['original_exception'] is None for case in report['cases']) == 175
    assert len(report['preceding_full_source_failures']) == 11
    return report


def verify(root, wanted):
    assert root.resolve().is_relative_to(D) and digest(root / 'source-pins.json') == wanted
    pins = json.loads((root / 'source-pins.json').read_text())
    for name, expected in pins['source_files'].items():
        assert digest(root / 'source' / name) == expected, name
    for name, expected in pins['toolchain'].items():
        assert digest(Path(name)) == expected, name
    oracle_verify(Path(pins['oracle']), pins['oracle_manifest_sha256'])
    assert prior_sources(root / 'source') == pins['prior_qualification']
    return pins


def prepare(args):
    producing = native()
    assert args.repo.resolve().is_relative_to(Y)
    oracle = oracle_verify(args.oracle, args.oracle_manifest_sha256)
    prior = prior_sources(args.repo)
    baseline = subprocess.run(['git', '-C', str(args.repo), 'show', COMMIT + ':scripts/acmg_scoring.py'], check=True, capture_output=True).stdout
    assert hashlib.sha256(baseline).hexdigest() == SOURCE_SHA
    for path in [args.work_root, args.build_root]:
        assert path.resolve().is_relative_to(D) and not path.exists()
        path.mkdir(parents=True)
    files = [args.repo / 'rust' / name for name in ['Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml', 'source_functions.tsv']]
    files += sorted((args.repo / 'rust/src').rglob('*.rs'))
    files += sorted((args.repo / 'rust/tests').glob('*.rs'))
    files += sorted((args.repo / 'rust/tests').glob('*.py'))
    files += [args.repo / 'rust/src/acmg_variant_mechanism/PANDAS_LICENSE']
    files += [args.repo / 'scripts' / name for name in ['acmg_consequence.py', 'acmg_variant_mechanism.py', 'acmg_scoring.py']]
    pins = dict(created_utc=utc(), prepare_native=producing, repo_path=str(args.repo), build_root=str(args.build_root),
        source_files={}, toolchain={}, baseline_commit=COMMIT, prior_qualification=prior,
        oracle=str(args.oracle), oracle_manifest_sha256=args.oracle_manifest_sha256, oracle_native=oracle['native_environment'],
        scope='R5 then R6 L599-616 only; literal parser subset, native float64/current R4 mask and exact ordered assignments; stage-only atomic Python fallback; no production/main/family claim',
        production_backend_selected=False, scientific_acceptance=False)
    for path in files:
        assert not path.is_symlink(), path
        before, wanted = identity(path), digest(path)
        name = str(path.relative_to(args.repo))
        destination = args.work_root / 'source' / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
        assert digest(destination) == wanted and identity(path) == before and digest(path) == wanted, name
        pins['source_files'][name] = wanted
    assert pins['source_files']['scripts/acmg_scoring.py'] == SOURCE_SHA
    assert pins['source_files']['rust/tests/ranking_gene_weights_oracle.py'] == oracle['harness_sha256']
    for name in ['cargo', 'rustc', 'rustdoc', 'rustfmt']:
        path = TOOLCHAIN / name
        pins['toolchain'][str(path)] = digest(path)
    with (args.work_root / 'source-pins.json').open('x') as stream:
        json.dump(pins, stream, indent=2)
        stream.write('\n')
    wanted = digest(args.work_root / 'source-pins.json')
    verify(args.work_root, wanted)
    (args.work_root / 'logs').mkdir()
    print(json.dumps(dict(pins_sha256=wanted, source_files=len(pins['source_files']), native=producing)), flush=True)


def execute(args):
    producing = native()
    pins = verify(args.work_root, args.pins_sha256)
    build, source, logs = Path(pins['build_root']), args.work_root / 'source', args.work_root / 'logs'
    cargo_home, temporary = build / 'cargo-home', build / 'tmp'
    cargo_home.mkdir()
    (cargo_home / 'registry').symlink_to(REGISTRY)
    temporary.mkdir()
    env = dict(os.environ, CARGO_HOME=str(cargo_home), CARGO_TARGET_DIR=str(build / 'target'), CARGO_BUILD_JOBS='2',
        RUSTC=str(TOOLCHAIN / 'rustc'), RUSTDOC=str(TOOLCHAIN / 'rustdoc'),
        PATH=str(TOOLCHAIN) + ':' + str(PYTHON.parent) + ':/usr/bin:/bin', POLARS_MAX_THREADS='1',
        OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', TMPDIR=str(temporary), RUSTFLAGS='-Dwarnings')
    report = dict(started_utc=utc(), native=producing, source_pins_sha256=args.pins_sha256,
        oracle_manifest_sha256=pins['oracle_manifest_sha256'], completed=False, prior_qualification=pins['prior_qualification'],
        scope=pins['scope'], production_backend_selected=False, scientific_acceptance=False)
    try:
        binary = args.work_root / 'priva'
        reuse = [args.reuse_built_work_root, args.reuse_built_pins_sha256, args.binary_sha256]
        if any(value is not None for value in reuse):
            assert all(value is not None for value in reuse)
            old = verify(args.reuse_built_work_root, args.reuse_built_pins_sha256)
            assert old['source_files'].keys() == pins['source_files'].keys()
            assert old['toolchain'] == pins['toolchain']
            assert old['oracle'] == pins['oracle'] and old['oracle_manifest_sha256'] == pins['oracle_manifest_sha256']
            differences = {name for name in old['source_files'] if old['source_files'][name] != pins['source_files'][name]}
            assert differences <= {'rust/tests/ranking_gene_weights_compare.py', 'rust/tests/ranking_gene_weights_validate.py'}
            prior = json.loads((args.reuse_built_work_root / 'validation.json').read_text())
            assert prior['binary_sha256'] == args.binary_sha256 == digest(args.reuse_built_work_root / 'priva')
            cargo_log = args.reuse_built_work_root / 'logs/cargo-tests.log'
            counts = re.findall(r'test result: ok\. (\d+) passed; 0 failed;', cargo_log.read_text())
            assert sum(map(int, counts)) == prior['rust_tests_passed'] == 147
            assert not re.search(r'test result: FAILED|error:', cargo_log.read_text())
            shutil.copy2(args.reuse_built_work_root / 'priva', binary)
            report['rust_tests_passed'] = 147
            report['reused_completed_build_and_rust_tests'] = dict(work_root=str(args.reuse_built_work_root),
                source_pins_sha256=args.reuse_built_pins_sha256, actual_native=prior['native'],
                previous_whole_validation_completed=prior['completed'], previous_failure=prior.get('failure'),
                changed_test_harness_files=sorted(differences), unchanged_Rust_and_prior_test_sources=True,
                binary_sha256=args.binary_sha256, cargo_log_sha256=digest(cargo_log))
        else:
            invoke([TOOLCHAIN / 'rustc', '--version', '--verbose'], logs / 'toolchain.log', env)
            invoke([TOOLCHAIN / 'cargo', 'fmt', '--manifest-path', source / 'rust/Cargo.toml', '--check'], logs / 'rustfmt.log', env)
            invoke([TOOLCHAIN / 'cargo', 'test', '--locked', '--offline', '--manifest-path', source / 'rust/Cargo.toml', '-j', '2'], logs / 'cargo-tests.log', env)
            rust_counts = re.findall(r'test result: ok\. (\d+) passed; 0 failed;', (logs / 'cargo-tests.log').read_text())
            report['rust_tests_passed'] = sum(map(int, rust_counts))
            assert report['rust_tests_passed'] == 147
            shutil.copy2(build / 'target/debug/priva', binary)
        report['binary_sha256'] = digest(binary)
        invoke([PYTHON, '-B', source / 'rust/tests/ranking_gene_weights_compare.py', 'compare', '--binary', binary,
            '--binary-sha256', report['binary_sha256'], '--oracle', pins['oracle'],
            '--oracle-manifest-sha256', pins['oracle_manifest_sha256'], '--output', args.work_root / 'comparison'], logs / 'comparison.log', env)
        result = json.loads((args.work_root / 'comparison/comparison.json').read_text())
        assert result['successful']
        report['differential_counts'] = {key: result[key] for key in ['actual_outputs', 'original_failures', 'deferred_to_python', 'deferred_original_failures']}
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
    commands = parser.add_subparsers(dest='command', required=True)
    command = commands.add_parser('format')
    for name in ['repo', 'output']:
        command.add_argument('--' + name, type=Path, required=True)
    command = commands.add_parser('prepare')
    for name in ['repo', 'oracle', 'work-root', 'build-root']:
        command.add_argument('--' + name, type=Path, required=True)
    command.add_argument('--oracle-manifest-sha256', required=True)
    command = commands.add_parser('execute')
    command.add_argument('--work-root', type=Path, required=True)
    command.add_argument('--pins-sha256', required=True)
    command.add_argument('--reuse-built-work-root', type=Path)
    command.add_argument('--reuse-built-pins-sha256')
    command.add_argument('--binary-sha256')
    args = parser.parse_args()
    {'format': format_sources, 'prepare': prepare, 'execute': execute}[args.command](args)


if __name__ == '__main__':
    main()
