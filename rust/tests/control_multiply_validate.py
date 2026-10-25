#!/usr/bin/env python3
"""Freeze/build/compare the opt-in original R3 multiplication/assignment in bounded Slurm steps.

prepare --repo REPO --oracle ROOT --oracle-manifest-sha256 SHA
        --pedigree-oracle ROOT --pedigree-oracle-manifest-sha256 SHA
        --setup-oracle ROOT --setup-oracle-manifest-sha256 SHA
        --terminal-oracle ROOT --terminal-oracle-manifest-sha256 SHA
        --integer-oracle ROOT --integer-oracle-manifest-sha256 SHA
        --work-root NEW_DATA_VALIDATION --build-root NEW_DATA_BUILD
execute --work-root DATA_VALIDATION --pins-sha256 SHA
        [--reuse-built-work-root CLOSED_ROOT --reuse-built-pins-sha256 SHA
         --binary-sha256 SHA]
The actual-Python L528 oracle predates this Rust implementation and binds normal
original C1-C5 returns to full-ranker R4 checkpoints. C1-C5/SciPy and all prior
stages remain unchanged Python. Fallback runs only L528 on current buffers, or
whole-Python dispatch must occur before any earlier stage in later integration.
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


def oracle_verify(root, wanted, count):
    assert root.resolve().is_relative_to(D) and digest(root / 'manifest.json') == wanted
    for name, expected in json.loads((root / 'manifest.json').read_text()).items():
        path = root / name
        assert path.resolve().is_relative_to(root) and not path.is_symlink()
        assert digest(path) == expected, name
    oracle = json.loads((root / 'oracle.json').read_text())
    assert oracle['source_commit'] == COMMIT and len(oracle['cases']) == count
    if count == 278:
        witness = json.loads((root / 'sex-predicate-unicode-witness.json').read_text())
        assert witness['nonascii_lower_entirely_recognized_ascii_letters'] == []
    return oracle


def verify(root, wanted):
    assert root.resolve().is_relative_to(D) and digest(root / 'source-pins.json') == wanted
    pins = json.loads((root / 'source-pins.json').read_text())
    for name, expected in pins['source_files'].items():
        assert digest(root / 'source' / name) == expected, name
    for name, expected in pins['toolchain'].items():
        assert digest(Path(name)) == expected, name
    for key, count in [('oracle', 116), ('pedigree_oracle', 278), ('setup_oracle', 129), ('terminal_oracle', 91), ('integer_oracle', 15)]:
        oracle_verify(Path(pins[key]), pins[key + '_manifest_sha256'], count)
    return pins


def prepare(args):
    producing = native()
    assert args.repo.resolve().is_relative_to(Y)
    oracles = {}
    for key, count in [('oracle', 116), ('pedigree_oracle', 278), ('setup_oracle', 129), ('terminal_oracle', 91), ('integer_oracle', 15)]:
        oracles[key] = oracle_verify(getattr(args, key), getattr(args, key + '_manifest_sha256'), count)
    baseline = subprocess.run(['git', '-C', str(args.repo), 'show', COMMIT + ':scripts/acmg_scoring.py'],
                              check=True, capture_output=True).stdout
    assert hashlib.sha256(baseline).hexdigest() == SOURCE_SHA
    for path in [args.work_root, args.build_root]:
        assert path.resolve().is_relative_to(D) and not path.exists(), 'fresh contained paths required'
    args.work_root.mkdir(parents=True)
    args.build_root.mkdir(parents=True)
    files = [args.repo / 'rust' / name for name in ['Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml', 'source_functions.tsv']]
    files += sorted((args.repo / 'rust/src').rglob('*.rs'))
    files += sorted((args.repo / 'rust/tests').glob('*.rs'))
    files += sorted((args.repo / 'rust/tests').glob('*.py'))
    files += [args.repo / 'rust/src/acmg_variant_mechanism/PANDAS_LICENSE']
    files += [args.repo / 'scripts' / name for name in ['acmg_consequence.py', 'acmg_variant_mechanism.py', 'acmg_scoring.py']]
    pins = dict(created_utc=utc(), prepare_native=producing, repo_path=str(args.repo),
                build_root=str(args.build_root), source_files={}, toolchain={}, baseline_commit=COMMIT,
                scope='R3 L528 only after unchanged Python C1-C5; in-memory composition with current-buffer stage fallback; separate from main1524420 frozen sources',
                production_backend_selected=False, scientific_acceptance=False)
    for key in oracles:
        pins[key] = str(getattr(args, key))
        pins[key + '_manifest_sha256'] = getattr(args, key + '_manifest_sha256')
    pins['oracle_native'] = oracles['oracle']['native_environment']
    for path in files:
        assert not path.is_symlink(), path
        before, wanted = identity(path), digest(path)
        name = str(path.relative_to(args.repo))
        target = args.work_root / 'source' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        assert digest(target) == wanted and identity(path) == before and digest(path) == wanted, name
        pins['source_files'][name] = wanted
    assert pins['source_files']['scripts/acmg_scoring.py'] == SOURCE_SHA
    for source, key in [('control_multiply_oracle.py', 'oracle'), ('pedigree_context_oracle.py', 'pedigree_oracle'), ('ranking_setup_oracle.py', 'setup_oracle'), ('terminal_order_oracle.py', 'terminal_oracle'),
                        ('terminal_order_integer_oracle.py', 'integer_oracle')]:
        assert pins['source_files']['rust/tests/' + source] == oracles[key]['harness_sha256'], source
    for name in ['cargo', 'rustc', 'rustdoc', 'rustfmt']:
        path = TOOLCHAIN / name
        pins['toolchain'][str(path)] = digest(path)
    with (args.work_root / 'source-pins.json').open('x') as stream:
        json.dump(pins, stream, indent=2); stream.write('\n')
    wanted = digest(args.work_root / 'source-pins.json')
    verify(args.work_root, wanted)
    (args.work_root / 'logs').mkdir()
    print(json.dumps(dict(pins_sha256=wanted, source_files=len(pins['source_files']), native=producing)), flush=True)


def execute(args):
    producing = native()
    pins = verify(args.work_root, args.pins_sha256)
    build, source, logs = Path(pins['build_root']), args.work_root / 'source', args.work_root / 'logs'
    cargo_home, temporary = build / 'cargo-home', build / 'tmp'
    cargo_home.mkdir(); (cargo_home / 'registry').symlink_to(REGISTRY); temporary.mkdir()
    env = dict(os.environ, CARGO_HOME=str(cargo_home), CARGO_TARGET_DIR=str(build / 'target'),
               CARGO_BUILD_JOBS='2', RUSTC=str(TOOLCHAIN / 'rustc'), RUSTDOC=str(TOOLCHAIN / 'rustdoc'),
               PATH=str(TOOLCHAIN) + ':' + str(PYTHON.parent) + ':/usr/bin:/bin',
               POLARS_MAX_THREADS='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               PYTHONDONTWRITEBYTECODE='1', TMPDIR=str(temporary), RUSTFLAGS='-Dwarnings')
    report = dict(started_utc=utc(), native=producing, source_pins_sha256=args.pins_sha256,
                  oracle_manifest_sha256=pins['oracle_manifest_sha256'], completed=False,
                  scope=pins['scope'], production_backend_selected=False, scientific_acceptance=False)
    try:
        binary = args.work_root / 'priva'
        reuse = [args.reuse_built_work_root, args.reuse_built_pins_sha256, args.binary_sha256]
        if any(value is not None for value in reuse):
            assert all(value is not None for value in reuse)
            previous = verify(args.reuse_built_work_root, args.reuse_built_pins_sha256)
            assert previous['source_files'].keys() == pins['source_files'].keys()
            for key in ['oracle', 'pedigree_oracle', 'setup_oracle', 'terminal_oracle', 'integer_oracle']:
                assert previous[key] == pins[key] and previous[key + '_manifest_sha256'] == pins[key + '_manifest_sha256']
            differences = {name for name in previous['source_files'] if previous['source_files'][name] != pins['source_files'][name]}
            assert differences <= {'rust/tests/control_multiply_compare.py', 'rust/tests/control_multiply_validate.py'}
            old_report = json.loads((args.reuse_built_work_root / 'validation.json').read_text())
            assert old_report['binary_sha256'] == args.binary_sha256
            assert digest(args.reuse_built_work_root / 'priva') == args.binary_sha256
            cargo_log = args.reuse_built_work_root / 'logs/cargo-tests.log'
            cargo_results = re.findall(r'test result: ok\. (\d+) passed; (\d+) failed;', cargo_log.read_text())
            assert cargo_results and sum(int(passed) for passed, _ in cargo_results) == 127
            assert all(int(failed) == 0 for _, failed in cargo_results)
            assert not re.search(r'test result: FAILED|error:', cargo_log.read_text())
            shutil.copy2(args.reuse_built_work_root / 'priva', binary)
            report['reused_completed_build_and_rust_tests'] = dict(work_root=str(args.reuse_built_work_root),
                source_pins_sha256=args.reuse_built_pins_sha256, rust_tests_passed=127,
                validation_sha256=digest(args.reuse_built_work_root / 'validation.json'), actual_native=old_report['native'],
                previous_whole_validation_completed=old_report['completed'],
                previous_failure=old_report.get('failure'),
                changed_test_harness_files=sorted(differences), unchanged_rust_scientific_and_prior_test_sources=True,
                log_sha256={path.name: digest(path) for path in sorted((args.reuse_built_work_root / 'logs').glob('*.log'))})
        else:
            invoke([TOOLCHAIN / 'rustc', '--version', '--verbose'], logs / 'toolchain.log', env)
            invoke([TOOLCHAIN / 'cargo', 'fmt', '--manifest-path', source / 'rust/Cargo.toml', '--check'], logs / 'rustfmt.log', env)
            invoke([TOOLCHAIN / 'cargo', 'test', '--locked', '--offline', '--manifest-path',
                    source / 'rust/Cargo.toml', '-j', '2'], logs / 'cargo-tests.log', env)
            shutil.copy2(build / 'target/debug/priva', binary)
        report['binary_sha256'] = digest(binary)
        for script, key, directory in [('control_multiply_compare.py', 'oracle', 'comparison'),
                ('pedigree_context_compare.py', 'pedigree_oracle', 'pedigree-comparison'),
                ('ranking_setup_compare.py', 'setup_oracle', 'setup-comparison'),
                ('terminal_order_compare.py', 'terminal_oracle', 'terminal-comparison'),
                ('terminal_order_compare.py', 'integer_oracle', 'integer-comparison')]:
            invoke([PYTHON, '-B', source / 'rust/tests' / script, 'compare', '--binary', binary,
                    '--binary-sha256', report['binary_sha256'], '--oracle', pins[key],
                    '--oracle-manifest-sha256', pins[key + '_manifest_sha256'], '--output', args.work_root / directory],
                   logs / (directory + '.log'), env)
        for prefix, script, scientific in [('coerce', 'coerce_pext', 'acmg_consequence.py'),
                                           ('condition', 'condition_masks', 'acmg_variant_mechanism.py')]:
            oracle = args.work_root / (prefix + '-oracle')
            invoke([PYTHON, '-B', source / 'rust/tests' / (script + '_oracle.py'), 'capture', '--source-file',
                    source / 'scripts' / scientific, '--output', oracle], logs / (prefix + '-oracle.log'), env)
            invoke([PYTHON, '-B', source / 'rust/tests' / (script + '_compare.py'), 'compare', '--binary', binary,
                    '--binary-sha256', report['binary_sha256'], '--oracle', oracle,
                    '--oracle-manifest-sha256', digest(oracle / 'manifest.json'), '--output', args.work_root / (prefix + '-comparison')],
                   logs / (prefix + '-comparison.log'), env)
        mechanism = args.work_root / 'mechanism-oracle'
        invoke([PYTHON, '-B', source / 'rust/tests/mechanism_masks_oracle.py', 'capture', '--output', mechanism],
               logs / 'mechanism-oracle.log', env)
        manifest = args.work_root / 'mechanism-oracle-manifest.json'
        manifest.write_text(json.dumps({p.name: digest(p) for p in sorted(mechanism.iterdir()) if p.is_file()}, indent=2) + '\n')
        invoke([PYTHON, '-B', source / 'rust/tests/mechanism_masks_compare.py', 'compare', '--binary', binary,
                '--oracle', mechanism, '--oracle-manifest', manifest, '--output', args.work_root / 'mechanism-comparison'],
               logs / 'mechanism-comparison.log', env)
        invoke([PYTHON, '-B', source / 'rust/tests/differential_scoring.py', binary], logs / 'scoring-comparison.log', env)
        invoke([PYTHON, '-B', source / 'rust/tests/differential_final_order.py', binary], logs / 'final-order-comparison.log', env)
        verify(args.work_root, args.pins_sha256)
        assert digest(binary) == report['binary_sha256']
        report['completed'] = True
    except Exception as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error)); raise
    finally:
        report['finished_utc'] = utc()
        with (args.work_root / 'validation.json').open('x') as stream:
            json.dump(report, stream, indent=2); stream.write('\n')
        print(json.dumps(report), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    prepare_parser = commands.add_parser('prepare')
    for key in ['repo', 'oracle', 'pedigree-oracle', 'setup-oracle', 'terminal-oracle', 'integer-oracle', 'work-root', 'build-root']:
        prepare_parser.add_argument('--' + key, type=Path, required=True)
    for key in ['oracle', 'pedigree-oracle', 'setup-oracle', 'terminal-oracle', 'integer-oracle']:
        prepare_parser.add_argument('--' + key + '-manifest-sha256', required=True)
    execute_parser = commands.add_parser('execute')
    execute_parser.add_argument('--work-root', type=Path, required=True)
    execute_parser.add_argument('--pins-sha256', required=True)
    execute_parser.add_argument('--reuse-built-work-root', type=Path)
    execute_parser.add_argument('--reuse-built-pins-sha256')
    execute_parser.add_argument('--binary-sha256')
    args = parser.parse_args()
    (prepare if args.command == 'prepare' else execute)(args)


if __name__ == '__main__':
    main()
