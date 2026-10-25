#!/usr/bin/env python3
"""Freeze/build/test the opt-in R7 boundary, without scientific activation.

format --repo REPO --output NEW_DATA_ROOT
prepare --repo REPO --oracle ROOT --oracle-manifest-sha256 SHA
        --work-root NEW_DATA_VALIDATION --build-root NEW_DATA_BUILD
execute --work-root DATA_VALIDATION --pins-sha256 SHA
        [--reuse-built-work-root CLOSED_ROOT --reuse-built-pins-sha256 SHA
         --binary-sha256 SHA]
All previous source bodies/tests must remain byte-identical to the closed R5
proof. Fresh Cargo tests cover all modules; only the new R7 differential runs.
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

ORACLE = '739a3c6fdcbad01314237b5cd5e0d7210dbd5e8d82c7e977ebc40041a7f27b04'
PREVIOUS = D / 'validation/rust-ranking-gene-weights-20261024-v2/validation-v2'
PREVIOUS_PINS = '6763d46bd8bf1969f6437482e280b4c21bfd942243cf5de2695ddd5d287cefa4'
NUMERIC = D / 'validation/ranking-pext-numeric-20261025-v1'
DISPATCH = D / 'validation/ranking-pext-dispatch-20261025-v1'
DISPATCH_MANIFEST = '6fa9c7b144790c045895c04fa56666667a100bed502d068a69e394fbd2a1f883'
COMPILER_SHA = 'e6abf55ab1859e7c990be77fd593f51661ea99fb2fde626ae9f4108aaa34e7fa'
ADDITIONS = {
    'rust/src/acmg_scoring/mod.rs': b'pub mod ranking_pext;\npub(crate) mod ranking_pext_transport;\n',
    'rust/src/cli.rs': b'    if args.get(1).map(String::as_str) == Some("ranking-pext") {\n        return crate::acmg_scoring::ranking_pext_transport::dispatch(&args[2..]);\n    }\n',
}


def format_sources(args):
    producing = native()
    assert args.repo.resolve().is_relative_to(Y) and args.output.resolve().is_relative_to(D)
    args.output.mkdir(parents=True, exist_ok=False)
    files = ['rust/src/acmg_scoring/ranking_pext.rs', 'rust/src/acmg_scoring/ranking_pext_transport.rs', 'rust/tests/ranking_pext.rs']
    before = {}
    for name in files:
        path = args.repo / name
        before[name] = digest(path)
        target = args.output / 'source' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())
    command = [TOOLCHAIN / 'rustfmt', '--edition', '2024'] + [args.repo / name for name in files]
    result = subprocess.run(command, capture_output=True, text=True)
    (args.output / 'rustfmt.log').write_text(result.stdout + result.stderr)
    (args.output / 'format.json').write_text(json.dumps(dict(native=producing, returncode=result.returncode,
        source_before=before, source_after={name: digest(args.repo / name) for name in files}), indent=2) + '\n')
    result.check_returncode()
    print(json.dumps(dict(native=producing, returncode=0)))


def prior_sources(repo):
    assert digest(PREVIOUS / 'source-pins.json') == PREVIOUS_PINS
    previous = json.loads((PREVIOUS / 'source-pins.json').read_text())
    result = json.loads((PREVIOUS / 'validation.json').read_text())
    comparison = json.loads((PREVIOUS / 'comparison/comparison.json').read_text())
    assert result['completed'] and result['rust_tests_passed'] == 147
    assert comparison['successful'] and (comparison['actual_outputs'], comparison['deferred_to_python'], comparison['deferred_original_failures']) == (120, 83, 28)
    for name, wanted in previous['source_files'].items():
        old, current = PREVIOUS / 'source' / name, repo / name
        assert digest(old) == wanted
        if name in ADDITIONS:
            extra = ADDITIONS[name]
            assert current.read_bytes().count(extra) == 1 and current.read_bytes().replace(extra, b'') == old.read_bytes(), name
        elif name == 'rust/source_functions.tsv':
            assert current.read_bytes().startswith(old.read_bytes()), name
        else:
            assert digest(current) == wanted, name
    return dict(work_root=str(PREVIOUS), source_pins_sha256=PREVIOUS_PINS,
        validation_sha256=digest(PREVIOUS / 'validation.json'), native=result['native'],
        all_previous_bodies_and_tests_unchanged=True, prior_Python_suites_not_repeated=True)


def oracle_verify(root, wanted):
    assert root.resolve().is_relative_to(D) and digest(root / 'manifest.json') == wanted == ORACLE
    for name, expected in json.loads((root / 'manifest.json').read_text()).items():
        path = root / name
        assert not path.is_symlink() and path.resolve().is_relative_to(root) and digest(path) == expected, name
    result = json.loads((root / 'oracle.json').read_text())
    assert len(result['cases']) == 201 and sum(case['original_exception'] is None for case in result['cases']) == 168
    assert len(result['preceding_full_source_failures']) == 13 and result['source_commit'] == COMMIT
    return result


def numeric_proof():
    result = json.loads((NUMERIC / 'result.json').read_text())
    assert result['success'] and (result['values_compared'], result['powf_mismatches'], result['explicit_half_mismatches']) == (807012, 0, 32)
    assert result['numpy_version'] == '1.26.4' and result['source_sha256']['compiler'] == COMPILER_SHA == digest(TOOLCHAIN / 'rustc')
    assert not result['full_R7_qualification'] and not result['production_backend_changed']
    assert digest(DISPATCH / 'manifest.json') == DISPATCH_MANIFEST
    for name, wanted in json.loads((DISPATCH / 'manifest.json').read_text()).items():
        assert digest(DISPATCH / name) == wanted, name
    dispatch = json.loads((DISPATCH / 'result.json').read_text())
    assert (dispatch['cases'], dispatch['values'], dispatch['series_vs_direct_ufunc_mismatches']) == (252, 807012, 32)
    assert dispatch['numexpr_min_elements'] == 1000000
    for profile in dispatch['profiles']:
        assert profile['compiler_sha256'] == COMPILER_SHA
        assert (profile['generic_powf_mismatches'], profile['source_scalar_dispatch_mismatches']) == (32, 0)
    assert {profile['opt_level'] for profile in dispatch['profiles']} == {0, 2}
    return dict(root=str(NUMERIC), result_sha256=digest(NUMERIC / 'result.json'),
        compiler_sha256=COMPILER_SHA, ndarray_only_generic_powf_exact_values=807012,
        ndarray_only_half_power_candidate_mismatches=32, probe_profile=result['compile_command'],
        actual_Series_dispatch_root=str(DISPATCH), actual_Series_dispatch_manifest_sha256=DISPATCH_MANIFEST,
        actual_Series_dispatch_result_sha256=digest(DISPATCH / 'result.json'),
        actual_Series_generic_powf_mismatches=32, source_scalar_dispatch_mismatches=0,
        source_scalar_dispatch_profiles=[0, 2], probe_is_not_full_R7_proof=True)


def verify(root, wanted):
    assert root.resolve().is_relative_to(D) and digest(root / 'source-pins.json') == wanted
    pins = json.loads((root / 'source-pins.json').read_text())
    for name, expected in pins['source_files'].items():
        assert digest(root / 'source' / name) == expected, name
    for path, expected in pins['toolchain'].items():
        assert digest(Path(path)) == expected, path
    oracle_verify(Path(pins['oracle']), pins['oracle_manifest_sha256'])
    assert prior_sources(root / 'source') == pins['prior_qualification']
    assert numeric_proof() == pins['numerical_groundwork']
    return pins


def prepare(args):
    producing = native()
    assert args.repo.resolve().is_relative_to(Y)
    oracle = oracle_verify(args.oracle, args.oracle_manifest_sha256)
    prior = prior_sources(args.repo)
    proof = numeric_proof()
    for path in (args.work_root, args.build_root):
        assert path.resolve().is_relative_to(D) and not path.exists()
        path.mkdir(parents=True)
    files = [args.repo / 'rust' / name for name in ['Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml', 'source_functions.tsv']]
    files += sorted((args.repo / 'rust/src').rglob('*.rs'))
    files += sorted((args.repo / 'rust/tests').glob('*.rs')) + sorted((args.repo / 'rust/tests').glob('*.py'))
    files += [args.repo / 'rust/src/acmg_variant_mechanism/PANDAS_LICENSE']
    files += [args.repo / 'scripts' / name for name in ['acmg_consequence.py', 'acmg_variant_mechanism.py', 'acmg_scoring.py']]
    pins = dict(created_utc=utc(), prepare_native=producing, repo_path=str(args.repo), build_root=str(args.build_root),
        source_files={}, toolchain={}, baseline_commit=COMMIT, prior_qualification=prior,
        oracle=str(args.oracle), oracle_manifest_sha256=args.oracle_manifest_sha256, oracle_native=oracle['native'],
        numerical_groundwork=proof, actual_Cargo_profile='test/dev, opt-level=0, debug=0; no fast-math',
        scope='R7 L621-665 only; existing scalar parser, actual boxed primitives, finite valid f64 parameters/native scores/aligned R4 mask; atomic current-stage Python fallback',
        production_backend_selected=False, scientific_acceptance=False)
    for path in files:
        assert not path.is_symlink()
        raw, wanted = identity(path), digest(path)
        name = str(path.relative_to(args.repo))
        target = args.work_root / 'source' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        assert identity(path) == raw and digest(path) == digest(target) == wanted
        pins['source_files'][name] = wanted
    assert pins['source_files']['scripts/acmg_scoring.py'] == SOURCE_SHA
    assert pins['source_files']['rust/tests/ranking_pext_oracle.py'] == oracle['harness_sha256']
    for name in ['cargo', 'rustc', 'rustdoc', 'rustfmt']:
        path = TOOLCHAIN / name
        pins['toolchain'][str(path)] = digest(path)
    with (args.work_root / 'source-pins.json').open('x') as stream:
        json.dump(pins, stream, indent=2); stream.write('\n')
    wanted = digest(args.work_root / 'source-pins.json')
    verify(args.work_root, wanted)
    (args.work_root / 'logs').mkdir()
    print(json.dumps(dict(pins_sha256=wanted, source_files=len(pins['source_files']), native=producing)))


def execute(args):
    pins = verify(args.work_root, args.pins_sha256)
    producing = native()
    build, source, logs = Path(pins['build_root']), args.work_root / 'source', args.work_root / 'logs'
    cargo_home, temporary = build / 'cargo-home', build / 'tmp'
    cargo_home.mkdir(); (cargo_home / 'registry').symlink_to(REGISTRY); temporary.mkdir()
    env = dict(os.environ, CARGO_HOME=str(cargo_home), CARGO_TARGET_DIR=str(build / 'target'), CARGO_BUILD_JOBS='2',
        RUSTC=str(TOOLCHAIN / 'rustc'), RUSTDOC=str(TOOLCHAIN / 'rustdoc'),
        PATH=str(TOOLCHAIN) + ':' + str(PYTHON.parent) + ':/usr/bin:/bin', POLARS_MAX_THREADS='1',
        OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', TMPDIR=str(temporary), RUSTFLAGS='-Dwarnings')
    report = dict(started_utc=utc(), native=producing, source_pins_sha256=args.pins_sha256,
        completed=False, prior_qualification=pins['prior_qualification'], numerical_groundwork=pins['numerical_groundwork'],
        scope=pins['scope'], production_backend_selected=False, scientific_acceptance=False)
    try:
        binary = args.work_root / 'priva'
        reuse = [args.reuse_built_work_root, args.reuse_built_pins_sha256, args.binary_sha256]
        if any(value is not None for value in reuse):
            assert all(value is not None for value in reuse)
            old = verify(args.reuse_built_work_root, args.reuse_built_pins_sha256)
            assert old['source_files'].keys() == pins['source_files'].keys() and old['toolchain'] == pins['toolchain']
            assert old['oracle'] == pins['oracle'] and old['oracle_manifest_sha256'] == pins['oracle_manifest_sha256']
            changed = {name for name in old['source_files'] if old['source_files'][name] != pins['source_files'][name]}
            assert changed <= {'rust/tests/ranking_pext_compare.py', 'rust/tests/ranking_pext_validate.py'}
            previous = json.loads((args.reuse_built_work_root / 'validation.json').read_text())
            assert previous['rust_tests_passed'] == 162 and previous['binary_sha256'] == args.binary_sha256 == digest(args.reuse_built_work_root / 'priva')
            cargo_log = args.reuse_built_work_root / 'logs/cargo-tests.log'
            counts = re.findall(r'test result: ok\. (\d+) passed; 0 failed;', cargo_log.read_text())
            assert sum(map(int, counts)) == 162 and not re.search(r'test result: FAILED|error:', cargo_log.read_text())
            shutil.copy2(args.reuse_built_work_root / 'priva', binary)
            report['rust_tests_passed'] = 162
            report['reused_completed_build_and_rust_tests'] = dict(work_root=str(args.reuse_built_work_root),
                source_pins_sha256=args.reuse_built_pins_sha256, actual_native=previous['native'],
                previous_whole_validation_completed=previous['completed'], previous_failure=previous.get('failure'),
                changed_harness_files=sorted(changed), unchanged_Rust_tests_toolchain_and_profile=True,
                binary_sha256=args.binary_sha256, cargo_log_sha256=digest(cargo_log))
        else:
            invoke([TOOLCHAIN / 'rustc', '--version', '--verbose'], logs / 'toolchain.log', env)
            invoke([TOOLCHAIN / 'cargo', 'fmt', '--manifest-path', source / 'rust/Cargo.toml', '--check'], logs / 'rustfmt.log', env)
            invoke([TOOLCHAIN / 'cargo', 'test', '--locked', '--offline', '--manifest-path', source / 'rust/Cargo.toml', '-j', '2'], logs / 'cargo-tests.log', env)
            counts = re.findall(r'test result: ok\. (\d+) passed; 0 failed;', (logs / 'cargo-tests.log').read_text())
            report['rust_tests_passed'] = sum(map(int, counts))
            assert report['rust_tests_passed'] == 162
            shutil.copy2(build / 'target/debug/priva', binary)
        report['binary_sha256'] = digest(binary)
        invoke([PYTHON, '-B', source / 'rust/tests/ranking_pext_compare.py', 'compare', '--binary', binary,
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
            json.dump(report, stream, indent=2); stream.write('\n')
        print(json.dumps(report), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='procedure', required=True)
    command = commands.add_parser('format')
    for name in ('repo', 'output'):
        command.add_argument('--' + name, type=Path, required=True)
    command = commands.add_parser('prepare')
    for name in ('repo', 'oracle', 'work-root', 'build-root'):
        command.add_argument('--' + name, type=Path, required=True)
    command.add_argument('--oracle-manifest-sha256', required=True)
    command = commands.add_parser('execute')
    command.add_argument('--work-root', type=Path, required=True)
    command.add_argument('--pins-sha256', required=True)
    command.add_argument('--reuse-built-work-root', type=Path)
    command.add_argument('--reuse-built-pins-sha256')
    command.add_argument('--binary-sha256')
    args = parser.parse_args()
    {'format': format_sources, 'prepare': prepare, 'execute': execute}[args.procedure](args)


if __name__ == '__main__':
    main()
