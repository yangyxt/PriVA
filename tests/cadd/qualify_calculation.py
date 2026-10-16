#!/usr/bin/env python3
"""Real-model fixture qualification with complete differences and runtime identity.

This does not qualify annotation, public-VCF processing or family prioritization.
Run inside the allocation, with outputs under the durable validation data root.
"""
import argparse
import contextlib
import gzip
import hashlib
import importlib.util
import io
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('cadd_task', ROOT / 'bin/cadd_task.py')
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)


def invoke(runtime, script, args, stdin=None, stdout=None):
    env = dict(os.environ)
    env.pop('PYTHONPATH', None)
    env.pop('PYTHONHOME', None)
    env.update(PYTHONDONTWRITEBYTECODE='1', PYTHONNOUSERSITE='1',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    subprocess.run([runtime, str(script)] + args, stdin=stdin, stdout=stdout,
                   env=env, check=True)


def file_identity(path):
    result = {'path': str(path), 'sha256': C.digest(path), 'size': Path(path).stat().st_size}
    if str(path).endswith('.gz'):
        digest = hashlib.sha256()
        with gzip.open(path, 'rb') as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(block)
        result['uncompressed_sha256'] = digest.hexdigest()
    return result


def compare_tables(actual, expected, columns, report, exact_bytes=False):
    """Report every row/field difference, including missing/extra/order changes."""
    differences = 0
    with open(report, 'x') as output:
        try:
            for row_index, (left, right) in enumerate(itertools.zip_longest(
                    C.rows(actual, columns), C.rows(expected, columns)), 1):
                if left is None or right is None:
                    output.write(json.dumps({'row': row_index, 'kind': 'missing_actual_row' if left is None else 'extra_actual_row',
                                             'actual': left, 'expected': right}) + '\n')
                    differences += 1
                else:
                    for column, a, b in zip(columns, left, right):
                        if a != b:
                            output.write(json.dumps({'row': row_index, 'kind': 'field_or_order_difference',
                                                     'column': column, 'actual': a, 'expected': b,
                                                     'actual_variant': left[:5], 'expected_variant': right[:5]}) + '\n')
                            differences += 1
            if exact_bytes and differences == 0 and C.digest(actual) != C.digest(expected):
                output.write(json.dumps({'kind': 'serialization_difference', 'actual_sha256': C.digest(actual), 'expected_sha256': C.digest(expected)}) + '\n')
                differences += 1
        except (ValueError, OSError) as error:
            output.write(json.dumps({'kind': 'table_parse_failure', 'error': str(error)}) + '\n')
            differences += 1
    return differences


def direct_calculation(runtime, source, env, annotations, out):
    scripts = source / 'src/scripts'
    features = out / 'direct.features.csv.gz'
    plain = out / 'direct.anno.tsv'
    with gzip.open(annotations, 'rb') as handle:
        plain.write_bytes(handle.read())
    with open(plain, 'rb') as handle:
        invoke(runtime, scripts / 'trackTransformation.py', ['-b', '-c', env['CADD_IMPUTE_CONFIG'],
               '-o', str(features), '--noheader'], stdin=handle)
    raw = out / 'direct.raw.tsv'
    with open(raw, 'wb') as handle:
        invoke(runtime, scripts / 'predictSKmodel.py', ['-i', str(features), '-m', env['CADD_MODEL'],
               '-a', str(annotations)], stdout=handle)
    hierarchy = out / 'direct.hierarchy.tsv'
    with open(raw, 'rb') as src, open(hierarchy, 'wb') as dest:
        invoke(runtime, scripts / 'max_line_hierarchy.py', ['--all'], stdin=src, stdout=dest)
    direct = out / 'direct.scores.tsv'
    with open(hierarchy, 'rb') as src, open(direct, 'wb') as dest:
        invoke(runtime, scripts / 'appendPHREDscore.py', ['-t', env['CADD_PHRED_TABLE']], stdin=src, stdout=dest)
    return {'features': features, 'raw': raw, 'hierarchy': hierarchy, 'scores': direct}


def qualify_assembly(args, assembly, source, runtime, out, runtime_identity):
    profile = json.loads((ROOT / 'conf/cadd' / (assembly + '.reference.json')).read_text())
    profile.pop('annotation')
    profile['software_id'] = args.software_id + ';sif-sha256=' + runtime_identity['container']['sha256']
    profile['resource_manifest_path'] = str(out / 'resources.json')
    env = profile['calculation']['env']
    env.update(CADD_SOURCE_ROOT=str(source), CADD_PYTHON2=runtime,
               CADD_IMPUTE_CONFIG=str(source / 'config' / ('impute_' + assembly + '_v1.7.cfg')),
               CADD_MODEL=str(source / 'data/models' / assembly / ('CADDv1.7-' + assembly + '.mod')),
               CADD_PHRED_TABLE=str(source / 'data/models' / assembly / ('conversionTable_CADDv1.7-' + assembly + '.txt')))
    for item, key in zip(profile['calculation']['verify_files'], ['CADD_IMPUTE_CONFIG', 'CADD_MODEL', 'CADD_PHRED_TABLE']):
        item['path'] = env[key]
    profile['calculation']['resource_roots'] = [item['path'] for item in profile['calculation']['verify_files']]
    if args.standalone_pin:
        profile['calculation']['standalone'] = C.read_json(args.standalone_pin)
    ppath = out / 'profile.json'
    ppath.write_text(json.dumps(profile, indent=2) + '\n')
    with contextlib.redirect_stdout(io.StringIO()):
        C.lock_resources(argparse.Namespace(profile=str(ppath), output=profile['resource_manifest_path'], checksum=True))
    profile['contract']['resource_manifest_sha256'] = C.digest(profile['resource_manifest_path'])
    ppath.write_text(json.dumps(profile, indent=2) + '\n')
    fixture = source / 'test' / ('output_inclAnno_' + assembly + '.tsv.gz')
    columns = profile['contract']['columns']
    annotations = out / 'fixture.anno.tsv.gz'
    count = 0
    with gzip.open(annotations, 'wt') as handle:
        handle.write('#' + '\t'.join(columns) + '\n')
        for fields in C.rows(fixture, columns + ['RawScore', 'PHRED']):
            handle.write('\t'.join(fields[:-2]) + '\n')
            count += 1
    bundle = out / 'annotation_bundle'
    C.execute(argparse.Namespace(command='import-annotations', profile=str(ppath), assembly=assembly,
                                 input=str(annotations), output=str(bundle), stub=False))
    scored = out / 'calculation_bundle'
    C.execute(argparse.Namespace(command='calculate', profile=str(ppath), assembly=assembly,
                                 input=str(bundle), output=str(scored), stub=False))
    direct = direct_calculation(runtime, source, env, bundle / 'annotations.tsv.gz', out)
    actual = scored / 'scores.tsv'
    direct_differences = compare_tables(actual, direct['scores'], columns + ['RawScore', 'PHRED'], out / 'direct.mismatches.jsonl', exact_bytes=True)
    published_differences = compare_tables(actual, fixture, columns + ['RawScore', 'PHRED'], out / 'published-fixture.mismatches.jsonl')
    intermediate = {}
    names = {'features': 'features.csv.gz', 'raw': 'raw.tsv', 'hierarchy': 'hierarchy.tsv'}
    for stage, name in names.items():
        matches = list((out / 'calculation_bundle.work').rglob(name))
        if len(matches) != 1:
            raise AssertionError('Expected one retained original ' + stage + ' artifact; got ' + str(matches))
        left, right = file_identity(matches[0]), file_identity(direct[stage])
        key = 'uncompressed_sha256' if stage == 'features' else 'sha256'
        intermediate[stage] = {'actual': left, 'direct': right, 'exact': left[key] == right[key]}
        if stage != 'features':
            intermediate[stage]['field_differences'] = compare_tables(matches[0], direct[stage], columns + ['RawScore'], out / (stage + '.mismatches.jsonl'), exact_bytes=True)
        elif left[key] != right[key]:
            # Feature matrices have no header; compare every numeric field by its
            # ordered feature index, retaining rows and complete differences.
            with gzip.open(matches[0], 'rt') as left_file, gzip.open(direct[stage], 'rt') as right_file, open(out / 'features.mismatches.jsonl', 'x') as report:
                for number, (a, b) in enumerate(itertools.zip_longest(left_file, right_file), 1):
                    if a != b:
                        report.write(json.dumps({'row': number, 'actual_csv': a, 'expected_csv': b}) + '\n')
    exact = actual.read_bytes() == direct['scores'].read_bytes() and all(v['exact'] for v in intermediate.values())
    result = {'assembly': assembly, 'rows': count, 'direct_original_exact': exact, 'direct_field_differences': direct_differences,
              'published_fixture_identical': published_differences == 0, 'published_fixture_field_differences': published_differences,
              'published_fixture_interpretation': 'separate reference-profile provenance comparison; not hidden by direct-wrapper equality',
              'fixture': file_identity(fixture), 'annotation_input': file_identity(bundle / 'annotations.tsv.gz'),
              'output': file_identity(actual), 'direct_output': file_identity(direct['scores']), 'intermediates': intermediate,
              'profile_sha256': C.digest(ppath), 'software_id': profile['software_id'], 'upstream_commit': C.PIN,
              'runtime': runtime_identity, 'model_assets': [file_identity(item['path']) for item in profile['calculation']['verify_files']]}
    C.write_json(out / 'result.json', result)
    return result


def qualify(args):
    source = Path(args.source).resolve(strict=True)
    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    # Preserve argv[0]: installed container entry points select their environment
    # from the symlink name. Resolving that symlink would lose the Python2 route.
    runtime = str(Path(args.python2).absolute())
    runtime_identity = {}
    try:
        runtime_identity = {'container': file_identity(Path(args.runtime_image).resolve(strict=True)),
                            'entrypoint': file_identity(runtime)}
        version = subprocess.run([runtime, '-c', 'import json,sys,numpy,pandas,sklearn; print(json.dumps({"python":sys.version,"numpy":numpy.__version__,"pandas":pandas.__version__,"sklearn":sklearn.__version__}))'],
                                 capture_output=True, text=True, check=True)
        runtime_identity['versions'] = json.loads(version.stdout.strip().splitlines()[-1])
    except Exception as error:
        C.write_json(root / 'runtime.failure.json', {
            'error_type': type(error).__name__, 'error': str(error),
            'command': getattr(error, 'cmd', None), 'returncode': getattr(error, 'returncode', None),
            'stdout': getattr(error, 'stdout', None), 'stderr': getattr(error, 'stderr', None),
            'identities_observed': runtime_identity})
        raise SystemExit('Runtime qualification failed; inspect runtime.failure.json')
    C.write_json(root / 'runtime.json', runtime_identity)
    results = []
    failed = False
    for assembly in ['GRCh37', 'GRCh38']:
        out = root / assembly
        out.mkdir()
        try:
            result = qualify_assembly(args, assembly, source, runtime, out, runtime_identity)
            failed |= not result['direct_original_exact']
            print('{assembly}: {rows} rows; direct original exact={direct_original_exact}; published fixture identical={published_fixture_identical}'.format(**result))
        except Exception as error:
            failed = True
            result = {'assembly': assembly, 'status': 'failed', 'error_type': type(error).__name__, 'error': str(error),
                      'command': getattr(error, 'cmd', None), 'returncode': getattr(error, 'returncode', None),
                      'completed_artifacts': [file_identity(p) for p in out.rglob('*') if p.is_file()]}
            C.write_json(out / 'failure.json', result)
            print(assembly + ': failed; all available artifacts and details retained', file=sys.stderr)
        results.append(result)
    C.write_json(root / 'qualification.json', {'schema_version': 1, 'scope': 'upstream-fixture-real-model-calculation-only',
                                             'family_acceptance': False, 'annotation_qualification': False,
                                             'direct_original_exact': not failed, 'results': results})
    if failed:
        raise SystemExit('Qualification failed; inspect every assembly result/failure and mismatch artifact')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['source', 'python2', 'software-id', 'runtime-image', 'output']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--standalone-pin', help='Optional qualified standalone executable/source pin JSON')
    qualify(parser.parse_args())


if __name__ == '__main__':
    main()
