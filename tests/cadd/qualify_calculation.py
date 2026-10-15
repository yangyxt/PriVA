#!/usr/bin/env python3
"""Real installed-model qualification on upstream fixtures, not family acceptance.

Run in the allocation with an installed original CADD Python2 runtime (or pinned
container shim). This uses bundled models/annotations and requires no annotation
cache downloads. Generated files belong in the caller's durable validation area.
"""
import argparse
import contextlib
import gzip
import importlib.util
import io
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


def qualify(args):
    source = Path(args.source).resolve(strict=True)
    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    runtime = str(Path(args.python2).absolute())
    if not Path(runtime).is_file():
        raise FileNotFoundError(runtime)
    results = []
    for assembly in ['GRCh37', 'GRCh38']:
        out = root / assembly
        out.mkdir()
        profile = json.loads((ROOT / 'conf/cadd' / (assembly + '.reference.json')).read_text())
        profile.pop('annotation')
        profile['software_id'] = args.software_id
        profile['resource_manifest_path'] = str(out / 'resources.json')
        env = profile['calculation']['env']
        env.update(CADD_SOURCE_ROOT=str(source), CADD_PYTHON2=runtime,
                   CADD_IMPUTE_CONFIG=str(source / 'config' / ('impute_' + assembly + '_v1.7.cfg')),
                   CADD_MODEL=str(source / 'data/models' / assembly / ('CADDv1.7-' + assembly + '.mod')),
                   CADD_PHRED_TABLE=str(source / 'data/models' / assembly / ('conversionTable_CADDv1.7-' + assembly + '.txt')))
        for item, key in zip(profile['calculation']['verify_files'],
                             ['CADD_IMPUTE_CONFIG', 'CADD_MODEL', 'CADD_PHRED_TABLE']):
            item['path'] = env[key]
        profile['calculation']['resource_roots'] = [item['path'] for item in profile['calculation']['verify_files']]
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
        # Independent direct invocation of the four original scripts, not the
        # adapter or its Bash wrapper. Intermediate files preserve stream bytes.
        scripts = source / 'src/scripts'
        features = out / 'direct.features.csv.gz'
        with gzip.open(bundle / 'annotations.tsv.gz', 'rb') as handle:
            # GzipFile.fileno() addresses compressed bytes; decompress explicitly
            # into a normal file before passing it as subprocess standard input.
            plain = out / 'direct.anno.tsv'
            plain.write_bytes(handle.read())
        with open(plain, 'rb') as handle:
            invoke(runtime, scripts / 'trackTransformation.py', ['-b', '-c', env['CADD_IMPUTE_CONFIG'],
                   '-o', str(features), '--noheader'], stdin=handle)
        raw = out / 'direct.raw.tsv'
        with open(raw, 'wb') as handle:
            invoke(runtime, scripts / 'predictSKmodel.py', ['-i', str(features), '-m', env['CADD_MODEL'],
                   '-a', str(bundle / 'annotations.tsv.gz')], stdout=handle)
        hierarchy = out / 'direct.hierarchy.tsv'
        with open(raw, 'rb') as src, open(hierarchy, 'wb') as dest:
            invoke(runtime, scripts / 'max_line_hierarchy.py', ['--all'], stdin=src, stdout=dest)
        direct = out / 'direct.scores.tsv'
        with open(hierarchy, 'rb') as src, open(direct, 'wb') as dest:
            invoke(runtime, scripts / 'appendPHREDscore.py', ['-t', env['CADD_PHRED_TABLE']], stdin=src, stdout=dest)
        actual = scored / 'scores.tsv'
        if actual.read_bytes() != direct.read_bytes():
            raise AssertionError('Adapter/direct original calculation differ for ' + assembly)
        output_rows = list(C.rows(actual, columns + ['RawScore', 'PHRED']))
        published_rows = list(C.rows(fixture, columns + ['RawScore', 'PHRED']))
        differing = sum(a != b for a, b in zip(output_rows, published_rows)) + abs(len(output_rows) - len(published_rows))
        results.append({'assembly': assembly, 'rows': count, 'direct_original_exact': True,
                        'published_fixture_identical': differing == 0, 'published_fixture_differing_rows': differing,
                        'fixture_sha256': C.digest(fixture), 'output_sha256': C.digest(actual),
                        'profile_sha256': C.digest(ppath), 'software_id': args.software_id,
                        'upstream_commit': C.PIN})
    C.write_json(root / 'qualification.json', {'schema_version': 1, 'scope': 'upstream-fixture-real-model-calculation-only',
                                             'family_acceptance': False, 'annotation_qualification': False, 'results': results})
    for result in results:
        print('{assembly}: {rows} rows; direct original exact={direct_original_exact}; published fixture identical={published_fixture_identical}'.format(**result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--python2', required=True)
    parser.add_argument('--software-id', required=True)
    parser.add_argument('--output', required=True)
    qualify(parser.parse_args())


if __name__ == '__main__':
    main()
