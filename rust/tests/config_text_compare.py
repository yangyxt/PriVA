#!/usr/bin/env python3
"""Compare shared configuration text output with both unchanged Python writers.

python config_text_compare.py compare --repo FROZEN_SOURCE --probe PROBE \
    --work-root RESULTS --random-cases 20000

Only synthetic text enters this test. Imports call no installer procedure,
download or file publication. The shared JSON parser's existing integer,
surrogate and nesting bounds are outside this helper's supported JSON domain.
All returned output bytes and error class/messages inside that domain match
exactly; mismatches are retained without normalization or tolerance.
"""
import argparse
import hashlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import time

sys.dont_write_bytecode = True


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def cases(random_cases):
    keys = ['a', '', 'a.b', 'a+b?', 'é', 'a:b', 'a\nb', 'a\t', '#key']
    spaces = ['', ' ', '\t', '  ', '\u001c', '\u0085', '\u00a0', '\u2028']
    scalars = ['', 'old', 'value#same', 'v # c', '"x #y"', '"x\\"y"',
               "'x''#y'", '"unfinished # c', "'unfinished # c", '[]', 'null',
               '"closed" trailing', "''", '""', "'''", '"x\\\\"']
    for key, space, scalar, ending in itertools.product(keys, spaces, scalars, ['', '\n', '\r\n']):
        yield {'text': key + ' \t: ' + scalar + space + '# comment' + ending,
               'updates': {key: 'new/λ "file"'}}
    for key, count, indent in itertools.product(keys, [1, 2, 3], ['', ' ', '\t']):
        yield {'text': (indent + key + ': old\n') * count, 'updates': {key: 'new'}}
    values = [None, True, False, 0, 2**64-1, -(2**63), -0.0, 1.25, 1e-7,
              float('nan'), float('inf'), -float('inf'), 'λé😀', '\u001c\n\t',
              [1, 'λ', None], {'b': 1, 'a': 'é'}]
    originals = ['{}', ' { }\r\n', '{"z":1,"a":2,"z":3}', '{"λ": [1, null]}',
                 '{"old":NaN,"n":-0.0}', '{"a": 1,}', '{not json}', '{',
                 '\u00a0{}', '\u001c{}', '\ufeff{}']
    for original, value in itertools.product(originals, values):
        yield {'text': original, 'updates': {'z': value, 'a': 'é'}}
    for original in ['', '\n', '\r\n', '# hi\n', 'a: old', 'a: old\r\n\n', '{}']:
        yield {'text': original, 'updates': {}}
        yield {'text': original, 'updates': {'z': 'last', 'a': 'first'}}
    rng = random.Random(475012)
    alphabet = 'abc:#\'"\\ \t\r\u001c\u00a0\u2028λ😀'
    for _ in range(random_cases):
        original = ''.join(rng.choice(alphabet + '\n') for _ in range(rng.randrange(0, 100)))
        if rng.randrange(3):
            key = rng.choice(keys)
            original = original + '\n' + key + ': ' + ''.join(rng.choice(alphabet) for _ in range(rng.randrange(25)))
        else:
            key = rng.choice(keys)
        updates = {key: rng.choice(values)}
        if rng.randrange(2):
            updates[rng.choice(keys)] = rng.choice(values)
        yield {'text': original, 'updates': updates}


def compare(args):
    args.work_root.mkdir(parents=True, exist_ok=False)
    files = [args.repo / 'scripts/installer_config_paths.py', args.repo / 'scripts/reference_resources.py',
             args.repo / 'rust/src/config_text.rs', args.repo / 'rust/src/json_text.rs',
             args.repo / 'rust/src/python_text.rs', args.probe]
    pins = {str(path): digest(path) for path in files}
    paths = load(files[0], 'config_paths_reference')
    reference = load(files[1], 'reference_resources_reference')
    writers = {'paths': paths.render_config, 'reference': reference.render_config_values}
    totals = {'cases': 0, 'successes': 0, 'failures': 0, 'mismatches': 0}
    before = time.monotonic()
    with (args.work_root / 'probe.stderr.log').open('x') as stderr, \
         (args.work_root / 'mismatches.jsonl').open('x') as mismatches:
        probe = subprocess.Popen([str(args.probe), 'serve'], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                 stderr=stderr, text=True, encoding='utf-8', bufsize=1)
        try:
            for row in cases(args.random_cases):
                for mode, writer in writers.items():
                    request = dict(row, mode=mode)
                    try:
                        expected = {'ok': True, 'output': writer(row['text'], row['updates'])}
                    except Exception as error:
                        expected = {'ok': False, 'kind': type(error).__name__, 'error': str(error)}
                    probe.stdin.write(json.dumps(request, ensure_ascii=True) + '\n')
                    probe.stdin.flush()
                    line = probe.stdout.readline()
                    if not line:
                        raise RuntimeError('native probe ended before its response')
                    actual = json.loads(line)
                    totals['cases'] += 1
                    totals['successes' if expected['ok'] else 'failures'] += 1
                    if expected != actual:
                        totals['mismatches'] += 1
                        mismatches.write(json.dumps({'request': request, 'expected': expected, 'actual': actual},
                                                    ensure_ascii=True) + '\n')
            probe.stdin.close()
            returncode = probe.wait(timeout=30)
        finally:
            if probe.poll() is None:
                probe.kill()
                probe.wait()
    stable = all(digest(Path(path)) == sha for path, sha in pins.items())
    report = dict(totals, job=os.environ.get('SLURM_JOB_ID'), elapsed_s=time.monotonic()-before,
                  source_hashes=pins, source_stable=stable, probe_returncode=returncode,
                  python_version=sys.version, random_seed=475012,
                  scope='exact text/error comparison; no installer caller activated',
                  passed=totals['mismatches'] == 0 and returncode == 0 and stable)
    (args.work_root / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    if not report['passed']:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    command = sub.add_parser('compare')
    for key in ('repo', 'probe', 'work-root'):
        command.add_argument('--' + key, type=Path, required=True)
    command.add_argument('--random-cases', type=int, default=20000)
    command.set_defaults(function=compare)
    args = parser.parse_args()
    args.function(args)


if __name__ == '__main__':
    main()
