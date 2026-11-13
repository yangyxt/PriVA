#!/usr/bin/env python3
"""Materialize the explicitly selected domain-sorted-rescue-v1 Python oracle.

python das_sorted_reference.py render --source SCRIPT --output NEW_SCRIPT

The installer keeps the original Python source and its completion identities.
This test/reference procedure changes only the two rescue candidate iterators,
preserving every other source byte and verifying the complete AST. It never
executes DAS or opens an installed cache. Generated references belong outside
the source checkout. Pin their reported source/output SHA256 in the run profile.
"""
from __future__ import annotations
import argparse
import ast
import copy
import hashlib
import json
from pathlib import Path

PROFILE = 'domain-sorted-rescue-v1'
FUNCTIONS = ('rescue_bimodal_go_mf_domains', 'rescue_non_gomf_domains')


def render_source(raw: bytes) -> bytes:
    tree = ast.parse(raw)
    expected = copy.deepcopy(tree)
    lines = raw.splitlines(keepends=True)
    offsets = [0]
    for line in lines:
        offsets.append(offsets[-1] + len(line))
    changes = []
    for name in FUNCTIONS:
        functions = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name]
        if len(functions) != 1:
            raise ValueError('Expected one original rescue procedure: ' + name)
        loops = [n for n in ast.walk(functions[0]) if isinstance(n, ast.For)
                 and isinstance(n.target, ast.Name) and n.target.id == 'domain_path'
                 and isinstance(n.iter, ast.Name) and n.iter.id == 'tolerant_domains']
        if len(loops) != 1:
            raise ValueError('Expected one unchanged candidate iterator: ' + name)
        iterator = loops[0].iter
        start = offsets[iterator.lineno-1] + iterator.col_offset
        end = offsets[iterator.end_lineno-1] + iterator.end_col_offset
        if raw[start:end] != b'tolerant_domains':
            raise ValueError('Unexpected candidate iterator source bytes')
        changes.append((start, end, b'sorted(tolerant_domains)'))
        expected_function = next(n for n in expected.body if isinstance(n, ast.FunctionDef) and n.name == name)
        expected_loop = next(n for n in ast.walk(expected_function) if isinstance(n, ast.For)
                             and isinstance(n.target, ast.Name) and n.target.id == 'domain_path'
                             and isinstance(n.iter, ast.Name) and n.iter.id == 'tolerant_domains')
        expected_loop.iter = ast.Call(func=ast.Name(id='sorted', ctx=ast.Load()), args=[expected_loop.iter], keywords=[])
    output = raw
    for start, end, replacement in sorted(changes, reverse=True):
        output = output[:start] + replacement + output[end:]
    if ast.dump(ast.parse(output), include_attributes=False) != ast.dump(expected, include_attributes=False):
        raise ValueError('Sorted reference changed more than the two approved iterators')
    return output


def render(args):
    before = args.source.read_bytes()
    output = render_source(before)
    if args.source.read_bytes() != before:
        raise RuntimeError('Reference source changed while preparing the profile')
    with args.output.open('xb') as stream:
        stream.write(output)
    print(json.dumps({'profile': PROFILE, 'source_sha256': hashlib.sha256(before).hexdigest(),
                      'output_sha256': hashlib.sha256(output).hexdigest(), 'approved_iterator_changes': 2}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='procedure', required=True)
    p = sub.add_parser('render')
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.set_defaults(function=render)
    args = parser.parse_args(); args.function(args)


if __name__ == '__main__':
    main()
