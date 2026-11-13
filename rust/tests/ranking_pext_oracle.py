#!/usr/bin/env python3
"""Capture exact R7 boxed scalar inputs and original current-buffer outcomes.

capture --oracle CLOSED_ROOT_R7 --full-oracle CLOSED_FULL_RANKER --output NEW_ROOT
inspect --oracle CLOSED_TRANSPORT_ROOT

The original root-owned oracle and scientific source are immutable. A Python
profile observer records arguments to the actual coerce_pext_value calls;
neither pandas boxing nor the scientific R7 statements are reimplemented.
Runtime, binary transports and pickles belong on the data side, inside Slurm.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import logging
import os
from pathlib import Path
import pickle
import struct
import sys

import numpy as np
import pandas as pd

D = Path('/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade')
ROOT_MANIFEST = 'f1cb80551e98c4961a20318e1b1e8a4128999f88bc32f8378cc5c2d6200e76a5'
CONTEXT_SOURCE = Path('/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/rust/tests/ranking_gene_weights_compare.py')
CONTEXT_SHA = 'eb47e2dcfd747bde6ff313d87c4ab366d7a4e9fe4fad110653b05233399f29c1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def cell(value, scalar_origin=False):
    if value is None:
        return ['null', '']
    if value is pd.NA:
        return ['pd_na', '']
    if scalar_origin and isinstance(value, np.bool_):
        return ['numpy_bool', str(int(value))]
    if scalar_origin and isinstance(value, np.float32):
        return ['float32', struct.pack('>f', value).hex()]
    if scalar_origin and isinstance(value, np.floating) and not isinstance(value, np.float64):
        return ['unsupported', type(value).__name__]
    if isinstance(value, (bool, np.bool_)):
        return ['bool', str(int(value))]
    if isinstance(value, (int, np.integer)):
        if value < -(2**63) or value > 2**64 - 1:
            return ['bigint', str(value)]
        return ['uint' if value > 2**63 - 1 else 'int', str(value)]
    if isinstance(value, (float, np.floating)):
        return ['float', struct.pack('>d', value).hex()]
    if isinstance(value, str):
        try:
            return ['str', value.encode().hex()]
        except UnicodeEncodeError:
            pass
    return ['unsupported', type(value).__module__ + '.' + type(value).__qualname__]


def observed_run(inputs, original, scoring, helper, harness):
    boxed, selection = [], []
    parser = scoring.coerce_pext_value
    def profile(frame, event, arg):
        if event == 'call' and frame.f_code is parser.__code__:
            boxed.append(frame.f_locals['value'])
        if event == 'return' and frame.f_code is original.__code__:
            selection.append(frame.f_locals.get('pext_col_to_use'))
    previous = sys.getprofile()
    try:
        sys.setprofile(profile)
        actual = harness.run_current(inputs, original, scoring, helper)
    finally:
        sys.setprofile(previous)
    assert len(selection) == 1
    return actual, boxed, selection[0]


def export(root, name, prepared, original, scoring, helper, harness, origin, expected=None):
    folder = root / name
    folder.mkdir()
    (folder / 'prepared-input.pkl').write_bytes(prepared)
    actual, boxed, selection = observed_run(pickle.loads(prepared), original, scoring, helper, harness)
    repeated, _, _ = observed_run(pickle.loads(prepared), original, scoring, helper, harness)
    assert actual == repeated, name
    result, before, after, output = actual
    if expected is not None:
        assert result == expected, name
    for filename, payload in [('input.pkl', before), ('input-after.pkl', after), ('expected.pkl', output)]:
        if payload is not None:
            (folder / filename).write_bytes(payload)
    write(folder / 'result.json', result)
    frame, tissues, cutoff, floor, shape, mask = pickle.loads(before)
    token = lambda value: ['opaque', helper.canonical(helper.encode(value)).hex()]
    with (folder / 'input.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value'] + [item for label in frame.columns for item in [str(label) + '_kind', str(label)]])
        for row in range(len(frame)):
            writer.writerow(token(frame.index[row]) + [item for col in range(frame.shape[1]) for item in cell(frame.iloc[row, col])])
    with (folder / 'mask.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value', 'value_kind', 'value', 'raw_data', 'is_null'])
        if isinstance(mask, pd.Series):
            array = mask.array
            for row, value in enumerate(mask):
                writer.writerow(token(mask.index[row]) + cell(value) + [
                    str(int(array._data[row])) if hasattr(array, '_mask') else '',
                    str(int(array._mask[row])) if hasattr(array, '_mask') else '0'])
    with (folder / 'options.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['name', 'kind', 'value'])
        for key, value in zip(['tissues', 'cutoff', 'floor', 'shape'], [tissues, cutoff, floor, shape]):
            writer.writerow([key] + cell(value))
    with (folder / 'coercion-input.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['source_position', 'value_kind', 'value'])
        for position, value in enumerate(boxed):
            writer.writerow([position] + cell(value, scalar_origin=True))
    report = dict(name=name, origin=origin, rows=len(frame), columns=[str(x) for x in frame.columns],
        dtypes=[str(x) for x in frame.dtypes], mask_dtype=str(mask.dtype) if isinstance(mask, pd.Series) else 'unsupported',
        index_kind='multi' if isinstance(frame.index, pd.MultiIndex) else 'flat', selected_column=selection,
        coerce_arguments_observed=len(boxed), scalar_types=[type(x).__module__ + '.' + type(x).__qualname__ for x in boxed],
        scalar_kinds=[cell(x, scalar_origin=True)[0] for x in boxed], original_exception=result['exception'],
        warnings=result['warnings'], input_before_sha256=result['input_before_sha256'],
        input_after_sha256=result['input_after_sha256'], source_logs=result['source_logs'])
    return report


def extra_cases(harness):
    # Actual pandas map boxing distinguishes object numpy_bool from native bool.
    for dtype in ('object', 'bool', 'boolean', 'float32', 'float64', 'Float64'):
        values = [True, False, True, False] if dtype in ('bool', 'boolean') else [.02, np.nan, -0., .1]
        inputs = list(harness.base())
        inputs[0]['PEXT_MEAN'] = pd.Series(values, dtype=dtype)
        yield 'boxing_' + dtype, tuple(inputs)
    inputs = list(harness.base())
    inputs[0]['PEXT_MEAN'] = pd.Series([np.bool_(False), False, np.float32(.02), np.int64(-1)], dtype=object)
    yield 'boxing_numpy_objects', tuple(inputs)
    for text in ('\u001c Liver \u001f, Brain', '\u00a0Brain-Cortex\u2003', 'Missing,Liver', ' , ,'):
        inputs = list(harness.base())
        inputs[0]['PEXT_Liver'] = [.01] * 4
        inputs[0]['PEXT_Brain_Cortex'] = [.03] * 4
        inputs[1] = text
        yield 'tissue_' + text.encode().hex(), tuple(inputs)
    vals = ['-0', '0', '-0&0', '0&-0', '\u0660.\u0660\u0662', '0.0_2', '-NAN&.02', '.02&-NAN',
            'INF', '-INF', 'NA', '.', 'NAN', np.float64(-0.)]
    inputs = list(harness.base(vals))
    inputs[0]['sort_index'] = np.array([-.0, .0, np.inf, -np.inf, .95, .475, -.1] * 2)
    yield 'signed_zero_unicode_parser', tuple(inputs)
    for shape in (.5, np.nextafter(.5, 0.), np.nextafter(.5, 1.), 1., 2., .3, 1.5, np.finfo(float).tiny):
        values = np.array([0., -0., np.nextafter(0., 1.), .0000001, .02,
            np.nextafter(.1, 0.), .1, np.nextafter(.1, 1.), .5, 1.])
        inputs = list(harness.base(values))
        inputs[0]['sort_index'] = np.linspace(-.1, .95, len(values))
        inputs[4] = float(shape)
        yield 'power_' + struct.pack('>d', shape).hex(), tuple(inputs)
    for floor in (0., -0., 1., np.nextafter(1., 0.)):
        inputs = list(harness.base([0., -0., .02, .1]))
        inputs[0]['sort_index'] = [-0., 0., np.inf, -np.inf]
        inputs[3] = float(floor)
        yield 'floor_' + struct.pack('>d', floor).hex(), tuple(inputs)
    for cutoff in (np.nextafter(0., 1.), np.finfo(float).tiny, 1., np.finfo(float).max):
        inputs = list(harness.base([0., -0., .02, .1]))
        inputs[2] = float(cutoff)
        yield 'cutoff_' + struct.pack('>d', cutoff).hex(), tuple(inputs)
    for selected in (False, True):
        inputs = list(harness.base())
        inputs[-1] = pd.Series([selected] * 4, name='mask')
        inputs[0].attrs = {'scope': 'R7 attrs', 'nested': [1, 2]}
        inputs[0]['PEXT_MEAN'].attrs = {'scope': 'cached source-column attrs'}
        yield 'attrs_' + str(selected), tuple(inputs)
    inputs = list(harness.base())
    inputs[0].drop(columns=['PEXT_MEAN'], inplace=True)
    inputs[2:5] = [None, 'invalid', -1.]
    inputs[-1] = pd.Series([1, 2, 3, 4], dtype='int64')
    yield 'no_column_skips_parameter_mask_score_validation', tuple(inputs)


def capture(args):
    assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_STEP_ID')
    for path in (args.oracle, args.full_oracle, args.output):
        assert path.resolve().is_relative_to(D)
    assert sha(args.oracle / 'manifest.json') == ROOT_MANIFEST
    for name, wanted in json.loads((args.oracle / 'manifest.json').read_text()).items():
        assert sha(args.oracle / name) == wanted, name
    args.output.mkdir(parents=True, exist_ok=False)
    os.chdir(args.output)
    (args.output / 'capture-harness.py').write_bytes(Path(__file__).read_bytes())
    write(args.output / 'attempt.json', dict(native={k: os.environ.get(k) for k in ['SLURM_JOB_ID', 'SLURM_STEP_ID']}, source_sha256=sha(Path(__file__))))
    harness = load('original_root_R7', args.oracle / 'capture-harness.py')
    assert sha(args.full_oracle / 'source-pins.json') == harness.FULL_PINS
    assert sha(args.full_oracle / 'capture.json') == harness.FULL_CAPTURE
    helper = load('full_source_R7_helper', args.full_oracle / 'source/ranking_control_oracle.py')
    helper.verify_sources(args.full_oracle, harness.FULL_PINS)
    scoring, _, _, _ = helper.original_modules(args.full_oracle / 'source/scripts')
    scoring.logger.setLevel(logging.INFO)
    original, source = harness.excerpt(args.full_oracle / 'source/scripts/acmg_scoring.py', scoring)
    assert sha(CONTEXT_SOURCE) == CONTEXT_SHA
    context = load('qualified_source_context', CONTEXT_SOURCE)
    context.source_context(helper, scoring, args.full_oracle, args.output)
    reports = []
    root_oracle = json.loads((args.oracle / 'oracle.json').read_text())
    for case in root_oracle['cases']:
        folder = args.oracle / case['name']
        reports.append(export(args.output, case['name'], (folder / 'prepared-input.pkl').read_bytes(),
            original, scoring, helper, harness, dict(kind='verified-root-original-R7', case=case['name']),
            json.loads((folder / 'result.json').read_text())))
    for name, inputs in extra_cases(harness):
        reports.append(export(args.output, 'edge_' + name, pickle.dumps(inputs, protocol=5),
            original, scoring, helper, harness, dict(kind='additional-original-R7-before-Rust')))
    write(args.output / 'runtime-identities.json', helper.runtime_identity())
    write(args.output / 'oracle.json', dict(source_commit=helper.COMMIT, source=source,
        root_oracle=str(args.oracle), root_manifest_sha256=ROOT_MANIFEST, full_oracle=str(args.full_oracle),
        context_source=str(CONTEXT_SOURCE), context_source_sha256=CONTEXT_SHA, harness_sha256=sha(Path(__file__)),
        original_harness_sha256=sha(args.oracle / 'capture-harness.py'), cases=reports,
        preceding_full_source_failures=root_oracle['preceding_full_source_failures'],
        native={k: os.environ.get(k) for k in ['SLURM_JOB_ID', 'SLURM_STEP_ID', 'SLURMD_NODENAME']},
        production_backend_selected=False, scientific_acceptance=False))
    write(args.output / 'manifest.json', {str(p.relative_to(args.output)): sha(p) for p in sorted(args.output.rglob('*')) if p.is_file()})
    print(json.dumps(dict(cases=len(reports), original_outputs=sum(c['original_exception'] is None for c in reports),
        original_failures=sum(c['original_exception'] is not None for c in reports), manifest_sha256=sha(args.output / 'manifest.json'))))


def inspect(args):
    assert args.oracle.resolve().is_relative_to(D)
    for name, wanted in json.loads((args.oracle / 'manifest.json').read_text()).items():
        assert sha(args.oracle / name) == wanted, name
    oracle = json.loads((args.oracle / 'oracle.json').read_text())
    print(json.dumps(dict(cases=len(oracle['cases']), manifest_sha256=sha(args.oracle / 'manifest.json'), all_bytes_verified=True)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='procedure', required=True)
    command = commands.add_parser('capture')
    for name in ('oracle', 'full-oracle', 'output'):
        command.add_argument('--' + name, type=Path, required=True)
    commands.add_parser('inspect').add_argument('--oracle', type=Path, required=True)
    args = parser.parse_args()
    {'capture': capture, 'inspect': inspect}[args.procedure](args)


if __name__ == '__main__':
    main()
