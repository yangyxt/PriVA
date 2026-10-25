#!/usr/bin/env python3
"""Exactly compare opt-in R7 plans and genuine current-buffer stage fallback.

compare --binary FILE --binary-sha256 SHA --oracle ROOT
        --oracle-manifest-sha256 SHA --output NEW_DATA_ROOT
Retain source locals, logs, failures, warnings, immediate pickles and raw bits.
Preflight uses only closed transports; paired pandas buffers are constructed
after the native result, so unsupported inputs cannot be warmed/normalized.
"""
import argparse
import ast
import copy
import csv
import hashlib
import importlib.util
import json
import logging
import os
from pathlib import Path
import pickle
import subprocess
import warnings

import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def compare(args):
    assert sha(args.binary) == args.binary_sha256
    assert sha(args.oracle / 'manifest.json') == args.oracle_manifest_sha256
    for name, wanted in json.loads((args.oracle / 'manifest.json').read_text()).items():
        assert sha(args.oracle / name) == wanted, name
    oracle = json.loads((args.oracle / 'oracle.json').read_text())
    root = Path(oracle['root_oracle'])
    assert sha(root / 'manifest.json') == oracle['root_manifest_sha256']
    assert sha(root / 'capture-harness.py') == oracle['original_harness_sha256']
    harness = load('unchanged_root_R7', root / 'capture-harness.py')
    full = Path(oracle['full_oracle'])
    helper = load('full_R7_helper', full / 'source/ranking_control_oracle.py')
    helper.verify_sources(full, harness.FULL_PINS)
    scoring, _, _, _ = helper.original_modules(full / 'source/scripts')
    scoring.logger.setLevel(logging.INFO)
    original, source = harness.excerpt(full / 'source/scripts/acmg_scoring.py', scoring)
    assert source == oracle['source']
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'comparison-harness.py').write_bytes(Path(__file__).read_bytes())
    # Same separate genuine full-source context as the producing captures.
    # This never replays any earlier stage on the paired R7-entry buffers.
    context_path = Path(oracle['context_source'])
    assert sha(context_path) == oracle['context_source_sha256']
    context = load('qualified_R7_source_context', context_path)
    context.source_context(helper, scoring, full, args.output)
    expressions = ast.parse((full / 'source/scripts/acmg_scoring.py').read_text())
    power = next(node for node in ast.walk(expressions) if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute) and isinstance(node.func.value, ast.Name)
        and node.func.value.id == 'np' and node.func.attr == 'power' and node.lineno == 651)
    parent = {id(child): node for node in ast.walk(expressions) for child in ast.iter_child_nodes(node)}
    product, total = parent[id(power)], parent[id(parent[id(power)])]
    assert isinstance(product, ast.BinOp) and isinstance(total, ast.BinOp)
    math = [compile(ast.Expression(copy.deepcopy(node)), str(full / 'source/scripts/acmg_scoring.py'), 'eval')
            for node in [power, product, total]]
    reports, failures = [], []
    try:
        for case in oracle['cases']:
            source_folder = args.oracle / case['name']
            expected = json.loads((source_folder / 'result.json').read_text())
            output = args.output / (case['name'] + '.csv')
            metadata = args.output / (case['name'] + '.metadata.csv')
            selected = case['selected_column']
            source_position = case['columns'].index(selected) if selected in case['columns'] else 0
            axis = sha(source_folder / 'input.pkl')
            argv = [str(args.binary), 'ranking-pext', '--input', str(source_folder / 'input.csv'),
                '--options-input', str(source_folder / 'options.csv'), '--coercion-input', str(source_folder / 'coercion-input.csv'),
                '--coercion-column-position', str(source_position), '--mask-input', str(source_folder / 'mask.csv'),
                '--mask-dtype', case['mask_dtype'], '--row-index-kind', case['index_kind'], '--axis-metadata', axis,
                '--output', str(output), '--metadata-output', str(metadata)]
            for i, dtype in enumerate(case['dtypes']):
                argv += ['--column-dtype', str(i) + '=' + dtype]
            files = ['input.csv', 'options.csv', 'coercion-input.csv', 'mask.csv']
            transports_before = [sha(source_folder / name) for name in files]
            process = subprocess.run(argv, capture_output=True, text=True)
            (args.output / (case['name'] + '.stderr')).write_text(process.stderr)
            report = dict(name=case['name'], returncode=process.returncode, successful=False)
            reports.append(report)
            assert [sha(source_folder / name) for name in files] == transports_before
            # No inspection/coercion of actual current buffers before preflight.
            inputs = pickle.loads((source_folder / 'prepared-input.pkl').read_bytes())
            frame, tissues, cutoff, floor, shape, mask = inputs
            if process.returncode:
                assert 'retain Python backend' in process.stderr and any(text in process.stderr for text in
                    ['UnsupportedInput', 'unsupported typed cell']), (case['name'], process.stderr)
                assert not output.exists() and not metadata.exists()
                identity = id(frame)
                result, before, after, payload = harness.run_current(inputs, original, scoring, helper)
                assert id(inputs[0]) == identity
                assert result == expected, case['name']
                assert before == (source_folder / 'input.pkl').read_bytes() and after == (source_folder / 'input-after.pkl').read_bytes(), case['name']
                if payload is not None:
                    assert payload == (source_folder / 'expected.pkl').read_bytes(), case['name']
                report.update(successful=True, deferred_to_python=True, source_failed=result['exception'] is not None,
                    original_stage_once_on_current_buffers=True, prior_stages_not_replayed=True,
                    no_preflight_buffer_observation=True, original_failure_phase_warning_partial_mutation_exact=True)
                continue
            assert case['original_exception'] is None, (case['name'], 'source failures must defer')
            entry = helper.encode(inputs)
            before = pickle.dumps(inputs, protocol=5)
            assert hashlib.sha256(before).hexdigest() == expected['input_before_sha256'], case['name']
            with output.open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            with metadata.open(newline='') as stream:
                meta = {row['name']: row['value'] for row in csv.DictReader(stream)}
            assert meta['axis'] == axis
            applied = meta['applied'] == 'true'
            assert applied == (selected is not None)
            actual = dict(zip(harness.INPUTS, inputs))
            actual['pext_col_to_use'] = None
            if meta['tissue_list_present'] == 'true':
                # Reify source string identities on the original arguments.
                # split/strip can reuse an unchanged token/whole argument.
                # Selection itself remains the checked native plan.
                actual['pext_tissue_list'] = [t.strip() for t in tissues.split(',') if t.strip()]
                assert actual['pext_tissue_list'] == [meta['tissue.' + str(i)] for i in range(int(meta['tissue_count']))]
            if meta['first_tissue_column_present'] == 'true':
                actual['first_tissue_col'] = f"PEXT_{actual['pext_tissue_list'][0].replace('-', '_')}"
                assert actual['first_tissue_col'] == meta['first_tissue_column']
            if applied:
                actual['pext_col_to_use'] = actual['first_tissue_col'] if meta['used_tissue'] == 'true' else next(
                    value for value in original.__code__.co_consts if isinstance(value, str) and value == 'PEXT_MEAN')
                assert actual['pext_col_to_use'] == meta['column']
            logs = harness.Logs(helper)
            scoring.logger.addHandler(logs)
            def floats(field):
                return np.array([int(row[field], 16) for row in rows], dtype=np.uint64).view(np.float64)
            def series(values, source_series, name):
                result = pd.Series(values, index=source_series.index, name=name)
                result.attrs = copy.deepcopy(source_series.attrs)
                result.flags.allows_duplicate_labels = source_series.flags.allows_duplicate_labels
                return result
            try:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    if applied:
                        assert len(rows) == len(frame)
                        assert [int(row['source_position']) for row in rows] == list(range(len(frame)))
                        with (source_folder / 'input.csv').open(newline='') as stream:
                            tokens = [row[:2] for row in list(csv.reader(stream))[1:]]
                        assert [[row['index_kind'], row['index_value']] for row in rows] == tokens
                        if meta['used_tissue'] == 'true':
                            scoring.logger.info(f'Using tissue-specific pext column for sorting: {selected}')
                        else:
                            scoring.logger.info('Using PEXT_MEAN column for expression-aware sorting')
                        assert int(meta['source_column_position']) == list(frame.columns).index(selected)
                        assert int(meta['score_column_position']) == list(frame.columns).index('sort_index')
                        column = frame[selected]
                        actual['pext_numeric'] = series(floats('pext_numeric_bits'), column, column.name)
                        actual['scaled_pext'] = series(floats('scaled_pext_bits'), actual['pext_numeric'], column.name)
                        assert actual['pext_numeric'].name is column.name and actual['scaled_pext'].name is column.name
                        target_position = list(frame.columns).index('pext_sort_index') if 'pext_sort_index' in frame.columns else len(frame.columns)
                        assert int(meta['target_column_position']) == target_position
                        frame['pext_sort_index'] = floats('pext_sort_index_bits')
                        # Retain source-equivalent loc row construction/metadata;
                        # all numeric values and mask selections come from Rust.
                        left = frame.loc[mask, 'sort_index']
                        right = frame.loc[mask, 'pext_sort_index']
                        chosen = np.array([bool(int(row['selected'])) for row in rows], dtype=bool)
                        assert len(left) == len(right) == int(chosen.sum())
                        assert helper.encode(left.index) == helper.encode(frame.index[chosen])
                        value = series(floats('score_bits')[chosen], left, None)
                        frame.loc[mask, 'sort_index'] = value
                        scoring.logger.info('Applied pext-based sort modulation using %s (cutoff=%s, floor=%s, shape=%s)',
                            selected, cutoff, floor, shape)
                    else:
                        assert not rows
            finally:
                scoring.logger.removeHandler(logs)
            after = pickle.dumps(inputs, protocol=5)
            actual = harness.selected(actual)
            payload = pickle.dumps(actual, protocol=5)
            assert before == (source_folder / 'input.pkl').read_bytes(), case['name']
            assert after == (source_folder / 'input-after.pkl').read_bytes(), case['name']
            assert payload == (source_folder / 'expected.pkl').read_bytes(), case['name']
            assert helper.encode(actual) == expected['output'], case['name']
            assert helper.encode(inputs) == expected['input_after'] and entry == expected['input_before'], case['name']
            assert logs.records == expected['source_logs'], case['name']
            assert [dict(category=x.category.__name__, message=str(x.message), source_line=x.lineno) for x in caught] == expected['warnings'], case['name']
            assert [dict(event='return', line=666, locals=helper.encode(actual))] == expected['source_events'], case['name']
            if applied:
                closed_locals = pickle.loads((source_folder / 'expected.pkl').read_bytes())
                for expression, field in zip(math, ['power_bits', 'product_bits', 'sum_bits']):
                    values = eval(expression, dict(scoring.__dict__), closed_locals)
                    assert np.asarray(values, dtype=np.float64).view(np.uint64).tolist() == [int(row[field], 16) for row in rows], (case['name'], field)
            report.update(successful=True, actual_Rust_output=True, exact_current_buffer_pickles=True,
                exact_locals_raw_bits_dtypes_axes_attrs_flags=True, exact_logs_warnings_and_return_phase=True,
                scalar_parser_and_operation_order_exact=True, numerical_tolerance=False)
    except Exception as error:
        failures.append(dict(type=type(error).__name__, message=str(error)))
        raise
    finally:
        result = dict(successful=not failures and len(reports) == len(oracle['cases']) and all(r['successful'] for r in reports),
            cases=reports, failures=failures, oracle_manifest_sha256=args.oracle_manifest_sha256,
            source_code_sha256=sha(Path(__file__)), binary_sha256=args.binary_sha256,
            actual_outputs=sum(r.get('actual_Rust_output', False) for r in reports),
            original_failures=0, deferred_to_python=sum(r.get('deferred_to_python', False) for r in reports),
            deferred_original_failures=sum(r.get('source_failed', False) for r in reports),
            native={k: os.environ.get(k) for k in ['SLURM_JOB_ID', 'SLURM_STEP_ID', 'SLURMD_NODENAME']},
            no_production_activation=True, scientific_acceptance=False)
        (args.output / 'comparison.json').write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps({k: v for k, v in result.items() if k != 'cases'}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest='procedure', required=True).add_parser('compare')
    for name in ['binary', 'oracle', 'output']:
        command.add_argument('--' + name, type=Path, required=True)
    for name in ['binary-sha256', 'oracle-manifest-sha256']:
        command.add_argument('--' + name, required=True)
    compare(parser.parse_args())


if __name__ == '__main__':
    main()
