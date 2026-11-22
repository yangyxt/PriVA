#!/usr/bin/env python3
"""Compare typed R5->R6 plans against actual original current-buffer stages.

compare --binary FILE --binary-sha256 SHA --oracle ROOT
        --oracle-manifest-sha256 SHA --output NEW_DATA_ROOT
No parser/numeric substitute is used for unsupported inputs: run the unchanged
source block once on the original current buffers, without earlier-stage replay.
"""
import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import subprocess
import warnings

import numpy as np
import pandas as pd
from pandas.core.arrays.floating import NUMPY_FLOAT_TO_DTYPE


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source_context(helper, scoring, full_oracle, output):
    """Run one real separate source fixture before any paired current buffer.

    The captured oracle had already run full-source nullable cases. Their R8
    sort populates a process-global pandas Float64 dtype cache. Reproduce that
    native context through its actual producer, never a cache reset/write or a
    replay of earlier stages on the paired R5-entry buffers.
    """
    registry = NUMPY_FLOAT_TO_DTYPE[np.dtype('float64')]
    def cache():
        return {key: dict(type=type(value).__name__, repr=repr(value))
                for key, value in getattr(registry, '_cache', {}).items()}
    case = next(case for case in json.loads((full_oracle / 'capture.json').read_text())['cases']
                if case['name'] == 'rank_nullable_Float64' and case['kind'] == 'ranking')
    folder = full_oracle / case['folder']
    before = cache()
    data, pedigree = pickle.loads((folder / 'prepared-input.pkl').read_bytes())
    helper.encode(data)
    helper.encode(pedigree)
    bound = json.loads((folder / 'traced/arguments.json').read_text())['arguments']
    kwargs = {key['value']: helper.decode_atom(value) for key, value in bound['items']
              if key['value'] not in {'df', 'ped_df'}}
    calls = []
    def trace(frame, event, arg):
        if event == 'call' and frame.f_code.co_name == 'index_class' and frame.f_locals.get('self') is registry:
            stack, current = [], frame
            while current is not None:
                stack.append(dict(function=current.f_code.co_name, file=current.f_code.co_filename, line=current.f_lineno))
                current = current.f_back
            calls.append(stack)
        return trace
    previous = __import__('sys').gettrace()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            __import__('sys').settrace(trace)
            value = scoring.sort_and_rank_variants(data, pedigree, **kwargs)
        finally:
            __import__('sys').settrace(previous)
    assert helper.encode(value) == json.loads((folder / 'traced/outcome.json').read_text())['output']
    assert 'index_class' in cache()
    assert 'index_class' in before or calls
    result = dict(case=case['name'], arguments_sha256=sha(folder / 'traced/arguments.json'),
        prepared_input_sha256=sha(folder / 'prepared-input.pkl'), outcome_sha256=sha(folder / 'traced/outcome.json'),
        registry_before=before, registry_after=cache(), actual_source_registry_calls=calls,
        warnings=[dict(category=item.category.__name__, message=str(item.message)) for item in caught],
        distinct_fixture_buffers=True, paired_R5_buffers_not_created_yet=True,
        before_any_Rust_case=True, caches_not_reset=True, no_manual_cache_warming=True)
    (output / 'source-context-bootstrap.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def compare(args):
    assert sha(args.binary) == args.binary_sha256 and sha(args.oracle / 'manifest.json') == args.oracle_manifest_sha256
    manifest = json.loads((args.oracle / 'manifest.json').read_text())
    for name, wanted in manifest.items():
        assert sha(args.oracle / name) == wanted, name
    oracle = json.loads((args.oracle / 'oracle.json').read_text())
    harness = load('frozen_R5_R6_oracle', args.oracle / 'capture-harness.py')
    helper = load('full_ranker_R5_R6_helper', Path(oracle['full_oracle']) / 'source/ranking_control_oracle.py')
    scoring, _, _, _ = helper.original_modules(Path(oracle['full_oracle']) / 'source/scripts')
    original, _ = harness.excerpt(Path(oracle['full_oracle']) / 'source/scripts/acmg_scoring.py', scoring)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'comparison-harness.py').write_bytes(Path(__file__).read_bytes())
    reports, failures = [], []
    context = None
    try:
        context = source_context(helper, scoring, Path(oracle['full_oracle']), args.output)
        for case in oracle['cases']:
            source = args.oracle / case['name']
            output, metadata = [args.output / (case['name'] + suffix) for suffix in ['.csv', '.metadata.csv']]
            # Construct the genuine prepared graph, then perform the producing
            # entry observation exactly once. No nullable-array preflight or
            # cache reset/rehydration occurs before an unsupported-stage replay.
            inputs = pickle.loads((source / 'prepared-input.pkl').read_bytes())
            entry = helper.encode(inputs)
            before = pickle.dumps(inputs, protocol=5)
            assert hashlib.sha256(before).hexdigest() == case['input_before_sha256'], case['name']
            frame, _, _, mask = inputs
            axis = sha(source / 'input.pkl')
            argv = [str(args.binary), 'ranking-gene-weights', '--input', str(source / 'input.csv'),
                '--penalty-input', str(source / 'penalty.csv'), '--penalty-dtype', case['penalty_dtype'] or 'unsupported',
                '--row-index-kind', case['index_kind'], '--axis-metadata', axis,
                '--output', str(output), '--metadata-output', str(metadata)]
            for role, argument in zip(['dispensable', 'relevant'], case['gene_arguments']):
                if argument['kind'] == 'none':
                    argv += ['--' + role + '-kind', 'none']
                elif argument['kind'] in ['str', 'path']:
                    argv += ['--' + role + '-kind', 'file', '--' + role + '-file', argument['value']]
                else:
                    argv += ['--' + role + '-kind', 'unsupported']
            for i, dtype in enumerate(case['dtypes']):
                argv += ['--column-dtype', f'{i}={dtype}']
            transport_before = [sha(source / name) for name in ['input.csv', 'penalty.csv']]
            process = subprocess.run(argv, capture_output=True, text=True)
            (args.output / (case['name'] + '.stderr')).write_text(process.stderr)
            report = dict(name=case['name'], successful=False, returncode=process.returncode)
            reports.append(report)
            assert [sha(source / name) for name in ['input.csv', 'penalty.csv']] == transport_before
            assert pickle.dumps(inputs, protocol=5) == before, case['name']
            expected_record = json.loads((source / 'result.json').read_text())
            if process.returncode:
                assert 'retain Python backend' in process.stderr and (
                    'UnsupportedInput' in process.stderr or 'unsupported typed cell' in process.stderr), (case['name'], process.stderr)
                assert not output.exists() and not metadata.exists()
                original_frame = id(frame)
                result, fb, fa, _ = harness.run_current(inputs, original, scoring, helper, entry_observation=entry)
                assert id(inputs[0]) == original_frame
                assert result == expected_record, case['name']
                assert fb == before and hashlib.sha256(fa).hexdigest() == case['input_after_sha256'], case['name']
                report.update(successful=True, deferred_to_python=True, original_current_buffer_stage_only_exact=True,
                    earlier_stages_not_replayed=True, preflight_current_buffers_unchanged=True,
                    fallback_source_failed=result['exception'] is not None,
                    original_failure_type_message_phase_and_partial_mutation_exact=True,
                    original_parser_and_warning_behavior_retained=True)
                continue
            assert case['original_exception'] is None, (case['name'], 'Rust must defer the original source failure')
            with output.open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            with metadata.open(newline='') as stream:
                meta = {row['name']: row['value'] for row in csv.DictReader(stream)}
            with (source / 'input.csv').open(newline='') as stream:
                source_rows = list(csv.reader(stream))[1:]
            assert meta['axis'] == axis and int(meta['score_column_position']) == case['columns'].index('sort_index')
            assert len(rows) == case['rows'] * int(meta['step_count'])
            actual = dict(zip(harness.INPUTS, inputs))
            events = []

            def snapshot(event, line):
                events.append(dict(event=event, line=line, snapshot=pickle.dumps(actual, protocol=5)))

            def floats(items, field):
                return np.array([int(row[field], 16) for row in items], dtype=np.uint64).view(np.float64)

            def series(values, index, name=None, attrs=None, duplicate=None):
                result = pd.Series(values, index=index, name=name)
                result.attrs = copy.deepcopy(frame.attrs if attrs is None else attrs)
                result.flags.allows_duplicate_labels = frame.flags.allows_duplicate_labels if duplicate is None else duplicate
                return result

            def stage(i, name):
                assert meta[f'step.{i}.name'] == name
                selected_rows = [row for row in rows if int(row['step']) == i]
                assert [int(row['source_position']) for row in selected_rows] == list(range(case['rows']))
                assert [[row['index_kind'], row['index_value']] for row in selected_rows] == [row[:2] for row in source_rows]
                return selected_rows

            position = 0
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                for role, target, factor, line in [('dispensable', 'dispensable_gene_index', .9, 604),
                                                    ('relevant', 'relevant_gene_index', 1.2, 611)]:
                    enabled = meta[role + '.enabled'] == 'true'
                    assert enabled == (actual[role + '_gene_list'] is not None)
                    if enabled:
                        values = [meta[f'{role}.gene.{i}'] for i in range(int(meta[role + '.count']))]
                        # Preserve CPython set construction at the typed boundary.
                        # Parsing/membership and all numeric values come from Rust.
                        actual[role + '_genes'] = set(values)
                        items = stage(position, target)
                        assert int(meta[f'step.{position}.line']) == line
                        target_position = list(frame.columns).index(target) if target in frame.columns else len(frame.columns)
                        assert int(meta[f'step.{position}.multiplier_column_position']) == target_position
                        frame[target] = 1.0
                        source_series = frame['SYMBOL']
                        membership = series(np.array([bool(int(row['selected'])) for row in items], dtype=bool),
                            frame.index, name='SYMBOL', attrs=source_series.attrs,
                            duplicate=source_series.flags.allows_duplicate_labels)
                        frame.loc[membership, target] = factor
                        assert np.array(frame[target], dtype=np.float64).view(np.uint64).tolist() == [int(row['multiplier_bits'], 16) for row in items]
                        frame.loc[:, 'sort_index'] = series(floats(items, 'score_bits'), frame.index)
                        position += 1
                    snapshot('line', 606 if role == 'dispensable' else 613)
                items = stage(position, 'zygosity_inheritance_mechanism_compatibility')
                assert int(meta[f'step.{position}.line']) == 613
                assert int(meta[f'step.{position}.multiplier_column_position']) == list(frame.columns).index('zygosity_inheritance_mechanism_compatibility')
                selected = np.array([bool(int(row['selected'])) for row in items], dtype=bool)
                product = series(floats(items, 'score_bits')[selected], frame.index[selected])
                frame.loc[mask, 'sort_index'] = product
                snapshot('return', 617)
                assert position + 1 == int(meta['step_count'])
            after = pickle.dumps(inputs, protocol=5)
            assert hashlib.sha256(after).hexdigest() == case['input_after_sha256'], case['name']
            assert helper.encode(actual) == expected_record['output'], case['name']
            assert helper.encode(inputs) == expected_record['input_after'], case['name']
            assert [dict(event=item['event'], line=item['line'], locals=helper.encode(pickle.loads(item['snapshot'])))
                    for item in events] == expected_record['source_events'], case['name']
            assert [dict(category=item.category.__name__, message=str(item.message)) for item in caught] == [
                {key: item[key] for key in ['category', 'message']} for item in case['warnings']], case['name']
            assert not case['source_logs'], case['name']
            report.update(successful=True, full_typed_output_exact=True, raw_float_bits_and_nulls_exact=True,
                ordered_intermediate_states_exact=True, opaque_axes_metadata_and_column_positions_exact=True,
                source_input_mutation_and_pickle_exact=True, warnings_and_stage_order_exact=True,
                preflight_current_buffers_unchanged=True, prior_stages_not_replayed=True)
        # Supported coverage cannot be replaced by an all-Python result.
        exact_names = {item['name'] for item in reports if item.get('full_typed_output_exact')}
        assert {'extra_none', 'extra_both', 'extra_same_file', 'extra_parse_empty', 'extra_parse_comments',
                'extra_axis_duplicate', 'extra_axis_gapped', 'extra_penalty_masked_raw_0', 'extra_penalty_masked_raw_1'} <= exact_names
        for name, wanted in manifest.items():
            assert sha(args.oracle / name) == wanted, name
    except Exception as error:
        failures.append(dict(type=type(error).__name__, message=str(error), case=reports[-1]['name'] if reports else None))
        raise
    finally:
        final = dict(binary_sha256=args.binary_sha256, oracle_manifest_sha256=args.oracle_manifest_sha256,
            harness_sha256=sha(Path(__file__)), native_environment={key: os.environ.get(key) for key in ['SLURM_JOB_ID', 'SLURM_STEP_ID', 'SLURMD_NODENAME']},
            successful=not failures and len(reports) == len(oracle['cases']), cases=reports, failures=failures,
            actual_outputs=sum(case.get('full_typed_output_exact', False) for case in reports),
            original_failures=0, deferred_to_python=sum(case.get('deferred_to_python', False) for case in reports),
            deferred_original_failures=sum(case.get('fallback_source_failed', False) for case in reports),
            comparison='R5 then R6 only; exact sequential score bits/intermediate states/column plans and genuine current-buffer parser/failure/partial-mutation fallback',
            separate_original_source_context_bootstrap=context,
            scientific_acceptance=False, production_backend_selected=False)
        (args.output / 'comparison.json').write_text(json.dumps(final, indent=2) + '\n')
        print(json.dumps({key: value for key, value in final.items() if key != 'cases'}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest='command', required=True).add_parser('compare')
    for name in ['binary', 'oracle', 'output']:
        command.add_argument('--' + name, type=Path, required=True)
    for name in ['binary-sha256', 'oracle-manifest-sha256']:
        command.add_argument('--' + name, required=True)
    compare(parser.parse_args())


if __name__ == '__main__':
    main()
