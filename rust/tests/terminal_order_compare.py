#!/usr/bin/env python3
"""Exact differential check of the opt-in R8 permutation/rank boundary.

compare --binary FILE --binary-sha256 SHA --oracle ROOT
        --oracle-manifest-sha256 SHA --output NEW_DATA_ROOT
The original payload/axes are gathered using Rust's permutation, never inferred
from TSV. Full output values/dtypes/null masks/float bits/attrs are compared.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import subprocess

import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(args):
    assert sha(args.binary) == args.binary_sha256
    assert sha(args.oracle / 'manifest.json') == args.oracle_manifest_sha256
    manifest = json.loads((args.oracle / 'manifest.json').read_text())
    for name, wanted in manifest.items():
        assert sha(args.oracle / name) == wanted, name
    oracle = json.loads((args.oracle / 'oracle.json').read_text())
    harness = Path(oracle['full_oracle']) / 'source/ranking_control_oracle.py'
    spec = importlib.util.spec_from_file_location('frozen_full_oracle', harness)
    helper = importlib.util.module_from_spec(spec); spec.loader.exec_module(helper)
    spec = importlib.util.spec_from_file_location('captured_R8_oracle', args.oracle / 'capture-harness.py')
    transport = importlib.util.module_from_spec(spec); spec.loader.exec_module(transport)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'comparison-harness.py').write_bytes(Path(__file__).read_bytes())
    cases, failures = [], []
    try:
        for case in oracle['cases']:
            source = args.oracle / case['name']
            output = args.output / (case['name'] + '.csv')
            metadata = args.output / (case['name'] + '.metadata.csv')
            axis = sha(source / 'entry.json')
            argv = [str(args.binary), 'terminal-order', '--input', str(source / 'input.csv'),
                    '--output', str(output), '--metadata-output', str(metadata), '--axis-metadata', axis]
            for position, dtype in enumerate(case['dtypes']):
                argv += ['--column-dtype', f'{position}={dtype}']
            before = sha(source / 'input.csv')
            process = subprocess.run(argv, text=True, capture_output=True)
            (args.output / (case['name'] + '.stderr')).write_text(process.stderr)
            report = dict(name=case['name'], returncode=process.returncode, successful=False)
            cases.append(report)
            assert sha(source / 'input.csv') == before
            if case['name'] == 'full_rank_chrom_dtype_mixed':
                assert process.returncode and 'UnsupportedInput' in process.stderr and 'retain Python' in process.stderr
                assert not output.exists() and not metadata.exists()
                report.update(successful=True, deferred_to_python=True)
                continue
            if case['original_exception'] is not None:
                expected = case['original_exception']
                assert process.returncode and expected['type'] + ': ' + expected['message'] in process.stderr, (case['name'], process.stderr, expected)
                assert not output.exists() and not metadata.exists()
                report.update(successful=True, original_failure_preserved=True)
                continue
            assert process.returncode == 0, (case['name'], process.stderr)
            with output.open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            positions = [int(row['source_position']) for row in rows]
            ranks = [int(row['variant_rank']) for row in rows]
            assert positions == case['permutation'] and ranks == case['ranks'], case['name']
            with (source / 'input.csv').open(newline='') as stream:
                original = list(csv.reader(stream))[1:]
            for row, position in zip(rows, positions):
                assert [row['index_kind'], row['index_value']] == original[position][:2]
            maxima = json.loads((source / 'maxima.json').read_text())[0]
            for row, position in zip(rows, positions):
                expected = transport.cell(transport.atom(maxima['values'][position]))
                assert [row['max_kind'], row['max_value']] == expected, (case['name'], position, row, expected)
            with metadata.open(newline='') as stream:
                records = list(csv.DictReader(stream))
            columns = [record for record in records if record['record'] == 'column']
            assert [record['name'] for record in columns] == case['output_columns']
            assert [record['dtype'] for record in columns] == case['output_dtypes']
            assert [int(record['position']) for record in columns] == list(range(len(columns)))
            assert next(record for record in records if record['record'] == 'axis')['name'] == axis
            entry = pickle.loads((source / 'entry.pkl').read_bytes())
            temporary = next(record for record in records if record['record'] == 'temporary')
            assert temporary['dtype'] == maxima['dtype']['text']
            expected_position = list(entry.columns).index('max_variant_score') if 'max_variant_score' in entry else len(entry.columns)
            assert int(temporary['position']) == expected_position
            # Original sort_values returns an empty copy before taking rows;
            # iloc[[]] would expand its RangeIndex and lose source metadata.
            gathered = entry.copy(deep=True) if len(entry) == 0 else entry.iloc[positions].copy(deep=True)
            gathered['variant_rank'] = np.asarray(ranks, dtype='int64')
            if 'max_variant_score' in gathered:
                gathered = gathered.drop(columns='max_variant_score')
            expected = pickle.loads((source / 'expected.pkl').read_bytes())
            pd.testing.assert_frame_equal(gathered, expected, check_exact=True)
            assert helper.encode(gathered) == helper.encode(expected), case['name']
            assert sha(source / 'entry.pkl') == case['entry_sha256']
            report.update(successful=True, full_typed_output_exact=True, index_tokens_exact=True,
                          temporary_max_bits_exact=True, input_unchanged=True)
    except Exception as error:
        failures.append(dict(type=type(error).__name__, message=str(error), case=cases[-1]['name'] if cases else None))
        raise
    finally:
        result = dict(binary_sha256=args.binary_sha256, oracle_manifest_sha256=args.oracle_manifest_sha256,
            harness_sha256=sha(Path(__file__)), native_environment={name: os.environ.get(name) for name in
                ['SLURM_JOB_ID', 'SLURM_STEP_ID', 'SLURMD_NODENAME']},
            cases=cases, failures=failures, successful=not failures and len(cases) == len(oracle['cases']),
            actual_outputs=sum(x.get('full_typed_output_exact', False) for x in cases),
            original_failures=sum(x.get('original_failure_preserved', False) for x in cases),
            deferred_to_python=sum(x.get('deferred_to_python', False) for x in cases),
            scientific_acceptance=False, production_backend_selected=False)
        (args.output / 'comparison.json').write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps({key: value for key, value in result.items() if key != 'cases'}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest='command', required=True)
    procedure = command.add_parser('compare')
    for name in ['binary', 'oracle', 'output']:
        procedure.add_argument('--' + name, type=Path, required=True)
    for name in ['binary-sha256', 'oracle-manifest-sha256']:
        procedure.add_argument('--' + name, required=True)
    compare(parser.parse_args())


if __name__ == '__main__':
    main()
