#!/usr/bin/env python3
"""Capture native-integer R8 behavior before extending the Rust score domain.

capture --base-oracle CLOSED_R8_ROOT --output NEW_DATA_ROOT
One real full-ranker Int64 checkpoint is retained. Extra integer cases execute
only the original verbatim R8 statements, explicitly not earlier ranking math.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

D = Path('/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade')
BASE_MANIFEST = '579b02b4afa4677b247cc446eed3e01d6749929a6bc3a5cc932fd3483fb0420c'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def capture(args):
    assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_STEP_ID')
    assert args.base_oracle.resolve().is_relative_to(D) and args.output.resolve().is_relative_to(D)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'integer-capture.py').write_bytes(Path(__file__).read_bytes())
    native = {name: os.environ.get(name) for name in ['SLURM_JOB_ID', 'SLURM_STEP_ID', 'SLURMD_NODENAME']}
    (args.output / 'attempt.json').write_text(json.dumps(dict(native_environment=native, source_sha256=sha(Path(__file__))), indent=2)+'\n')
    assert sha(args.base_oracle / 'manifest.json') == BASE_MANIFEST
    base = json.loads((args.base_oracle / 'oracle.json').read_text())
    original = load('r8_preimplementation_oracle', args.base_oracle / 'capture-harness.py')
    full = Path(base['full_oracle'])
    helper = load('unchanged_full_oracle', full / 'source/ranking_control_oracle.py')
    helper.verify_sources(full, original.FULL_PINS)
    excerpt, details = original.r8_function(full / 'source/scripts/acmg_scoring.py')
    (args.output / 'capture-harness.py').write_bytes((args.base_oracle / 'capture-harness.py').read_bytes())
    original.write(args.output / 'runtime-identities.json', helper.runtime_identity())
    reports = []
    source = args.base_oracle / 'full_rank_nullable_Int64'
    reports.append(original.export_case(args.output, 'original_full_rank_nullable_Int64',
        pickle.loads((source / 'entry.pkl').read_bytes()), pickle.loads((source / 'expected.pkl').read_bytes()),
        helper, excerpt, dict(kind='actual-closed-full-ranker-R8-checkpoint', base_oracle_manifest_sha256=BASE_MANIFEST)))
    for dtype in ['int64', 'Int64', 'uint64', 'UInt64']:
        signed = dtype.lower() == 'int64'
        values = ([2**53+1, 2**53, 2**63-1, -2**63, -1, 0, 2**53+2, 2**53+3] if signed else
                  [2**53+1, 2**53, 2**64-1, 2**63, 1, 0, 2**53+2, 2**53+3])
        for missing_coordinate in [False, True]:
            frame = helper.base_frame(len(values))
            frame['sort_index'] = pd.Series(values, dtype=dtype)
            frame['ACMG_quant_score'] = pd.Series(list(reversed(values)), dtype=dtype)
            frame['pos'] = pd.Series([1, 1, 2, 2, 3, 3, 4, 5], dtype='Int64')
            if missing_coordinate:
                frame.loc[4, 'pos'] = pd.NA
            frame.index = pd.Index([7, 7, -2, 40, -2, 9, 9, 100], name='integer checkpoint')
            name = f'integer_{dtype}_missing_coordinate_{missing_coordinate}'
            reports.append(original.export_case(args.output, name, frame, None, helper, excerpt,
                dict(kind='verbatim-R8-excerpt-only', source_promotion_observed_not_repaired=True)))
        if dtype in ['Int64', 'UInt64']:
            frame = helper.base_frame(6)
            frame['sort_index'] = pd.Series([pd.NA, pd.NA, 2**53+1, pd.NA, 0, pd.NA], dtype=dtype)
            frame['ACMG_quant_score'] = pd.Series([1, pd.NA, 3, 3, 0, pd.NA], dtype=dtype)
            frame['pos'] = pd.Series([1, 1, 2, 2, pd.NA, 3], dtype='Int64')
            reports.append(original.export_case(args.output, 'integer_' + dtype + '_all_missing_groups', frame, None,
                helper, excerpt, dict(kind='verbatim-R8-excerpt-only', nullable_mask_preserved=True)))
        frame = helper.base_frame(0)
        frame['sort_index'] = pd.Series([], dtype=dtype)
        frame['ACMG_quant_score'] = pd.Series([], dtype=dtype)
        reports.append(original.export_case(args.output, 'integer_' + dtype + '_empty', frame, None,
            helper, excerpt, dict(kind='verbatim-R8-excerpt-only')))
    original.write(args.output / 'oracle.json', dict(native_environment=native, created_utc=helper.utc(),
        source_commit=helper.COMMIT, full_oracle=str(full), base_oracle_manifest_sha256=BASE_MANIFEST,
        harness_sha256=sha(Path(__file__)), excerpt=details, cases=reports,
        scope='native-integer R8 supplement; one actual full-ranker checkpoint, additional original R8 excerpt cases only',
        scientific_acceptance=False, production_backend_selected=False))
    original.write(args.output / 'manifest.json', {str(path.relative_to(args.output)): sha(path) for path in sorted(args.output.rglob('*')) if path.is_file()})
    print(json.dumps(dict(cases=len(reports), oracle_sha256=sha(args.output / 'oracle.json'),
        manifest_sha256=sha(args.output / 'manifest.json'), native_environment=native,
        maxima_dtypes={item['name']: json.loads((args.output / item['name'] / 'maxima.json').read_text())[0]['dtype']['text'] for item in reports})))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest='command', required=True).add_parser('capture')
    command.add_argument('--base-oracle', type=Path, required=True)
    command.add_argument('--output', type=Path, required=True)
    capture(parser.parse_args())


if __name__ == '__main__':
    capture_module = None
    main()
