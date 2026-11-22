#!/usr/bin/env python3
"""Compare the R1 plan with actual-source outputs, failures and input identity.

compare --binary FILE --binary-sha256 SHA --oracle ROOT
        --oracle-manifest-sha256 SHA --output NEW_ROOT
Applies Rust's plan to original payload buffers; no Python BIOTYPE comparison or
score arithmetic is performed by the adapter.
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
    assert sha(args.binary)==args.binary_sha256 and sha(args.oracle/'manifest.json')==args.oracle_manifest_sha256
    for name,wanted in json.loads((args.oracle/'manifest.json').read_text()).items(): assert sha(args.oracle/name)==wanted,name
    oracle=json.loads((args.oracle/'oracle.json').read_text())
    spec=importlib.util.spec_from_file_location('frozen_full_oracle',Path(oracle['full_oracle'])/'source/ranking_control_oracle.py')
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'comparison-harness.py').write_bytes(Path(__file__).read_bytes())
    reports,failures=[],[]
    try:
        for case in oracle['cases']:
            source=args.oracle/case['name'];output=args.output/(case['name']+'.csv');metadata=args.output/(case['name']+'.metadata.csv')
            axis=sha(source/'input.json')
            argv=[str(args.binary),'rank-setup','--input',str(source/'input.csv'),'--output',str(output),
                '--metadata-output',str(metadata),'--axis-metadata',axis]
            for position,dtype in enumerate(case['dtypes']):argv+=['--column-dtype',f'{position}={dtype}']
            input_sha=sha(source/'input.csv')
            process=subprocess.run(argv,capture_output=True,text=True)
            (args.output/(case['name']+'.stderr')).write_text(process.stderr)
            report=dict(name=case['name'],returncode=process.returncode,successful=False);reports.append(report)
            assert sha(source/'input.csv')==input_sha
            affected_duplicates=any(case['columns'].count(name)>1 for name in ['ACMG_quant_score','BIOTYPE','sort_index'])
            if affected_duplicates:
                assert process.returncode and 'UnsupportedInput' in process.stderr and 'retain Python' in process.stderr
                assert not output.exists() and not metadata.exists()
                report.update(successful=True,deferred_to_python=True);continue
            if case['original_exception'] is not None:
                error=case['original_exception']
                assert process.returncode and error['type']+': '+error['message'] in process.stderr,(case['name'],process.stderr,error)
                assert not output.exists() and not metadata.exists()
                report.update(successful=True,original_failure_preserved=True);continue
            assert process.returncode==0,(case['name'],process.stderr)
            with output.open(newline='') as stream:rows=list(csv.DictReader(stream))
            positions=[int(row['source_position']) for row in rows]
            assert positions==case['permutation'],case['name']
            with (source/'input.csv').open(newline='') as stream:input_rows=list(csv.reader(stream))[1:]
            for row,position in zip(rows,positions):assert [row['index_kind'],row['index_value']]==input_rows[position][:2]
            with metadata.open(newline='') as stream:actions=list(csv.DictReader(stream))
            columns=[row for row in actions if row['record']=='column']
            assert [row['name'] for row in columns]==case['output_columns']
            assert [row['dtype'] for row in columns]==case['output_dtypes']
            assert [int(row['position']) for row in columns]==list(range(len(columns)))
            assert next(row for row in actions if row['record']=='axis')['name']==axis
            original=pickle.loads((source/'input.pkl').read_bytes())
            before=pickle.dumps(original,protocol=5)
            assert hashlib.sha256(before).hexdigest()==case['input_before_sha256']
            dropped=[row for row in actions if row['record']=='drop']
            for row in dropped:
                position=int(row['position'])
                assert position==int(row['source_column'])
                assert original.columns[position]==row['name']
                assert row['name'] in ['haplo_insuf_index','zygosity_inheritance_mechanism_index']
            copy_actions=[row for row in columns if row['operation']=='copy_score']
            assert len(copy_actions)==1
            action=copy_actions[0];source_name=original.columns[int(action['source_column'])]
            assert source_name=='ACMG_quant_score' and action['name']=='sort_index'
            source_columns=[int(row['source_column']) for row in columns]
            for row,position in zip(columns,source_columns):
                assert row['dtype']==str(original.dtypes.iloc[position])
                if row['operation']=='retain':assert row['name']==original.columns[position]
                else:assert row['operation']=='copy_score'
            # Apply every emitted column position to original buffers, including
            # unrelated duplicate labels and the raw score copy. No Python R1
            # drop/assignment/filter implementation supplies the plan here.
            gathered=original.copy(deep=True).iloc[:,source_columns].copy(deep=True)
            gathered.columns=pd.Index([row['name'] for row in columns],dtype=original.columns.dtype,name=original.columns.name)
            keep=np.zeros(len(original),dtype=bool);keep[positions]=True
            gathered=gathered.loc[keep,:]
            after=pickle.dumps(original,protocol=5)
            assert before==after and case['input_pickle_unchanged'],case['name']
            expected=pickle.loads((source/'expected.pkl').read_bytes())
            pd.testing.assert_frame_equal(gathered,expected,check_exact=True)
            assert helper.encode(gathered)==helper.encode(expected),case['name']
            report.update(successful=True,full_typed_output_exact=True,input_pickle_unchanged=True,index_tokens_exact=True,
                          emitted_column_plan_applied=True)
    except Exception as error:
        failures.append(dict(type=type(error).__name__,message=str(error),case=reports[-1]['name'] if reports else None));raise
    finally:
        result=dict(binary_sha256=args.binary_sha256,oracle_manifest_sha256=args.oracle_manifest_sha256,
            harness_sha256=sha(Path(__file__)),native_environment={name:os.environ.get(name) for name in ['SLURM_JOB_ID','SLURM_STEP_ID','SLURMD_NODENAME']},
            cases=reports,failures=failures,successful=not failures and len(reports)==len(oracle['cases']),
            actual_outputs=sum(x.get('full_typed_output_exact',False) for x in reports),
            original_failures=sum(x.get('original_failure_preserved',False) for x in reports),
            deferred_to_python=sum(x.get('deferred_to_python',False) for x in reports),
            scientific_acceptance=False,production_backend_selected=False)
        (args.output/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({key:value for key,value in result.items() if key!='cases'}),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    command=parser.add_subparsers(dest='command',required=True).add_parser('compare')
    for name in ['binary','oracle','output']:command.add_argument('--'+name,type=Path,required=True)
    for name in ['binary-sha256','oracle-manifest-sha256']:command.add_argument('--'+name,required=True)
    compare(parser.parse_args())


if __name__=='__main__':main()
