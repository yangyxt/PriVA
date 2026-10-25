#!/usr/bin/env python3
"""Compare opt-in R3 bits, in-place assignment and current-buffer fallback.

compare --binary FILE --binary-sha256 SHA --oracle ROOT
        --oracle-manifest-sha256 SHA --output NEW_ROOT
Fallback executes only original L528 on post-C1-C5 buffers; earlier stages,
warnings and numerical-library calls are never replayed by the fallback.
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


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def compare(args):
    assert sha(args.binary)==args.binary_sha256 and sha(args.oracle/'manifest.json')==args.oracle_manifest_sha256
    for name,wanted in json.loads((args.oracle/'manifest.json').read_text()).items():assert sha(args.oracle/name)==wanted
    oracle=json.loads((args.oracle/'oracle.json').read_text())
    harness=load('actual_R3_oracle',args.oracle/'capture-harness.py')
    helper=load('actual_full_oracle',Path(oracle['full_oracle'])/'source/ranking_control_oracle.py')
    original,_=harness.original(Path(oracle['full_oracle'])/'source/scripts/acmg_scoring.py')
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'comparison-harness.py').write_bytes(Path(__file__).read_bytes())
    reports,failures=[],[]
    try:
        for case in oracle['cases']:
            source=args.oracle/case['name'];output=args.output/(case['name']+'.csv');metadata=args.output/(case['name']+'.metadata.csv')
            # Recreate the producing observer's actual entry, rather than
            # rehydrate its later snapshot. Pickling can sever dtype-cache key
            # aliases to module strings even when entry bytes compare equal.
            # No cache reset, dtype coercion or earlier source stage is run.
            data=pickle.loads((source/'prepared-input.pkl').read_bytes())
            helper.encode(data)
            before=pickle.dumps(data,protocol=5)
            assert hashlib.sha256(before).hexdigest()==case['input_before_sha256']
            index_kind='multi' if isinstance(data.index,pd.MultiIndex) else 'flat'
            axis=sha(source/'input.json')
            argv=[str(args.binary),'control-multiply','--input',str(source/'input.csv'),'--output',str(output),
                '--metadata-output',str(metadata),'--row-index-kind',index_kind,'--axis-metadata',axis]
            for i,dtype in enumerate(case['dtypes']):argv+=['--column-dtype',f'{i}={dtype}']
            input_sha=sha(source/'input.csv')
            process=subprocess.run(argv,capture_output=True,text=True)
            (args.output/(case['name']+'.stderr')).write_text(process.stderr)
            report=dict(name=case['name'],successful=False,returncode=process.returncode,
                genuine_prepared_source_entry_reconstructed=True);reports.append(report)
            assert sha(source/'input.csv')==input_sha
            reason=None
            if index_kind=='multi':reason='MultiIndex assignment warning route'
            elif any(case['columns'].count(name)>1 for name in ['sort_index','control_common_index']):reason='ambiguous required labels'
            elif all(name in case['columns'] for name in ['sort_index','control_common_index']):
                if any(case['dtypes'][case['columns'].index(name)]!='float64' for name in ['sort_index','control_common_index']):
                    reason='unqualified source dtype/promotion route'
                elif np.any(np.isnan(data['sort_index'].to_numpy(copy=False)) & np.isnan(data['control_common_index'].to_numpy(copy=False))):
                    reason='two-NaN operand payload route'
            assert pickle.dumps(data,protocol=5)==before,case['name']
            report['fallback_preflight_input_pickle_unchanged']=True
            if reason:
                assert process.returncode and 'UnsupportedInput' in process.stderr and 'retain Python backend' in process.stderr,(case['name'],process.stderr)
                assert not output.exists() and not metadata.exists()
                # Stage-only fallback on current post-C1-C5 buffers; never replay
                # control/PED/scoring stages or their warnings/mutations.
                result,fallback_before,fallback_after=harness.run(data,original,helper)
                assert result==json.loads((source/'result.json').read_text()),case['name']
                assert hashlib.sha256(fallback_before).hexdigest()==case['input_before_sha256']
                assert hashlib.sha256(fallback_after).hexdigest()==case['input_after_sha256']
                report.update(successful=True,deferred_to_python=True,reason=reason,
                    original_current_buffer_stage_only_exact=True,earlier_stages_not_replayed=True,
                    fallback_source_failed=result['exception'] is not None)
                continue
            if case['original_exception'] is not None:
                error=case['original_exception']
                assert process.returncode and error['type']+': '+error['message'] in process.stderr,(case['name'],process.stderr,error)
                assert 'source_failure_line=528' in process.stderr
                assert not output.exists() and not metadata.exists()
                assert pickle.dumps(data,protocol=5)==before
                report.update(successful=True,original_failure_preserved=True,source_failure_line=528)
                continue
            assert process.returncode==0,(case['name'],process.stderr)
            with output.open(newline='') as stream:rows=list(csv.DictReader(stream))
            assert [int(row['source_position']) for row in rows]==list(range(case['rows']))
            with (source/'input.csv').open(newline='') as stream:source_rows=list(csv.reader(stream))[1:]
            assert [[row['index_kind'],row['index_value']] for row in rows]==[row[:2] for row in source_rows]
            with metadata.open(newline='') as stream:meta={row['name']:row['value'] for row in csv.DictReader(stream)}
            target=int(meta['target_source_column_position'])
            assert target==case['columns'].index('sort_index') and meta['dtype']=='float64' and meta['axis']==axis
            values=np.array([int(row['float_bits'],16) for row in rows],dtype=np.uint64).view(np.float64)
            product=pd.Series(values,index=data.index,name=None)
            product.attrs=copy.deepcopy(data.attrs);product.flags.allows_duplicate_labels=data.flags.allows_duplicate_labels
            same_object=id(data)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                data.loc[:,'sort_index']=product
            assert id(data)==same_object
            actual_warnings=[dict(category=item.category.__name__,message=str(item.message)) for item in caught]
            assert actual_warnings==case['warnings'],case['name']
            # The oracle snapshots immediately after the statement, before
            # subsequent semantic observers can warm any pandas dtype caches.
            after=pickle.dumps(data,protocol=5)
            assert hashlib.sha256(after).hexdigest()==case['input_after_sha256'],case['name']
            expected=pickle.loads((source/'expected.pkl').read_bytes())
            pd.testing.assert_frame_equal(data,expected,check_exact=True)
            assert helper.encode(data)==helper.encode(expected),case['name']
            assert (before==after)==case['source_input_pickle_unchanged']
            report.update(successful=True,full_typed_output_exact=True,raw_float_bits_exact=True,
                source_input_mutation_and_pickle_exact=True,source_column_position_exact=True,
                original_same_frame_assignment_and_warnings_exact=True)
    except Exception as error:
        failures.append(dict(type=type(error).__name__,message=str(error),case=reports[-1]['name'] if reports else None));raise
    finally:
        final=dict(binary_sha256=args.binary_sha256,oracle_manifest_sha256=args.oracle_manifest_sha256,
            harness_sha256=sha(Path(__file__)),native_environment={key:os.environ.get(key) for key in ['SLURM_JOB_ID','SLURM_STEP_ID','SLURMD_NODENAME']},
            successful=not failures and len(reports)==len(oracle['cases']),cases=reports,failures=failures,
            actual_outputs=sum(case.get('full_typed_output_exact',False) for case in reports),
            original_failures=sum(case.get('original_failure_preserved',False) for case in reports),
            deferred_to_python=sum(case.get('deferred_to_python',False) for case in reports),
            deferred_original_failures=sum(case.get('fallback_source_failed',False) for case in reports),
            comparison='strict native-float64 L528 bits/in-place assignment/metadata/mutation; current-buffer original-stage fallback counted separately',
            scientific_acceptance=False,production_backend_selected=False)
        (args.output/'comparison.json').write_text(json.dumps(final,indent=2)+'\n')
        print(json.dumps({key:value for key,value in final.items() if key!='cases'}),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    command=parser.add_subparsers(dest='command',required=True).add_parser('compare')
    for name in ['binary','oracle','output']:command.add_argument('--'+name,type=Path,required=True)
    for name in ['binary-sha256','oracle-manifest-sha256']:command.add_argument('--'+name,required=True)
    compare(parser.parse_args())


if __name__=='__main__':main()
