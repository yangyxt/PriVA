#!/usr/bin/env python3
"""Compare typed R2 context/masks, raw masked bits, warnings and source failures.

compare --binary FILE --binary-sha256 SHA --oracle ROOT
        --oracle-manifest-sha256 SHA --output NEW_ROOT
Original source locals/logs/pickle effects stay retained. Python dtype-cache
internals, arbitrary objects and declared fallback cases are never counted exact.
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
import re
import subprocess

import numpy as np
import pandas as pd


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(args):
    assert sha(args.binary)==args.binary_sha256 and sha(args.oracle/'manifest.json')==args.oracle_manifest_sha256
    for name,wanted in json.loads((args.oracle/'manifest.json').read_text()).items():assert sha(args.oracle/name)==wanted
    oracle=json.loads((args.oracle/'oracle.json').read_text())
    spec=importlib.util.spec_from_file_location('frozen_full_oracle',Path(oracle['full_oracle'])/'source/ranking_control_oracle.py')
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'comparison-harness.py').write_bytes(Path(__file__).read_bytes())
    reports,failures=[],[]
    try:
        for case in oracle['cases']:
            source=args.oracle/case['name'];result=args.output/(case['name']+'.csv');metadata=args.output/(case['name']+'.metadata.csv')
            axis=sha(source/'input.json')
            argv=[str(args.binary),'pedigree-context','--input',str(source/'input.csv'),'--family-input',str(source/'family.csv'),
                '--output',str(result),'--metadata-output',str(metadata),'--axis-metadata',axis]
            for i,dtype in enumerate(case['dtypes']):argv+=['--column-dtype',f'{i}={dtype}']
            if case['pedigree_columns'] is not None:
                argv+=['--pedigree-input',str(source/'pedigree.csv')]
                for i,dtype in enumerate(case['pedigree_transport_dtypes']):argv+=['--pedigree-column-dtype',f'{i}={dtype}']
            hashes={name:sha(source/name) for name in ['input.csv','family.csv']+(['pedigree.csv'] if case['pedigree_columns'] is not None else [])}
            process=subprocess.run(argv,capture_output=True,text=True)
            (args.output/(case['name']+'.stderr')).write_text(process.stderr)
            report=dict(name=case['name'],successful=False,returncode=process.returncode);reports.append(report)
            assert all(sha(source/name)==wanted for name,wanted in hashes.items())
            events=json.loads((source/'source-events.json').read_text())
            actual_proband=next((event['locals']['proband'] for event in reversed(events) if 'proband' in event['locals']),None)
            active=case['pedigree_columns'] is not None and case['family']['kind']!='none'
            reason=None
            if active:
                if case['unsupported_family_scalar_identity'] or case['unsupported_family_object_scalar_identity']:reason='NumPy family scalar/object identity'
                elif '#FamilyID' in case['pedigree_columns'] and case['pedigree_dtypes'][case['pedigree_columns'].index('#FamilyID')] not in ['object','string']:reason='native numeric family-column NumPy comparison'
                elif any(case['pedigree_columns'].count(name)>1 for name in ['#FamilyID','Phenotype','IndividualID','Sex']):reason='duplicate pedigree labels'
                elif actual_proband is not None and actual_proband['kind']!='str':reason='nontext selected proband'
                elif actual_proband is not None:
                    name=actual_proband['value']
                    if case['columns'].count(name)==1 and case['dtypes'][case['columns'].index(name)] in ['Int64','UInt64','Float64','boolean']:
                        reason='nullable numeric sample accessor/cache effects'
                    elif name not in case['columns'] and not re.fullmatch('[A-Za-z0-9_.-]*',name):reason='unported missing-sample repr'
            if reason:
                assert process.returncode and 'UnsupportedInput' in process.stderr and 'retain Python backend' in process.stderr,(case['name'],process.stderr)
                assert not result.exists() and not metadata.exists()
                report.update(successful=True,deferred_to_python=True,reason=reason,source_pickle_unchanged=case['input_pickle_unchanged']);continue
            source_warning='none'
            if case['source_log']:
                assert len(case['source_log'])==1
                source_warning='no_pedigree' if case['source_log'][0]['template'].startswith('No ped_table') else 'unknown_sex'
            if case['original_exception'] is not None:
                error=case['original_exception']
                assert process.returncode and error['type']+': '+error['message'] in process.stderr,(case['name'],process.stderr,error)
                assert f"source_failure_line={case['first_failure_line']}" in process.stderr,(case['name'],process.stderr)
                assert not result.exists() and not metadata.exists()
                if source_warning!='none':
                    assert f'source_warning={source_warning}' in process.stderr
                    assert process.stderr.index('source_warning=')<process.stderr.index(error['type']+': ')
                else:assert 'source_warning=' not in process.stderr
                assert case['input_pickle_unchanged'],'source object effect requires explicit Python fallback'
                report.update(successful=True,original_failure_preserved=True,source_failure_line=case['first_failure_line'],warning_descriptor_exact=True);continue
            assert process.returncode==0,(case['name'],process.stderr)
            with result.open(newline='') as stream:rows=list(csv.DictReader(stream))
            assert [int(row['source_position']) for row in rows]==list(range(case['rows']))
            with (source/'input.csv').open(newline='') as stream:source_rows=list(csv.reader(stream))[1:]
            assert [[row['index_kind'],row['index_value']] for row in rows]==[row[:2] for row in source_rows]
            with metadata.open(newline='') as stream:meta=dict((row['name'],row['value']) for row in csv.DictReader(stream))
            assert meta['axis']==axis and meta['source_warning']==source_warning
            data,pedigree,family=pickle.loads((source/'input.pkl').read_bytes())
            helper.encode((data,pedigree,family))
            before=pickle.dumps((data,pedigree,family),protocol=5)
            assert hashlib.sha256(before).hexdigest()==case['input_before_sha256']
            expected=pickle.loads((source/'expected.pkl').read_bytes())
            proband=meta['proband'] if meta['proband_present']=='1' else None
            assert helper.encode(proband)==helper.encode(expected['proband'])
            assert (meta['is_male']=='1')==expected['proband_is_male'] and (meta['is_female']=='1')==expected['proband_is_female']
            if meta['selected_pedigree_position']:
                # Source .loc[mask] takes a DataFrame before extracting a row;
                # preserve that take's index/scalar construction and metadata.
                row=int(meta['selected_pedigree_position']);selected=pedigree.iloc[[row]].iloc[0]
                assert helper.encode(selected)==helper.encode(expected['proband_row'])
                if source_warning=='unknown_sex':
                    raw=(selected.get('Sex',''),selected['IndividualID'])
                    assert helper.encode(raw)==case['source_log'][0]['args']
                    assert case['source_log'][0]['template']%raw==case['source_log'][0]['message']
            else:assert 'proband_row' not in expected
            values=np.array([row['value']=='1' for row in rows],dtype=bool)
            nulls=np.array([row['is_null']=='1' for row in rows],dtype=bool)
            dtype=meta['mask_dtype'];assert dtype==case['mask_dtype']
            if dtype=='boolean':array=pd.arrays.BooleanArray(values,nulls)
            else:
                assert not nulls.any()
                array=values.astype(object) if dtype=='object' else values
            mask=pd.Series(array,index=data.index,name=proband)
            if meta['inherit_frame_metadata']=='1':
                mask.attrs=copy.deepcopy(data.attrs);mask.flags.allows_duplicate_labels=data.flags.allows_duplicate_labels
            pd.testing.assert_series_equal(mask,expected['proband_het'],check_exact=True)
            assert helper.encode(mask)==helper.encode(expected['proband_het']),case['name']
            after=pickle.dumps((data,pedigree,family),protocol=5)
            assert before==after and case['input_pickle_unchanged'],case['name']
            report.update(successful=True,typed_context_and_mask_exact=True,raw_masked_bits_exact=True,source_index_metadata_exact=True,
                original_input_pickle_unchanged=True,warning_descriptor_and_original_args_exact=True,
                internal_normalized_sex_text_not_exported=True)
    except Exception as error:
        failures.append(dict(type=type(error).__name__,message=str(error),case=reports[-1]['name'] if reports else None));raise
    finally:
        final=dict(binary_sha256=args.binary_sha256,oracle_manifest_sha256=args.oracle_manifest_sha256,
            harness_sha256=sha(Path(__file__)),native_environment={key:os.environ.get(key) for key in ['SLURM_JOB_ID','SLURM_STEP_ID','SLURMD_NODENAME']},
            successful=not failures and len(reports)==len(oracle['cases']),cases=reports,failures=failures,
            actual_outputs=sum(case.get('typed_context_and_mask_exact',False) for case in reports),
            original_failures=sum(case.get('original_failure_preserved',False) for case in reports),
            deferred_to_python=sum(case.get('deferred_to_python',False) for case in reports),
            comparison='exact exported primitive context, mask data/null bits/index/dtype/attrs/flags, structured source warnings and supported original failure phase/type/message; internal normalized text stays source evidence',
            scientific_acceptance=False,production_backend_selected=False)
        (args.output/'comparison.json').write_text(json.dumps(final,indent=2)+'\n')
        print(json.dumps({key:value for key,value in final.items() if key!='cases'}),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    command=parser.add_subparsers(dest='command',required=True).add_parser('compare')
    for key in ['binary','oracle','output']:command.add_argument('--'+key,type=Path,required=True)
    for key in ['binary-sha256','oracle-manifest-sha256']:command.add_argument('--'+key,required=True)
    compare(parser.parse_args())


if __name__=='__main__':main()
