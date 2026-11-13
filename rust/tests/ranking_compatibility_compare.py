#!/usr/bin/env python3
"""Compare optional R4 masks/metadata/mutation and stage-only Python fallback.

compare --binary FILE --binary-sha256 SHA --oracle ROOT
        --oracle-manifest-sha256 SHA --output NEW_ROOT
No prior C1-C5/R3, R2 warnings or R6 multiplication is replayed at this boundary.
"""
import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import pickle
import struct
import subprocess
import warnings

import numpy as np
import pandas as pd


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def reason(case,source_rows,het_rows):
    if not (case['heterozygosity_dtype'] in {'bool','boolean'}
            or case['heterozygosity_dtype']=='object' and not case['rows']):return 'R2 mask representation'
    if [row[:2] for row in source_rows]!=[row[:2] for row in het_rows]:return 'unaligned R2 mask'
    used=['var_plausible_patho_mechs','variant_effect','variant_lof_score','variant_gof_score','variant_dn_score']
    native={'object','string','float64','int64','uint64','bool'}
    if any(dtype not in native for label,dtype in zip(case['columns'],case['dtypes']) if label in used):return 'mechanism dtype/cache representation'
    if case['first_failure_line']==537:return None # Existing helper errors precede later guards.
    if any(case['columns'].count(name)>1 for name in ['chrom','ACMG_quant_score','zygosity_inheritance_mechanism_compatibility']):return 'ambiguous R4 label'
    if 'chrom' in case['columns']:
        pos=case['columns'].index('chrom');dtype=case['dtypes'][pos]
        if dtype not in {'object','string','int64','uint64','bool'}:return 'chromosome dtype'
        for row in source_rows:
            kind,value=row[2+2*pos:4+2*pos]
            if dtype in {'object','string'} and kind in {'str','null','pd_na'}:continue
            if dtype=='object' and kind=='float' and math.isnan(struct.unpack('>d',bytes.fromhex(value))[0]):continue
            if dtype=='int64' and kind=='int':continue
            if dtype=='uint64' and kind in {'uint','int'} and int(value)>=0:continue
            if dtype=='bool' and kind=='bool':continue
            return 'chromosome cell representation'
    if 'ACMG_quant_score' in case['columns']:
        pos=case['columns'].index('ACMG_quant_score');dtype=case['dtypes'][pos]
        if dtype not in {'float64','int64','uint64','bool'}:return 'raw score dtype/promotion'
    return None


def compare(args):
    assert sha(args.binary)==args.binary_sha256 and sha(args.oracle/'manifest.json')==args.oracle_manifest_sha256
    for name,wanted in json.loads((args.oracle/'manifest.json').read_text()).items():assert sha(args.oracle/name)==wanted
    oracle=json.loads((args.oracle/'oracle.json').read_text())
    harness=load('frozen_R4_oracle',args.oracle/'capture-harness.py')
    helper=load('full_ranker_R4_helper',Path(oracle['full_oracle'])/'source/ranking_control_oracle.py')
    scoring,_,_,_=helper.original_modules(Path(oracle['full_oracle'])/'source/scripts')
    scoring.logger.setLevel(oracle['source_logging_configuration']['captured_source_level'])
    original,_=harness.excerpt(Path(oracle['full_oracle'])/'source/scripts/acmg_scoring.py',scoring)
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'comparison-harness.py').write_bytes(Path(__file__).read_bytes())
    reports,failures=[],[]
    try:
        for case in oracle['cases']:
            source=args.oracle/case['name'];output=args.output/(case['name']+'.csv');metadata=args.output/(case['name']+'.metadata.csv')
            # Reconstruct the producing prepared state; snapshot rehydration can
            # sever hidden cache aliases. Never normalize/reset source caches.
            inputs=pickle.loads((source/'prepared-input.pkl').read_bytes())
            entry=helper.encode(inputs);before=pickle.dumps(inputs,protocol=5)
            assert hashlib.sha256(before).hexdigest()==case['input_before_sha256'],case['name']
            with (source/'input.csv').open(newline='') as stream:source_rows=list(csv.reader(stream))[1:]
            with (source/'heterozygosity.csv').open(newline='') as stream:het_rows=list(csv.reader(stream))[1:]
            why=reason(case,source_rows,het_rows)
            frame,male,female,het=inputs;axis=sha(source/'input.pkl')
            argv=[str(args.binary),'ranking-compatibility','--input',str(source/'input.csv'),
                '--heterozygosity-input',str(source/'heterozygosity.csv'),'--heterozygosity-dtype',case['heterozygosity_dtype'],
                '--proband-is-male',str(male).lower(),'--proband-is-female',str(female).lower(),
                '--output',str(output),'--metadata-output',str(metadata),'--axis-metadata',axis]
            for i,dtype in enumerate(case['dtypes']):argv+=['--column-dtype',f'{i}={dtype}']
            input_shas=[sha(source/name) for name in ['input.csv','heterozygosity.csv']]
            process=subprocess.run(argv,capture_output=True,text=True)
            (args.output/(case['name']+'.stderr')).write_text(process.stderr)
            report=dict(name=case['name'],successful=False,returncode=process.returncode);reports.append(report)
            assert [sha(source/name) for name in ['input.csv','heterozygosity.csv']]==input_shas
            assert pickle.dumps(inputs,protocol=5)==before
            if why:
                assert process.returncode and 'UnsupportedInput' in process.stderr and 'retain Python backend' in process.stderr,(case['name'],why,process.stderr)
                assert not output.exists() and not metadata.exists()
                same_frame=id(frame)
                result,fb,fa,_=harness.run_current(inputs,original,scoring,helper,entry_observation=entry)
                assert id(inputs[0])==same_frame
                assert result==json.loads((source/'result.json').read_text()),case['name']
                assert fb==before and hashlib.sha256(fa).hexdigest()==case['input_after_sha256']
                report.update(successful=True,deferred_to_python=True,reason=why,
                    original_current_buffer_stage_only_exact=True,earlier_stages_not_replayed=True,
                    fallback_source_failed=result['exception'] is not None,preflight_current_buffers_unchanged=True)
                continue
            if case['original_exception'] is not None:
                error=case['original_exception']
                assert process.returncode and error['type']+': '+error['message'] in process.stderr,(case['name'],process.stderr,error)
                assert 'source_failure_line='+str(case['first_failure_line']) in process.stderr
                assert not output.exists() and not metadata.exists()
                assert case['source_input_pickle_unchanged'] and pickle.dumps(inputs,protocol=5)==before,case['name']
                report.update(successful=True,original_failure_preserved=True,source_failure_line=case['first_failure_line'])
                continue
            assert process.returncode==0,(case['name'],process.stderr)
            with output.open(newline='') as stream:
                reader=csv.DictReader(stream);headers=reader.fieldnames;rows=list(reader)
            assert [int(row['source_position']) for row in rows]==list(range(case['rows']))
            assert [[row['index_kind'],row['index_value']] for row in rows]==[row[:2] for row in source_rows]
            with metadata.open(newline='') as stream:meta={row['name']:row['value'] for row in csv.DictReader(stream)}
            assert meta['axis']==axis
            def series(name,values,dtype):
                origin=meta[name+'.metadata_origin']
                if origin=='fresh':attrs={};duplicate=True
                else:
                    assert origin.startswith('column:')
                    original_column=frame.iloc[:,int(origin.split(':')[1])]
                    attrs=copy.deepcopy(original_column.attrs);duplicate=original_column.flags.allows_duplicate_labels
                label=None if meta[name+'.series_name_kind']=='none' else meta[name+'.series_name']
                result=pd.Series(values,index=frame.index,name=label,dtype=dtype)
                result.attrs=attrs;result.flags.allows_duplicate_labels=duplicate
                return result
            masks={name:series(name,[bool(int(row[name])) for row in rows],'bool') for name in headers[5:32]}
            assert len(masks)==27
            actual=dict(zip(harness.INPUTS,inputs));actual['mechanism_masks']=masks
            actual['chrom']=series('chrom',[bytes.fromhex(row['chrom_hex']).decode() for row in rows],'object')
            for name in ['x_locus','has_single_alt_lof_compatible','dn_without_hi_or_ambiguous',
                         'heterozygous_lof_incompatible','high_acmg_lof']:
                actual[name]=series(name,[bool(int(row[name])) for row in rows],'bool')
            actual['exact_gof']=masks['is_exact_gof'];actual['asserted_lof_effect']=masks['is_predicted_lof']
            data=np.array([bool(int(row['compatibility_penalty_rows'])) for row in rows],dtype=bool)
            null=np.array([bool(int(row['compatibility_penalty_rows_null'])) for row in rows],dtype=bool)
            penalty_values=pd.arrays.BooleanArray(data,null) if meta['penalty_dtype']=='boolean' else data
            assert meta['penalty_dtype']=='boolean' or not null.any()
            actual['compatibility_penalty_rows']=series('compatibility_penalty_rows',penalty_values,meta['penalty_dtype'])
            actual['pext_modulation_rows']=actual['high_acmg_lof']
            expected_record=json.loads((source/'result.json').read_text())
            log=dict(level='INFO',message=meta['log_message'],template=meta['log_message'],args=helper.encode(()))
            assert [log]==case['source_logs'],case['name']
            target='zygosity_inheritance_mechanism_compatibility'
            position=case['columns'].index(target) if target in case['columns'] else len(case['columns'])
            assert int(meta['target_source_column_position'])==position
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                frame[target]=1.0
                frame.loc[actual['compatibility_penalty_rows'],target]=0.9
            after=pickle.dumps(inputs,protocol=5)
            assert hashlib.sha256(after).hexdigest()==case['input_after_sha256'],case['name']
            assert [dict(category=item.category.__name__,message=str(item.message)) for item in caught]==[
                {key:item[key] for key in ['category','message']} for item in case['warnings']],case['name']
            assert [struct.pack('>d',value).hex() for value in frame.iloc[:,position]]==[row['compatibility_bits'] for row in rows]
            actual={name:actual[name] for name in harness.INPUTS+harness.LOCALS}
            assert helper.encode(actual)==expected_record['output'],case['name']
            assert harness.aliases(actual)==expected_record['output_aliases']
            assert helper.encode(inputs)==expected_record['input_after']
            report.update(successful=True,full_typed_output_exact=True,raw_boolean_null_bits_and_float_bits_exact=True,
                opaque_metadata_and_aliases_exact=True,source_input_mutation_and_pickle_exact=True,
                source_logs_warnings_column_position_exact=True,R6_multiplication_not_run=True,
                preflight_current_buffers_unchanged=True)
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
            comparison='R4 only; exact masks/raw null data/metadata/aliases/INFO/warnings/column position/current-buffer mutation and fallback; R6 unchanged',
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
