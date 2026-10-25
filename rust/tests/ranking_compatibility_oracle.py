#!/usr/bin/env python3
"""Capture original R4 L537-597 before its optional Rust implementation.

capture --full-oracle CLOSED_ROOT --output NEW_DATA_ROOT --source-log-level INFO
inspect --oracle CLOSED_ROOT
Executes unchanged source statements and the real mechanism helper at their
original position. Full-source entry/output checkpoints and isolated edge cases
retain ordered masks, logs, warnings, failure phase and immediate object state.
No production dispatch, family data, scientific source or dataset index changes.
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
import struct
import sys
import warnings

import numpy as np
import pandas as pd

D=Path('/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade')
FULL_PINS='9e63fb73e49c9f8ce700b0abd3e588702d2489c5430059e76c5d5b723eabb46a'
FULL_CAPTURE='db919b02f8738da60f85d16e0a7ed2bb53f8eaaf6047ed8f4db912fe4938508c'
FUNCTION_SHA='8aa29d9a4f83ff5629bbdf5e6279a43a374d71ca87fe04fe9db88253cc57d2ca'
INPUTS=['df','proband_is_male','proband_is_female','proband_het']
LOCALS=['mechanism_masks','chrom','x_locus','has_single_alt_lof_compatible',
    'exact_gof','asserted_lof_effect','dn_without_hi_or_ambiguous',
    'heterozygous_lof_incompatible','high_acmg_lof','compatibility_penalty_rows','pext_modulation_rows']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path,value):
    with path.open('x') as stream:json.dump(value,stream,indent=2,ensure_ascii=True,allow_nan=False);stream.write('\n')


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def excerpt(path,scoring):
    text=path.read_text()
    source=next(node for node in ast.parse(text).body if isinstance(node,ast.FunctionDef) and node.name=='sort_and_rank_variants')
    assert hashlib.sha256(ast.get_source_segment(text,source).encode()).hexdigest()==FUNCTION_SHA
    body=[copy.deepcopy(node) for node in source.body if 537<=node.lineno<=597]
    assert body[0].lineno==537 and body[-1].end_lineno==597
    result=ast.Return(value=ast.Call(func=ast.Name(id='locals',ctx=ast.Load()),args=[],keywords=[]),
        lineno=598,col_offset=4,end_lineno=598,end_col_offset=19)
    function=ast.FunctionDef(name='original_R4_excerpt',args=ast.arguments(posonlyargs=[],
        args=[ast.arg(arg=name) for name in INPUTS],kwonlyargs=[],kw_defaults=[],defaults=[]),
        body=body+[result],decorator_list=[],lineno=537,col_offset=0,end_lineno=598,end_col_offset=19)
    namespace=dict(scoring.__dict__)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[])),str(path),'exec'),namespace)
    return namespace['original_R4_excerpt'],dict(lines=[537,597],full_function_sha256=FUNCTION_SHA,
        statement_ast_sha256=hashlib.sha256(ast.dump(ast.Module(body=body,type_ignores=[])).encode()).hexdigest(),
        mechanism_helper_module=scoring._variant_mechanism_masks.__module__,
        mechanism_helper_bytecode_sha256=hashlib.sha256(scoring._variant_mechanism_masks.__code__.co_code).hexdigest())


class Logs(logging.Handler):
    def __init__(self,helper):super().__init__();self.helper=helper;self.records=[]
    def emit(self,record):
        self.records.append(dict(level=record.levelname,message=record.getMessage(),
            template=record.msg,args=self.helper.encode(record.args)))


def selected(values):return {name:values[name] for name in INPUTS+LOCALS if name in values}


def aliases(values):
    pairs=[]
    for left,right in [('pext_modulation_rows','high_acmg_lof')]:
        if left in values and right in values:pairs.append(dict(left=left,right=right,same_object=values[left] is values[right]))
    for left,right in [('exact_gof','is_exact_gof'),('asserted_lof_effect','is_predicted_lof')]:
        if left in values and 'mechanism_masks' in values:
            pairs.append(dict(left=left,right='mechanism_masks.'+right,same_object=values[left] is values['mechanism_masks'][right]))
    return pairs


def run_current(inputs,original,scoring,helper,entry_observation=None):
    """Run only original R4 on these buffers; never rehydrate a fallback input."""
    entry=helper.encode(inputs) if entry_observation is None else entry_observation
    before=pickle.dumps(inputs,protocol=5)
    events=[]
    def trace(frame,event,arg):
        if frame.f_code is original.__code__ and event in ['exception','return']:
            events.append(dict(event=event,line=frame.f_lineno,locals=selected(frame.f_locals)))
        return trace
    observed=Logs(helper);scoring.logger.addHandler(observed)
    previous=sys.gettrace()
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            try:
                sys.settrace(trace);result=original(*inputs)
            except Exception as error:
                exception=dict(type=type(error).__name__,message=str(error));result=None
            else:exception=None
            finally:sys.settrace(previous)
    finally:scoring.logger.removeHandler(observed)
    # Match the producing full-oracle convention before post-observation.
    after=pickle.dumps(inputs,protocol=5)
    output=selected(result) if result is not None else None
    output_pickle=pickle.dumps(output,protocol=5) if output is not None else None
    info=dict(exception=exception,input_before_sha256=hashlib.sha256(before).hexdigest(),
        input_after_sha256=hashlib.sha256(after).hexdigest(),input_pickle_unchanged=before==after,
        input_before=entry,input_after=helper.encode(inputs),output=helper.encode(output),
        output_aliases=aliases(output or {}),source_logs=observed.records,
        warnings=[dict(category=item.category.__name__,message=str(item.message),source_line=item.lineno) for item in caught],
        source_events=[dict(event=item['event'],line=item['line'],locals=helper.encode(item['locals'])) for item in events],
        first_failure_line=next((item['line'] for item in events if item['event']=='exception'),None),
        returned_same_frame=output is not None and output['df'] is inputs[0])
    return info,before,after,output_pickle


def run(prepared,original,scoring,helper):
    return run_current(pickle.loads(prepared),original,scoring,helper)


def cell(value):
    if value is None:return ['null','']
    if value is pd.NA:return ['pd_na','']
    if isinstance(value,(bool,np.bool_)):return ['bool',str(int(value))]
    if isinstance(value,(int,np.integer)):return ['uint' if value>2**63-1 else 'int',str(value)]
    if isinstance(value,(float,np.floating)):return ['float',struct.pack('>d',value).hex()]
    if isinstance(value,str):return ['str',value.encode().hex()]
    return ['unsupported',type(value).__module__+'.'+type(value).__qualname__]


def export(root,name,inputs,original,scoring,helper,origin,expected=None,expected_error=None):
    folder=root/name;folder.mkdir()
    prepared=pickle.dumps(inputs,protocol=5);(folder/'prepared-input.pkl').write_bytes(prepared)
    result,before,after,output=run(prepared,original,scoring,helper)
    repeated,rb,ra,ro=run(prepared,original,scoring,helper)
    assert result==repeated and (before,after,output)==(rb,ra,ro),name
    (folder/'input.pkl').write_bytes(before);(folder/'input-after.pkl').write_bytes(after)
    write(folder/'result.json',result)
    if output is not None:(folder/'expected.pkl').write_bytes(output)
    frame,male,female,het=pickle.loads(before)
    with (folder/'input.csv').open('x',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['index_kind','index_value']+[item for label in frame.columns for item in [str(label)+'_kind',str(label)]])
        for row in range(len(frame)):
            writer.writerow(['opaque',helper.canonical(helper.encode(frame.index[row])).hex()]+[
                item for col in range(frame.shape[1]) for item in cell(frame.iloc[row,col])])
    with (folder/'heterozygosity.csv').open('x',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['index_kind','index_value','value_kind','value','raw_data','is_null'])
        for row,value in enumerate(het):
            array=het.array
            raw=bool(array._data[row]) if hasattr(array,'_mask') else None
            missing=bool(array._mask[row]) if hasattr(array,'_mask') else False
            writer.writerow(['opaque',helper.canonical(helper.encode(het.index[row])).hex(),*cell(value),
                '' if raw is None else int(raw),int(missing)])
    report=dict(name=name,origin=origin,rows=len(frame),columns=list(frame.columns),dtypes=[str(value) for value in frame.dtypes],
        index_kind='multi' if isinstance(frame.index,pd.MultiIndex) else 'flat',heterozygosity_dtype=str(het.dtype),
        heterozygosity_index_equals_frame=het.index.equals(frame.index),heterozygosity_name=helper.encode(het.name),
        proband_is_male=helper.encode(male),proband_is_female=helper.encode(female),
        original_exception=result['exception'],first_failure_line=result['first_failure_line'],
        input_before_sha256=result['input_before_sha256'],input_after_sha256=result['input_after_sha256'],
        source_input_pickle_unchanged=result['input_pickle_unchanged'],warnings=result['warnings'],source_logs=result['source_logs'])
    if expected is not None:
        assert result['exception'] is None
        actual=pickle.loads(output)
        assert {key:helper.encode(actual[key]) for key in LOCALS+['df']}=={
            key:expected[key] for key in LOCALS+['df']},name
        report['full_source_R5_entry_exact']=True
    if expected_error is not None:
        assert result['exception']==expected_error,(name,result['exception'],expected_error)
        report['full_source_R4_failure_exact']=True
    return report


def base(n=4):
    return pd.DataFrame({'var_plausible_patho_mechs':['','dominant_LOF','x_linked_recessive_LOF','dominant_DN'][:n],
        'variant_effect':['predicted_LOF_high_confidence']*n,'variant_lof_score':[1]*n,
        'ACMG_quant_score':[.95]*n,'sort_index':[.475]*n,'chrom':['X']*n,'untouched':['u']*n})


def extras():
    profiles=['','dominant','dominant_LOF','dominant_GOF','dominant_DN','recessive_LOF',
        'x_linked_recessive_LOF','x_linked_recessive','x_linked_dominant_LOF','x_linked_dominant',
        'x_linked_dominant_DN','x_linked_unspecified_LOF',' dominant_LOF','dominant_LOF\n','dominant_LOF\n\n',
        'dominant_DN;dominant_LOF','x_linked_dominant_DN;dominant','dominant_LOF_extra']
    for male,female,label in [(False,False,'unknown'),(True,False,'male'),(False,True,'female'),(True,True,'both')]:
        frame=pd.concat([base(1)]*len(profiles),ignore_index=True);frame['var_plausible_patho_mechs']=profiles
        yield 'histories_'+label,(frame,male,female,pd.Series([True]*len(frame),index=frame.index,name='PROBAND'))
    chroms=['X','x','chrX','CHRX','cHrX','chrchrX','CHR X',' X','X ','chrX\n','X\n','XX','23','chr23','Χ','ⅹ','ｘ','',None,pd.NA,np.nan,1,True]
    for dtype,values in [('object',chroms),('string',chroms[:-2]),('int64',[1,23,0,-1]),
        ('float64',[1.,23.,np.nan,np.inf]),('Int64',[1,None,23,0]),('Float64',[1.,None,23.,0.]),
        ('boolean',[True,None,False,True]),('category',['X',None,'chrX','1'])]:
        frame=pd.concat([base(1)]*len(values),ignore_index=True);frame['chrom']=pd.Series(values,dtype=dtype)
        yield 'chrom_dtype_'+dtype,(frame,True,False,pd.Series([True]*len(frame),name='PROBAND'))
    edges=[np.nextafter(.89,-np.inf),.89,np.nextafter(.89,np.inf),0.,-0.,np.inf,-np.inf,np.nan]
    for dtype,values in [('float64',edges),('float32',edges),('int64',[-1,0,1,2]),
        ('uint64',[0,1,2,2**64-1]),('bool',[True,False,True,False]),('Float64',[.95,None,-0.,np.inf]),
        ('Int64',[1,None,-1,0]),('UInt64',[1,None,0,2**64-1]),('boolean',[True,None,False,True]),
        ('object',[.95,None,'0.95',True]),('string',['.95',None,'1','0'])]:
        frame=pd.concat([base(1)]*len(values),ignore_index=True);frame['ACMG_quant_score']=pd.Series(values,dtype=dtype)
        yield 'score_dtype_'+dtype,(frame,False,False,pd.Series([True]*len(frame),name='PROBAND'))
    for dtype,values in [('bool',[True,False,True,False]),('boolean',[True,None,False,None]),
        ('object',[True,False,None,pd.NA]),('int64',[1,0,1,0]),('string',['1',None,'0','True'])]:
        frame=base();het=pd.Series(values,dtype=dtype,name='PROBAND')
        yield 'het_dtype_'+dtype,(frame,True,False,het)
    for bits in [[False,False,False,False],[True,True,True,True]]:
        frame=base();het=pd.Series(pd.arrays.BooleanArray(np.array(bits),np.array([False,True,False,True])),name='PROBAND')
        yield 'het_masked_raw_'+str(int(bits[0])),(frame,True,False,het)
    for name,index in [('duplicate',pd.Index([7,7,-2,7],name='duplicate')),
        ('gapped',pd.Index([8,-2,90,4],name='gapped')),
        ('multi',pd.MultiIndex.from_tuples([('a',1),('a',1),('b',2),('a',3)],names=['left','right']))]:
        frame=base();frame.index=index;frame.attrs={'scope':'R4 fixture'};frame.columns.name='column-axis'
        het=pd.Series([True,False,True,True],index=index,name='PROBAND');het.attrs={'context':'R2 mask'}
        yield 'axis_'+name,(frame,True,False,het)
    for index,label in [(pd.Index([3,2,1,0]),'reordered'),(pd.Index([0,1,2,7]),'missing_extra'),
        (pd.Index([0,0,1,2]),'different_duplicates')]:
        yield 'unaligned_het_'+label,(base(),True,False,pd.Series([True]*4,index=index,name='PROBAND'))
    for n in [0,1,4]:
        for dtype in ['bool','boolean','object']:
            frame=base(n);het=pd.Series([True]*n,dtype=dtype,index=frame.index,name='PROBAND')
            if not n:het.name=None
            yield 'size_'+str(n)+'_het_'+dtype,(frame,False,False,het)
    for name in ['var_plausible_patho_mechs','variant_effect','ACMG_quant_score','chrom']:
        yield 'missing_'+name,(base().drop(columns=[name]),True,False,pd.Series([True]*4,name='PROBAND'))
    yield 'missing_both_mechanism_and_score',(base().drop(columns=['var_plausible_patho_mechs','variant_effect','ACMG_quant_score']),True,False,pd.Series([True]*4))
    for name in ['var_plausible_patho_mechs','variant_effect','variant_lof_score','ACMG_quant_score','chrom','untouched']:
        frame=base();frame=pd.concat([frame,frame[[name]]],axis=1)
        yield 'duplicate_'+name,(frame,True,False,pd.Series([True]*4,name='PROBAND'))
    target='zygosity_inheritance_mechanism_compatibility'
    for dtype,values in [('float64',[7.]*4),('int64',[7]*4),('string',['old']*4),('Float64',[7.,None,7.,7.])]:
        frame=base();frame.insert(1,target,pd.Series(values,dtype=dtype))
        yield 'existing_compatibility_'+dtype,(frame,False,False,pd.Series([True]*4,name='PROBAND'))
    frame=base();frame.insert(1,target,7.);frame=pd.concat([frame,frame[[target]]],axis=1)
    yield 'duplicate_compatibility',(frame,False,False,pd.Series([True]*4,name='PROBAND'))
    for name,dtype in [('var_plausible_patho_mechs','Int64'),('variant_effect','boolean'),('variant_lof_score','boolean')]:
        frame=base();frame[name]=pd.Series([None,1,0,1],dtype=dtype)
        yield 'mechanism_nullable_'+name,(frame,True,False,pd.Series([True]*4,name='PROBAND'))
    for sample_dtype in ['object','string']:
        frame=base();sample=pd.Series(['0/1','0/1:10:20',None,'0/0:1'],dtype=sample_dtype,name='PROBAND')
        yield 'actual_R2_het_'+sample_dtype,(frame,True,False,sample.str.count('1')==1)
    frame=base();frame.attrs={'scope':'attrs and flags'};frame.flags.allows_duplicate_labels=False
    het=pd.Series([True]*4,name='PROBAND');het.attrs={'scope':'R2 attrs'};het.flags.allows_duplicate_labels=False
    yield 'attrs_flags',(frame,True,False,het)


def capture(args):
    native={key:os.environ.get(key) for key in ['SLURM_JOB_ID','SLURM_STEP_ID','SLURMD_NODENAME']}
    assert native['SLURM_JOB_ID'] and native['SLURM_STEP_ID']
    assert args.full_oracle.resolve().is_relative_to(D) and args.output.resolve().is_relative_to(D)
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'capture-harness.py').write_bytes(Path(__file__).read_bytes())
    write(args.output/'attempt.json',dict(native_environment=native,harness_sha256=sha(Path(__file__))))
    assert sha(args.full_oracle/'source-pins.json')==FULL_PINS and sha(args.full_oracle/'capture.json')==FULL_CAPTURE
    helper=load('frozen_full_R4_oracle',args.full_oracle/'source/ranking_control_oracle.py')
    helper.verify_sources(args.full_oracle,FULL_PINS)
    scoring,_,_,_=helper.original_modules(args.full_oracle/'source/scripts')
    logging_configuration=dict(original_logger_level=scoring.logger.level,
        original_effective_level=scoring.logger.getEffectiveLevel(),
        captured_source_level=args.source_log_level,
        scope='explicit fixture caller logging configuration; scientific source and arguments unchanged')
    scoring.logger.setLevel(getattr(logging,args.source_log_level))
    original,source=excerpt(args.full_oracle/'source/scripts/acmg_scoring.py',scoring)
    write(args.output/'runtime-identities.json',helper.runtime_identity())
    write(args.output/'chromosome-uppercase-witness.json',dict(
        nonascii_upper_entirely_CHRX=[point for point in range(128,0x110000)
            if chr(point).upper() and set(chr(point).upper())<=set('CHRX')],
        scope='ASCII X and optional leading CHR predicate only, not general Unicode uppercase parity'))
    captured=json.loads((args.full_oracle/'capture.json').read_text())
    reports=[];preceding=[]
    for case in captured['cases']:
        if case['kind']!='ranking':continue
        folder=args.full_oracle/case['folder']
        events=json.loads((folder/'traced/trace.json').read_text())['events']
        def state(event):return json.loads((args.full_oracle/'states'/(event['state']+'.json')).read_text())
        entry=next((event for event in events if event['function']=='sort_and_rank_variants' and event['phase']=='R4'),None)
        if entry is None:
            preceding.append(dict(name=case['name'],first_failure=case['first_uncaught_failure']));continue
        # The closed observer deduplicates unchanged states across phase labels.
        # A skipped R5 can therefore have no explicit R5 event. Its first later
        # state still carries the R4 locals; prove these against a fresh real
        # L599 boundary below rather than infer completion from a phase label.
        following=next((event for event in events if event['function']=='sort_and_rank_variants'
                        and event['line']>597),None)
        expected=state(following) if following else None
        expected_error=None
        if following is None:
            assert case['first_uncaught_failure']['phase']=='R4'
            expected_error=dict(type=case['exception']['type'].rsplit('.',1)[-1],message=case['exception']['message'])
        # Replay the actual full function to recover genuine entry buffers and
        # their cross-object aliases; never synthesize a frame from JSON values.
        data,pedigree=pickle.loads((folder/'prepared-input.pkl').read_bytes())
        helper.encode(data);helper.encode(pedigree)
        bound=json.loads((folder/'traced/arguments.json').read_text())['arguments']
        kwargs={key['value']:helper.decode_atom(value) for key,value in bound['items']
                if key['value'] not in {'df','ped_df'}}
        entries=[];boundaries=[]
        def trace(frame,event,arg):
            if frame.f_code is scoring.sort_and_rank_variants.__code__ and event=='line' and frame.f_lineno==537:
                values=tuple(frame.f_locals[name] for name in INPUTS)
                entries.append(pickle.dumps(values,protocol=5))
                assert {name:helper.encode(value) for name,value in zip(INPUTS,values)}=={
                    name:state(entry)[name] for name in INPUTS},case['name']
            if frame.f_code is scoring.sort_and_rank_variants.__code__ and event=='line' and frame.f_lineno==599:
                boundaries.append(pickle.dumps(selected(frame.f_locals),protocol=5))
            return trace
        previous=sys.gettrace()
        with warnings.catch_warnings(record=True):
            warnings.simplefilter('always')
            try:
                sys.settrace(trace);full_output=scoring.sort_and_rank_variants(data,pedigree,**kwargs)
            except Exception as error:
                assert case['status']=='error'
                assert dict(type=type(error).__name__,message=str(error))==dict(
                    type=case['exception']['type'].rsplit('.',1)[-1],message=case['exception']['message'])
            else:
                assert case['status']=='value'
                assert helper.encode(full_output)==json.loads((folder/'traced/outcome.json').read_text())['output'],case['name']
            finally:sys.settrace(previous)
        assert len(entries)==1,case['name']
        if expected is not None:
            assert len(boundaries)==1
            boundary=pickle.loads(boundaries[0])
            assert {key:helper.encode(boundary[key]) for key in LOCALS+['df']}=={
                key:expected[key] for key in LOCALS+['df']},case['name']
        else:assert not boundaries
        reports.append(export(args.output,'full_'+case['name'],pickle.loads(entries[0]),original,scoring,helper,
            dict(kind='actual-full-source-R4-entry',full_case=case['folder'],scientific_arguments_sha256=sha(folder/'traced/arguments.json')),
            expected,expected_error))
        if boundaries:(args.output/('full_'+case['name'])/'actual-full-L599-boundary.pkl').write_bytes(boundaries[0])
    for name,inputs in extras():
        reports.append(export(args.output,'extra_'+name,inputs,original,scoring,helper,dict(kind='unchanged-original-R4-statements')))
    write(args.output/'oracle.json',dict(source_commit=helper.COMMIT,source=source,harness_sha256=sha(Path(__file__)),
        native_environment=native,full_oracle=str(args.full_oracle),full_source_pins_sha256=FULL_PINS,
        full_capture_sha256=FULL_CAPTURE,cases=reports,preceding_full_source_failures=preceding,
        source_logging_configuration=logging_configuration,
        scope='R4 compatibility masks/logs/column assignment only; R6 multiplication remains later; no scientific or production activation',
        scientific_acceptance=False,production_backend_selected=False))
    write(args.output/'manifest.json',{str(path.relative_to(args.output)):sha(path) for path in sorted(args.output.rglob('*')) if path.is_file()})
    print(json.dumps(dict(cases=len(reports),outputs=sum(case['original_exception'] is None for case in reports),
        source_failures=sum(case['original_exception'] is not None for case in reports),
        preceding_failures=len(preceding),manifest_sha256=sha(args.output/'manifest.json'),native_environment=native)),flush=True)


def inspect(args):
    assert args.oracle.resolve().is_relative_to(D)
    data=json.loads((args.oracle/'oracle.json').read_text())
    for name,wanted in json.loads((args.oracle/'manifest.json').read_text()).items():assert sha(args.oracle/name)==wanted
    print(json.dumps(dict(cases=len(data['cases']),full_source_cases=sum(case['origin']['kind']=='actual-full-source-R4-entry' for case in data['cases']),
        outputs=sum(case['original_exception'] is None for case in data['cases']),
        source_failures=[dict(name=case['name'],exception=case['original_exception'],line=case['first_failure_line'])
            for case in data['cases'] if case['original_exception'] is not None],
        manifest_sha256=sha(args.oracle/'manifest.json'))),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    command=commands.add_parser('capture')
    for name in ['full-oracle','output']:command.add_argument('--'+name,type=Path,required=True)
    command.add_argument('--source-log-level',choices=['INFO','WARNING'],default='INFO')
    commands.add_parser('inspect').add_argument('--oracle',type=Path,required=True)
    args=parser.parse_args();(capture if args.command=='capture' else inspect)(args)


if __name__=='__main__':main()
