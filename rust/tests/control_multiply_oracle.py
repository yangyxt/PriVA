#!/usr/bin/env python3
"""Capture unchanged R3 multiplication/assignment before its optional Rust port.

capture --full-oracle CLOSED_ROOT --pedigree-oracle CLOSED_R2_ROOT
        --source-audit CLOSED_PROBES --output NEW_DATA_ROOT
Runs the unchanged L528 AST statement and binds actual full-ranker checkpoints.
For original fixtures, the preceding C1-C5 calls remain actual pinned Python/SciPy.
No scientific source, production dispatch, installation or family data is changed.
"""
import argparse
import ast
import copy
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pickle
import struct
import warnings

import numpy as np
import pandas as pd

D=Path('/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade')
FULL_PINS='9e63fb73e49c9f8ce700b0abd3e588702d2489c5430059e76c5d5b723eabb46a'
FULL_CAPTURE='db919b02f8738da60f85d16e0a7ed2bb53f8eaaf6047ed8f4db912fe4938508c'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path,value):
    with path.open('x') as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def original(path):
    text=path.read_text()
    function=next(node for node in ast.parse(text).body if isinstance(node,ast.FunctionDef) and node.name=='sort_and_rank_variants')
    body=[copy.deepcopy(node) for node in function.body if node.lineno==528]
    assert len(body)==1 and isinstance(body[0],ast.Assign)
    result=ast.Return(value=ast.Name(id='df',ctx=ast.Load()),lineno=529,col_offset=4,end_lineno=529,end_col_offset=13)
    wrapper=ast.FunctionDef(name='original_R3_multiply',args=ast.arguments(posonlyargs=[],args=[ast.arg(arg='df')],
        kwonlyargs=[],kw_defaults=[],defaults=[]),body=body+[result],decorator_list=[],lineno=528,col_offset=0,end_lineno=529,end_col_offset=13)
    scope={}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[wrapper],type_ignores=[])),str(path),'exec'),scope)
    return scope['original_R3_multiply'],dict(lines=[528,528],source_statement=ast.get_source_segment(text,body[0]),
        full_function_sha256=hashlib.sha256(ast.get_source_segment(text,function).encode()).hexdigest(),
        statement_ast_sha256=hashlib.sha256(ast.dump(ast.Module(body=body,type_ignores=[])).encode()).hexdigest())


def f64(bits):return struct.unpack('>d',struct.pack('>Q',bits))[0]


def examples():
    edges=[0.,-0.,1.,-1.,np.nextafter(1.,2.),np.finfo(float).max,np.nextafter(0.,1.),np.inf,-np.inf,
           f64(0x7ff8000000000123),f64(0xfff8000000000456),f64(0x7ff0000000000789),f64(0xfff0000000000abc)]
    left=[a for a in edges for _ in edges];right=edges*len(edges)
    for n in [0,1,17,33,len(left),1024]:
        frame=pd.DataFrame({'sort_index':pd.Series((left*8)[:n],dtype='float64'),
            'control_common_index':pd.Series((right*8)[:n],dtype='float64'),'untouched':['x']*n})
        yield 'float64_edges_'+str(n),frame
    for name,index in [('gapped',pd.Index([7,-2,100,3],name='gapped')),
        ('duplicate',pd.Index([7,7,-2,7],name='duplicate')),
        ('object',pd.Index(['x',None,'x',''],dtype=object,name='opaque')),
        ('multi',pd.MultiIndex.from_tuples([('x',1),('x',1),('a',2),('b',3)],names=['one','two']))]:
        frame=pd.DataFrame({'sort_index':[0.,-0.,np.inf,np.nan],'control_common_index':[.5,.5,0.,1.],'untouched':['x']*4})
        frame.index=index;frame.columns.name='column-axis';frame.attrs={'scope':'synthetic source audit'}
        yield 'axis_'+name,frame
    entries=[('float64',[1.,-0.,np.inf,np.nan]),('int64',[1,-1,0,2]),('uint64',[1,0,2,2**64-1]),
        ('bool',[True,False,True,False]),('Float64',[1.,None,-0.,np.inf]),('Int64',[1,None,-1,0]),
        ('UInt64',[1,None,0,2**64-1]),('boolean',[True,None,False,True]),
        ('object',[1.,None,'x',True]),('string',['1',None,'x','2'])]
    for role in ['sort_index','control_common_index']:
        for dtype,values in entries:
            frame=pd.DataFrame({'sort_index':[1.,2.,3.,4.],'control_common_index':[.5,.5,1.,1.],'untouched':['x']*4})
            frame[role]=pd.Series(values,dtype=dtype)
            yield role+'_dtype_'+dtype,frame
    for drop in [['sort_index'],['control_common_index'],['sort_index','control_common_index']]:
        yield 'missing_'+'_'.join(drop),pd.DataFrame({'sort_index':[1.],'control_common_index':[.5]}).drop(columns=drop)
    for name in ['sort_index','control_common_index','untouched']:
        frame=pd.DataFrame({'sort_index':[1.,2.],'control_common_index':[.5,1.],'untouched':['a','b']})
        yield 'duplicate_'+name,pd.concat([frame,frame[[name]]],axis=1)


def run(frame,function,helper):
    helper.encode(frame)
    before=pickle.dumps(frame,protocol=5)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:result=function(frame)
        except Exception as error:
            exception=dict(type=type(error).__name__,message=str(error))
            lines=[];trace=error.__traceback__
            while trace is not None:
                if trace.tb_frame.f_code is function.__code__:lines.append(trace.tb_lineno)
                trace=trace.tb_next
            exception['source_lines']=lines;result=None
        else:exception=None
    after=pickle.dumps(frame,protocol=5)
    return dict(exception=exception,returned_same_object=result is frame,
        input_before_sha256=hashlib.sha256(before).hexdigest(),input_after_sha256=hashlib.sha256(after).hexdigest(),
        input_pickle_unchanged=before==after,output=helper.encode(frame),
        warnings=[dict(category=item.category.__name__,message=str(item.message)) for item in caught]),before,after


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def verified(root,wanted):
    assert root.resolve().is_relative_to(D) and sha(root/'manifest.json')==wanted
    for name,digest in json.loads((root/'manifest.json').read_text()).items():
        path=root/name
        assert not path.is_symlink() and path.resolve().is_relative_to(root) and sha(path)==digest,name


def cell(value):
    if value is None:return ['null','']
    if value is pd.NA:return ['pd_na','']
    if isinstance(value,(bool,np.bool_)):return ['bool',str(int(value))]
    if isinstance(value,(int,np.integer)):return ['uint' if value>2**63-1 else 'int',str(value)]
    if isinstance(value,(float,np.floating)):return ['float',struct.pack('>d',value).hex()]
    if isinstance(value,str):return ['str',value.encode().hex()]
    raise TypeError('outside declared primitive transport: '+str(type(value)))


def export(output,name,frame,function,helper,origin,expected=None,expected_error=None):
    folder=output/name;folder.mkdir()
    prepared=pickle.dumps(frame,protocol=5)
    (folder/'prepared-input.pkl').write_bytes(prepared)
    result,before,after=run(pickle.loads(prepared),function,helper)
    repeated,_,repeated_after=run(pickle.loads(prepared),function,helper)
    assert result==repeated and after==repeated_after,name
    (folder/'input.pkl').write_bytes(before);(folder/'input-after.pkl').write_bytes(after)
    input_frame=pickle.loads(before)
    write(folder/'input.json',helper.encode(input_frame));write(folder/'result.json',result)
    with (folder/'input.csv').open('x',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['index_kind','index_value']+[item for label in input_frame.columns for item in [label+'_kind',label]])
        for i in range(len(input_frame)):
            writer.writerow(['opaque',helper.canonical(helper.encode(input_frame.index[i])).hex()]+[
                item for j in range(input_frame.shape[1]) for item in cell(input_frame.iloc[i,j])])
    report=dict(name=name,origin=origin,rows=len(input_frame),columns=list(input_frame.columns),
        dtypes=[str(value) for value in input_frame.dtypes],index_type=type(input_frame.index).__name__,
        original_exception=result['exception'],warnings=result['warnings'],
        input_before_sha256=result['input_before_sha256'],input_after_sha256=result['input_after_sha256'],
        source_input_pickle_unchanged=result['input_pickle_unchanged'],source_returned_same_object=result['returned_same_object'])
    if result['exception'] is None:
        (folder/'expected.pkl').write_bytes(after)
        if expected is not None:
            assert result['output']==expected,name
            report['actual_full_ranker_R4_entry_exact']=True
    elif expected_error is not None:
        assert result['exception']==expected_error,(name,result['exception'],expected_error)
        report['actual_full_ranker_R3_failure_exact']=True
    return report


def capture(args):
    native={key:os.environ.get(key) for key in ['SLURM_JOB_ID','SLURM_STEP_ID','SLURMD_NODENAME']}
    assert native['SLURM_JOB_ID'] and native['SLURM_STEP_ID']
    assert args.output.resolve().is_relative_to(D)
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'capture-harness.py').write_bytes(Path(__file__).read_bytes())
    write(args.output/'attempt.json',dict(native_environment=native,source_sha256=sha(Path(__file__))))
    assert sha(args.full_oracle/'source-pins.json')==FULL_PINS and sha(args.full_oracle/'capture.json')==FULL_CAPTURE
    helper=load('closed_full_ranking_oracle',args.full_oracle/'source/ranking_control_oracle.py')
    helper.verify_sources(args.full_oracle,FULL_PINS)
    pedigree_pin='88674165b68a09fca27d2b34fe70542cb1210278e7cd81471fa825edcafb6018'
    probe_pin='60abc4dcc10caf7eb9791ecb722781f2d1fdc673da41a8d2e930f0b76fbf16c8'
    verified(args.pedigree_oracle,pedigree_pin);verified(args.source_audit,probe_pin)
    function,source=original(args.full_oracle/'source/scripts/acmg_scoring.py')
    assert source['full_function_sha256']=='8aa29d9a4f83ff5629bbdf5e6279a43a374d71ca87fe04fe9db88253cc57d2ca'
    scoring,_,_,_=helper.original_modules(args.full_oracle/'source/scripts')
    write(args.output/'runtime-identities.json',helper.runtime_identity())
    full={case['folder']:case for case in json.loads((args.full_oracle/'capture.json').read_text())['cases']}
    reports,preceding_failures=[],[]
    pedigree=json.loads((args.pedigree_oracle/'oracle.json').read_text())
    for case in pedigree['cases']:
        if case['origin']['kind']!='closed-full-ranker-R3-locals-or-R2-failure' or case['original_exception'] is not None:continue
        actual=full[case['origin']['full_case']]
        trace=json.loads((args.full_oracle/actual['folder']/'traced/trace.json').read_text())['events']
        def state(event):return json.loads((args.full_oracle/'states'/(event['state']+'.json')).read_text())
        normal=[event for event in trace if event['function']=='compute_control_common' and event['event']=='return'
                and state(event).get('return_value',{}).get('kind')=='dataframe']
        frame,_,_=pickle.loads((args.pedigree_oracle/case['name']/'input.pkl').read_bytes())
        proband=pickle.loads((args.pedigree_oracle/case['name']/'expected.pkl').read_bytes())['proband']
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            try:control=scoring.compute_control_common(frame,proband)
            except Exception as error:
                exception=dict(type=type(error).__name__,message=str(error))
                assert not normal and actual['first_uncaught_failure']['phase'].startswith('C')
                assert exception==dict(type=actual['exception']['type'].rsplit('.',1)[-1],message=actual['exception']['message'])
                preceding_failures.append(dict(name=actual['name'],exception=exception,
                    full_case=actual['folder'],first_failure=actual['first_uncaught_failure'],not_a_Rust_boundary_case=True))
                continue
        assert len(normal)==1 and control is frame
        assert helper.encode(control)==state(normal[0])['return_value'],actual['name']
        following=[event for event in trace if event['function']=='sort_and_rank_variants' and event['phase']=='R4']
        expected=state(following[0])['df'] if following else None
        error=None
        if expected is None:
            assert actual['first_uncaught_failure']['phase']=='R3'
            error=dict(type=actual['exception']['type'].rsplit('.',1)[-1],message=actual['exception']['message'],source_lines=[528])
        reports.append(export(args.output,'full_'+actual['name'],control,function,helper,
            dict(kind='actual-normal-Python-C1-C5-return-to-R4-entry-or-R3-failure',full_case=actual['folder'],
                actual_preceding_control_return_state=normal[0]['state'],preceding_python_control_warnings=[
                    dict(category=item.category.__name__,message=str(item.message)) for item in caught]),expected,error))
    probes=json.loads((args.source_audit/'audit.json').read_text())
    for case in probes['cases']:
        frame=pickle.loads((args.source_audit/case['name']/'prepared-input.pkl').read_bytes())
        reports.append(export(args.output,'extra_'+case['name'],frame,function,helper,
            dict(kind='actual-verbatim-L528-source-probes',preceding_probe_manifest_sha256=probe_pin)))
    write(args.output/'oracle.json',dict(native_environment=native,source_commit=helper.COMMIT,source=source,
        full_oracle=str(args.full_oracle),full_source_pins_sha256=FULL_PINS,full_capture_sha256=FULL_CAPTURE,
        pedigree_oracle=str(args.pedigree_oracle),pedigree_manifest_sha256=pedigree_pin,
        source_audit=str(args.source_audit),source_audit_manifest_sha256=probe_pin,harness_sha256=sha(Path(__file__)),
        cases=reports,actual_python_control_failures_before_multiplication=preceding_failures,
        preceding_C1_C5_reexecuted_with_actual_frozen_Python_SciPy=True,Rust_not_implemented=True,
        scientific_acceptance=False,production_backend_selected=False))
    write(args.output/'manifest.json',{str(path.relative_to(args.output)):sha(path) for path in sorted(args.output.rglob('*')) if path.is_file()})
    print(json.dumps(dict(native_environment=native,cases=len(reports),source_outputs=sum(case['original_exception'] is None for case in reports),
        source_failures=sum(case['original_exception'] is not None for case in reports),preceding_Python_control_failures=len(preceding_failures),
        manifest_sha256=sha(args.output/'manifest.json'),harness_sha256=sha(Path(__file__)),Rust_not_implemented=True)),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    command=parser.add_subparsers(dest='command',required=True).add_parser('capture')
    for name in ['full-oracle','pedigree-oracle','source-audit','output']:command.add_argument('--'+name,type=Path,required=True)
    capture(parser.parse_args())


if __name__=='__main__':main()

