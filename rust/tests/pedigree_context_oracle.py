#!/usr/bin/env python3
"""Capture actual original R2 locals before its independent optional Rust port.

capture --full-oracle CLOSED_FULL_ROOT --setup-oracle CLOSED_R1_ROOT
        --output NEW_DATA_ROOT
Runs the one original R2 if-statement, retaining source locations, pandas types,
warnings, partial failure locals and immediate pre/post input serialization.
No C1-C5, numerical scoring, family data or production dispatch is called.
"""
import argparse
import ast
import copy
import csv
import hashlib
import importlib.util
import itertools
import json
import logging
import os
from pathlib import Path
import pickle
import sys
import warnings

import numpy as np
import pandas as pd

D = Path('/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade')
FULL_PINS = '9e63fb73e49c9f8ce700b0abd3e588702d2489c5430059e76c5d5b723eabb46a'
FULL_CAPTURE = 'db919b02f8738da60f85d16e0a7ed2bb53f8eaaf6047ed8f4db912fe4938508c'
SETUP_MANIFEST = '09a077d520e7d6f71716539df4cff1aee8816d34245590df63a35b06fe1cd612'
FUNCTION_SHA = '8aa29d9a4f83ff5629bbdf5e6279a43a374d71ca87fe04fe9db88253cc57d2ca'
LOCALS = ['proband_row', 'proband', 'proband_sex', 'proband_is_male', 'proband_is_female', 'proband_het']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True, allow_nan=False); stream.write('\n')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def excerpt(path, logger):
    text = path.read_text()
    source = next(node for node in ast.parse(text).body if isinstance(node, ast.FunctionDef) and node.name == 'sort_and_rank_variants')
    assert hashlib.sha256(ast.get_source_segment(text, source).encode()).hexdigest() == FUNCTION_SHA
    body = [copy.deepcopy(node) for node in source.body if 501 <= node.lineno <= 524]
    assert len(body) == 1 and isinstance(body[0], ast.If) and body[0].lineno == 501 and body[0].end_lineno == 524
    # Returning locals is observation after the unchanged one original statement.
    end = ast.Return(value=ast.Call(func=ast.Name(id='locals', ctx=ast.Load()), args=[], keywords=[]),
        lineno=525, col_offset=4, end_lineno=525, end_col_offset=19)
    function = ast.FunctionDef(name='original_R2_excerpt', args=ast.arguments(posonlyargs=[],
        args=[ast.arg(arg=name) for name in ['df', 'ped_df', 'fam_name']], kwonlyargs=[], kw_defaults=[], defaults=[]),
        body=body+[end], decorator_list=[], lineno=501, col_offset=0, end_lineno=525, end_col_offset=19)
    namespace = dict(pd=pd, logger=logger)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[])), str(path), 'exec'), namespace)
    return namespace['original_R2_excerpt'], dict(lines=[501, 524], full_function_sha256=FUNCTION_SHA,
        statement_ast_sha256=hashlib.sha256(ast.dump(ast.Module(body=body, type_ignores=[])).encode()).hexdigest())


def frame(values, dtype='object', index=None):
    result = pd.DataFrame({'PROBAND': pd.Series(values, dtype=dtype)})
    result['__oracle_row_position'] = np.arange(len(result), dtype=np.int64)
    if index is not None: result.index = index
    return result


def ped(**changes):
    values = {'#FamilyID': ['F', 'F'], 'IndividualID': ['PROBAND', 'OTHER'], 'Phenotype': ['2', '2'], 'Sex': ['1', '2']}
    values.update(changes)
    return pd.DataFrame(values)


def extras():
    samples = ['0/1', '0/1:10:20', '0/1:8,1:9', '1/1', '1|1', '1:20', '0|1:20',
               '0/0:1', './.:1', '', '.', '１', '1\n1', '\x001', None, pd.NA, np.nan, 1, True, np.inf]
    for dtype, values in [('object', samples), ('string', samples[:14]+[None,pd.NA,np.nan])]:
        yield 'sample_whole_string_' + dtype, frame(values, dtype), ped(), 'F'
    groups = [('empty', []), ('all_null', [None,pd.NA,np.nan]), ('int', [1,0,-1]), ('float', [1.,2.5,np.nan]),
        ('bool', [True,False,None]), ('int_float', [1,2.5,None]), ('bool_int', [True,1,None]),
        ('bool_float', [True,2.5,None]), ('str_int', ['1',1,None]), ('str_bool', ['1',True,None]),
        ('str_float', ['1',1.,np.nan]), ('only_empty', ['',None])]
    for name, values in groups:
        yield 'object_inference_' + name, frame(values), ped(), 'F'
    for dtype, values in [('int64',[1,0]), ('uint64',[2**64-1,1]), ('float64',[1.,np.nan]), ('bool',[True,False]),
        ('Int64',[1,None]), ('UInt64',[1,None]), ('Float64',[1.,None]), ('boolean',[True,None]),
        ('string',['1',None])]:
        yield 'sample_dtype_' + dtype, frame(values,dtype), ped(), 'F'
        yield 'sample_empty_dtype_' + dtype, frame([],dtype), ped(), 'F'
    for n in [0, 2]:
        for pedigree, family, label in [(None,None,'both_none'),(ped(),None,'family_none'),(None,'F','ped_none'),
                                      (ped(),pd.NA,'family_pdNA'),(ped(),np.nan,'family_nan')]:
            yield 'missing_branch_' + label + '_' + str(n), frame(['1']*n), pedigree, family
    for sex in ['1','1.0','M',' MaLe\n','2','2.0','F',' FeMaLe ', '', '0', None, pd.NA, np.nan,
                True, False, 1,2,1.,2.,-0.,np.inf,-np.inf,3.25, '\x1cM\x1f', '\u2003F\u00a0', 'İ','Σ','ΟΣ','K']:
        name = str(len(repr(sex))) + '_' + hashlib.sha256(pickle.dumps(sex,protocol=5)).hexdigest()[:12]
        yield 'sex_' + name, frame(['1',None]), ped(Sex=[sex,'2']), 'F'
    for name, values, family in [('text',['F','F'],'F'),('numbers',[1,1.0],1),('numeric_bool',[True,1],True),
        ('string_number',['1',1],1),('nullable',[pd.NA,'F'],'F'),('missing',[None,'F'],None)]:
        pedigree=ped();pedigree['#FamilyID']=pd.Series(values,dtype=object)
        yield 'family_' + name, frame(['1']), pedigree, family
    phenotypes=['2',2,2.,'2.0',' 2',True,None,pd.NA,np.nan]
    for i, phenotype in enumerate(phenotypes):
        yield 'phenotype_' + str(i),frame(['1']),ped(Phenotype=[phenotype,'1']), 'F'
    for name, values, dtype in [('string',['1','2'],'string'),('Int64',[None,2],'Int64'),
        ('Float64',[None,2.],'Float64'),('UInt64',[0,2],'UInt64'),('boolean',[None,True],'boolean')]:
        pedigree=ped();pedigree['Phenotype']=pd.Series(values,dtype=dtype)
        yield 'phenotype_dtype_' + name,frame(['1']),pedigree,'F'
    for name in ['#FamilyID','Phenotype','IndividualID','Sex']:
        yield 'missing_ped_' + name.replace('#',''),frame(['1']),ped().drop(columns=[name]),'F'
    yield 'missing_two_ped_fields',frame(['1']),ped().drop(columns=['#FamilyID','Phenotype']),'F'
    yield 'first_affected_is_missing_sample',frame(['1']),ped(IndividualID=['ABSENT','PROBAND']),'F'
    yield 'no_affected_even_empty_df',frame([]),ped(Phenotype=['1','1']),'F'
    yield 'empty_pedigree',frame(['1']),ped().iloc[:0],'F'
    for i, identifier in enumerate([None,pd.NA,np.nan,1,True,1.]):
        yield 'nontext_proband_' + str(i),frame(['1']),ped(IndividualID=[identifier,'OTHER']),'F'
    for label in ['PROBAND', '#FamilyID','Phenotype','IndividualID','Sex']:
        if label=='PROBAND':
            data=frame(['1']);data=pd.concat([data,data[['PROBAND']]],axis=1);pedigree=ped()
        else:
            data=frame(['1']);pedigree=ped();pedigree=pd.concat([pedigree,pedigree[[label]]],axis=1)
        yield 'duplicate_' + label.replace('#',''),data,pedigree,'F'
    for name,index in [('duplicate',pd.Index([7,7,-2],name='duplicate')),('range',pd.RangeIndex(10,16,2,name='range')),
        ('multi',pd.MultiIndex.from_tuples([('a',1),('a',1),('z',-2)],names=['left','right']))]:
        data=frame(['1',None,'11'],'string',index);data.attrs['context']='opaque';data.columns.name='sample columns'
        pedigree=ped();pedigree.index=pd.Index([9,9],name='duplicate pedigree')
        yield 'axes_' + name,data,pedigree,'F'
    for count in range(1,6):
        for positions in itertools.combinations(range(5),count):
            values=[True,1,2**64-1,2.5,'1']
            yield 'object_primitive_subset_' + ''.join(map(str,positions)),frame([values[i] for i in positions]+[None,pd.NA,np.nan]),ped(),'F'
    yield 'warning_before_numeric_failure',frame([1],'int64'),ped(Sex=['unknown','2']),'F'
    for n in [0,2]:
        yield 'disabled_branch_ignores_missing_ped_' + str(n),frame([1]*n,'UInt64'),pd.DataFrame(index=[0]),None
    for family,values,label in [(2**53+1,[float(2**53),2**53+1],'int_above_2p53'),
        (float(2**53),[2**53+1,2**53],'float_above_2p53'),(2**64-1,[float(2**64),2**64-1],'uint_max')]:
        pedigree=ped();pedigree['#FamilyID']=pd.Series(values,dtype=object)
        yield 'family_exact_numeric_' + label,frame(['1']),pedigree,family
    for point in range(0x110000):
        if chr(point).isspace():
            yield 'sex_whitespace_' + hex(point),frame(['1']),ped(Sex=[chr(point)+'MaLe'+chr(point),'2']),'F'
    for n in [0,3]:
        for present in [False,True]:
            data=frame((['1',None,'11'] if n else []),'string')
            data.attrs['scope']='caller frame attributes';data.flags.allows_duplicate_labels=False
            pedigree=ped();pedigree.attrs['scope']='pedigree attributes';pedigree.flags.allows_duplicate_labels=False
            yield 'attrs_flags_' + str(n) + '_' + str(present),data,pedigree if present else None,'F' if present else None
    for dtype,values,family in [('float64',[float(2**53),0.],2**53+1),('int64',[2**53+1,0],float(2**53)),
        ('uint64',[2**64-1,0],float(2**64)),('Float64',[float(2**53),None],2**53+1)]:
        data=frame(['1']);data['OTHER']='11';pedigree=ped();pedigree['#FamilyID']=pd.Series(values,dtype=dtype)
        yield 'unported_family_numeric_dtype_' + dtype,data,pedigree,family
    for label,value,family in [('np_int',np.int64(2**53+1),float(2**53)),('np_float',np.float64(2**53),2**53+1),('np_bool',np.bool_(True),1)]:
        data=frame(['1']);data['OTHER']='11';pedigree=ped();pedigree['#FamilyID']=pd.Series([value,0],dtype=object)
        yield 'unported_family_object_' + label,data,pedigree,family
    for label,value,key in [('np_int',np.int64(2**53+1),float(2**53)),('np_float',np.float64(2**53),2**53+1),('np_bool',np.bool_(True),1)]:
        data=frame(['1']);data['OTHER']='11';pedigree=ped();pedigree['#FamilyID']=pd.Series([key,0],dtype=object)
        yield 'unported_family_scalar_' + label,data,pedigree,value


class Logs(logging.Handler):
    def __init__(self, helper):
        super().__init__();self.helper=helper;self.records=[]
    def emit(self, record):
        self.records.append(dict(level=record.levelname,message=record.getMessage(),template=record.msg,args=self.helper.encode(record.args)))


def export(output, name, data, pedigree, family, original, logger, helper, cell, origin, expected=None, expected_error=None):
    folder=output/name;folder.mkdir()
    prepared=pickle.dumps((data,pedigree,family),protocol=5);(folder/'input.pkl').write_bytes(prepared)
    data,pedigree,family=pickle.loads(prepared)
    before_semantics=helper.encode((data,pedigree,family));write(folder/'input.json',before_semantics)
    before=pickle.dumps((data,pedigree,family),protocol=5)
    partial=[]
    def trace(pyframe,event,arg):
        if pyframe.f_code is original.__code__ and event in ['exception','return']:
            partial.append(dict(event=event,line=pyframe.f_lineno,locals={key:pyframe.f_locals[key] for key in LOCALS if key in pyframe.f_locals}))
        return trace
    observer=Logs(helper);logger.addHandler(observer)
    previous=sys.gettrace()
    try:
        with warnings.catch_warnings(record=True) as observed:
            try:
                sys.settrace(trace); result=original(data,pedigree,family)
            except Exception as error:
                result=None;exception=dict(type=type(error).__name__,message=str(error))
            else:exception=None
            finally:sys.settrace(previous)
    finally:logger.removeHandler(observer)
    after=pickle.dumps((data,pedigree,family),protocol=5)
    (folder/'input-after.pkl').write_bytes(after)
    after_semantics=helper.encode((data,pedigree,family));write(folder/'input-after.json',after_semantics)
    write(folder/'source-events.json',[dict(event=event['event'],line=event['line'],
        locals={key:helper.encode(value) for key,value in event['locals'].items()}) for event in partial])
    write(folder/'source-log.json',observer.records)
    report=dict(name=name,origin=origin,rows=len(data),columns=list(data.columns),dtypes=[str(x) for x in data.dtypes],
        pedigree_columns=list(pedigree.columns) if pedigree is not None else None,
        pedigree_dtypes=[str(x) for x in pedigree.dtypes] if pedigree is not None else None,
        family=helper.encode(family),original_exception=exception,input_pickle_unchanged=before==after,
        input_semantically_unchanged=before_semantics==after_semantics,input_before_sha256=hashlib.sha256(before).hexdigest(),
        input_after_sha256=hashlib.sha256(after).hexdigest(),source_log=observer.records,
        first_failure_line=next((event['line'] for event in partial if event['event']=='exception'),None),
        python_warnings=[dict(type=type(item.message).__name__,message=str(item.message),line=item.lineno) for item in observed])
    report['unsupported_family_scalar_identity']=isinstance(family,np.generic)
    report['pedigree_transport_dtypes']=list(report['pedigree_dtypes']) if pedigree is not None else None
    report['unsupported_family_object_scalar_identity']=False
    if pedigree is not None:
        for position,name in enumerate(pedigree.columns):
            if name=='#FamilyID' and str(pedigree.dtypes.iloc[position])=='object' and any(isinstance(value,np.generic) for value in pedigree.iloc[:,position]):
                report['unsupported_family_object_scalar_identity']=True
                report['pedigree_transport_dtypes'][position]='unsupported_numpy_object_family'
    assert report['input_semantically_unchanged'],name
    for label,table in [('input',data),('pedigree',pedigree)]:
        if table is None:continue
        with (folder/(label+'.csv')).open('x',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['index_kind','index_value']+[item for col in table.columns for item in [col+'_kind',col]])
            for row in range(len(table)):
                writer.writerow(['opaque',helper.canonical(helper.encode(table.index[row])).hex()]+[
                    item for col in range(table.shape[1]) for item in cell(table.iloc[row,col])])
    with (folder/'family.csv').open('x',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['kind','value'])
        kind,value=cell(family);writer.writerow(['unsupported_numpy_family_'+kind if isinstance(family,np.generic) else kind,value])
    if result is not None:
        selected={key:result[key] for key in LOCALS if key in result}
        (folder/'expected.pkl').write_bytes(pickle.dumps(selected,protocol=5));write(folder/'expected.json',helper.encode(selected))
        report.update(mask_dtype=str(selected['proband_het'].dtype),mask_name=helper.encode(selected['proband_het'].name))
        if expected is not None:
            assert {key:helper.encode(value) for key,value in selected.items()}=={key:expected[key] for key in LOCALS if key in expected},name
            report['actual_full_ranker_R3_locals_exact']=True
    if expected_error is not None:
        assert exception==expected_error,(name,exception,expected_error)
        report['actual_full_ranker_R2_failure_exact']=True
    return report


def capture(args):
    native={key:os.environ.get(key) for key in ['SLURM_JOB_ID','SLURM_STEP_ID','SLURMD_NODENAME']}
    assert native['SLURM_JOB_ID'] and native['SLURM_STEP_ID']
    for path in [args.full_oracle,args.setup_oracle,args.output]:assert path.resolve().is_relative_to(D)
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'capture-harness.py').write_bytes(Path(__file__).read_bytes())
    write(args.output/'attempt.json',dict(native_environment=native,harness_sha256=sha(Path(__file__))))
    assert sha(args.full_oracle/'source-pins.json')==FULL_PINS and sha(args.full_oracle/'capture.json')==FULL_CAPTURE
    assert sha(args.setup_oracle/'manifest.json')==SETUP_MANIFEST
    helper=load('frozen_full_oracle',args.full_oracle/'source/ranking_control_oracle.py');helper.verify_sources(args.full_oracle,FULL_PINS)
    setup=load('frozen_setup_oracle',args.setup_oracle/'capture-harness.py')
    for name,wanted in json.loads((args.setup_oracle/'manifest.json').read_text()).items():assert sha(args.setup_oracle/name)==wanted
    scoring,_,_,_=helper.original_modules(args.full_oracle/'source/scripts')
    original,source=excerpt(args.full_oracle/'source/scripts/acmg_scoring.py',scoring.logger)
    write(args.output/'runtime-identities.json',helper.runtime_identity())
    # A predicate-only sex contract need not manufacture an exported Python
    # spelling. Bind the Unicode fact it relies on before any Rust code exists.
    write(args.output/'sex-predicate-unicode-witness.json',dict(python_version=sys.version,
        whitespace_codepoints=[point for point in range(0x110000) if chr(point).isspace()],
        nonascii_lower_entirely_recognized_ascii_letters=[point for point in range(128,0x110000)
            if chr(point).lower() and set(chr(point).lower())<=set('malef')],
        semantics='all Unicode scalar codepoints under producing CPython; only male/female predicate letters, not a general lowercase implementation'))
    reports=[];full={case['folder']:case for case in json.loads((args.full_oracle/'capture.json').read_text())['cases']}
    for case in json.loads((args.setup_oracle/'oracle.json').read_text())['cases']:
        if case['origin']['kind']!='actual-closed-full-ranker-R2-entry-or-R1-failure' or case['original_exception'] is not None:continue
        actual=full[case['origin']['full_case']];folder=args.full_oracle/actual['folder']
        _,pedigree=pickle.loads((folder/'prepared-input.pkl').read_bytes())
        bound=json.loads((folder/'traced/arguments.json').read_text())['arguments']
        saved=next(value for key,value in bound['items'] if key==helper.encode('fam_name'))
        family=helper.decode_atom(saved)
        assert helper.encode(family)==saved,'exact recorded family argument required'
        trace=json.loads((folder/'traced/trace.json').read_text())
        events=[event for event in trace['events'] if event['function']=='sort_and_rank_variants' and event['phase']=='R3']
        expected,error=None,None
        if events:expected=json.loads((args.full_oracle/'states'/(events[0]['state']+'.json')).read_text())
        elif actual['first_uncaught_failure'] and actual['first_uncaught_failure']['phase']=='R2':
            error=dict(type=actual['exception']['type'].rsplit('.',1)[-1],message=actual['exception']['message'])
        else:continue
        data=pickle.loads((args.setup_oracle/case['name']/'expected.pkl').read_bytes())
        reports.append(export(args.output,'full_'+actual['name'],data,pedigree,family,original,scoring.logger,helper,setup.cell,
            dict(kind='closed-full-ranker-R3-locals-or-R2-failure',full_case=actual['folder']),expected,error))
    for name,data,pedigree,family in extras():
        reports.append(export(args.output,'extra_'+name,data,pedigree,family,original,scoring.logger,helper,setup.cell,
            dict(kind='verbatim-original-R2-statement')))
    write(args.output/'oracle.json',dict(source_commit=helper.COMMIT,source=source,harness_sha256=sha(Path(__file__)),
        native_environment=native,full_oracle=str(args.full_oracle),full_source_pins_sha256=FULL_PINS,
        full_capture_sha256=FULL_CAPTURE,setup_oracle=str(args.setup_oracle),setup_manifest_sha256=SETUP_MANIFEST,
        cases=reports,scope='R2 context only; no numerical/control/scoring execution',scientific_acceptance=False,production_backend_selected=False))
    write(args.output/'manifest.json',{str(path.relative_to(args.output)):sha(path) for path in sorted(args.output.rglob('*')) if path.is_file()})
    print(json.dumps(dict(cases=len(reports),outputs=sum(case['original_exception'] is None for case in reports),
        source_failures=sum(case['original_exception'] is not None for case in reports),
        pickle_mutations=sum(not case['input_pickle_unchanged'] for case in reports),
        oracle_sha256=sha(args.output/'oracle.json'),manifest_sha256=sha(args.output/'manifest.json'),native_environment=native)),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    capture_parser=commands.add_parser('capture')
    for key in ['full-oracle','setup-oracle','output']:capture_parser.add_argument('--'+key,type=Path,required=True)
    capture(parser.parse_args())


if __name__=='__main__':
    main()
