#!/usr/bin/env python3
"""Capture unchanged09223ef mask behavior; typed CSV is an explicit test transport."""
import argparse, csv, hashlib, importlib, json, math, platform, random, struct, sys
from pathlib import Path
import pandas as pd
import numpy as np
Y=Path('/stornext/snfs170/next-gen/scratch/yangyxt')
SOURCE=Y/'softwares/priva-reference-09223ef/scripts'
COMMIT='09223ef5a8b0ceb4c6d61055c6eefe62b51ba903'
COLUMNS=['var_plausible_patho_mechs','variant_effect','variant_lof_score','variant_gof_score','variant_dn_score']

def encode(v):
    if v is None or v is pd.NA: return ['null','']
    if isinstance(v,(bool,np.bool_)): return ['bool',str(int(v))]
    if isinstance(v,(int,np.integer)): return ['uint' if v > 2**63-1 else 'int',str(v)]
    if isinstance(v,(float,np.floating)): return ['float',struct.pack('>d',v).hex()]
    if isinstance(v,str): return ['str',v.encode().hex()]
    raise TypeError(type(v))

def cases():
    histories=['',' ',None,float('nan'),'\u001c\u001f','\u00a0\u2003','\u200b']
    for base in ['recessive','dominant','x_linked_recessive','x_linked_dominant','x_linked_unspecified']:
        for suffix in ['', '_LOF','_GOF','_DN','_','_other','_\n','_LOF\n','_LOF\n\n']:
            for outer in ['{}',';{};',' {}','{} ','{}\n','{}\n\n','x;{};y','x{}']:
                histories.append(outer.format(base+suffix))
    effects=['','uncertain',' uncertain','uncertain\n','predicted_LOF_high_confidence','predicted_LOF_high_confidence\n',None]
    for token in ['GOF','DOMINANT_NEGATIVE','DN','LOF','gof']:
        for template in ['{}','exact_known_{}','exact_known_{}\n','exact_known_{}\n\n','x+{}','+{}+x','x{}','exact_known_X+{}+Y']:
            effects.append(template.format(token))
    values=[None,pd.NA,float('nan'),float('inf'),-float('inf'),True,False,-1,0,1,2,3,0.9999999999999999,1.0000000000000002,1.9999999999999998,2.0000000000000004]
    values += ['', ' ','\t2\n','\u00a02','2\u00a0','2_0','2,0','+2','-0','2.','02','2e0','2E+00','0x2','NaN','nan','NAN','inf','Infinity','-INF','+infinity','1e309','-1e309','1e-999','2\x00junk','２','٢','1.0000000000000001','0.99999999999999999','1.9999999999999999','2.0000000000000001']
    rng=random.Random(19223)
    for n in range(1,35):
        values.extend(['0.'+'9'*n,'1.'+'0'*n+'1','1.'+'9'*n,'2.'+'0'*n+'1'])
    for _ in range(1500):
        values.append(rng.choice(['0.','1.','2.'])+''.join(str(rng.randrange(10)) for _ in range(rng.randrange(1,28))))
    def frame(profiles,effects,lof,gof,dn):
        d=pd.DataFrame(dict(zip(COLUMNS,[profiles,effects,lof,gof,dn])),dtype=object)
        d.index=pd.Index([7 if i%3 else -2 for i in range(len(d))],name='gapped_duplicate')
        return d
    yield 'profiles',frame(histories,['uncertain']*len(histories),[1]*len(histories),[0]*len(histories),[0]*len(histories))
    yield 'effects',frame(['dominant_GOF']*len(effects),effects,[1]*len(effects),[0]*len(effects),[0]*len(effects))
    yield 'numeric',frame(['']*len(values),['']*len(values),values,values,values)
    for label,values2 in [('signed_unsigned_mix',[2**63, -1, '2', None]),('integer_invalid_mix',[2**63-1, 'bad', '1.9999999999999999', float('nan')]),('score_neighbors',[np.nextafter(1.,0.),1.,np.nextafter(1.,2.),np.nextafter(2.,1.),2.,np.nextafter(2.,3.)])]:
        yield label,frame(['']*len(values2),['']*len(values2),values2,values2,values2)
    for label,vals in [('integer_only',['0'*17+'2']),('integer_bad',['0'*17+'2','bad']),('integer_null',['0'*17+'2',None]),('integer_float',['0'*17+'2',0.0]),('nul_forms',['2\0junk','2.0\0junk','2e0\0junk','2 \0junk','inf\0junk'])]:
        yield label,frame(['']*len(vals),['']*len(vals),vals,vals,vals)
    vals=['0'*17+'2',pd.NA]
    d=frame(['',''],['',''],vals,vals,vals)
    for key in COLUMNS[2:]:d[key]=pd.array(vals,dtype='string')
    yield 'nullable_integer_string',d
    for label,vals,dtype in [('nullable_integer_empty',['0'*17+'2',''],'string'),('nullable_integer_bad',['0'*17+'2','bad'],'string'),('nullable_bool_null',[True,pd.NA],'boolean'),('nullable_bool_valid',[True,False],'boolean'),('nullable_int_score',[2,pd.NA],'Int64'),('nullable_float_score',[2.0,pd.NA],'Float64')]:
        d=frame(['']*len(vals),['']*len(vals),vals,vals,vals)
        for key in COLUMNS[2:]:d[key]=pd.array(vals,dtype=dtype)
        yield label,d
    # Each score phase independently exercises actual column-wide dtype/null rules.
    matrix={
      'object':[['0'*17+'2'],['0'*17+'2',None],['0'*17+'2',''],['0'*17+'2',' '],['0'*17+'2','bad'],['0'*17+'2',float('nan')],[True,None],[2**63,-1,'2'],['2\0junk','2.0\0junk']],
      'string':[['0'*17+'2'],['0'*17+'2',pd.NA],['0'*17+'2',''],['0'*17+'2',' '],['0'*17+'2','bad'],['0'*17+'2','nan'],['0'*17+'2','NaN'],['2\0junk','2.0\0junk']],
      'Int64':[[2],[2,pd.NA],[0,1,2,-1],[pd.NA,pd.NA]],
      'Float64':[[2.0],[2.0,pd.NA],[0.9999999999999999,1.0,1.9999999999999998,2.0],[float('inf'),-float('inf')],[pd.NA,pd.NA]],
      'boolean':[[True],[True,False],[True,pd.NA],[pd.NA,pd.NA]],
    }
    for dtype,groups in matrix.items():
        for number,vals in enumerate(groups):
            for phase,key in enumerate(COLUMNS[2:]):
                d=pd.DataFrame({COLUMNS[0]:['']*len(vals),COLUMNS[1]:['']*len(vals)})
                d[key]=pd.Series(vals,dtype=dtype)
                yield f'matrix_{dtype}_{number}_{phase}',d
    for earlier in range(2):
        d=pd.DataFrame({COLUMNS[0]:[''],COLUMNS[1]:['']})
        for key in COLUMNS[2+earlier:]:d[key]=pd.Series([pd.NA],dtype='boolean')
        yield f'matrix_score_failure_order_{earlier}',d
    yield 'text_scalars',frame([0,1,True,False,2.5,None],[None,True,False,1,2.5,float('nan')],[None]*6,[None]*6,[None]*6)
    for bits in range(8):
        d=frame(['dominant'],['exact_known_GOF'],[2],[2],[2])
        for i,col in enumerate(COLUMNS[2:]):
            if not bits&(1<<i):del d[col]
        yield f'optional_{bits}',d
    for name,cols in [('missing_both',[]),('missing_profile',[COLUMNS[1]]),('missing_effect',[COLUMNS[0]])]:
        yield name,pd.DataFrame({k:[''] for k in cols})
    for name,cols in [('duplicate_profile',[COLUMNS[0],COLUMNS[0],COLUMNS[1]]),('duplicate_effect',[COLUMNS[0],COLUMNS[1],COLUMNS[1]]),('duplicate_lof',[COLUMNS[0],COLUMNS[1],COLUMNS[2],COLUMNS[2]]),('duplicate_effect_lof',[COLUMNS[0],COLUMNS[1],COLUMNS[1],COLUMNS[2],COLUMNS[2]])]:
        yield name,pd.DataFrame([['']*len(cols)],columns=cols)
    d=pd.concat([pd.Series([pd.NA],dtype='Int64'),pd.Series([pd.NA],dtype='Int64'),pd.Series([''])],axis=1)
    d.columns=[COLUMNS[0],COLUMNS[0],COLUMNS[1]]
    yield 'duplicate_nullable_profile',d
    d=pd.concat([pd.Series(['']),pd.Series([pd.NA],dtype='Int64'),pd.Series([pd.NA],dtype='Int64'),pd.Series(['']),pd.Series([''])],axis=1)
    d.columns=[COLUMNS[0],COLUMNS[1],COLUMNS[1],COLUMNS[2],COLUMNS[2]]
    yield 'duplicate_nullable_effect_lof',d
    yield 'empty',pd.DataFrame({k:pd.Series(dtype=object) for k in COLUMNS})
    for dtype in ['string','Float64','Int64','boolean']:
        d=pd.DataFrame({COLUMNS[0]:pd.Series([None,None],dtype=dtype),COLUMNS[1]:['uncertain','']})
        yield 'nullable_text_'+dtype,d

def capture(args):
    args.output.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(SOURCE)); module=importlib.import_module('acmg_variant_mechanism')
    assert hashlib.sha256((SOURCE/'acmg_variant_mechanism.py').read_bytes()).hexdigest() == 'c7e8da80ca0ec47fb38a5fd1d965641c1791189cb0d25a651253f44508d0bf26'
    result={'harness_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source_commit':COMMIT,'source_path':str(SOURCE/'acmg_variant_mechanism.py'),'source_sha256':hashlib.sha256((SOURCE/'acmg_variant_mechanism.py').read_bytes()).hexdigest(),'runtime':{'python':platform.python_version(),'pandas':pd.__version__,'numpy':np.__version__},'cases':[],'scientific_acceptance':False}
    for name,df in cases():
        before=df.copy(deep=True); entry={'name':name,'rows':len(df),'dtypes':{k:str(v) for k,v in df.dtypes.items()},'index_name':df.index.name,'index_dtype':str(df.index.dtype)}
        with (args.output/(name+'.input.csv')).open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['index_kind','index_value']+[v for c in df.columns for v in (c+'_kind',c)])
            for i,row in enumerate(df.itertuples(index=False,name=None)):
                w.writerow(encode(df.index[i])+[x for v in row for x in encode(v)])
        try:
            masks=module._variant_mechanism_masks(df)
            entry['keys']=list(masks)
            entry['series_metadata']={k:{'name':v.name,'dtype':str(v.dtype),'index_names':list(v.index.names),'index_dtype':str(v.index.dtype)} for k,v in masks.items()}
            with (args.output/(name+'.expected.csv')).open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['index_kind','index_value',*masks])
                for i in range(len(df)):w.writerow(encode(df.index[i])+[int(m.iloc[i]) for m in masks.values()])
            assert all(m.index.equals(df.index) and str(m.dtype)=='bool' for m in masks.values())
        except Exception as e:entry['exception']={'type':type(e).__name__,'message':str(e)}
        pd.testing.assert_frame_equal(before,df,check_exact=True);entry['input_unchanged']=True
        result['cases'].append(entry)
    (args.output/'oracle.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True);c=s.add_parser('capture');c.add_argument('--output',type=Path,required=True)
    a=p.parse_args();capture(a)
if __name__=='__main__':main()
