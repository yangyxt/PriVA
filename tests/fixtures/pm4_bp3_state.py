"""Exact typed state encoding without pandas column reads or object callbacks."""
import base64
import hashlib
import json
import struct
import numpy as np
import pandas as pd


def kind(value):
    cls = type(value)
    return cls.__module__ + '.' + cls.__qualname__


def atom(value):
    label = kind(value)
    if value is pd.NA: return dict(type=label, value='pd.NA')
    if value is None: return dict(type=label, value=None)
    if type(value) in (str, bool): return dict(type=label, value=value)
    if type(value) is int: return dict(type=label, decimal=str(value))
    if type(value) is float: return dict(type=label, f64_be=struct.pack('>d', value).hex())
    if isinstance(value, np.generic):
        return dict(type=label, dtype=value.dtype.str, bytes=value.tobytes().hex())
    if type(value) in (list, tuple): return dict(type=label, items=[atom(item) for item in value])
    if type(value) is dict: return dict(type=label, items=[[atom(key), atom(item)] for key, item in value.items()])
    if type(value) in (set, frozenset):
        items=[atom(item) for item in value]
        return dict(type=label, items=sorted(items, key=lambda item: json.dumps(item, sort_keys=True)))
    # No repr(), str(), len(), dtype or custom protocol is invoked here.
    return dict(type=label, unsupported_object_identity_not_serialized=True)


def array(value):
    label = kind(value)
    if type(value) is np.ndarray:
        result = dict(type=label, dtype=value.dtype.str, shape=list(value.shape), strides=list(value.strides),
                      writeable=bool(value.flags.writeable), c_contiguous=bool(value.flags.c_contiguous),
                      f_contiguous=bool(value.flags.f_contiguous))
        if value.dtype.kind == 'O':
            digest=hashlib.sha256(); small=[]
            for item in value.flat:
                encoded=atom(item)
                raw=json.dumps(encoded, sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode()+b'\n'
                digest.update(raw)
                if value.size <= 32: small.append(encoded)
            result.update(typed_ordered_element_sha256=digest.hexdigest(), elements=small if value.size <= 32 else None)
        else:
            raw=value.tobytes(order='C')
            result.update(raw_C_order_base64=base64.b64encode(raw).decode(), raw_sha256=hashlib.sha256(raw).hexdigest())
        return result
    if type(value) in (pd.arrays.BooleanArray, pd.arrays.IntegerArray, pd.arrays.FloatingArray):
        return dict(type=label, dtype=str(value.dtype), data=array(value._data), mask=array(value._mask))
    if type(value) is pd.arrays.StringArray:
        return dict(type=label, dtype=str(value.dtype), ndarray=array(value._ndarray))
    if type(value) is pd.Categorical:
        return dict(type=label, codes=array(value._codes), categories=index(value.categories), ordered=value.ordered)
    return dict(type=label, unsupported_storage_not_dereferenced=True)


def index(value):
    result=dict(type=kind(value), names=atom(list(value.names)), values=array(value._values))
    if type(value) is pd.RangeIndex:
        result['range']=[value.start, value.stop, value.step]
    if type(value) is pd.MultiIndex:
        result['levels']=[index(level) for level in value.levels]
        result['codes']=[array(code) for code in value.codes]
        result['sortorder']=value.sortorder
    return result


def series(value):
    return dict(type=kind(value), name=atom(value.name), index=index(value.index),
        values=array(value._mgr.blocks[0].values), attrs=atom(value.attrs),
        allows_duplicate_labels=value.flags.allows_duplicate_labels)


def frame(value):
    cache=value._item_cache
    before=[(key,id(item)) for key,item in cache.items()]
    result=dict(type=kind(value), index=index(value.index), columns=index(value.columns),
        arrays=[array(value._mgr.iget_values(position)) for position in range(len(value.columns))],
        attrs=atom(value.attrs), allows_duplicate_labels=value.flags.allows_duplicate_labels,
        cache=[dict(key=atom(key), value=series(item) if isinstance(item,pd.Series) else atom(item)) for key,item in cache.items()])
    assert before==[(key,id(item)) for key,item in cache.items()], 'Observer altered the input cache'
    return result


def output(value):
    if isinstance(value, pd.Series): return series(value)
    if isinstance(value, pd.DataFrame): return frame(value)
    if isinstance(value, np.ndarray): return array(value)
    if type(value) is tuple: return dict(type=kind(value), items=[output(item) for item in value])
    return atom(value)
