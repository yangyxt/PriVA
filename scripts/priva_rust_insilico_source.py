"""Admission for the complete, optional PP3/BP4/BP7 column implementation.

The original file, callable code/default objects and mature package profile are
checked before the first getter. Unknown axes/storage/objects execute the whole
original function. This is a fixed-profile bridge, not a pandas replacement.
"""
import ast
import builtins
import gc
import hashlib
import json
import logging
from pathlib import Path
import stat
from types import CodeType, FunctionType

import numpy as np
import pandas as pd

REFERENCE_SHA = 'b8a093a97d0ad9c4730540b2a57a60a302794bf84a42ba304a22e782053a81a2'


def regular_bytes(path):
    path = Path(path)
    before = path.lstat()
    if not stat.S_ISREG(before.st_mode):
        raise ValueError('nonregular source/executable')
    raw = path.read_bytes()
    after = path.lstat()
    fields = ('st_dev', 'st_ino', 'st_mode', 'st_size', 'st_mtime_ns', 'st_ctime_ns', 'st_nlink')
    if any(getattr(before, key) != getattr(after, key) for key in fields):
        raise ValueError('source/executable changed while reading')
    return raw


def reference_guard(reference):
    raw = regular_bytes(reference.__file__)
    if hashlib.sha256(raw).hexdigest() != REFERENCE_SHA:
        raise ValueError('original insilico source changed')
    compiled = compile(raw, str(reference.__file__), 'exec')
    for name in ('PP3_BP4_criteria', 'BP7_criteria'):
        function = getattr(reference, name)
        codes = [code for code in compiled.co_consts if type(code) is CodeType and code.co_name == name]
        if (type(function) is not FunctionType or len(codes) != 1
                or function.__code__ != codes[0] or function.__globals__ is not vars(reference)
                or function.__closure__ is not None or function.__kwdefaults__ is not None):
            raise ValueError('actual original callable changed')
        if name == 'BP7_criteria' and function.__defaults__ is not None:
            raise ValueError('BP7 defaults changed')
        if name == 'PP3_BP4_criteria':
            defaults = function.__defaults__
            if type(defaults) is not tuple or len(defaults) != 2 or defaults[0] is not None:
                raise ValueError('PP3 defaults changed')
            # The original mutable mapping is an actual argument, never replaced
            # by a freshly invented confidence policy. Its contents are handled
            # by the exact same mapping phase as an explicitly supplied dict.
        for label in ('len', 'sorted', 'int', 'str'):
            if function.__builtins__.get(label) is not getattr(builtins, label):
                raise ValueError('original builtin changed')
    if reference.pd is not pd or reference.np is not np:
        raise ValueError('original scientific namespace changed')
    for label in ('len', 'sorted', 'int', 'str'):
        if vars(reference).get(label, getattr(builtins, label)) is not getattr(builtins, label):
            raise ValueError('original builtin shadowed')
    # Reuse the existing COMPLETE, pinned Conda/PyPI mature ownership profile;
    # do not weaken it to a pandas/NumPy version-only assertion. No PM4 code or
    # interval procedure is executed or changed by these read-only checks.
    from priva_rust_pm4_bp3 import _package_profile
    import priva_rust_pm4_bp3_packages as packages
    binding, digest = _package_profile()
    packages.guard(binding, digest)
    if type(reference.logger) is not logging.Logger or 'warning' in vars(reference.logger):
        raise ValueError('custom original logger binding')
    # The original logger calls execute inline, including actual handler and
    # formatter callbacks, their mutations and exceptions. Those mature effects
    # are not represented by Rust, so a Handler/Formatter subclass is not an
    # unsupported calculation profile. Retain only the ordinary logger/method
    # ownership checks; after_log rechecks bindings after each real callback.
    methods = {label.split(':')[1].rsplit('.', 1)[-1]
               for label, _declaration, _mode in packages.SPECS
               if label.startswith('logging:Logger.')} | {'warning'}
    logger = reference.logger
    while logger is not None:
        if type(logger) not in (logging.Logger, logging.RootLogger):
            raise ValueError('custom logger')
        data = object.__getattribute__(logger, '__dict__')
        if any(name in data for name in methods) or type(data.get('name')) is not str:
            raise ValueError('per-instance logger method/name changed')
        if (type(data.get('level')) is not int or type(data.get('disabled')) is not bool
                or type(data.get('handlers')) is not list):
            raise ValueError('custom logger state')
        if not data.get('propagate'): break
        logger = data.get('parent')
    if (np.logical_not is not np.core.umath.logical_not
            or np.asarray is not np.core.multiarray.asarray
            or np.zeros is not np.core.multiarray.zeros):
        raise ValueError('NumPy operation binding changed')
    # Additional slots used by this unit have FIXED canonical declarations;
    # checking an arbitrary replacement against its own source would be unsafe.
    # The two additional pandas source files are pinned beside the unchanged
    # complete ownership profile. Other files retain that profile's hashes.
    expected = json.loads(regular_bytes(binding))['source_files']
    expected.update({'pandas.core.tools.numeric': '7fc1ca5272933483ef96b15ae1b6cbcbaa4c1f750b4a071ed38a91ccae6a57f6',
                     'pandas.core.reshape.concat': 'ab05ec02523da672e5775a63f6ea877f6ce29d777e7e3f06137919f9735569c5'})
    extras = (
        (pd.to_numeric, 'pandas.core.tools.numeric:to_numeric', 'plain'),
        (pd.concat, 'pandas.core.reshape.concat:concat', 'plain'),
        (pd.Series.astype, 'pandas.core.generic:NDFrame.astype', 'plain'),
        (pd.Series.notna, 'pandas.core.series:Series.notna', 'plain'),
        (pd.Series.isna, 'pandas.core.series:Series.isna', 'plain'),
        (pd.Series.map, 'pandas.core.series:Series.map', 'plain'),
        (pd.Series.unique, 'pandas.core.series:Series.unique', 'plain'),
        (pd.Series.__getitem__, 'pandas.core.series:Series.__getitem__', 'plain'),
        (pd.Series.__le__, 'pandas.core.arraylike:OpsMixin.__le__', 'operator'),
        (pd.Series.__ne__, 'pandas.core.arraylike:OpsMixin.__ne__', 'operator'),
        (logging.Logger.warning, 'logging:Logger.warning', 'plain'))
    from pandas.core.strings.accessor import StringMethods
    extras += ((StringMethods.strip, 'pandas.core.strings.accessor:StringMethods.strip', 'contains'),)
    names = {declaration.split(':')[0] for _function, declaration, _mode in extras}
    names.add('pandas.core.ops.common')
    context = packages.source_context(names, expected)
    for function, declaration, mode in extras:
        if mode == 'plain': packages.plain_function(function, declaration, context)
        else: packages.wrapper(function, declaration, mode, context)
    return raw


def placements(placement, count):
    """Read only existing Cython references; never request lazy .as_array."""
    existing = []
    references = gc.get_referents(placement)
    for value in references:
        if value is None or value is type(placement):
            continue
        if type(value) is slice:
            if any(item is not None and type(item) is not int for item in (value.start, value.stop, value.step)):
                raise ValueError('custom placement slice')
            existing.append(list(range(*value.indices(count))))
        elif type(value) is np.ndarray and value.ndim == 1 and value.dtype.kind == 'i':
            existing.append([int(item) for item in value])
        else:
            raise ValueError('unknown placement representation')
    if not 1 <= len(existing) <= 2 or len(references) > 4 or any(item != existing[0] for item in existing):
        raise ValueError('ambiguous placement representation')
    return existing[0]


def raw_columns(frame):
    from pandas.core.internals.managers import BlockManager
    from pandas.core.internals.blocks import NumpyBlock, ExtensionBlock
    if type(frame._mgr) is not BlockManager:
        raise ValueError('custom frame manager')
    count = len(frame._mgr.axes[0])
    result = [None] * count
    for block in frame._mgr.blocks:
        if type(block) not in (NumpyBlock, ExtensionBlock):
            raise ValueError('custom frame block')
        locations = placements(block._mgr_locs, count)
        values = block.values
        if type(values) is np.ndarray:
            if values.ndim != 2 or values.shape[0] != len(locations):
                raise ValueError('custom block dimensions')
            columns = [values[i] for i in range(len(locations))]
        else:
            if len(locations) != 1:
                raise ValueError('ambiguous extension block')
            columns = [values]
        for location, column in zip(locations, columns):
            if not 0 <= location < count or result[location] is not None:
                raise ValueError('invalid/overlapping block placement')
            result[location] = column
    if any(value is None for value in result):
        raise ValueError('incomplete block placement')
    return result


def axis_values(axis):
    """Ordinary unique labels only, without index engine/cache materialization."""
    if type(axis) is pd.RangeIndex:
        return list(vars(axis)['_range'])
    if type(axis) is not pd.Index:
        raise ValueError('unsupported axis type')
    values = vars(axis)['_data']
    if type(values) is not np.ndarray or values.ndim != 1 or values.dtype.kind not in 'OiU':
        raise ValueError('unsupported axis storage')
    result = []
    for value in values:
        if type(value) in (str, int): result.append(value)
        elif type(value) in (np.int64, np.int32): result.append(int(value))
        elif type(value) is np.str_: result.append(str(value))
        else: raise ValueError('unsupported axis label')
    return result


def plain_array(values):
    if type(values) is np.ndarray:
        if values.ndim != 1 or values.dtype.kind not in 'ObifU':
            raise ValueError('unsupported raw array')
        if values.dtype.kind == 'f' and values.dtype != np.dtype('float64'):
            # NumPy float32 scalar promotion changes exact threshold behavior;
            # do not silently substitute f64 comparison for that profile.
            raise ValueError('unsupported floating precision')
        if values.dtype.kind == 'i' and values.dtype != np.dtype('int64'):
            raise ValueError('unsupported integer storage width')
        if values.dtype.kind == 'O':
            for value in values:
                if (type(value) not in (str, bool, int, float, np.int64, np.int32, np.float64, np.float32, np.bool_)
                        and value is not None and value is not pd.NA):
                    raise ValueError('unsupported object scalar')
        return
    if type(values) in (pd.arrays.StringArray, pd.arrays.BooleanArray, pd.arrays.IntegerArray, pd.arrays.FloatingArray):
        if type(values) is pd.arrays.FloatingArray and str(values.dtype) != 'Float64':
            raise ValueError('unsupported nullable floating precision')
        if type(values) is pd.arrays.IntegerArray and str(values.dtype) != 'Int64':
            raise ValueError('unsupported nullable integer width')
        return
    raise ValueError('unsupported extension storage')


def frame_guard(frame, pvs1):
    from priva_rust_pm4_bp3_packages import instance_guard
    if type(frame) is not pd.DataFrame:
        raise ValueError('custom DataFrame has an original-only getter contract')
    instance_guard(frame)
    columns = axis_values(frame._mgr.axes[0])
    labels = axis_values(frame._mgr.axes[1])
    if len(set(columns)) != len(columns) or any(type(name) is not str for name in columns):
        raise ValueError('ambiguous/non-string column axis')
    if len(set(labels)) != len(labels) or len({type(label) for label in labels}) > 1:
        raise ValueError('ambiguous/mixed row index')
    for values in raw_columns(frame): plain_array(values)
    for cached in frame._item_cache.values():
        if type(cached) is not pd.Series: raise ValueError('custom cached Series')
        instance_guard(cached, series=True)
    if pvs1 is not None:
        if type(pvs1) is pd.Series:
            instance_guard(pvs1, series=True)
            pvs_labels = axis_values(pvs1._mgr.axes[0])
            if len(set(pvs_labels)) != len(pvs_labels) or set(pvs_labels) != set(labels):
                raise ValueError('unsupported PVS1 alignment shape')
            plain_array(pvs1._mgr.blocks[0].values)
        elif type(pvs1) is np.ndarray:
            plain_array(pvs1)
            if pvs1.shape != (len(labels),): raise ValueError('unsupported PVS1 shape')
        else:
            raise ValueError('unsupported PVS1 object')
    return columns, labels


def source_tree(raw):
    tree = ast.parse(raw)
    return [node for node in tree.body if type(node) is ast.FunctionDef
            and node.name in ('PP3_BP4_criteria', 'BP7_criteria')]
