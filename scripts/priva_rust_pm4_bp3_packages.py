"""Canonical source-derived package admission for optional PM4/BP3 procedures.

No mutable API is trusted because it happened to be present at import. Each
Python callable is checked against its named installed-source definition, with
source-derived defaults and the two recognized pandas decorator closures. C
constructors/types are checked against the pinned extension and immutable C
metadata. No dataframe/accessor/logging method is invoked during this proof.
"""
import argparse
import ast
import builtins
from functools import lru_cache
import hashlib
import importlib
import inspect
import json
import logging
from pathlib import Path
import stat
import struct
import sys
from types import BuiltinFunctionType, CodeType, FunctionType, ModuleType

import numpy as np
import pandas as pd


# Public slot -> actual source declaration. These names are a reviewed source
# ledger, not values or code identities captured from the current live slot.
SPECS = (
    ('pandas.core.frame:DataFrame.__init__', 'pandas.core.frame:DataFrame.__init__', 'plain'),
    ('pandas.core.frame:DataFrame.__getitem__', 'pandas.core.frame:DataFrame.__getitem__', 'plain'),
    ('pandas.core.frame:DataFrame._get_item_cache', 'pandas.core.frame:DataFrame._get_item_cache', 'plain'),
    ('pandas.core.frame:DataFrame._ixs', 'pandas.core.frame:DataFrame._ixs', 'plain'),
    ('pandas.core.frame:DataFrame.__len__', 'pandas.core.frame:DataFrame.__len__', 'plain'),
    ('pandas.core.frame:DataFrame.assign', 'pandas.core.frame:DataFrame.assign', 'plain'),
    ('pandas.core.generic:NDFrame.get', 'pandas.core.generic:NDFrame.get', 'plain'),
    ('pandas.core.generic:NDFrame.__finalize__', 'pandas.core.generic:NDFrame.__finalize__', 'plain'),
    ('pandas.core.generic:NDFrame.attrs', 'pandas.core.generic:NDFrame.attrs', 'property'),
    ('pandas.core.generic:NDFrame.flags', 'pandas.core.generic:NDFrame.flags', 'property'),
    ('pandas.core.generic:NDFrame.__getattr__', 'pandas.core.generic:NDFrame.__getattr__', 'plain'),
    ('pandas.core.generic:NDFrame.__setattr__', 'pandas.core.generic:NDFrame.__setattr__', 'plain'),
    ('pandas.core.series:Series.__init__', 'pandas.core.series:Series.__init__', 'plain'),
    ('pandas.core.series:Series._constructor', 'pandas.core.series:Series._constructor', 'property'),
    ('pandas.core.series:Series.dtype', 'pandas.core.series:Series.dtype', 'property'),
    ('pandas.core.series:Series.name', 'pandas.core.series:Series.name', 'property'),
    ('pandas.core.series:Series._values', 'pandas.core.series:Series._values', 'property'),
    ('pandas.core.frame:DataFrame._constructor', 'pandas.core.frame:DataFrame._constructor', 'property'),
    ('pandas.core.series:Series.fillna', 'pandas.core.generic:NDFrame.fillna', 'plain'),
    ('pandas.core.series:Series.isin', 'pandas.core.series:Series.isin', 'plain'),
    ('pandas.core.series:Series.__and__', 'pandas.core.arraylike:OpsMixin.__and__', 'operator'),
    ('pandas.core.series:Series.__rand__', 'pandas.core.arraylike:OpsMixin.__rand__', 'operator'),
    ('pandas.core.series:Series.__or__', 'pandas.core.arraylike:OpsMixin.__or__', 'operator'),
    ('pandas.core.series:Series.__ror__', 'pandas.core.arraylike:OpsMixin.__ror__', 'operator'),
    ('pandas.core.series:Series.__invert__', 'pandas.core.generic:NDFrame.__invert__', 'plain'),
    ('pandas.core.series:Series.__lt__', 'pandas.core.arraylike:OpsMixin.__lt__', 'operator'),
    ('pandas.core.series:Series.__gt__', 'pandas.core.arraylike:OpsMixin.__gt__', 'operator'),
    ('pandas.core.series:Series.__ge__', 'pandas.core.arraylike:OpsMixin.__ge__', 'operator'),
    ('pandas.core.series:Series.__eq__', 'pandas.core.arraylike:OpsMixin.__eq__', 'operator'),
    ('pandas.core.series:Series.__add__', 'pandas.core.arraylike:OpsMixin.__add__', 'operator'),
    ('pandas.core.series:Series.__sub__', 'pandas.core.arraylike:OpsMixin.__sub__', 'operator'),
    ('pandas.core.series:Series.sum', 'pandas.core.series:Series.sum', 'plain'),
    ('pandas.core.series:Series.any', 'pandas.core.series:Series.any', 'plain'),
    ('pandas.core.series:Series.__array__', 'pandas.core.series:Series.__array__', 'plain'),
    ('pandas.core.series:Series._logical_method', 'pandas.core.series:Series._logical_method', 'plain'),
    ('pandas.core.series:Series._cmp_method', 'pandas.core.series:Series._cmp_method', 'plain'),
    ('pandas.core.series:Series._arith_method', 'pandas.core.series:Series._arith_method', 'plain'),
    ('pandas.core.series:Series._align_for_op', 'pandas.core.series:Series._align_for_op', 'plain'),
    ('pandas.core.generic:NDFrame.sum', 'pandas.core.generic:NDFrame.sum', 'plain'),
    ('pandas.core.generic:NDFrame._min_count_stat_function', 'pandas.core.generic:NDFrame._min_count_stat_function', 'plain'),
    ('pandas.core.strings.accessor:StringMethods.__init__', 'pandas.core.strings.accessor:StringMethods.__init__', 'plain'),
    ('pandas.core.strings.accessor:StringMethods._validate', 'pandas.core.strings.accessor:StringMethods._validate', 'staticmethod'),
    ('pandas.core.strings.accessor:StringMethods.len', 'pandas.core.strings.accessor:StringMethods.len', 'plain'),
    ('pandas.core.strings.accessor:StringMethods.contains', 'pandas.core.strings.accessor:StringMethods.contains', 'contains'),
    ('pandas.core.strings.accessor:StringMethods._wrap_result', 'pandas.core.strings.accessor:StringMethods._wrap_result', 'plain'),
    ('pandas.core.strings.object_array:ObjectStringArrayMixin._str_len', 'pandas.core.strings.object_array:ObjectStringArrayMixin._str_len', 'plain'),
    ('pandas.core.strings.object_array:ObjectStringArrayMixin._str_contains', 'pandas.core.strings.object_array:ObjectStringArrayMixin._str_contains', 'plain'),
    ('pandas.core.strings.object_array:ObjectStringArrayMixin._str_map', 'pandas.core.strings.object_array:ObjectStringArrayMixin._str_map', 'plain'),
    ('pandas.core.accessor:CachedAccessor.__init__', 'pandas.core.accessor:CachedAccessor.__init__', 'plain'),
    ('pandas.core.accessor:CachedAccessor.__get__', 'pandas.core.accessor:CachedAccessor.__get__', 'plain'),
    ('pandas.core.arrays.masked:BaseMaskedArray.__init__', 'pandas.core.arrays.masked:BaseMaskedArray.__init__', 'plain'),
    ('pandas.core.arrays.masked:BaseMaskedArray.fillna', 'pandas.core.arrays.masked:BaseMaskedArray.fillna', 'plain'),
    ('pandas.core.arrays.masked:BaseMaskedArray.__array__', 'pandas.core.arrays.masked:BaseMaskedArray.__array__', 'plain'),
    ('pandas.core.arrays.masked:BaseMaskedArray.__invert__', 'pandas.core.arrays.masked:BaseMaskedArray.__invert__', 'plain'),
    ('pandas.core.arrays.boolean:BooleanArray.__init__', 'pandas.core.arrays.boolean:BooleanArray.__init__', 'plain'),
    ('pandas.core.arrays.boolean:BooleanArray._logical_method', 'pandas.core.arrays.boolean:BooleanArray._logical_method', 'plain'),
    ('pandas.core.arrays.integer:IntegerArray.__init__', 'pandas.core.arrays.numeric:NumericArray.__init__', 'plain'),
    ('pandas.core.arrays.string_:StringArray.__init__', 'pandas.core.arrays.string_:StringArray.__init__', 'plain'),
    ('pandas.core.arrays.string_:StringArray.fillna', 'pandas.core.arrays._mixins:NDArrayBackedExtensionArray.fillna', 'plain'),
    ('pandas.core.arrays.string_:StringArray._str_map', 'pandas.core.arrays.string_:StringArray._str_map', 'plain'),
    ('pandas.core.arrays.numpy_:NumpyExtensionArray.__init__', 'pandas.core.arrays.numpy_:NumpyExtensionArray.__init__', 'plain'),
    ('pandas.core.arrays.numpy_:NumpyExtensionArray.fillna', 'pandas.core.arrays._mixins:NDArrayBackedExtensionArray.fillna', 'plain'),
    ('pandas.core.arrays.numpy_:NumpyExtensionArray.__array__', 'pandas.core.arrays.numpy_:NumpyExtensionArray.__array__', 'plain'),
    ('pandas.core.internals.managers:BlockManager.iget_values', 'pandas.core.internals.managers:BlockManager.iget_values', 'plain'),
    ('pandas.core.internals.blocks:Block.iget', 'pandas.core.internals.blocks:Block.iget', 'plain'),
    ('pandas.core.ops.array_ops:logical_op', 'pandas.core.ops.array_ops:logical_op', 'plain'),
    ('pandas.core.ops.array_ops:arithmetic_op', 'pandas.core.ops.array_ops:arithmetic_op', 'plain'),
    ('pandas.core.ops.array_ops:comparison_op', 'pandas.core.ops.array_ops:comparison_op', 'plain'),
    ('pandas.core.roperator:rand_', 'pandas.core.roperator:rand_', 'plain'),
    ('pandas.core.roperator:ror_', 'pandas.core.roperator:ror_', 'plain'),
    ('logging:Logger.info', 'logging:Logger.info', 'plain'),
    ('logging:Logger.error', 'logging:Logger.error', 'plain'),
    ('logging:Logger.isEnabledFor', 'logging:Logger.isEnabledFor', 'plain'),
    ('logging:Logger.getEffectiveLevel', 'logging:Logger.getEffectiveLevel', 'plain'),
    ('logging:Logger._log', 'logging:Logger._log', 'plain'),
    ('logging:Logger.findCaller', 'logging:Logger.findCaller', 'plain'),
    ('logging:Logger.makeRecord', 'logging:Logger.makeRecord', 'plain'),
    ('logging:Logger.handle', 'logging:Logger.handle', 'plain'),
    ('logging:Logger.callHandlers', 'logging:Logger.callHandlers', 'plain'),
    ('logging:Handler.handle', 'logging:Handler.handle', 'plain'),
    ('logging:Handler.filter', 'logging:Filterer.filter', 'plain'),
    ('logging:Handler.format', 'logging:Handler.format', 'plain'),
    ('logging:StreamHandler.emit', 'logging:StreamHandler.emit', 'plain'),
    ('logging:Formatter.format', 'logging:Formatter.format', 'plain'),
    ('logging:Formatter.formatMessage', 'logging:Formatter.formatMessage', 'plain'),
    ('logging:Formatter.formatTime', 'logging:Formatter.formatTime', 'plain'),
    ('logging:LogRecord.__init__', 'logging:LogRecord.__init__', 'plain'),
    ('logging:LogRecord.getMessage', 'logging:LogRecord.getMessage', 'plain'),
    ('logging:_is_internal_frame', 'logging:_is_internal_frame', 'plain'),
)


def regular_bytes(path):
    path = Path(path); before = path.lstat()
    if not stat.S_ISREG(before.st_mode): raise ValueError('nonregular package source')
    raw = path.read_bytes(); after = path.lstat()
    fields = ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns', 'st_mode')
    if any(getattr(before, name) != getattr(after, name) for name in fields): raise ValueError('package source race')
    return raw


def resolve(label):
    module_name, name = label.split(':'); owner = importlib.import_module(module_name)
    for part in name.split('.')[:-1]: owner = inspect.getattr_static(owner, part)
    return owner, name.split('.')[-1]


def slot(label):
    owner, name = resolve(label); return inspect.getattr_static(owner, name)


def executable(code):
    names = ('co_argcount', 'co_posonlyargcount', 'co_kwonlyargcount', 'co_nlocals',
        'co_stacksize', 'co_flags', 'co_code', 'co_names', 'co_varnames', 'co_filename',
        'co_name', 'co_qualname', 'co_firstlineno', 'co_exceptiontable', 'co_freevars', 'co_cellvars')
    return tuple(getattr(code, name) for name in names), tuple(executable(value) if type(value) is CodeType else value for value in code.co_consts)


def codes(code):
    yield code
    for value in code.co_consts:
        if type(value) is CodeType: yield from codes(value)


@lru_cache(maxsize=64)
def parsed(raw, filename): return ast.parse(raw), compile(raw, filename, 'exec')


def definitions(tree):
    result = {}
    def walk(body, prefix=()):
        for node in body:
            if type(node) is ast.ClassDef: walk(node.body, prefix + (node.name,))
            elif type(node) in (ast.FunctionDef, ast.AsyncFunctionDef):
                label = '.'.join(prefix + (node.name,)); result.setdefault(label, []).append(node)
                walk(node.body, prefix + (node.name, '<locals>'))
            elif type(node) in (ast.If, ast.Try):
                walk(node.body, prefix); walk(node.orelse, prefix)
                if type(node) is ast.Try:
                    walk(node.finalbody, prefix)
                    for handler in node.handlers: walk(handler.body, prefix)
    walk(tree.body); return result


def typed(value):
    from pandas._libs import lib
    if type(value) is lib._NoDefault and value is lib._NoDefault.no_default: return ['pandas_no_default']
    if type(value) in (type(None), str, bool, int): return [type(value).__name__, value]
    if type(value) is float: return ['float', struct.pack('>d', value).hex()]
    if type(value) in (tuple, list, set, frozenset):
        values = [typed(item) for item in value]
        if type(value) in (set, frozenset): values.sort(key=lambda item: json.dumps(item, sort_keys=True))
        return [type(value).__name__, values]
    if type(value) is dict: return ['dict', [[typed(key), typed(item)] for key, item in value.items()]]
    if type(value) is type and value.__module__ == 'builtins' and value.__name__ in ('dict', 'list', 'bool', 'int', 'float', 'str', 'object'):
        if value is not getattr(builtins, value.__name__): raise ValueError('builtin class default changed')
        return ['builtin_class', value.__name__]
    raise ValueError('unsupported package default or closure value')


def source_default(node):
    if type(node) is ast.Constant: return typed(node.value)
    if type(node) is ast.UnaryOp and type(node.op) in (ast.USub, ast.UAdd) and type(node.operand) is ast.Constant:
        value = node.operand.value
        if type(value) not in (int, float): raise ValueError('nonliteral signed default')
        return typed(-value if type(node.op) is ast.USub else value)
    if type(node) in (ast.Tuple, ast.List):
        return ['tuple' if type(node) is ast.Tuple else 'list', [source_default(item) for item in node.elts]]
    if type(node) is ast.Dict:
        return ['dict', [[source_default(key), source_default(value)] for key, value in zip(node.keys, node.values)]]
    if type(node) is ast.Name and node.id in ('dict', 'list', 'bool', 'int', 'float', 'str', 'object'):
        return ['builtin_class', node.id]
    if type(node) is ast.Attribute and type(node.value) is ast.Name:
        if (node.value.id, node.attr) == ('lib', 'no_default'): return ['pandas_no_default']
        if (node.value.id, node.attr) == ('np', 'nan'): return ['float', '7ff8000000000000']
    raise ValueError('unrecognized source default expression')


def source_context(names, expected=None):
    result = {}
    for name in sorted(names):
        module = importlib.import_module(name)
        if type(module) is not ModuleType: raise ValueError('package module changed')
        raw = regular_bytes(module.__file__); digest = hashlib.sha256(raw).hexdigest()
        if expected is not None and expected.get(name) != digest: raise ValueError('package source changed: ' + name)
        tree, code = parsed(raw, str(module.__file__))
        result[name] = dict(module=module, raw=raw, code=code, definitions=definitions(tree), sha256=digest)
    return result


def plain_function(function, declaration, context, *, property_getter=False):
    name, qualname = declaration.split(':'); data = context[name]
    nodes = data['definitions'].get(qualname, [])
    nodes = [node for node in nodes if not any(type(value) is ast.Name and value.id == 'overload' for value in node.decorator_list)]
    if property_getter:
        nodes = [node for node in nodes if any(type(value) is ast.Name and value.id == 'property' for value in node.decorator_list)]
    if len(nodes) != 1 or type(nodes[0]) is not ast.FunctionDef:
        raise ValueError('ambiguous canonical declaration: ' + declaration)
    node = nodes[0]
    firstline = min([node.lineno] + [value.lineno for value in node.decorator_list])
    expected = [code for code in codes(data['code']) if code.co_qualname == qualname and code.co_firstlineno == firstline]
    if len(expected) != 1: raise ValueError('ambiguous canonical source code: ' + declaration)
    if (type(function) is not FunctionType or function.__globals__ is not data['module'].__dict__
            or executable(function.__code__) != executable(expected[0])):
        raise ValueError('package executable function changed: ' + declaration)
    positional = ['NoneType', None] if not node.args.defaults else ['tuple', [source_default(value) for value in node.args.defaults]]
    keyword = [(arg.arg, value) for arg, value in zip(node.args.kwonlyargs, node.args.kw_defaults) if value is not None]
    keywords = ['NoneType', None] if not keyword else ['dict', [[typed(key), source_default(value)] for key, value in keyword]]
    if typed(function.__defaults__) != positional or typed(function.__kwdefaults__) != keywords:
        raise ValueError('package source defaults changed: ' + declaration)
    closure = dict(zip(function.__code__.co_freevars, function.__closure__ or ()))
    if closure:
        if set(closure) != {'__class__'}: raise ValueError('unexpected package closure: ' + declaration)
        parent = qualname.rsplit('.', 1)[0]
        if '<locals>' in parent or closure['__class__'].cell_contents is not slot(name + ':' + parent):
            raise ValueError('package class closure changed: ' + declaration)
    return expected[0]


def wrapper(function, declaration, mode, context):
    if mode == 'operator':
        wrapper_declaration = 'pandas.core.ops.common:_unpack_zerodim_and_defer.<locals>.new_method'
        expected_names = {'is_cmp', 'method'}; body_name = 'method'
    else:
        wrapper_declaration = 'pandas.core.strings.accessor:forbid_nonstring_types.<locals>._forbid_nonstring_types.<locals>.wrapper'
        expected_names = {'allowed_types', 'func', 'func_name'}; body_name = 'func'
    # The wrapper itself has no defaults; its captured function is validated
    # against the PUBLIC SLOT's named source declaration, never its own name.
    name, qualname = wrapper_declaration.split(':'); data = context[name]
    matches = [code for code in codes(data['code']) if code.co_qualname == qualname]
    if (type(function) is not FunctionType or function.__globals__ is not data['module'].__dict__ or len(matches) != 1
            or executable(function.__code__) != executable(matches[0])
            or function.__defaults__ is not None or function.__kwdefaults__ is not None):
        raise ValueError('package decorator code/defaults changed: ' + declaration)
    cells = dict(zip(function.__code__.co_freevars, function.__closure__ or ()))
    if set(cells) != expected_names: raise ValueError('package decorator closure changed')
    inner = cells[body_name].cell_contents
    plain_function(inner, declaration, context)
    if function.__dict__.get('__wrapped__') is not inner: raise ValueError('package wrapped binding changed')
    if mode == 'operator':
        operator_name = declaration.rsplit('.', 1)[-1].removeprefix('__').removesuffix('__')
        if typed(cells['is_cmp'].cell_contents) != typed(operator_name in ('eq', 'ne', 'lt', 'le', 'gt', 'ge')):
            raise ValueError('package operator closure flag changed')
    else:
        node = context[declaration.split(':')[0]]['definitions'][declaration.split(':')[1]][0]
        decorators = [value for value in node.decorator_list if type(value) is ast.Call
            and type(value.func) is ast.Name and value.func.id == 'forbid_nonstring_types']
        if len(decorators) != 1 or len(decorators[0].args) != 1: raise ValueError('string decorator source changed')
        forbidden = ast.literal_eval(decorators[0].args[0])
        if type(forbidden) is not list or any(type(value) is not str for value in forbidden): raise ValueError('string decorator literal changed')
        allowed = {'string', 'empty', 'bytes', 'mixed', 'mixed-integer'} - set(forbidden)
        if typed(cells['allowed_types'].cell_contents) != typed(allowed) or typed(cells['func_name'].cell_contents) != typed(node.name):
            raise ValueError('string decorator values changed')


def _numpy(expected_extension=None):
    from numpy.core import _multiarray_umath as native
    raw = regular_bytes(native.__file__); digest = hashlib.sha256(raw).hexdigest()
    if expected_extension is not None and expected_extension != digest: raise ValueError('NumPy extension changed')
    for name in ('array', 'asarray', 'zeros', 'empty'):
        value, canonical = getattr(np, name), getattr(native, name)
        if (type(value) is not BuiltinFunctionType or value is not canonical
                or value.__name__ != name or value.__module__ != 'numpy'):
            raise ValueError('NumPy C constructor changed: ' + name)
    for name in ('ndarray', 'dtype'):
        value, canonical = getattr(np, name), getattr(native, name)
        if (value is not canonical or type.__getattribute__(value, '__name__') != name
                or type.__getattribute__(value, '__module__') != 'numpy'
                or not type.__getattribute__(value, '__flags__') & (1 << 8)):
            raise ValueError('NumPy immutable C class changed: ' + name)
    for name, dtype in (('bool_', 'bool'), ('int64', 'int64')):
        if getattr(np, name) is not native.dtype(dtype).type: raise ValueError('NumPy scalar class changed')
    # NumPy1.26 exports the ufunc objects here, but not the ufunc TYPE name.
    # Derive that immutable C type only from a type/name-checked native object;
    # a Python replacement cannot supply a C ufunc type or its readonly name.
    canonical_ufunc = type(native.bitwise_and)
    if (np.ufunc is not canonical_ufunc or type.__getattribute__(canonical_ufunc, '__name__') != 'ufunc'
            or type.__getattribute__(canonical_ufunc, '__module__') != 'numpy'
            or not type.__getattribute__(canonical_ufunc, '__flags__') & (1 << 8)):
        raise ValueError('NumPy immutable ufunc class changed')
    for name in ('bitwise_and', 'bitwise_or', 'invert', 'less', 'greater', 'greater_equal', 'equal'):
        value = getattr(np, name)
        if type(value) is not canonical_ufunc or value is not getattr(native, name) or value.__name__ != name:
            raise ValueError('NumPy ufunc changed: ' + name)
    return digest


def builtin_guard(namespace, name):
    value = getattr(builtins, name)
    if name == 'len':
        valid = type(value) is BuiltinFunctionType and value.__module__ == 'builtins' and value.__name__ == 'len'
    else:
        valid = type(value) is type and value.__module__ == 'builtins' and value.__name__ == name and bool(value.__flags__ & (1 << 8))
    if not valid or namespace.get(name, value) is not value: raise ValueError('original builtin binding changed')


def special_bindings(context):
    from pandas.core import accessor, arraylike, frame, generic, roperator, series
    from pandas.core.strings import accessor as string_module
    if pd.DataFrame is not frame.DataFrame or pd.Series is not series.Series: raise ValueError('public pandas class changed')
    if inspect.getattr_static(frame.DataFrame, '_constructor_sliced') is not series.Series:
        raise ValueError('pandas sliced constructor changed')
    from pandas._libs.properties import AxisProperty
    for owner, name, axis in ((frame.DataFrame, 'columns', 0), (frame.DataFrame, 'index', 1), (series.Series, 'index', 0)):
        value = inspect.getattr_static(owner, name)
        if type(value) is not AxisProperty or value.axis != axis: raise ValueError('pandas axis descriptor changed')
    for owner in (frame.DataFrame, series.Series, string_module.StringMethods, logging.Logger, logging.Handler, logging.Formatter):
        if inspect.getattr_static(owner, '__getattribute__') is not object.__getattribute__:
            raise ValueError('custom package attribute access')
    descriptor = inspect.getattr_static(series.Series, 'str')
    if (type(descriptor) is not accessor.CachedAccessor or type(descriptor.__dict__) is not dict
            or set(descriptor.__dict__) != {'_name', '_accessor'}
            or typed(descriptor.__dict__['_name']) != typed('str')
            or descriptor.__dict__['_accessor'] is not string_module.StringMethods):
        raise ValueError('Series.str descriptor/accessor changed')
    import operator
    import _operator
    for name in ('and_', 'or_', 'invert', 'eq', 'gt', 'ge', 'lt', 'add', 'sub'):
        value = getattr(operator, name)
        if (value is not getattr(_operator, name) or type(value) is not BuiltinFunctionType
                or value.__name__ != name or value.__module__ != '_operator'):
            raise ValueError('canonical operator binding changed')
    for module, name, expected in ((arraylike, 'operator', operator), (arraylike, 'roperator', roperator),
                                  (generic, 'operator', operator)):
        if module.__dict__.get(name) is not expected: raise ValueError('package operator global changed')
    for name in ('BooleanArray', 'IntegerArray', 'StringArray', 'NumpyExtensionArray'):
        module_name = {'BooleanArray': 'boolean', 'IntegerArray': 'integer', 'StringArray': 'string_', 'NumpyExtensionArray': 'numpy_'}[name]
        module = importlib.import_module('pandas.core.arrays.' + module_name)
        if getattr(pd.arrays, name) is not getattr(module, name): raise ValueError('public pandas array class changed')
    if logging._logRecordFactory is not logging.LogRecord: raise ValueError('custom log record factory')


def source_names():
    return {declaration.split(':')[0] for _label, declaration, _mode in SPECS} | {
        label.split(':')[0] for label, _declaration, _mode in SPECS} | {
        'pandas.core.ops.common', 'pandas.core.arraylike', 'pandas.core.accessor', 'pandas.core.strings.accessor'}


def validate(context):
    for label, declaration, mode in SPECS:
        value = slot(label)
        if mode == 'property':
            if type(value) is not property: raise ValueError('package property changed: ' + label)
            value = value.fget
        elif mode == 'staticmethod':
            if type(value) is not staticmethod: raise ValueError('package staticmethod changed: ' + label)
            value = value.__func__
        if mode in ('operator', 'contains'): wrapper(value, declaration, mode, context)
        else: plain_function(value, declaration, context, property_getter=mode == 'property')
    special_bindings(context)


def freeze(path):
    if type(np.__version__) is not str or type(pd.__version__) is not str or np.__version__ != '1.26.4' or pd.__version__ != '2.2.3':
        raise ValueError('unqualified package versions')
    context = source_context(source_names()); validate(context)
    value = dict(schema=2, numpy=np.__version__, pandas=pd.__version__, python=list(sys.version_info[:3]),
        source_files={name: item['sha256'] for name, item in context.items()},
        canonical_slot_declarations=[list(row) for row in SPECS], numpy_extension_sha256=_numpy(),
        identity='named canonical source/code/defaults and recognized closure values; no post-import API identity snapshot',
        package_code_comparison='all recursive executable fields excluding only co_linetable debugger positions')
    with Path(path).open('x') as stream: json.dump(value, stream, sort_keys=True, indent=2); stream.write('\n')
    return value


def guard(binding, digest):
    if type(np.__version__) is not str or type(pd.__version__) is not str or np.__version__ != '1.26.4' or pd.__version__ != '2.2.3':
        raise ValueError('unqualified package versions')
    raw = regular_bytes(binding)
    if hashlib.sha256(raw).hexdigest() != digest: raise ValueError('package proof changed')
    trusted = json.loads(raw)
    if (trusted.get('schema') != 2 or trusted['canonical_slot_declarations'] != [list(row) for row in SPECS]
            or trusted['numpy'] != np.__version__ or trusted['pandas'] != pd.__version__
            or trusted['python'] != list(sys.version_info[:3])):
        raise ValueError('package profile changed')
    if pd.options.mode.copy_on_write is not False: raise ValueError('CoW profile requires original')
    try:
        _numpy(trusted['numpy_extension_sha256'])
        context = source_context(source_names(), trusted['source_files']); validate(context)
    except ImportError as error: raise ValueError('qualified package API unavailable') from error


def instance_guard(value, *, series=False):
    # No object method or descriptor is invoked. Cached StringMethods is allowed
    # only with its ordinary constructor state and the actual owning Series.
    data = object.__getattribute__(value, '__dict__')
    names = {label.split(':')[1].rsplit('.', 1)[-1] for label, _declaration, _mode in SPECS
             if label.startswith('pandas.core.series:Series.')}
    names |= {'get', 'assign'}
    if any(type(key) is not str or key in names for key in data): raise ValueError('per-instance pandas method override')
    if series and 'str' in data:
        from pandas.core.strings.accessor import StringMethods
        accessor = data['str']
        if type(accessor) is not StringMethods: raise ValueError('custom cached string accessor')
        fields = object.__getattribute__(accessor, '__dict__')
        expected = {'_data', '_orig', '_parent', '_index', '_name', '_inferred_dtype', '_is_categorical', '_is_string', '__frozen'}
        if set(fields) != expected or any(fields[key] is not value for key in ('_data', '_orig', '_parent')):
            raise ValueError('cached string accessor state changed')
        raw = value._mgr.blocks[0].values
        if (fields['_index'] is not value.index or type(fields['_name']) is not str or fields['_name'] != value.name
                or typed(fields['_inferred_dtype']) != typed('string') or fields['_is_categorical'] is not False
                or fields['_is_string'] is not (type(raw) is pd.arrays.StringArray) or fields['__frozen'] is not True):
            raise ValueError('cached string accessor context changed')


def logger_instance_guard(logger):
    methods = {label.split(':')[1].rsplit('.', 1)[-1] for label, _declaration, _mode in SPECS if label.startswith('logging:Logger.')}
    current = logger
    while current is not None:
        if type(current) not in (logging.Logger, logging.RootLogger): raise ValueError('custom logger')
        data = object.__getattribute__(current, '__dict__')
        if any(name in data for name in methods) or type(data.get('name')) is not str:
            raise ValueError('per-instance logger method/name changed')
        if type(data.get('level')) is not int or type(data.get('disabled')) is not bool:
            raise ValueError('custom logger state')
        handlers = data.get('handlers')
        if type(handlers) is not list: raise ValueError('custom logger handlers')
        shadow = {'handle', 'filter', 'emit', 'format', 'formatMessage', 'formatTime', 'getMessage',
                  'flush', 'acquire', 'release', 'handleError'}
        for handler in handlers:
            if type(handler) not in (logging.StreamHandler, logging.NullHandler): raise ValueError('custom logging handler')
            state = object.__getattribute__(handler, '__dict__')
            if any(key in state for key in shadow): raise ValueError('per-instance handler method override')
            formatter = state.get('formatter')
            if formatter is not None:
                if type(formatter) is not logging.Formatter: raise ValueError('custom formatter')
                fields = object.__getattribute__(formatter, '__dict__')
                if any(key in fields for key in shadow | {'converter'}): raise ValueError('per-instance formatter override')
        if not data.get('propagate'): break
        current = data.get('parent')


def main():
    parser = argparse.ArgumentParser(description=__doc__); commands = parser.add_subparsers(dest='procedure', required=True)
    command = commands.add_parser('freeze'); command.add_argument('--output', type=Path, required=True)
    command = commands.add_parser('check'); command.add_argument('--binding', type=Path, required=True); command.add_argument('--sha256', required=True)
    args = parser.parse_args()
    if args.procedure == 'freeze':
        result = freeze(args.output); print(json.dumps(dict(passed=True, APIs=len(SPECS), source_modules=len(result['source_files'])), sort_keys=True))
    else: guard(args.binding, args.sha256); print(json.dumps(dict(passed=True, APIs=len(SPECS)), sort_keys=True))


if __name__ == '__main__': main()
