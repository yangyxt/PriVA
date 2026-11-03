"""Explicit package admission for the optional frequency procedure.

This module does not freeze or approve a live environment. Native callers must
supply a separately reviewed binding. SciPy's actual binomial methods and the
mature gzip/pickle decoders remain the scientific dependency boundary.
"""
import ast
import copy
import hashlib
import importlib
import inspect
import json
from pathlib import Path
import sys
from types import BuiltinFunctionType, FunctionType, MethodType

import numpy as np
import pandas as pd
import scipy
import scipy.stats

import priva_rust_pm4_bp3_packages as canonical
import priva_rust_frequency_source as source

SCHEMA = 'priva-frequency-mature-packages-v1'
REFERENCE_VERIFICATION = '4f47b8c187b7f28ff1a4f3d7fa58aad90d823dfddb351707dd0f0e00f607aee2'
VERSIONS = {'numpy': '1.26.4', 'pandas': '2.2.3', 'scipy': '1.14.1'}
# Public slots and their named, unmodified source definitions. Additional
# pandas slots are admitted by the existing reusable canonical-package guard.
SPECS = (
    ('scipy.stats._distn_infrastructure:rv_discrete.sf',
     'scipy.stats._distn_infrastructure:rv_discrete.sf'),
    ('scipy.stats._distn_infrastructure:rv_discrete.cdf',
     'scipy.stats._distn_infrastructure:rv_discrete.cdf'),
    ('scipy.stats._discrete_distns:binom_gen._sf',
     'scipy.stats._discrete_distns:binom_gen._sf'),
    ('scipy.stats._discrete_distns:binom_gen._cdf',
     'scipy.stats._discrete_distns:binom_gen._cdf'),
    ('gzip:open', 'gzip:open'),
)
METHOD_SPECS = (
    ('pandas.core.series:Series.__add__', 'pandas.core.arraylike:OpsMixin.__add__', 'operator'),
    ('pandas.core.series:Series.__radd__', 'pandas.core.arraylike:OpsMixin.__radd__', 'operator'),
    ('pandas.core.series:Series.__sub__', 'pandas.core.arraylike:OpsMixin.__sub__', 'operator'),
    ('pandas.core.series:Series.__rsub__', 'pandas.core.arraylike:OpsMixin.__rsub__', 'operator'),
    ('pandas.core.series:Series.__mul__', 'pandas.core.arraylike:OpsMixin.__mul__', 'operator'),
    ('pandas.core.series:Series.__rmul__', 'pandas.core.arraylike:OpsMixin.__rmul__', 'operator'),
    ('pandas.core.series:Series.__truediv__', 'pandas.core.arraylike:OpsMixin.__truediv__', 'operator'),
    ('pandas.core.series:Series.__rtruediv__', 'pandas.core.arraylike:OpsMixin.__rtruediv__', 'operator'),
    ('pandas.core.series:Series.__and__', 'pandas.core.arraylike:OpsMixin.__and__', 'operator'),
    ('pandas.core.series:Series.__rand__', 'pandas.core.arraylike:OpsMixin.__rand__', 'operator'),
    ('pandas.core.series:Series.__or__', 'pandas.core.arraylike:OpsMixin.__or__', 'operator'),
    ('pandas.core.series:Series.__ror__', 'pandas.core.arraylike:OpsMixin.__ror__', 'operator'),
    ('pandas.core.series:Series.__eq__', 'pandas.core.arraylike:OpsMixin.__eq__', 'operator'),
    ('pandas.core.series:Series.__ne__', 'pandas.core.arraylike:OpsMixin.__ne__', 'operator'),
    ('pandas.core.series:Series.__lt__', 'pandas.core.arraylike:OpsMixin.__lt__', 'operator'),
    ('pandas.core.series:Series.__le__', 'pandas.core.arraylike:OpsMixin.__le__', 'operator'),
    ('pandas.core.series:Series.__gt__', 'pandas.core.arraylike:OpsMixin.__gt__', 'operator'),
    ('pandas.core.series:Series.__ge__', 'pandas.core.arraylike:OpsMixin.__ge__', 'operator'),
    ('pandas.core.series:Series.__invert__', 'pandas.core.generic:NDFrame.__invert__', 'plain'),
    ('pandas.core.series:Series.fillna', 'pandas.core.generic:NDFrame.fillna', 'plain'),
    ('pandas.core.series:Series.astype', 'pandas.core.generic:NDFrame.astype', 'plain'),
    ('pandas.core.series:Series.clip', 'pandas.core.generic:NDFrame.clip', 'plain'),
    ('pandas.core.series:Series.mask', 'pandas.core.generic:NDFrame.mask', 'plain'),
    ('pandas.core.series:Series.replace', 'pandas.core.generic:NDFrame.replace', 'plain'),
    ('pandas.core.series:Series.map', 'pandas.core.series:Series.map', 'plain'),
    ('pandas.core.series:Series.groupby', 'pandas.core.series:Series.groupby', 'plain'),
    ('pandas.core.series:Series.isna', 'pandas.core.series:Series.isna', 'plain'),
    ('pandas.core.series:Series.isnull', 'pandas.core.series:Series.isnull', 'plain'),
    ('pandas.core.series:Series.notna', 'pandas.core.series:Series.notna', 'plain'),
    ('pandas.core.series:Series.isin', 'pandas.core.series:Series.isin', 'plain'),
    ('pandas.core.series:Series.sum', 'pandas.core.series:Series.sum', 'plain'),
    ('pandas.core.series:Series.any', 'pandas.core.series:Series.any', 'plain'),
    ('pandas.core.series:Series.max', 'pandas.core.series:Series.max', 'plain'),
    ('pandas.core.series:Series.to_dict', 'pandas.core.series:Series.to_dict', 'nonkeyword'),
    ('pandas.core.series:Series.eq', 'pandas.core.series:Series.eq', 'plain'),
    ('pandas.core.series:Series.ge', 'pandas.core.series:Series.ge', 'plain'),
    ('pandas.core.frame:DataFrame.apply', 'pandas.core.frame:DataFrame.apply', 'plain'),
    ('pandas.core.groupby.generic:SeriesGroupBy.max', 'pandas.core.groupby.groupby:GroupBy.max', 'plain'),
    ('pandas.core.strings.accessor:StringMethods.lower', 'pandas.core.strings.accessor:StringMethods.lower', 'contains'),
    ('pandas.core.strings.accessor:StringMethods.upper', 'pandas.core.strings.accessor:StringMethods.upper', 'contains'),
    ('pandas.core.strings.accessor:StringMethods.removeprefix', 'pandas.core.strings.accessor:StringMethods.removeprefix', 'contains'),
    ('pandas.core.strings.accessor:StringMethods.contains', 'pandas.core.strings.accessor:StringMethods.contains', 'contains'),
)
NDARRAY_METHODS = ('astype', 'clip', 'sum', 'any', 'max')
ROLES = ('numpy', 'pandas', 'scipy', 'scipy_distribution_methods',
         'scipy_discrete_distribution', 'scipy_special_ufuncs',
         'numpy_multiarray', 'pandas_frame', 'pandas_series',
         'pandas_internals', 'pickle', 'pickle_extension', 'gzip')


def _sha(value):
    return type(value) is str and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def _identity(path):
    st = Path(path).lstat()
    return tuple(getattr(st, name) for name in
                 ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns', 'st_mode', 'st_nlink'))


def _role_paths():
    modules = {role: importlib.import_module(name) for role, name in (
        ('numpy', 'numpy'), ('pandas', 'pandas'), ('scipy', 'scipy'),
        ('scipy_distribution_methods', 'scipy.stats._distn_infrastructure'),
        ('scipy_discrete_distribution', 'scipy.stats._discrete_distns'),
        ('scipy_special_ufuncs', 'scipy.special._ufuncs'),
        ('numpy_multiarray', 'numpy.core._multiarray_umath'),
        ('pandas_frame', 'pandas.core.frame'), ('pandas_series', 'pandas.core.series'),
        ('pandas_internals', 'pandas._libs.internals'),
        ('pickle', 'pickle'), ('pickle_extension', '_pickle'), ('gzip', 'gzip'))}
    return {role: str(module.__file__) for role, module in modules.items()}


def _declared_plain_function(function, declaration, context):
    """Resolve only gzip.open's pinned literal named default; never freeze live values."""
    if declaration != 'gzip:open':
        return canonical.plain_function(function, declaration, context)
    data = context['gzip']
    assignments = [node for node in ast.parse(data['raw']).body
                   if type(node) is ast.Assign and len(node.targets) == 1
                   and type(node.targets[0]) is ast.Name
                   and node.targets[0].id == '_COMPRESS_LEVEL_BEST']
    if (len(assignments) != 1 or type(assignments[0].value) is not ast.Constant
            or type(assignments[0].value.value) is not int
            or assignments[0].value.value != 9):
        raise ValueError('unreviewed gzip compression-level source declaration')
    literal = assignments[0].value.value
    live = data['module'].__dict__.get('_COMPRESS_LEVEL_BEST')
    if type(live) is not int or live != literal:
        raise ValueError('gzip compression-level global differs from source')
    nodes = data['definitions'].get('open', [])
    if len(nodes) != 1 or type(nodes[0]) is not ast.FunctionDef:
        raise ValueError('ambiguous gzip.open declaration')
    node = nodes[0]
    if ([item.arg for item in node.args.args] !=
            ['filename', 'mode', 'compresslevel', 'encoding', 'errors', 'newline']
            or len(node.args.defaults) != 5
            or type(node.args.defaults[1]) is not ast.Name
            or node.args.defaults[1].id != '_COMPRESS_LEVEL_BEST'):
        raise ValueError('unreviewed gzip.open named-default position')
    # The original compiled code/raw source/globals remain untouched. Only
    # this local AST view for default validation substitutes the source's
    # exact literal assignment. All canonical code/closure/default checks run.
    resolved = copy.copy(node)
    resolved.args = copy.copy(node.args)
    resolved.args.defaults = list(node.args.defaults)
    resolved.args.defaults[1] = ast.copy_location(ast.Constant(value=literal), node.args.defaults[1])
    definitions = dict(data['definitions'], open=[resolved])
    checked_context = dict(context, gzip=dict(data, definitions=definitions))
    return canonical.plain_function(function, declaration, checked_context)


def _nonkeyword(function, declaration, context):
    """Validate the actual Series.to_dict decorator and captured source body."""
    name = 'pandas.util._decorators'
    qualname = 'deprecate_nonkeyword_arguments.<locals>.decorate.<locals>.wrapper'
    data = context[name]
    expected = [code for code in canonical.codes(data['code']) if code.co_qualname == qualname]
    if (type(function) is not FunctionType or function.__globals__ is not data['module'].__dict__
            or len(expected) != 1
            or canonical.executable(function.__code__) != canonical.executable(expected[0])
            or function.__defaults__ is not None or function.__kwdefaults__ is not None):
        raise ValueError('nonkeyword decorator code/defaults changed')
    cells = dict(zip(function.__code__.co_freevars, function.__closure__ or ()))
    if set(cells) != {'allow_args', 'func', 'msg', 'num_allow_args'}:
        raise ValueError('nonkeyword decorator closure changed')
    inner = cells['func'].cell_contents
    canonical.plain_function(inner, declaration, context)
    if function.__dict__.get('__wrapped__') is not inner:
        raise ValueError('nonkeyword wrapped source changed')
    module, qualname = declaration.split(':')
    nodes = [node for node in context[module]['definitions'][qualname]
             if not any(type(value) is ast.Name and value.id == 'overload' for value in node.decorator_list)]
    decorators = [value for value in nodes[0].decorator_list if type(value) is ast.Call
                  and type(value.func) is ast.Name and value.func.id == 'deprecate_nonkeyword_arguments']
    if len(decorators) != 1 or decorators[0].args:
        raise ValueError('nonkeyword source decorator changed')
    values = {value.arg: ast.literal_eval(value.value) for value in decorators[0].keywords}
    if values != {'version': '3.0', 'allowed_args': ['self'], 'name': 'to_dict'}:
        raise ValueError('unreviewed nonkeyword source decorator arguments')
    expected_message = 'Starting with pandas version 3.0 all arguments of to_dict{arguments} will be keyword-only.'
    for key, value in (('allow_args', ['self']), ('num_allow_args', 1), ('msg', expected_message)):
        if canonical.typed(cells[key].cell_contents) != canonical.typed(value):
            raise ValueError('nonkeyword source closure value changed: ' + key)


class Guard:
    """Use immutable source declarations, never a snapshot of mutable API slots."""
    def __init__(self, binding, digest, reference):
        if not _sha(digest):
            raise ValueError('canonical package binding SHA is required')
        raw = source.regular_bytes(binding)
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('frequency package binding changed')
        self.binding = json.loads(raw)
        value = self.binding
        if (value.get('schema') != SCHEMA
                or value.get('reference_verification_sha256') != REFERENCE_VERIFICATION
                or value.get('versions') != VERSIONS
                or value.get('python') != list(sys.version_info[:3])
                or value.get('canonical_slot_declarations') != [list(row) for row in SPECS]
                or value.get('frequency_method_declarations') != [list(row) for row in METHOD_SPECS]
                or value.get('ndarray_method_names') != list(NDARRAY_METHODS)):
            raise ValueError('frequency mature-package profile changed')
        if {name: module.__version__ for name, module in
                (('numpy', np), ('pandas', pd), ('scipy', scipy))} != VERSIONS:
            raise ValueError('unqualified mature-package version')
        if pd.options.mode.copy_on_write is not False:
            raise ValueError('frequency profile requires original pandas ownership')
        records = value.get('package_files')
        if type(records) is not dict or set(records) != set(ROLES):
            raise ValueError('exact thirteen-file original dependency frontier required')
        paths = _role_paths()
        self.identities = {str(binding): _identity(binding)}
        for role, path in paths.items():
            record = records[role]
            if (record.get('path') != path or not _sha(record.get('sha256'))
                    or type(record.get('bytes')) is not int):
                raise ValueError('package file role/path/pin differs: ' + role)
            content = source.regular_bytes(path)
            if len(content) != record['bytes'] or hashlib.sha256(content).hexdigest() != record['sha256']:
                raise ValueError('package bytes differ: ' + role)
            self.identities[path] = _identity(path)
        # Reuse the separately pinned canonical pandas/NumPy guard. Its binding
        # is an explicit input, not generated or refreshed by this adapter.
        base = value['canonical_pandas_binding']
        canonical.guard(base['path'], base['sha256'])
        self.base = base
        base_raw = source.regular_bytes(base['path'])
        if hashlib.sha256(base_raw).hexdigest() != base['sha256']:
            raise ValueError('canonical pandas binding changed during admission')
        self.identities[base['path']] = _identity(base['path'])
        base_document = json.loads(base_raw)
        self.base_context = canonical.source_context(canonical.source_names(), base_document['source_files'])
        for item in self.base_context.values():
            self.identities[str(item['module'].__file__)] = _identity(item['module'].__file__)
        names = {item.split(':')[0] for row in SPECS for item in row}
        names |= {item.split(':')[0] for row in METHOD_SPECS for item in row[:2]}
        names |= {'pandas.core.tools.numeric', 'pandas.util._decorators', 'pandas.core.ops.common',
                  'numpy.core.fromnumeric', 'numpy.core.numeric', 'numpy.core.multiarray'}
        if set(value.get('source_files', {})) != names:
            raise ValueError('canonical frequency dependency source set differs')
        self.context = canonical.source_context(names, value['source_files'])
        for item in self.context.values():
            self.identities[str(item['module'].__file__)] = _identity(item['module'].__file__)
        self.reference = reference
        self.check()

    def check(self, *, complete=False):
        for path, identity in self.identities.items():
            if _identity(path) != identity:
                raise ValueError('mature package metadata changed during procedure')
        for slot, declaration in SPECS:
            _declared_plain_function(canonical.slot(slot), declaration, self.context)
        for slot, _declaration, _mode in METHOD_SPECS:
            self.declared_method(slot)
        canonical.plain_function(pd.to_numeric, 'pandas.core.tools.numeric:to_numeric', self.context)
        canonical.plain_function(np.full, 'numpy.core.numeric:full', self.context)
        from numpy.core import _multiarray_umath as native
        from numpy.core import fromnumeric, multiarray
        for name in ('minimum', 'maximum', 'isnan', 'logical_not'):
            value = getattr(np, name)
            if type(value) is not np.ufunc or value is not getattr(native, name) or value.__name__ != name:
                raise ValueError('canonical frequency ufunc changed: ' + name)
        dispatcher = native._ArrayFunctionDispatcher
        if not type.__getattribute__(dispatcher, '__flags__') & (1 << 8):
            raise ValueError('mutable NumPy dispatcher class')
        if not type.__getattribute__(np.ndarray, '__flags__') & (1 << 8):
            raise ValueError('mutable NumPy ndarray class')
        if type(np.round) is not dispatcher or np.round is not fromnumeric.round:
            raise ValueError('canonical NumPy round dispatcher changed')
        canonical.plain_function(np.round.__wrapped__, 'numpy.core.fromnumeric:round', self.context)
        if (type(np.where) is not dispatcher or np.where is not multiarray.where
                or np.where.__wrapped__ is not native.where
                or type(native.where) is not BuiltinFunctionType or native.where.__name__ != 'where'):
            raise ValueError('canonical NumPy where dispatcher changed')
        import gzip
        import pickle
        import _pickle
        from scipy.stats._discrete_distns import binom_gen
        if (self.reference.binom is not scipy.stats.binom
                or type(self.reference.binom) is not binom_gen
                or self.reference.pickle is not pickle or self.reference.gzip is not gzip
                or self.reference.np is not np or self.reference.pd is not pd):
            raise ValueError('original scientific dependency binding changed')
        load = pickle.load
        if (load is not _pickle.load or type(load) is not BuiltinFunctionType
                or load.__module__ != '_pickle' or load.__name__ != 'load'
                or load.__self__ is not _pickle):
            raise ValueError('mature pickle decoder binding changed')
        for name in ('sf', 'cdf'):
            method = getattr(self.reference.binom, name)
            if (type(method) is not MethodType or method.__self__ is not self.reference.binom
                    or method.__func__ is not canonical.slot('scipy.stats._distn_infrastructure:rv_discrete.' + name)):
                raise ValueError('mature binomial public method changed')
        if complete:
            # Read-only metadata and the cached source-derived code/default
            # proof suffice between callbacks; do not hash package payloads
            # again at every scientific log line.
            canonical.validate(self.base_context)

    def binomial(self, name, method):
        self.check()
        if (name not in ('sf', 'cdf') or type(method) is not MethodType
                or method.__self__ is not self.reference.binom
                or method.__func__ is not canonical.slot('scipy.stats._distn_infrastructure:rv_discrete.' + name)):
            raise ValueError('unapproved numerical dependency call')
        return method

    def declared_method(self, label, implementation=None):
        """Check an invoked function against its named source declaration.

        Comparing two references to the same mutable live class slot does not
        establish its source. A previously captured bound method keeps its
        lookup-before-arguments semantics; its actual implementation is checked.
        """
        rows = [row for row in METHOD_SPECS if row[0] == label]
        if len(rows) != 1:
            raise ValueError('method has no approved source declaration: ' + label)
        slot, declaration, mode = rows[0]
        value = canonical.slot(slot) if implementation is None else implementation
        if mode in ('contains', 'operator'): canonical.wrapper(value, declaration, mode, self.context)
        elif mode == 'nonkeyword': _nonkeyword(value, declaration, self.context)
        else: canonical.plain_function(value, declaration, self.context)
        return value

    def operator(self, name, *operands):
        operations = {'add': ('__add__', '__radd__'), 'sub': ('__sub__', '__rsub__'),
            'mul': ('__mul__', '__rmul__'), 'div': ('__truediv__', '__rtruediv__'),
            'and': ('__and__', '__rand__'), 'or': ('__or__', '__ror__'),
            'eq': ('__eq__',), 'ne': ('__ne__',), 'lt': ('__lt__', '__gt__'),
            'le': ('__le__', '__ge__'), 'gt': ('__gt__', '__lt__'),
            'ge': ('__ge__', '__le__'), 'invert': ('__invert__',)}
        if name not in operations: raise ValueError('undeclared source operator')
        if any(type(value) is pd.Series for value in operands):
            for method in operations[name]:
                self.declared_method('pandas.core.series:Series.' + method)

    def numeric(self, operation):
        canonical.plain_function(operation, 'pandas.core.tools.numeric:to_numeric', self.context)

    def array_function(self, name, operation):
        # Compare the actual callable already captured by the source AST as
        # well as its canonical slot; restoring a live slot after an argument
        # callback cannot authorize a different captured callable.
        self.check()
        if name not in ('where', 'round', 'minimum', 'maximum', 'isnan', 'logical_not', 'zeros', 'full'):
            raise ValueError('undeclared NumPy source call')
        if operation is not getattr(np, name):
            raise ValueError('captured NumPy callable differs from the source binding')

    def method_is_canonical(self, operation, name):
        owner = getattr(operation, '__self__', None)
        from pandas.core.strings.accessor import StringMethods
        from pandas.core.groupby.generic import SeriesGroupBy
        prefixes = {pd.Series: 'pandas.core.series:Series.',
                    pd.DataFrame: 'pandas.core.frame:DataFrame.',
                    StringMethods: 'pandas.core.strings.accessor:StringMethods.',
                    SeriesGroupBy: 'pandas.core.groupby.generic:SeriesGroupBy.'}
        if type(owner) in prefixes:
            if type(operation) is not MethodType: return False
            if name in vars(owner): return False
            try: self.declared_method(prefixes[type(owner)] + name, operation.__func__)
            except (ValueError, TypeError, KeyError): return False
            return True
        if type(owner) is np.ndarray:
            # These named C slots belong to the pinned immutable ndarray type.
            return (name in NDARRAY_METHODS and type(operation) is BuiltinFunctionType
                    and operation == inspect.getattr_static(np.ndarray, name).__get__(owner, np.ndarray))
        return False
