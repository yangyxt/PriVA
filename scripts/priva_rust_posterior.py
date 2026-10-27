"""Optional default-domain posterior plan; pandas owns both publications.

No caller frame is changed here. Unsupported source, storage, values or native
results return None before the original current-buffer apply/zip expression.
"""
import ast
import builtins
import csv
import hashlib
from functools import lru_cache
import os
from pathlib import Path
import re
import stat
import struct
import subprocess
import tempfile
from types import CodeType, FunctionType

import numpy as np
import pandas as pd

from priva_rust_terminal_order import _plain_metadata


CRITERIA = ('PVS1', 'PS1', 'PS2', 'PS3', 'PS4', 'PM1', 'PM2', 'PM3', 'PM4',
            'PM5', 'PM6', 'PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'BA1', 'BS1',
            'BS2', 'BS3', 'BS4', 'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7')
AUDITS = ('audit_BS1_disease_incidence_gate', 'audit_BS1_gene_pathogenic_AF_gate',
          'audit_BS1_gene_pathogenic_AF_blocked_by_PM2',
          'audit_PM2_backdown_for_BS1_disease_incidence', 'frequency_conflict_PM2_BS1')
_SOURCE_SHA = '64ec17b5dd7d28a2404dedbe6c459709920d93e051a9e4bb6c273a7e7224320d'
_NUMPY_POWER = np.power
_SERIES_TYPE = pd.Series
_BUILTINS = {name: getattr(builtins, name) for name in
             ('isinstance', 'int', 'float', 'sum', 'len', 'round')}
_INTEGER_TYPES = frozenset((bool, int, np.bool_, np.int32, np.int64, np.uint32, np.uint64))
_FLOAT_TYPES = frozenset((float, np.float32, np.float64))
_DTYPES = frozenset(('object', 'bool', 'int32', 'int64', 'uint32', 'uint64',
                     'float32', 'float64', 'boolean', 'Int32', 'Int64',
                     'UInt32', 'UInt64', 'Float32', 'Float64'))
_CLASSES = frozenset(('Pathogenic', 'Likely Pathogenic', 'Uncertain Significance',
                      'Likely Benign', 'Benign'))
_FLOAT64_STRENGTHS = {struct.pack('>d', float(n)): n for n in range(6)}
_FLOAT64_STRENGTHS[bytes.fromhex('8000000000000000')] = 0
_FLOAT32_STRENGTHS = {struct.pack('>f', float(n)): n for n in range(6)}
_FLOAT32_STRENGTHS[bytes.fromhex('80000000')] = 0
_TIMEOUT_SECONDS = 300
_PACKAGE_SOURCE_SHA = {
    'core/frame.py': '9b813159606b62af03bfe76a7089da282d2d6091f308cf2f4d90dd0661a01105',
    'core/series.py': '0ba5405fc2a49ac5bed13e3274832ac3b2e8aee15ae58cb5e2172bbe64d1f0d0',
    'core/apply.py': '4dec5994c0e9f8b9b194245e42e1c1e20fe308cf821ca0fe351f45cbf7793547',
    'core/dtypes/cast.py': '527cbd79a643d631917d69a1ea9ebfeb91d27af98d9b5eaf3f406f4f3bc7e614',
    'core/internals/managers.py': 'b680e0a168699ce262c32b6ac91fd7e40989926a9eb58bc1aba29b311f53ebe5',
    'core/internals/blocks.py': '5bb2188a704537c296e51b74a3b049be6720a604ad41e2d4ab73f1a5512cad8f',
    'util/_decorators.py': '9f53b2291446f9d7024544a6c9e8e9293c9f3ff8aeda4545d1327ff7220991ea',
}


def _source_guard(function):
    """Bind the actual caller function and defaults, not only its disk text."""
    if (type(function) is not FunctionType or function.__closure__ is not None
            or function.__name__ != 'calculate_posterior_probability'):
        raise ValueError('custom posterior callable requires Python')
    namespace = function.__globals__
    if namespace.get(function.__name__) is not function:
        raise ValueError('posterior binding changed')
    path = Path(function.__code__.co_filename)
    before = path.lstat()
    if not stat.S_ISREG(before.st_mode):
        raise ValueError('posterior source is not regular')
    raw = path.read_bytes()
    after = path.lstat()
    fields = ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns', 'st_mode')
    if any(getattr(before, name) != getattr(after, name) for name in fields):
        raise ValueError('posterior source changed while reading')
    text = raw.decode('utf-8')
    nodes = [node for node in ast.parse(text).body
             if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
             and node.name == function.__name__]
    if len(nodes) != 1 or type(nodes[0]) is not ast.FunctionDef or nodes[0].decorator_list:
        raise ValueError('ambiguous posterior definition')
    if hashlib.sha256(ast.get_source_segment(text, nodes[0]).encode()).hexdigest() != _SOURCE_SHA:
        raise ValueError('posterior source changed')
    codes = [code for code in compile(raw, str(path), 'exec').co_consts
             if isinstance(code, CodeType) and code.co_name == function.__name__]
    if len(codes) != 1 or function.__code__ != codes[0]:
        raise ValueError('actual posterior code changed')
    defaults = function.__defaults__
    if (type(defaults) is not tuple or tuple(type(value) for value in defaults) != (float, int, int)
            or defaults != (0.1, 2, 350) or function.__kwdefaults__ is not None):
        raise ValueError('nondefault posterior model requires Python')
    if (namespace.get('pd') is not pd or namespace.get('np') is not np
            or pd.Series is not _SERIES_TYPE or np.power is not _NUMPY_POWER
            or type(_NUMPY_POWER) is not np.ufunc):
        raise ValueError('posterior numerical bindings changed')
    for name, expected in _BUILTINS.items():
        if namespace.get(name, expected) is not expected or function.__builtins__.get(name) is not expected:
            raise ValueError('posterior conversion bindings changed')


def _executable_code(code):
    """Package-code identity, excluding only debugger column/line ranges.

    Installed pandas bytecode can retain different co_linetable encodings from
    its build interpreter. Every executable field and nested constant remains
    protected. The PriVA scalar above still requires full CodeType equality.
    """
    constants = tuple(_executable_code(value) if type(value) is CodeType else value for value in code.co_consts)
    names = ('co_argcount', 'co_posonlyargcount', 'co_kwonlyargcount', 'co_nlocals',
             'co_stacksize', 'co_flags', 'co_code', 'co_names', 'co_varnames',
             'co_filename', 'co_name', 'co_qualname', 'co_firstlineno',
             'co_exceptiontable', 'co_freevars', 'co_cellvars')
    return tuple(getattr(code, name) for name in names), constants


@lru_cache(maxsize=7)
def _compiled_package(raw, filename):
    # Source bytes are still re-read, race-checked and hashed on every request.
    # Only compilation of identical trusted bytes is reused, never frame data.
    return compile(raw, filename, 'exec')


def _package_guard():
    """Bind the qualified apply/row-boxing APIs before accessing a frame."""
    from numpy.core import umath
    from pandas.core import apply, frame, series
    from pandas.core.dtypes import cast
    from pandas.core.internals import blocks, managers
    from pandas._libs import lib
    if np.power is not umath.power:
        raise ValueError('NumPy power ufunc changed')
    if pd.DataFrame is not frame.DataFrame or pd.Series is not series.Series:
        raise ValueError('pandas public classes changed')
    package = Path(pd.__file__).parent
    compiled = {}
    for name, digest in _PACKAGE_SOURCE_SHA.items():
        path = package / name
        before = path.lstat()
        if not stat.S_ISREG(before.st_mode):
            raise ValueError('pandas package source is not regular')
        raw = path.read_bytes()
        after = path.lstat()
        fields = ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns', 'st_mode')
        if (any(getattr(before, key) != getattr(after, key) for key in fields)
                or hashlib.sha256(raw).hexdigest() != digest):
            raise ValueError('pandas package source changed')
        compiled[name] = _compiled_package(raw, str(path))

    def codes(code):
        yield code
        for value in code.co_consts:
            if type(value) is CodeType:
                yield from codes(value)

    def api(function, module, path, qualname, line):
        if type(function) is not FunctionType or function.__globals__ is not module.__dict__:
            raise ValueError('pandas method/global binding changed')
        candidates = [code for code in codes(compiled[path])
                      if code.co_qualname == qualname and code.co_firstlineno == line]
        if len(candidates) != 1 or _executable_code(function.__code__) != _executable_code(candidates[0]):
            raise ValueError('pandas executable method code changed')

    api(pd.DataFrame.apply, frame, 'core/frame.py', 'DataFrame.apply', 10170)
    defaults = pd.DataFrame.apply.__defaults__
    expected = (0, False, None, (), 'compat', 'python', None)
    if (type(defaults) is not tuple or tuple(type(value) for value in defaults) != tuple(type(value) for value in expected)
            or defaults != expected or pd.DataFrame.apply.__kwdefaults__ is not None
            or pd.DataFrame.apply.__closure__ is not None):
        raise ValueError('pandas apply defaults changed')
    from pandas.util import _decorators
    wrapper = pd.Series.to_dict
    api(wrapper, _decorators, 'util/_decorators.py',
        'deprecate_nonkeyword_arguments.<locals>.decorate.<locals>.wrapper', 325)
    closure = dict(zip(wrapper.__code__.co_freevars, wrapper.__closure__ or ()))
    if set(closure) != {'allow_args', 'func', 'msg', 'num_allow_args'}:
        raise ValueError('pandas to_dict wrapper changed')
    values = {name: cell.cell_contents for name, cell in closure.items()}
    if (type(values['allow_args']) is not list or len(values['allow_args']) != 1
            or type(values['allow_args'][0]) is not str or values['allow_args'][0] != 'self'
            or type(values['num_allow_args']) is not int or values['num_allow_args'] != 1
            or type(values['msg']) is not str
            or values['msg'] != 'Starting with pandas version 3.0 all arguments of to_dict{arguments} will be keyword-only.'
            or wrapper.__defaults__ is not None or wrapper.__kwdefaults__ is not None):
        raise ValueError('pandas to_dict wrapper defaults changed')
    original = values['func']
    if wrapper.__dict__.get('__wrapped__') is not original:
        raise ValueError('pandas to_dict wrapped binding changed')
    api(original, series, 'core/series.py', 'Series.to_dict', 2037)
    if (type(original.__defaults__) is not tuple or len(original.__defaults__) != 1
            or original.__defaults__[0] is not dict or original.__kwdefaults__ is not None
            or original.__closure__ is not None):
        raise ValueError('pandas to_dict defaults changed')
    generator = apply.FrameColumnApply.__dict__.get('series_generator')
    if type(generator) is not property:
        raise ValueError('pandas row generator changed')
    api(generator.fget, apply, 'core/apply.py', 'FrameColumnApply.series_generator', 1247)
    api(apply.frame_apply, apply, 'core/apply.py', 'frame_apply', 82)
    api(cast.maybe_box_native, cast, 'core/dtypes/cast.py', 'maybe_box_native', 182)
    api(managers.BlockManager.iget_values, managers, 'core/internals/managers.py', 'BlockManager.iget_values', 1026)
    api(blocks.Block.iget, blocks, 'core/internals/blocks.py', 'Block.iget', 1253)
    api(blocks.ExtensionBlock.iget, blocks, 'core/internals/blocks.py', 'ExtensionBlock.iget', 2372)
    if (series.__dict__.get('maybe_box_native') is not cast.maybe_box_native
            or any(cast.__dict__.get(name) is not getattr(lib, name) for name in ('is_float', 'is_integer', 'is_bool'))):
        raise ValueError('pandas row-boxing globals changed')


def _storage(matrix):
    """Reject unknown raw storage before any extension dtype/value callback."""
    # Private/public pandas APIs are inspected only in the qualified context.
    # Missing APIs also select the original function, before touching a frame.
    from pandas.core.internals.blocks import ExtensionBlock, NumpyBlock
    from pandas.core.internals.managers import BlockManager
    extension_types = (pd.arrays.IntegerArray, pd.arrays.FloatingArray, pd.arrays.BooleanArray)
    numpy_type = pd.arrays.NumpyExtensionArray
    if type(matrix._mgr) is not BlockManager:
        raise ValueError('custom matrix manager requires Python')
    for block in matrix._mgr.blocks:
        if type(block) not in (NumpyBlock, ExtensionBlock):
            raise ValueError('custom matrix block requires Python')
        values = block.values
        if type(values) is np.ndarray:
            continue
        if type(values) in extension_types:
            if (type(values._data) is np.ndarray and type(values._mask) is np.ndarray
                    and values._mask.dtype == np.dtype('bool')):
                continue
        if type(values) is numpy_type and type(values._ndarray) is np.ndarray:
            continue
        raise ValueError('custom backing storage requires Python')
    return extension_types, numpy_type


def _raw_columns(matrix):
    """Read original blocks without item-cache or lazy manager-map writes."""
    from pandas._libs.internals import BlockPlacement
    result = [None] * len(matrix.columns)
    for block in matrix._mgr.blocks:
        placement = block._mgr_locs
        if type(placement) is not BlockPlacement:
            raise ValueError('custom block placement requires Python')
        values = block.values
        if type(values) is np.ndarray:
            if values.ndim != 2 or len(placement) != len(values):
                raise ValueError('unexpected raw block shape')
            selected = ((position, values[offset]) for offset, position in enumerate(placement))
        else:
            if len(placement) != 1:
                raise ValueError('unexpected extension block shape')
            selected = ((placement[0], values),)
        for position, column in selected:
            if not 0 <= position < len(result) or result[position] is not None:
                raise ValueError('ambiguous raw block placement')
            result[position] = column
    if any(value is None for value in result):
        raise ValueError('missing raw column')
    return result


def _strength(value):
    kind = type(value)
    if kind in _INTEGER_TYPES:
        result = int(value)
        if 0 <= result <= 5:
            return result
    elif kind in _FLOAT_TYPES:
        # Inspect encodings before any arithmetic/comparison. In particular a
        # signaling NaN must not trigger an extra warning during admission.
        if kind is np.float32:
            encoded = value.tobytes()
            if np.little_endian:
                encoded = encoded[::-1]
            if encoded in _FLOAT32_STRENGTHS:
                return _FLOAT32_STRENGTHS[encoded]
        else:
            encoded = struct.pack('>d', value)
            if encoded in _FLOAT64_STRENGTHS:
                return _FLOAT64_STRENGTHS[encoded]
    raise ValueError('unsupported, missing or nonintegral strength requires Python')


def _admit(matrix, path):
    if (pd.__version__, np.__version__) != ('2.2.3', '1.26.4'):
        raise ValueError('unqualified pandas/NumPy context')
    _package_guard()
    if type(matrix) is not pd.DataFrame:
        raise ValueError('matrix subclass requires Python')
    extension_types, numpy_type = _storage(matrix)
    if (type(matrix.columns) is not pd.Index or any(type(name) is not str for name in matrix.columns)
            or tuple(matrix.columns) not in (CRITERIA, CRITERIA + ('variant_id',) + AUDITS)):
        raise ValueError('unsupported or duplicate columns require Python')
    if type(matrix.index) not in (pd.Index, pd.RangeIndex) or not len(matrix):
        raise ValueError('empty or custom row axis requires Python')
    if not (_plain_metadata(matrix.attrs) and _plain_metadata(tuple(matrix.index.names))
            and _plain_metadata(tuple(matrix.columns.names))
            and all(_plain_metadata(value) for value in matrix.index)):
        raise ValueError('custom metadata/row labels require Python')
    dtypes = (pd.Int32Dtype, pd.Int64Dtype, pd.UInt32Dtype, pd.UInt64Dtype,
              pd.Float32Dtype, pd.Float64Dtype, pd.BooleanDtype)
    columns = []
    for name, array in zip(matrix.columns, _raw_columns(matrix)):
        if type(array) is np.ndarray:
            values = array
        elif type(array) in extension_types:
            if array._mask.any():
                raise ValueError('missing matrix cell requires original row boxing')
            values = array._data
        elif type(array) is numpy_type and type(array._ndarray) is np.ndarray:
            values = array._ndarray
        else:
            raise ValueError('unsupported series array requires Python')
        dtype = array.dtype
        if not (isinstance(dtype, np.dtype) or type(dtype) in dtypes):
            raise ValueError('custom dtype requires Python')
        if name in CRITERIA:
            if str(dtype) not in _DTYPES:
                raise ValueError('unsupported criterion dtype requires Python')
            columns.append([_strength(value) for value in values])
        elif name == 'variant_id':
            if any(type(value) is not str for value in values):
                raise ValueError('nonstring ignored ID requires original row boxing')
        elif any(type(value) not in _INTEGER_TYPES or int(value) not in (0, 1) for value in values):
            # The caller publishes these flags with .astype(int), then copies
            # their arrays into this matrix. Keep those actual integer0/1
            # values as well as ordinary booleans; never coerce other cells.
            raise ValueError('nonstandard audit value requires original row boxing')
    with path.open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['source_position', *CRITERIA])
        for position in range(len(matrix)):
            writer.writerow([str(position), *(str(column[position]) for column in columns)])
    return len(matrix)


def _gather(output, count):
    with output.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['source_position', 'probability_kind', 'probability_value', 'ACMG_class']:
            raise ValueError('invalid posterior header')
        rows = list(reader)
    if len(rows) != count:
        raise ValueError('posterior row count changed')
    result = []
    for position, row in enumerate(rows):
        if set(row) != set(reader.fieldnames) or row['source_position'] != str(position):
            raise ValueError('posterior order/shape changed')
        if row['ACMG_class'] not in _CLASSES:
            raise ValueError('invalid posterior class')
        if row['probability_kind'] == 'int' and row['probability_value'] == '0':
            probability = 0
        elif row['probability_kind'] == 'float64' and re.fullmatch('[0-9a-f]{16}', row['probability_value']):
            bits = int(row['probability_value'], 16)
            if bits > 0x3ff0000000000000:
                raise ValueError('nonfinite/negative/out-of-range native posterior')
            probability = np.array([bits], dtype=np.uint64).view(np.float64)[0]
        else:
            raise ValueError('invalid native posterior scalar')
        result.append((probability, row['ACMG_class']))
    return tuple(result)


def try_posterior_probability(matrix, binary, *, source_function):
    """Return an atomic default-domain plan, or None without frame mutation."""
    try:
        _source_guard(source_function)
        if type(binary) not in (str, type(Path('.'))):
            return None
        executable = os.path.abspath(os.fspath(binary))
        with tempfile.TemporaryDirectory(prefix='priva-posterior-') as directory:
            root = Path(directory)
            count = _admit(matrix, root / 'input.csv')
            environment = os.environ.copy()
            environment.setdefault('POLARS_MAX_THREADS', '1')
            done = subprocess.run([executable, 'posterior-probability', '--input', str(root / 'input.csv'),
                '--output', str(root / 'output.csv')], capture_output=True, env=environment,
                check=False, timeout=_TIMEOUT_SECONDS)
            if done.returncode:
                return None
            return _gather(root / 'output.csv', count)
    except Exception:
        return None
