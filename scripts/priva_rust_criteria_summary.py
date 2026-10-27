"""Optional native criteria summaries, with original pandas publication seams.

The caller has constructed/validated its matrix and loaded the selected ClinGen
map. Only successful primitive summary/replacement plans cross this boundary.
No caller frame is changed here. Unsupported requests retain the original
Python block, including partial mutations and failures, at its current position.
"""
import ast
import copy
import csv
import hashlib
import os
from pathlib import Path
import subprocess
import stat
import tempfile
from types import CodeType, FunctionType, SimpleNamespace

import numpy as np
import pandas as pd

from priva_rust_terminal_order import _plain_metadata


CRITERIA = ('PVS1', 'PS1', 'PS2', 'PS3', 'PS4', 'PM1', 'PM2', 'PM3', 'PM4',
            'PM5', 'PM6', 'PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'BA1', 'BS1',
            'BS2', 'BS3', 'BS4', 'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7')
_DTYPES = frozenset(('object', 'int64', 'int32', 'uint64', 'uint32', 'float64',
                     'float32', 'bool', 'Int64', 'Int32', 'UInt64', 'UInt32',
                     'Float64', 'Float32', 'boolean'))
# Unsupported pandas APIs must retain import-time fallback. In the qualified
# context these are the same original class tuple expressions and values.
_EXTENSION_DTYPES = _ARRAY_TYPES = ()
if (pd.__version__, np.__version__) == ('2.2.3', '1.26.4'):
    try:
        _EXTENSION_DTYPES = (pd.Int64Dtype, pd.Int32Dtype, pd.UInt64Dtype, pd.UInt32Dtype,
                     pd.Float64Dtype, pd.Float32Dtype, pd.BooleanDtype)
        _ARRAY_TYPES = (pd.arrays.NumpyExtensionArray, pd.arrays.IntegerArray,
                pd.arrays.FloatingArray, pd.arrays.BooleanArray)
    except AttributeError:
        _EXTENSION_DTYPES = _ARRAY_TYPES = ()
_NUMERIC_TYPES = frozenset((bool, int, float, np.bool_, np.int32, np.int64,
                            np.uint32, np.uint64, np.float32, np.float64))
_SILENT_MAP_TYPES = frozenset((type(None), bool, int, float, complex, bytes,
                              list, tuple, dict, set, frozenset)) | _NUMERIC_TYPES
_TIMEOUT_SECONDS = 300
_FUNCTION_PINS = {
    'summarize_acmg_criteria': '7d29acab47945ec8970de94727dad060ff00d8ad1da5d8110a165739f8f401fd',
    'create_criteria_summary': 'f9e8e5a4ee6e9c8f73577e71e3a00872263e566684a50381b38397536351df24',
    'fix_clingen_term': 'cc637a3aa3cd9aa3b93ddeb1aed59fcd5ed726d263f35eac142e9ae80a4e61dc',
    'translate_strength_level': '80f7e9ca08b66eff8440862d41880167de76a5dc21046511a5495888e95cf5f7',
}


def _source_guard(functions):
    """Bind live callables and their publication seam, not just disk source.

    Compile (never execute) the pinned source in this interpreter and compare
    actual code objects. This catches monkeypatched code even when its filename
    and on-disk source still look original; no cross-version AST dump hash.
    """
    if type(functions) is not tuple or len(functions) != len(_FUNCTION_PINS):
        raise ValueError('missing actual summary function bindings')
    if any(type(f) is not FunctionType or f.__closure__ is not None for f in functions):
        raise ValueError('custom summary callable requires Python')
    namespace = functions[0].__globals__
    if any(f.__globals__ is not namespace or f.__name__ != name or namespace.get(name) is not f
           for name, f in zip(_FUNCTION_PINS, functions)):
        raise ValueError('summary helper bindings changed')
    filename = functions[0].__code__.co_filename
    if any(f.__code__.co_filename != filename for f in functions):
        raise ValueError('summary helpers have different sources')
    path = Path(filename)
    before = path.lstat()
    if not stat.S_ISREG(before.st_mode):
        raise ValueError('summary source is not a regular file')
    raw = path.read_bytes()
    after = path.lstat()
    fields = ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns', 'st_mode')
    if any(getattr(before, k) != getattr(after, k) for k in fields):
        raise ValueError('summary source changed while reading')
    text = raw.decode('utf-8'); tree = ast.parse(text)
    compiled = compile(raw, filename, 'exec')
    for name, function in zip(_FUNCTION_PINS, functions):
        nodes = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name == name]
        if len(nodes) != 1 or type(nodes[0]) is not ast.FunctionDef or nodes[0].decorator_list:
            raise ValueError('ambiguous summary helper definition')
        if hashlib.sha256(ast.get_source_segment(text, nodes[0]).encode()).hexdigest() != _FUNCTION_PINS[name]:
            raise ValueError('summary helper/publication source changed')
        codes = [c for c in compiled.co_consts if isinstance(c, CodeType) and c.co_name == name]
        if len(codes) != 1 or function.__code__ != codes[0]:
            raise ValueError('actual summary helper code changed')
        defaults = ((True,) if name == 'summarize_acmg_criteria' else
                    (None, None) if name == 'create_criteria_summary' else None)
        actual = function.__defaults__
        if defaults is None:
            if actual is not None: raise ValueError('summary helper defaults changed')
        elif (type(actual) is not tuple or len(actual) != len(defaults)
              or any(a is not b for a, b in zip(actual, defaults))):
            raise ValueError('summary helper defaults changed')
        kw = function.__kwdefaults__
        if name == 'summarize_acmg_criteria':
            if type(kw) is not dict or list(kw) != ['summary_binary'] or kw['summary_binary'] is not None:
                raise ValueError('summary keyword defaults changed')
        elif kw is not None:
            raise ValueError('summary helper keyword defaults changed')
    if namespace.get('pd') is not pd or namespace.get('np') is not np:
        raise ValueError('summary numerical namespace changed')
    for name, expected in [('str', str), ('int', int), ('bool', bool), ('isinstance', isinstance)]:
        if namespace.get(name, expected) is not expected:
            raise ValueError('summary conversion namespace changed')


def _ordinary_array(value):
    """Inspect storage types before any extension dtype or iterator callback."""
    kind = type(value)
    if kind is np.ndarray:
        return
    if kind in (pd.arrays.IntegerArray, pd.arrays.FloatingArray, pd.arrays.BooleanArray):
        if type(value._data) is np.ndarray and type(value._mask) is np.ndarray:
            return
    elif kind in (pd.arrays.NumpyExtensionArray, pd.arrays.StringArray):
        if type(value._ndarray) is np.ndarray:
            return
    raise ValueError('custom backing storage requires Python')


def _ordinary_matrix_storage(matrix):
    # Only reached in the qualified pandas context. Raw block values are stored
    # fields; Series.dtype would call an unknown extension array's dtype getter
    # before the later Series-array type gate could reject it.
    from pandas.core.internals.blocks import ExtensionBlock, NumpyBlock
    from pandas.core.internals.managers import BlockManager
    if type(matrix._mgr) is not BlockManager:
        raise ValueError('custom matrix manager requires Python')
    for block in matrix._mgr.blocks:
        if type(block) not in (NumpyBlock, ExtensionBlock):
            raise ValueError('custom matrix block requires Python')
        _ordinary_array(block.values)


def _admit(frame, matrix, clingen_map, path):
    if (pd.__version__, np.__version__) != ('2.2.3', '1.26.4'):
        raise ValueError('unqualified pandas/NumPy context')
    if type(frame) is not pd.DataFrame or type(matrix) is not pd.DataFrame:
        raise ValueError('frame subclass requires Python')
    _ordinary_matrix_storage(matrix)
    for current in (frame, matrix):
        if type(current.columns) is not pd.Index or any(type(k) is not str for k in current.columns):
            raise ValueError('custom column axis requires Python')
        if not current.columns.is_unique or type(current.index) not in (pd.Index, pd.RangeIndex):
            raise ValueError('duplicate columns or unsupported index require Python')
        if not (_plain_metadata(current.attrs) and
                _plain_metadata(tuple(current.index.names)) and
                _plain_metadata(tuple(current.columns.names)) and
                all(_plain_metadata(k) for k in current.index)):
            raise ValueError('custom metadata/labels require Python')
        if not current.index.is_unique:
            raise ValueError('duplicate-row alignment requires Python')
    if not len(matrix) or not frame.index.equals(matrix.index):
        raise ValueError('empty frame or unequal alignment requires Python')
    if list(matrix.columns) != list(CRITERIA) + ['variant_id']:
        raise ValueError('unexpected matrix columns require Python')
    if type(clingen_map) is not dict or any(type(k) is not str for k in clingen_map):
        raise ValueError('custom/nonstring map keys require Python')
    identifiers = list(matrix['variant_id'].array)
    if any(type(v) is not str for v in identifiers):
        # apply(axis=1) can box numeric IDs differently from a column iterator.
        # Keep that lookup and any resulting failure in the original function.
        raise ValueError('nonstring variant IDs require original row boxing')
    columns = []
    for name in CRITERIA:
        series = matrix[name]
        dtype = series.dtype
        if not (isinstance(dtype, np.dtype) or type(dtype) in _EXTENSION_DTYPES):
            raise ValueError('custom extension dtype requires Python')
        if str(dtype) not in _DTYPES or type(series.array) not in _ARRAY_TYPES:
            raise ValueError('unsupported criterion representation requires Python')
        if clingen_map and str(dtype) in ('bool', 'boolean'):
            # Original ordered .loc publication owns bool promotion/warnings
            # and nullable-boolean assignment failures.
            raise ValueError('Boolean replacement publication requires Python')
        values = list(series.array)
        if any(type(v) not in _NUMERIC_TYPES for v in values):
            raise ValueError('custom or missing strength requires Python')
        if any(v not in (0, 1, 2, 3, 4, 5) for v in values):
            raise ValueError('source strength validation must be retained')
        columns.append([int(v) for v in values])
    evidence = []
    for identifier in identifiers:
        value = clingen_map.get(identifier) if clingen_map else None
        if type(value) is not str and type(value) not in _SILENT_MAP_TYPES and value is not pd.NA:
            # isinstance can call a custom __class__ getter. Let the original
            # row/mask operations own every such callback, warning and failure.
            raise ValueError('custom map value requires original Python callbacks')
        # Ordinary nonstring values, including real map NaNs, are silence.
        evidence.append(value if type(value) is str else None)
    with path.open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['source_position', *CRITERIA, 'evidence_kind', 'evidence'])
        for i, value in enumerate(evidence):
            writer.writerow([str(i), *(c[i] for c in columns),
                             'none' if value is None else 'str', '' if value is None else value])
    return columns, [value is not None and bool(value.strip()) for value in evidence]


def _metadata(series, source):
    series.attrs = copy.deepcopy(source.attrs)
    series.flags.allows_duplicate_labels = source.flags.allows_duplicate_labels
    return series


def _gather(frame, matrix, columns, expected_override, output):
    with output.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['source_position', 'summary_kind', 'summary', 'overridden', *CRITERIA]:
            raise ValueError('invalid native summary header')
        rows = list(reader)
    if len(rows) != len(matrix):
        raise ValueError('native summary row count changed')
    summaries, mask = [], []
    replacements = {name: [] for name in CRITERIA}
    for i, row in enumerate(rows):
        if set(row) != set(reader.fieldnames) or row['source_position'] != str(i):
            raise ValueError('native summary source order changed')
        if row['overridden'] not in ('0', '1'):
            raise ValueError('invalid native override mask')
        overridden = row['overridden'] == '1'
        if overridden != expected_override[i]:
            raise ValueError('native plan changed the pre-resolved override mask')
        mask.append(overridden)
        if row['summary_kind'] == 'nan' and row['summary'] == '':
            summaries.append(np.nan)
        elif row['summary_kind'] == 'str' and row['summary']:
            summaries.append(row['summary'])
        else:
            raise ValueError('invalid native summary scalar')
        for column, name in enumerate(CRITERIA):
            if row[name] not in ('0', '1', '2', '3', '4', '5'):
                raise ValueError('invalid native replacement strength')
            value = int(row[name])
            if not overridden and value != columns[column][i]:
                raise ValueError('native plan changed an uncurated strength')
            if overridden:
                replacements[name].append(value)
    summary = _metadata(pd.Series(summaries, index=matrix.index), matrix)
    overridden = _metadata(pd.Series(mask, index=matrix.index, name='variant_id', dtype=bool), matrix)
    selected_index = frame.index[np.asarray(mask, dtype=bool)]
    # map on an empty Series retains its input dtype; a nonempty integer map
    # returns int64. Preserve that host seam, including assignment promotion.
    dtype = 'int64' if len(selected_index) else summary.dtype
    updates = {name: _metadata(pd.Series(values, index=selected_index,
                 name='ACMG_criteria', dtype=dtype), frame)
               for name, values in replacements.items()}
    return SimpleNamespace(summary=summary, overridden=overridden, updates=updates)


def try_criteria_summary(frame, matrix, clingen_map, binary, *, source_functions=None):
    """Return a complete native plan, or None before changing either frame."""
    try:
        if source_functions is None:
            import acmg_scoring
            source_functions = tuple(getattr(acmg_scoring, name) for name in _FUNCTION_PINS)
        _source_guard(source_functions)
        executable = os.path.abspath(os.fspath(binary))
        with tempfile.TemporaryDirectory(prefix='priva-criteria-summary-') as directory:
            root = Path(directory)
            columns, expected_override = _admit(frame, matrix, clingen_map, root / 'input.csv')
            env = os.environ.copy()
            env.setdefault('POLARS_MAX_THREADS', '1')
            process = subprocess.run([executable, 'criteria-summary', '--input', str(root / 'input.csv'),
                                      '--output', str(root / 'output.csv')],
                                     capture_output=True, env=env, check=False, timeout=_TIMEOUT_SECONDS)
            if process.returncode:
                return None
            return _gather(frame, matrix, columns, expected_override, root / 'output.csv')
    except Exception:
        return None
