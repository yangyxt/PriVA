"""Optional atomic R7 plan on the current post-gene-weight ranking buffers.

The ranker retains original tissue selection, conditional parameter validation,
logs and both ordered assignments. Unsupported input or transport keeps the
original PEXT block on these same buffers, without replaying prior stages.
"""
import copy
import csv
import os
from pathlib import Path
import re
import struct
import subprocess
import tempfile

import numpy as np
import pandas as pd

from priva_rust_terminal_order import _cell, _plain_metadata


_AXIS_TOKEN = 'caller-owned'
_TARGET = 'pext_sort_index'
_TIMEOUT_SECONDS = 300
_MAX_ROWS = 1_000_000
_FIELDS = ['source_position', 'index_kind', 'index_value', 'selected',
           'pext_numeric_bits', 'scaled_pext_bits', 'power_bits', 'product_bits',
           'sum_bits', 'pext_sort_index_bits', 'score_bits']


def _scalar(value):
    """Encode the actual boxed map argument, including NumPy scalar origin."""
    if value is None:
        return ['null', '']
    if value is pd.NA:
        return ['pd_na', '']
    if type(value) is np.bool_:
        return ['numpy_bool', str(int(value))]
    if type(value) is np.float32:
        return ['float32', struct.pack('>f', value).hex()]
    if type(value) in (bool, int, float, str, np.float64) or type(value) in (
            np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64):
        if type(value) is int and not -(2**63) <= value <= 2**64 - 1:
            return ['bigint', str(value)]
        return list(_cell(value))
    raise ValueError('custom/unsupported PEXT scalar requires Python')


def _input(frame, column, tissues, cutoff, floor, shape, mask, root):
    if (pd.__version__, np.__version__) != ('2.2.3', '1.26.4'):
        raise ValueError('unqualified pandas/NumPy context')
    if type(frame) is not pd.DataFrame or type(frame.columns) is not pd.Index:
        raise ValueError('frame/column-axis subclass requires Python')
    if type(frame.index) not in (pd.Index, pd.RangeIndex) or not 0 < len(frame) <= _MAX_ROWS:
        raise ValueError('axis/empty/expression-size route requires Python')
    names = list(frame.columns)
    if type(column) is not str or any(type(name) is not str for name in names) or any(
            names.count(name) != 1 for name in (column, 'sort_index')) or names.count(_TARGET) > 1:
        raise ValueError('ambiguous consumed/assigned labels require Python')
    if not (_plain_metadata(frame.attrs) and _plain_metadata(tuple(frame.index.names))
            and _plain_metadata(tuple(frame.columns.names))
            and all(_plain_metadata(label) for label in frame.index)):
        raise ValueError('custom metadata/axis labels require Python')
    if type(tissues) not in (str, type(None)) or any(type(value) is not float
            for value in (cutoff, floor, shape)) or not all(np.isfinite(value) for value in (cutoff, floor, shape)):
        raise ValueError('unsupported option scalar requires Python')
    if cutoff <= 0 or not 0 <= floor <= 1 or shape <= 0:
        raise ValueError('original parameter validation requires Python')
    if type(mask) is not pd.Series or mask.index is not frame.index or not (
            _plain_metadata(mask.name) and _plain_metadata(mask.attrs)):
        raise ValueError('unaligned/custom R4 PEXT mask requires Python')
    array = mask.array
    if type(array) not in (pd.arrays.BooleanArray, pd.arrays.NumpyExtensionArray):
        raise ValueError('custom R4 PEXT mask storage requires Python')
    dtype = mask.dtype
    if type(dtype) is pd.BooleanDtype and type(array) is pd.arrays.BooleanArray:
        mask_dtype = 'boolean'
        data, nulls = array._data, array._mask
        if type(data) is not np.ndarray or type(nulls) is not np.ndarray or (
                data.dtype != np.dtype('bool') or nulls.dtype != np.dtype('bool')):
            raise ValueError('raw Boolean mask storage requires Python')
    elif type(dtype) is type(np.dtype('bool')) and dtype == np.dtype('bool'):
        mask_dtype = 'bool'
        data, nulls = mask.to_numpy(copy=False), np.zeros(len(mask), dtype=bool)
    else:
        raise ValueError('R4 PEXT mask dtype requires Python')
    score = frame.iloc[:, names.index('sort_index')]
    if type(score.array) is not pd.arrays.NumpyExtensionArray or (
            type(score.dtype) is not type(np.dtype('float64')) or score.dtype != np.dtype('float64')):
        raise ValueError('native float64 scores required')
    source = frame[column]
    if type(source.array) not in (pd.arrays.NumpyExtensionArray, pd.arrays.StringArray,
            pd.arrays.BooleanArray, pd.arrays.FloatingArray, pd.arrays.IntegerArray):
        raise ValueError('custom PEXT storage requires Python')
    allowed_dtypes = tuple(type(np.dtype(value)) for value in (
        'object', 'bool', 'int8', 'int16', 'int32', 'int64', 'uint8', 'uint16', 'uint32', 'uint64', 'float32', 'float64'))
    if type(source.dtype) not in allowed_dtypes + (
            pd.StringDtype, pd.BooleanDtype, pd.Float32Dtype, pd.Float64Dtype,
            pd.Int8Dtype, pd.Int16Dtype, pd.Int32Dtype, pd.Int64Dtype,
            pd.UInt8Dtype, pd.UInt16Dtype, pd.UInt32Dtype, pd.UInt64Dtype):
        raise ValueError('custom PEXT dtype requires Python')
    if type(source.dtype) is pd.StringDtype and source.dtype.storage != 'python':
        raise ValueError('unqualified string storage requires Python')
    if not _plain_metadata(source.attrs):
        raise ValueError('custom source-Series metadata requires Python')
    # Reject objects before map can inspect them. Then use the original public
    # Series.map boxing, rather than guessing numeric/nullable scalar origins.
    for value in source.array:
        _scalar(value)
    boxed = []
    source.map(lambda value: boxed.append(_scalar(value)))
    if len(boxed) != len(frame):
        raise ValueError('non-rowwise source map requires Python')
    score_values = score.to_numpy(copy=False)
    with (root / 'input.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value'] +
                        [part for name in names for part in (name + '_kind', name)])
        for row in range(len(frame)):
            writer.writerow(['opaque', str(row)] + [part for name in names
                for part in (_cell(score_values[row]) if name == 'sort_index' else ('null', ''))])
    with (root / 'options.csv').open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['name', 'kind', 'value'])
        for name, value in zip(('tissues', 'cutoff', 'floor', 'shape'), (tissues, cutoff, floor, shape)):
            writer.writerow([name] + list(_cell(value)))
    with (root / 'coercion.csv').open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['source_position', 'value_kind', 'value'])
        for row, value in enumerate(boxed):
            writer.writerow([row] + value)
    with (root / 'mask.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value', 'value_kind', 'value', 'raw_data', 'is_null'])
        for row in range(len(frame)):
            writer.writerow(['opaque', str(row), 'pd_na' if nulls[row] else 'bool',
                '' if nulls[row] else str(int(data[row])), str(int(data[row])) if mask_dtype == 'boolean' else '',
                str(int(nulls[row]))])
    return names, source, mask_dtype, data & ~nulls, score_values


def _series(values, source, name):
    result = pd.Series(values, index=source.index, name=name, dtype='float64')
    result.attrs = copy.deepcopy(source.attrs)
    result.flags.allows_duplicate_labels = source.flags.allows_duplicate_labels
    return result


def _plan(frame, column, tissues, names, source, selected, scores, root):
    with (root / 'output.csv').open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != _FIELDS:
            raise ValueError('invalid native PEXT header')
        rows = list(reader)
    if len(rows) != len(frame) or any(set(row) != set(_FIELDS)
            or row['source_position'] != str(i) or row['index_kind'] != 'opaque'
            or row['index_value'] != str(i) or row['selected'] != str(int(selected[i]))
            or any(re.fullmatch('[0-9a-f]{16}', row[name]) is None for name in _FIELDS[4:])
            for i, row in enumerate(rows)):
        raise ValueError('native PEXT row/mask/bits identity changed')
    tissue_list = [value.strip() for value in tissues.split(',') if value.strip()] if tissues else None
    first = f"PEXT_{tissue_list[0].replace('-', '_')}" if tissue_list else None
    used_tissue = first is not None and first in names
    expected_column = first if used_tissue else ('PEXT_MEAN' if 'PEXT_MEAN' in names else None)
    if column != expected_column:
        raise ValueError('native original tissue selection differs')
    expected = [('axis', _AXIS_TOKEN), ('applied', 'true'), ('column', column),
        ('used_tissue', str(used_tissue).lower()), ('tissue_list_present', str(tissue_list is not None).lower()),
        ('tissue_count', str(len(tissue_list or []))), ('first_tissue_column_present', str(first is not None).lower()),
        ('first_tissue_column', first or ''), ('source_column_position', str(names.index(column))),
        ('target_column_position', str(names.index(_TARGET) if _TARGET in names else len(names))),
        ('score_column_position', str(names.index('sort_index')))]
    expected += [(f'tissue.{i}', value) for i, value in enumerate(tissue_list or [])]
    with (root / 'metadata.csv').open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['name', 'value'] or list(reader) != [dict(name=k, value=v) for k, v in expected]:
            raise ValueError('native complete ordered PEXT metadata changed')
    bits = {name: np.asarray([int(row[name], 16) for row in rows], dtype=np.uint64) for name in _FIELDS[4:]}
    if not np.array_equal(bits['score_bits'][~selected], scores.view(np.uint64)[~selected]):
        raise ValueError('native PEXT changed an unselected score')
    return dict(pext_numeric=_series(bits['pext_numeric_bits'].view(np.float64), source, source.name),
        scaled_pext=_series(bits['scaled_pext_bits'].view(np.float64), source, source.name),
        pext_sort_index=bits['pext_sort_index_bits'].view(np.float64),
        score_values=bits['score_bits'].view(np.float64), selected=selected)


def pext_score_product(frame, mask, plan):
    """Restore original loc operand metadata after the first target assignment."""
    left = frame.loc[mask, 'sort_index']
    right = frame.loc[mask, 'pext_sort_index']
    return _series(plan['score_values'][plan['selected']], left, None)


def try_pext(frame, column, tissues, cutoff, floor, shape, mask, binary):
    """Return one fully validated R7 plan, or None without live publication."""
    try:
        if type(binary) not in (str, type(Path())):
            return None
        executable = os.path.abspath(os.fspath(binary))
        with tempfile.TemporaryDirectory(prefix='priva-r7-') as directory:
            root = Path(directory)
            names, source, dtype, selected, scores = _input(
                frame, column, tissues, cutoff, floor, shape, mask, root)
            command = [executable, 'ranking-pext', '--input', str(root / 'input.csv'),
                '--options-input', str(root / 'options.csv'), '--coercion-input', str(root / 'coercion.csv'),
                '--coercion-column-position', str(names.index(column)), '--mask-input', str(root / 'mask.csv'),
                '--mask-dtype', dtype, '--row-index-kind', 'flat', '--axis-metadata', _AXIS_TOKEN,
                '--output', str(root / 'output.csv'), '--metadata-output', str(root / 'metadata.csv')]
            for position, name in enumerate(names):
                command += ['--column-dtype', f'{position}=' + ('float64' if name == 'sort_index' else 'object')]
            env = os.environ.copy(); env.setdefault('POLARS_MAX_THREADS', '1')
            process = subprocess.run(command, capture_output=True, env=env, check=False, timeout=_TIMEOUT_SECONDS)
            if process.returncode != 0:
                return None
            return _plan(frame, column, tissues, names, source, selected, scores, root)
    except (MemoryError, RecursionError):
        raise
    except Exception:
        return None
