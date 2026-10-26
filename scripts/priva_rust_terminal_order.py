"""Optional adapter for PriVA's existing Rust terminal-order procedure.

The caller has already executed the original control calculations and ranking
weights. Only the final group maximum, row ordering and variant rank are native.
Payload buffers and axes stay in pandas; only consumed primitive cells and
opaque row positions cross the boundary. Unsupported inputs and failed native
requests return None before the current frame is changed, so the caller can
execute its original Python statements on that same frame.
"""
import csv
import os
from pathlib import Path
import struct
import subprocess
import tempfile

import numpy as np
import pandas as pd


_KEYS = frozenset(('chrom', 'pos', 'ref', 'alt', 'sort_index',
                   'ACMG_quant_score', 'variant_id'))
_DTYPES = frozenset(('object', 'string', 'float64', 'Float64', 'int64',
                     'Int64', 'uint64', 'UInt64', 'bool'))
_AXIS_TOKEN = 'caller-owned'
# Operational bound only: expiry retains the original current-frame Python tail.
_TIMEOUT_SECONDS = 300
_NUMPY_SCALARS = frozenset(getattr(np, 'sctypeDict', {}).values())


def _plain_metadata(value, parents=frozenset()):
    """Exclude user callbacks during pandas' private metadata copies."""
    if value is None or value is pd.NA or type(value) in (str, bytes, bool, int, float):
        return True
    if type(value) in _NUMPY_SCALARS:
        return not value.dtype.hasobject
    if type(value) not in (dict, list, tuple) or id(value) in parents:
        return False
    parents = parents | {id(value)}
    values = [part for item in value.items() for part in item] if type(value) is dict else value
    return all(_plain_metadata(part, parents) for part in values)


def _cell(value):
    if value is pd.NA:
        return 'pd_na', ''
    if value is None:
        return 'null', ''
    if type(value) not in (bool, int, float, str) and type(value) not in _NUMPY_SCALARS:
        raise TypeError('custom consumed scalar requires Python')
    if isinstance(value, (bool, np.bool_)):
        return 'bool', str(int(value))
    if isinstance(value, (int, np.integer)):
        if not -(2**63) <= value < 2**64:
            raise ValueError('integer outside the typed transport')
        return ('uint' if value > 2**63 - 1 else 'int'), str(value)
    if isinstance(value, (float, np.floating)):
        if isinstance(value, np.floating) and value.dtype.itemsize > 8:
            raise ValueError('extended float requires Python')
        return 'float', struct.pack('>d', value).hex()
    if isinstance(value, str):
        return 'str', value.encode('utf-8').hex()
    raise TypeError('nonprimitive consumed cell requires Python')


def _input(frame, path):
    # This is the pandas/NumPy context qualified for the native R8 boundary.
    # Other versions retain their own original Python behavior automatically.
    if (pd.__version__, np.__version__) != ('2.2.3', '1.26.4'):
        raise ValueError('unqualified pandas/NumPy context')
    if type(frame) is not pd.DataFrame or type(frame.columns) is not pd.Index:
        raise ValueError('frame/column-axis subclass requires Python')
    names = list(frame.columns)
    if any(type(name) is not str for name in names):
        raise ValueError('nonstring column labels require Python')
    if not frame.columns.is_unique:
        raise ValueError('duplicate-column input requires Python')
    if type(frame.index) not in (pd.Index, pd.RangeIndex, pd.MultiIndex) or not (
            _plain_metadata(frame.attrs) and _plain_metadata(tuple(frame.index.names)) and
            _plain_metadata(tuple(frame.columns.names)) and
            all(_plain_metadata(label) for label in frame.index)):
        raise ValueError('custom metadata/axis labels require Python')
    if not {'chrom', 'pos', 'ref', 'alt', 'sort_index', 'ACMG_quant_score'} <= set(names):
        raise ValueError('original source must handle missing columns')
    dtypes, columns = [], []
    for position, name in enumerate(names):
        series = frame.iloc[:, position]
        dtype = str(series.dtype)
        dtypes.append(dtype)
        if name not in _KEYS:
            # The native algorithm never consumes payload values. Its metadata
            # plan still carries their original labels/dtypes; pandas gathers
            # the real values, including arbitrary Python objects, afterwards.
            columns.append(None)
            continue
        if dtype not in _DTYPES:
            raise ValueError('consumed dtype requires Python')
        if isinstance(series.dtype, pd.StringDtype) and series.dtype.storage != 'python':
            raise ValueError('Arrow string storage requires Python')
        columns.append([_cell(value) for value in series.array])
    with path.open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value'] +
                        [part for name in names for part in (name + '_kind', name)])
        for row in range(len(frame)):
            writer.writerow(['opaque', str(row)] + [part for column in columns
                for part in (('null', '') if column is None else column[row])])
    return names, dtypes


def _number(kind, value):
    if kind == 'null' and value == '':
        return pd.NA
    if kind == 'float' and len(value) == 16:
        return struct.unpack('>d', bytes.fromhex(value))[0]
    if kind in ('int', 'uint'):
        result = int(value)
        if (kind == 'int' and -(2**63) <= result < 2**63) or (
                kind == 'uint' and 0 <= result < 2**64):
            return result
    raise ValueError('invalid native maximum')


def _gather(frame, output, metadata, names, dtypes):
    with output.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['source_position', 'index_kind', 'index_value',
                                 'max_kind', 'max_value', 'variant_rank']:
            raise ValueError('invalid native row header')
        rows = list(reader)
    positions = [int(row['source_position']) for row in rows]
    ranks = [int(row['variant_rank']) for row in rows]
    if len(rows) != len(frame) or sorted(positions) != list(range(len(frame))):
        raise ValueError('native result is not a complete row permutation')
    if any(row['index_kind'] != 'opaque' or row['index_value'] != str(position)
           for row, position in zip(rows, positions)):
        raise ValueError('native row identity changed')
    if any(rank < 0 or rank > len(frame) for rank in ranks):
        raise ValueError('invalid native rank')
    with metadata.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['record', 'position', 'name', 'dtype']:
            raise ValueError('invalid native metadata header')
        records = list(reader)
    expected = [(name, 'int64' if name == 'variant_rank' else dtype)
                for name, dtype in zip(names, dtypes) if name != 'max_variant_score']
    if 'variant_rank' not in names:
        expected.append(('variant_rank', 'int64'))
    wanted = [dict(record='column', position=str(i), name=name, dtype=dtype)
              for i, (name, dtype) in enumerate(expected)]
    temporary_position = names.index('max_variant_score') if 'max_variant_score' in names else len(names)
    if len(records) != len(wanted) + 2 or records[:len(wanted)] != wanted:
        raise ValueError('native output column plan changed')
    temporary, axis = records[-2:]
    if (temporary['record'], temporary['position'], temporary['name']) != (
            'temporary', str(temporary_position), 'max_variant_score'):
        raise ValueError('native temporary column plan changed')
    if axis != dict(record='axis', position='', name=_AXIS_TOKEN, dtype='opaque'):
        raise ValueError('native axis token changed')
    maxima = [None] * len(frame)
    for row, position in zip(rows, positions):
        maxima[position] = _number(row['max_kind'], row['max_value'])
    # Rehearse the exact source mutation order on a private copy. A failure here
    # also leaves the live R8-entry frame available to the original Python tail.
    maximum_dtype = temporary['dtype']
    maxima = [np.nan if value is pd.NA and maximum_dtype == 'float64' else value
              for value in maxima]
    working = frame.copy(deep=True)
    working['max_variant_score'] = pd.Series(maxima, index=frame.index, dtype=maximum_dtype)
    # pandas sort_values keeps an empty RangeIndex; iloc[[]] would expand it.
    gathered = working.copy(deep=True) if not len(frame) else working.iloc[positions].copy(deep=True)
    gathered['variant_rank'] = np.asarray(ranks, dtype='int64')
    gathered = gathered.drop('max_variant_score', axis=1)
    if list(gathered.columns) != [name for name, _ in expected] or [
            str(dtype) for dtype in gathered.dtypes] != [dtype for _, dtype in expected]:
        raise ValueError('pandas gather disagrees with native metadata')
    return gathered


def try_terminal_order(frame, binary):
    """Return the native R8 output, or None for an atomic Python fallback.

    The complete ranker supplies this explicit opt-in binary path. This adapter
    does not select a production backend from the environment or recompute any
    preceding control, compatibility, gene or expression weight.
    """
    try:
        executable = os.path.abspath(os.fspath(binary))
        with tempfile.TemporaryDirectory(prefix='priva-r8-') as directory:
            root = Path(directory)
            names, dtypes = _input(frame, root / 'input.csv')
            command = [executable, 'terminal-order', '--input', str(root / 'input.csv'),
                       '--output', str(root / 'output.csv'), '--metadata-output', str(root / 'metadata.csv'),
                       '--axis-metadata', _AXIS_TOKEN]
            for position, dtype in enumerate(dtypes):
                command.extend(('--column-dtype', f'{position}={dtype}'))
            env = os.environ.copy()
            env.setdefault('POLARS_MAX_THREADS', '1')
            result = subprocess.run(command, capture_output=True, env=env, check=False,
                                    timeout=_TIMEOUT_SECONDS)
            if result.returncode != 0:
                return None
            return _gather(frame, root / 'output.csv', root / 'metadata.csv', names, dtypes)
    except (MemoryError, RecursionError):
        raise
    except Exception:
        # Native transport/unsupported-input errors are not new scientific
        # exceptions: the original Python source handles this same current df.
        return None
