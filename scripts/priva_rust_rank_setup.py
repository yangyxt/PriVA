"""Optional R1 adapter: native copy/drop/score-copy/filter plan, original buffers.

The Rust procedure decides the column plan and protein-coding row selection.
Pandas retains caller-owned values and axes for the subsequent original ranking
stages. A rejected request leaves the input unchanged for the original Python
statements, including their own failures. No score arithmetic occurs here.
"""
import csv
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import pandas as pd

from priva_rust_terminal_order import _cell, _plain_metadata


_STALE = ('haplo_insuf_index', 'zygosity_inheritance_mechanism_index')
_DTYPES = frozenset(('object', 'string', 'float64', 'Float64', 'int64',
                     'Int64', 'uint64', 'UInt64', 'bool', 'boolean'))
_AXIS_TOKEN = 'caller-owned'
_TIMEOUT_SECONDS = 300


def _input(frame, path):
    if (pd.__version__, np.__version__) != ('2.2.3', '1.26.4'):
        raise ValueError('unqualified pandas/NumPy context')
    if type(frame) is not pd.DataFrame or type(frame.columns) is not pd.Index:
        raise ValueError('frame/column-axis subclass requires Python')
    names = list(frame.columns)
    if any(type(name) is not str for name in names):
        raise ValueError('nonstring column labels require Python')
    if any(names.count(name) != 1 for name in ('ACMG_quant_score', 'BIOTYPE')) or names.count('sort_index') > 1:
        raise ValueError('original source must handle required/assigned labels')
    if type(frame.index) not in (pd.Index, pd.RangeIndex, pd.MultiIndex) or not (
            _plain_metadata(frame.attrs) and _plain_metadata(tuple(frame.index.names)) and
            _plain_metadata(tuple(frame.columns.names)) and
            all(_plain_metadata(label) for label in frame.index)):
        raise ValueError('custom metadata/axis labels require Python')
    dtypes, columns = [], []
    for position, name in enumerate(names):
        series = frame.iloc[:, position]
        dtype = str(series.dtype)
        dtypes.append(dtype)
        if name in _STALE:
            # The original source drops these columns before examining cells.
            # Retain their dtype/position, without evaluating dropped objects.
            columns.append(None)
            continue
        if dtype not in _DTYPES:
            raise ValueError('column dtype requires Python')
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


def _column_plan(names, dtypes):
    score = names.index('ACMG_quant_score')
    columns = []
    for source, name in enumerate(names):
        if name in _STALE:
            continue
        copy_score = name == 'sort_index'
        origin = score if copy_score else source
        columns.append(dict(record='column', position=str(len(columns)), name=name,
            dtype=dtypes[origin], source_column=str(origin),
            operation='copy_score' if copy_score else 'retain'))
    if 'sort_index' not in names:
        columns.append(dict(record='column', position=str(len(columns)), name='sort_index',
            dtype=dtypes[score], source_column=str(score), operation='copy_score'))
    drops = [dict(record='drop', position=str(position), name=name,
                  dtype=dtypes[position], source_column=str(position), operation='remove_legacy')
             for stale in _STALE for position, name in enumerate(names) if name == stale]
    return columns, drops


def _gather(frame, output, metadata, names, dtypes):
    with output.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['source_position', 'index_kind', 'index_value']:
            raise ValueError('invalid native row header')
        rows = list(reader)
    positions = [int(row['source_position']) for row in rows]
    if any(position < 0 or position >= len(frame) for position in positions) or any(
            left >= right for left, right in zip(positions, positions[1:])):
        raise ValueError('native row positions are not an ordered subset')
    if any(row['index_kind'] != 'opaque' or row['index_value'] != str(position)
           for row, position in zip(rows, positions)):
        raise ValueError('native row identity changed')
    with metadata.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['record', 'position', 'name', 'dtype', 'source_column', 'operation']:
            raise ValueError('invalid native metadata header')
        actions = list(reader)
    columns, drops = _column_plan(names, dtypes)
    axis = dict(record='axis', position='', name=_AXIS_TOKEN, dtype='opaque',
                source_column='', operation='')
    if actions != columns + drops + [axis]:
        raise ValueError('native column/drop/axis plan changed')
    # 26 October 2026 CT: keep the source's mutation order. The native column
    # plan supplies the one score-copy assignment; the emitted row mask alone
    # selects the final rows. Do not recalculate BIOTYPE or score arithmetic.
    working = frame.copy(deep=True)
    for stale in _STALE:
        if stale in working.columns:
            working = working.drop(columns=[stale])
    score_copy = next(column for column in columns if column['operation'] == 'copy_score')
    working[score_copy['name']] = frame.iloc[:, int(score_copy['source_column'])]
    keep = np.zeros(len(frame), dtype=bool)
    keep[positions] = True
    gathered = working.loc[keep, :]
    if list(gathered.columns) != [column['name'] for column in columns] or [
            str(dtype) for dtype in gathered.dtypes] != [column['dtype'] for column in columns]:
        raise ValueError('pandas gather disagrees with native metadata')
    return gathered


def try_rank_setup(frame, binary):
    """Return native R1 output, or None before mutating the caller's frame."""
    try:
        executable = os.path.abspath(os.fspath(binary))
        with tempfile.TemporaryDirectory(prefix='priva-r1-') as directory:
            root = Path(directory)
            names, dtypes = _input(frame, root / 'input.csv')
            command = [executable, 'rank-setup', '--input', str(root / 'input.csv'),
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
        return None
