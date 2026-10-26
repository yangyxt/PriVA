"""Optional R3 product on the ranker's current post-control buffers.

Only the existing native float64 multiplication runs in Rust. The original
control calculations have already returned. Return a complete product Series,
or None without changing the frame, so the caller performs its one original
.loc assignment or its original multiplication/assignment fallback. Never
replay the control calculations, their warnings or their numerical libraries.
"""
import copy
import csv
import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd

from priva_rust_terminal_order import _cell, _plain_metadata


_OPERANDS = ('sort_index', 'control_common_index')
_AXIS_TOKEN = 'caller-owned'
_TIMEOUT_SECONDS = 300
_MAX_ROWS = 1_000_000  # Qualified standard NumPy route; larger inputs keep Python.


def _input(frame, path):
    if (pd.__version__, np.__version__) != ('2.2.3', '1.26.4'):
        raise ValueError('unqualified pandas/NumPy context')
    if type(frame) is not pd.DataFrame or type(frame.columns) is not pd.Index:
        raise ValueError('frame/column-axis subclass requires Python')
    if type(frame.index) not in (pd.Index, pd.RangeIndex) or len(frame) > _MAX_ROWS:
        raise ValueError('axis warning or numerical execution route requires Python')
    names = list(frame.columns)
    if any(type(name) is not str for name in names) or any(names.count(name) != 1 for name in _OPERANDS):
        raise ValueError('original source must handle required labels')
    if not (_plain_metadata(frame.attrs) and _plain_metadata(tuple(frame.index.names))
            and _plain_metadata(tuple(frame.columns.names))
            and all(_plain_metadata(label) for label in frame.index)):
        raise ValueError('custom metadata/axis labels require Python')
    positions = [names.index(name) for name in _OPERANDS]
    values = []
    for position in positions:
        series = frame.iloc[:, position]
        if type(series.dtype) is not type(np.dtype('float64')) or series.dtype != np.dtype('float64'):
            raise ValueError('native float64 operands required')
        values.append(series.to_numpy(copy=False))
    if np.any(np.isnan(values[0]) & np.isnan(values[1])):
        raise ValueError('two-NaN payload route requires Python')
    # Unconsumed columns retain their source positions, without inspecting
    # payload cells or invoking a custom payload dtype. Their placeholders do
    # not leave pandas, and only the native product is returned.
    dtypes = ['float64' if position in positions else 'object' for position in range(len(names))]
    consumed = {position: [_cell(value) for value in column] for position, column in zip(positions, values)}
    with path.open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value'] +
                        [part for name in names for part in (name + '_kind', name)])
        for row in range(len(frame)):
            writer.writerow(['opaque', str(row)] + [part for position in range(len(names))
                for part in (consumed[position][row] if position in consumed else ('null', ''))])
    return positions[0], dtypes


def _product(frame, output, metadata, target):
    with output.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['source_position', 'index_kind', 'index_value', 'float_bits']:
            raise ValueError('invalid native product header')
        rows = list(reader)
    if len(rows) != len(frame) or any(
            set(row) != {'source_position', 'index_kind', 'index_value', 'float_bits'}
            or row['source_position'] != str(position) or row['index_kind'] != 'opaque'
            or row['index_value'] != str(position)
            or re.fullmatch('[0-9a-f]{16}', row['float_bits']) is None
            for position, row in enumerate(rows)):
        raise ValueError('native product row identity or bits changed')
    with metadata.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['name', 'value']:
            raise ValueError('invalid native product metadata header')
        records = list(reader)
    if records != [dict(name='target_source_column_position', value=str(target)),
                   dict(name='dtype', value='float64'), dict(name='axis', value=_AXIS_TOKEN)]:
        raise ValueError('native target/dtype/axis changed')
    values = np.asarray([int(row['float_bits'], 16) for row in rows], dtype=np.uint64).view(np.float64)
    product = pd.Series(values, index=frame.index, name=None)
    product.attrs = copy.deepcopy(frame.attrs)
    product.flags.allows_duplicate_labels = frame.flags.allows_duplicate_labels
    return product


def try_control_multiply(frame, binary):
    """Return an atomic R3 product, or None for the current-buffer Python RHS."""
    try:
        executable = os.path.abspath(os.fspath(binary))
        with tempfile.TemporaryDirectory(prefix='priva-r3-') as directory:
            root = Path(directory)
            target, dtypes = _input(frame, root / 'input.csv')
            command = [executable, 'control-multiply', '--input', str(root / 'input.csv'),
                       '--output', str(root / 'output.csv'), '--metadata-output', str(root / 'metadata.csv'),
                       '--row-index-kind', 'flat', '--axis-metadata', _AXIS_TOKEN]
            for position, dtype in enumerate(dtypes):
                command.extend(('--column-dtype', f'{position}={dtype}'))
            env = os.environ.copy()
            env.setdefault('POLARS_MAX_THREADS', '1')
            result = subprocess.run(command, capture_output=True, env=env, check=False,
                                    timeout=_TIMEOUT_SECONDS)
            if result.returncode != 0:
                return None
            return _product(frame, root / 'output.csv', root / 'metadata.csv', target)
    except (MemoryError, RecursionError):
        raise
    except Exception:
        return None
