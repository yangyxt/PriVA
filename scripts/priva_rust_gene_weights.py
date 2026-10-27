"""Optional complete R5/R6 plan on the ranker's actual post-R4 buffers.

Native parsing/membership and the three ordered products are validated before
the caller changes its frame. Rejection leaves the original two gene branches
and compatibility multiplication on that same frame, without earlier replay.
"""
import copy
import csv
import os
from pathlib import Path
import re
import stat
import subprocess
import tempfile

import numpy as np
import pandas as pd

from priva_rust_terminal_order import _cell, _plain_metadata


_SCORE = 'sort_index'
_COMPATIBILITY = 'zygosity_inheritance_mechanism_compatibility'
_BRANCHES = (('dispensable', 'dispensable_gene_index', 0.9, 604),
             ('relevant', 'relevant_gene_index', 1.2, 611))
_AXIS_TOKEN = 'caller-owned'
_TIMEOUT_SECONDS = 300
_MAX_ROWS = 1_000_000
_INFERENCE_TOKENS = frozenset(('inf', 'infinity', 'nan', 'true', 'false', 'null', 'none'))


def _identity(path):
    value = path.lstat()
    if not stat.S_ISREG(value.st_mode):
        raise ValueError('literal regular gene-list route required')
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns,
            value.st_ctime_ns, value.st_mode, value.st_nlink)


def _files(arguments):
    result = []
    for value in arguments:
        if value is None:
            result.append(None)
            continue
        if type(value) not in (str, type(Path())):
            raise ValueError('custom gene-file argument requires Python')
        path = Path(os.path.abspath(os.fspath(value)))
        # pandas infers compression from selected names even if their bytes
        # happen to be plain ASCII. Those parser routes remain original.
        if path.name.lower().endswith(('.gz', '.bz2', '.zip', '.xz', '.zst', '.tar', '.tgz')):
            raise ValueError('compressed/parser route requires Python')
        result.append((path, _identity(path)))
    return result


def _input(frame, mask, enabled, root):
    if (pd.__version__, np.__version__) != ('2.2.3', '1.26.4'):
        raise ValueError('unqualified pandas/NumPy context')
    if type(frame) is not pd.DataFrame or type(frame.columns) is not pd.Index:
        raise ValueError('frame/column-axis subclass requires Python')
    if type(frame.index) not in (pd.Index, pd.RangeIndex) or len(frame) > _MAX_ROWS:
        raise ValueError('axis warning or expression-size route requires Python')
    names = list(frame.columns)
    consumed = (_SCORE, _COMPATIBILITY, 'SYMBOL') if any(enabled) else (_SCORE, _COMPATIBILITY)
    targets = [target for active, (_, target, _, _) in zip(enabled, _BRANCHES) if active]
    if any(type(name) is not str for name in names) or any(names.count(name) != 1 for name in consumed):
        raise ValueError('original source must handle required labels')
    if any(names.count(name) > 1 for name in targets):
        raise ValueError('ambiguous assigned target requires Python')
    if not (_plain_metadata(frame.attrs) and _plain_metadata(tuple(frame.index.names))
            and _plain_metadata(tuple(frame.columns.names))
            and all(_plain_metadata(label) for label in frame.index)):
        raise ValueError('custom frame metadata requires Python')
    if type(mask) is not pd.Series or mask.index is not frame.index or not (
            _plain_metadata(mask.name) and _plain_metadata(mask.attrs)):
        raise ValueError('unaligned/custom R4 mask requires Python')
    dtype = str(mask.dtype)
    if dtype == 'boolean' and type(mask.array) is pd.arrays.BooleanArray:
        data, nulls = mask.array._data, mask.array._mask
        if type(data) is not np.ndarray or type(nulls) is not np.ndarray or (
                data.dtype != np.dtype('bool') or nulls.dtype != np.dtype('bool')):
            raise ValueError('raw R4 Boolean storage requires Python')
    elif dtype == 'bool' and mask.dtype == np.dtype('bool'):
        data, nulls = mask.to_numpy(copy=False), np.zeros(len(mask), dtype=bool)
    elif dtype == 'object' and not len(mask):
        data = nulls = np.empty(0, dtype=bool)
    else:
        raise ValueError('R4 mask dtype requires Python')
    dtypes, columns, source_values = ['object'] * len(names), {}, {}
    for name in consumed:
        position = names.index(name)
        source = frame.iloc[:, position]
        source_dtype = str(source.dtype)
        if name in (_SCORE, _COMPATIBILITY):
            if source.dtype != np.dtype('float64') or type(source.dtype) is not type(np.dtype('float64')):
                raise ValueError('native float64 operands required')
            source_values[name] = source.to_numpy(copy=False)
        elif source_dtype not in ('object', 'string', 'int64', 'uint64', 'float64', 'bool') or (
                isinstance(source.dtype, pd.StringDtype) and source.dtype.storage != 'python'):
            raise ValueError('SYMBOL representation requires Python')
        if name == 'SYMBOL' and not _plain_metadata(frame[name].attrs):
            raise ValueError('custom source-Series metadata requires Python')
        dtypes[position] = source_dtype
        columns[position] = [_cell(value) for value in source.array]
    with (root / 'input.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value'] +
                        [part for name in names for part in (name + '_kind', name)])
        for row in range(len(frame)):
            writer.writerow(['opaque', str(row)] + [part for position in range(len(names))
                for part in (columns[position][row] if position in columns else ('null', ''))])
    with (root / 'penalty.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value', 'value_kind', 'value', 'raw_data', 'is_null'])
        for row in range(len(frame)):
            writer.writerow(['opaque', str(row), 'pd_na' if nulls[row] else 'bool',
                '' if nulls[row] else str(int(data[row])), str(int(data[row])) if dtype == 'boolean' else '',
                str(int(nulls[row]))])
    return names, dtypes, dtype, data & ~nulls, source_values[_COMPATIBILITY]


def _series(frame, values, index=None):
    result = pd.Series(values, index=frame.index if index is None else index, name=None, dtype='float64')
    result.attrs = copy.deepcopy(frame.attrs)
    result.flags.allows_duplicate_labels = frame.flags.allows_duplicate_labels
    return result


def gene_membership(frame, values):
    """Restore membership metadata at the original post-column-creation point."""
    source = frame['SYMBOL']
    result = pd.Series(values, index=frame.index, name=source.name, dtype='bool')
    result.attrs = copy.deepcopy(source.attrs)
    result.flags.allows_duplicate_labels = source.flags.allows_duplicate_labels
    return result


def _plan(frame, names, enabled, selected, compatibility, root):
    steps, planned_names = [], list(names)
    for active, (role, target, factor, line) in zip(enabled, _BRANCHES):
        if active:
            position = planned_names.index(target) if target in planned_names else len(planned_names)
            if position == len(planned_names):
                planned_names.append(target)
            steps.append((role, target, factor, line, position))
    steps.append(('compatibility', _COMPATIBILITY, None, 613, names.index(_COMPATIBILITY)))
    with (root / 'output.csv').open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['step', 'source_position', 'index_kind', 'index_value', 'selected',
                                 'multiplier_bits', 'score_bits']:
            raise ValueError('invalid native weight header')
        rows = list(reader)
    fields = {'step', 'source_position', 'index_kind', 'index_value', 'selected', 'multiplier_bits', 'score_bits'}
    if len(rows) != len(frame) * len(steps) or any(set(row) != fields or row['selected'] not in ('0', '1')
            or re.fullmatch('[0-9a-f]{16}', row['multiplier_bits']) is None
            or re.fullmatch('[0-9a-f]{16}', row['score_bits']) is None for row in rows):
        raise ValueError('incomplete or malformed native weight rows')
    plan = {}
    for step, (role, target, factor, line, position) in enumerate(steps):
        items = rows[step * len(frame):(step + 1) * len(frame)]
        if any(row['step'] != str(step) or row['source_position'] != str(i)
                or row['index_kind'] != 'opaque' or row['index_value'] != str(i)
                for i, row in enumerate(items)):
            raise ValueError('native step/row identity changed')
        membership = np.asarray([row['selected'] == '1' for row in items], dtype=bool)
        multiplier_bits = np.asarray([int(row['multiplier_bits'], 16) for row in items], dtype=np.uint64)
        if factor is None:
            if not np.array_equal(membership, selected) or not np.array_equal(
                    multiplier_bits, compatibility.view(np.uint64)):
                raise ValueError('native current R4 mask/multiplier changed')
        else:
            expected = np.asarray([factor if value else 1.0 for value in membership], dtype='float64').view(np.uint64)
            if not np.array_equal(multiplier_bits, expected):
                raise ValueError('native gene multiplier changed')
        values = np.asarray([int(row['score_bits'], 16) for row in items], dtype=np.uint64).view(np.float64)
        plan[role + '_selected'] = membership
        plan[role + '_product'] = _series(frame, values if factor is not None else values[selected],
            None if factor is not None else frame.index[selected])
    expected = [('axis', _AXIS_TOKEN), ('score_column_position', str(names.index(_SCORE))), ('step_count', str(len(steps)))]
    for i, (_, target, _, line, position) in enumerate(steps):
        expected += [(f'step.{i}.name', target), (f'step.{i}.line', str(line)),
                     (f'step.{i}.multiplier_column_position', str(position))]
    with (root / 'metadata.csv').open(newline='') as stream:
        reader = csv.DictReader(stream)
        records = list(reader)
        if reader.fieldnames != ['name', 'value']:
            raise ValueError('invalid native weight metadata header')
    position = len(expected)
    for active, (role, _, _, _) in zip(enabled, _BRANCHES):
        expected.append((role + '.enabled', str(active).lower())); position += 1
        if active:
            if position >= len(records) or records[position]['name'] != role + '.count' or (
                    re.fullmatch('0|[1-9][0-9]*', records[position]['value']) is None):
                raise ValueError('invalid native gene-list count')
            count = int(records[position]['value'])
            if count > len(records):
                raise ValueError('incomplete native gene-list metadata')
            expected.append((role + '.count', str(count))); position += 1
            genes = []
            for i in range(count):
                if position >= len(records):
                    raise ValueError('incomplete native gene list')
                record = records[position]
                if record['name'] != f'{role}.gene.{i}' or re.fullmatch('[A-Za-z_][A-Za-z0-9_.-]*', record['value']) is None or (
                        record['value'] == 'NA' or record['value'].lower() in _INFERENCE_TOKENS):
                    raise ValueError('invalid native gene-list record')
                expected.append((record['name'], record['value'])); genes.append(record['value']); position += 1
            plan[role + '_genes'] = set(genes)
    if records != [dict(name=name, value=value) for name, value in expected]:
        raise ValueError('native complete ordered metadata changed')
    return plan


def try_gene_weights(frame, dispensable, relevant, mask, binary):
    """Return the complete ordered R5/R6 plan, or None before publication."""
    try:
        enabled = (dispensable is not None, relevant is not None)
        executable = os.path.abspath(os.fspath(binary))
        with tempfile.TemporaryDirectory(prefix='priva-r5-r6-') as directory:
            root = Path(directory)
            names, dtypes, dtype, selected, compatibility = _input(frame, mask, enabled, root)
            files = _files((dispensable, relevant))
            command = [executable, 'ranking-gene-weights', '--input', str(root / 'input.csv'),
                '--penalty-input', str(root / 'penalty.csv'), '--penalty-dtype', dtype,
                '--row-index-kind', 'flat', '--axis-metadata', _AXIS_TOKEN,
                '--output', str(root / 'output.csv'), '--metadata-output', str(root / 'metadata.csv')]
            for (role, _, _, _), record in zip(_BRANCHES, files):
                command += ['--' + role + '-kind', 'none' if record is None else 'file']
                if record is not None:
                    command += ['--' + role + '-file', str(record[0])]
            for position, dtype in enumerate(dtypes):
                command.extend(('--column-dtype', f'{position}={dtype}'))
            env = os.environ.copy(); env.setdefault('POLARS_MAX_THREADS', '1')
            result = subprocess.run(command, capture_output=True, env=env, check=False, timeout=_TIMEOUT_SECONDS)
            if result.returncode != 0:
                return None
            plan = _plan(frame, names, enabled, selected, compatibility, root)
            if any(_identity(path) != before for record in files if record is not None for path, before in [record]):
                return None
            return plan
    except (MemoryError, RecursionError):
        raise
    except Exception:
        return None
