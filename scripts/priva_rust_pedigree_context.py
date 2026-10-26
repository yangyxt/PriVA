"""Optional R2 adapter for the qualified native pedigree-context procedure.

Native code selects the first affected pedigree row, sex flags and the whole
sample-string count('1') mask. Python retains original scalar and axis objects
and constructs the emitted Boolean data/null buffers without recomputing that
rule. Rejected requests leave the original R2 statement available on the same
inputs. Warnings belong to the caller, after an entire result is accepted.
"""
import copy
import csv
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import pandas as pd

from priva_rust_terminal_order import _cell, _plain_metadata


_DTYPES = frozenset(('object', 'string', 'float64', 'Float64', 'int64',
                     'Int64', 'uint64', 'UInt64', 'bool', 'boolean'))
_AXIS_TOKEN = 'caller-owned'
_TIMEOUT_SECONDS = 300
_META_NAMES = ('proband_present', 'proband', 'selected_pedigree_position',
               'is_male', 'is_female', 'mask_dtype', 'inherit_frame_metadata',
               'source_warning', 'axis')


def _plain_frame(frame):
    if type(frame) is not pd.DataFrame or type(frame.columns) is not pd.Index:
        raise ValueError('frame/column-axis subclass requires Python')
    if any(type(name) is not str for name in frame.columns):
        raise ValueError('nonstring column labels require Python')
    if type(frame.index) not in (pd.Index, pd.RangeIndex, pd.MultiIndex) or not (
            _plain_metadata(frame.attrs) and _plain_metadata(tuple(frame.index.names)) and
            _plain_metadata(tuple(frame.columns.names)) and
            all(_plain_metadata(label) for label in frame.index)):
        raise ValueError('custom metadata/axis labels require Python')


def _write_frame(frame, path):
    # Use private pandas buffers while preparing the request. No caller column
    # accessor or arbitrary Python-object callback is needed by this transport.
    working = frame.copy(deep=True)
    names, dtypes, columns = list(working.columns), [], []
    for position in range(len(names)):
        series = working.iloc[:, position]
        dtype = str(series.dtype)
        if dtype not in _DTYPES:
            raise ValueError('column dtype requires Python')
        if isinstance(series.dtype, pd.StringDtype) and series.dtype.storage != 'python':
            raise ValueError('Arrow string storage requires Python')
        dtypes.append(dtype)
        columns.append([_cell(value) for value in series.array])
    with path.open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value'] +
                        [part for name in names for part in (name + '_kind', name)])
        for row in range(len(working)):
            writer.writerow(['opaque', str(row)] + [part for column in columns for part in column[row]])
    return dtypes


def _input(frame, pedigree, family, root):
    if (pd.__version__, np.__version__) != ('2.2.3', '1.26.4'):
        raise ValueError('unqualified pandas/NumPy context')
    _plain_frame(frame)
    if isinstance(family, np.generic):
        raise ValueError('NumPy family scalar identity requires Python')
    active = pedigree is not None and family is not None
    if active:
        _plain_frame(pedigree)
        names = list(pedigree.columns)
        if any(names.count(name) != 1 for name in ('#FamilyID', 'Phenotype', 'IndividualID')) or names.count('Sex') > 1:
            raise ValueError('original source must handle required/ambiguous PED labels')
        families = pedigree.iloc[:, names.index('#FamilyID')]
        if str(families.dtype) not in ('object', 'string') or any(
                isinstance(value, np.generic) for value in families.array):
            raise ValueError('family comparison scalar identity requires Python')
    dtypes = _write_frame(frame, root / 'input.csv')
    pedigree_dtypes = _write_frame(pedigree, root / 'pedigree.csv') if active else None
    with (root / 'family.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['kind', 'value'])
        writer.writerow(_cell(family) if active else ('null', ''))
    return dtypes, pedigree_dtypes


def _gather(frame, pedigree, family, output, metadata):
    with output.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['source_position', 'index_kind', 'index_value', 'value', 'is_null']:
            raise ValueError('invalid native mask header')
        rows = list(reader)
    if len(rows) != len(frame) or any(
            row != dict(source_position=str(position), index_kind='opaque', index_value=str(position),
                        value=row['value'], is_null=row['is_null'])
            or row['value'] not in ('0', '1') or row['is_null'] not in ('0', '1')
            for position, row in enumerate(rows)):
        raise ValueError('native mask positions/identity/bits changed')
    with metadata.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ['name', 'value']:
            raise ValueError('invalid native context header')
        records = list(reader)
    if [row['name'] for row in records] != list(_META_NAMES) or any(set(row) != {'name', 'value'} for row in records):
        raise ValueError('native context fields changed')
    meta = {row['name']: row['value'] for row in records}
    if meta['axis'] != _AXIS_TOKEN or any(meta[name] not in ('0', '1') for name in (
            'proband_present', 'is_male', 'is_female', 'inherit_frame_metadata')):
        raise ValueError('native axis/context flags changed')
    active = pedigree is not None and family is not None
    if active:
        position = int(meta['selected_pedigree_position'])
        if not 0 <= position < len(pedigree) or meta['selected_pedigree_position'] != str(position) or (
                meta['proband_present'] != '1' or meta['inherit_frame_metadata'] != '1'):
            raise ValueError('native selected PED position changed')
        # The original .loc[mask] first makes a DataFrame and then .iloc[0]
        # constructs its Series. Preserve that take/row scalar construction.
        selected = pedigree.iloc[[position]].iloc[0]
        proband = selected['IndividualID']
        if type(proband) is not str or meta['proband'] != proband:
            raise ValueError('native selected proband/scalar identity changed')
        if meta['is_male'] == meta['is_female'] == '1' or meta['source_warning'] != (
                'none' if '1' in (meta['is_male'], meta['is_female']) else 'unknown_sex'):
            raise ValueError('native sex/warning flags changed')
        names = list(frame.columns)
        if names.count(proband) != 1:
            raise ValueError('native sample label changed')
        sample_dtype = str(frame.iloc[:, names.index(proband)].dtype)
        if sample_dtype in ('Int64', 'UInt64', 'Float64', 'boolean') or meta['mask_dtype'] != (
                'boolean' if sample_dtype == 'string' else 'bool'):
            raise ValueError('native selected-sample mask dtype changed')
    else:
        if [meta[name] for name in _META_NAMES] != [
                '0', '', '', '0', '0', 'object' if not len(frame) else 'bool', '0', 'no_pedigree', _AXIS_TOKEN]:
            raise ValueError('native disabled-PED context changed')
        selected, proband = None, None
    values = np.asarray([row['value'] == '1' for row in rows], dtype=bool)
    nulls = np.asarray([row['is_null'] == '1' for row in rows], dtype=bool)
    if meta['mask_dtype'] == 'boolean':
        array = pd.arrays.BooleanArray(values, nulls)
    else:
        if nulls.any():
            raise ValueError('native nonnullable mask contains nulls')
        array = values.astype(object) if meta['mask_dtype'] == 'object' else values
    mask = pd.Series(array, index=frame.index, name=proband)
    if meta['inherit_frame_metadata'] == '1':
        mask.attrs = copy.deepcopy(frame.attrs)
        mask.flags.allows_duplicate_labels = frame.flags.allows_duplicate_labels
    return dict(proband=proband, proband_row=selected, proband_het=mask,
                proband_is_male=meta['is_male'] == '1', proband_is_female=meta['is_female'] == '1',
                source_warning=meta['source_warning'])


def try_pedigree_context(frame, pedigree, family, binary):
    """Return accepted R2 locals, or None before any caller warning/mutation."""
    try:
        executable = os.path.abspath(os.fspath(binary))
        with tempfile.TemporaryDirectory(prefix='priva-r2-') as directory:
            root = Path(directory)
            dtypes, pedigree_dtypes = _input(frame, pedigree, family, root)
            command = [executable, 'pedigree-context', '--input', str(root / 'input.csv'),
                       '--family-input', str(root / 'family.csv'), '--output', str(root / 'output.csv'),
                       '--metadata-output', str(root / 'metadata.csv'), '--axis-metadata', _AXIS_TOKEN]
            for position, dtype in enumerate(dtypes):
                command.extend(('--column-dtype', f'{position}={dtype}'))
            if pedigree_dtypes is not None:
                command.extend(('--pedigree-input', str(root / 'pedigree.csv')))
                for position, dtype in enumerate(pedigree_dtypes):
                    command.extend(('--pedigree-column-dtype', f'{position}={dtype}'))
            env = os.environ.copy()
            env.setdefault('POLARS_MAX_THREADS', '1')
            result = subprocess.run(command, capture_output=True, env=env, check=False,
                                    timeout=_TIMEOUT_SECONDS)
            if result.returncode != 0:
                return None
            return _gather(frame, pedigree, family, root / 'output.csv', root / 'metadata.csv')
    except (MemoryError, RecursionError):
        raise
    except Exception:
        return None
