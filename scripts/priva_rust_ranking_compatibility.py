"""Optional atomic R4 plan on current post-control/post-R3 buffers.

The existing native procedure returns the complete mechanism/compatibility
masks. Python restores their original axes, metadata and aliases before the
caller emits its unchanged INFO call and two compatibility assignments. A
rejected request leaves the original R4 block on the same current frame; prior
pedigree warnings, control calculations and R3 are never replayed here.
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


_MECHANISM = (
    'modern_profile', 'is_exact_gof', 'is_exact_dn', 'is_exact_nonlof',
    'is_predicted_lof', 'is_uncertain', 'has_recessive_compatible',
    'has_dominant_compatible', 'has_x_linked_recessive_compatible',
    'has_x_linked_dominant_compatible', 'has_x_linked_unspecified_compatible',
    'has_rec_lof_history', 'has_rec_gof_history', 'has_rec_dn_history',
    'has_rec_unresolved_history', 'has_dom_lof_history', 'has_dom_gof_history',
    'has_dom_dn_history', 'has_dom_unresolved_history', 'has_x_recessive_lof_history',
    'has_x_recessive_gof_history', 'has_x_recessive_dn_history',
    'has_x_recessive_unresolved_history', 'has_x_dominant_lof_history',
    'has_x_dominant_gof_history', 'has_x_dominant_dn_history',
    'has_x_dominant_unresolved_history',
)
_ROLES = ('chrom', 'x_locus', 'has_single_alt_lof_compatible', 'exact_gof',
          'asserted_lof_effect', 'dn_without_hi_or_ambiguous',
          'heterozygous_lof_incompatible', 'high_acmg_lof',
          'compatibility_penalty_rows', 'compatibility_penalty_rows_null',
          'pext_modulation_rows')
_CONSUMED = ('var_plausible_patho_mechs', 'variant_effect', 'variant_lof_score',
             'variant_gof_score', 'variant_dn_score', 'ACMG_quant_score', 'chrom')
_REQUIRED = ('var_plausible_patho_mechs', 'variant_effect', 'ACMG_quant_score')
_DTYPES = frozenset(('object', 'string', 'float64', 'int64', 'uint64', 'bool'))
_TARGET = 'zygosity_inheritance_mechanism_compatibility'
_AXIS_TOKEN = 'caller-owned'
_TIMEOUT_SECONDS = 300
_MAX_ROWS = 1_000_000


def _input(frame, male, female, heterozygosity, root):
    if (pd.__version__, np.__version__) != ('2.2.3', '1.26.4'):
        raise ValueError('unqualified pandas/NumPy context')
    if type(frame) is not pd.DataFrame or type(frame.columns) is not pd.Index:
        raise ValueError('frame/column-axis subclass requires Python')
    if type(frame.index) not in (pd.Index, pd.RangeIndex, pd.MultiIndex) or len(frame) > _MAX_ROWS:
        raise ValueError('axis/size requires Python')
    names = list(frame.columns)
    if any(type(name) is not str for name in names) or any(names.count(name) != 1 for name in _REQUIRED):
        raise ValueError('original source must handle required labels')
    if any(names.count(name) > 1 for name in (*_CONSUMED, _TARGET)):
        raise ValueError('ambiguous consumed label requires Python')
    if not (_plain_metadata(frame.attrs) and _plain_metadata(tuple(frame.index.names))
            and _plain_metadata(tuple(frame.columns.names))
            and all(_plain_metadata(label) for label in frame.index)):
        raise ValueError('custom metadata/axis labels require Python')
    if type(male) is not bool or type(female) is not bool or type(heterozygosity) is not pd.Series:
        raise ValueError('R2 scalar/mask subclass requires Python')
    if heterozygosity.index is not frame.index or not (
            _plain_metadata(heterozygosity.name) and _plain_metadata(heterozygosity.attrs)):
        raise ValueError('unaligned/custom R2 mask requires Python')
    mask_dtype = str(heterozygosity.dtype)
    if mask_dtype == 'boolean' and type(heterozygosity.array) is pd.arrays.BooleanArray:
        # The pinned pandas representation retains hidden data beneath nulls.
        # Both buffers are transported; do not reconstruct data from fillna.
        data, nulls = heterozygosity.array._data, heterozygosity.array._mask
        if type(data) is not np.ndarray or type(nulls) is not np.ndarray or (
                data.dtype != np.dtype('bool') or nulls.dtype != np.dtype('bool')):
            raise ValueError('R2 nullable storage requires Python')
    elif mask_dtype == 'bool' and heterozygosity.dtype == np.dtype('bool'):
        data = heterozygosity.to_numpy(copy=False)
        nulls = np.zeros(len(heterozygosity), dtype=bool)
    elif mask_dtype == 'object' and not len(heterozygosity):
        data = nulls = np.empty(0, dtype=bool)
    else:
        raise ValueError('R2 mask dtype requires Python')
    sources, columns, dtypes = {}, {}, ['object'] * len(names)
    for name in _CONSUMED:
        if name not in names:
            continue
        position = names.index(name)
        series = frame[name]
        dtype = str(series.dtype)
        if dtype not in _DTYPES or (isinstance(series.dtype, pd.StringDtype) and series.dtype.storage != 'python'):
            raise ValueError('consumed dtype/storage requires Python')
        if name == 'ACMG_quant_score' and dtype not in ('float64', 'int64', 'uint64', 'bool'):
            raise ValueError('raw score dtype requires Python')
        if name == 'chrom' and dtype not in ('object', 'string', 'int64', 'uint64', 'bool'):
            raise ValueError('chromosome dtype requires Python')
        if not _plain_metadata(series.attrs):
            raise ValueError('custom source-Series metadata requires Python')
        sources[position] = series
        columns[position] = [_cell(value) for value in series.array]
        dtypes[position] = dtype
    with (root / 'input.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value'] +
                        [part for name in names for part in (name + '_kind', name)])
        for row in range(len(frame)):
            writer.writerow(['opaque', str(row)] + [part for position in range(len(names))
                for part in (columns[position][row] if position in columns else ('null', ''))])
    with (root / 'heterozygosity.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['index_kind', 'index_value', 'value_kind', 'value', 'raw_data', 'is_null'])
        for row in range(len(frame)):
            writer.writerow(['opaque', str(row), 'pd_na' if nulls[row] else 'bool',
                '' if nulls[row] else str(int(data[row])), str(int(data[row])) if mask_dtype == 'boolean' else '',
                str(int(nulls[row]))])
    return dtypes, mask_dtype, sources


def _series_contract(name, names):
    if name in _MECHANISM:
        position = _MECHANISM.index(name)
        label = None if 1 <= position <= 4 else ('variant_effect' if position == 5 else 'var_plausible_patho_mechs')
    else:
        label = ('chrom' if name in ('chrom', 'x_locus') and 'chrom' in names else
                 'var_plausible_patho_mechs' if name == 'dn_without_hi_or_ambiguous' else
                 'ACMG_quant_score' if name in ('high_acmg_lof', 'pext_modulation_rows') else None)
    if name in ('chrom', 'x_locus'):
        origin = names.index('chrom') if 'chrom' in names else None
    elif name in ('is_exact_gof', 'is_exact_dn', 'is_exact_nonlof', 'exact_gof'):
        origin = None
    elif name in ('is_predicted_lof', 'is_uncertain', 'asserted_lof_effect', 'heterozygous_lof_incompatible'):
        origin = names.index('variant_effect')
    elif name in ('high_acmg_lof', 'pext_modulation_rows', 'compatibility_penalty_rows', 'compatibility_penalty_rows_null'):
        origin = names.index('ACMG_quant_score')
    else:
        origin = names.index('var_plausible_patho_mechs')
    return label, origin


def _plan(frame, heterozygosity, sources, root):
    fields = ['source_position', 'index_kind', 'index_value', 'chrom_hex', 'compatibility_bits',
              *_MECHANISM, *[name for name in _ROLES if name != 'chrom']]
    with (root / 'output.csv').open(newline='') as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != fields:
            raise ValueError('invalid native compatibility header')
        rows = list(reader)
    if len(rows) != len(frame) or any(set(row) != set(fields) or row['source_position'] != str(position)
            or row['index_kind'] != 'opaque' or row['index_value'] != str(position)
            or re.fullmatch('(?:[0-9a-f]{2})*', row['chrom_hex']) is None
            or any(row[name] not in ('0', '1') for name in fields[5:])
            or row['compatibility_bits'] not in ('3ff0000000000000', '3feccccccccccccd')
            for position, row in enumerate(rows)):
        raise ValueError('native row identity, mask or compatibility bits changed')
    names = list(frame.columns)
    penalty_dtype = 'boolean' if str(heterozygosity.dtype) == 'boolean' else 'bool'
    count = lambda name: sum(row[name] == '1' for row in rows)
    penalty_count = sum(row['compatibility_penalty_rows'] == '1' and row['compatibility_penalty_rows_null'] == '0' for row in rows)
    log_message = ('Zygosity/inheritance/mechanism compatibility signals for ranking: '
        f"hub_profile_rows={count('modern_profile')}, asserted_lof_effect={count('is_predicted_lof')}, "
        f"dominant_DN_without_HI_or_ambiguous={count('dn_without_hi_or_ambiguous')}, "
        f"heterozygous_lof_incompatible={count('heterozygous_lof_incompatible')}, "
        f'compatibility_penalty_rows={penalty_count}')
    expected = [('target_source_column_position', str(names.index(_TARGET) if _TARGET in names else len(names))),
                ('penalty_dtype', penalty_dtype), ('axis', _AXIS_TOKEN), ('log_message', log_message)]
    for name in (*_MECHANISM, *_ROLES):
        label, origin = _series_contract(name, names)
        expected.extend(((name + '.series_name_kind', 'none' if label is None else 'str'),
                         (name + '.series_name', '' if label is None else label),
                         (name + '.metadata_origin', 'fresh' if origin is None else 'column:' + str(origin))))
    with (root / 'metadata.csv').open(newline='') as stream:
        reader = csv.DictReader(stream)
        records = list(reader)
        if reader.fieldnames != ['name', 'value'] or records != [dict(name=name, value=value) for name, value in expected]:
            raise ValueError('native complete metadata/INFO contract changed')
    if any((penalty_dtype != 'boolean' and row['compatibility_penalty_rows_null'] != '0')
            or row['compatibility_bits'] != ('3feccccccccccccd' if row['compatibility_penalty_rows'] == '1'
                and row['compatibility_penalty_rows_null'] == '0' else '3ff0000000000000')
            or row['exact_gof'] != row['is_exact_gof'] or row['asserted_lof_effect'] != row['is_predicted_lof']
            or row['pext_modulation_rows'] != row['high_acmg_lof'] for row in rows):
        raise ValueError('native alias/null/publication contract changed')
    def series(name, values, dtype):
        label, origin = _series_contract(name, names)
        result = pd.Series(values, index=frame.index, name=label, dtype=dtype)
        if origin is not None:
            result.attrs = copy.deepcopy(sources[origin].attrs)
            result.flags.allows_duplicate_labels = sources[origin].flags.allows_duplicate_labels
        return result
    masks = {name: series(name, [row[name] == '1' for row in rows], 'bool') for name in _MECHANISM}
    plan = dict(mechanism_masks=masks,
                chrom=series('chrom', [bytes.fromhex(row['chrom_hex']).decode('utf-8') for row in rows], 'object'))
    for name in ('x_locus', 'has_single_alt_lof_compatible', 'dn_without_hi_or_ambiguous',
                 'heterozygous_lof_incompatible', 'high_acmg_lof'):
        plan[name] = series(name, [row[name] == '1' for row in rows], 'bool')
    plan['exact_gof'] = masks['is_exact_gof']
    plan['asserted_lof_effect'] = masks['is_predicted_lof']
    data = np.asarray([row['compatibility_penalty_rows'] == '1' for row in rows], dtype=bool)
    nulls = np.asarray([row['compatibility_penalty_rows_null'] == '1' for row in rows], dtype=bool)
    values = pd.arrays.BooleanArray(data, nulls) if penalty_dtype == 'boolean' else data
    plan['compatibility_penalty_rows'] = series('compatibility_penalty_rows', values, penalty_dtype)
    plan['pext_modulation_rows'] = plan['high_acmg_lof']
    return {'mechanism_masks': masks, **{name: plan[name] for name in _ROLES
                                       if name != 'compatibility_penalty_rows_null'}}


def try_ranking_compatibility(frame, male, female, heterozygosity, binary):
    """Return all R4 locals atomically, or None before caller publication."""
    try:
        executable = os.path.abspath(os.fspath(binary))
        with tempfile.TemporaryDirectory(prefix='priva-r4-') as directory:
            root = Path(directory)
            dtypes, mask_dtype, sources = _input(frame, male, female, heterozygosity, root)
            command = [executable, 'ranking-compatibility', '--input', str(root / 'input.csv'),
                '--heterozygosity-input', str(root / 'heterozygosity.csv'), '--heterozygosity-dtype', mask_dtype,
                '--proband-is-male', str(male).lower(), '--proband-is-female', str(female).lower(),
                '--output', str(root / 'output.csv'), '--metadata-output', str(root / 'metadata.csv'),
                '--axis-metadata', _AXIS_TOKEN]
            for position, dtype in enumerate(dtypes):
                command.extend(('--column-dtype', f'{position}={dtype}'))
            env = os.environ.copy()
            env.setdefault('POLARS_MAX_THREADS', '1')
            result = subprocess.run(command, capture_output=True, env=env, check=False, timeout=_TIMEOUT_SECONDS)
            if result.returncode != 0:
                return None
            return _plan(frame, heterozygosity, sources, root)
    except (MemoryError, RecursionError):
        raise
    except Exception:
        return None
