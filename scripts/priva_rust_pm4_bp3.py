"""Optional PM4/BP3 procedures with unchanged Python defaults and fallback.

Native code owns full-REF BED coordinates, ID membership and criterion masks.
Mature pandas publication and pybedtools/bedtools intersection remain in place.
Source-order logging and one actual interval intersection are preserved.
"""
import ast
import copy
import csv
import hashlib
import io
import logging
import os
from pathlib import Path
import stat
import struct
import subprocess
import sys
from types import FunctionType
import uuid

import numpy as np
import pandas as pd

import priva_rust_pm4_bp3_source as guards
import acmg_pm4_bp3_protein_length as _reference
import shutil
import tempfile

FLAGS = ('splicing_len_changing', 'splicing_frameshift', '5UTR_len_changing',
    '5UTR_frameshift', 'vep_consq_len_changing', 'splicing_span_intol_domain',
    '5UTR_span_intol_domain', 'span_functional_domains')
MASKS = ('in_repeat', 'low_complexity_domain', 'in_repetitive_region', 'inframe_indels',
    'splicing_inframe', 'utr_inframe', 'vep_inframe', 'all_inframe', 'stop_loss_variants',
    'nmd_escaping_trunc', 'has_functional_domain', 'frameshift_with_pvs1', 'strong_nmd',
    'exclude_from_pm4', 'pm4_inframe_eligible', 'pm4_eligible', 'bp3_eligible')


def sha(path): return hashlib.sha256(guards.regular_bytes(path)).hexdigest()


def raw_columns(frame):
    """Read ordinary block storage without changing item or manager caches."""
    from pandas.core.internals.managers import BlockManager
    from pandas.core.internals.blocks import NumpyBlock, ExtensionBlock
    if type(frame._mgr) is not BlockManager: raise ValueError('custom manager')
    result = [None] * len(frame.columns)
    for block in frame._mgr.blocks:
        if type(block) not in (NumpyBlock, ExtensionBlock): raise ValueError('custom block')
        values = block.values
        if type(values) is np.ndarray:
            for position, column in enumerate(block._mgr_locs.as_array): result[int(column)] = values[position]
        elif type(values) in (pd.arrays.BooleanArray, pd.arrays.StringArray):
            if len(block._mgr_locs) != 1: raise ValueError('ambiguous extension placement')
            result[int(block._mgr_locs.as_array[0])] = values
        else: raise ValueError('custom raw array')
    return result


def strings(values, *, default=None):
    if type(values) is np.ndarray and values.dtype == np.dtype('O'):
        result = []
        for value in values:
            if type(value) is str: result.append(value)
            elif default is not None and (value is None or value is pd.NA
                    or type(value) in (float, np.float64) and (
                        int.from_bytes(struct.pack('>d', value), 'big') & 0x7fffffffffffffff) > 0x7ff0000000000000):
                result.append(default)
            else: raise ValueError('unsupported text scalar')
        return result
    if type(values) is pd.arrays.StringArray and type(values._ndarray) is np.ndarray:
        return strings(values._ndarray, default=default)
    raise ValueError('unsupported string storage')


def booleans(values):
    if type(values) is np.ndarray and values.dtype == np.dtype('bool'):
        return values.tolist()
    if (type(values) is pd.arrays.BooleanArray and type(values._data) is np.ndarray
            and type(values._mask) is np.ndarray and values._data.dtype == np.dtype('bool')
            and values._mask.dtype == np.dtype('bool')):
        return [False if missing else bool(value) for value, missing in zip(values._data, values._mask)]
    raise ValueError('unsupported boolean storage')


def integer(values, count):
    if type(values) is not np.ndarray or values.dtype != np.dtype('int64') or values.shape != (count,):
        raise ValueError('unsupported integer vector')
    return values.tolist()


def logger_guard(logger):
    # Arbitrary handlers/streams can mutate inputs or throw at different phases.
    # They remain on the untouched original path; no custom method is called here.
    guards.logger_instance_guard(logger)
    if type(logger) is not logging.Logger or logger.disabled or logger.filters: raise ValueError('custom logger')
    current = logger
    while current is not None:
        if current.filters: raise ValueError('custom logger filter')
        for handler in current.handlers:
            if type(handler) is logging.NullHandler: continue
            if type(handler) is not logging.StreamHandler or handler.filters: raise ValueError('custom logging handler')
            if type(handler.stream) not in (io.StringIO, io.TextIOWrapper): raise ValueError('custom logging stream')
            if handler.formatter is not None and type(handler.formatter) is not logging.Formatter:
                raise ValueError('custom formatter')
        if not current.propagate: break
        current = current.parent


class _Bridge:
    def __init__(self, reference, *, binary, binary_sha256, scratch, package_binding, package_binding_sha256):
        self.reference = reference
        self.binary = Path(binary) if binary is not None else None
        self.binary_sha256 = binary_sha256
        self.scratch = Path(scratch)
        self.package_binding = Path(package_binding)
        self.package_binding_sha256 = package_binding_sha256
        self.events = []
        self.last_masks = None
        self._tail = None
        self._emitters = None

    def _original_tail(self):
        raw = guards.regular_bytes(self.reference.__file__)
        if hashlib.sha256(raw).hexdigest() != guards.REFERENCE_SHA: raise ValueError('original source changed')
        tree = ast.parse(raw)
        matches = [node for node in tree.body if type(node) is ast.FunctionDef and node.name == 'PM4_BP3_criteria']
        if len(matches) != 1: raise ValueError('ambiguous original tail')
        node = copy.deepcopy(matches[0])
        first = node.body[1]
        if (type(first) is not ast.Assign or len(first.targets) != 1 or type(first.targets[0]) is not ast.Name
                or first.targets[0].id != 'in_repeat_vars' or type(first.value) is not ast.Call
                or type(first.value.func) is not ast.Name or first.value.func.id != 'find_overlaps_bedtools_efficient'):
            raise ValueError('original first assignment changed')
        original_body = copy.deepcopy(node.body[2:])
        # The real caller formatter includes funcName/lineno. Retain the actual
        # original function name and source line attributes in this continuation.
        node.name = 'PM4_BP3_criteria'
        node.args.args.append(ast.arg(arg='in_repeat_vars'))
        node.body = [node.body[0]] + original_body
        ast.fix_missing_locations(node)
        assert ast.dump(ast.Module(body=node.body[1:], type_ignores=[]), include_attributes=True) == ast.dump(
            ast.Module(body=matches[0].body[2:], type_ignores=[]), include_attributes=True)
        compiled = compile(ast.Module(body=[node], type_ignores=[]), str(self.reference.__file__), 'exec')
        code = next(value for value in compiled.co_consts if hasattr(value, 'co_name') and value.co_name == node.name)
        return FunctionType(code, self.reference.__dict__, node.name)

    def _original_log_sites(self):
        raw = guards.regular_bytes(self.reference.__file__)
        if hashlib.sha256(raw).hexdigest() != guards.REFERENCE_SHA: raise ValueError('original source changed')
        tree = ast.parse(raw)
        function = next(node for node in tree.body if type(node) is ast.FunctionDef and node.name == 'PM4_BP3_criteria')
        sites = [node for node in function.body if type(node) is ast.Expr and type(node.value) is ast.Call
                 and type(node.value.func) is ast.Attribute and type(node.value.func.value) is ast.Name
                 and node.value.func.value.id == 'logger']
        branch = [node for node in function.body if type(node) is ast.If][-1]
        sites += [branch.orelse[0], branch.body[0]]
        if len(sites) != 9: raise ValueError('original logging site count changed')
        names = ('in_repetitive_region', 'all_inframe', 'stop_loss_variants', 'nmd_escaping_trunc',
                 'has_functional_domain', 'pm4_array', 'bp3_array', 'both_assigned')
        emitters = []
        for site in sites:
            node = copy.deepcopy(function)
            node.args = ast.arguments(posonlyargs=[], args=[ast.arg(arg=name) for name in names],
                                      kwonlyargs=[], kw_defaults=[], defaults=[])
            node.body = [copy.deepcopy(site)]
            ast.fix_missing_locations(node)
            assert ast.dump(node.body[0], include_attributes=True) == ast.dump(site, include_attributes=True)
            compiled = compile(ast.Module(body=[node], type_ignores=[]), str(self.reference.__file__), 'exec')
            code = next(item for item in compiled.co_consts if hasattr(item, 'co_name') and item.co_name == function.name)
            emitters.append((FunctionType(code, self.reference.__dict__, function.name), site.lineno))
        return emitters

    def _admit(self, df, *, whole=False, pvs=None, loc=None):
        if self.binary is None: raise ValueError('Python default')
        guards.reference_guard(self.reference)
        guards.package_guard(self.package_binding, self.package_binding_sha256)
        if type(df) is not pd.DataFrame or len(df) == 0: raise ValueError('nonordinary or empty frame')
        guards.frame_instance_guard(df)
        if type(df._item_cache) is not dict or not guards.plain(df.attrs): raise ValueError('custom frame metadata')
        if type(df.index) not in (pd.Index, pd.RangeIndex, pd.MultiIndex): raise ValueError('custom index')
        if type(df.columns) is not pd.Index or type(df.columns._values) is not np.ndarray: raise ValueError('custom columns')
        names = list(df.columns._values)
        if any(type(name) is not str for name in names) or len(set(names)) != len(names): raise ValueError('ambiguous columns')
        arrays = raw_columns(df)
        raw = dict(zip(names, arrays))
        cached = {}
        for name, value in df._item_cache.items():
            if type(value) is pd.Series: guards.series_instance_guard(value)
            if (type(name) is not str or type(value) is not pd.Series or value.index is not df.index
                    or type(value.name) is not str or value.name != name):
                raise ValueError('custom cached column')
            if not guards.plain(value.attrs): raise ValueError('custom cached metadata')
            if len(value._mgr.blocks) != 1: raise ValueError('custom cached manager')
            data = value._mgr.blocks[0].values
            if type(data) not in (np.ndarray, pd.arrays.StringArray, pd.arrays.BooleanArray): raise ValueError('custom cached storage')
            cached[name] = data
        effective = lambda name: cached.get(name, raw[name])
        count = len(df)
        # The BED list selection reads raw chrom/pos/id, whereas .get('ref')
        # and later Series expressions deliberately respect existing caches.
        payload = dict(chrom=strings(raw['chrom']), pos=integer(raw['pos'], count),
            variant_id=strings(raw['variant_id']), ref=strings(effective('ref')))
        for position, reference in zip(payload['pos'], payload['ref']):
            if position <= np.iinfo(np.int64).min or position - 1 + len(reference) > np.iinfo(np.int64).max:
                raise ValueError('BED int64 overflow requires original NumPy')
        if whole:
            logger_guard(self.reference.logger)
            payload['effective_id'] = strings(effective('variant_id'))
            payload['DOMAINS'] = strings(effective('DOMAINS'), default='.')
            payload['Consequence'] = strings(effective('Consequence'))
            payload['NMD'] = strings(effective('NMD'), default='.')
            payload['LoF_filter'] = strings(effective('LoF_filter'), default='.')
            for name in FLAGS: payload[name] = booleans(effective(name))
            payload['pvs1'] = integer(pvs, count)
            if type(loc) is not np.ndarray or loc.shape != (count,): raise ValueError('unsupported domain vector')
            payload['loc_intol_domain'] = booleans(loc)
        if sha(self.binary) != self.binary_sha256 or not os.access(self.binary, os.X_OK): raise ValueError('native binary identity changed')
        return payload

    def _run(self, procedure, columns, *, overlap=None):
        # UUID scratch does not consume the deterministic pybedtools fixture RNG.
        target = self.scratch / ('pm4-' + uuid.uuid4().hex)
        target.mkdir()
        source, output = target / 'input.tsv', target / 'output.tsv'
        with source.open('x', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(columns)
            count = len(next(iter(columns.values())))
            writer.writerows(zip(*columns.values()))
        command = [str(self.binary), 'pm4-bp3', 'prepare-bed' if procedure == 'prepare-bed' else 'assign', '--input', str(source), '--output', str(output)]
        if overlap is not None:
            path = target / 'overlap.tsv'
            with path.open('x', newline='') as handle:
                writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
                writer.writerow(['variant_id']); writer.writerows((name,) for name in sorted(overlap))
            command += ['--overlap', str(path)]
        self.events.append(dict(event='native_attempt', procedure=procedure, argv=command))
        completed = subprocess.run(command, capture_output=True, timeout=120, check=False)
        (target / 'stdout.log').write_bytes(completed.stdout)
        (target / 'stderr.log').write_bytes(completed.stderr)
        if completed.returncode != 0: raise ValueError('native transport failure')
        if sha(self.binary) != self.binary_sha256: raise ValueError('native binary changed during execution')
        with output.open(newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            rows = list(reader); names = reader.fieldnames
        if len(rows) != count or any(row.get('ordinal') != str(i) or None in row for i, row in enumerate(rows)):
            raise ValueError('native row/ordinal contract failure')
        self.events.append(dict(event='native_result', procedure=procedure, count=count))
        return names, rows

    def _overlap(self, df, regions_file, method, payload):
        import pybedtools
        # Preserve the original selection and original-frame ref cache timing.
        selected = df[['chrom', 'pos', 'variant_id']]
        df.get('ref', 'A')
        fields = dict(ordinal=list(range(len(df))), chrom=payload['chrom'], pos=payload['pos'],
                      variant_id=payload['variant_id'], ref=payload['ref'])
        try:
            names, rows = self._run('prepare-bed', fields)
            if names != ['ordinal', 'chrom', 'start', 'end', 'variant_id']: raise ValueError('native BED schema')
            if any(row['chrom'] != payload['chrom'][i] or row['variant_id'] != payload['variant_id'][i] for i, row in enumerate(rows)):
                raise ValueError('native BED identity')
            starts = np.array([int(row['start']) for row in rows], dtype=np.int64)
            ends = self._bed_end(df, np.array([int(row['end']) for row in rows], dtype=np.int64))
        except (ValueError, OSError, subprocess.SubprocessError, OverflowError):
            self.events.append(dict(event='original_overlap_before_bedtools'))
            return self.reference.find_overlaps_bedtools_efficient(df, regions_file, method=method)
        bed = selected.assign(start=starts, end=ends)[['chrom', 'start', 'end', 'variant_id']]
        variants_bed = pybedtools.BedTool.from_dataframe(bed)
        regions_bed = pybedtools.BedTool(regions_file)
        self.events.append(dict(event='intersection_started', method=method))
        if method == 'any': intersected = variants_bed.intersect(regions_bed, wa=True, u=True)
        else: intersected = variants_bed.intersect(regions_bed, wa=True, f=1, u=True)
        overlapping = set()
        for feature in intersected: overlapping.add(feature[3])
        self.events.append(dict(event='intersection_completed', unique_string_ids=len(overlapping)))
        return overlapping

    def find_overlaps_bedtools_efficient(self, variants_df, regions_file, method='any'):
        try:
            if type(method) is not str or method not in ('any', 'all'): raise ValueError('unsupported method')
            if type(regions_file) is not str: raise ValueError('custom regions path')
            payload = self._admit(variants_df)
        except (ValueError, TypeError, KeyError, AttributeError, OSError):
            self.events.append(dict(event='original_overlap_pre_admission'))
            return self.reference.find_overlaps_bedtools_efficient(variants_df, regions_file, method=method)
        return self._overlap(variants_df, regions_file, method, payload)

    def PM4_BP3_criteria(self, df, pvs1_criteria, repeat_regions_file, loc_intol_domain):
        try:
            if type(repeat_regions_file) is not str: raise ValueError('custom regions path')
            payload = self._admit(df, whole=True, pvs=pvs1_criteria, loc=loc_intol_domain)
            # Prepare the exact continuation only after admission, but before
            # starting any mature-tool effect. The default path reads no source.
            if self._tail is None: self._tail = self._original_tail()
            if self._emitters is None: self._emitters = self._original_log_sites()
        except (ValueError, TypeError, KeyError, AttributeError, OSError):
            self.events.append(dict(event='original_whole_pre_admission'))
            return self.reference.PM4_BP3_criteria(df, pvs1_criteria, repeat_regions_file, loc_intol_domain)
        overlapping = self._overlap(df, repeat_regions_file, 'all', payload)
        fields = dict(ordinal=list(range(len(df))), variant_id=payload['effective_id'], DOMAINS=payload['DOMAINS'],
            Consequence=payload['Consequence'], NMD=payload['NMD'], LoF_filter=payload['LoF_filter'])
        for name in FLAGS: fields[name] = [int(value) for value in payload[name]]
        fields['pvs1'] = payload['pvs1']
        fields['loc_intol_domain'] = [int(value) for value in payload['loc_intol_domain']]
        try:
            names, rows = self._run('pm4-bp3', fields, overlap=overlapping)
            if names != ['ordinal', *MASKS, 'pm4_array', 'bp3_array']: raise ValueError('native mask schema')
            masks = {name: np.array([row[name] == '1' for row in rows], dtype=np.bool_) for name in MASKS}
            if any(row[name] not in ('0', '1') for row in rows for name in MASKS): raise ValueError('native mask value')
            pm4 = np.array([int(row['pm4_array']) for row in rows], dtype=np.int64)
            bp3 = np.array([int(row['bp3_array']) for row in rows], dtype=np.int64)
            if any(value not in (0, 2) for value in pm4) or any(value not in (0, 1) for value in bp3): raise ValueError('native strength value')
        except (ValueError, OSError, subprocess.SubprocessError, OverflowError):
            self.events.append(dict(event='original_tail_after_one_intersection'))
            return self._tail(df, pvs1_criteria, repeat_regions_file, loc_intol_domain, overlapping)
        self.last_masks = masks
        both = (pm4 > 0) & (bp3 > 0)
        log_arguments = [masks[name] for name in ('in_repetitive_region', 'all_inframe', 'stop_loss_variants',
            'nmd_escaping_trunc', 'has_functional_domain')] + [pm4, bp3, both]
        def emit(position): self._emitters[position][0](*log_arguments)
        ids = df['variant_id']
        # This is the exact original Series.isin publication constructor/finalizer;
        # Rust supplies its boolean result, preserving name/index/attrs/flags.
        in_repeat = ids._constructor(masks['in_repeat'], index=ids.index, copy=False).__finalize__(ids, method='isin')
        df['DOMAINS']
        emit(0)
        for name in ('Consequence', 'Consequence', 'splicing_len_changing', 'splicing_frameshift',
                     '5UTR_len_changing', '5UTR_frameshift', 'vep_consq_len_changing', 'Consequence', 'Consequence'): df[name]
        emit(1)
        df['Consequence']
        emit(2)
        for name in ('Consequence', 'Consequence', 'NMD', 'LoF_filter'): df[name]
        emit(3)
        for name in ('splicing_span_intol_domain', '5UTR_span_intol_domain', 'span_functional_domains'): df[name]
        emit(4)
        df['Consequence']
        emit(5)
        emit(6)
        emit(8 if both.any() else 7)
        self.events.append(dict(event='native_whole_published', rows=len(df)))
        return pm4, bp3, in_repeat

    def _bed_end(self, df, values):
        # The original StringArray.str.len() returns nullable Int64 even without
        # missing values. Preserve that mature-parser input representation.
        source = df._item_cache['ref']._mgr.blocks[0].values
        if type(source) is pd.arrays.StringArray:
            from pandas.core.arrays.integer import IntegerArray
            if pd.arrays.IntegerArray is not IntegerArray:
                raise ValueError('nullable integer class binding changed')
            return pd.arrays.IntegerArray(values, np.zeros(len(values), dtype=np.bool_))
        return values


_PACKAGE_BINDING = Path(__file__).with_name('priva_rust_pm4_bp3_bindings.json')
_PACKAGE_SHA256 = '8c8989440cb2af506fd64318d5c67c93d853db20394fef194f092bcd4b7f366f'
_PACKAGE_PROFILES = (
    ('6fe3b164ac2de7ba79cdb49a8ef8808d7714128dd2733cb0099ce1e70df77662',
     _PACKAGE_BINDING, _PACKAGE_SHA256),
    ('6a650810a8735cbdb988f74729973740a337645dd1585fef696e73e2dbc65beb',
     Path(__file__).with_name('priva_rust_pm4_bp3_bindings_pypi.json'),
     '12e5da0d390fed119d3312ea14ed192c00e508e0003f5031eacc5c5c0ab2c479'),
)


def _package_profile():
    """Select one complete, pinned profile from the actual NumPy extension.

    Conda and the CPython 3.12 Linux x86_64 PyPI wheel have different extension
    bytes. Each admits only its own full package/source/closure proof; versions
    alone, a newly observed binary or a mixture of profiles are insufficient.
    The unchanged package guard validates the selected proof before input use.
    """
    from numpy.core import _multiarray_umath as native
    extension = sha(native.__file__)
    matches = [(path, digest) for expected, path, digest in _PACKAGE_PROFILES
               if expected == extension]
    if len(matches) != 1:
        raise ValueError('unqualified NumPy binary profile')
    return matches[0]


def _optional_bridge(ranking_binary):
    """Set up private transport only; failure leaves the original call untouched."""
    if type(ranking_binary) is not str:
        return None, None
    scratch = None
    try:
        binary = Path(ranking_binary).expanduser().resolve(strict=True)
        digest = sha(binary)
        if not os.access(binary, os.X_OK):
            return None, None
        binding, binding_digest = _package_profile()
        # UUID names do not consume pybedtools' tempfile candidate stream.
        scratch = Path(tempfile.gettempdir()) / ('priva-pm4-bp3-' + uuid.uuid4().hex)
        scratch.mkdir(mode=0o700)
        return _Bridge(_reference, binary=binary, binary_sha256=digest,
            scratch=scratch, package_binding=binding,
            package_binding_sha256=binding_digest), scratch
    except (OSError, ValueError, TypeError):
        if scratch is not None:
            shutil.rmtree(scratch, ignore_errors=True)
        return None, None


def find_overlaps_bedtools_efficient(variants_df, regions_file, method='any', *, ranking_binary=None):
    """Find source-faithful overlaps; default and unsupported inputs use Python."""
    if ranking_binary is None:
        return _reference.find_overlaps_bedtools_efficient(variants_df, regions_file, method=method)
    bridge, scratch = _optional_bridge(ranking_binary)
    if bridge is None:
        return _reference.find_overlaps_bedtools_efficient(variants_df, regions_file, method=method)
    try:
        return bridge.find_overlaps_bedtools_efficient(variants_df, regions_file, method=method)
    finally:
        # Cleanup must not replace an actual mature-tool/source exception.
        shutil.rmtree(scratch, ignore_errors=True)


def PM4_BP3_criteria(df, pvs1_criteria, repeat_regions_file, loc_intol_domain, *, ranking_binary=None):
    """Assign PM4/BP3 with an optional tool and exact current-buffer fallback."""
    if ranking_binary is None:
        return _reference.PM4_BP3_criteria(df, pvs1_criteria, repeat_regions_file, loc_intol_domain)
    bridge, scratch = _optional_bridge(ranking_binary)
    if bridge is None:
        return _reference.PM4_BP3_criteria(df, pvs1_criteria, repeat_regions_file, loc_intol_domain)
    try:
        return bridge.PM4_BP3_criteria(df, pvs1_criteria, repeat_regions_file, loc_intol_domain)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
