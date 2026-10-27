"""Portable, source-derived PM4/BP3 inputs; no expected results are guessed."""
import itertools
import logging
import numpy as np
import pandas as pd

CALLBACKS = []
FLAGS = ('splicing_len_changing', 'splicing_frameshift', '5UTR_len_changing',
         '5UTR_frameshift', 'vep_consq_len_changing', 'splicing_span_intol_domain',
         '5UTR_span_intol_domain', 'span_functional_domains')
REGIONS = b'chr1\t10\t20\nchr1\t30\t34\nchr1\t40\t41\nchr2\t0\t3\n'


class LengthProbe:
    def __init__(self, raises=False): self.raises = raises
    def __len__(self):
        CALLBACKS.append('ref.__len__')
        if self.raises: raise RuntimeError('deliberate ref length callback')
        return 3
    def __repr__(self): return 'LengthProbe(raises=%r)' % self.raises


class AccessFrame(pd.DataFrame):
    @property
    def _constructor(self): return AccessFrame
    def __getitem__(self, key):
        CALLBACKS.append(('frame.__getitem__', key if isinstance(key, str) else tuple(key)))
        return super().__getitem__(key)


class RaiseAtLog(logging.Handler):
    def __init__(self, at): super().__init__(); self.at = at; self.seen = 0
    def emit(self, record):
        self.seen += 1
        CALLBACKS.append(('log', self.seen, record.getMessage()))
        if self.seen == self.at: raise RuntimeError('deliberate logging failure at %d' % self.at)


def base(rows=4):
    positions = [11, 21, 31, 41]
    frame = pd.DataFrame(dict(chrom=['chr1'] * rows,
        pos=np.array([positions[i % 4] for i in range(rows)], dtype=np.int64),
        variant_id=['v%d' % i for i in range(rows)], ref=['A', 'AAA', 'AAAA', 'AA'] * (rows // 4) + ['A', 'AAA', 'AAAA', 'AA'][:rows % 4],
        DOMAINS=['.'] * rows, Consequence=['inframe_deletion'] * rows,
        NMD=['.'] * rows, LoF_filter=['.'] * rows))
    for name in FLAGS: frame[name] = np.zeros(rows, dtype=np.bool_)
    frame.attrs = dict(prototype='PM4/BP3', nested=dict(values=[1, 'x']))
    return frame, np.zeros(rows, dtype=np.int64), np.zeros(rows, dtype=np.bool_)


def descriptors():
    helper = ['overlap-any', 'overlap-all', 'overlap-unknown-method', 'overlap-none-method',
              'overlap-duplicate-ids', 'overlap-integer-ids', 'overlap-null-id', 'overlap-unicode-ref',
              'overlap-empty-ref', 'overlap-missing-ref', 'overlap-pos-zero', 'overlap-pos-negative',
              'overlap-pos-float', 'overlap-missing-regions', 'overlap-invalid-regions',
              'overlap-full-ref-boundary', 'overlap-empty-frame']
    whole = ['whole-base', 'whole-truth-grid', 'whole-duplicate-ids', 'whole-integer-ids',
             'whole-null-id', 'whole-duplicate-index', 'whole-nonmonotonic-index', 'whole-multi-index',
             'whole-nullable-booleans', 'whole-string-python', 'whole-categorical-text',
             'whole-object-bools', 'whole-float-flags', 'whole-int-flags', 'whole-null-domains',
             'whole-numeric-domains', 'whole-mixed-consequence', 'whole-null-consequence',
             'whole-pvs-float', 'whole-pvs-nan', 'whole-pvs-snan-raise', 'whole-pvs-infinity',
             'whole-pvs-short', 'whole-pvs-long', 'whole-pvs-series', 'whole-loc-short', 'whole-loc-int',
             'whole-ref-callback', 'whole-ref-callback-raises', 'whole-frame-getter',
             'whole-warm-cache', 'whole-stale-ref-cache', 'whole-stale-id-cache', 'whole-stale-consequence-cache',
             'whole-missing-ref', 'whole-missing-domains', 'whole-missing-consequence', 'whole-missing-nmd',
             'whole-missing-late-flag', 'whole-duplicate-column', 'whole-empty', 'whole-extra-object',
             'whole-log-fail-1', 'whole-log-fail-3', 'whole-log-fail-6',
             'whole-stoploss-pvs-precedence', 'whole-nmd-pvs-precedence', 'whole-substring-semantics',
             'whole-flags-no-duplicate-labels', 'performance-10000', 'performance-100000']
    return [dict(name=name, procedure='overlap' if name in helper else 'whole',
                 performance=name.startswith('performance-')) for name in helper + whole]


def build(name):
    CALLBACKS.clear()
    frame, pvs, loc = base()
    method, regions, np_errors, log_failure = 'any', REGIONS, None, None
    if name == 'whole-truth-grid':
        rows = []
        consequences = ('inframe_deletion', 'inframe_insertion', 'stop_lost', 'frameshift_variant',
                        'stop_gained', 'inframe_deletion&stop_lost', 'frameshift_variant&stop_lost', 'missense_variant')
        for serial, (consequence, repeated, low_complexity, domain, strength, splice, utr, vep, escape) in enumerate(
            itertools.product(consequences, (False, True), (False, True), (False, True), range(5),
                              (False, True), (False, True), (False, True), (False, True))):
            row = dict(chrom='chr1', pos=11 if repeated else 21, variant_id='grid%d' % serial, ref='A',
                DOMAINS='Low_complexity' if low_complexity else '.', Consequence=consequence,
                NMD='NMD_escaping' if escape else '.', LoF_filter='.',
                splicing_len_changing=splice, splicing_frameshift=False, **{'5UTR_len_changing':utr,
                '5UTR_frameshift':False}, vep_consq_len_changing=vep,
                splicing_span_intol_domain=False, **{'5UTR_span_intol_domain':False}, span_functional_domains=False,
                _strength=strength, _domain=domain)
            rows.append(row)
        frame = pd.DataFrame(rows)
        pvs = frame.pop('_strength').to_numpy(dtype=np.int64)
        loc = frame.pop('_domain').to_numpy(dtype=np.bool_)
    elif name.startswith('performance-'):
        count = int(name.split('-')[1]); frame, pvs, loc = base(count)
        frame['Consequence'] = np.resize(np.array(['inframe_deletion','missense_variant','stop_lost','frameshift_variant'], dtype=object), count)
        frame['DOMAINS'] = np.resize(np.array(['.', 'Low_complexity', '.', '.'], dtype=object), count)
        pvs[:] = np.arange(count) % 5
        loc[:] = (np.arange(count) % 7 == 0)
    elif name.endswith('duplicate-ids'):
        frame['variant_id'] = ['shared', 'shared', 'other', 'other']
    elif name.endswith('integer-ids'): frame['variant_id'] = np.array([1, 2, 3, 4], dtype=np.int64)
    elif name.endswith('null-id'): frame['variant_id'] = ['.', None, 'x', 'x']
    elif name == 'overlap-unicode-ref': frame['ref'] = ['é', '😀AA', 'A\u0301', 'AC']
    elif name == 'overlap-empty-ref': frame['ref'] = ['', 'A', '', 'A']
    elif name in ('overlap-missing-ref', 'whole-missing-ref'): frame = frame.drop(columns='ref')
    elif name == 'overlap-pos-zero': frame.loc[0, 'pos'] = 0
    elif name == 'overlap-pos-negative': frame.loc[0, 'pos'] = -1
    elif name == 'overlap-pos-float': frame['pos'] = frame['pos'].astype(np.float64)
    elif name == 'overlap-missing-regions': regions = None
    elif name == 'overlap-invalid-regions': regions = b'chr1\tbad\t20\n'
    elif name == 'overlap-full-ref-boundary':
        frame['pos'] = [10, 11, 19, 20]; frame['ref'] = ['AAA', 'A' * 10, 'AA', 'AA']
    elif name in ('overlap-empty-frame', 'whole-empty'): frame, pvs, loc = base(0)
    elif name == 'whole-duplicate-index': frame.index = pd.Index(['a', 'a', 'b', 'b'], name='duplicate')
    elif name == 'whole-nonmonotonic-index': frame.index = pd.Index([9, -1, 4, 0], name='source_order')
    elif name == 'whole-multi-index': frame.index = pd.MultiIndex.from_tuples([('b',2),('a',1),('b',2),('x',0)], names=['first','second'])
    elif name == 'whole-nullable-booleans':
        for field in FLAGS: frame[field] = pd.array([False, True, pd.NA, False], dtype='boolean')
    elif name == 'whole-string-python':
        for field in ('chrom','variant_id','ref','DOMAINS','Consequence','NMD','LoF_filter'):
            frame[field] = frame[field].astype('string[python]')
    elif name == 'whole-categorical-text': frame['Consequence'] = pd.Categorical(frame['Consequence'])
    elif name == 'whole-object-bools':
        for field in FLAGS: frame[field] = pd.Series([True, False, None, True], dtype=object)
    elif name == 'whole-float-flags': frame['splicing_frameshift'] = np.array([0., 1., np.nan, 0.])
    elif name == 'whole-int-flags': frame['splicing_len_changing'] = np.array([0, 1, 2, -1], dtype=np.int64)
    elif name == 'whole-null-domains': frame['DOMAINS'] = [None, '.', np.nan, 'Low_complexity']
    elif name == 'whole-numeric-domains': frame['DOMAINS'] = [1, 2, 3, 4]
    elif name == 'whole-mixed-consequence': frame['Consequence'] = ['inframe_deletion', 3, None, 'stop_lost']
    elif name == 'whole-null-consequence': frame['Consequence'] = [None] * 4
    elif name == 'whole-pvs-float': pvs = np.array([0., 1., 2., 3.], dtype=np.float64)
    elif name == 'whole-pvs-nan': pvs = np.array([np.nan, 0., 1., 3.], dtype=np.float64)
    elif name == 'whole-pvs-snan-raise':
        pvs = np.array([0x7ff0000000000001,0,0,0],dtype=np.uint64).view(np.float64); np_errors='raise'
    elif name == 'whole-pvs-infinity': pvs = np.array([np.inf, -np.inf, 0., 1.], dtype=np.float64)
    elif name == 'whole-pvs-short': pvs = np.zeros(2, dtype=np.int64)
    elif name == 'whole-pvs-long': pvs = np.zeros(5, dtype=np.int64)
    elif name == 'whole-pvs-series': pvs = pd.Series([0, 1, 2, 3], index=[3, 2, 1, 0], name='pvs')
    elif name == 'whole-loc-short': loc = np.zeros(2, dtype=np.bool_)
    elif name == 'whole-loc-int': loc = np.array([0, 1, 2, -1], dtype=np.int64)
    elif name in ('whole-ref-callback','whole-ref-callback-raises'):
        frame['ref'] = [LengthProbe(name.endswith('-raises')), 'A', 'AA', 'AAA']
    elif name == 'whole-frame-getter': frame = AccessFrame(frame)
    elif name == 'whole-warm-cache':
        for field in ('variant_id','ref','Consequence','DOMAINS'): frame[field]
    elif name == 'whole-stale-ref-cache':
        frame['ref']; frame._item_cache['ref'] = pd.Series(['A' * 20] * 4, index=frame.index, name='ref')
    elif name == 'whole-stale-id-cache':
        frame['variant_id']; frame._item_cache['variant_id'] = pd.Series(['v1','v0','v3','v2'], index=frame.index, name='variant_id')
    elif name == 'whole-stale-consequence-cache':
        frame['Consequence']; frame._item_cache['Consequence'] = pd.Series(['stop_lost'] * 4, index=frame.index, name='Consequence')
    elif name == 'whole-missing-domains': frame = frame.drop(columns='DOMAINS')
    elif name == 'whole-missing-consequence': frame = frame.drop(columns='Consequence')
    elif name == 'whole-missing-nmd': frame = frame.drop(columns='NMD')
    elif name == 'whole-missing-late-flag': frame = frame.drop(columns='span_functional_domains')
    elif name == 'whole-duplicate-column':
        frame = pd.concat([frame, frame[['Consequence']]], axis=1)
    elif name == 'whole-extra-object': frame['unused'] = [LengthProbe(True)] * 4
    elif name.startswith('whole-log-fail-'): log_failure = int(name.rsplit('-', 1)[1])
    elif name == 'whole-stoploss-pvs-precedence':
        frame['Consequence'] = ['inframe_deletion&stop_lost'] * 4; pvs[:] = [0, 1, 2, 4]
    elif name == 'whole-nmd-pvs-precedence':
        frame['Consequence'] = ['inframe_deletion&stop_gained'] * 4
        frame['NMD'] = ['escaping'] * 4; pvs[:] = [0, 1, 2, 4]
    elif name == 'whole-substring-semantics':
        frame['Consequence'] = ['not_inframe_deletion', 'not_stop_lost', 'frameshift_variant', 'stop_gained']
        frame['NMD'] = ['.', '.', 'not_escaping', '.']; frame['LoF_filter'] = ['.', '.', '.', 'xEND_TRUNCx']
        frame['DOMAINS'] = ['low_complexity','xLow_complexityx','.','.']
    elif name == 'whole-flags-no-duplicate-labels': frame.flags.allows_duplicate_labels = False
    elif name not in ('whole-base','overlap-any','overlap-all','overlap-unknown-method','overlap-none-method'):
        raise ValueError('Unknown source-derived case: ' + name)
    if name == 'overlap-all': method = 'all'
    elif name == 'overlap-unknown-method': method = 'ALL'
    elif name == 'overlap-none-method': method = None
    elif name == 'overlap-full-ref-boundary': method = 'all'
    return dict(frame=frame, pvs=pvs, loc=loc, method=method, regions=regions,
                np_errors=np_errors, log_failure=log_failure)
