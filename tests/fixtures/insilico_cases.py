"""Finite, synthetic original-reference cases for complete PP3/BP4/BP7.

27 October 2026 CT: capture actual Python behavior before replacement. No
expected biological result is encoded here. ``list`` enumerates descriptors
without importing pandas/NumPy or invoking the scientific reference. ``build``
is an importable fixture constructor for the separately frozen Slurm capture.
All values describe fabricated annotation rows, never family/model data.
"""
import argparse
from collections import Counter
import json

SCORES = ('PrimateAI', 'am_pathogenicity', 'CADD_phred', 'CADD_reg_phred')
DS = ('SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL')
PP3_COLUMNS = ('Consequence', 'LoF', *SCORES, 'splicing_lof', *DS,
               'am_class', 'vep_consq_lof', '5UTR_lof', 'CLNSIG', 'CLNREVSTAT')
BP7_COLUMNS = ('Consequence', 'splicing_len_changing', 'splicing_lof', 'Conservation')
OMITTED = object()


def specifications():
    result = []
    def add(group, procedure='PP3_BP4_criteria', **changes):
        result.append(dict(id=f'{len(result):04d}-{group}', group=group,
                           procedure=procedure, **changes))
    add('baseline')
    add('baseline', 'BP7_criteria')
    thresholds = [('PrimateAI', .8), ('am_pathogenicity', .564),
                  ('CADD_phred', 15.), ('CADD_phred', 20.), ('CADD_reg_phred', 15.)]
    thresholds += [(name, .1) for name in DS]
    for column, threshold in thresholds:
        for direction in ('below', 'at', 'above'):
            add('threshold', column=column, threshold=threshold, direction=direction)
    for direction in ('below', 'at', 'above'):
        add('threshold', 'BP7_criteria', column='Conservation', threshold=5., direction=direction)
    coercions = ('bad', 'NaN', 'nan', ' NA ', '1,2', '0.8&0.9', 'none',
                 'pdna', 'nan_float', 'empty', 'space', 'spaced_numeric',
                 'overflow_text', 'true', 'positive_infinity', 'negative_infinity')
    for column in SCORES:
        for value in coercions:
            add('score-coercion', column=column, value=value)
    for column in DS:
        for value in ('bad', 'pdna', 'none', 'nan_float', 'empty'):
            add('splice-score-coercion', column=column, value=value)
    for value in ('mixed-substrings', 'stop-and-splice', 'frameshift-and-missense',
                  'capitalized', 'none', 'pdna', 'nan_float', 'integer',
                  'numpy_string', 'string_subclass'):
        add('consequence', value=value)
    for value in ('HC', 'OS', 'LC', 'hc', 'os', 'HC&LC', 'HC,OS', 'OS;LC',
                  '|HC|', ' HC ', 'HCC', 'XHC', 'HC\n', 'HC\n\n', 'none',
                  'pdna', 'integer', 'empty'):
        add('loftee-token', value=value)
    for column in ('am_class', 'CLNSIG', 'CLNREVSTAT'):
        for value in ('none', 'pdna', 'integer', 'nan_float', 'empty', 'mixed-labels'):
            add('text-label', column=column, value=value)
    for value in ('omitted', 'none', 'all-zero', 'all-two', 'all-three', 'all-four',
                  'all-five', 'boolean', 'float-nan', 'float-infinity',
                  'series-aligned', 'series-reversed', 'series-disjoint',
                  'series-short', 'series-duplicate', 'nullable-series',
                  'short-array', 'long-array', 'matrix', 'python-list',
                  'python-tuple', 'object-string', 'scalar-int'):
        add('pvs1', value=value)
    axes = ('empty', 'singleton', 'duplicate-index', 'named-index', 'reversed-index',
            'multi-index', 'integer-columns', 'multi-columns', 'unused-duplicate',
            'input-attrs', 'populated-item-cache', 'nullable-strings', 'categorical-strings')
    for procedure in ('PP3_BP4_criteria', 'BP7_criteria'):
        for value in axes:
            add('axes-storage', procedure, value=value)
    bool_values = ('plain-true', 'plain-false', 'nullable-missing', 'object-none',
                   'object-nan', 'integer-zero', 'integer-one', 'string-false',
                   'string-empty', 'string-x', 'float-nan')
    for procedure, columns in (('PP3_BP4_criteria', ('splicing_lof', 'vep_consq_lof', '5UTR_lof')),
                               ('BP7_criteria', ('splicing_lof', 'splicing_len_changing'))):
        for column in columns:
            for value in bool_values:
                add('boolean-domain', procedure, column=column, value=value)
    for column in PP3_COLUMNS:
        add('missing-column', column=column)
    for column in BP7_COLUMNS:
        add('missing-column', 'BP7_criteria', column=column)
    for procedure, columns in (('PP3_BP4_criteria', ('Consequence', 'LoF', 'PrimateAI', 'CADD_reg_phred', DS[0])),
                               ('BP7_criteria', ('Consequence', 'splicing_lof', 'Conservation'))):
        for column in columns:
            add('duplicate-column', procedure, column=column)
    for length in ('false', 'true', 'missing'):
        for lof in ('false', 'true', 'missing'):
            add('bp7-literal-or', 'BP7_criteria', length=length, lof=lof)
    for value in ('default', 'empty', 'numeric-labels', 'strength-one', 'strength-two',
                  'strength-three', 'string-strengths', 'negative',
                  'none', 'callable', 'mapping-subclass'):
        add('confidence-mapping', value=value)
    for ordinal in (1, 2, 3, 4):
        add('warning-prefix', ordinal=ordinal)
    for ordinal in (1, 2):
        add('info-prefix', ordinal=ordinal)
    for column in ('Consequence', 'LoF', *SCORES, *DS, 'CLNSIG', 'CLNREVSTAT'):
        add('getter-prefix', column=column)
    for column in BP7_COLUMNS:
        add('getter-prefix', 'BP7_criteria', column=column)
    add('warning-mutation')
    add('legacy-os-test-missing-DS')
    return result


def scalar(name, np, pd):
    values = {'none': None, 'pdna': pd.NA, 'nan_float': float('nan'),
              'empty': '', 'space': ' ', 'spaced_numeric': ' 0.81 ',
              'overflow_text': '1e309', 'true': True,
              'positive_infinity': float('inf'), 'negative_infinity': -float('inf'),
              'integer': 7, 'numpy_string': np.str_('missense_variant')}
    if name == 'string_subclass':
        class Text(str): pass
        return Text('missense_variant')
    return values.get(name, name)


def build(spec, event_callback=None):
    """Return fresh mutable arguments and observer directives; invoke no reference."""
    import numpy as np
    import pandas as pd
    size = 7
    data = {
        'Consequence': ['missense_variant', 'splice_region_variant', '5_prime_UTR_variant',
                        'intron_variant', 'stop_gained', 'frameshift_variant', 'synonymous_variant'],
        'LoF': ['', 'OS', '', '', 'HC', 'LC', ''],
        'PrimateAI': [.9, .1, .1, .1, .1, .1, .1],
        'am_pathogenicity': [.9, .1, .1, .1, .1, .1, .1],
        'CADD_phred': [21., 10., 10., 21., 21., 10., 10.],
        'CADD_reg_phred': [10.] * size,
        'splicing_lof': [False, True, False, False, False, False, False],
        'splicing_len_changing': [False, True, False, False, False, False, False],
        'am_class': ['likely_pathogenic', '', '', '', '', '', ''],
        'vep_consq_lof': [False, False, False, False, True, False, False],
        '5UTR_lof': [False, False, True, False, False, False, False],
        'CLNSIG': ['', '', '', '', '', '', ''],
        'CLNREVSTAT': ['', '', '', '', '', '', ''],
        'Conservation': [0.] * size,
    }
    data.update({name: [0., .5, 0., 0., 0., 0., 0.] for name in DS})
    frame = pd.DataFrame(data)
    pvs, confidence, directive = np.zeros(size, dtype=np.int64), OMITTED, {}
    group, value, column = spec['group'], spec.get('value'), spec.get('column')
    # This child fixes fixture sensitivity, not original rules: isolate the
    # scientific signal being varied so a companion signal cannot mask it.
    def neutral(consequence):
        frame['Consequence'] = consequence
        frame['LoF'] = ''
        frame['PrimateAI'] = .1
        frame['am_pathogenicity'] = .1
        frame['am_class'] = ''
        frame['CADD_phred'] = float('nan')
        frame['CADD_reg_phred'] = 10.
        for name in ('splicing_lof', 'vep_consq_lof', '5UTR_lof'):
            frame[name] = False
        for name in DS: frame[name] = 0.
        frame['CLNSIG'], frame['CLNREVSTAT'] = '', ''
    if group == 'threshold':
        if column == 'Conservation':
            neutral('synonymous_variant'); frame['splicing_len_changing'] = False
        elif column in ('PrimateAI', 'am_pathogenicity'):
            neutral('missense_variant')
        elif column in DS:
            neutral('splice_region_variant')
        else:
            neutral('intron_variant')
            if column == 'CADD_reg_phred': frame['CADD_phred'] = 10.
    elif group == 'splice-score-coercion':
        neutral('splice_region_variant')
        for name in DS: frame[name] = float('nan')
    elif group == 'loftee-token':
        neutral('splice_region_variant')
        for name in DS: frame[name] = float('nan')
    elif group == 'consequence':
        neutral(''); frame['PrimateAI'] = .9
    elif group == 'text-label':
        neutral('missense_variant')
        if column == 'am_class':
            frame['am_pathogenicity'] = .9
        else:
            frame['PrimateAI'] = [.9, .1, .9, .1, .9, .1, .1]
            frame['CLNSIG'] = 'benign/pathogenic'
            frame['CLNREVSTAT'] = 'reviewed_by_expert_panel'
    elif group == 'boolean-domain':
        if spec['procedure'] == 'BP7_criteria':
            neutral('synonymous_variant')
            frame['splicing_len_changing'], frame['splicing_lof'] = True, True
        elif column == 'splicing_lof':
            neutral('splice_region_variant')
        elif column == '5UTR_lof':
            neutral('5_prime_UTR_variant')
        else:
            neutral('intron_variant')
    elif group == 'confidence-mapping':
        neutral('missense_variant')
        frame['PrimateAI'] = [.9, .1, .9, .1, .9, .1, .1]
    if group == 'threshold':
        threshold = spec['threshold']
        numeric = threshold if spec['direction'] == 'at' else float(np.nextafter(threshold,
            -np.inf if spec['direction'] == 'below' else np.inf))
        frame[column] = numeric
    elif group in ('score-coercion', 'splice-score-coercion', 'loftee-token'):
        frame[column or 'LoF'] = pd.Series([scalar(value, np, pd)] * size, dtype=object)
    elif group == 'consequence':
        special = {'mixed-substrings': 'not_missense_variant_extra', 'stop-and-splice': 'stop_gained&splice_region_variant',
                   'frameshift-and-missense': 'missense_variant&frameshift_variant', 'capitalized': 'Missense_variant'}
        frame['Consequence'] = pd.Series([special.get(value, scalar(value, np, pd))] * size, dtype=object)
    elif group == 'text-label':
        actual = ('benign/pathogenic' if column == 'CLNSIG' else
                  'criteria_provided,_multiple_submitters,_no_conflicts' if column == 'CLNREVSTAT' else
                  'benign_pathogenic') if value == 'mixed-labels' else scalar(value, np, pd)
        frame[column] = pd.Series([actual] * size, dtype=object)
    elif group == 'pvs1':
        if value == 'omitted': pvs = OMITTED
        elif value == 'none': pvs = None
        elif value.startswith('all-'):
            pvs = np.full(size, dict(zero=0, two=2, three=3, four=4, five=5)[value[4:]], dtype=np.int64)
        elif value == 'boolean': pvs = np.array([True, False, True, False, True, False, True])
        elif value == 'float-nan': pvs = np.array([0., np.nan, 3., 4., -np.inf, np.inf, 2.])
        elif value == 'float-infinity': pvs = np.full(size, np.inf)
        elif value.startswith('series-'):
            indexes = {'aligned': frame.index, 'reversed': frame.index[::-1], 'disjoint': range(10, 17),
                       'short': range(6), 'duplicate': [0, 0, 1, 2, 3, 4, 5]}
            index = indexes[value[7:]]
            pvs = pd.Series(([0, 3, 2, 0, 4, 1, 5])[:len(index)], index=index)
        elif value == 'nullable-series': pvs = pd.Series([0, 3, pd.NA, 0, 4, 1, 5], dtype='Int64')
        elif value == 'short-array': pvs = np.zeros(6, dtype=np.int64)
        elif value == 'long-array': pvs = np.zeros(8, dtype=np.int64)
        elif value == 'matrix': pvs = np.zeros((size, 1), dtype=np.int64)
        elif value == 'python-list': pvs = [0] * size
        elif value == 'python-tuple': pvs = (0,) * size
        elif value == 'object-string': pvs = np.array(['0'] * size, dtype=object)
        elif value == 'scalar-int': pvs = 3
        else: raise ValueError('unknown PVS fixture')
    elif group == 'axes-storage':
        if value == 'empty': frame = frame.iloc[:0].copy(); pvs = np.zeros(0, dtype=np.int64)
        elif value == 'singleton': frame = frame.iloc[:1].copy(); pvs = np.zeros(1, dtype=np.int64)
        elif value == 'duplicate-index': frame.index = [0, 0, 1, 1, 2, 2, 3]
        elif value == 'named-index': frame.index = pd.Index(['a', 'b', 'c', 'd', 'e', 'f', 'g'], name='synthetic')
        elif value == 'reversed-index': frame.index = frame.index[::-1]
        elif value == 'multi-index': frame.index = pd.MultiIndex.from_tuples([('x', i) for i in range(size)], names=['a', 'b'])
        elif value == 'integer-columns': frame.columns = list(range(len(frame.columns)))
        elif value == 'multi-columns': frame.columns = pd.MultiIndex.from_tuples([(name, 'x') for name in frame.columns])
        elif value == 'unused-duplicate': frame = pd.concat([frame, pd.DataFrame({'unused': range(size)}), pd.DataFrame({'unused': range(size)})], axis=1)
        elif value == 'input-attrs': frame.attrs = {'synthetic': ['retained', {'phase': 1}]}
        elif value == 'populated-item-cache':
            for name in PP3_COLUMNS: frame[name]
        elif value == 'nullable-strings':
            for name in ('Consequence', 'LoF', 'am_class', 'CLNSIG', 'CLNREVSTAT'): frame[name] = frame[name].astype('string')
        elif value == 'categorical-strings': frame['Consequence'] = frame['Consequence'].astype('category')
        else: raise ValueError('unknown axis fixture')
    elif group == 'boolean-domain':
        entries = {'plain-true': (True, bool), 'plain-false': (False, bool),
                   'nullable-missing': (pd.NA, 'boolean'), 'object-none': (None, object),
                   'object-nan': (np.nan, object), 'integer-zero': (0, np.int64),
                   'integer-one': (1, np.int64), 'string-false': ('False', object),
                   'string-empty': ('', object), 'string-x': ('x', object), 'float-nan': (np.nan, float)}
        item, dtype = entries[value]
        frame[column] = pd.Series([item] * size, dtype=dtype)
    elif group == 'missing-column': frame = frame.drop(columns=[column])
    elif group == 'duplicate-column': frame = pd.concat([frame, frame[[column]].copy()], axis=1)
    elif group == 'bp7-literal-or':
        # Four source operands, NINE realizable combinations: neither field can
        # simultaneously equal False and be missing. Do not invent sixteen
        # independent truth assignments to dependent source predicates.
        frame['Consequence'] = 'synonymous_variant'
        states = dict(false=False, true=True, missing=None)
        frame['splicing_len_changing'] = states[spec['length']]
        frame['splicing_lof'] = states[spec['lof']]
    elif group == 'confidence-mapping':
        frame['CLNSIG'] = ['benign', 'pathogenic', 'benign/pathogenic', '', 'benign', 'pathogenic', '']
        frame['CLNREVSTAT'] = 'reviewed_by_expert_panel'
        variants = {'default': OMITTED, 'empty': {}, 'numeric-labels': {0: 3},
                    'strength-one': {'reviewed_by_expert_panel': 1},
                    'strength-two': {'reviewed_by_expert_panel': 2},
                    'strength-three': {'reviewed_by_expert_panel': 3},
                    'string-strengths': {'reviewed_by_expert_panel': '3'},
                    'negative': {'reviewed_by_expert_panel': -3}, 'none': None}
        if value == 'callable': confidence = lambda item: 3
        elif value == 'mapping-subclass':
            class Mapping(dict):
                def __missing__(self, key):
                    self.calls.append(key)
                    if event_callback is not None:
                        event_callback(dict(event='confidence-missing', key=key))
                    return 2
            confidence = Mapping(reviewed_by_expert_panel=3)
            confidence.calls = []
            frame['CLNREVSTAT'] = ['reviewed_by_expert_panel', 'unlisted-review',
                'reviewed_by_expert_panel', 'unlisted-review', 'reviewed_by_expert_panel',
                'unlisted-review', 'unlisted-review']
        else: confidence = variants[value]
    elif group in ('warning-prefix', 'info-prefix', 'getter-prefix', 'warning-mutation'):
        directive = dict(spec)
        if group in ('warning-prefix', 'warning-mutation'):
            for name in SCORES: frame[name] = ['bad'] * size
        if group == 'info-prefix': pvs = np.full(size, 3, dtype=np.int64)
    elif group == 'legacy-os-test-missing-DS': frame = frame.drop(columns=list(DS))
    elif group != 'baseline': raise ValueError('unknown case group')
    return dict(frame=frame, pvs=pvs, confidence=confidence, directive=directive)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    listing = sub.add_parser('list', help='Describe finite fixtures without scientific execution')
    listing.add_argument('--include-descriptors', action='store_true')
    args = parser.parse_args()
    rows = specifications()
    report = dict(cases=len(rows), by_group=dict(Counter(row['group'] for row in rows)),
                  by_procedure=dict(Counter(row['procedure'] for row in rows)),
                  expectations='Actual reference outcomes only; no predicted return/error is asserted')
    if args.include_descriptors: report['specifications'] = rows
    print(json.dumps(report, sort_keys=True))


if __name__ == '__main__': main()
