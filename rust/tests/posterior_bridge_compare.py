#!/usr/bin/env python3
"""Compare the optional posterior procedure and the literal caller tail.

prepare --candidate ROOT --reference ROOT --output NEW_DIR
numerical --candidate ROOT --reference ROOT --inputs DIR --binary PATH --output NEW_DIR
worker --candidate ROOT --reference ROOT --inputs DIR --binary PATH --old-binary PATH
       --mode default|native|fallback|old-command --context fresh|source-bootstrap --output NEW_DIR
compare --inputs DIR --workers DIR --output NEW.json

Only already-assigned criteria and the post-matrix-write caller seam execute.
The original scalar, default model, pandas publications and failure behavior
remain the oracle. No source file, external cache, model or runtime config is
modified. Every path is explicit; generated fixtures belong outside the repo.
"""
import argparse
import ast
import contextlib
import copy
import csv
import dis
from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
import logging
import math
import os
from pathlib import Path
import struct
import subprocess
import sys
import traceback
from types import CodeType
from unittest.mock import patch
import warnings

import numpy as np
import pandas as pd


CRITERIA = ('PVS1', 'PS1', 'PS2', 'PS3', 'PS4', 'PM1', 'PM2', 'PM3', 'PM4',
            'PM5', 'PM6', 'PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'BA1', 'BS1',
            'BS2', 'BS3', 'BS4', 'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7')
AUDITS = ('audit_BS1_disease_incidence_gate', 'audit_BS1_gene_pathogenic_AF_gate',
          'audit_BS1_gene_pathogenic_AF_blocked_by_PM2',
          'audit_PM2_backdown_for_BS1_disease_incidence', 'frequency_conflict_PM2_BS1')
SCALAR_SHA = '64ec17b5dd7d28a2404dedbe6c459709920d93e051a9e4bb6c273a7e7224320d'
RUST_SCALAR_SHA = '24367a7460ccbddbd7475861acbd881b77d048654fe3026f315946d1998f0d04'
REFERENCE_SHA = '9111a76e3fc94a064e14d97b535bb1dd2132f7c61836f8e8892477d4c5e69646'
CALLER_AST = 'fba4729cf524acc543df0d45a493007d7dc7f1e39ffb1bb8be4c474fac8330be'
MODES = ('default', 'native', 'fallback', 'old-command')
CONTEXTS = ('fresh', 'source-bootstrap')
CALLBACK_ACTIVE = False
CALLBACK_EVENTS = []


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def write_json(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False); stream.write('\n')


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def definition(text, name):
    nodes = [node for node in ast.parse(text).body
             if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name]
    if len(nodes) != 1 or type(nodes[0]) is not ast.FunctionDef or nodes[0].decorator_list:
        raise AssertionError('Expected one undecorated synchronous ' + name)
    return nodes[0]


def source_proof(reference, candidate):
    original = (reference / 'scripts/acmg_scoring.py').read_bytes()
    current = (candidate / 'scripts/acmg_scoring.py').read_bytes()
    assert digest(original) == REFERENCE_SHA
    left = ast.get_source_segment(original.decode(), definition(original, 'calculate_posterior_probability'))
    right = ast.get_source_segment(current.decode(), definition(current, 'calculate_posterior_probability'))
    assert left == right and digest(right.encode()) == SCALAR_SHA
    assert digest((candidate / 'rust/src/acmg_scoring/calculate_posterior_probability.rs').read_bytes()) == RUST_SCALAR_SHA
    contract = load_module('posterior_caller_source_contract', candidate / 'tests/test_rust_ranking_caller.py')
    caller = contract.caller_function(contract.acmg_tree())
    restored = contract.historical_caller(caller)
    assert contract.ast_digest(restored) == CALLER_AST
    original_caller = definition((reference / 'scripts/acmg_criteria_assign.py').read_text(), 'ACMG_criteria_assign')
    assert ast.dump(original_caller, include_attributes=False) == ast.dump(restored, include_attributes=False)
    return dict(original_whole_scoring_sha256=REFERENCE_SHA, current_whole_scoring_sha256=digest(current),
        identical_original_posterior_segment_sha256=SCALAR_SHA, unchanged_Rust_scalar_sha256=RUST_SCALAR_SHA,
        restored_whole_caller_AST=CALLER_AST,
        current_whole_caller_sha256=digest((candidate / 'scripts/acmg_criteria_assign.py').read_bytes()),
        no_earlier_stage_or_SciPy_execution=True)


def count_vectors(slots, categories=4):
    if not categories:
        yield ()
    else:
        for count in range(slots + 1):
            for rest in count_vectors(slots - count, categories - 1):
                yield (count,) + rest


def multiplicity(counts, slots):
    ignored = slots - sum(counts)
    result = math.factorial(slots) * 2 ** ignored
    for count in (*counts, ignored):
        result //= math.factorial(count)
    return result


def reduction_ledger(source):
    """Derive the complete count reduction from the pinned scalar AST.

    Category members have identical dyadic weights. All intermediate sums are
    exact multiples of 1/8 below 24, so category order/slot permutations cannot
    change arithmetic. Ignored categories retain both distinct source values.
    Multiplicities account for every one of the 6**28 canonical assignments.
    """
    text = source.read_text(); node = definition(text, 'calculate_posterior_probability')
    assert digest(ast.get_source_segment(text, node).encode()) == SCALAR_SHA
    assignments = {item.targets[0].id: item for item in node.body
                   if isinstance(item, ast.Assign) and len(item.targets) == 1 and isinstance(item.targets[0], ast.Name)}
    groups = {}
    for name in ('patho_matches', 'benign_matches'):
        expression = assignments[name].value
        expected_prefix = 'P' if name == 'patho_matches' else 'B'
        expected = ast.parse('[n for c,n in criteria_assignment.items() if c.startswith(%r) and int(n) > 0]' % expected_prefix,
                             mode='eval').body
        assert ast.dump(expression) == ast.dump(expected)
        groups[name] = tuple(key for key in CRITERIA if key.startswith(expected_prefix))
    selectors = {}
    for name, parent in (('ps_matches', 'patho_matches'), ('pm_matches', 'patho_matches'),
                         ('pp_matches', 'patho_matches'), ('pvs_matches', 'patho_matches'),
                         ('bs_matches', 'benign_matches'), ('bp_matches', 'benign_matches'),
                         ('ba_matches', 'benign_matches'), ('bm_matches', 'benign_matches')):
        expression = assignments[name].value
        assert isinstance(expression, ast.ListComp) and len(expression.generators) == 1
        generator = expression.generators[0]
        assert isinstance(generator.iter, ast.Name) and generator.iter.id == parent
        test = generator.ifs[0]
        assert (len(generator.ifs) == 1 and isinstance(test, ast.Compare)
                and len(test.ops) == 1 and isinstance(test.ops[0], ast.Eq))
        assert isinstance(test.comparators[0], ast.Constant)
        value = test.comparators[0].value
        assert type(value) is int and value in range(1, 6)
        selectors[name] = value
    assert selectors == dict(ps_matches=3, pm_matches=2, pp_matches=1, pvs_matches=4,
                             bs_matches=3, bp_matches=1, ba_matches=5, bm_matches=2)
    namespace = dict(exp_base=2)
    weights = copy.deepcopy(assignments['evidence_weights'])
    exec(compile(ast.Module(body=[weights], type_ignores=[]), str(source), 'exec'), namespace)
    dyadic = {key: Fraction(value) for key, value in namespace['evidence_weights'].items()}
    assert dyadic == {4: Fraction(1), 3: Fraction(1, 2), 2: Fraction(1, 4), 1: Fraction(1, 8), 5: Fraction(2)}
    assert len(groups['patho_matches']) == 16 and len(groups['benign_matches']) == 12
    assert set(groups['patho_matches']).isdisjoint(groups['benign_matches'])
    assert set(groups['patho_matches']) | set(groups['benign_matches']) == set(CRITERIA)
    pathogenic, benign = list(count_vectors(16)), list(count_vectors(12))
    assert len(pathogenic) == 4845 and len(benign) == 1820
    assert sum(multiplicity(counts, 16) for counts in pathogenic) == 6 ** 16
    assert sum(multiplicity(counts, 12) for counts in benign) == 6 ** 12
    no_ba = [counts for counts in benign if counts[3] == 0]
    assert len(no_ba) == 455
    assert sum(multiplicity(counts, 12) for counts in no_ba) == 5 ** 12
    by_exponent = {}
    for counts in no_ba:
        exponent = counts[0] + 2 * counts[1] + 4 * counts[2]
        by_exponent.setdefault(exponent, []).append(counts)
    assert sorted(by_exponent) == list(range(47)) + [48]
    # Every category tuple is emitted, not inferred solely from binomial counts.
    records = dict(pathogenic=[dict(counts=list(counts), multiplicity=str(multiplicity(counts, 16))) for counts in pathogenic],
        benign=[dict(counts=list(counts), multiplicity=str(multiplicity(counts, 12))) for counts in benign],
        no_BA_exponents=[dict(eighths=exponent, representative=list(values[0]),
                              category_vectors=[list(counts) for counts in values])
                         for exponent, values in sorted(by_exponent.items())])
    return dict(source_sha256=SCALAR_SHA, source_groups={key: list(value) for key, value in groups.items()},
        source_selectors=selectors, source_dyadic_weights={str(key): str(value) for key, value in dyadic.items()},
        category_order=dict(pathogenic=[1, 2, 3, 4], benign=[1, 2, 3, 5]),
        ignored_source_values=dict(pathogenic=[0, 5], benign=[0, 4]),
        pathogenic_category_vectors=4845, benign_category_vectors=1820,
        joint_category_vectors=4845 * 1820, BA_containing_category_vectors=4845 * (1820 - 455),
        no_BA_category_vectors=4845 * 455, no_BA_numeric_representatives=4845 * 48,
        canonical_strength_assignments=str(6 ** 28), no_BA_strength_assignments=str(6 ** 16 * 5 ** 12),
        exact_dyadic_sums_bound='0..24 in increments of1/8, each exactly representable in binary64',
        BA_intermediate_finiteness_requires_separate_executed_component_and_interval_proof=True, records=records)


class CallbackStrength:
    def __init__(self, fail=False): self.fail = fail

    def __repr__(self):
        if CALLBACK_ACTIVE: CALLBACK_EVENTS.append('string')
        return 'CallbackStrength(fail=' + str(self.fail) + ')'

    __str__ = __repr__

    def __int__(self):
        if CALLBACK_ACTIVE:
            CALLBACK_EVENTS.append('int')
            warnings.warn('posterior strength conversion', UserWarning)
            if self.fail: raise RuntimeError('posterior strength conversion failed')
        return 1


class CallbackIntegerArray(pd.arrays.IntegerArray):
    @property
    def dtype(self):
        if CALLBACK_ACTIVE:
            CALLBACK_EVENTS.append('extension dtype')
            if object.__getattribute__(self, '__dict__').get('posterior_fail', False):
                warnings.warn('posterior extension dtype', UserWarning)
                raise RuntimeError('posterior extension dtype failed')
        return super().dtype


class PublicationFrame(pd.DataFrame):
    def __setitem__(self, name, value):
        target = object.__getattribute__(self, '__dict__').get('posterior_fail_column')
        if target is not None and name == target:
            warnings.warn('posterior publication ' + name, UserWarning)
            raise RuntimeError('posterior publication failed: ' + name)
        return super().__setitem__(name, value)


class WriteFailureFrame(pd.DataFrame):
    def to_csv(self, *args, **kwargs):
        raise RuntimeError('posterior matrix write failed')


def cases():
    result = [dict(name='integer_canonical', dtype='int64', native=True),
              dict(name='actual_production_int_audits', dtype='int64', production=True, native=True),
              dict(name='ignored_P5_B4', dtype='int64', pattern='ignored', native=True),
              dict(name='all_BA_builtin_integer', dtype='int64', pattern='all_BA', native=True),
              dict(name='mixed_BA_float_publication', dtype='int64', pattern='mixed_BA', native=True),
              dict(name='PS_and_PVS_PM_branch_order', dtype='int64', pattern='branches', native=True)]
    for dtype in ('object', 'bool', 'boolean', 'int32', 'uint64', 'float32', 'float64',
                  'Int32', 'Int64', 'UInt64', 'Float32', 'Float64'):
        result.append(dict(name='ordinary_' + dtype, dtype=dtype, native=True))
    for name in ('duplicate_row_index', 'nonmonotonic_row_index', 'negative_zero', 'warm_column_cache', 'stale_column_cache',
                 'false_duplicate_labels_flag', 'score_publication_failure', 'class_publication_failure', 'frame_row_count_mismatch'):
        result.append(dict(name=name, dtype='object' if name.endswith('cache') else 'float64', native=True))
    for name in ('empty', 'no_columns', 'missing_criterion', 'arbitrary_ignored_column', 'arbitrary_pathogenic_column',
                 'duplicate_column', 'reordered_columns', 'nonstring_column', 'string_strength', 'fractional_strength',
                 'out_of_domain_strength', 'negative_strength', 'none_strength', 'nullable_missing', 'quiet_NaN',
                 'signaling_NaN', 'positive_infinity', 'negative_infinity', 'callback_strength', 'callback_strength_failure',
                 'extension_callback', 'extension_callback_failure', 'matrix_write_failure', 'nondefault_prior',
                 'changed_apply_method', 'changed_to_dict_method', 'changed_helper_global'):
        result.append(dict(name=name, dtype='object', native=False))
    return result


def build_case(spec):
    global CALLBACK_ACTIVE
    CALLBACK_ACTIVE = False
    name, count = spec['name'], 3
    values = {key: pd.array([0] * count, dtype=spec['dtype']) for key in CRITERIA}
    matrix = pd.DataFrame(values)
    data = PublicationFrame({'untouched': ['beta\ttab\nquote"', 'second', 'third'],
                            'ACMG_quant_score': [7.0, -0.0, float('nan')], 'ACMG_class': ['retained'] * count})
    if spec['dtype'] not in ('bool', 'boolean'):
        matrix.loc[0, 'PS1'] = 3; matrix.loc[1, 'PVS1'] = 4; matrix.loc[1, 'PM1'] = 2
    else:
        matrix.loc[0, 'PS1'] = True; matrix.loc[1, 'PM1'] = True
    if spec.get('production'):
        matrix['variant_id'] = ['a', 'b', 'c']
        for key in AUDITS: matrix[key] = np.array([0, 1, 0], dtype=np.int64)
    if spec.get('pattern') == 'ignored':
        matrix.iloc[:, :16] = 5; matrix.iloc[:, 16:] = 4
    if spec.get('pattern') in ('all_BA', 'mixed_BA'):
        matrix.loc[0, 'BS1'] = 5
        if spec['pattern'] == 'all_BA': matrix['BA1'] = 5
    if spec.get('pattern') == 'branches':
        matrix.loc[0, ['PS1', 'PS2', 'PS3']] = 3
        matrix.loc[1, ['PVS1', 'PS4', 'PM1', 'PM2']] = [4, 4, 2, 2]
        matrix.loc[2, ['PVS1', 'PM1', 'PP1']] = [4, 2, 1]
    if name == 'empty': matrix, data = matrix.iloc[:0], data.iloc[:0]
    elif name == 'no_columns': matrix = matrix.iloc[:, :0]
    elif name == 'missing_criterion': matrix = matrix.drop(columns='PS1')
    elif name == 'arbitrary_ignored_column': matrix['unrelated'] = ['a', 'b', 'c']
    elif name == 'arbitrary_pathogenic_column': matrix['PS_extra'] = [3, 1, 0]
    elif name == 'duplicate_column': matrix = pd.concat([matrix, matrix[['PS1']]], axis=1)
    elif name == 'reordered_columns': matrix = matrix.iloc[:, ::-1]
    elif name == 'nonstring_column': matrix.columns = [17, *CRITERIA[1:]]
    elif name == 'nullable_missing': matrix = matrix.astype('Int64'); matrix.loc[1, 'PS1'] = pd.NA
    elif name in ('quiet_NaN', 'signaling_NaN', 'positive_infinity', 'negative_infinity', 'negative_zero'):
        bits = {'quiet_NaN': 0x7ff8000000000001, 'signaling_NaN': 0x7ff0000000000001,
                'positive_infinity': 0x7ff0000000000000, 'negative_infinity': 0xfff0000000000000,
                'negative_zero': 0x8000000000000000}[name]
        matrix.loc[1, 'PS1'] = np.array([bits], dtype=np.uint64).view(np.float64)[0]
    elif name in ('string_strength', 'fractional_strength', 'out_of_domain_strength', 'negative_strength', 'none_strength'):
        matrix.loc[1, 'PS1'] = {'string_strength': '1', 'fractional_strength': 1.5,
            'out_of_domain_strength': 6, 'negative_strength': -1, 'none_strength': None}[name]
    elif name.startswith('callback_strength'):
        matrix.loc[1, 'PS1'] = CallbackStrength(name.endswith('failure'))
    elif name.startswith('extension_callback'):
        array = CallbackIntegerArray(np.array([0, 1, 0], dtype=np.int64), np.zeros(count, dtype=bool))
        matrix['PS1'] = array
        retained = [block.values for block in matrix._mgr.blocks if type(block.values) is CallbackIntegerArray]
        assert len(retained) == 1
        object.__getattribute__(retained[0], '__dict__')['posterior_fail'] = name.endswith('failure')
    elif name == 'matrix_write_failure': matrix = WriteFailureFrame(matrix)
    elif name.endswith('publication_failure'):
        object.__getattribute__(data, '__dict__')['posterior_fail_column'] = 'ACMG_quant_score' if name.startswith('score') else 'ACMG_class'
    elif name == 'frame_row_count_mismatch': data = data.iloc[:2]
    if name in ('duplicate_row_index', 'nonmonotonic_row_index'):
        matrix.index = pd.Index([7, 7, 2] if name.startswith('duplicate') else [8, 3, 5])
        data.index = matrix.index.copy()
    if name == 'false_duplicate_labels_flag':
        matrix.flags.allows_duplicate_labels = False; data.flags.allows_duplicate_labels = False
    for frame in (matrix, data):
        frame.index.name = 'original rows'; frame.columns.name = 'original columns'
        frame.attrs = {'nested': {'source': ['retain', 7]}}
    if name == 'warm_column_cache': matrix['PS1']
    elif name == 'stale_column_cache': matrix._item_cache['PS1'] = pd.Series([5, 5, 5], index=matrix.index, name='PS1')
    return data, matrix


def scalar(value):
    kind = type(value)
    if value is None: return dict(kind='None')
    if value is pd.NA: return dict(kind='pd.NA')
    if value is pd.NaT: return dict(kind='pd.NaT')
    if kind in (str, bool, int): return dict(kind=kind.__name__, value=value)
    if kind is float: return dict(kind='float', bits=struct.pack('>d', value).hex())
    if kind is bytes: return dict(kind='bytes', hex=value.hex())
    if isinstance(value, np.generic): return dict(kind=kind.__module__ + '.' + kind.__qualname__, dtype=value.dtype.str, bits=value.tobytes().hex())
    if kind in (tuple, list): return dict(kind=kind.__name__, values=[scalar(item) for item in value])
    if kind is dict: return dict(kind='dict', items=[[scalar(key), scalar(item)] for key, item in value.items()])
    return dict(kind='unsupported', type=kind.__module__ + '.' + kind.__qualname__)


def array_state(value):
    if type(value) is np.ndarray:
        if value.dtype.kind == 'O': return dict(dtype=value.dtype.str, shape=list(value.shape), values=[scalar(item) for item in value.flat])
        return dict(dtype=value.dtype.str, shape=list(value.shape), bits=value.tobytes().hex())
    if isinstance(value, (pd.arrays.IntegerArray, pd.arrays.FloatingArray, pd.arrays.BooleanArray)):
        return dict(type=type(value).__qualname__, data=array_state(value._data), mask=array_state(value._mask))
    return dict(type=type(value).__qualname__)


def frame_state(frame):
    """Inspect raw buffers, without item-cache or dtype/conversion callbacks."""
    return dict(type=type(frame).__qualname__, shape=list(frame.shape),
        index=dict(type=type(frame.index).__qualname__, values=[scalar(item) for item in frame.index], names=scalar(tuple(frame.index.names))),
        columns=dict(type=type(frame.columns).__qualname__, values=[scalar(item) for item in frame.columns], names=scalar(tuple(frame.columns.names))),
        attrs=scalar(frame._attrs), allows_duplicate_labels=frame.flags.allows_duplicate_labels,
        blocks=[dict(placement=[int(value) for value in block._mgr_locs], array=array_state(block.values)) for block in frame._mgr.blocks])


def frame_fingerprint(frame):
    return digest(json.dumps(frame_state(frame), sort_keys=True, separators=(',', ':')).encode())


def tail_code(candidate):
    contract = load_module('posterior_current_caller_contract', candidate / 'tests/test_rust_ranking_caller.py')
    caller = contract.caller_function(contract.acmg_tree())
    assert contract.ast_digest(contract.historical_caller(caller)) == CALLER_AST
    dispatch = contract.posterior_dispatch(caller)
    start = caller.body.index(dispatch) - 1
    original_write = ast.parse('criteria_matrix.to_csv(output_matrix, sep="\\t", index=False)').body[0]
    assert ast.dump(caller.body[start]) == ast.dump(original_write)
    body = caller.body[start:start + 4]
    assert [ast.unparse(node.targets[0]) for node in body[2:]] == ["anno_df['ACMG_quant_score']", "anno_df['ACMG_class']"]
    module = ast.Module(body=copy.deepcopy(body), type_ignores=[])
    code = compile(module, str(candidate / 'scripts/acmg_criteria_assign.py'), 'exec')
    instructions = list(dis.get_instructions(code)); stores = {}
    for label, node in zip(('score', 'class'), body[2:]):
        target = node.targets[0]
        span = target.lineno, target.end_lineno, target.col_offset, target.end_col_offset
        selected = [item for item in instructions if item.opname == 'STORE_SUBSCR' and tuple(item.positions) == span]
        assert len(selected) == 1
        stores[selected[0].offset] = label
    assert len(stores) == 2
    return code, stores, digest(ast.dump(module, include_attributes=False).encode())


def prepare(args):
    args.output.mkdir()
    proof = source_proof(args.reference, args.candidate)
    ledger = reduction_ledger(args.candidate / 'scripts/acmg_scoring.py')
    _, stores, tail = tail_code(args.candidate)
    inventory = cases()
    assert len({item['name'] for item in inventory}) == len(inventory)
    write_json(args.output / 'ledger.json', ledger)
    write_json(args.output / 'cases.json', inventory)
    write_json(args.output / 'inputs.json', dict(source_proof=proof, literal_tail_AST=tail, STORE_SUBSCR=stores,
        ledger_sha256=digest((args.output / 'ledger.json').read_bytes()),
        cases_sha256=digest((args.output / 'cases.json').read_bytes()), cases=len(inventory),
        modes=list(MODES), contexts=list(CONTEXTS), outcomes=len(inventory) * 8, comparisons=len(inventory) * 6,
        numerical_representatives=232560, finite_intermediate_source_component_calls=4845 + 1820,
        whole_caller_scope='literal matrix write, current-buffer posterior dispatch/apply-zip, original two separate assignments; earlier stages and ranker not executed'))


def modules(candidate):
    sys.path.insert(0, str(candidate / 'scripts'))
    import acmg_scoring as scoring
    import priva_rust_posterior as bridge
    assert Path(scoring.__file__).resolve() == (candidate / 'scripts/acmg_scoring.py').resolve()
    assert Path(bridge.__file__).resolve() == (candidate / 'scripts/priva_rust_posterior.py').resolve()
    assert scoring.__name__ == 'acmg_scoring' and scoring.__package__ == ''
    assert scoring.logger.name == 'acmg_scoring'
    bridge._source_guard(scoring.calculate_posterior_probability)
    return scoring, bridge


def read_inputs(args):
    proof = source_proof(args.reference, args.candidate)
    inputs = json.loads((args.inputs / 'inputs.json').read_bytes())
    assert inputs['source_proof'] == proof
    for name in ('ledger', 'cases'):
        assert digest((args.inputs / (name + '.json')).read_bytes()) == inputs[name + '_sha256']
    return inputs


def raw_identity(frame):
    """Within-process mutation check, including pre-existing stale cache aliases."""
    result = dict(manager=id(frame._mgr), index=id(frame.index), columns=id(frame.columns),
        attrs=id(frame._attrs), flags=id(frame._flags), cache=id(frame._item_cache),
        frame=frame_fingerprint(frame),
        blocks=[(id(block), id(block.values)) for block in frame._mgr.blocks],
        maps=[None if value is None else (id(value), value.tobytes())
              for value in (frame._mgr._blknos, frame._mgr._blklocs)],
        cache_entries=[(scalar(name), id(value), id(value._mgr), id(value.index),
                        tuple((id(block.values), array_state(block.values)) for block in value._mgr.blocks))
                       for name, value in frame._item_cache.items()])
    return result


def execute_case(scoring, bridge, code, stores, spec, mode, output):
    global CALLBACK_ACTIVE
    data, matrix = build_case(spec)
    before = dict(frame=frame_state(data), matrix=frame_state(matrix))
    name = spec['name']; path = output / (name + '.matrix.tsv')
    selected = None if mode == 'default' else output.binary if mode != 'old-command' else output.old_binary
    logs, publications, calls, plans, processes, boxed_rows, plan_rows = [], [], [], [], [], [], []
    pending = [None]; attempts = [0]
    function = scoring.calculate_posterior_probability
    real_try, real_run = bridge.try_posterior_probability, bridge.subprocess.run
    real_apply, real_dict = pd.DataFrame.apply, pd.Series.to_dict
    actual_write_code = matrix.to_csv.__func__.__code__
    CALLBACK_EVENTS.clear()

    def invoke(command, *args, **kwargs):
        if len(command) > 1 and command[1] == 'posterior-probability':
            if mode == 'fallback':
                result = subprocess.CompletedProcess(command, 73, stdout=b'', stderr=b'explicit posterior fixture rejection')
            else:
                result = real_run(command, *args, **kwargs)
            processes.append(dict(executable=str(command[0]), subcommand=command[1], returncode=result.returncode,
                                  stdout=scalar(result.stdout), stderr=scalar(result.stderr)))
            return result
        return real_run(command, *args, **kwargs)

    def dispatch(current, *args, **kwargs):
        attempts[0] += 1
        original = raw_identity(current)
        plan = real_try(current, *args, **kwargs)
        assert raw_identity(current) == original, 'posterior admission/plan mutated current buffers, caches or aliases'
        plans.append(plan is not None)
        if plan is not None: plan_rows.append([scalar(pair) for pair in plan])
        return plan

    def changed_apply(current, *args, **kwargs):
        CALLBACK_EVENTS.append('changed apply')
        return real_apply(current, *args, **kwargs)

    def changed_dict(current, *args, **kwargs):
        CALLBACK_EVENTS.append('changed to_dict')
        return real_dict(current, *args, **kwargs)

    def changed_sum(values):
        CALLBACK_EVENTS.append('changed scalar sum')
        return sum(values)

    class Log(logging.Handler):
        def emit(self, record): logs.append([record.levelname, record.name, record.getMessage()])

    def local_trace(frame, event, arg):
        if event == 'call':
            frame.f_trace_opcodes = True; frame.f_trace_lines = False
        elif event == 'opcode':
            # An opcode event is BEFORE its instruction; the next same-frame
            # opcode proves the pending STORE completed. Exception clears it.
            if pending[0] is not None:
                publications.append(dict(phase=pending[0], frame=frame_state(data), matrix=frame_state(matrix)))
                pending[0] = None
            if frame.f_lasti in stores: pending[0] = stores[frame.f_lasti]
        elif event == 'exception': pending[0] = None
        return local_trace

    def trace(frame, event, arg):
        return local_trace(frame, event, arg) if frame.f_code is code else None

    def profile(frame, event, arg):
        if event == 'call':
            if frame.f_code is function.__code__:
                calls.append('scalar')
            elif frame.f_code is real_apply.__code__ and frame.f_locals.get('self') is matrix:
                calls.append('apply')
            elif frame.f_code is actual_write_code and (frame.f_locals.get('self') is matrix
                    or (frame.f_locals.get('args') and frame.f_locals['args'][0] is matrix)):
                calls.append('matrix_write')
        elif event == 'return' and frame.f_code is function.__code__:
            boxed_rows.append(dict(row_type=type(frame.f_locals['row']).__qualname__,
                criteria=scalar(frame.f_locals.get('criteria_assignment')), returned=scalar(arg)))

    namespace = dict(anno_df=data, criteria_matrix=matrix, output_matrix=str(path),
                     calculate_posterior_probability=function, ranking_binary=selected)
    error, raw_traceback = None, None
    old_trace, old_profile = sys.gettrace(), sys.getprofile()
    old_err, old_callback = np.geterr(), np.geterrcall()
    callback = lambda message, flag: CALLBACK_EVENTS.append(['numpy seterr', message, flag])
    errmode = 'call' if name in ('signaling_NaN', 'positive_infinity', 'negative_infinity') else 'warn'
    observer_frame = sys._getframe(); previous_opcodes = observer_frame.f_trace_opcodes
    logger = logging.getLogger(); previous_level = logger.level; handler = Log()
    logger.addHandler(handler); logger.setLevel(logging.INFO)
    try:
        with warnings.catch_warnings(record=True) as caught, contextlib.ExitStack() as stack:
            warnings.simplefilter('always')
            stack.enter_context(patch.object(bridge, 'try_posterior_probability', dispatch))
            stack.enter_context(patch.object(bridge.subprocess, 'run', invoke))
            if name == 'nondefault_prior': stack.enter_context(patch.object(function, '__defaults__', (0.2, 2, 350)))
            elif name == 'changed_apply_method': stack.enter_context(patch.object(pd.DataFrame, 'apply', changed_apply))
            elif name == 'changed_to_dict_method': stack.enter_context(patch.object(pd.Series, 'to_dict', changed_dict))
            elif name == 'changed_helper_global': stack.enter_context(patch.dict(scoring.__dict__, {'sum': changed_sum}))
            np.seterr(all=errmode); np.seterrcall(callback)
            expected_err = np.geterr()
            CALLBACK_ACTIVE = True
            observer_frame.f_trace_opcodes = True
            sys.setprofile(profile); sys.settrace(trace)
            try:
                exec(code, namespace)
            except Exception as failure:
                error = dict(type=type(failure).__name__, message=str(failure), args=scalar(failure.args))
                raw_traceback = traceback.format_exc()
            finally:
                observer_frame.f_trace_opcodes = previous_opcodes
                sys.settrace(old_trace); sys.setprofile(old_profile)
                CALLBACK_ACTIVE = False
            warning_events = [[type(item.message).__name__, str(item.message)] for item in caught]
            assert np.geterr() == expected_err and np.geterrcall() is callback
    finally:
        CALLBACK_ACTIVE = False
        np.seterr(**old_err); np.seterrcall(old_callback)
        logger.removeHandler(handler); logger.setLevel(previous_level)
    outcome = dict(error=error, warnings=warning_events, logs=logs, calls=calls,
        boxed_rows=boxed_rows, publications=publications,
        initial=before, final_frame=frame_state(data), final_matrix=frame_state(matrix),
        numpy_error_state_restored=True, callback_events=copy.deepcopy(CALLBACK_EVENTS),
        matrix_write=scalar(path.read_bytes()) if path.exists() else None)
    # Native output intentionally skips the original apply/scalar calls. These
    # invocation counters are checked separately; final row kinds/publications,
    # callbacks, warnings and errors remain full-field comparison outcomes.
    outcome_without_dispatch_counts = dict(outcome)
    outcome_without_dispatch_counts.pop('calls'); outcome_without_dispatch_counts.pop('boxed_rows')
    report = dict(name=name, expected_native=spec['native'], outcome=outcome_without_dispatch_counts,
        actual_original_calls=calls, original_boxed_rows=boxed_rows, actual_native_plan_rows=plan_rows,
        attempts=attempts[0], plans=plans, processes=processes)
    write_json(output / (name + '.json'), report)
    if raw_traceback is not None:
        with (output / (name + '.traceback.txt')).open('x') as stream: stream.write(raw_traceback)
    return report


class WorkerOutput:
    """Keep execution arguments out of each immutable outcome comparison."""
    def __init__(self, root, binary, old_binary): self.root, self.binary, self.old_binary = root, binary, old_binary
    def __truediv__(self, name): return self.root / name


def worker(args):
    inputs = read_inputs(args)
    args.output.mkdir()
    scoring, bridge = modules(args.candidate)
    code, stores, tail = tail_code(args.candidate)
    assert tail == inputs['literal_tail_AST']
    rows = []
    for spec in json.loads((args.inputs / 'cases.json').read_bytes()):
        if args.context == 'source-bootstrap':
            bootstrap = args.output / (spec['name'] + '.bootstrap'); bootstrap.mkdir()
            execute_case(scoring, bridge, code, stores, spec, 'default', WorkerOutput(bootstrap, args.binary, args.old_binary))
        rows.append(execute_case(scoring, bridge, code, stores, spec, args.mode,
                                 WorkerOutput(args.output, args.binary, args.old_binary)))
    read_inputs(args)
    write_json(args.output / 'result.json', dict(mode=args.mode, context=args.context, cases=rows,
        input_manifest_sha256=digest((args.inputs / 'inputs.json').read_bytes()),
        canonical_namespace=dict(module='acmg_scoring', package='', logger='acmg_scoring')))


def compare(args):
    inputs = json.loads((args.inputs / 'inputs.json').read_bytes())
    specs = json.loads((args.inputs / 'cases.json').read_bytes())
    mismatches, comparisons, outcomes, coverage = [], 0, 0, {}
    for context in CONTEXTS:
        reports = {mode: json.loads((args.workers / (context + '-' + mode) / 'result.json').read_bytes()) for mode in MODES}
        baseline = reports['default']['cases']
        for mode, report in reports.items():
            assert report['mode'] == mode and report['context'] == context
            assert len(report['cases']) == len(specs)
            assert report['canonical_namespace'] == dict(module='acmg_scoring', package='', logger='acmg_scoring')
            totals = dict(plans=0, original_apply_calls=0, actual_processes=0, actual_old_command_rejections=0, errors=0)
            for spec, reference, current in zip(specs, baseline, report['cases']):
                assert spec['name'] == reference['name'] == current['name']
                outcomes += 1
                totals['plans'] += sum(current['plans'])
                totals['original_apply_calls'] += current['actual_original_calls'].count('apply')
                totals['actual_processes'] += len(current['processes'])
                totals['errors'] += current['outcome']['error'] is not None
                phases = [item['phase'] for item in current['outcome']['publications']]
                assert phases in ([], ['score'], ['score', 'class'])
                if spec['native']:
                    wanted_phases = [] if spec['name'] in ('score_publication_failure', 'frame_row_count_mismatch') else \
                        ['score'] if spec['name'] == 'class_publication_failure' else ['score', 'class']
                    assert phases == wanted_phases, (spec['name'], mode, phases)
                    assert current['actual_original_calls'][0] == 'matrix_write'
                    if mode != 'native':
                        assert current['actual_original_calls'] == ['matrix_write', 'apply', 'scalar', 'scalar', 'scalar']
                if mode == 'default':
                    assert current['attempts'] == 0 and not current['plans'] and not current['processes']
                elif mode == 'native' and spec['native']:
                    assert current['plans'] == [True], spec['name']
                    assert current['actual_original_calls'].count('apply') == 0
                    assert current['original_boxed_rows'] == []
                    assert current['actual_native_plan_rows'] == [[item['returned'] for item in reference['original_boxed_rows']]]
                else:
                    assert not any(current['plans'])
                    assert current['actual_original_calls'] == reference['actual_original_calls']
                    assert current['original_boxed_rows'] == reference['original_boxed_rows']
                if mode == 'old-command' and spec['native']:
                    assert len(current['processes']) == 1 and current['processes'][0]['returncode'] != 0
                    assert current['processes'][0]['returncode'] != 73
                    totals['actual_old_command_rejections'] += 1
                if mode == 'fallback' and spec['native']:
                    assert len(current['processes']) == 1 and current['processes'][0]['returncode'] == 73
                if mode != 'default':
                    comparisons += 1
                    if current['outcome'] != reference['outcome']:
                        fields = sorted(key for key in set(current['outcome']) | set(reference['outcome'])
                                        if current['outcome'].get(key) != reference['outcome'].get(key))
                        mismatches.append(dict(context=context, mode=mode, case=spec['name'], fields=fields,
                                               expected=reference['outcome'], actual=current['outcome']))
            coverage[context + '-' + mode] = totals
    assert comparisons == inputs['comparisons'] and outcomes == inputs['outcomes']
    result = dict(passed=not mismatches, comparisons=comparisons, outcomes=outcomes,
                  mismatch_count=len(mismatches), mismatches=mismatches, coverage=coverage)
    write_json(args.output, result)
    if mismatches: raise AssertionError('all posterior caller mismatches retained in ' + str(args.output))


def assigned_row(pathogenic, benign):
    left = list(itertools.chain.from_iterable([value] * count for value, count in zip((1, 2, 3, 4), pathogenic)))
    right = list(itertools.chain.from_iterable([value] * count for value, count in zip((1, 2, 3, 5), benign)))
    assert len(left) <= 16 and len(right) <= 12
    return left + [0] * (16 - len(left)) + right + [0] * (12 - len(right))


def typed_probability(value):
    if type(value) is int and value == 0: return 'int', '0'
    assert type(value) is np.float64
    return 'float64', struct.pack('>d', value).hex()


def finite_interval_proof(scoring, ledger):
    """Execute every source component, then bound every BA Cartesian pair.

    This is not BA output sampling. All4845 P and1820 B source count states
    execute; exact rational interval composition covers every intermediate
    operation of all6613425 BA-containing category pairs before int0 override.
    Values are many orders inside binary64's normal finite range. Category
    sums are exact dyadics and every prefix/ignored-value case is in the ledger.
    """
    function = scoring.calculate_posterior_probability
    records, captured = dict(pathogenic=[], benign=[]), [None]
    fields = ('odds_benign', 'ps_odds', 'pv_pm_odds', 'pp_odds', 'odds_path', 'posterior_probability')
    def profile(frame, event, arg):
        if event == 'return' and frame.f_code is function.__code__:
            captured[0] = {key: frame.f_locals[key] for key in fields}
    old_profile, old_err, old_callback = sys.getprofile(), np.geterr(), np.geterrcall()
    try:
        np.seterr(all='raise')
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            for side in ('pathogenic', 'benign'):
                for item in ledger['records'][side]:
                    counts = tuple(item['counts'])
                    strengths = assigned_row(counts, (0, 0, 0, 0)) if side == 'pathogenic' else assigned_row((0, 0, 0, 0), counts)
                    row = pd.Series(strengths, index=CRITERIA, dtype='int64')
                    captured[0] = None
                    sys.setprofile(profile)
                    try: probability, class_ = function(row)
                    finally: sys.setprofile(old_profile)
                    values = captured[0]; assert values is not None
                    for name in fields[:-1]:
                        assert np.isfinite(values[name]) and float(values[name]) >= np.finfo(np.float64).tiny, (side, counts, name)
                    kind, bits = typed_probability(probability)
                    if side == 'benign' and counts[3]: assert (kind, bits, class_) == ('int', '0', 'Benign')
                    else: assert kind == 'float64'
                    records[side].append(dict(counts=list(counts),
                        components={key: scalar(value) for key, value in values.items()},
                        numeric={key: float(values[key]) for key in fields[:-1]}, final_kind=kind, final_bits=bits, final_class=class_))
            assert caught == []
            assert np.geterr() == dict.fromkeys(old_err, 'raise') and np.geterrcall() is old_callback
    finally:
        sys.setprofile(old_profile); np.seterr(**old_err); np.seterrcall(old_callback)
    pvalues = [item['numeric']['odds_path'] for item in records['pathogenic']]
    bvalues = [item['numeric']['odds_benign'] for item in records['benign']]
    pbound = min(pvalues), max(pvalues); bbound = min(bvalues), max(bvalues)
    normal = Fraction(float(np.finfo(np.float64).tiny))
    maximum = Fraction(float(np.finfo(np.float64).max))
    def outward(low, high):
        assert -maximum < low <= high < maximum
        result = Fraction(math.nextafter(float(low), -math.inf)), Fraction(math.nextafter(float(high), math.inf))
        assert -maximum < result[0] <= low <= high <= result[1] < maximum
        return result
    pbound, bbound = tuple(map(Fraction, pbound)), tuple(map(Fraction, bbound))
    x = outward(pbound[0] * bbound[0], pbound[1] * bbound[1])
    assert normal < x[0] <= x[1] < maximum
    tenth = Fraction(0.1)
    numerator = outward(x[0] * tenth, x[1] * tenth)
    subtract = outward(x[0] - 1, x[1] - 1)
    multiply = outward(subtract[0] * tenth, subtract[1] * tenth)
    denominator = outward(multiply[0] + 1, multiply[1] + 1)
    assert denominator[0] > 0
    quotient = outward(numerator[0] / denominator[1], numerator[1] / denominator[0])
    scaled = outward(quotient[0] * 10000, quotient[1] * 10000)
    assert normal < numerator[0] and normal < quotient[0] and normal < scaled[0]
    assert max(numerator[1], quotient[1], scaled[1]) < maximum
    # A finite binary64 x near1 has spacing at least2**-53. Thus x-1 is
    # either exact0 or has magnitude >=2**-53 (Sterbenz's exact subtraction
    # applies on[0.5,2]); outside that interval its magnitude is >=0.5.
    # Its following multiplication by0.1 therefore cannot underflow either.
    subtraction_nonzero = Fraction(1, 2 ** 53)
    assert subtraction_nonzero * tenth > normal
    # Scalar.round(4) scales, applies rint and unscales. A rounded nonzero
    # integer is >=1, so its division by10000 is also normal; zero is exact.
    assert Fraction(1, 10000) > normal
    bounds = {name: [str(low), str(high)] for name, (low, high) in
              dict(pathogenic=pbound, benign=bbound, product=x, numerator=numerator,
                   subtraction=subtract, multiplication=multiply, denominator=denominator,
                   quotient=quotient, round_scale=scaled).items()}
    return dict(passed=True, original_source_component_calls=4845 + 1820,
        original_component_count=dict(pathogenic=len(pvalues), benign=len(bvalues)),
        BA_category_pairs_covered_by_complete_interval_composition=6613425,
        source_component_records=records, outward_exact_rational_bounds=bounds,
        no_numpy_warnings_or_error_state_callbacks=True,
        subtraction_nonzero_lower_bound=str(subtraction_nonzero),
        subtraction_zero_is_exact_and_following_nonzero_product_is_normal=True,
        rounded_nonzero_unscale_lower_bound='1/10000; zero remains exact',
        original_BA_override='builtin int0 after all finite source operations; class Benign')


def numerical(args):
    inputs = read_inputs(args)
    args.output.mkdir()
    scoring, bridge = modules(args.candidate)
    bridge._package_guard()
    ledger = json.loads((args.inputs / 'ledger.json').read_bytes())
    assert reduction_ledger(args.candidate / 'scripts/acmg_scoring.py') == ledger
    rows, positions = [], []
    for pindex, pathogenic in enumerate(ledger['records']['pathogenic']):
        for benign in ledger['records']['no_BA_exponents']:
            rows.append(assigned_row(pathogenic['counts'], benign['representative']))
            positions.append((pindex, benign['eighths']))
    assert len(rows) == inputs['numerical_representatives'] == 232560
    frame = pd.DataFrame(np.asarray(rows, dtype=np.int64), columns=CRITERIA)
    source_bytes_before = frame_fingerprint(frame)
    input_path, output_path = args.output / 'request.csv', args.output / 'native.csv'
    with input_path.open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['source_position', *CRITERIA])
        writer.writerows((position, *row) for position, row in enumerate(rows))
    command = [str(args.binary), 'posterior-probability', '--input', str(input_path), '--output', str(output_path)]
    process = subprocess.run(command, capture_output=True)
    (args.output / 'native.stdout').write_bytes(process.stdout); (args.output / 'native.stderr').write_bytes(process.stderr)
    assert process.returncode == 0
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        expected = frame.apply(scoring.calculate_posterior_probability, axis=1)
    assert caught == [] and frame_fingerprint(frame) == source_bytes_before
    with output_path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == ['source_position', 'probability_kind', 'probability_value', 'ACMG_class']
        actual = list(reader)
    mismatches = []
    with (args.output / 'expected.csv').open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['source_position', 'probability_kind', 'probability_value', 'ACMG_class'])
        for position, (probability, class_) in enumerate(expected):
            kind, value = typed_probability(probability)
            wanted = dict(source_position=str(position), probability_kind=kind, probability_value=value, ACMG_class=class_)
            writer.writerow(wanted.values())
            current = actual[position] if position < len(actual) else None
            if current != wanted:
                mismatches.append(dict(position=position, pathogenic_category_index=positions[position][0],
                    benign_exponent_eighths=positions[position][1], expected=wanted, actual=current))
        if len(actual) != len(expected): mismatches.append(dict(expected_rows=len(expected), actual_rows=len(actual)))
    interval = finite_interval_proof(scoring, ledger)
    write_json(args.output / 'finite-intermediates.json', interval)
    read_inputs(args)
    result = dict(passed=not mismatches, comparisons=len(expected), actual_rows=len(actual),
        scalar_kind_raw_f64_bits_class_and_ordinal_all_protected=True,
        original_source_component_calls=interval['original_source_component_calls'],
        BA_category_pairs_proved_finite_before_int0=6613425,
        input_sha256=digest(input_path.read_bytes()), native_sha256=digest(output_path.read_bytes()),
        expected_sha256=digest((args.output / 'expected.csv').read_bytes()),
        finite_intermediate_proof_sha256=digest((args.output / 'finite-intermediates.json').read_bytes()),
        mismatch_count=len(mismatches), mismatches=mismatches,
        native_argv=command, native_returncode=process.returncode)
    write_json(args.output / 'result.json', result)
    if mismatches: raise AssertionError('all default-domain posterior mismatches retained')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='procedure', required=True)
    for name in ('prepare', 'numerical', 'worker'):
        command = commands.add_parser(name)
        command.add_argument('--candidate', type=Path, required=True)
        command.add_argument('--reference', type=Path, required=True)
        command.add_argument('--output', type=Path, required=True)
        if name != 'prepare':
            command.add_argument('--inputs', type=Path, required=True)
            command.add_argument('--binary', type=Path, required=True)
        if name == 'worker':
            command.add_argument('--old-binary', type=Path, required=True)
            command.add_argument('--mode', choices=MODES, required=True)
            command.add_argument('--context', choices=CONTEXTS, required=True)
    command = commands.add_parser('compare')
    command.add_argument('--inputs', type=Path, required=True)
    command.add_argument('--workers', type=Path, required=True)
    command.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    globals()[args.procedure.replace('-', '_')](args)


if __name__ == '__main__': main()
