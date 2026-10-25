#!/usr/bin/env python3
"""Capture the unchanged 09223ef full-ranking/control oracle before integration.

prepare --repo REPO --output NEW_DATA_ROOT
capture --root DATA_ROOT --pins-sha256 SHA
validate --root DATA_ROOT --capture-sha256 SHA

Use real bounded Slurm steps. All fixtures, pickles, typed states, native-package
identities and logs belong under the supplied data root. No Rust is invoked,
no scientific source is rewritten, and no family/cache data is read. Tracing is
checked against a separate untraced execution for every case. Production parser
excerpts are original AST statements, explicitly not criterion assignment.
"""
import argparse
import ast
import datetime
import hashlib
import importlib
import inspect
import io
import json
import logging
import os
from pathlib import Path
import pickle
import platform
import struct
import subprocess
import sys
import traceback
import warnings

import numpy as np
import pandas as pd
import scipy
from scipy.stats import beta, binom

Y = Path('/stornext/snfs170/next-gen/scratch/yangyxt')
D = Path('/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade')
COMMIT = '09223ef5a8b0ceb4c6d61055c6eefe62b51ba903'
SOURCE_FILES = [
    'acmg_scoring.py', 'acmg_consequence.py', 'acmg_variant_mechanism.py',
    'gene_mechanism_hub.py', 'gene_mechanism_common.py', 'gene_mechanism_conditions.py',
    'gene_mechanism_variants.py', 'gene_mechanism_resources.py', 'gene_mechanism_annotation.py',
    'acmg_inheritance.py', 'hpo_penetrance.py', 'clinvar_vcv.py',
    'resort_step3_only.py', 'acmg_criteria_assign.py',
]
SCIENCE_PINS = {
    'acmg_scoring.py': '9111a76e3fc94a064e14d97b535bb1dd2132f7c61836f8e8892477d4c5e69646',
    'acmg_consequence.py': 'a9efd0964f25e84dc01ad276d4f0905a10492a5d1f94c2b277573ac94244d354',
    'acmg_variant_mechanism.py': 'c7e8da80ca0ec47fb38a5fd1d965641c1791189cb0d25a651253f44508d0bf26',
}
FUNCTIONS = {
    'acmg_scoring.py': ['compute_control_common', 'sort_and_rank_variants'],
    'acmg_consequence.py': ['coerce_pext_value'],
    'acmg_variant_mechanism.py': ['_variant_mechanism_masks'],
    'resort_step3_only.py': ['resort_family'],
    'acmg_criteria_assign.py': ['ACMG_criteria_assign'],
}
CONTROL_COLUMNS = ['control_AC', 'control_AN', 'control_AF', 'control_nhomalt']
CONTROL_LOCALS = ['df', 'an', 'ac', 'af', 'nhomo', 'valid_an', 'gt', 'is_hom',
                  'use_nhomo', 'use_ac', 'use_af', 'common', 'm', 'ac_safe',
                  'an_safe', 'confident_af', 'valid_for_beta', 'n_flagged']
RANK_LOCALS = ['df', 'df_sorted', 'proband_row', 'proband', 'proband_sex',
               'proband_is_male', 'proband_is_female', 'proband_het', 'mechanism_masks',
               'chrom', 'x_locus', 'has_single_alt_lof_compatible', 'exact_gof',
               'asserted_lof_effect', 'dn_without_hi_or_ambiguous', 'heterozygous_lof_incompatible',
               'high_acmg_lof', 'compatibility_penalty_rows', 'pext_modulation_rows',
               'dispensable_genes', 'relevant_genes', 'pext_col_to_use', 'pext_tissue_list',
               'first_tissue_col', 'pext_numeric', 'scaled_pext', 'max_scores', 'codes', 'var_key']


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def native_environment():
    assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_STEP_ID'), 'real Slurm execution required'
    return {name: os.environ.get(name) for name in ['SLURM_JOB_ID', 'SLURM_STEP_ID', 'SLURMD_NODENAME']}


def write_json(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True, allow_nan=False)
        stream.write('\n')


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True, allow_nan=False).encode()


def typename(value):
    return type(value).__module__ + '.' + type(value).__qualname__


def scalar(value):
    """Lossless declared scalars, rather than JSON's conflated null/NaN/zero."""
    if value is None:
        return dict(kind='none')
    if value is pd.NA:
        return dict(kind='pd.NA')
    if value is pd.NaT:
        return dict(kind='pd.NaT')
    if isinstance(value, np.generic):
        return dict(kind='numpy_scalar', type=typename(value), dtype=value.dtype.str,
                    storage_hex=value.tobytes().hex())
    if isinstance(value, bool):
        return dict(kind='bool', value=value)
    if isinstance(value, int):
        return dict(kind='int', decimal=str(value))
    if isinstance(value, float):
        return dict(kind='float', ieee_f64_be=struct.pack('>d', value).hex())
    if isinstance(value, str):
        return dict(kind='str', value=value)
    if isinstance(value, bytes):
        return dict(kind='bytes', hex=value.hex())
    if isinstance(value, Path):
        return dict(kind='path', value=str(value))
    if isinstance(value, pd.Timestamp):
        return dict(kind='timestamp', iso=value.isoformat(), ticks=str(value.value), timezone=str(value.tz))
    return None


def dtype(value):
    result = dict(type=typename(value), text=str(value), repr=repr(value))
    if isinstance(value, pd.CategoricalDtype):
        result.update(categories=encode(value.categories), ordered=value.ordered)
    return result


def array_storage(value):
    if isinstance(value, np.ndarray) and value.dtype.kind != 'O':
        return dict(dtype=value.dtype.str, shape=list(value.shape), c_order_hex=value.tobytes(order='C').hex())
    return None


def missing(value):
    try:
        result = pd.isna(value)
        return bool(result) if np.ndim(result) == 0 else encode(result)
    except Exception as error:
        return dict(unavailable=typename(error), message=str(error))


def encode(value):
    atom = scalar(value)
    if atom is not None:
        return atom
    if isinstance(value, pd.DataFrame):
        return dict(kind='dataframe', type=typename(value), columns=encode(value.columns), index=encode(value.index),
                    shape=list(value.shape), attrs=encode(value.attrs),
                    allows_duplicate_labels=value.flags.allows_duplicate_labels,
                    data=[encode(value.iloc[:, position]) for position in range(value.shape[1])])
    if isinstance(value, pd.Series):
        array = value.array
        storage = array_storage(value.to_numpy(copy=False))
        if hasattr(array, '_data') and hasattr(array, '_mask'):
            storage = dict(data=array_storage(array._data), mask=array_storage(array._mask))
        return dict(kind='series', type=typename(value), name=encode(value.name), index=encode(value.index),
                    dtype=dtype(value.dtype), values=[encode(x) for x in array], null_mask=[missing(x) for x in array],
                    storage=storage, attrs=encode(value.attrs), allows_duplicate_labels=value.flags.allows_duplicate_labels)
    if isinstance(value, pd.Index):
        result = dict(kind='index', type=typename(value), dtype=dtype(value.dtype), names=encode(value.names),
                      values=[encode(x) for x in value], null_mask=[missing(x) for x in value])
        if isinstance(value, pd.RangeIndex):
            result.update(start=value.start, stop=value.stop, step=value.step)
        if isinstance(value, pd.MultiIndex):
            result.update(levels=[encode(x) for x in value.levels], codes=[encode(x) for x in value.codes])
        return result
    if isinstance(value, np.ndarray):
        return dict(kind='ndarray', dtype=dtype(value.dtype), shape=list(value.shape), storage=array_storage(value),
                    values=[encode(x) for x in value.flat])
    if isinstance(value, dict):
        return dict(kind='dict', items=[[encode(key), encode(item)] for key, item in value.items()])
    if isinstance(value, (list, tuple)):
        return dict(kind=type(value).__name__, values=[encode(x) for x in value])
    if isinstance(value, (set, frozenset)):
        items = [encode(x) for x in value]
        return dict(kind=type(value).__name__, values=sorted(items, key=canonical))
    return dict(kind='opaque_diagnostic', type=typename(value), repr=repr(value))


def phase(function, line):
    if function == 'compute_control_common':
        for start, end, name in [(390, 411, 'C1'), (413, 425, 'C2'), (427, 447, 'C3'),
                                 (449, 462, 'C4'), (463, 471, 'C5')]:
            if start <= line <= end:
                return name
        return 'C-boundary'
    if function == 'sort_and_rank_variants':
        for start, end, name in [(474, 498, 'R1'), (501, 524, 'R2'), (527, 528, 'R3'),
                                 (537, 597, 'R4'), (599, 611, 'R5'), (613, 616, 'R6'),
                                 (621, 665, 'R7'), (668, 697, 'R8')]:
            if start <= line <= end:
                return name
        return 'R-boundary'
    return {'coerce_pext_value': 'PEXT-scalar', '_variant_mechanism_masks': 'mechanism-helper',
            'resort_family': 'resort-caller', '<module>': 'production-parser-excerpt'}.get(function, 'SciPy')


class StateStore:
    def __init__(self, root):
        self.root = root
        root.mkdir()

    def put(self, value):
        raw = canonical(value)
        key = hashlib.sha256(raw).hexdigest()
        path = self.root / (key + '.json')
        if not path.exists():
            path.write_bytes(raw + b'\n')
        return key


class Trace:
    """Observe original Python frames; never replace SciPy or source functions."""
    def __init__(self, source, store):
        self.source = source
        self.store = store
        self.events = []
        self.previous = {}
        self.exception_events = []

    def __call__(self, frame, event, arg):
        name = frame.f_code.co_name
        path = Path(frame.f_code.co_filename)
        source = path.parent == self.source and (name in {
            'compute_control_common', 'sort_and_rank_variants', 'coerce_pext_value',
            '_variant_mechanism_masks', 'resort_family', '<module>'})
        scipy_call = name in {'cdf', '_cdf', 'ppf', '_ppf'} and getattr(frame.f_locals.get('self'), 'name', None) in {'binom', 'beta'}
        if event == 'call' and not (source or scipy_call):
            return None
        if not (source or scipy_call):
            return None
        if scipy_call and event == 'line':
            return self
        values = frame.f_locals
        if name == 'compute_control_common':
            selected = CONTROL_LOCALS
        elif name == 'sort_and_rank_variants':
            selected = RANK_LOCALS
        elif scipy_call:
            selected = ['k', 'x', 'q', 'n', 'p', 'a', 'b', 'args', 'kwds', 'loc', 'scale']
        elif name == 'coerce_pext_value':
            if event == 'line':
                return self
            selected = ['value', 'item', 'values']
        elif name == '_variant_mechanism_masks':
            if event == 'line':
                return self
            selected = ['df']
        else:
            selected = ['anno_df', 'ped_df', 'fam_df', 'proband', 'sorted_df', 'drop_cols']
        state = {key: encode(values[key]) for key in selected if key in values}
        if event == 'return':
            state['return_value'] = encode(arg)
        if scipy_call:
            state['distribution'] = getattr(values['self'], 'name')
        key = self.store.put(state)
        frame_id = id(frame)
        if event == 'line' and self.previous.get(frame_id) == key:
            return self
        self.previous[frame_id] = key
        item = dict(event=event, function=name, file=str(path), line=frame.f_lineno,
                    phase=phase(name, frame.f_lineno), before_source_line=event == 'line', state=key)
        if scipy_call or event == 'exception':
            stack = []
            current = frame
            while current is not None:
                if Path(current.f_code.co_filename).parent == self.source:
                    stack.append(dict(function=current.f_code.co_name, line=current.f_lineno,
                                      phase=phase(current.f_code.co_name, current.f_lineno)))
                current = current.f_back
            item['source_stack'] = stack
        if event == 'exception':
            error = arg[1]
            item['exception'] = dict(type=typename(error), message=str(error))
            # Keep the actual object alive. A handled exception's integer id can
            # otherwise be reused by a later, unrelated uncaught exception.
            self.exception_events.append((error, len(self.events)))
        self.events.append(item)
        if event == 'return':
            self.previous.pop(frame_id, None)
        return self


def verify_sources(root, wanted):
    assert root.resolve().is_relative_to(D)
    assert digest(root / 'source-pins.json') == wanted
    pins = json.loads((root / 'source-pins.json').read_text())
    for name, expected in pins['files'].items():
        assert digest(root / 'source' / name) == expected, name
    return pins


def prepare(args):
    native = native_environment()
    assert args.repo.resolve().is_relative_to(Y)
    assert args.output.resolve().is_relative_to(D) and not args.output.exists()
    args.output.mkdir(parents=True)
    source = args.output / 'source'
    (source / 'scripts').mkdir(parents=True)
    pins = dict(schema='priva.full-ranking-oracle-source.v1', commit=COMMIT, repo=str(args.repo),
                created_utc=utc(), native_environment=native, files={}, functions={}, current_content_matches={},
                scientific_acceptance=False, rust_qualification=False, production_backend_selected=False)
    for name in SOURCE_FILES:
        raw = subprocess.run(['git', '-C', str(args.repo), 'show', COMMIT + ':scripts/' + name],
                             check=True, capture_output=True).stdout
        sha = hashlib.sha256(raw).hexdigest()
        if name in SCIENCE_PINS:
            assert sha == SCIENCE_PINS[name], name
        target = source / 'scripts' / name
        target.write_bytes(raw)
        pins['files']['scripts/' + name] = sha
        pins['current_content_matches'][name] = digest(args.repo / 'scripts' / name) == sha
        text = raw.decode()
        nodes = {node.name: node for node in ast.parse(text).body if isinstance(node, ast.FunctionDef)}
        for function in FUNCTIONS.get(name, []):
            node = nodes[function]
            pins['functions'][function] = dict(file=name, lines=[node.lineno, node.end_lineno],
                source_sha256=hashlib.sha256(ast.get_source_segment(text, node).encode()).hexdigest(),
                ast_sha256=hashlib.sha256(ast.dump(node).encode()).hexdigest())
    target = source / 'ranking_control_oracle.py'
    target.write_bytes(Path(__file__).read_bytes())
    pins['files'][target.name] = digest(target)
    write_json(args.output / 'source-pins.json', pins)
    wanted = digest(args.output / 'source-pins.json')
    verify_sources(args.output, wanted)
    print(json.dumps(dict(source_pins_sha256=wanted, source_files=len(pins['files']),
                          functions=pins['functions'], native_environment=native)))


def runtime_identity():
    """Installed build and actually loaded numerical-library identities, not versions alone."""
    import scipy.special._ufuncs as ufuncs
    import scipy.stats._continuous_distns as continuous
    import scipy.stats._discrete_distns as discrete
    import scipy.stats._distn_infrastructure as infrastructure
    import threadpoolctl
    expected = dict(python='3.12.3', numpy='1.26.4', pandas='2.2.3', scipy='1.14.1')
    actual = dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__)
    assert actual == expected, actual
    result = dict(versions=actual, python_version=sys.version, executable=str(Path(sys.executable).resolve()),
                  byteorder=sys.byteorder, machine=platform.machine(), libc=platform.libc_ver(),
                  numerical_threads={name: os.environ.get(name) for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']},
                  packages=[], native_files={}, wrapper_files={}, functions={}, ufuncs={}, threadpools=threadpoolctl.threadpool_info())
    prefix = Path(sys.prefix)
    for path in sorted((prefix / 'conda-meta').glob('*.json')):
        metadata = json.loads(path.read_text())
        result['packages'].append({**{key: metadata.get(key) for key in ['name', 'version', 'build', 'build_number', 'subdir', 'depends']},
                                   'metadata_path': str(path), 'metadata_sha256': digest(path)})
    for module in [continuous, discrete, infrastructure]:
        path = Path(module.__file__)
        result['wrapper_files'][str(path)] = digest(path)
    for name, function in [('binom_public_cdf', binom.cdf), ('binom_private_cdf', binom._cdf),
                           ('beta_public_ppf', beta.ppf), ('beta_private_ppf', beta._ppf)]:
        result['functions'][name] = dict(module=function.__module__, qualname=function.__qualname__,
            signature=str(inspect.signature(function)), source_file=inspect.getsourcefile(function),
            first_line=function.__code__.co_firstlineno,
            source_sha256=hashlib.sha256(inspect.getsource(function).encode()).hexdigest())
    for name in ['_binom_cdf', '_beta_ppf']:
        function = getattr(ufuncs, name)
        result['ufuncs'][name] = dict(type=typename(function), nin=function.nin, nout=function.nout,
                                     types=function.types, module_file=ufuncs.__file__)
    for name, module in list(sys.modules.items()):
        filename = getattr(module, '__file__', None)
        if name.startswith(('numpy.', 'scipy.', 'pandas.')) and filename and filename.endswith('.so'):
            path = Path(filename)
            result['native_files'][str(path)] = dict(sha256=digest(path), bytes=path.stat().st_size)
    for pool in result['threadpools']:
        path = Path(pool['filepath'])
        result['native_files'][str(path)] = dict(sha256=digest(path), bytes=path.stat().st_size)
    executable = Path(result['executable'])
    result['native_files'][str(executable)] = dict(sha256=digest(executable), bytes=executable.stat().st_size)
    return result


def base_frame(count=1, **changes):
    frame = pd.DataFrame(dict(
        __oracle_row_position=[f'input-row-{i}' for i in range(count)],
        ACMG_quant_score=[0.95] * count, BIOTYPE=['protein_coding'] * count,
        chrom=['1'] * count, pos=list(range(100, 100 + count)), ref=['A'] * count, alt=['G'] * count,
        variant_id=[f'v{i}' for i in range(count)], SYMBOL=['BOTH'] * count,
        var_plausible_patho_mechs=[''] * count, variant_effect=['predicted_LOF_high_confidence'] * count,
        PROBAND=['0/1'] * count, SECOND=['1/1'] * count,
    ))
    for name, value in changes.items():
        frame[name] = value
    frame.attrs = dict(fixture='synthetic ranking-entry table; no criteria assignment or family data')
    return frame


def ped(sex='1', **changes):
    frame = pd.DataFrame({'#FamilyID': ['F', 'F'], 'IndividualID': ['PROBAND', 'SECOND'],
                          'PaternalID': ['0', '0'], 'MaternalID': ['0', '0'],
                          'Sex': [sex, '2'], 'Phenotype': ['2', 2]})
    for name, value in changes.items():
        frame[name] = value
    return frame


def case_matrix(root):
    fixtures = root / 'fixtures'
    fixtures.mkdir()
    lists = {
        'dispensable': '# synthetic exact symbols\nBOTH\nDROP\nNA\n spaced \ncase\n',
        'relevant': '# synthetic exact symbols\nBOTH\nKEEP\nNA\n spaced \nCASE\n',
        'empty': '', 'comments': '# comment only\n# another comment\n',
    }
    paths = {}
    for name, content in lists.items():
        path = fixtures / (name + '.txt')
        path.write_text(content)
        paths[name] = str(path)
    paths['missing'] = str(fixtures / 'missing-list.txt')
    cases = []

    def add(name, kind, frame=None, pedigree=None, tags=(), expected=None, **kwargs):
        assert name not in {case['name'] for case in cases}
        cases.append(dict(name=name, kind=kind, frame=frame, ped=pedigree, tags=list(tags),
                          expected=expected, kwargs=kwargs))

    def rank(name, frame=None, pedigree=None, tags=(), expected=None, **kwargs):
        options = dict(fam_name='F' if pedigree is not None else None, dispensable_gene_list=None)
        options.update(kwargs)
        add(name, 'ranking', base_frame() if frame is None else frame, pedigree, tags, expected, **options)

    # Exhaust schema presence, including exact original partial-schema failures.
    for bits in range(16):
        columns = {name: [value] for bit, (name, value) in enumerate(zip(CONTROL_COLUMNS, [10, 100, .9, 4])) if bits & (1 << bit)}
        frame = pd.DataFrame(dict(PROBAND=['1/1'], **columns), index=pd.Index([7], name='control-label'))
        expected = dict(status='value') if bits in (0, 15) else dict(status='error', phase='C1' if not bits & 2 else 'C2')
        add(f'control_schema_{bits:04b}', 'control', frame, tags=['16-control-schema-patterns'], expected=expected, proband='PROBAND')

    def control(name, rows, tags=(), expected=None, **kwargs):
        frame = pd.DataFrame(rows)
        for column in CONTROL_COLUMNS:
            if column not in frame:
                frame[column] = np.nan
        add(name, 'control', frame, tags=tags, expected=expected, proband='PROBAND', **kwargs)

    control('control_route_priority', [
        dict(PROBAND='1/1', control_AN=100, control_nhomalt=4, control_AC=0, control_AF=.9),
        dict(PROBAND='1|1', control_AN=101, control_nhomalt=5, control_AC=1, control_AF=.9),
        dict(PROBAND='1', control_AN=99, control_nhomalt=6, control_AC=2, control_AF=.9),
        dict(PROBAND='0/1', control_AN=100, control_nhomalt=7, control_AC=3, control_AF=.9),
        dict(PROBAND='1/1:DP:1', control_AN=100, control_nhomalt=8, control_AC=4, control_AF=.9),
        dict(PROBAND='1/1', control_AN=100, control_nhomalt=np.nan, control_AC=5, control_AF=.9),
        dict(PROBAND='1/1', control_AN=100, control_nhomalt=np.nan, control_AC=np.nan, control_AF=.9),
        dict(PROBAND='0/1', control_AN=0, control_nhomalt=100, control_AC=100, control_AF=.9),
        dict(PROBAND=' 1/1', control_AN=100, control_nhomalt=10, control_AC=0, control_AF=.9),
        dict(PROBAND='1/1 ', control_AN=100, control_nhomalt=10, control_AC=0, control_AF=.9),
    ], tags=['C2-priority', 'GT-only-control-parser', 'actual-SciPy', 'original-mixed-AF-mask-failure'],
            expected=dict(status='error', phase='C4'))
    control('control_fractional_rounding_clipping', [
        dict(PROBAND='1/1', control_AN=an, control_nhomalt=k, control_AC=0, control_AF=.9)
        for an, k in [(3, .49), (5, .5), (7, 1.99), (9, -2.4), (.49, 0)]
    ] + [dict(PROBAND='0/1', control_AN=an, control_AC=k, control_AF=.9)
         for an, k in [(2.5, 1.9), (3.5, -2.9), (.49, .99), (100, 101)]], tags=['fractional', 'round-half-even', 'negative-counts'])
    control('control_missing_bad_text', [
        dict(PROBAND='0/1', control_AN=an, control_AC=ac, control_AF=af, control_nhomalt=hom)
        for an, ac, af, hom in [('bad', 10, .9, 1), (None, 10, .9, 1), (-1, 10, .9, 1),
                               (100, 'bad', '.9', 'bad'), ('100', '10', 'bad', None), (100, None, None, None)]
    ], tags=['missing', 'numeric-coercion'])
    for column in CONTROL_COLUMNS:
        for value, label in [(np.inf, 'positive_inf'), (-np.inf, 'negative_inf'), (np.nan, 'nan')]:
            row = dict(PROBAND='1/1', control_AN=100, control_AC=4, control_AF=.9, control_nhomalt=3)
            row[column] = value
            control(f'control_{column}_{label}', [row], tags=['nonfinite', 'original-failure-or-value'])
    for column in ['control_AC', 'control_AF']:
        for value, label in [(np.inf, 'positive_inf'), (-np.inf, 'negative_inf'), (np.nan, 'nan')]:
            row = dict(PROBAND='0/1', control_AN=100, control_AC=np.nan, control_AF=np.nan, control_nhomalt=np.nan)
            row[column] = value
            control(f'control_routed_{column}_{label}', [row], tags=['nonfinite', 'explicit-count-or-AF-route'])
    control('control_beta_original_ac_an', [
        dict(PROBAND='0/1', control_AN=100, control_AC=np.nan, control_AF=.9),
        dict(PROBAND='1/1', control_AN=100, control_AC=0, control_nhomalt=30, control_AF=0),
        dict(PROBAND='0/1', control_AN=2, control_AC=4, control_AF=.9),
        dict(PROBAND='0/1', control_AN=100.5, control_AC=95.7, control_AF=0),
        dict(PROBAND='0/1', control_AN=100, control_AC=-1.5, control_AF=.9),
    ], tags=['Beta-uses-original-counts', 'AC-greater-than-AN', 'original-mixed-AF-mask-failure'],
            expected=dict(status='error', phase='C4'))
    for name, row in [
        ('AF_missing_AC', dict(PROBAND='0/1', control_AN=100, control_AC=np.nan, control_AF=.9)),
        ('hom_AC_zero', dict(PROBAND='1/1', control_AN=100, control_AC=0, control_nhomalt=30, control_AF=0)),
        ('invalid_shape', dict(PROBAND='0/1', control_AN=2, control_AC=4, control_AF=.9)),
        ('fractional_original_counts', dict(PROBAND='0/1', control_AN=100.5, control_AC=95.7, control_AF=0)),
    ]:
        control('control_beta_single_' + name, [row], tags=['actual-Beta-PPF', 'original-counts', name],
                expected=dict(status='value'))
    control('control_large_counts', [dict(PROBAND='0/1', control_AN=10**12, control_AC=2*10**10, control_AF=.02),
                                     dict(PROBAND='1/1', control_AN=10**12+1, control_nhomalt=10**10, control_AC=2*10**10)], tags=['large-counts'])
    strict_cdf = float(binom.cdf(4, 100, .01))
    for value, label in [(np.nextafter(strict_cdf, -np.inf), 'below'), (strict_cdf, 'equal'), (np.nextafter(strict_cdf, np.inf), 'above')]:
        control('control_cdf_confidence_' + label, [dict(PROBAND='0/1', control_AN=100, control_AC=4)],
                tags=['strict-CDF-neighbor', label], expected=dict(status='value'), conf_level=value)
    control('control_af_cutoff_neighbors', [dict(PROBAND='0/1', control_AN=100, control_AF=value)
            for value in [np.nextafter(.01, -np.inf), .01, np.nextafter(.01, np.inf)]],
            tags=['strict-AF-neighbor'], expected=dict(status='value'))
    for cutoff in [0.0, 1.0, -0.01, np.nan, np.inf]:
        control('control_cutoff_' + repr(cutoff), [dict(PROBAND='0/1', control_AN=100, control_AC=4, control_AF=.9)],
                tags=['parameter-domain'], cutoff=cutoff)
    for confidence in [-1.0, 0.0, 1.0, np.nan, np.inf]:
        control('control_confidence_' + repr(confidence), [dict(PROBAND='0/1', control_AN=100, control_AC=4, control_AF=.9)],
                tags=['parameter-domain'], conf_level=confidence)
    control('control_Beta_penalty_floor_and_ceiling', [
        dict(PROBAND='0/1', control_AN=100000, control_AC=99999, control_AF=0),
        dict(PROBAND='0/1', control_AN=100, control_AC=10, control_AF=0),
    ], tags=['Beta-penalty-clipping'])
    for proband in [None, 'ABSENT']:
        frame = pd.DataFrame(dict(control_AN=[100], control_AC=[0], control_AF=[.9], control_nhomalt=[10]))
        add('control_proband_' + str(proband), 'control', frame, tags=['missing-proband-AC-route'], proband=proband)
    for dtype_name in ['Int64', 'Float64', 'string']:
        frame = pd.DataFrame({name: pd.Series([None, value], dtype=dtype_name) for name, value in zip(CONTROL_COLUMNS, [10, 100, 1, 4])})
        frame['PROBAND'] = ['0/1', '1/1']
        add('control_nullable_' + dtype_name, 'control', frame, tags=['extension-dtypes'], proband='PROBAND')
    duplicate = pd.DataFrame(dict(control_AN=[100, 100, 100], control_AC=[0, 10, 0], control_AF=[.9, .9, .9],
                                  control_nhomalt=[4, np.nan, np.nan], PROBAND=['1/1', '0/1', '0/1']), index=[7, 7, -2])
    add('control_duplicate_index', 'control', duplicate, tags=['duplicate-index', 'positional-vs-label'], proband='PROBAND')

    frame = base_frame(5, BIOTYPE=['protein_coding', 'lncRNA', 'protein_coding', 'protein_coding', 'other'])
    frame.index = pd.Index([7, 100, 7, -3, 41], name='duplicate-gapped-index')
    rank('rank_filter_duplicate_gapped', frame, tags=['R1', 'duplicate-index', 'gapped-index'])
    rank('rank_empty', base_frame(0), tags=['empty-frame'])
    rank('rank_all_filtered', base_frame(2, BIOTYPE=['lncRNA', '']), tags=['R1', 'all-filtered'])
    for sex, label in [('1', 'male1'), ('male', 'male_word'), ('1.0', 'male_float_text'), ('2', 'female2'),
                       ('female', 'female_word'), ('2.0', 'female_float_text'), ('0', 'unknown'), (None, 'missing'), (' M ', 'strip_case')]:
        profile = ['', 'dominant_LOF', 'dominant', 'dominant_DN', 'dominant_GOF', 'recessive_LOF',
                   'x_linked_recessive_LOF', 'x_linked_dominant_LOF', 'x_linked_recessive',
                   'x_linked_dominant', ' dominant_LOF', 'dominant_LOF\n']
        frame = base_frame(len(profile), chrom=['X'] * len(profile), var_plausible_patho_mechs=profile)
        rank('rank_sex_' + label, frame, ped(sex), tags=['R2', 'R4', 'X-sex-compatibility'])
    rank('rank_no_sex_column', base_frame(), ped().drop(columns=['Sex']), tags=['R2'])
    sample_strings = ['0/1', '0/1:10:20', '0/1:8,1:9', '1/1', '1|1', '1:20', '0|1:20', '0/0:1', './.:1', None]
    rank('rank_entire_sample_string_count', base_frame(len(sample_strings), PROBAND=sample_strings), ped(), tags=['R2', 'full-sample-string'])
    rank('rank_first_affected_row', base_frame(), ped(IndividualID=['SECOND', 'PROBAND']), tags=['PED-first-affected'])
    rank('rank_numeric_sample_str_failure', base_frame(PROBAND=[1]), ped(), tags=['original-exception'],
         expected=dict(status='error', phase='R2'))
    rank('rank_missing_affected', base_frame(), ped(Phenotype=['1', '1']), tags=['PED-exception'], expected=dict(status='error', phase='R2'))
    rank('rank_missing_sample', base_frame().drop(columns=['PROBAND']), ped(), tags=['PED-exception'], expected=dict(status='error', phase='R2'))
    rank('rank_family_type_mismatch', base_frame(), ped(), fam_name=1, tags=['PED-type-semantics'], expected=dict(status='error', phase='R2'))
    for missing_column in ['ACMG_quant_score', 'BIOTYPE', 'var_plausible_patho_mechs', 'variant_effect', 'chrom', 'pos']:
        rank('rank_missing_' + missing_column, base_frame().drop(columns=[missing_column]), tags=['required-column-failure-order'])
    frame = base_frame().drop(columns=['var_plausible_patho_mechs', 'variant_effect'])
    frame['control_AF'] = .9
    rank('rank_control_failure_before_mechanism', frame, tags=['first-failure-order'], expected=dict(status='error', phase='C1'))
    effects = ['predicted_LOF_high_confidence', 'uncertain', 'exact_known_GOF', 'exact_known_GOF+DOMINANT_NEGATIVE',
               'exact_known_DOMINANT_NEGATIVE', 'LOF', 'predicted_lof_high_confidence', 'predicted_LOF_high_confidence\n']
    rank('rank_exact_effect_and_scores', base_frame(len(effects), variant_effect=effects, variant_lof_score=[1]*len(effects),
                                                 variant_gof_score=[0, 0, 0, 2, 0, 'bad', np.nan, 0]), ped(), tags=['R4', 'exact-tags'])
    rank('rank_compatibility_raw_threshold', base_frame(3, ACMG_quant_score=[np.nextafter(.89, -np.inf), .89, np.nextafter(.89, np.inf)]),
         ped(), tags=['strict-raw-score-neighbor'])
    symbols = ['BOTH', 'DROP', 'KEEP', 'NA', None, ' spaced ', 'spaced', 'case', 'CASE']
    frame = base_frame(len(symbols), SYMBOL=symbols, ACMG_quant_score=[.9500000000000001]*len(symbols))
    rank('rank_gene_overlap_sequential', frame, ped(), dispensable_gene_list=paths['dispensable'],
         relevant_gene_list=paths['relevant'], tags=['R5', 'R6', 'sequential-floating-multiplication', 'NA-gene-symbol'])
    frame = base_frame(3, control_AN=[100, 2, 100], control_AC=[95, 4, 10], control_AF=[0, .9, .9],
                       control_nhomalt=[np.nan]*3, PEXT_MEAN=[.02, .02, .02])
    rank('rank_control_gene_compatibility_pext_sequence', frame, ped(),
         dispensable_gene_list=paths['dispensable'], relevant_gene_list=paths['relevant'],
         tags=['R3', 'R5', 'R6', 'R7', 'all-multiplications-in-source-order', 'Beta-NaN'])
    for option in ['dispensable_gene_list', 'relevant_gene_list']:
        for name in ['empty', 'comments', 'missing']:
            rank('rank_' + option + '_' + name, **{option: paths[name]}, tags=['gene-list-file-contract'])
    frame = base_frame()
    frame['dispensable_gene_index'] = 123.0; frame['relevant_gene_index'] = 456.0; frame['pext_sort_index'] = 789.0
    rank('rank_skipped_branches_retain_stale', frame, tags=['existing-column-position', 'skipped-column-retention'])
    original = base_frame(3, PEXT_MEAN=[.02, .03, .04])
    original.insert(1, 'variant_rank', [-10]*3); original.insert(2, 'sort_index', [-20.0]*3)
    original.insert(3, 'pext_sort_index', [-30.0]*3); original.insert(4, 'control_common_index', [-40.0]*3)
    original.insert(5, 'haplo_insuf_index', [-50.0]*3); original.insert(6, 'zygosity_inheritance_mechanism_index', [-60.0]*3)
    original.insert(7, 'max_variant_score', [-70.0]*3)
    rank('rank_existing_columns_positions', original, tags=['R1', 'R8', 'existing-column-position'])
    pext_values = ['NAN&0.02', '0.02&NAN', 'inf&0.2', '-inf', None, np.nan, 'NaN', '-nan', .02, 'bad&.02', '0&-0']
    rank('rank_pext_nan_inf_and_tokens', base_frame(len(pext_values), PEXT_MEAN=pext_values), tags=['R7', 'PEXT-order', 'PEXT-nonfinite'])
    rank('rank_pext_cutoff_neighbors', base_frame(3, PEXT_MEAN=[np.nextafter(.1, -np.inf), .1, np.nextafter(.1, np.inf)]),
         tags=['R7', 'strict-PEXT-cutoff-neighbor'])
    rank('rank_handled_pext_error_then_missing_coordinate', base_frame(PEXT_MEAN=['bad&.02']).drop(columns=['pos']),
         tags=['handled-parser-error-before-terminal-exception'], expected=dict(status='error', phase='R8'))
    frame = base_frame(3, PEXT_MEAN=[.02]*3, PEXT_Brain_Cortex=[.04]*3, PEXT_Liver=[.08]*3)
    for tissues, label in [('', 'mean'), ('Brain-Cortex,Liver', 'first_present'), ('Missing,Liver', 'first_missing_fallback_mean'),
                           (', ,Brain-Cortex,Liver', 'empty_tokens'), ('Missing', 'missing_fallback_mean')]:
        rank('rank_tissue_' + label, frame.copy(), pext_tissues=tissues, tags=['R7', 'first-tissue-fallback'])
    rank('rank_tissue_second_present_without_mean', frame.drop(columns=['PEXT_MEAN', 'PEXT_Brain_Cortex']),
         pext_tissues='Missing,Liver', pext_low_expression_cutoff=0, tags=['R7', 'does-not-use-second-tissue'])
    for name, value in [('pext_low_expression_cutoff', 0), ('pext_low_expression_cutoff', -1),
                        ('pext_penalty_floor', -.01), ('pext_penalty_floor', 1.01), ('pext_penalty_shape', 0),
                        ('pext_low_expression_cutoff', np.nan), ('pext_penalty_floor', np.nan), ('pext_penalty_shape', np.nan),
                        ('pext_low_expression_cutoff', np.inf), ('pext_penalty_shape', np.inf)]:
        for selected in [False, True]:
            frame = base_frame(PEXT_MEAN=[.02]) if selected else base_frame()
            rank(f'rank_parameter_{name}_{value}_{selected}', frame, **{name: value}, tags=['R7', 'parameter-validation-gated'])
    rank('rank_pext_raw_threshold_not_lof', base_frame(3, PEXT_MEAN=[.02]*3, variant_effect=['uncertain']*3,
         ACMG_quant_score=[np.nextafter(.89, -np.inf), .89, np.nextafter(.89, np.inf)]), tags=['R7', 'strict-raw-score-neighbor', 'not-LOF'])
    grouped = base_frame(6, chrom=['2', '10', '2', '2', None, '1'], pos=[100, 100, 100, 101, 100, None],
                         ACMG_quant_score=[.8, .95, .9, np.nan, .99, .99],
                         variant_id=['shared', 'shared', None, np.nan, pd.NA, 'last'])
    rank('rank_group_nulls_missing_shared_ids', grouped, tags=['R8', 'null-coordinate-group', 'missing-ID-rank-zero', 'shared-IDs'])
    rank('rank_absent_ids_null_coordinate_spellings', grouped.drop(columns=['variant_id']), tags=['R8', 'absent-ID-stringification'])
    rank('rank_all_nan_group', base_frame(3, pos=[100, 100, 101], ACMG_quant_score=[np.nan, np.nan, .9]), tags=['R8', 'all-null-group'])
    rank('rank_nonfinite_scores', base_frame(5, ACMG_quant_score=[np.inf, -np.inf, np.nan, .95, -0.0]), tags=['R8', 'nonfinite-score'])
    rank('rank_ties_duplicate_labels', base_frame(4, chrom=['1']*4, pos=[100]*4, variant_id=['b', 'a', 'b', None]), tags=['R8', 'equal-six-keys'])
    for field in ['chrom', 'pos', 'ref', 'alt']:
        frame = base_frame(3); frame.loc[1, field] = None
        rank('rank_null_key_' + field, frame, tags=['R8', 'null-coordinate-group'])
    for values, label in [([2, 10, 1], 'integer'), (['2', '10', '1'], 'string'), ([2, '10', 1], 'mixed')]:
        rank('rank_chrom_dtype_' + label, base_frame(3, chrom=values), tags=['R8', 'typed-sort'])
    for values, label in [([0, -1, 1.5], 'zero_negative_fractional'), (['100', '2', '11'], 'string')]:
        rank('rank_pos_dtype_' + label, base_frame(3, pos=values), tags=['R8', 'typed-sort'])
    for dtype_name in ['Float64', 'Int64']:
        frame = base_frame(3)
        frame['ACMG_quant_score'] = pd.Series([1, None, 0], dtype=dtype_name)
        frame['pos'] = pd.Series([100, None, 102], dtype='Int64')
        rank('rank_nullable_' + dtype_name, frame, tags=['R8', 'extension-dtypes'])
    for field in ['ACMG_quant_score', 'variant_id', 'chrom', 'PEXT_MEAN']:
        frame = base_frame(PEXT_MEAN=[.02])
        frame = pd.concat([frame, frame[[field]]], axis=1)
        rank('rank_duplicate_column_' + field, frame, tags=['duplicate-column', 'original-exception'])
    # Caller parsing is source-executed, separately identified from the library.
    annotation = fixtures / 'annotation.tsv'
    parser_frame = base_frame(3, SYMBOL=['NA', 'BOTH', 'BOTH'], variant_id=['NA', 'v', 'v'])
    parser_frame.iloc[2] = parser_frame.iloc[1]
    parser_bytes = parser_frame.to_csv(sep='\t', index=False).encode()
    annotation.write_bytes(parser_bytes)
    ped_bytes = {
        'headered': '#FamilyID\tIndividualID\tPaternalID\tMaternalID\tSex\tPhenotype\nF\tPROBAND\t0\t0\t1\t2\n',
        'headerless': 'F\tPROBAND\t0\t0\t1\t2\n',
        'leading_zero_family': '#FamilyID\tIndividualID\tPaternalID\tMaternalID\tSex\tPhenotype\n001\tPROBAND\t0\t0\t1\t2\n',
        'NA_sample': '#FamilyID\tIndividualID\tPaternalID\tMaternalID\tSex\tPhenotype\nF\tNA\t0\t0\t0\t2\n',
    }
    for name, content in ped_bytes.items():
        path = fixtures / ('ped_' + name + '.tsv'); path.write_text(content)
        add('production_parser_' + name, 'production_parser', tags=['caller-parser-only', 'production-string-PED'],
            anno_table=str(annotation), ped_table=str(path), fam_name='001' if name == 'leading_zero_family' else 'F')
        target = fixtures / ('resort_' + name + '.tsv')
        add('resort_parser_' + name, 'resort', tags=['actual-resort-function', 'inferred-PED'],
            filtered_tsv=str(target), original_bytes=parser_bytes, ped_table=str(path),
            fam_name='001' if name == 'leading_zero_family' else 'F', pext_tissues='', relevant_gene_list=None,
            dispensable_gene_list=None, pext_low_expression_cutoff=.1, pext_penalty_floor=.8, pext_penalty_shape=.5)
    for label, raw in [('zero_bytes', b''), ('whitespace', b' \t\n'), ('header_only', b'ACMG_quant_score\tBIOTYPE\n'),
                       ('missing_score', b'BIOTYPE\nprotein_coding\n'), ('missing_file', None)]:
        target = fixtures / ('resort_' + label + '.tsv')
        add('resort_input_' + label, 'resort', tags=['actual-resort-input-contract'], filtered_tsv=str(target),
            original_bytes=raw, ped_table=None, fam_name=None, pext_tissues='', relevant_gene_list=None,
            dispensable_gene_list=None, pext_low_expression_cutoff=.1, pext_penalty_floor=.8, pext_penalty_shape=.5)
    return cases, paths


def production_parser_code(source):
    """Execute verbatim original parser statements, excluding scientific assignment."""
    path = source / 'acmg_criteria_assign.py'
    function = next(node for node in ast.parse(path.read_text()).body
                    if isinstance(node, ast.FunctionDef) and node.name == 'ACMG_criteria_assign')
    nodes = []
    for node in function.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if (name == 'anno_df' and node.lineno == 196) or (
                    name in {'ped_df', 'fam_df', 'proband'} and isinstance(node.value, ast.Constant) and node.value.value is None):
                nodes.append(node)
        if isinstance(node, ast.If) and isinstance(node.test, ast.Name) and node.test.id == 'ped_table':
            nodes.append(node)
    assert len(nodes) == 5, [node.lineno for node in nodes]
    details = [dict(lines=[node.lineno, node.end_lineno], source_sha256=hashlib.sha256(
        ast.get_source_segment(path.read_text(), node).encode()).hexdigest()) for node in nodes]
    return compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), details


def original_modules(source):
    for name in ['acmg_scoring', 'acmg_consequence', 'acmg_variant_mechanism', 'resort_step3_only']:
        assert name not in sys.modules, 'fresh process required for frozen imports'
    sys.path.insert(0, str(source))
    scoring = importlib.import_module('acmg_scoring')
    resort = importlib.import_module('resort_step3_only')
    for module in [scoring, resort, sys.modules['acmg_consequence'], sys.modules['acmg_variant_mechanism']]:
        assert Path(module.__file__).parent == source
    parser_code, parser_details = production_parser_code(source)
    return scoring, resort, parser_code, parser_details


def exception_record(error, source):
    records = []
    for item in traceback.extract_tb(error.__traceback__):
        records.append(dict(file=item.filename, line=item.lineno, function=item.name, source_line=item.line,
                            frozen_scientific_source=Path(item.filename).parent == source))
    return dict(type=typename(error), message=str(error), args=encode(error.args), traceback=records)


def execute_case(case, scoring, resort, parser_code, source, store, folder, traced, prepared_input):
    folder.mkdir()
    # pandas deep copies share mutable dtype caches. Independent loads from the
    # SAME prepared bytes isolate those caches without normalizing source data.
    frame, pedigree = pickle.loads(prepared_input)
    kwargs = dict(case['kwargs'])
    file_path, original_bytes = None, None
    if case['kind'] == 'resort':
        file_path = Path(kwargs['filtered_tsv'])
        original_bytes = kwargs.pop('original_bytes')
        if original_bytes is None:
            assert not file_path.exists()
        else:
            file_path.write_bytes(original_bytes)
    before_frame, before_ped = encode(frame), encode(pedigree)
    write_json(folder / 'input.json', dict(frame=before_frame, pedigree=before_ped))
    if case['kind'] == 'control':
        function, positional = scoring.compute_control_common, [frame]
    elif case['kind'] == 'ranking':
        function, positional = scoring.sort_and_rank_variants, [frame, pedigree]
    elif case['kind'] == 'resort':
        function, positional = resort.resort_family, []
    else:
        function, positional = None, []
    if function is not None:
        bound = inspect.signature(function).bind(*positional, **kwargs)
        bound.apply_defaults()
        scientific_arguments = dict(function=function.__name__, signature=str(inspect.signature(function)),
                                    arguments=encode(dict(bound.arguments)))
    else:
        scientific_arguments = dict(function='original-production-parser-excerpts-only', arguments=encode(kwargs))
    write_json(folder / 'arguments.json', scientific_arguments)
    before_pickle = pickle.dumps((frame, pedigree), protocol=5)
    (folder / 'input.pkl').write_bytes(before_pickle)
    observer = Trace(source, store) if traced else None
    previous = sys.gettrace()
    log = io.StringIO()
    handler = logging.StreamHandler(log)
    handler.setFormatter(logging.Formatter('%(levelname)s %(name)s %(message)s'))
    logging.getLogger().addHandler(handler)
    result, caught = None, None
    try:
        with warnings.catch_warnings(record=True) as observed_warnings:
            if observer:
                sys.settrace(observer)
            try:
                if function is not None:
                    result = function(*positional, **kwargs)
                else:
                    namespace = dict(pd=pd, **kwargs)
                    exec(parser_code, namespace)
                    result = {name: namespace[name] for name in ['anno_df', 'ped_df', 'fam_df', 'proband']}
            except Exception as error:
                caught = error
            finally:
                sys.settrace(previous)
    finally:
        logging.getLogger().removeHandler(handler)
    (folder / 'source.log').write_text(log.getvalue())
    # Capture serialization immediately after execution, before post-observation
    # can populate another lazy pandas metadata property.
    after_pickle = pickle.dumps((frame, pedigree), protocol=5)
    after_frame, after_ped = encode(frame), encode(pedigree)
    (folder / 'input-after.pkl').write_bytes(after_pickle)
    write_json(folder / 'input-after.json', dict(frame=after_frame, pedigree=after_ped))
    record = dict(status='error' if caught else 'value', output=encode(result),
                  returned_input_frame=result is frame and frame is not None,
                  input_frame_semantically_unchanged=before_frame == after_frame,
                  pedigree_semantically_unchanged=before_ped == after_ped,
                  input_pickle_bytes_unchanged=before_pickle == after_pickle,
                  prepared_input_sha256=hashlib.sha256(prepared_input).hexdigest(),
                  input_pickle_before_sha256=hashlib.sha256(before_pickle).hexdigest(),
                  input_pickle_after_sha256=hashlib.sha256(after_pickle).hexdigest(),
                  warnings=[dict(type=typename(item.message), message=str(item.message), file=item.filename, line=item.lineno)
                            for item in observed_warnings])
    if caught:
        record['exception'] = exception_record(caught, source)
    if file_path is not None:
        after_bytes = file_path.read_bytes() if file_path.exists() else None
        record['file_effect'] = dict(existed_before=original_bytes is not None, exists_after=after_bytes is not None,
                                    bytes_unchanged=after_bytes == original_bytes,
                                    before_sha256=hashlib.sha256(original_bytes).hexdigest() if original_bytes is not None else None,
                                    after_sha256=hashlib.sha256(after_bytes).hexdigest() if after_bytes is not None else None)
        if original_bytes is not None:
            (folder / 'file-before.tsv').write_bytes(original_bytes)
        if after_bytes is not None:
            (folder / 'file-after.tsv').write_bytes(after_bytes)
    if result is not None:
        (folder / 'output.pkl').write_bytes(pickle.dumps(result, protocol=5))
    write_json(folder / 'outcome.json', record)
    failure = None
    if observer:
        if caught:
            matches = [index for error, index in observer.exception_events if error is caught]
            assert matches, 'uncaught exception must have an observed source propagation'
            event = observer.events[matches[0]]
            assert event['exception'] == dict(type=typename(caught), message=str(caught))
            core = [item for item in event.get('source_stack', []) if item['phase'][:1] in {'C', 'R'} and item['phase'][1:2].isdigit()]
            failure = dict(event_index=matches[0], leaf_function=event['function'], leaf_line=event['line'],
                           phase=core[0]['phase'] if core else event['phase'], source_stack=event.get('source_stack', []))
        write_json(folder / 'trace.json', dict(schema='priva.original-source-trace.v1', events=observer.events,
            first_uncaught_failure=failure, semantics='line states are observed before that original source line; '
            'preceding statements have executed; return states contain actual returned values; no source or SciPy replacement'))
    return record, failure


def self_tests():
    """Check that the observer can expose the distinctions this oracle must retain."""
    checks = []
    def check(name, condition):
        assert condition, name
        checks.append(name)
    check('signed-zero-bits', encode(0.0) != encode(-0.0))
    check('nan-payload-bits', encode(float('nan')) != encode(struct.unpack('>d', bytes.fromhex('7ff8000000000123'))[0]))
    check('numpy-vs-python-bool', encode(True) != encode(np.bool_(True)))
    check('none-pdNA-nan-distinct', len({canonical(encode(x)) for x in [None, pd.NA, float('nan')]}) == 3)
    check('dtype-distinct', encode(pd.Series([1], dtype='Int64')) != encode(pd.Series([1.0], dtype='float64')))
    check('null-mask-distinct', encode(pd.Series([pd.NA], dtype='Float64')) != encode(pd.Series([np.nan], dtype='float64')))
    check('column-order-distinct', encode(pd.DataFrame({'a': [1], 'b': [2]})) != encode(pd.DataFrame({'b': [2], 'a': [1]})))
    check('duplicate-label-order', encode(pd.Index([7, 7, -2], name='i')) != encode(pd.Index([7, -2, 7], name='i')))
    check('index-name-distinct', encode(pd.Index([1], name='a')) != encode(pd.Index([1], name='b')))
    check('input-value-mutation-detected', encode(pd.DataFrame({'a': [1.0]})) != encode(pd.DataFrame({'a': [np.nextafter(1., 2.)]})))
    nullable = pd.Series([1.0, None], dtype='Float64')
    _ = nullable.dtype.numpy_dtype
    prepared = pickle.dumps(nullable, protocol=5)
    first, second = pickle.loads(prepared), pickle.loads(prepared)
    check('independent-dtype-cache-objects', first.dtype is not second.dtype and first.dtype._cache is not second.dtype._cache)
    _ = first.dtype.index_class
    check('warming-one-dtype-does-not-warm-peer', 'index_class' not in second.dtype._cache)
    return checks


def capture(args):
    environment = native_environment()
    pins = verify_sources(args.root, args.pins_sha256)
    assert digest(Path(__file__)) == pins['files']['ranking_control_oracle.py']
    assert not (args.root / 'capture.json').exists()
    source = args.root / 'source/scripts'
    scoring, resort, parser_code, parser_details = original_modules(source)
    identities = runtime_identity()
    write_json(args.root / 'runtime-identities.json', identities)
    functions_before = {name: hashlib.sha256(function.__code__.co_code).hexdigest() for name, function in [
        ('compute_control_common', scoring.compute_control_common), ('sort_and_rank_variants', scoring.sort_and_rank_variants),
        ('coerce_pext_value', scoring.coerce_pext_value), ('_variant_mechanism_masks', scoring._variant_mechanism_masks)]}
    cases, lists = case_matrix(args.root)
    fixed_inputs = {str(path.relative_to(args.root)): digest(path) for path in (args.root / 'fixtures').iterdir() if path.is_file()}
    store = StateStore(args.root / 'states')
    (args.root / 'cases').mkdir()
    checks = self_tests()
    reports = []
    for position, case in enumerate(cases):
        folder = args.root / 'cases' / f"{position:03d}_{case['name']}"
        folder.mkdir()
        prepared_input = pickle.dumps((case['frame'], case['ped']), protocol=5)
        (folder / 'prepared-input.pkl').write_bytes(prepared_input)
        untraced, _ = execute_case(case, scoring, resort, parser_code, source, store, folder / 'untraced', False, prepared_input)
        traced, failure = execute_case(case, scoring, resort, parser_code, source, store, folder / 'traced', True, prepared_input)
        # Warnings may be emitted once by pandas/SciPy; scientific results,
        # exceptions/traceback, inputs and source file effects must still match.
        scientific_untraced = {key: value for key, value in untraced.items() if key != 'warnings'}
        scientific_traced = {key: value for key, value in traced.items() if key != 'warnings'}
        equivalent = scientific_untraced == scientific_traced
        report = dict(name=case['name'], kind=case['kind'], tags=case['tags'], expected=case['expected'],
                      folder=str(folder.relative_to(args.root)), status=traced['status'], observer_exactly_equivalent=equivalent,
                      input_frame_semantically_unchanged=traced['input_frame_semantically_unchanged'],
                      returned_input_frame=traced['returned_input_frame'], first_uncaught_failure=failure)
        if traced['status'] == 'error':
            report['exception'] = {key: traced['exception'][key] for key in ['type', 'message']}
        reports.append(report)
        write_json(folder / 'case.json', report)
        if (position + 1) % 25 == 0:
            print(json.dumps(dict(captured=position + 1, total=len(cases), equivalent=all(item['observer_exactly_equivalent'] for item in reports))), flush=True)
    verify_sources(args.root, args.pins_sha256)
    for name, wanted in fixed_inputs.items():
        assert digest(args.root / name) == wanted, 'source input fixture mutated: ' + name
    functions_after = {name: hashlib.sha256(function.__code__.co_code).hexdigest() for name, function in [
        ('compute_control_common', scoring.compute_control_common), ('sort_and_rank_variants', scoring.sort_and_rank_variants),
        ('coerce_pext_value', scoring.coerce_pext_value), ('_variant_mechanism_masks', scoring._variant_mechanism_masks)]}
    assert functions_before == functions_after
    report = dict(schema='priva.full-ranking-control-python-oracle.v1', started_source_commit=COMMIT,
                  finished_utc=utc(), native_environment=environment, source_pins_sha256=args.pins_sha256,
                  runtime_identities_sha256=digest(args.root / 'runtime-identities.json'),
                  function_bytecode_unchanged=functions_after, production_parser_original_statement_pins=parser_details,
                  self_tests=checks, fixed_input_files=fixed_inputs, gene_list_paths=lists, cases=reports,
                  counts=dict(cases=len(reports), values=sum(item['status'] == 'value' for item in reports),
                              original_exceptions=sum(item['status'] == 'error' for item in reports),
                              observer_mismatches=sum(not item['observer_exactly_equivalent'] for item in reports)),
                  scope='synthetic unchanged-Python ranking/control and distinct caller-parser groundwork; '
                        'no criterion assignment, Rust integration, family or installation acceptance',
                  scientific_acceptance=False, rust_qualification=False, production_backend_selected=False)
    report['files'] = {str(path.relative_to(args.root)): digest(path) for path in sorted(args.root.rglob('*')) if path.is_file()}
    write_json(args.root / 'capture.json', report)
    print(json.dumps(dict(capture_sha256=digest(args.root / 'capture.json'), **report['counts'])), flush=True)
    assert not report['counts']['observer_mismatches'], 'preserved observer mismatch requires review'


def decode_atom(value):
    if value['kind'] == 'numpy_scalar':
        return np.frombuffer(bytes.fromhex(value['storage_hex']), dtype=np.dtype(value['dtype']))[0]
    if value['kind'] == 'float':
        return struct.unpack('>d', bytes.fromhex(value['ieee_f64_be']))[0]
    if value['kind'] == 'int':
        return int(value['decimal'])
    if value['kind'] == 'str':
        return value['value']
    if value['kind'] == 'bool':
        return value['value']
    if value['kind'] == 'none':
        return None
    raise AssertionError(value)


def validate_core(args, native):
    assert args.root.resolve().is_relative_to(D)
    assert digest(args.root / 'capture.json') == args.capture_sha256
    captured = json.loads((args.root / 'capture.json').read_text())
    verify_sources(args.root, captured['source_pins_sha256'])
    for name, wanted in captured['files'].items():
        assert digest(args.root / name) == wanted, name
    reports = {item['name']: item for item in captured['cases']}
    assert len(reports) == len(captured['cases'])
    checks = list(captured['self_tests'])
    failures = []

    def check(name, condition):
        checks.append(name)
        if not condition:
            failures.append(name)

    def output(name):
        with (args.root / reports[name]['folder'] / 'traced/output.pkl').open('rb') as stream:
            return pickle.load(stream)

    def trace(name):
        return json.loads((args.root / reports[name]['folder'] / 'traced/trace.json').read_text())['events']

    def states(name, function=None, event=None):
        return [(item, json.loads((args.root / 'states' / (item['state'] + '.json')).read_text()))
                for item in trace(name) if (function is None or item['function'] == function)
                and (event is None or item['event'] == event)]

    check('all-traced-untraced-exact', all(item['observer_exactly_equivalent'] for item in reports.values()))
    for name, item in reports.items():
        wanted = item['expected']
        if wanted:
            check(name + '-status', item['status'] == wanted['status'])
            if 'phase' in wanted:
                check(name + '-first-failure-phase', item['first_uncaught_failure']['phase'] == wanted['phase'])
        if item['status'] == 'error':
            check(name + '-uncaught-failure-retained', item['first_uncaught_failure'] is not None)
        if item['kind'] == 'ranking':
            check(name + '-ranker-input-nonmutation', item['input_frame_semantically_unchanged'])
            if item['status'] == 'value':
                check(name + '-ranker-returns-copy', not item['returned_input_frame'])
        if item['kind'] == 'control' and item['status'] == 'value':
            check(name + '-control-returns-same-frame', item['returned_input_frame'])
    check('all-16-control-presence-patterns', len([name for name in reports if name.startswith('control_schema_')]) == 16)
    check('control-no-column-inplace-mutation', not reports['control_schema_0000']['input_frame_semantically_unchanged'])
    calls = [state for _, state in states('control_route_priority', '_cdf', 'call')]
    check('real-SciPy-private-CDF-both-count-routes', len(calls) == 2)
    # SciPy1.14.1 names the private binomial count argument x, not k.
    check('nhomalt-priority-exact-k', len(calls) >= 1 and [decode_atom(x) for x in calls[0]['x']['values']] == [4, 5, 6, 8])
    check('nhomalt-AN-half-ties-even-n', len(calls) >= 1 and [decode_atom(x) for x in calls[0]['n']['values']] == [50, 50, 50, 50])
    check('AC-route-excludes-nhomalt-rows', len(calls) >= 2 and [decode_atom(x) for x in calls[1]['x']['values']] == [3, 5, 0, 0])
    check('real-SciPy-Beta-public-PPF', bool(states('control_beta_single_fractional_original_counts', 'ppf', 'call')))
    check('real-SciPy-Beta-native-wrapper', bool(states('control_beta_single_fractional_original_counts', '_ppf', 'call')))
    check('Beta-invalid-shape-remains-NaN', np.isnan(output('control_beta_single_invalid_shape').loc[0, 'control_common_index']))
    check('AF-route-Beta-retains-missing-original-AC-as-zero', output('control_beta_single_AF_missing_AC').loc[0, 'control_common_index'] == .5)
    check('hom-route-Beta-retains-original-AC-zero', output('control_beta_single_hom_AC_zero').loc[0, 'control_common_index'] == .5)
    for name, wanted in [('below', True), ('equal', False), ('above', False)]:
        check('strict-CDF-' + name, bool(output('control_cdf_confidence_' + name).loc[0, 'control_common']) == wanted)
    check('strict-AF-neighbors', output('control_af_cutoff_neighbors')['control_common'].tolist() == [False, False, True])
    overlap = states('rank_gene_overlap_sequential', 'sort_and_rank_variants')
    sort_values = []
    for event, state in overlap:
        if event['phase'] in ['R5', 'R6'] and 'df' in state:
            column = next(column for column in state['df']['data'] if column['name'] == encode('sort_index'))
            value = column['values'][0]
            if not sort_values or sort_values[-1] != value:
                sort_values.append(value)
    check('sequential-gene-compatibility-floating-states', len(sort_values) >= 3)
    positions = output('rank_existing_columns_positions').columns.tolist()
    check('overwritten-columns-retain-position', [positions.index(name) for name in
          ['variant_rank', 'sort_index', 'pext_sort_index', 'control_common_index']] == [1, 2, 3, 4])
    check('only-two-named-stale-columns-removed', 'haplo_insuf_index' not in positions and 'zygosity_inheritance_mechanism_index' not in positions)
    pext = output('rank_pext_nan_inf_and_tokens').set_index('__oracle_row_position')
    check('PEXT-NaN-order-observed', pext.loc['input-row-0', 'pext_sort_index'] == 1.0 and pext.loc['input-row-1', 'pext_sort_index'] < 1.0)
    check('PEXT-infinities-clipped-later', pext.loc['input-row-2', 'pext_sort_index'] == 1.0 and pext.loc['input-row-3', 'pext_sort_index'] == .8)
    chosen = [state['pext_col_to_use'] for _, state in states('rank_tissue_first_missing_fallback_mean', 'sort_and_rank_variants')
              if 'pext_col_to_use' in state and state['pext_col_to_use'] != encode(None)]
    check('first-missing-tissue-falls-back-to-mean', bool(chosen) and all(item == encode('PEXT_MEAN') for item in chosen))
    check('does-not-try-second-tissue-or-validate-unused-cutoff', reports['rank_tissue_second_present_without_mean']['status'] == 'value')
    check('missing-ID-rank-zero', 0 in output('rank_group_nulls_missing_shared_ids')['variant_rank'].tolist())
    check('absent-ID-coordinate-stringification-positive', (output('rank_absent_ids_null_coordinate_spellings')['variant_rank'] > 0).all())
    parsed = output('production_parser_headered')
    check('production-parser-complete-row-dedup', len(parsed['anno_df']) == 2)
    check('production-PED-preserves-leading-zero-family', output('production_parser_leading_zero_family')['ped_df'].loc[0, '#FamilyID'] == '001')
    check('production-PED-preserves-NA-like-ID', output('production_parser_NA_sample')['ped_df'].loc[0, 'IndividualID'] == 'NA')
    check('resort-inferred-PED-retains-original-family-mismatch', reports['resort_parser_leading_zero_family']['status'] == 'error')
    check('resort-headerless-PED-retains-original-failure', reports['resort_parser_headerless']['status'] == 'error')
    check('resort-zero-byte-original-parser-error', reports['resort_input_zero_bytes'].get('exception', {}).get('type') == 'pandas.errors.EmptyDataError')
    result = dict(schema='priva.full-ranking-oracle-validation.v1', finished_utc=utc(), native_environment=native,
                  capture_sha256=args.capture_sha256, validation_harness_sha256=digest(Path(__file__)),
                  checks=checks, passed=len(checks) - len(failures), failed=failures, successful=not failures,
                  counts=captured['counts'], observer_exact=all(item['observer_exactly_equivalent'] for item in reports.values()),
                  qualification='oracle groundwork only; no Rust/full-family/criterion-assignment/installation acceptance',
                  scientific_acceptance=False, rust_qualification=False, production_backend_selected=False)
    (args.root / 'validation-harness.py').write_bytes(Path(__file__).read_bytes())
    write_json(args.root / 'validation.json', result)
    print(json.dumps({key: value for key, value in result.items() if key != 'checks'}))
    if failures:
        raise SystemExit(1)


def validate(args):
    native = native_environment()
    assert args.root.resolve().is_relative_to(D)
    attempt = args.root / 'validation-attempts' / (native['SLURM_JOB_ID'] + '.' + native['SLURM_STEP_ID'])
    attempt.mkdir(parents=True, exist_ok=False)
    (attempt / 'harness.py').write_bytes(Path(__file__).read_bytes())
    report = dict(started_utc=utc(), native_environment=native,
                  harness_sha256=digest(attempt / 'harness.py'), capture_sha256=args.capture_sha256,
                  successful=False)
    try:
        validate_core(args, native)
        report['successful'] = True
    except BaseException as error:
        report['failure'] = dict(type=typename(error), message=str(error))
        raise
    finally:
        report['finished_utc'] = utc()
        write_json(attempt / 'attempt.json', report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest='command', required=True)
    prepare_parser = command.add_parser('prepare')
    prepare_parser.add_argument('--repo', type=Path, required=True)
    prepare_parser.add_argument('--output', type=Path, required=True)
    capture_parser = command.add_parser('capture')
    capture_parser.add_argument('--root', type=Path, required=True)
    capture_parser.add_argument('--pins-sha256', required=True)
    validate_parser = command.add_parser('validate')
    validate_parser.add_argument('--root', type=Path, required=True)
    validate_parser.add_argument('--capture-sha256', required=True)
    args = parser.parse_args()
    {'prepare': prepare, 'capture': capture, 'validate': validate}[args.command](args)


if __name__ == '__main__':
    sys.dont_write_bytecode = True
    main()
