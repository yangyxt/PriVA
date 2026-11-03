"""Exact source/callable proof for the optional reported-evidence facade.

The four historical public functions and the complete interpreter are pinned.
The host compiles only a finite named inventory of mature expressions. Original
functions are never rewritten; unrelated changes, duplicates, decorators and
asynchronous definitions are rejected.
"""
import ast
import builtins
import hashlib
import inspect
from pathlib import Path
import stat
from types import CodeType, FunctionType, ModuleType


SOURCE_SHA = {
    'acmg_ps3_bs3_pp5_bp6_reported': '6a60267f56378f30d4a916fd139d953d1629429237769651f6017a7d37e3eb4f',
    'mavedb_interpreter': 'ae50157d3b9ab36526b26be1e333a4d802ccd3733e4b3dfae3ee48b45a155c17',
}
PUBLIC = ('mavedb_interpretation_per_row', 'mavedb_score_interpretation',
          'PS3_BS3_criteria', 'PP5_BP6_criteria')
HELPERS = ('_to_log2', 'validate_spec', 'interpret_score',
           'parse_neutral_band_from_metadata', 'pick_interpreter_from_metadata')
FAMILIES = {
    'Stability_K50': ('StabilityK50Subtype', 'StabilityK50Interpreter', ('LOG10_K50_STANDARDIZED', 'LOG10_K50_RELATIVE', 'K50_RAW')),
    'Stability_DG': ('StabilityDGSubtype', 'StabilityDGInterpreter', ('DDG_KCAL_MOL', 'DG_KCAL_MOL')),
    'Abundance': ('AbundanceSubtype', 'AbundanceInterpreter', ('VAMP_SEQ', 'CELL_SURFACE_TRAFF', 'TOTAL_EXPRESSION')),
    'Activity': ('ActivitySubtype', 'ActivityInterpreter', ('ENZYME_ACTIVITY_REL', 'REPORTER_ACTIVITY_REL')),
    'Binding': ('BindingSubtype', 'BindingInterpreter', ('KD', 'PKD', 'KA', 'KD_APPARENT')),
    'Fitness': ('FitnessSubtype', 'FitnessInterpreter', ('LOG2_ENRICHMENT',)),
    'Splicing': ('SplicingSubtype', 'SplicingInterpreter', ('DELTA_PSI_PERCENT', 'DELTA_PSI_FRACTION')),
    'Electrophysiology': ('EphysSubtype', 'EphysInterpreter', ('V1_2', 'CURRENT_DENSITY')),
    'Complementation': ('ComplementationSubtype', 'ComplementationInterpreter', ('YEAST', 'MAMMALIAN_CELL')),
    'Localization': ('LocalizationSubtype', 'LocalizationInterpreter', ('CELLULAR_LOCALIZATION',)),
    'Aggregation': ('AggregationSubtype', 'AggregationInterpreter', ('AGGREGATION_METRIC',)),
    'Interaction': ('InteractionSubtype', 'InteractionInterpreter', ('TWO_HYBRID',)),
    'MPRA': ('MPRASubtype', 'MPRAInterpreter', ('REPORTER_ACTIVITY',)),
}


def package_effect_guard():
    """Canonical additional mature materializers used by native mask values.

    The shared fixed-profile guard already protects the operations bypassed by
    Rust (str.contains, comparison/logical/invert), constructors and closures.
    These named source declarations protect result labelling and alignment;
    no mutable runtime callable is captured as its own expected identity.
    """
    from priva_rust_pm4_bp3_packages import plain_function, slot, source_context
    hashes = {
        'pandas.core.series': '0ba5405fc2a49ac5bed13e3274832ac3b2e8aee15ae58cb5e2172bbe64d1f0d0',
        'pandas.core.ops.common': '8d57ff2ffa0de9b29c50eb87e8581a28ec75f2c7bd0b71e5d893ddd14a23e2de',
        'pandas.core.indexes.base': '2d50f8000608294d074d8bf14a9cdc489f4bdd32ba5b81f233bcd00c44c3ff24',
        'pandas.core.indexes.range': 'aade484b66dab639ce1def7450ae63112ee9661825a695bffb5c2278b959cac0',
        'pandas.core.dtypes.missing': '04fcdb9abeceee28668cb284f40df5724e7800d8c0b6ba7cf9d553db4311e5f4',
    }
    context = source_context(set(hashes), hashes)
    import pandas as pd
    from pandas.core.indexes.base import Index
    from pandas.core.indexes.range import RangeIndex
    from pandas.core.arrays.floating import FloatingArray
    if pd.Index is not Index or pd.RangeIndex is not RangeIndex or pd.arrays.FloatingArray is not FloatingArray:
        raise ValueError('additional public pandas type binding changed')
    for cls in (Index, RangeIndex, FloatingArray):
        if inspect.getattr_static(cls, '__getattribute__') is not object.__getattribute__:
            raise ValueError('custom pandas storage/index getter')
    for label in ('pandas.core.series:Series._construct_result',
                  'pandas.core.series:Series.map', 'pandas.core.series:Series.__setitem__',
                  'pandas.core.ops.common:get_op_result_name', 'pandas.core.ops.common:_maybe_match_name',
                  'pandas.core.indexes.base:Index.equals', 'pandas.core.indexes.range:RangeIndex.equals'):
        plain_function(slot(label), label, context)
    from functools import partial
    import numpy as np
    from numpy.core import _multiarray_umath
    from pandas.core.dtypes import missing
    plain_function(missing.isna, 'pandas.core.dtypes.missing:isna', context)
    if pd.isna is not missing.isna or missing.INF_AS_NA is not False:
        raise ValueError('nondefault pandas missingness binding')
    actual = missing._isna
    if type(actual) is partial:
        if actual.args or actual.keywords != {'inf_as_na': False} or actual.__dict__:
            raise ValueError('nondefault pandas missingness closure')
        actual = actual.func
    plain_function(actual, 'pandas.core.dtypes.missing:_isna', context)
    if (missing.nan_checker is not _multiarray_umath.isnan or np.isnan is not _multiarray_umath.isnan
            or type(np.isnan) is not np.ufunc or np.isnan.__name__ != 'isnan'):
        raise ValueError('pandas/NumPy missingness operator changed')


def source_bytes(module):
    if type(module) is not ModuleType or module.__name__ not in SOURCE_SHA:
        raise ValueError('noncanonical reported-evidence module')
    path = Path(module.__file__)
    before = path.lstat()
    if not stat.S_ISREG(before.st_mode):
        raise ValueError('reported-evidence source is not regular')
    raw = path.read_bytes()
    after = path.lstat()
    fields = ('st_dev', 'st_ino', 'st_mode', 'st_size', 'st_mtime_ns', 'st_ctime_ns')
    if (any(getattr(before, name) != getattr(after, name) for name in fields)
            or hashlib.sha256(raw).hexdigest() != SOURCE_SHA[module.__name__]):
        raise ValueError('reported-evidence source changed')
    return raw


def _codes(code):
    yield code
    for value in code.co_consts:
        if type(value) is CodeType:
            yield from _codes(value)


def _function(function, module, code):
    if (type(function) is not FunctionType or function.__globals__ is not module.__dict__
            or function.__closure__ is not None or function.__code__ != code
            or function.__kwdefaults__ is not None):
        raise ValueError('actual reported-evidence callable changed')
    if function.__builtins__ is not builtins.__dict__:
        raise ValueError('custom reported-evidence builtins')


def _definitions(raw):
    tree = ast.parse(raw)
    result = {}
    for node in tree.body:
        if type(node) in (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef):
            if node.name in result:
                raise ValueError('duplicate reported-evidence definition')
            if type(node) is ast.AsyncFunctionDef or node.decorator_list:
                # AssaySpec's single source dataclass decorator is recognized
                # separately; no public/helper decorator is normalized away.
                if not (type(node) is ast.ClassDef and node.name == 'AssaySpec'
                        and len(node.decorator_list) == 1
                        and type(node.decorator_list[0]) is ast.Name
                        and node.decorator_list[0].id == 'dataclass'):
                    raise ValueError('unexpected reported-evidence decorator/async definition')
            result[node.name] = node
    return tree, result


def _enum_members(enum):
    import enum as library
    if type(enum) is not library.EnumType:
        raise ValueError('custom subtype metaclass')
    fields = type.__getattribute__(enum, '__dict__')
    if '_value_' in fields or '_name_' in fields:
        raise ValueError('custom enum value descriptor')
    members = fields.get('_member_map_')
    if type(members) is not dict or any(type(key) is not str for key in members):
        raise ValueError('custom enum member mapping')
    if inspect.getattr_static(enum, '__getattribute__') is not object.__getattribute__:
        raise ValueError('custom enum attribute access')
    for key, item in members.items():
        if type(item) is not enum:
            raise ValueError('custom enum member type')
        data = object.__getattribute__(item, '__dict__')
        if type(data.get('_value_')) is not str or data['_value_'] != key or data.get('_name_') != key:
            raise ValueError('enum name/value changed')
    return members


def _spec_constructor(module):
    """Bind the exact generated initializer and ordinary attribute effects.

    This is the source-derived dataclasses initializer body for the declared
    declared fields. It is compiled only; no AssaySpec instance is constructed.
    """
    cls = module.__dict__.get('AssaySpec')
    if type(cls) is not type:
        raise ValueError('custom AssaySpec metaclass')
    fields = type.__getattribute__(cls, '__dict__')
    names = ('assay_family', 'subtype', 'directionality', 'scale', 'baseline',
             'units', 'neutral_band', 'flexibility_requirement')
    if any(name in fields for name in names[:-1]) or fields.get(names[-1]) is not False:
        raise ValueError('AssaySpec field descriptor/default changed')
    if (inspect.getattr_static(cls, '__getattribute__') is not object.__getattribute__
            or inspect.getattr_static(cls, '__setattr__') is not object.__setattr__):
        raise ValueError('AssaySpec getter/setter changed')
    parameters = ','.join(('self', *names[:-1], names[-1] + '=False'))
    text = 'def __create_fn__():\n def __init__(' + parameters + '):\n'
    text += ''.join('  self.' + name + '=' + name + '\n' for name in names)
    text += ' return __init__\n'
    expected = [code for code in _codes(compile(text, '<string>', 'exec')) if code.co_name == '__init__']
    actual = fields.get('__init__')
    if (len(expected) != 1 or type(actual) is not FunctionType or actual.__code__ != expected[0]
            or actual.__globals__ is not module.__dict__ or actual.__builtins__ is not builtins.__dict__
            or actual.__defaults__ != (False,) or actual.__kwdefaults__ is not None or actual.__closure__ is not None):
        raise ValueError('AssaySpec generated initializer changed')


def guard(reported, interpreter):
    """Verify code, globals, enum/registry bindings and the original defaults.

    No frame, Series, AssaySpec getter, logger method, or scientific helper is
    called. The source code is compiled, never imported into a replacement
    namespace as an oracle.
    """
    import math
    import numpy as np
    import pandas as pd
    import logging
    raws = {module.__name__: source_bytes(module) for module in (reported, interpreter)}
    from types import BuiltinFunctionType
    for name in ('log', 'isnan'):
        function = getattr(math, name)
        if (type(function) is not BuiltinFunctionType or function.__self__ is not math
                or function.__name__ != name or function.__module__ != 'math'):
            raise ValueError('interpreter libm function changed')
    import struct
    if type(np.nan) is not float or struct.pack('>d', np.nan).hex() != '7ff8000000000000':
        raise ValueError('source NumPy NaN binding changed')
    import enum as enum_library
    enum_path = Path(enum_library.__file__)
    enum_raw = enum_path.read_bytes()
    if hashlib.sha256(enum_raw).hexdigest() != 'c8ead615c159598370295649eb296819ad4b40d50b200c4fec2d4269bf7af9ae':
        raise ValueError('enum package source changed')
    value_property = inspect.getattr_static(enum_library.Enum, 'value')
    if type(value_property) is not enum_library.property:
        raise ValueError('enum value descriptor changed')
    getter = object.__getattribute__(value_property, '__dict__').get('fget')
    expected = [code for code in _codes(compile(enum_raw, str(enum_path), 'exec')) if code.co_qualname == 'Enum.value']
    if len(expected) != 1 or type(getter) is not FunctionType or getter.__code__ != expected[0]:
        raise ValueError('enum value API changed')
    # CPython3.12's exact pinned EnumType has no __getattr__; absence is part
    # of its profile, not a request for an invented compatibility method.
    if enum_library.EnumType.__dict__.get('__getattr__') is not None or inspect.getattr_static(enum_library.EnumType, '__getattribute__') is not type.__getattribute__:
        raise ValueError('enum attribute lookup API changed')
    for name in ('__getitem__',):
        method = enum_library.EnumType.__dict__.get(name)
        expected = [code for code in _codes(compile(enum_raw, str(enum_path), 'exec')) if code.co_qualname == 'EnumType.' + name]
        if len(expected) != 1 or type(method) is not FunctionType or method.__code__ != expected[0]:
            raise ValueError('enum parsing API changed')
    builtin_types = ('str', 'float', 'bool', 'int', 'set', 'list', 'dict', 'tuple', 'enumerate',
                     'Exception', 'ValueError', 'TypeError', 'IndexError', 'KeyError', 'AttributeError', 'SyntaxError')
    for name in builtin_types:
        value = getattr(builtins, name)
        if (type(value) is not type or value.__module__ != 'builtins' or value.__name__ != name
                or not value.__flags__ & (1 << 8)):
            raise ValueError('noncanonical immutable builtin class')
    for name in ('len', 'isinstance', 'hasattr'):
        value = getattr(builtins, name)
        if (type(value) is not BuiltinFunctionType or value.__module__ != 'builtins'
                or value.__name__ != name or value.__self__ is not builtins):
            raise ValueError('noncanonical builtin function')
    _spec_constructor(interpreter)
    for module, names in ((reported, PUBLIC), (interpreter, HELPERS)):
        raw = raws[module.__name__]
        _tree, definitions = _definitions(raw)
        compiled = compile(raw, module.__file__, 'exec')
        codes = {code.co_qualname: code for code in _codes(compiled)}
        for name in names:
            if name not in definitions or type(definitions[name]) is not ast.FunctionDef:
                raise ValueError('missing unique synchronous source function')
            _function(module.__dict__.get(name), module, codes[name])
        for name in (*builtin_types, 'len', 'isinstance', 'hasattr'):
            expected = getattr(builtins, name)
            if module.__dict__.get(name, expected) is not expected:
                raise ValueError('reported-evidence builtin global changed')
        if module.__dict__.get('np') is not np:
            raise ValueError('reported-evidence NumPy binding changed')
        if module is reported:
            if (module.__dict__.get('pd') is not pd or module.__dict__.get('logging') is not logging
                    or module.__dict__.get('pick_interpreter_from_metadata') is not interpreter.pick_interpreter_from_metadata
                    or type(module.__dict__.get('logger')) is not logging.Logger
                    or object.__getattribute__(module.logger, '__dict__').get('name') != module.__name__):
                raise ValueError('reported-evidence pandas/interpreter/logger binding changed')
        else:
            if module.__dict__.get('math') is not math:
                raise ValueError('interpreter math binding changed')
            if (type(module.__dict__.get('REGISTRY')) is not dict
                    or type(module.__dict__.get('FAMILY_TO_SUBTYPE_ENUM')) is not dict):
                raise ValueError('custom interpreter registry/family mapping')
            for family, (enum_name, class_name, members) in FAMILIES.items():
                enum = module.__dict__.get(enum_name)
                actual_members = _enum_members(enum)
                if (module.FAMILY_TO_SUBTYPE_ENUM.get(family) is not enum
                        or list(actual_members) != list(members)):
                    raise ValueError('interpreter subtype enum changed')
                cls = module.__dict__.get(class_name)
                value = module.REGISTRY.get(family)
                if (type(cls) is not type or type(value) is not cls
                        or type(object.__getattribute__(value, '__dict__')) is not dict
                        or object.__getattribute__(value, '__dict__')
                        or type.__getattribute__(cls, '__bases__') != (module.BaseInterpreter,)
                        or inspect.getattr_static(cls, '__getattribute__') is not object.__getattribute__):
                    raise ValueError('interpreter registry instance changed')
                method = cls.__dict__.get('interpret')
                _function(method, module, codes[class_name + '.interpret'])
                expected_defaults = (None,) if family == 'Electrophysiology' else None
                if method.__defaults__ != expected_defaults:
                    raise ValueError('interpreter method defaults changed')
            if set(module.REGISTRY) != set(FAMILIES) or set(module.FAMILY_TO_SUBTYPE_ENUM) != set(FAMILIES):
                raise ValueError('interpreter registry/family set changed')
            for name in ('standardize', 'call_from_band', 'interpret'):
                method = module.BaseInterpreter.__dict__.get(name)
                _function(method, module, codes['BaseInterpreter.' + name])
                if method.__defaults__ is not None:
                    raise ValueError('base interpreter defaults changed')
            for name, values in (('Scale', ('CENTERED_0', 'WT_EQ_1', 'RAW', 'PERCENT', 'FRACTION')),
                                 ('Directionality', ('HIGHER_IS_FUNCTION', 'LOWER_IS_FUNCTION', 'NON_MONOTONIC'))):
                enum = module.__dict__.get(name)
                if list(_enum_members(enum)) != list(values):
                    raise ValueError('interpreter direction/scale enum changed')
    defaults = reported.PS3_BS3_criteria.__defaults__
    expected = {'practice_guideline': 4, 'reviewed_by_expert_panel': 3,
                'criteria_provided,_multiple_submitters,_no_conflicts': 2}
    if (type(defaults) is not tuple or len(defaults) != 2 or type(defaults[0]) is not str or defaults[0] != ''
            or type(defaults[1]) is not dict or any(type(key) is not str for key in defaults[1])
            or list(defaults[1].items()) != list(expected.items())
            or any(type(value) is not int for value in defaults[1].values())):
        raise ValueError('PS3/BS3 source defaults changed')
    for name in PUBLIC:
        if name != 'PS3_BS3_criteria' and getattr(reported, name).__defaults__ is not None:
            raise ValueError('reported public defaults changed')
    for name in HELPERS:
        expected = (None,) if name == 'interpret_score' else None
        if getattr(interpreter, name).__defaults__ != expected:
            raise ValueError('interpreter helper defaults changed')
    return raws


# Each retained expression is a mature operation, not an executable copy of a
# PriVA function/branch. Rust selects every next effect and publication. Exact
# whole-file pins above bind all original control and the source positions;
# this finite table additionally rejects an expression absent from that source.
EFFECTS = {
    'row.missing': (PUBLIC[0], "pd.isna(row.get('MaveDB_urn', np.nan))"),
    'row.text.MaveDB_high_conf': (PUBLIC[0], "str(row.get('MaveDB_high_conf', ''))"),
    'row.text.MaveDB_pvalue': (PUBLIC[0], "str(row.get('MaveDB_pvalue', ''))"),
    'row.text.MaveDB_score': (PUBLIC[0], "str(row.get('MaveDB_score', ''))"),
    'row.text.MaveDB_urn': (PUBLIC[0], "str(row.get('MaveDB_urn'))"),
    'row.urn-values': (PUBLIC[0], "metadata_df['urn'].values"),
    'row.metadata-column': (PUBLIC[0], "metadata_df['urn']"),
    'row.rationale': (PUBLIC[0], "meta_row['rationale']"),
    'row.family': (PUBLIC[0], "meta_row.get('assay_family', 'Unknown')"),
    'row.subtype': (PUBLIC[0], "meta_row.get('subtype_value', 'Unknown')"),
    'row.float-score': (PUBLIC[0], 'float(score)'),
    'row.float-pvalue': (PUBLIC[0], 'float(pvalue)'),
    'row.call': (PUBLIC[0], "result.get('call', 'Unknown')"),
    'row.interpreter-ps3': (PUBLIC[0], "result['MaveDB_PS3']"),
    'row.interpreter-bs3': (PUBLIC[0], "result['MaveDB_BS3']"),
    'score.select.PS3': (PUBLIC[1], "df.loc[df['MaveDB_PS3'], ['MaveDB_pvalue', 'MaveDB_score', 'MaveDB_high_conf', 'MaveDB_score_interpretation']]"),
    'score.select.BS3': (PUBLIC[1], "df.loc[df['MaveDB_BS3'], ['MaveDB_pvalue', 'MaveDB_score', 'MaveDB_high_conf', 'MaveDB_score_interpretation']]"),
    'criteria.read-table': (PUBLIC[2], 'pd.read_table(mavedb_metadata_tsv, low_memory=False)'),
    'criteria.text.Pathogenic': (PUBLIC[2], "df['CLNSIG'].fillna('')"),
    'criteria.text.Benign': (PUBLIC[2], "df['CLNSIG'].fillna('')"),
    'criteria.ranks.Pathogenic': (PUBLIC[2], "df['CLNREVSTAT'].map(high_confidence_status)"),
    'criteria.ranks.Benign': (PUBLIC[2], "df['CLNREVSTAT'].map(high_confidence_status)"),
    'criteria.score.PS3': (PUBLIC[2], "df['MaveDB_PS3']"),
    'criteria.score.BS3': (PUBLIC[2], "df['MaveDB_BS3']"),
    'reported.text.athogenic': (PUBLIC[3], "df['CLNSIG'].fillna('')"),
    'reported.text.enign': (PUBLIC[3], "df['CLNSIG'].fillna('')"),
    'pick.urn-column': ('pick_interpreter_from_metadata', "metadata_df['urn']"),
    'pick.dict-lookup': ('pick_interpreter_from_metadata', 'metadata_df.get(urn)'),
    'band.initial-text': ('parse_neutral_band_from_metadata', "str(meta.get('neutral_bands', '[NA, NA]')).strip(' []')"),
}

# Repeated source expressions have distinct operation sites. Whole-file hashes
# protect these physical locations; never silently take the first occurrence.
EFFECT_LINES = {
    'row.metadata-column': 61,
    'criteria.text.Pathogenic': 218, 'criteria.ranks.Pathogenic': 218,
    'criteria.text.Benign': 222, 'criteria.ranks.Benign': 222,
    'reported.text.athogenic': 251, 'reported.text.enign': 252,
}


def compile_effects(raws, filenames):
    """Source-only compilation of the finite mature-expression inventory."""
    modules = {}
    for name, raw in raws.items():
        if hashlib.sha256(raw).hexdigest() != SOURCE_SHA.get(name):
            raise ValueError('effect source changed')
        _tree, definitions = _definitions(raw)
        modules[name] = definitions
    if set(modules) != set(SOURCE_SHA):
        raise ValueError('incomplete reported source unit')
    compiled = {}
    for key, (function_name, text) in EFFECTS.items():
        module_name = 'acmg_ps3_bs3_pp5_bp6_reported' if function_name in PUBLIC else 'mavedb_interpreter'
        node = modules[module_name].get(function_name)
        if type(node) is not ast.FunctionDef:
            raise ValueError('missing effect source function')
        expression = ast.parse(text, mode='eval').body
        signature = ast.dump(expression, include_attributes=False)
        matches = [part for part in ast.walk(node) if isinstance(part, ast.expr)
                   and ast.dump(part, include_attributes=False) == signature]
        if key in EFFECT_LINES:
            matches = [part for part in matches if part.lineno == EFFECT_LINES[key]]
        if not matches:
            raise ValueError('mature effect absent from original source: ' + key)
        if len(matches) != 1:
            raise ValueError('ambiguous mature effect source site: ' + key)
        compiled[key] = compile(ast.Expression(copy_location(matches[0])), filenames[module_name], 'eval')
    return compiled


def copy_location(expression):
    # Compile an independent expression node: never mutate the original AST or
    # a measured callable to observe its publications.
    import copy
    return ast.fix_missing_locations(copy.deepcopy(expression))
