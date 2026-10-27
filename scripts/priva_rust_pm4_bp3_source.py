"""Strict source and package admission for optional PM4/BP3 procedures."""
import ast
import builtins
from functools import lru_cache
import hashlib
import importlib
import inspect
import json
from pathlib import Path
import stat
from types import CodeType, FunctionType

import numpy as np
import pandas as pd
import priva_rust_pm4_bp3_packages as package_guards

REFERENCE_SHA = '08ef5af382801363727c6aec4915ca5751267067b43ef7996e3715e3fd1bc515'


def regular_bytes(path):
    path = Path(path)
    before = path.lstat()
    if not stat.S_ISREG(before.st_mode): raise ValueError('nonregular source')
    raw = path.read_bytes()
    after = path.lstat()
    fields = ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns', 'st_mode')
    if any(getattr(before, key) != getattr(after, key) for key in fields): raise ValueError('source race')
    return raw


def reference_guard(module):
    raw = regular_bytes(module.__file__)
    if hashlib.sha256(raw).hexdigest() != REFERENCE_SHA: raise ValueError('original source changed')
    compiled = compile(raw, str(module.__file__), 'exec')
    for name, defaults in (('find_overlaps_bedtools_efficient', ('any',)), ('PM4_BP3_criteria', None)):
        definitions = [node for node in ast.parse(raw).body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name]
        if len(definitions) != 1 or type(definitions[0]) is not ast.FunctionDef or definitions[0].decorator_list:
            raise ValueError('ambiguous original function')
        codes = [code for code in compiled.co_consts if type(code) is CodeType and code.co_name == name]
        function = getattr(module, name)
        if (type(function) is not FunctionType or len(codes) != 1 or function.__code__ != codes[0]
                or function.__globals__ is not module.__dict__ or function.__closure__ is not None
                or function.__kwdefaults__ is not None):
            raise ValueError('actual original callable changed')
        actual_defaults = function.__defaults__
        if defaults is None:
            if actual_defaults is not None: raise ValueError('original defaults changed')
        elif (type(actual_defaults) is not tuple or len(actual_defaults) != 1
                or type(actual_defaults[0]) is not str or actual_defaults[0] != 'any'):
            raise ValueError('original defaults changed')
        for name in ('len', 'set', 'int'):
            package_guards.builtin_guard(function.__builtins__, name)
    if module.__dict__.get('np') is not np or module.__dict__.get('pd') is not pd:
        raise ValueError('original library binding changed')
    if np.zeros is not np.core.multiarray.zeros: raise ValueError('NumPy zeros binding changed')
    for name in ('len', 'set', 'int'):
        package_guards.builtin_guard(module.__dict__, name)


def executable(code):
    # Match all recursive executable fields. Only debugger line/column ranges
    # are excluded for installed pandas build-interpreter compatibility. PriVA
    # functions above still require full CodeType equality.
    names = ('co_argcount', 'co_posonlyargcount', 'co_kwonlyargcount', 'co_nlocals',
        'co_stacksize', 'co_flags', 'co_code', 'co_names', 'co_varnames', 'co_filename',
        'co_name', 'co_qualname', 'co_firstlineno', 'co_exceptiontable', 'co_freevars', 'co_cellvars')
    return tuple(getattr(code, name) for name in names), tuple(executable(item) if type(item) is CodeType else item for item in code.co_consts)


def codes(code):
    yield code
    for item in code.co_consts:
        if type(item) is CodeType: yield from codes(item)


@lru_cache(maxsize=16)
def compiled(raw, filename): return compile(raw, filename, 'exec')


def defaults(value):
    from pandas._libs import lib
    if value is lib.no_default and value is lib._NoDefault.no_default:
        return ['pandas._libs.lib._NoDefault.no_default']
    if type(value) in (type(None), str, bool, int): return [type(value).__name__, value]
    if type(value) is tuple: return ['tuple', [defaults(item) for item in value]]
    if type(value) is dict: return ['dict', [[defaults(key), defaults(item)] for key, item in value.items()]]
    raise ValueError('unsupported package default')


def apis():
    return [(label, declaration, mode) for label, declaration, mode in package_guards.SPECS]


def freeze_package_bindings(path):
    return package_guards.freeze(path)


def package_guard(binding, digest):
    package_guards.guard(binding, digest)


def frame_instance_guard(frame):
    package_guards.instance_guard(frame)


def series_instance_guard(series):
    package_guards.instance_guard(series, series=True)


def logger_instance_guard(logger):
    package_guards.logger_instance_guard(logger)


def plain(value, active=None):
    if type(value) in (type(None), str, bool, int, float): return True
    if type(value) not in (list, tuple, dict): return False
    active = set() if active is None else active
    if id(value) in active: return False
    active.add(id(value))
    items = value.items() if type(value) is dict else enumerate(value)
    result = all(plain(key, active) and plain(item, active) for key, item in items)
    active.remove(id(value))
    return result
