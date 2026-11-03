"""Explicit, unregistered Rust-owned PP3/BP4/BP7 procedure profile.

Rust owns the complete named procedures and their scientific control flow.
This host is a sequential object/effect service: actual getters, bound methods,
annotation definition, formatting/log callbacks and ownership happen only when
visited. Unrepresented operands execute ONE actual mature operation at that
frontier and are counted as mixed; the original criterion is never executed.

The earlier priva_rust_insilico.Bridge/insilico-columns-v1 profile is retained
as a separate transitional control. This class does not register or activate
itself in the caller. The new profile needs separate native qualification.
"""
import ast
import builtins
from collections import Counter
import copy
import csv
import hashlib
import inspect
import io
import operator
import os
from pathlib import Path
import struct
import subprocess
import time
from types import CodeType, FunctionType

import numpy as np
import pandas as pd

import acmg_pp3_bp4_bp7_insilico as _reference
import priva_rust_insilico as columns
import priva_rust_insilico_source as guards

PROFILE = 'insilico-procedures-v2'
TransportError = columns.TransportError
Unsupported = columns.Unsupported
_MISSING = object()
_FUNCTIONS = frozenset(('PP3_BP4_criteria', 'BP7_criteria', '_coerce_score_column'))
_ATTRIBUTES = frozenset(('str', 'contains', 'fillna', 'astype', 'isna', 'notna',
    'strip', 'map', 'max', 'sum', 'any', 'unique', 'get', 'index', 'warning', 'info'))
_GLOBALS = frozenset((('np', 'zeros'), ('np', 'logical_not'), ('np', 'asarray'),
    ('pd', 'Series'), ('pd', 'to_numeric'), ('pd', 'concat'),
    ('', 'int'), ('', 'str'), ('', 'sorted'), ('', 'logger')))
_OPERATIONS = frozenset((*columns._BINARY, 'invert', 'not', 'contains', 'strip',
    'fillna', 'isna', 'notna', 'as_str', 'map', 'max_columns', 'numeric', 'sum', 'any', 'unique',
    'sorted', 'zeros', 'assign', 'series_default', 'get', 'concat', 'asarray',
    'int', 'log'))


def _code_tree(code):
    """Explicit full fields; CodeType.__eq__ alone ignores some provenance."""
    def constant(value):
        if type(value) is CodeType: return ('code', _code_tree(value))
        if type(value) is tuple: return ('tuple', tuple(constant(item) for item in value))
        if type(value) is frozenset: return ('frozenset', frozenset(constant(item) for item in value))
        if type(value) is float: return ('float', struct.pack('>d', value))
        return (type(value).__name__, value)
    names = ('co_argcount', 'co_posonlyargcount', 'co_kwonlyargcount', 'co_nlocals',
        'co_stacksize', 'co_flags', 'co_code', 'co_names', 'co_varnames',
        'co_filename', 'co_name', 'co_qualname', 'co_firstlineno', 'co_linetable',
        'co_exceptiontable', 'co_freevars', 'co_cellvars')
    return (tuple((name, getattr(code, name)) for name in names),
            tuple(constant(value) for value in code.co_consts), tuple(code.co_positions()))


def _source_codes(reference):
    raw = guards.reference_guard(reference)
    compiled = compile(raw, str(reference.__file__), 'exec', dont_inherit=True)
    result = {code.co_name: code for code in compiled.co_consts if type(code) is CodeType}
    for name in ('PP3_BP4_criteria', 'BP7_criteria'):
        if _code_tree(getattr(reference, name).__code__) != _code_tree(result[name]):
            raise Unsupported('full original CodeType/source identity differs')
    helper = [code for code in result['PP3_BP4_criteria'].co_consts
              if type(code) is CodeType and code.co_name == '_coerce_score_column']
    if len(helper) != 1: raise Unsupported('original score helper declaration differs')
    result['_coerce_score_column'] = helper[0]
    return raw, result


def _rows(raw):
    try:
        result = list(csv.reader(io.StringIO(raw.decode('utf-8'), newline=''),
                                 delimiter='\t', strict=True))
    except (UnicodeError, csv.Error) as error:
        raise TransportError('corrupt procedure CSV/UTF-8') from error
    if not result or any(len(row) != 6 for row in result):
        raise TransportError('procedure frame schema differs')
    return result


def _bytes(rows):
    output = io.StringIO(newline='')
    writer = csv.writer(output, delimiter='\t', lineterminator='\n')
    writer.writerows(rows)
    result = output.getvalue().encode('utf-8')
    if not 0 < len(result) <= columns.MAX_FRAME:
        raise TransportError('procedure frame length differs')
    return result


def _number(value):
    if type(value) is not str or not value.isascii() or not value.isdigit():
        raise TransportError('invalid procedure ordinal/handle')
    number = int(value)
    if number >= 1 << 64: raise TransportError('procedure ordinal/handle overflow')
    return number


def _parts(rows, kind, *, allowed=True):
    result = []
    for row in rows:
        if row[0] != kind: continue
        if not allowed or row[1] != str(len(result)) or any(row[3:]):
            raise TransportError('nonconsecutive procedure detail')
        result.append(row[2])
    return result


def _binding_snapshot():
    from pandas.core.strings.accessor import StringMethods
    owners = ((pd.Series, ('__and__', '__or__', '__invert__', '__eq__', '__ne__',
        '__lt__', '__le__', '__gt__', '__ge__', 'fillna', 'astype', 'notna', 'isna',
        'map', 'sum', 'any', 'unique', '__array__')),
        (StringMethods, ('contains', 'strip')), (pd.DataFrame, ('max',)),
        (pd, ('to_numeric', 'concat', 'Series')), (np, ('logical_not', 'zeros', 'asarray')))
    snapshot = []
    for owner, names in owners:
        for name in names:
            function = getattr(owner, name)
            snapshot.append((owner, name, function, getattr(function, '__code__', None)))
    return snapshot


class _Host:
    """No policy branch, original-function call, or speculative operand read."""
    def __init__(self, reference, raw, codes, session):
        self.reference, self.raw, self.codes, self.session = reference, raw, codes, session
        self.objects, self.next_object = {}, 0
        self.attributes, self.globals = {}, {}
        self.sequence = 0
        self.pending = None
        self.source_error = None
        self.shims = {}
        self.bindings = _binding_snapshot()
        self.force_mature = False
        self.builtins = reference.PP3_BP4_criteria.__builtins__

    def put(self, value):
        handle = self.next_object
        self.next_object += 1
        self.objects[handle] = value
        return handle

    def get(self, handle):
        number = _number(handle)
        if number not in self.objects: raise TransportError('foreign mature object handle')
        return self.objects[number]

    def _at(self, site, function, args=(), kwargs=None):
        # The shim is explicitly NOT a genuine original CodeType. It only
        # preserves the actual warning/logger source frontier for one mature
        # effect on the same already evaluated objects.
        if site not in self.shims:
            name, line = site
            tree = ast.parse('def call(__function, __args, __kwargs):\n return __function(*__args, **__kwargs)\n')
            node = tree.body[0]
            node.name = name
            ast.increment_lineno(node, line - 2)
            compiled = compile(tree, str(self.reference.__file__), 'exec', dont_inherit=True)
            code = next(item for item in compiled.co_consts if type(item) is CodeType)
            self.shims[site] = FunctionType(code, vars(self.reference), name)
        return self.shims[site](function, args, {} if kwargs is None else kwargs)

    def _constant(self, site, value):
        def constants(code):
            for item in code.co_consts:
                if type(item) is tuple:
                    yield from item
                elif type(item) is not CodeType:
                    yield item
        def equal(item):
            if type(item) is not type(value): return False
            if type(item) is float: return struct.pack('>d', item) == struct.pack('>d', value)
            return item == value
        matches = [item for item in constants(self.codes[site[0]]) if equal(item)]
        if not matches:
            raise TransportError('Rust literal absent from pinned source constant pool')
        return matches[0]

    def _global(self, module, name):
        namespace = vars(self.reference)
        if module: return getattr(namespace[module], name)
        if name in namespace: return namespace[name]
        return self.builtins[name]

    def _same_bindings(self):
        return (self.reference.pd is pd and self.reference.np is np
            and all(vars(self.reference).get(name, self.builtins.get(name)) is getattr(builtins, name)
                    for name in ('len', 'sorted', 'int', 'str', 'float'))
            and all(getattr(owner, name) is function
                    and getattr(function, '__code__', None) is code
                    for owner, name, function, code in self.bindings))

    def _effect(self, site, action, args, text):
        if action == 'getitem' and len(args) == 2 and not text:
            return self._at(site, operator.getitem, args)
        if action == 'list' and not text: return list(args)
        if action == 'slice5' and len(args) == 1 and not text:
            return self._at(site, operator.getitem, (args[0], slice(None, 5)))
        if action in ('format', 'join_text') and len(text) == len(args) + 1:
            pieces = [text[0]]
            for index, value in enumerate(args):
                pieces.append(self._at(site, format, (value, '')) if action == 'format' else value)
                pieces.append(text[index + 1])
            # BUILD_STRING only accepts already formatted exact strings.
            if any(type(piece) is not str for piece in pieces):
                raise TransportError('formatted source message contains non-text')
            return ''.join(pieces)
        if action == 'define_score_helper' and len(args) == 1 and not text:
            tree = ast.parse(self.raw)
            outer = next(node for node in tree.body if type(node) is ast.FunctionDef
                         and node.name == 'PP3_BP4_criteria')
            helper = next(node for node in outer.body if type(node) is ast.FunctionDef)
            # Only the original definition/annotations are evaluated. The
            # original helper body is retained as an uncalled code object.
            wrapper = ast.parse('def _priva_definition(df):\n return _coerce_score_column\n').body[0]
            wrapper.body.insert(0, copy.deepcopy(helper))
            namespace = dict(vars(self.reference))
            exec(compile(ast.fix_missing_locations(ast.Module(body=[wrapper], type_ignores=[])),
                         str(self.reference.__file__), 'exec', dont_inherit=True), namespace)
            return namespace['_priva_definition'](args[0])
        if action == 'after_log' and not args and not text:
            try: guards.reference_guard(self.reference)
            except (OSError, ValueError, TypeError, AttributeError) as error:
                self.force_mature = True
                self.session.metrics['binding_changes'].append(type(error).__name__ + ': ' + str(error))
            return None
        raise TransportError('unknown mature-effect action/arity')

    def _call(self, details):
        if not details or details[0][0] != 'call': raise TransportError('missing operation descriptor')
        head = details[0]
        operation, target, function, receiver = head[1:5]
        if operation not in _OPERATIONS or head[5]: raise TransportError('unknown source operation')
        args, keywords, native = [], [], None
        section = 'argument'
        for record in details[1:]:
            kind = record[0]
            if kind == 'argument' and section == 'argument':
                if record[1] != str(len(args)) or any(record[3:]): raise TransportError('bad operation arguments')
                self.get(record[2]); args.append(record[2])
            elif kind == 'keyword' and section in ('argument', 'keyword'):
                section = 'keyword'
                if record[1] != str(len(keywords)) or any(record[4:]) or record[2] in dict(keywords):
                    raise TransportError('bad operation keywords')
                self.get(record[3]); keywords.append((record[2], record[3]))
            elif kind == 'native' and native is None and record[1] == '0' and not any(record[3:]):
                section = 'native'; native = record[2]
            else: raise TransportError('unknown/out-of-order operation detail')
        if target == 'bound':
            self.get(function)
            if receiver:
                obj = self.get(receiver)
                if not args or self.get(args[0]) is not obj:
                    raise TransportError('bound receiver/operands differ')
        elif target in (*columns._BINARY, 'invert', 'assignment'):
            if function or receiver or keywords: raise TransportError('foreign operator descriptor')
            if target in columns._BINARY and (operation != target or len(args) != 2):
                raise TransportError('binary operation arity differs')
            if target == 'invert' and (operation != 'invert' or len(args) != 1):
                raise TransportError('invert arity differs')
            if target == 'assignment' and (operation != 'assign' or len(args) != 3):
                raise TransportError('assignment arity differs')
        else: raise TransportError('foreign operation target')
        return (operation, target, function, receiver, tuple(args), tuple(keywords)), native

    def _mature(self, site, call):
        operation, target, function, receiver, handles, keyword_handles = call
        args = [self.get(handle) for handle in handles]
        kwargs = {name: self.get(handle) for name, handle in keyword_handles}
        if target in columns._BINARY:
            result = self._at(site, columns._BINARY[target], args)
        elif target == 'invert': result = self._at(site, operator.invert, args)
        elif target == 'assignment':
            self._at(site, operator.setitem, args)
            result = args[0]
        else:
            result = self._at(site, self.get(function), args[1:] if receiver else args, kwargs)
        self.session.metrics['mature_operations'][operation] += 1
        return result

    def _safe_operands(self, operands, operation):
        # This is operation-local encoding; it never scans the original frame,
        # future columns, PVS1 or confidence mapping at function entry.
        from priva_rust_pm4_bp3_packages import instance_guard
        for value in operands:
            if type(value) is pd.Series:
                instance_guard(value, series=True)
                axis = guards.axis_values(value._mgr.axes[0])
                if len(set(axis)) != len(axis) or len({type(item) for item in axis}) > 1:
                    raise Unsupported('mature duplicate/mixed Series axis operation')
            columns._value(value)
        if operation in ('and', 'or') and len(operands) == 2 and all(type(value) is pd.Series for value in operands):
            axes = [guards.axis_values(value._mgr.axes[0]) for value in operands]
            if axes[0] != axes[1]: raise Unsupported('mature logical Series alignment')
        return operands

    def _operands(self, call):
        operation, target, function, receiver, handles, keyword_handles = call
        args = [self.get(handle) for handle in handles]
        kwargs = {name: self.get(handle) for name, handle in keyword_handles}
        if self.force_mature or not self._same_bindings():
            raise Unsupported('mature binding changed at visited frontier')
        if target in columns._BINARY:
            # Python lists/tuples have their own comparison semantics; do not
            # turn PVS1 list < 3 into a vector computation that would succeed.
            if any(type(value) in (list, tuple) for value in args):
                raise Unsupported('mature sequence operator')
            return self._safe_operands(args, operation)
        if target == 'invert': return self._safe_operands(args, operation)
        if target == 'assignment':
            if (type(args[0]) is not np.ndarray or args[0].dtype != np.dtype('int64')
                    or args[0].ndim != 1 or not args[0].flags.writeable):
                raise Unsupported('mature output assignment ownership')
            self._safe_operands(args, operation)
            return args
        actual = self.get(function)
        owner = self.get(receiver) if receiver else None
        if operation in ('get', 'series_default', 'concat', 'asarray', 'int', 'log'):
            raise Unsupported('retained mature effect: ' + operation)
        if receiver:
            attribute = self.attributes.get(_number(function))
            if attribute is None or attribute[0] != _number(receiver):
                raise TransportError('bound method lacks actual attribute frontier')
            name = attribute[1]
            expected_name = {'as_str': 'astype', 'max_columns': 'max'}.get(operation, operation)
            if name != expected_name: raise TransportError('bound method/source operation differ')
            if type(owner) is pd.Series:
                # Refuse unknown extension/instance storage before inspecting
                # dtype equality or invoking a compatibility adapter. Unknown
                # dtype objects must not acquire an extra __eq__ callback.
                self._safe_operands([owner], operation)
            from pandas.core.strings.accessor import StringMethods
            if operation in ('contains', 'strip'):
                if type(owner) is not StringMethods or getattr(actual, '__func__', None) is not getattr(StringMethods, name):
                    raise Unsupported('custom mature string binding')
                if kwargs and not (operation == 'contains' and kwargs == {'regex': True}):
                    raise TransportError('source string options differ')
                operands = [owner._orig, *args[1:]]
            elif operation == 'as_str':
                if len(args) != 2 or args[1] is not str or kwargs:
                    raise Unsupported('mature string dtype binding')
                operands = [owner]
            elif operation == 'map':
                if (type(owner) is not pd.Series or len(args) != 2 or type(args[1]) is not dict
                        or kwargs != {'na_action': 'ignore'}):
                    raise Unsupported('mature callable/custom confidence map')
                operands = [owner, list(args[1]), list(args[1].values())]
            elif operation == 'max_columns':
                if type(owner) is not pd.DataFrame or len(args) != 1 or kwargs != {'axis': 1}:
                    raise Unsupported('mature horizontal maximum')
                if getattr(actual, '__func__', None) is not pd.DataFrame.max:
                    raise Unsupported('custom mature DataFrame maximum')
                from priva_rust_pm4_bp3 import _package_profile
                import priva_rust_pm4_bp3_packages as packages
                import json
                binding, _digest = _package_profile()
                expected = json.loads(guards.regular_bytes(binding))['source_files']
                context = packages.source_context({'pandas.core.frame'}, expected)
                packages.plain_function(actual.__func__, 'pandas.core.frame:DataFrame.max', context)
                labels = guards.axis_values(owner._mgr.axes[1])
                if len(set(labels)) != len(labels): raise Unsupported('mature horizontal duplicate index')
                index = owner._mgr.axes[1]
                operands = [pd.Series(value, index=index, copy=False)
                            for value in guards.raw_columns(owner)]
                if not operands: raise Unsupported('empty mature horizontal maximum')
            else:
                if type(owner) is not pd.Series or getattr(actual, '__func__', None) is not getattr(pd.Series, name):
                    raise Unsupported('custom mature Series method')
                if kwargs: raise Unsupported('mature method keyword profile')
                operands = args
                if operation == 'fillna':
                    if len(args) != 2: raise TransportError('source fillna arity')
                    dtype = owner.dtype
                    if type(args[1]) in (bool, np.bool_) and (dtype == np.dtype('float64') or type(dtype) is pd.Float64Dtype):
                        raise Unsupported('mature float/Boolean fill dtype/warning')
                    if dtype == np.dtype('O'):
                        raw = owner._mgr.blocks[0].values
                        missing = [value is None or value is pd.NA or type(value) in (float, np.float64) and np.isnan(value) for value in raw]
                        filled = [args[1] if absent else value for value, absent in zip(raw, missing)]
                        if any(missing) and all(type(value) in (bool, int, float, np.bool_, np.int64, np.float64) for value in filled):
                            raise Unsupported('mature object fill downcast warning')
        else:
            declaration = self.globals.get(_number(function))
            expected = {'numeric': ('pd', 'to_numeric'), 'not': ('np', 'logical_not'),
                        'zeros': ('np', 'zeros'), 'sorted': ('', 'sorted')}.get(operation)
            if expected is None or declaration != expected:
                raise TransportError('native global operation has foreign binding')
            if actual is not self._global(*expected): raise Unsupported('changed mature global binding')
            if operation == 'numeric':
                if kwargs != {'errors': 'coerce'} or len(args) != 1: raise TransportError('numeric source options')
            elif operation == 'zeros':
                if kwargs != {'dtype': int} or len(args) != 1: raise Unsupported('mature output dtype')
            elif kwargs: raise TransportError('unexpected source global keywords')
            operands = args
        return self._safe_operands(operands, operation)

    def _details(self, rows):
        objects = _parts(rows, 'object')
        text = _parts(rows, 'text')
        if any(row[0] not in ('object', 'text') for row in rows):
            raise TransportError('unknown effect detail')
        return [self.get(handle) for handle in objects], objects, text

    def service(self, rows):
        header, details = rows[0], rows[1:]
        if (header[0] != 'host' or header[1] != str(self.sequence) or header[5] != PROFILE
                or header[3] not in _FUNCTIONS):
            raise TransportError('host request sequence/profile/function differs')
        site = (header[3], _number(header[4]))
        if not 32 <= site[1] <= 138: raise TransportError('host source line outside procedure')
        action = header[2]
        self.session.metrics['host_effects'][action] += 1
        if action == 'release':
            if self.pending is not None or any(row[0] != 'object' for row in details):
                raise TransportError('release interleaved with operation/foreign detail')
            handles = [_number(value) for value in _parts(details, 'object')]
            if len(set(handles)) != len(handles) or any(value not in self.objects for value in handles):
                raise TransportError('foreign/duplicate released handle')
            # Do not first build a list of the referents: that would keep a
            # source temporary alive across its explicit drop frontier.
            for handle in handles:
                del self.objects[handle]
                self.attributes.pop(handle, None); self.globals.pop(handle, None)
            return 'released', []
        status, payload = 'object', []
        if action in ('prepare', 'commit', 'mature'):
            call, native = self._call(details)
            if action == 'prepare':
                if self.pending is not None or native is not None: raise TransportError('replayed prepare/native payload')
                self.session.metrics['source_operations'][call[0]] += 1
                reason, operands, encoded = None, None, None
                try:
                    operands = self._operands(call)
                    encoded = columns._request(self.sequence, call[0], operands).decode('utf-8')
                except (Unsupported, ValueError, TypeError, AttributeError) as error:
                    reason = str(error)
                self.pending = (site, call, operands, reason)
                if encoded is None:
                    self.session.metrics['mature_admissions'][call[0]] += 1
                    self.session.metrics['mature_sites'].append((site, call[0], reason))
                    return 'mature', []
                return 'prepared', [('operands', '0', encoded, '', '', '')]
            if self.pending is None or self.pending[:2] != (site, call):
                raise TransportError('foreign/replayed operation completion')
            _site, _call, operands, reason = self.pending
            self.pending = None
            if action == 'mature':
                if native is not None: raise TransportError('native payload on mature effect')
                result = self._mature(site, call)
            else:
                if native is None or operands is None or reason is not None:
                    raise TransportError('native commit without prepared operands')
                try:
                    result = columns._result(native.encode('utf-8'), self.sequence, operands, call[0])
                except columns._NativeUnsupported as error:
                    # The Rust driver sends a mature request after a kernel
                    # rejection. A commit MUST contain an ordinary result;
                    # legacy unsupported status here is a corrupt transport,
                    # not an original scientific exception/admission result.
                    raise TransportError('impossible unsupported native commit') from error
                if call[0] == 'assign':
                    original = self.get(call[4][0])
                    mask = self.get(call[4][1])
                    if type(mask) is pd.Series:
                        # This is the actual ndarray assignment's ownership
                        # conversion. Its result is used to verify this same
                        # operation, never a ceremonial getter/cache replay.
                        converted = self._at(site, np.asarray, (mask,))
                        expected = operands[1]._mgr.blocks[0].values
                        if type(expected) is pd.arrays.BooleanArray:
                            # A completed native assignment has a fully valid
                            # Boolean mask. Its canonical extension storage is
                            # two ndarray buffers, not an ndarray with tobytes.
                            # Keep the real source conversion above and verify
                            # its result against these already-admitted buffers.
                            if expected._mask.any():
                                raise TransportError('native assignment retained a missing Boolean mask')
                            expected = expected._data
                        if (type(expected) is not np.ndarray or expected.dtype != np.dtype('bool')
                                or type(converted) is not np.ndarray or converted.dtype != np.dtype('bool')
                                or converted.shape != expected.shape or converted.tobytes() != expected.tobytes()):
                            raise TransportError('assignment ownership conversion changed operand')
                    self._at(site, operator.setitem, (original, slice(None), result))
                    result = original
                self.session.metrics['native_operations'][call[0]] += 1
            return 'object', [('object', '0', str(self.put(result)), '', '', '')]
        if self.pending is not None: raise TransportError('effect interleaved with incomplete operation')
        if action == 'literal':
            if len(details) != 1 or details[0][0] != 'literal' or any(details[0][3:]):
                raise TransportError('literal detail schema')
            value = columns._uncell(details[0][1], details[0][2])
            result = self._constant(site, value)
        elif action == 'effect':
            if not details or details[0][0] != 'effect' or any(details[0][2:]):
                raise TransportError('effect action schema')
            args, _handles, text = self._details(details[1:])
            result = self._effect(site, details[0][1], args, text)
        else:
            args, handles, text = self._details(details)
            if action == 'column' and len(args) == len(text) == 1:
                key = self._constant(site, text[0])
                result = self._at(site, operator.getitem, (args[0], key))
            elif action == 'attribute' and len(args) == len(text) == 1 and text[0] in _ATTRIBUTES:
                result = self._at(site, getattr, (args[0], text[0]))
                handle = self.put(result)
                self.attributes[handle] = (_number(handles[0]), text[0])
                return 'object', [('object', '0', str(handle), '', '', '')]
            elif action == 'global' and not args and len(text) == 2 and tuple(text) in _GLOBALS:
                result = self._global(*text)
                handle = self.put(result); self.globals[handle] = tuple(text)
                return 'object', [('object', '0', str(handle), '', '', '')]
            elif action == 'length' and len(args) == 1 and not text:
                length = vars(self.reference).get('len', self.builtins['len'])
                result = self._at(site, length, args)
            elif action in ('truth', 'is_none') and len(args) == 1 and not text:
                result = self._at(site, bool, args) if action == 'truth' else args[0] is None
                return 'boolean', [('boolean', '0', str(int(result)), '', '', '')]
            else: raise TransportError('unknown source effect/arity')
        payload.append(('object', '0', str(self.put(result)), '', '', ''))
        return status, payload


class ProcedureSession(columns.NativeSession):
    """One binary child drives one whole function; no row child or memo cache."""
    def __init__(self, binary, digest, *, timeout=120):
        super().__init__(binary, digest, timeout=timeout)
        self.metrics = dict(processes=0, wire_requests=0, native_operations=Counter(),
            mature_operations=Counter(), mature_admissions=Counter(), mature_sites=[],
            source_operations=Counter(), host_effects=Counter(), binding_changes=[],
            transport_profile=PROFILE, cache_hits=0, original_entries=0,
            scientific_control='named Rust procedures', completed=False)

    def start(self):
        if self.process is not None: raise TransportError('procedure child already started')
        if (hashlib.sha256(guards.regular_bytes(self.binary)).hexdigest() != self.digest
                or not os.access(self.binary, os.X_OK)):
            raise TransportError('procedure executable identity differs')
        self.process = subprocess.Popen([str(self.binary), 'insilico-criteria', 'procedure-session'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env=dict(os.environ, POLARS_MAX_THREADS='1', RAYON_NUM_THREADS='1'), bufsize=0)
        os.set_blocking(self.process.stdin.fileno(), False)
        self.metrics['processes'] += 1

    def send(self, rows, deadline):
        raw = _bytes(rows)
        self._write(struct.pack('>I', len(raw)) + raw, deadline)

    def receive(self, deadline):
        length = struct.unpack('>I', self._read(4, deadline))[0]
        if not 0 < length <= columns.MAX_FRAME: raise TransportError('procedure response size')
        return _rows(self._read(length, deadline))

    def run(self, host, procedure, handles):
        self.start()
        deadline = time.monotonic() + self.timeout
        self.send([('run', procedure, *map(str, handles), PROFILE)], deadline)
        while True:
            rows = self.receive(deadline)
            if rows[0][0] == 'done':
                header = rows[0]
                if (header[1] != str(host.sequence) or header[3] != PROFILE
                        or header[4] != guards.REFERENCE_SHA or header[5] or host.pending is not None):
                    raise TransportError('procedure completion profile/sequence')
                if header[2] == 'source-error':
                    if len(rows) != 1 or host.source_error is None:
                        raise TransportError('foreign source exception acknowledgment')
                    self.finish()
                    raise host.source_error
                expected = 1 if procedure == 'BP7_criteria' else 2
                if header[2] != 'ok' or host.source_error is not None or len(rows) != expected + 1:
                    raise TransportError('procedure output/completion denominator')
                handles = _parts(rows[1:], 'object')
                if len(handles) != expected: raise TransportError('procedure output schema')
                result = tuple(host.get(handle) for handle in handles)
                self.finish()
                self.metrics['completed'] = True
                return result[0] if expected == 1 else result
            if host.source_error is not None: raise TransportError('operation requested after source error')
            try:
                status, payload = host.service(rows)
            except TransportError:
                raise
            except BaseException as error:
                host.source_error = error
                host.pending = None
                status, payload = 'source-error', []
            self.send([('host-result', str(host.sequence), status, PROFILE, '', ''), *payload], deadline)
            host.sequence += 1
            self.metrics['wire_requests'] += 1
            # Existing individual native-fault timeout remains 120 seconds;
            # this is not a whole-case/scientific time budget.
            deadline = time.monotonic() + self.timeout

    def finish(self):
        process = self.process
        process.stdin.close()
        try: code = process.wait(timeout=5)
        except subprocess.TimeoutExpired as error: raise TransportError('procedure child close timeout') from error
        remainder, diagnostic = process.stdout.read(), process.stderr.read(65537)
        process.stdout.close(); process.stderr.close()
        if code or remainder or diagnostic: raise TransportError('procedure child failed/extra output')
        if hashlib.sha256(guards.regular_bytes(self.binary)).hexdigest() != self.digest:
            raise TransportError('procedure executable changed during call')


class ProcedureBridge:
    """Explicit complete-function Rust control; no original-function fallback."""
    def __init__(self, *, binary, binary_sha256, reference=_reference):
        self.binary, self.binary_sha256, self.reference = Path(binary), binary_sha256, reference
        self.last_metrics = None

    def call(self, procedure, df, *args, **kwargs):
        if procedure not in ('PP3_BP4_criteria', 'BP7_criteria'): raise ValueError('unknown insilico procedure')
        with columns._LOCK:
            raw, codes = _source_codes(self.reference)
            # Binding has no DataFrame/argument-value observation. No eager
            # frame/PVS/confidence admission or target-function execution.
            bound = inspect.signature(getattr(self.reference, procedure)).bind(df, *args, **kwargs)
            bound.apply_defaults()
            session = ProcedureSession(self.binary, self.binary_sha256)
            host = _Host(self.reference, raw, codes, session)
            handles = (host.put(df), host.put(bound.arguments.get('pvs1_criteria')),
                       host.put(bound.arguments.get('high_confidence_status')))
            try:
                return session.run(host, procedure, handles)
            finally:
                if session.process is not None and session.process.poll() is None:
                    session.close(failed=True)
                self.last_metrics = copy.deepcopy(session.metrics)
