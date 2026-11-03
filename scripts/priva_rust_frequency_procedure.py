"""Demand-driven host effects for named native frequency procedures.

Rust chooses every source branch, getter, scalar/table operation and return.
This module transports existing objects and performs only the requested mature
effect. It never lowers, invokes or replays an original scientific body.

The first complete native program is gnomAD_rare_AF. This module is deliberately
unregistered. A separately pinned binary and package binding are required;
adding this source does not grant scientific qualification or default use.
"""
import csv
import hashlib
import io
import operator
import os
from pathlib import Path
import re
import struct
import subprocess
import time
from types import BuiltinFunctionType, CodeType, MethodType

import numpy as np
import pandas as pd

import priva_rust_frequency as columns
import priva_rust_frequency_packages as packages
import priva_rust_pm4_bp3_packages as canonical

PROFILE = 'frequency-native-procedure-ef2e701-v1'
PROCEDURE = 'gnomAD_rare_AF'
SOURCE_SHA256 = 'ef2e70115ddfe7087b4d24002ce4af95c4f4b91fbf38a06e8a52275eecca4fb3'
TransportError = columns.TransportError
Unsupported = columns.Unsupported
_OPERATORS = {'div': operator.truediv, 'ge': operator.ge, 'le': operator.le}
_TERMINAL = {'return', 'raise', 'host-raised', 'unsupported'}
_EFFECTS = {'isinstance', 'bind-where', 'getitem', 'getattr', 'prepare-kernel',
            'materialize', 'inspect-float', 'mature'}


def _code_fields(code):
    """Complete immutable code identity, including nested source locations."""
    def constant(value):
        if type(value) is CodeType: return ('code', _code_fields(value))
        if type(value) in (tuple, frozenset):
            items = [constant(item) for item in value]
            if type(value) is frozenset: items.sort(key=repr)
            return (type(value).__name__, tuple(items))
        if type(value) is float: return ('float', struct.pack('>d', value).hex())
        if type(value) in (type(None), bool, int, str, bytes): return (type(value).__name__, value)
        if value is Ellipsis: return ('ellipsis',)
        raise ValueError('unknown package code constant')
    names = ('co_argcount', 'co_posonlyargcount', 'co_kwonlyargcount', 'co_nlocals',
        'co_stacksize', 'co_flags', 'co_code', 'co_names', 'co_varnames', 'co_filename',
        'co_name', 'co_qualname', 'co_firstlineno', 'co_linetable', 'co_exceptiontable',
        'co_freevars', 'co_cellvars')
    return (tuple(getattr(code, name) for name in names), tuple(constant(x) for x in code.co_consts),
            tuple(code.co_lines()), tuple(code.co_positions()))


def _full_plain(function, declaration, context):
    expected = canonical.plain_function(function, declaration, context)
    if _code_fields(function.__code__) != _code_fields(expected):
        raise ValueError('complete visited package code identity differs: ' + declaration)


class _HostRaised(BaseException):
    """Keep actual operation exceptions distinct from our transport errors."""
    def __init__(self, error):
        self.error = error


def _natural(value):
    if type(value) is not str or re.fullmatch('0|[1-9][0-9]*', value) is None:
        raise TransportError('noncanonical procedure integer')
    result = int(value)
    if result >= 1 << 64: raise TransportError('procedure integer exceeds u64')
    return result


def _bits(value):
    if type(value) is not str or re.fullmatch('[0-9a-f]{16}', value) is None:
        raise TransportError('invalid procedure float bits')
    return struct.unpack('>d', bytes.fromhex(value))[0]


def _rows(raw):
    try:
        rows = list(csv.reader(io.StringIO(raw.decode('utf-8'), newline=''),
                               delimiter='\t', strict=True))
    except (UnicodeError, csv.Error) as error:
        raise TransportError('invalid procedure TSV') from error
    if not rows or any(len(row) != 6 for row in rows):
        raise TransportError('procedure record width differs')
    return rows


def _encode(rows):
    stream = io.StringIO(newline='')
    writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
    for row in rows:
        if len(row) != 6: raise TransportError('internal procedure record width differs')
        writer.writerow(row)
    return stream.getvalue().encode('utf-8')


class Registry:
    """Identity-only registration; no repr, equality, iteration or coercion."""
    def __init__(self):
        self.values = {}

    def add(self, value):
        key = len(self.values) + 1
        self.values[key] = value
        return key

    def get(self, key):
        if key == 0 or key not in self.values:
            raise TransportError('foreign procedure object handle')
        return self.values[key]

    def object(self, tag, value):
        if tag == 'handle': return self.get(_natural(value))
        if tag == 'none' and value == '': return None
        if tag == 'str': return value
        if tag == 'float': return _bits(value)
        if tag == 'int' and re.fullmatch('0|-?[1-9][0-9]*', value):
            number = int(value)
            if -(1 << 63) <= number < 1 << 63: return number
        raise TransportError('invalid procedure object encoding')


class Action:
    """Strict wire admission; it never selects a scientific operation."""
    def __init__(self, raw, sequence):
        rows = _rows(raw)
        header = rows.pop(0)
        if (header[0] != 'action' or header[1] != str(sequence)
                or header[5] != PROFILE):
            raise TransportError('procedure action profile/sequence differs')
        self.kind, self.site = header[2], header[4]
        if self.kind not in _EFFECTS | _TERMINAL | {'closed'}:
            raise TransportError('unknown procedure action')
        if self.kind in _EFFECTS:
            self.token = _natural(header[3])
        else:
            if header[3]: raise TransportError('terminal action has an effect token')
            self.token = None
        if self.kind in _EFFECTS | {'unsupported'}:
            match = re.fullmatch('gnomAD_rare_AF@(26[2-9]):(0|[1-9][0-9]*)', self.site)
            if match is None or int(match[2]) > 32:
                raise TransportError('unknown native source operation')
            self.line = int(match[1])
        elif self.site:
            raise TransportError('unexpected terminal source site')
        self.meta, self.objects, self.native = {}, {}, []
        for row in rows:
            if row[0] == 'meta':
                if row[1] in self.meta or any(row[3:]):
                    raise TransportError('duplicate/invalid procedure metadata')
                self.meta[row[1]] = row[2]
            elif row[0] == 'object':
                key = row[1], _natural(row[2])
                if key in self.objects or row[5]:
                    raise TransportError('duplicate/invalid procedure object')
                self.objects[key] = tuple(row[3:5])
            elif row[0] in ('native', 'index', 'cell'):
                self.native.append(row)
            else:
                raise TransportError('unknown procedure payload record')
        if self.kind == 'closed':
            self._expect(set(), set())
            return
        for name in ('native-operations', 'mature-operations'):
            if name not in self.meta: raise TransportError('missing native operation count')
            _natural(self.meta[name])
        counters = {'native-operations', 'mature-operations'}
        self._counters = counters
        if self.kind == 'isinstance':
            self._expect(counters | {'class'}, {('object', 0)})
            if self.meta['class'] not in ('float', 'str'):
                raise TransportError('unrecognized source type test')
        elif self.kind == 'bind-where': self._expect(counters, set())
        elif self.kind == 'getitem': self._expect(counters, {('owner', 0), ('key', 0)})
        elif self.kind == 'getattr':
            self._expect(counters | {'name'}, {('owner', 0)})
            if self.meta['name'] != 'fillna': raise TransportError('unrecognized source method')
        elif self.kind == 'inspect-float': self._expect(counters, {('object', 0)})
        elif self.kind == 'return': self._expect(counters, {('result', 0)})
        elif self.kind == 'raise':
            self._expect(counters | {'class', 'message'}, set())
            if (self.meta['class'], self.meta['message']) != ('ZeroDivisionError', 'float division by zero'):
                raise TransportError('unknown native scalar exception')
        elif self.kind == 'host-raised':
            self._expect(counters | {'exception'}, set())
            _natural(self.meta['exception'])
        elif self.kind == 'unsupported':
            self._expect(counters | {'reason'}, set())
            self._reason()
        else:
            self._call_schema()

    def _reason(self):
        if not 0 < len(self.meta['reason'].encode('utf-8')) <= 4096:
            raise TransportError('invalid explicit mature-operation reason')

    def _expect(self, metadata, objects, *, native=False):
        if (set(self.meta) != metadata or set(self.objects) != objects
                or bool(self.native) != native):
            raise TransportError('procedure action payload differs from named scope')

    def _call_schema(self):
        metadata = self._counters | {'invoke'}
        invoke = self.meta.get('invoke')
        if invoke == 'bound':
            keys = {('callable', 0)}
            count = sum(role == 'argument' for role, _ in self.objects)
            if count not in (1, 3): raise TransportError('invalid captured argument count')
            keys |= {('argument', ordinal) for ordinal in range(count)}
        elif invoke == 'binary':
            metadata |= {'operator'}
            if self.meta.get('operator') not in _OPERATORS:
                raise TransportError('unknown native operator')
            keys = {('left', 0), ('right', 0)}
        else:
            raise TransportError('invalid captured invocation')
        if self.kind == 'mature':
            metadata |= {'reason'}
            self._expect(metadata, keys)
            self._reason()
            return
        metadata |= {'kernel'}
        kernel = self.meta.get('kernel')
        count = {'fillna': 2, 'div': 2, 'ge': 2, 'le': 2, 'where': 3}.get(kernel)
        if count is None: raise TransportError('unknown native table kernel')
        if (kernel in ('fillna', 'where')) != (invoke == 'bound'):
            raise TransportError('kernel invocation class differs')
        if invoke == 'binary' and self.meta['operator'] != kernel:
            raise TransportError('kernel operator differs')
        if invoke == 'bound' and sum(role == 'argument' for role, _ in keys) != (1 if kernel == 'fillna' else 3):
            raise TransportError('kernel captured argument count differs')
        keys |= {('operand', ordinal) for ordinal in range(count)}
        self._expect(metadata, keys, native=self.kind == 'materialize')
        bindings = ([('left', 0, 0), ('right', 0, 1)] if invoke == 'binary' else
                    [('argument', i, i) for i in range(3)] if kernel == 'where' else
                    [('argument', 0, 1)])
        if any(self.objects[role, ordinal] != self.objects['operand', operand]
               for role, ordinal, operand in bindings):
            raise TransportError('native kernel operands differ from the captured invocation')
        if self.kind == 'materialize':
            first = self.native[0]
            if first[0] != 'native' or any(first[4:]) or _natural(first[3]) > columns.MAX_FRAME:
                raise TransportError('invalid materialization header')
            if any(row[0] == 'native' for row in self.native[1:]):
                raise TransportError('duplicate materialization header')

    def object(self, registry, role, ordinal=0):
        try: tag, value = self.objects[role, ordinal]
        except KeyError as error: raise TransportError('missing action object') from error
        return registry.object(tag, value)

    def invocation(self, registry):
        if self.meta['invoke'] == 'binary':
            return _OPERATORS[self.meta['operator']], (self.object(registry, 'left'), self.object(registry, 'right'))
        count = sum(role == 'argument' for role, _ in self.objects)
        return self.object(registry, 'callable'), tuple(self.object(registry, 'argument', i) for i in range(count))

    def operands(self, registry):
        count = sum(role == 'operand' for role, _ in self.objects)
        return tuple(self.object(registry, 'operand', i) for i in range(count))


class ProcedureSession(columns.NativeSession):
    """Reuse the bounded owned-child cleanup; only the named protocol differs."""
    def start(self):
        if self.closed: raise TransportError('native procedure session already closed')
        if self.process is not None: return
        if (hashlib.sha256(columns.guards.regular_bytes(self.binary)).hexdigest() != self.digest
                or not os.access(self.binary, os.X_OK)):
            raise TransportError('native procedure executable identity differs')
        environment = dict(os.environ, POLARS_MAX_THREADS='1', RAYON_NUM_THREADS='1')
        self.process = subprocess.Popen([str(self.binary), 'frequency-procedure-control', 'session'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env=environment, bufsize=0)
        os.set_blocking(self.process.stdin.fileno(), False)
        self.metrics['processes'] += 1
        self.metrics['transport_profile'] = PROFILE

    def exchange(self, command, body=(), *, token=''):
        raw = _encode([('request', self.sequence, command, token, PROFILE, ''), *body])
        if not 0 < len(raw) <= columns.MAX_FRAME:
            raise TransportError('procedure frame exceeds fixed bound')
        self.start()
        if self.process.poll() is not None: raise TransportError('native procedure process unavailable')
        deadline = time.monotonic() + self.timeout
        try:
            self._write(struct.pack('>I', len(raw)) + raw, deadline)
            self.metrics['wire_requests'] += 1
            size = struct.unpack('>I', self._read(4, deadline))[0]
            if not 0 < size <= columns.MAX_FRAME: raise TransportError('invalid procedure response size')
            response = self._read(size, deadline)
        except (OSError, BrokenPipeError) as error:
            raise TransportError('native procedure stream failed') from error
        action = Action(response, self.sequence)
        self.sequence += 1
        self.metrics['procedure_actions_received'] = self.metrics.get('procedure_actions_received', 0) + 1
        if action.kind != 'closed':
            self.metrics['native_math_operations_attempted'] = _natural(action.meta['native-operations'])
            self.metrics['mature_operations_requested'] = _natural(action.meta['mature-operations'])
        return action

    def request(self, operation, args, *, site):
        # Inherited close() owns all drain/wait/exception handling. Its exact
        # quit request is translated, never another scientific procedure.
        if operation != 'quit' or args != [] or site != 'close@0':
            raise TransportError('column request sent to procedure session')
        action = self.exchange('close')
        if action.kind != 'closed': raise TransportError('procedure close acknowledgment differs')
        return 'closed'


class VisitedPackages:
    """Admit only a reached operation against an immutable supplied binding.

    Construction does not inspect the dataframe or validate unvisited package
    methods, cache files, SciPy kernels, or later procedure fields. Source-file
    and immutable-extension checks are demand driven, without a live baseline.
    """
    def __init__(self, path, digest):
        if not packages._sha(digest): raise ValueError('explicit package binding SHA required')
        self.path, self.digest = Path(path), digest
        self.binding, self.context, self.identities = None, {}, {}

    def _load(self):
        import json
        if self.binding is None:
            raw = columns.guards.regular_bytes(self.path)
            if hashlib.sha256(raw).hexdigest() != self.digest:
                raise TransportError('immutable frequency package binding changed')
            document = json.loads(raw)
            if (document.get('schema') != packages.SCHEMA
                    or document.get('reference_verification_sha256') != packages.REFERENCE_VERIFICATION
                    or document.get('frequency_method_declarations') != [list(row) for row in packages.METHOD_SPECS]):
                raise TransportError('unknown frequency package declaration scope')
            self.binding = document
            self.identities[str(self.path)] = packages._identity(self.path)
        for path, identity in self.identities.items():
            if packages._identity(path) != identity:
                raise TransportError('visited package/source identity changed')

    def _role(self, role, module):
        self._load()
        path = str(module.__file__)
        record = self.binding.get('package_files', {}).get(role)
        if (type(record) is not dict or record.get('path') != path
                or not packages._sha(record.get('sha256')) or type(record.get('bytes')) is not int):
            raise TransportError('visited package role/path differs')
        if path not in self.identities:
            raw = columns.guards.regular_bytes(path)
            if len(raw) != record['bytes'] or hashlib.sha256(raw).hexdigest() != record['sha256']:
                raise TransportError('visited package source/extension differs')
            self.identities[path] = packages._identity(path)

    def _contexts(self, names):
        self._load()
        missing = set(names) - set(self.context)
        if missing:
            declared = self.binding.get('source_files', {})
            if any(not packages._sha(declared.get(name)) for name in missing):
                raise TransportError('visited source lacks an immutable declared pin')
            expected = {name: declared[name] for name in missing}
            try: added = canonical.source_context(missing, expected)
            except (ValueError, TypeError, KeyError) as error:
                raise TransportError('visited package source differs') from error
            self.context.update(added)
            for item in added.values():
                self.identities[str(item['module'].__file__)] = packages._identity(item['module'].__file__)

    def _pandas(self):
        from pandas.core import series
        from pandas._libs import internals
        self._role('pandas', pd)
        self._role('pandas_series', series)
        self._role('pandas_internals', internals)
        if pd.__version__ != self.binding['versions']['pandas'] or pd.options.mode.copy_on_write is not False:
            raise TransportError('visited pandas ownership/version differs')
        if pd.Series is not series.Series:
            raise TransportError('visited pandas Series constructor binding differs')

    def transport(self):
        from numpy.core import _multiarray_umath as native
        self._role('numpy', np)
        self._role('numpy_multiarray', native)
        if (np.__version__ != self.binding['versions']['numpy']
                or np.array is not native.array or type(np.array) is not BuiltinFunctionType
                or np.array.__name__ != 'array'
                or not type.__getattribute__(np.ndarray, '__flags__') & (1 << 8)):
            raise TransportError('visited native array transport binding differs')

    def fillna(self, operation, owner):
        if type(owner) is not pd.Series or type(operation) is not MethodType or operation.__self__ is not owner:
            return False
        if 'fillna' in vars(owner): return False
        self._pandas()
        self._contexts({'pandas.core.generic'})
        try:
            _full_plain(operation.__func__, 'pandas.core.generic:NDFrame.fillna', self.context)
        except (ValueError, TypeError, KeyError): return False
        return True

    def binary(self, name, values):
        if any(type(value) is pd.Series for value in values):
            self._pandas()
            self._contexts({'pandas.core.arraylike', 'pandas.core.ops.common'})
            # These exact source declarations cover normal/reflected dispatch;
            # comparing a class slot with itself is not a source proof.
            for method in {'div': ('__truediv__', '__rtruediv__'),
                           'ge': ('__ge__', '__le__'), 'le': ('__le__', '__ge__')}[name]:
                function = canonical.slot('pandas.core.series:Series.' + method)
                declaration = 'pandas.core.arraylike:OpsMixin.' + method
                canonical.wrapper(function, declaration, 'operator', self.context)
                expected = [code for code in canonical.codes(self.context['pandas.core.ops.common']['code'])
                    if code.co_qualname == '_unpack_zerodim_and_defer.<locals>.new_method']
                if len(expected) != 1 or _code_fields(function.__code__) != _code_fields(expected[0]):
                    raise ValueError('complete visited package decorator identity differs')
                closure = dict(zip(function.__code__.co_freevars, function.__closure__ or ()))
                _full_plain(closure['method'].cell_contents, declaration, self.context)

    def where(self, operation):
        from numpy.core import _multiarray_umath as native, multiarray
        self._role('numpy', np)
        self._role('numpy_multiarray', native)
        self._contexts({'numpy.core.multiarray'})
        if np.__version__ != self.binding['versions']['numpy']:
            raise TransportError('visited NumPy version differs')
        dispatcher = native._ArrayFunctionDispatcher
        if not type.__getattribute__(dispatcher, '__flags__') & (1 << 8):
            raise TransportError('mutable NumPy dispatcher class')
        return (type(operation) is dispatcher and operation is np.where
                and operation is multiarray.where and operation.__wrapped__ is native.where
                and type(native.where) is BuiltinFunctionType and native.where.__name__ == 'where')


class HostEffects:
    """Execute a single requested effect; no frequency branch lives here."""
    def __init__(self, registry, package_guard, source_path):
        self.registry, self.packages = registry, package_guard
        self.source_path = str(source_path)
        self._callers = {}
        self.events = []
        self._prepared = None

    def _filename(self, function):
        if function != PROCEDURE: raise TransportError('unknown source warning context')
        return self.source_path

    # Reuse the same small, source-positioned mature-call shim. It only calls
    # the supplied operation; none of the old lowered control flow is invoked.
    _call = columns._Runtime._call

    def _invoke(self, action, operation, args):
        try: return self._call(action.line, PROCEDURE, operation, args, {})
        except BaseException as error: raise _HostRaised(error) from error

    def _object_reply(self, value):
        return [('reply', 'object', self.registry.add(value), '', '', '')]

    def _prepare(self, action):
        values = action.operands(self.registry)
        operation, args = action.invocation(self.registry)
        kernel = action.meta['kernel']
        self._prepared = (tuple(sorted(action.objects.items())), kernel, values)
        try:
            if kernel == 'fillna':
                if not self.packages.fillna(operation, values[0]):
                    raise Unsupported('captured fillna owner/function is outside the native method scope')
                if args[0] is not values[1] and columns._cell(args[0]) != columns._cell(values[1]):
                    raise TransportError('fillna invocation operand differs')
            elif kernel == 'where':
                if not self.packages.where(operation):
                    raise Unsupported('captured NumPy where callback is outside the native scope')
                if any(left is not right for left, right in zip(args, values)):
                    raise TransportError('where invocation operand identity differs')
            else:
                self.packages.binary(kernel, values)
                # ndarray arithmetic can emit NumPy warnings. Its single
                # mature division remains explicit until that effect contract
                # is represented; ordinary Series arithmetic remains native.
                if kernel == 'div' and any(type(value) is np.ndarray for value in values):
                    raise Unsupported('NumPy ndarray division warning effect remains mature')
            self.packages.transport()
            raw = columns._request(0, kernel, values, action.site)
            operands = _rows(raw)[1:]
            return [('reply', 'prepared', len(values), '', '', ''), *operands]
        except Unsupported as error:
            return [('reply', 'not-represented', str(error), '', '', '')]
        except (ValueError, TypeError, KeyError) as error:
            return [('reply', 'not-represented', str(error), '', '', '')]

    def _materialize(self, action):
        signature = tuple(sorted(action.objects.items())), action.meta['kernel']
        if self._prepared is None or signature != self._prepared[:2]:
            raise TransportError('native result does not bind its prepared operands')
        values = self._prepared[2]
        header = action.native[0]
        raw = _encode([('result', '0', header[1], header[2], header[3], columns.PROFILE), *action.native[1:]])
        value = columns._result(raw, 0, values, action.meta['kernel'])
        self._prepared = None
        return self._object_reply(value)

    def execute(self, action):
        self.events.append({'token': action.token, 'site': action.site, 'effect': action.kind})
        try:
            if action.kind == 'isinstance':
                value = action.object(self.registry, 'object')
                result = self._invoke(action, isinstance, (value, float if action.meta['class'] == 'float' else str))
                return [('reply', 'truth', '1' if result else '0', '', '', '')]
            if action.kind == 'bind-where':
                return self._object_reply(columns._reference.np.where)
            if action.kind == 'getitem':
                return self._object_reply(self._invoke(action, operator.getitem,
                    (action.object(self.registry, 'owner'), action.object(self.registry, 'key'))))
            if action.kind == 'getattr':
                return self._object_reply(self._invoke(action, getattr,
                    (action.object(self.registry, 'owner'), action.meta['name'])))
            if action.kind == 'inspect-float':
                value = action.object(self.registry, 'object')
                return [('reply', 'float', struct.pack('>d', value).hex(), '', '', '')] if type(value) is float else [
                    ('reply', 'other-scalar', '', '', '', '')]
            if action.kind == 'prepare-kernel': return self._prepare(action)
            if action.kind == 'materialize': return self._materialize(action)
            if action.kind == 'mature':
                self._prepared = None
                operation, args = action.invocation(self.registry)
                return self._object_reply(self._invoke(action, operation, args))
            raise TransportError('non-effect sent to host dispatcher')
        except _HostRaised as wrapper:
            return [('reply', 'raised', self.registry.add(wrapper.error), '', '', '')]
        except (TransportError, Unsupported):
            raise
        except BaseException as error:
            # Keep the actual object/args/cause; native acknowledges this exact
            # visited failure before the outer caller re-raises it. No replay.
            return [('reply', 'raised', self.registry.add(error), '', '', '')]


def gnomad_rare_af(df, cutoff, *, binary, binary_sha256, package_binding,
                   package_binding_sha256, mode, timeout=120, evidence=None):
    """Run the whole named native procedure; no original-body fallback.

    Explicit mixed mode retains individual unsupported package operations.
    It is recorded separately from fully native math; neither means all mature
    packages were ported. An optional caller-owned evidence dict receives the
    actual effect/lifecycle tape even when the procedure raises.
    """
    if mode not in ('native_math', 'mixed'): raise ValueError('explicit native_math or mixed mode required')
    if hashlib.sha256(columns.guards.regular_bytes(columns._reference.__file__)).hexdigest() != SOURCE_SHA256:
        raise TransportError('frequency source context differs')
    registry = Registry()
    df_handle, cutoff_handle = registry.add(df), registry.add(cutoff)
    host = HostEffects(registry, VisitedPackages(package_binding, package_binding_sha256),
                       columns._reference.__file__)
    session = ProcedureSession(binary, binary_sha256, timeout=timeout)
    terminal, primary = False, None
    try:
        action = session.exchange('start', [('start', PROCEDURE, mode, df_handle, cutoff_handle, '')])
        while action.kind in _EFFECTS:
            body = host.execute(action)
            action = session.exchange('resume', body, token=action.token)
        terminal = True
        if action.kind == 'return': return action.object(registry, 'result')
        if action.kind == 'raise': raise ZeroDivisionError(action.meta['message'])
        if action.kind == 'host-raised':
            error = registry.get(_natural(action.meta['exception']))
            if not isinstance(error, BaseException): raise TransportError('foreign host exception handle')
            raise error
        if action.kind == 'unsupported': raise Unsupported(action.site + ': ' + action.meta['reason'])
        raise TransportError('unexpected procedure terminal action')
    except BaseException as error:
        primary = error
        raise
    finally:
        try:
            session.close(failed=not terminal)
        except BaseException as cleanup_error:
            if primary is None: raise
            session.metrics['cleanup_exception_preserving_primary'] = {
                'type': type(cleanup_error).__qualname__, 'message': str(cleanup_error)[:512]}
        finally:
            if evidence is not None:
                evidence.update(procedure=PROCEDURE, mode=mode, native_control_requested=PROCEDURE,
                    direct_original_body_dispatch_implemented=False,
                    original_entry_observer_required=True, prefix_replay_implemented=False, effects=host.events,
                    terminal_action=action.kind if terminal else None, lifecycle=session.metrics)
