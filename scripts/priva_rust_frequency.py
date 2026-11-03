"""Complete optional frequency table/rule transport.

The original nine scientific bodies and their caller segment determine phase
order. Rust/Polars evaluates admitted table operations; original mature SciPy,
gzip/pickle, object ownership, callbacks and explicit mixed operations remain
hosted. No production dispatch or fallback replays a completed function prefix.
Native qualification has not been performed by adding this source.
"""
import ast
from collections import Counter
import copy
import csv
import hashlib
import io
import math
import operator
import os
from pathlib import Path
import select
import re
import struct
import subprocess
import threading
import time
from types import CodeType, FunctionType

import numpy as np
import pandas as pd

import acmg_pm2_bs1_ba1_frequency as _reference
import priva_rust_frequency_source as guards
import priva_rust_frequency_packages as packages

PROFILE = 'frequency-columns-ef2e701-v1'
MAX_FRAME = 256 * 1024 * 1024
CLOSE_PHASE_SECONDS = 5.0
CLOSE_DIAGNOSTIC_BYTES = 65536
_GLOBAL = guards.RUNTIME
_LOCK = threading.RLock()
_MISSING = object()


class TransportError(RuntimeError):
    """Native transport fault; the scientific original has no such fault."""


class Unsupported(ValueError):
    """An explicit mature-operation fallback, never a transport waiver."""


class _NativeUnsupported(Unsupported):
    """Only a complete, schema-valid native unsupported response creates this."""


def _cell(value):
    if value is None: return 'none', ''
    if value is pd.NA: return 'na', ''
    if type(value) in (bool, np.bool_): return 'bool', str(int(value))
    if type(value) in (int, np.int64, np.int32):
        number = int(value)
        if not -(1 << 53) <= number <= 1 << 53:
            raise Unsupported('integer outside exact native-float domain')
        return 'int', str(number)
    if type(value) in (float, np.float64, np.float32):
        return 'float', struct.pack('>d', float(value)).hex()
    if type(value) in (str, np.str_): return 'str', str(value)
    raise Unsupported('unsupported scalar type')


def _decode_cell(tag, text):
    if tag == 'none' and not text: return None
    if tag == 'na' and not text: return pd.NA
    if tag == 'bool' and text in ('0', '1'): return text == '1'
    if tag == 'int': return int(text)
    if tag == 'float' and len(text) == 16: return struct.unpack('>d', bytes.fromhex(text))[0]
    if tag == 'str': return text
    raise TransportError('invalid native scalar encoding')


def _uncell(tag, text):
    try: return _decode_cell(tag, text)
    except (ValueError, TypeError, OverflowError, UnicodeError) as error:
        raise TransportError('invalid native scalar encoding') from error


def _label(value):
    # Sort-preserving, lossless ordinary homogeneous index labels. Signed
    # integers use an offset, not lexicographic decimal ordering (10 before 2).
    if type(value) is int and -(1 << 63) <= value < (1 << 63):
        return 'i' + f'{value + (1 << 63):020d}'
    if type(value) is str: return 's' + value
    if type(value) is tuple:
        import json
        return 't' + json.dumps([_label(item) for item in value], separators=(',', ':'))
    raise Unsupported('unsupported index label')


def _unlabel(value):
    if value.startswith('i') and len(value) == 21 and value[1:].isdigit():
        number = int(value[1:])
        if number < 1 << 64: return number - (1 << 63)
    if value.startswith('s'): return value[1:]
    if value.startswith('t'): raise TransportError('tuple labels must reuse the exact source axis')
    raise TransportError('invalid native index label')


def _array(value):
    guards.plain_array(value)
    if type(value) is np.ndarray:
        if value.dtype.kind == 'U':
            raise Unsupported('fixed-width Unicode ndarray dtype remains mature')
        if value.dtype.kind == 'f' and any(math.isnan(x) and struct.pack('>d', x).hex() != '7ff8000000000000' for x in value):
            raise Unsupported('noncanonical NaN payload is preserved by mature arithmetic')
        return str(value.dtype), list(value)
    if type(value) is pd.arrays.StringArray:
        if value.dtype.storage != 'python': raise Unsupported('non-Python string storage')
        return 'string', list(value._ndarray)
    if type(value) in (pd.arrays.BooleanArray, pd.arrays.IntegerArray, pd.arrays.FloatingArray):
        if type(value) is pd.arrays.FloatingArray and any(
                not missing and math.isnan(item) for item, missing in zip(value._data, value._mask)):
            raise Unsupported('unmasked nullable NaN is distinct from pandas.NA')
        return str(value.dtype), [pd.NA if missing else item for item, missing in zip(value._data, value._mask)]
    raise Unsupported('unsupported raw storage')


def _value(value):
    if type(value) is pd.Series:
        if len(value._mgr.blocks) != 1: raise Unsupported('non-single-block Series')
        dtype, cells = _array(value._mgr.blocks[0].values)
        labels = [_label(item) for item in guards.axis_values(value._mgr.axes[0])]
        return 'series', dtype, labels, cells
    if type(value) is np.ndarray:
        dtype, cells = _array(value)
        return 'array', dtype, [], cells
    if type(value) in (tuple, list):
        return 'list', 'object', [], list(value)
    _cell(value)
    return 'scalar', 'scalar', [], [value]


def _encode_request(sequence, operation, args, site):
    stream = io.StringIO(newline='')
    writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
    writer.writerow(('request', sequence, operation, len(args), site, PROFILE))
    for ordinal, value in enumerate(args):
        kind, dtype, labels, cells = _value(value)
        writer.writerow(('arg', ordinal, kind, dtype, '', ''))
        for i, label in enumerate(labels):
            writer.writerow(('index', ordinal, i, 'label', label, ''))
        for i, cell in enumerate(cells):
            tag, raw = _cell(cell)
            writer.writerow(('cell', ordinal, i, tag, raw, ''))
    return stream.getvalue().encode('utf-8')


def _request(sequence, operation, args, site):
    try: return _encode_request(sequence, operation, args, site)
    except Unsupported: raise
    except (ValueError, TypeError, OverflowError, UnicodeError) as error:
        raise Unsupported('unrepresented input encoding') from error


def _decode_result(raw, sequence, args, operation):
    rows = list(csv.reader(io.StringIO(raw.decode('utf-8'), newline=''), delimiter='\t', strict=True))
    if not rows or any(len(row) != 6 for row in rows): raise TransportError('invalid native result schema')
    header = rows.pop(0)
    if header[1] != str(sequence) or header[5] != PROFILE: raise TransportError('native result sequence/profile differs')
    if header[0] == 'unsupported':
        if rows or header[3:5] != ['', ''] or not 0 < len(header[2]) <= 4096:
            raise TransportError('invalid unsupported result')
        raise _NativeUnsupported(header[2])
    if header[0] != 'result' or header[2] not in ('scalar', 'series', 'array', 'list', 'pairs'):
        raise TransportError('invalid result kind')
    kind, dtype = header[2:4]
    if dtype not in ('scalar', 'bool', 'boolean', 'int64', 'Int64', 'float64', 'Float64', 'object', 'string'):
        raise TransportError('unknown native result dtype')
    if not header[4].isdigit(): raise TransportError('invalid result denominator')
    size = int(header[4]); labels, cells = [], []
    for row in rows:
        if row[1] != '0' or row[5]: raise TransportError('invalid result reserved field')
        if row[0] == 'index':
            if cells or row[2] != str(len(labels)) or row[3] != 'label': raise TransportError('invalid result index sequence')
            labels.append(row[4])
        elif row[0] == 'cell':
            if row[2] != str(len(cells)): raise TransportError('invalid result cell sequence')
            cells.append(_uncell(row[3], row[4]))
        else: raise TransportError('unknown result record')
    if len(cells) != size or (kind == 'series' and len(labels) != size) or (kind != 'series' and labels):
        raise TransportError('native result denominator mismatch')
    shapes = [('series', len(value._mgr.axes[0])) if type(value) is pd.Series else
              ('array', value.shape[0]) if type(value) is np.ndarray else
              ('list', len(value)) if type(value) in (list, tuple) else ('scalar', 1)
              for value in args]
    prototype_kind = 'series' if any(kind == 'series' for kind, _size in shapes) else next((kind for kind, _size in shapes if kind != 'scalar'), 'scalar')
    expected_kind = ('pairs' if operation == 'group_max_pairs' else
                     'array' if operation in ('where', 'full', 'zeros_bool', 'assign_array') else
                     'scalar' if operation in ('quit', 'sum', 'any') else
                     'list' if operation == 'sorted' else
                     'array' if operation in ('unique', 'zeros', 'assign') else prototype_kind)
    if kind != expected_kind: raise TransportError('native result kind differs from operation')
    if operation == 'quit' and cells != ['closed']:
        raise TransportError('native quit acknowledgment differs')
    if operation in ('zeros', 'zeros_bool', 'full'):
        if size != args[0]: raise TransportError('native zero-array denominator differs')
    elif operation not in ('quit', 'sum', 'any', 'group_max_pairs', 'full', 'zeros_bool'):
        expected_size = next((count for kind, count in shapes if kind != 'scalar'), 1)
        if (operation == 'unique' and size > expected_size
                or operation != 'unique' and size != expected_size):
            raise TransportError('native column denominator differs from operands')
    if kind == 'scalar':
        if size != 1: raise TransportError('nonscalar result')
        return cells[0]
    if kind == 'pairs':
        if dtype != 'object' or len(args) != 2 or len(cells) % 2:
            raise TransportError('invalid grouped pair schema')
        # This validates the wire result's key domain, complete membership and
        # ordering. It does not recompute any grouped scientific value.
        expected = []
        try: source_keys = _value(args[1])[3]
        except Unsupported as error:
            raise TransportError('grouping operands changed before native reply') from error
        for key in source_keys:
            if key is None or key is pd.NA or type(key) in (float, np.float64) and math.isnan(key):
                continue
            try: normalized = _uncell(*_cell(key))
            except Unsupported as error:
                raise TransportError('unapproved source group key in native reply') from error
            if type(normalized) not in (str, int):
                raise TransportError('native grouping reply for an unapproved key type')
            expected.append(normalized)
        if len({type(key) for key in expected}) > 1:
            raise TransportError('native grouping reply for mixed key types')
        expected = sorted(set(expected))
        keys, values = cells[::2], cells[1::2]
        if (any(type(key) not in (str, int) for key in keys)
                or len({(type(key), key) for key in keys}) != len(keys)
                or len(keys) != len(expected)
                or [(type(key), key) for key in keys] != [(type(key), key) for key in expected]
                or any(type(value) is not float for value in values)):
            raise TransportError('grouped result keys/count/order/value types differ')
        return dict(zip(keys, values))
    if kind == 'list': return cells
    try:
        if dtype in ('boolean', 'Int64', 'Float64', 'string'):
            values = pd.array(cells, dtype=pd.StringDtype(storage='python') if dtype == 'string' else dtype)
        else:
            values = np.array(cells, dtype=dtype)
    except (ValueError, TypeError) as error:
        raise TransportError('invalid native dtype/cell combination') from error
    if kind == 'array': return values
    prototypes = [value for value in args if type(value) is pd.Series]
    if (not prototypes or size != len(prototypes[0])
            or labels != [_label(item) for item in guards.axis_values(prototypes[0]._mgr.axes[0])]):
        raise TransportError('native result Series axis is foreign/incomplete')
    prototype = prototypes[0] if prototypes else None
    axis = None
    for value in prototypes:
        if [_label(item) for item in guards.axis_values(value._mgr.axes[0])] == labels:
            axis = value._mgr.axes[0]
            break
    if axis is None:
        axis = pd.Index([_unlabel(label) for label in labels])
    name = prototype.name if prototype is not None else None
    if len(prototypes) > 1 and any(value.name != name for value in prototypes[1:]): name = None
    if operation == 'max_columns': name = None
    result = pd.Series(values, index=axis, name=name)
    if prototype is not None and operation not in ('numeric', 'max_columns'):
        result.attrs = copy.deepcopy(vars(prototype).get('_attrs', {}))
    return result


def _result(raw, sequence, args, operation):
    try: return _decode_result(raw, sequence, args, operation)
    except (_NativeUnsupported, TransportError): raise
    except (ValueError, TypeError, OverflowError, UnicodeError, csv.Error,
            KeyError, IndexError, MemoryError) as error:
        raise TransportError('corrupt native response') from error


class NativeSession:
    """One child per complete procedure; no row cache or speculative reads."""
    def __init__(self, binary, digest, *, timeout=120):
        self.binary = Path(binary)
        self.digest = digest
        self.timeout = timeout
        self.process = None
        self.closed = False
        self.sequence = 0
        self.metrics = dict(processes=0, wire_requests=0, native_operations=Counter(),
                            mature_fallbacks=Counter(), source_operations=Counter(),
                            transport_profile=PROFILE, cache_hits=0)

    def start(self):
        if self.closed: raise TransportError('native session already closed')
        if self.process is not None: return
        if (hashlib.sha256(guards.regular_bytes(self.binary)).hexdigest() != self.digest
                or not os.access(self.binary, os.X_OK)):
            raise TransportError('native executable identity differs')
        environment = dict(os.environ, POLARS_MAX_THREADS='1', RAYON_NUM_THREADS='1')
        # No Python process, environment creation or shell interpretation.
        self.process = subprocess.Popen([str(self.binary), 'frequency-criteria', 'session'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env=environment, bufsize=0)
        os.set_blocking(self.process.stdin.fileno(), False)
        self.metrics['processes'] += 1

    def _read(self, count, deadline):
        result = bytearray()
        while len(result) < count:
            remaining = deadline - time.monotonic()
            if remaining <= 0: raise TransportError('native response timeout')
            ready, _, _ = select.select([self.process.stdout], [], [], remaining)
            if not ready: raise TransportError('native response timeout')
            chunk = os.read(self.process.stdout.fileno(), count - len(result))
            if not chunk: raise TransportError('native response truncated')
            result.extend(chunk)
        return bytes(result)

    def _write(self, frame, deadline):
        position = 0
        while position < len(frame):
            remaining = deadline - time.monotonic()
            if remaining <= 0: raise TransportError('native request write timeout')
            if self.process.poll() is not None: raise TransportError('native process unavailable during write')
            _, ready, _ = select.select([], [self.process.stdin], [], remaining)
            if not ready: raise TransportError('native request write timeout')
            try:
                count = os.write(self.process.stdin.fileno(), memoryview(frame)[position:position + 65536])
            except BlockingIOError:
                continue
            if not 0 < count <= min(65536, len(frame) - position):
                raise TransportError('native request write made invalid progress')
            position += count

    def request(self, operation, args, *, site):
        body = _request(self.sequence, operation, args, site)  # admission before child IO
        if not 0 < len(body) <= MAX_FRAME: raise Unsupported('column frame exceeds fixed bound')
        self.start()
        if self.process.poll() is not None: raise TransportError('native process unavailable')
        deadline = time.monotonic() + self.timeout
        try:
            self._write(struct.pack('>I', len(body)) + body, deadline)
            self.metrics['wire_requests'] += 1
            size = struct.unpack('>I', self._read(4, deadline))[0]
            if not 0 < size <= MAX_FRAME: raise TransportError('invalid native response length')
            raw = self._read(size, deadline)
        except (OSError, BrokenPipeError) as error:
            raise TransportError('native stream failed') from error
        sequence = self.sequence
        self.sequence += 1
        value = _result(raw, sequence, args, operation)
        if operation != 'quit': self.metrics['native_operations'][operation] += 1
        return value

    def _finish_child(self, process):
        """Boundedly drain and reap this exact owned Popen, on every close path.

        No process-group lookup, external PID lookup or scientific retry is
        involved. A child surviving bounded SIGKILL/wait remains an explicit
        unreaped transport failure, never a successful cleanup claim.
        """
        errors, signals, pipes = [], [], {}
        retained = {name: bytearray() for name in ('stdout', 'stderr')}
        counts = {name: 0 for name in retained}

        def error(phase, exception):
            errors.append(dict(phase=phase, type=type(exception).__qualname__,
                               message=str(exception)[:512]))

        def exited():
            try: return process.poll() is not None
            except BaseException as exception:
                error('poll', exception)
                return process.returncode is not None

        try: process.stdin.close()
        except BaseException as exception: error('stdin-close', exception)
        for name in retained:
            stream = getattr(process, name)
            try:
                descriptor = stream.fileno()
                os.set_blocking(descriptor, False)
                pipes[descriptor] = name
            except BaseException as exception: error(name + '-admission', exception)

        def drain(deadline):
            while pipes or not exited():
                remaining = deadline - time.monotonic()
                if remaining <= 0: return
                try:
                    ready, _, _ = select.select(list(pipes), [], [], min(.05, remaining))
                except BaseException as exception:
                    error('pipe-select', exception)
                    return
                for descriptor in ready:
                    name = pipes[descriptor]
                    try: block = os.read(descriptor, 65536)
                    except BlockingIOError: continue
                    except BaseException as exception:
                        error(name + '-read', exception)
                        pipes.pop(descriptor)
                        continue
                    if not block:
                        pipes.pop(descriptor)
                        continue
                    counts[name] += len(block)
                    retained[name].extend(block[:max(0, CLOSE_DIAGNOSTIC_BYTES - len(retained[name]))])

        # Drain while the child runs: waiting before reading stderr can block
        # an otherwise successful child on a full diagnostic pipe.
        for phase in ('grace', 'terminate', 'kill'):
            if phase != 'grace' and not exited():
                signals.append(phase)
                try: getattr(process, phase)()
                except BaseException as exception: error(phase, exception)
            drain(time.monotonic() + CLOSE_PHASE_SECONDS)
            if exited() and not pipes: break
        try: process.wait(timeout=CLOSE_PHASE_SECONDS)
        except BaseException as exception: error('final-wait', exception)
        for name in retained:
            try: getattr(process, name).close()
            except BaseException as exception: error(name + '-close', exception)
        result = dict(pid=process.pid, returncode=process.returncode,
            reaped=process.returncode is not None, termination_attempts=signals,
            errors=errors, diagnostic_limit_per_stream=CLOSE_DIAGNOSTIC_BYTES,
            remaining_open_pipes=len(pipes),
            streams={name: dict(bytes=counts[name], retained_hex=bytes(retained[name]).hex(),
                               truncated=counts[name] > len(retained[name])) for name in retained})
        self.metrics['native_close'] = result
        return result

    def close(self, *, failed=False):
        if self.process is None or self.closed: return
        process = self.process
        primary = None
        try:
            if not failed: self.request('quit', [], site='close@0')
        except BaseException as exception:
            primary = exception
        finally:
            try:
                cleanup = self._finish_child(process)
                self.closed = cleanup['reaped']
            except BaseException as cleanup_error:
                # Preserve the first fault even if unexpected cleanup code
                # fails. The final bounded stop/wait still targets only the
                # Popen created by this session, never a discovered PID.
                failures = [dict(phase='cleanup', type=type(cleanup_error).__qualname__,
                                 message=str(cleanup_error)[:512])]
                for name, operation in (('kill', process.kill),
                                        ('wait', lambda: process.wait(timeout=CLOSE_PHASE_SECONDS))):
                    try: operation()
                    except BaseException as error:
                        failures.append(dict(phase=name, type=type(error).__qualname__, message=str(error)[:512]))
                snapshots = {}
                for name in ('stdout', 'stderr'):
                    try:
                        descriptor = getattr(process, name).fileno()
                        os.set_blocking(descriptor, False)
                        block = os.read(descriptor, CLOSE_DIAGNOSTIC_BYTES + 1)
                        snapshots[name] = dict(retained_hex=block[:CLOSE_DIAGNOSTIC_BYTES].hex(),
                            truncated=len(block) > CLOSE_DIAGNOSTIC_BYTES, complete_stream_not_claimed=True)
                    except BaseException as error:
                        failures.append(dict(phase=name + '-snapshot', type=type(error).__qualname__, message=str(error)[:512]))
                for stream in (process.stdin, process.stdout, process.stderr):
                    try: stream.close()
                    except BaseException: pass
                self.closed = process.returncode is not None
                self.metrics['native_close_emergency'] = dict(pid=process.pid,
                    returncode=process.returncode, reaped=self.closed, errors=failures, snapshots=snapshots)
                if primary is not None: raise primary from cleanup_error
                raise TransportError('unexpected native cleanup failure; exact owned child stop/wait attempted') from cleanup_error
        if primary is not None: raise primary
        if not cleanup['reaped'] or cleanup['errors'] or cleanup['remaining_open_pipes']:
            raise TransportError('native child cleanup failed; bounded diagnostics retained')
        if failed: return
        if (cleanup['returncode'] or cleanup['termination_attempts']
                or any(row['bytes'] for row in cleanup['streams'].values())):
            raise TransportError('native session closed with failure/extra output; bounded diagnostics retained')
        if hashlib.sha256(guards.regular_bytes(self.binary)).hexdigest() != self.digest:
            raise TransportError('native executable changed during procedure')


_OPERATORS = dict(add=operator.add, sub=operator.sub, mul=operator.mul,
    div=operator.truediv, **{'and': operator.and_, 'or': operator.or_},
    eq=operator.eq, ne=operator.ne, lt=operator.lt, le=operator.le,
    gt=operator.gt, ge=operator.ge)


def _vector(value):
    return type(value) in (pd.Series, np.ndarray)


class _Grouped:
    """Owned original groupby construction with a native aggregate phase."""
    def __init__(self, runtime, line, function, actual, values, keys):
        self.runtime, self.line, self.function = runtime, line, function
        self.actual, self.values, self.keys = actual, values, keys

    def max(self):
        runtime = self.runtime
        actual = self.actual.max
        if not runtime.packages.method_is_canonical(actual, 'max'):
            return runtime._mature(self.line, self.function, 'group_max_pairs',
                actual, (), {}, 'group maximum source declaration differs')
        try: runtime.packages.declared_method('pandas.core.series:Series.to_dict')
        except (ValueError, TypeError, KeyError) as error:
            return runtime._mature(self.line, self.function, 'group_max_pairs',
                actual, (), {}, str(error))
        try:
            result = runtime.session.request('group_max_pairs', (self.values, self.keys),
                site=runtime._site(self.line, self.function))
        except Unsupported as error:
            return runtime._mature(self.line, self.function, 'group_max_pairs',
                actual, (), {}, str(error))
        return _GroupedDictionary(result, runtime)


class _GroupedDictionary:
    def __init__(self, value, runtime): self.value, self.runtime = value, runtime
    def to_dict(self):
        # A drift after max must not replay the already completed aggregate.
        self.runtime.packages.declared_method('pandas.core.series:Series.to_dict')
        return self.value


class _Runtime:
    """Execute one source operation, preserving the already-evaluated inputs.

    A native unsupported response permits exactly that mature operation once
    in explicit mixed mode. Transport corruption never invokes a fallback.
    Mature SciPy/codecs/ownership are separately identified dependencies; they
    are not hidden inside a claim of whole-native execution.
    """
    def __init__(self, session, package_guard, caller_path, *, mode):
        if mode not in ('native_math', 'mixed'):
            raise ValueError('explicit native_math or mixed execution mode required')
        self.session = session
        self.packages = package_guard
        self.caller_path = str(caller_path)
        self.mode = mode
        self.metrics = session.metrics
        self.metrics.update(mode=mode, whole_native=False, host_control=Counter(),
            mature_dependencies=Counter(), source_fallbacks=[], numerical_calls=[],
            ownership='original pandas/NumPy objects, indexing and writes')
        self._callers = {}

    def _site(self, line, function):
        site = function + '@' + str(line)
        self.metrics['source_operations'][site] += 1
        return site

    def _filename(self, function):
        if function == 'ACMG_criteria_assign': return self.caller_path
        if function in ('_variant_condition_masks', 'has_token'):
            return _reference._variant_condition_masks.__code__.co_filename
        return _reference.__file__

    def _call(self, line, function, operation, args, kwargs):
        # Calling the actual mature operation at the original filename/line
        # preserves its warning location. This is not an original-function
        # retry and does not re-evaluate any operand or earlier statement.
        key = function, line
        code = self._callers.get(key)
        if code is None:
            tree = ast.parse('def operation_at_source(_operation, _args, _kwargs):\n return _operation(*_args, **_kwargs)\n')
            ast.increment_lineno(tree, line - 2)
            tree.body[0].lineno = max(1, line - 1)
            code = next(item for item in compile(tree, self._filename(function), 'exec').co_consts
                        if type(item) is CodeType)
            self._callers[key] = code
        return FunctionType(code, {})(operation, args, kwargs)

    def _mature(self, line, function, name, operation, args, kwargs, reason):
        site = function + '@' + str(line)
        self.metrics['source_fallbacks'].append(dict(site=site, operation=name, reason=reason))
        self.metrics['mature_fallbacks'][name] += 1
        if self.mode == 'native_math':
            raise Unsupported(site + ': ' + name + ': ' + reason)
        return self._call(line, function, operation, args, kwargs)

    def _native(self, line, function, name, values, operation, args, kwargs, validate=None):
        site = self._site(line, function)
        try:
            result = self.session.request(name, values, site=site)
            if validate is not None: validate(result)
        except Unsupported as error:
            return self._mature(line, function, name, operation, args, kwargs, str(error))
        return result

    def binary(self, line, function, name, left, right):
        operation = _OPERATORS[name]
        if not _vector(left) and not _vector(right):
            self.metrics['host_control'][function + '@' + str(line)] += 1
            return self._call(line, function, operation, (left, right), {})
        try: self.packages.operator(name, left, right)
        except (ValueError, TypeError, KeyError) as error:
            return self._mature(line, function, name, operation, (left, right), {}, str(error))
        # NumPy array floating-warning behavior is observable. Keep warning
        # cases mature instead of silently replacing the original ufunc.
        if name in ('add', 'sub', 'mul', 'div') and any(type(x) is np.ndarray for x in (left, right)):
            try:
                cells = [cell for x in (left, right) for cell in _value(x)[3]]
                if any(type(x) in (float, np.float64) and not math.isfinite(x) for x in cells):
                    raise Unsupported('nonfinite NumPy arithmetic warning frontier')
                if name == 'div' and any(type(x) in (int, float, np.int64, np.float64) and x == 0 for x in _value(right)[3]):
                    raise Unsupported('NumPy division warning frontier')
            except (Unsupported, ValueError) as error:
                return self._mature(line, function, name, operation, (left, right), {}, str(error))
        def validate(result):
            # Inspect only the pure native result, before any mature operation.
            # Never retry a mature result, including a nonfinite one.
            if name in ('add', 'sub', 'mul', 'div') and type(result) is np.ndarray and result.dtype == np.dtype('float64'):
                if any(not math.isfinite(x) for x in result):
                    raise Unsupported('NumPy result warning frontier')
        return self._native(line, function, name, (left, right), operation, (left, right), {}, validate)

    def invert(self, line, function, value):
        if not _vector(value):
            self.metrics['host_control'][function + '@' + str(line)] += 1
            return self._call(line, function, operator.invert, (value,), {})
        try: self.packages.operator('invert', value)
        except (ValueError, TypeError, KeyError) as error:
            return self._mature(line, function, 'invert', operator.invert, (value,), {}, str(error))
        return self._native(line, function, 'invert', (value,), operator.invert, (value,), {})

    def numeric(self, line, function, operation, *args, **kwargs):
        try: self.packages.numeric(operation)
        except (ValueError, TypeError, KeyError) as error:
            return self._mature(line, function, 'numeric', operation, args, kwargs, str(error))
        if len(args) != 1 or kwargs != {'errors': 'coerce'} or not _vector(args[0]):
            return self._mature(line, function, 'numeric', operation, args, kwargs, 'noncanonical numeric signature')
        return self._native(line, function, 'numeric', args, operation, args, kwargs)

    def numpy(self, line, function, name, operation, *args, **kwargs):
        try: self.packages.array_function(name, operation)
        except (ValueError, TypeError, KeyError) as error:
            return self._mature(line, function, name, operation, args, kwargs, str(error))
        route = name
        values = args
        if name == 'zeros' and len(args) == 1 and set(kwargs) == {'dtype'}:
            if kwargs['dtype'] is bool: route = 'zeros_bool'
            elif kwargs['dtype'] is int: route = 'zeros'
            else: route = None
        elif name == 'full' and len(args) == 2 and not kwargs: pass
        elif name in ('where', 'minimum', 'maximum', 'round', 'isnan', 'logical_not') and not kwargs:
            if name == 'logical_not': route = 'not'
            if name == 'round' and len(args) != 1: route = None
        else: route = None
        if route is None:
            return self._mature(line, function, name, operation, args, kwargs, 'unrepresented NumPy signature')
        return self._native(line, function, route, values, operation, args, kwargs)

    def kernel(self, line, function, name, operation, *args, **kwargs):
        # The exact original bound method, original argument objects, keyword
        # order and eager call order survive. No CDF/SF approximation is used.
        self.packages.binomial(name, operation)
        self.metrics['numerical_calls'].append(dict(site=function + '@' + str(line),
            method=name, positional_count=len(args), keyword_order=list(kwargs)))
        self.metrics['mature_dependencies']['scipy.stats.binom.' + name] += 1
        value = self._call(line, function, operation, args, kwargs)
        self.packages.check()
        return value

    def after_callback(self, line, function, value):
        # The original logger call was evaluated before this method, including
        # its arguments and callbacks. Check the next phase's binding afterward.
        self.metrics['mature_dependencies']['logging'] += 1
        self.packages.check(complete=True)
        guards.reference_guard(_reference, self.caller_path)
        return value

    def filter_cache(self, line, function, mapping):
        def original_filter(): return {key: value for key, value in mapping.items() if value is not None}
        if type(mapping) is not dict or any(type(key) not in (str, int) for key in mapping):
            return self._mature(line, function, 'filter_none', original_filter, (), {}, 'custom mapping iteration/key hashing')
        # Nested cache objects remain owned by the mature decoder. Only the
        # literal None-vs-present predicate crosses this native interface.
        keys = list(mapping)
        predicates = [None if value is None else True for value in mapping.values()]
        mask = self._native(line, function, 'filter_none', (predicates,),
                            original_filter, (), {})
        if type(mask) is dict: return mask  # exactly one explicit mixed result
        if type(mask) is not list or len(mask) != len(keys) or any(type(x) is not bool for x in mask):
            raise TransportError('cache filter mask differs')
        return {key: mapping[key] for key, keep in zip(keys, mask) if keep}

    def _exon_operands(self, frame, callback, offset):
        if type(frame) is not pd.DataFrame or type(callback) is not FunctionType:
            raise Unsupported('custom exon apply owner/callback')
        if callback.__code__.co_filename != str(_reference.__file__):
            raise Unsupported('foreign exon callback')
        closure = dict(zip(callback.__code__.co_freevars, callback.__closure__ or ()))
        if set(closure) != {'clinvar_patho_exon_af_stat_dict'}:
            raise Unsupported('exon callback closure differs')
        cache = closure['clinvar_patho_exon_af_stat_dict'].cell_contents
        if type(cache) is not dict:
            raise Unsupported('custom exon cache mapping')
        features, exons, values = [], [], []
        for feature, entry in cache.items():
            if type(feature) is not str or type(entry) is not dict:
                raise Unsupported('malformed exon cache remains mature')
            for exon, summary in entry.items():
                if (type(exon) is not str or type(summary) not in (tuple, list)
                        or len(summary) <= offset
                        or type(summary[offset]) not in (int, float)):
                    raise Unsupported('malformed exon summary remains mature')
                features.append(feature); exons.append(exon); values.append(summary[offset])
        labels = guards.axis_values(frame._mgr.axes[0])
        if labels.count('Feature') != 1 or labels.count('EXON') != 1 or len(frame) == 0:
            raise Unsupported('missing/duplicate/empty exon frame remains mature')
        columns = guards.raw_columns(frame)
        a, b = columns[labels.index('Feature')], columns[labels.index('EXON')]
        # Row apply's mixed numeric coercion and extension boxing have their
        # own mature semantics. The native two-key lookup admits plain text.
        for column in (a, b):
            guards.plain_array(column)
            if any(type(cell) is not str for cell in column):
                raise Unsupported('exon row key coercion remains mature')
        index = frame._mgr.axes[1]
        return (pd.Series(a, index=index), pd.Series(b, index=index), features, exons, values)

    def exon_apply(self, line, function, operation, frame, args, kwargs):
        try:
            values = self._exon_operands(frame, args[0], 0 if line == 362 else 2)
        except (Unsupported, ValueError, TypeError) as error:
            return self._mature(line, function, 'exon_map', operation, args, kwargs, str(error))
        # The original ordinary row lookup makes this column-index engine
        # observable. Reproduce its two ordered lookups, after all no-effect
        # admission and before the native operation; never synthesize an engine.
        frame._mgr.axes[0].get_loc('Feature')
        frame._mgr.axes[0].get_loc('EXON')
        result = self._native(line, function, 'exon_map', values, operation, args, kwargs)
        if type(result) is not pd.Series:
            raise TransportError('exon result ownership differs')
        result.name = None
        result.attrs = copy.deepcopy(vars(frame).get('_attrs', {}))
        return result

    def _map(self, line, function, operation, owner, args, kwargs):
        mapping = args[0]
        default = 'cache_map_nan'
        if type(mapping) is FunctionType:
            # Only the exact source lambdas below have a native map contract.
            allowed = {('PM2_criteria', 356): 'clinvar_patho_af_dict',
                       ('PM2_criteria', 359): 'gene_max_patho_hom_freq_dict',
                       ('BS1_criteria', 598): 'hom_dict'}
            name = allowed.get((function, line))
            if name is None or mapping.__code__.co_filename != str(_reference.__file__):
                raise Unsupported('unrepresented/custom mapping callback')
            closure = dict(zip(mapping.__code__.co_freevars, mapping.__closure__ or ()))
            if set(closure) != {name}: raise Unsupported('mapping closure differs')
            mapping = closure[name].cell_contents
            default = 'cache_map_zero'
            if type(mapping) is not dict: raise Unsupported('custom cache mapping')
            if name == 'clinvar_patho_af_dict':
                # Validate without calling custom .get or consuming a partial
                # scientific prefix. Malformed entries stay with the original
                # lambda, which fails at its original first offending row.
                if any(type(value) is not dict for value in mapping.values()):
                    raise Unsupported('malformed gene AF cache remains mature')
                mapping = {key: value.get('af', 0) for key, value in mapping.items()}
        if type(mapping) is not dict or not (not kwargs or kwargs == {'na_action': 'ignore'}):
            raise Unsupported('custom map/na_action semantics')
        if any(type(key) not in (str, int) or type(key) is bool for key in mapping):
            raise Unsupported('mapping key equality outside native domain')
        if kwargs.get('na_action') == 'ignore' and any(cell is pd.NA or cell is None for cell in _value(owner)[3]):
            raise Unsupported('map NA sentinel identity remains mature')
        return default, (owner, list(mapping), list(mapping.values()))

    def method(self, line, function, operation, name, *args, **kwargs):
        owner = getattr(operation, '__self__', None)
        if type(owner) in (_Grouped, _GroupedDictionary):
            if kwargs or args or name not in ('max', 'to_dict'):
                raise TransportError('private grouped phase signature differs')
            return operation()
        if not self.packages.method_is_canonical(operation, name):
            return self._mature(line, function, name, operation, args, kwargs, 'noncanonical method/owner')
        if (name == 'apply' and function == 'PM2_criteria' and line in (362, 363)
                and len(args) == 1 and kwargs == {'axis': 1}):
            return self.exon_apply(line, function, operation, owner, args, kwargs)
        if (name == 'groupby' and function == '_gene_max_patho_hom_freq' and line == 65
                and type(owner) is pd.Series and len(args) == 1
                and type(args[0]) is pd.Series and not kwargs):
            actual = self._call(line, function, operation, args, kwargs)
            self.metrics['mature_dependencies']['groupby_ownership'] += 1
            return _Grouped(self, line, function, actual, owner, args[0])
        from pandas.core.strings.accessor import StringMethods
        if type(owner) is StringMethods:
            parent = vars(owner).get('_parent')
            if type(parent) is not pd.Series:
                return self._mature(line, function, name, operation, args, kwargs, 'custom string accessor')
            owner = parent
        if not _vector(owner):
            return self._mature(line, function, name, operation, args, kwargs,
                                'grouping/row callback or custom ownership remains mature')
        route, values = None, None
        if name in ('isna', 'isnull', 'notna', 'sum', 'any', 'lower', 'upper') and not args and not kwargs:
            route = 'isna' if name == 'isnull' else name
            values = (owner,)
        elif name == 'fillna' and len(args) == 1 and not kwargs:
            route, values = 'fillna', (owner, args[0])
        elif name == 'astype' and len(args) == 1 and not kwargs:
            choices = {int: 'cast_int64', bool: 'cast_bool', str: 'as_str',
                       'boolean': 'cast_boolean', 'bool': 'cast_bool', 'int64': 'cast_int64'}
            try: route = choices.get(args[0])
            except TypeError: route = None
            if type(owner) is np.ndarray and args[0] is str: route = None
            values = (owner,)
        elif name == 'clip' and not kwargs and len(args) == 2:
            route, values = 'clip', (owner, *args)
        elif name == 'clip' and not args and set(kwargs) <= {'lower', 'upper'}:
            route, values = 'clip', (owner, kwargs.get('lower'), kwargs.get('upper'))
        elif name == 'mask' and not kwargs and len(args) in (1, 2):
            route, values = ('mask_missing' if len(args) == 1 else 'mask_value'), (owner, *args)
        elif name == 'replace' and not kwargs and len(args) == 2 and type(args[0]) is int and args[0] == 0 and type(args[1]) is float and math.isnan(args[1]):
            route, values = 'replace_zero_nan', (owner,)
        elif name in ('eq', 'ge', 'removeprefix') and len(args) == 1 and not kwargs:
            route, values = name, (owner, args[0])
        elif name == 'isin' and len(args) == 1 and not kwargs and type(args[0]) in (set, list, tuple):
            if all(type(item) is str for item in args[0]):
                route, values = 'isin', (owner, list(args[0]))
        elif name == 'contains' and len(args) == 1 and type(args[0]) is str:
            if kwargs == {'regex': False}:
                route, values = 'contains_literal', (owner, args[0])
            elif kwargs == {'regex': True}:
                match = re.fullmatch(r'\(\?:\^\|;\)([a-z_]+)\(\?:;\|\$\)', args[0])
                if match: route, values = 'condition_token', (owner, match.group(1))
        elif name == 'map' and len(args) == 1:
            try: route, values = self._map(line, function, operation, owner, args, kwargs)
            except (Unsupported, ValueError) as error:
                return self._mature(line, function, name, operation, args, kwargs, str(error))
        if route is None:
            return self._mature(line, function, name, operation, args, kwargs, 'unrepresented method/dtype signature')
        return self._native(line, function, route, values, operation, args, kwargs)

    def assign(self, line, function, value, array, mask):
        # RHS was evaluated FIRST by the lowered AST. Only an ordinary writable
        # ndarray uses a native assignment. Publish into that same owned array;
        # returning a replacement would break aliases and the caller copy rule.
        def original_assignment():
            array[mask] = value
        if type(array) is not np.ndarray or not array.flags.writeable:
            self.metrics['mature_dependencies']['owned_assignment'] += 1
            return self._call(line, function, original_assignment, (), {})
        try:
            result = self.session.request('assign_array', (array, mask, value),
                                          site=self._site(line, function))
        except Unsupported as error:
            return self._mature(line, function, 'assign_array', original_assignment, (), {}, str(error))
        if type(result) is not np.ndarray or result.dtype != array.dtype or result.shape != array.shape:
            raise TransportError('native assignment publication shape/dtype differs')
        array[...] = result
        return None


class FrequencyBridge:
    """Explicit complete-unit API; no caller or default activation is installed.

    Construct only with reviewed executable and package bindings. ``run``
    executes one selected original procedure's private lowered body. In mixed
    mode, every host table operation is reported by source site; no original
    procedure is rerun after any native or original prefix has executed.
    """
    def __init__(self, *, binary, binary_sha256, package_binding, package_sha256,
                 caller_path=None, mode='native_math', timeout=120):
        if not packages._sha(binary_sha256):
            raise ValueError('pinned native executable required')
        self.caller_path = Path(caller_path or Path(_reference.__file__).with_name('acmg_criteria_assign.py'))
        guards.reference_guard(_reference, self.caller_path)
        self.package_guard = packages.Guard(package_binding, package_sha256, _reference)
        self.session = NativeSession(binary, binary_sha256, timeout=timeout)
        self.runtime = _Runtime(self.session, self.package_guard, self.caller_path, mode=mode)
        self.functions, self.source_sites = guards.lowered(_reference, self.caller_path, self.runtime)
        self.used = False

    @property
    def metrics(self): return self.runtime.metrics

    def run(self, procedure, *args, **kwargs):
        if procedure not in (*guards.PROCEDURES, 'ACMG_criteria_assign'):
            raise ValueError('unknown complete frequency procedure')
        if self.used: raise ValueError('one procedure per owned native session')
        self.used = True
        with _LOCK:
            try:
                value = self.functions[procedure](*args, **kwargs)
                self.package_guard.check(complete=True)
                guards.reference_guard(_reference, self.caller_path)
            except BaseException:
                try: self.session.close(failed=True)
                except BaseException as error:
                    self.metrics['close_failure'] = dict(type=type(error).__qualname__, message=str(error))
                raise
            self.session.close(failed=False)
            return value


def run(procedure, *args, binary=None, binary_sha256=None, package_binding=None,
        package_sha256=None, caller_path=None, mode='native_math', metrics=None, **kwargs):
    """Call an original function by default; opt in with every native binding.

    The captured caller segment is a separate qualification-only procedure,
    never an alternative to the full ACMG caller.
    """
    if binary is None:
        if procedure not in guards.PROCEDURES:
            raise ValueError('caller segment requires explicit qualification bridge')
        return getattr(_reference, procedure)(*args, **kwargs)
    bridge = FrequencyBridge(binary=binary, binary_sha256=binary_sha256,
        package_binding=package_binding, package_sha256=package_sha256,
        caller_path=caller_path, mode=mode)
    try: return bridge.run(procedure, *args, **kwargs)
    finally:
        if metrics is not None: metrics.update(bridge.metrics)


