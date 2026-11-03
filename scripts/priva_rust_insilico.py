"""Optional complete PP3/BP4/BP7 native columns; Python stays the default.

The unchanged source AST is the orchestration authority. Its statement,
DataFrame/Series getter, accessor-construction, logger and branch order is
retained. Native column procedures replace table calculations, including
coercion and reductions, not merely the final criterion masks. Mature boolean
Series selection and raw ndarray/Series ownership remain hosted explicitly.

Unknown whole-input profiles use the original before observable getters. A
later unsupported mature operation (for example, a logging callback introduced
an arbitrary Python object) executes that original operation ONCE at its source
frontier and is reported as mixed execution. Transport failures never trigger
replay/fallback. No source, config, model, threshold or default is changed.
"""
import ast
from collections import Counter
import copy
import csv
import hashlib
import io
import operator
import os
from pathlib import Path
import select
import struct
import subprocess
import threading
import time
from types import CodeType, FunctionType

import numpy as np
import pandas as pd

import acmg_pp3_bp4_bp7_insilico as _reference
import priva_rust_insilico_source as guards

PROFILE = 'insilico-columns-v1'
MAX_FRAME = 256 * 1024 * 1024
_GLOBAL = '_priva_insilico_column_runtime'
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
    raise Unsupported('unsupported index label')


def _unlabel(value):
    if value.startswith('i') and len(value) == 21 and value[1:].isdigit():
        number = int(value[1:])
        if number < 1 << 64: return number - (1 << 63)
    if value.startswith('s'): return value[1:]
    raise TransportError('invalid native index label')


def _array(value):
    guards.plain_array(value)
    if type(value) is np.ndarray:
        return str(value.dtype), list(value)
    if type(value) is pd.arrays.StringArray:
        if value.dtype.storage != 'python': raise Unsupported('non-Python string storage')
        return 'string', list(value._ndarray)
    if type(value) in (pd.arrays.BooleanArray, pd.arrays.IntegerArray, pd.arrays.FloatingArray):
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


def _encode_request(sequence, operation, args):
    stream = io.StringIO(newline='')
    writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
    writer.writerow(('request', sequence, operation, len(args), '', PROFILE))
    for ordinal, value in enumerate(args):
        kind, dtype, labels, cells = _value(value)
        writer.writerow(('arg', ordinal, kind, dtype, '', ''))
        for i, label in enumerate(labels):
            writer.writerow(('index', ordinal, i, 'label', label, ''))
        for i, cell in enumerate(cells):
            tag, raw = _cell(cell)
            writer.writerow(('cell', ordinal, i, tag, raw, ''))
    return stream.getvalue().encode('utf-8')


def _request(sequence, operation, args):
    try: return _encode_request(sequence, operation, args)
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
    if header[0] != 'result' or header[2] not in ('scalar', 'series', 'array', 'list'):
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
    prototype_kind = next((kind for kind, _size in shapes if kind != 'scalar'), 'scalar')
    expected_kind = ('scalar' if operation in ('quit', 'sum', 'any') else
                     'list' if operation == 'sorted' else
                     'array' if operation in ('unique', 'zeros', 'assign') else prototype_kind)
    if kind != expected_kind: raise TransportError('native result kind differs from operation')
    if operation == 'quit' and cells != ['closed']:
        raise TransportError('native quit acknowledgment differs')
    if operation == 'zeros':
        if size != args[0]: raise TransportError('native zero-array denominator differs')
    elif operation not in ('quit', 'sum', 'any'):
        expected_size = next((count for kind, count in shapes if kind != 'scalar'), 1)
        if (operation == 'unique' and size > expected_size
                or operation != 'unique' and size != expected_size):
            raise TransportError('native column denominator differs from operands')
    if kind == 'scalar':
        if size != 1: raise TransportError('nonscalar result')
        return cells[0]
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
            or len(set(labels)) != len(labels)
            or set(labels) != {_label(item) for item in guards.axis_values(prototypes[0]._mgr.axes[0])}):
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
        self.sequence = 0
        self.metrics = dict(processes=0, wire_requests=0, native_operations=Counter(),
                            mature_fallbacks=Counter(), source_operations=Counter(),
                            transport_profile=PROFILE, cache_hits=0)

    def start(self):
        if self.process is not None: return
        if (hashlib.sha256(guards.regular_bytes(self.binary)).hexdigest() != self.digest
                or not os.access(self.binary, os.X_OK)):
            raise TransportError('native executable identity differs')
        environment = dict(os.environ, POLARS_MAX_THREADS='1', RAYON_NUM_THREADS='1')
        # No Python process, environment creation or shell interpretation.
        self.process = subprocess.Popen([str(self.binary), 'insilico-criteria', 'session'],
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

    def request(self, operation, args):
        body = _request(self.sequence, operation, args)  # admission before child IO
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

    def close(self, *, failed=False):
        if self.process is None: return
        process = self.process
        if failed:
            # Close only this owned transport. No original work is replayed.
            process.stdin.close()
            try: process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.terminate(); process.wait(timeout=5)
            process.stdout.close(); process.stderr.close()
            return
        self.request('quit', [])
        process.stdin.close()
        try: code = process.wait(timeout=5)
        except subprocess.TimeoutExpired as error: raise TransportError('native close timeout') from error
        remainder = process.stdout.read()
        diagnostic = process.stderr.read(65537)
        process.stdout.close(); process.stderr.close()
        if code or remainder or diagnostic: raise TransportError('native session closed with failure/extra output')
        if hashlib.sha256(guards.regular_bytes(self.binary)).hexdigest() != self.digest:
            raise TransportError('native executable changed during procedure')


class _Matrix:
    def __init__(self, columns): self.columns = columns


_BINARY = {'and': operator.and_, 'or': operator.or_, 'eq': operator.eq,
    'ne': operator.ne, 'lt': operator.lt, 'le': operator.le, 'gt': operator.gt, 'ge': operator.ge}


class _Runtime:
    def __init__(self, session, reference):
        self.session = session
        self.reference = reference
        self.fallback_functions = {}
        self.force_mature = False
        self.binding_failures = []

    def after_log(self, value):
        # Handlers may mutate a later column, the confidence mapping or a
        # scientific binding. Revalidate AFTER the actual original logger call;
        # never speculate across that callback. A changed binding makes the
        # remaining exact source expressions mature, without replaying a prefix.
        try: guards.reference_guard(self.reference)
        except (OSError, ValueError, TypeError, AttributeError) as failure:
            self.force_mature = True
            self.binding_failures.append(str(failure))
        return value

    def fallback(self, site, function_name, operation, function, args, kwargs):
        self.session.metrics['mature_fallbacks'][operation] += 1
        # Warnings emitted by mature pandas must still identify the actual
        # original call site, not this adapter's file. No original call/getter
        # is repeated; the already-evaluated operands enter this source site.
        key = (site, function_name)
        if key not in self.fallback_functions:
            expression = ast.parse('def call(__function, __args, __kwargs):\n return __function(*__args, **__kwargs)\n').body[0]
            expression.name = function_name
            ast.increment_lineno(expression, site - 2)
            compiled = compile(ast.Module(body=[expression], type_ignores=[]), str(self.reference.__file__), 'exec')
            code = next(value for value in compiled.co_consts if type(value) is CodeType)
            self.fallback_functions[key] = FunctionType(code, vars(self.reference), function_name)
        return self.fallback_functions[key](function, args, kwargs)

    def calculate(self, site, name, operation, args, original, kwargs=None):
        self.session.metrics['source_operations'][operation] += 1
        if self.force_mature:
            return self.fallback(site, name, operation, original, args, kwargs or {})
        try:
            if (operation in ('and', 'or') and len(args) == 2
                    and all(type(value) is pd.Series for value in args)):
                # 27 October 2026 CT: pandas owns label alignment, including
                # its actual Index.join/cache/error effects. Until that whole
                # boundary is represented, differing ordered axes execute the
                # original operator once here, with these same operands.
                try:
                    labels = [tuple(_label(item) for item in guards.axis_values(value._mgr.axes[0]))
                              for value in args]
                except (ValueError, TypeError, AttributeError) as error:
                    raise Unsupported('unrepresented logical axis ownership') from error
                if labels[0] != labels[1]:
                    raise Unsupported('mature logical Series alignment')
            return self.session.request(operation, args)
        except Unsupported:
            return self.fallback(site, name, operation, original, args, kwargs or {})

    def binary(self, site, name, operation, left, right):
        return self.calculate(site, name, operation, [left, right], _BINARY[operation])

    def unary(self, site, name, operation, value):
        return self.calculate(site, name, operation, [value], operator.invert)

    def unary_call(self, site, name, operation, function, value):
        return self.calculate(site, name, operation, [value], function)

    def method(self, site, name, receiver, method, *args, **kwargs):
        if type(receiver) is _Matrix and method == 'max':
            if args or kwargs != {'axis': 1}: raise TransportError('unexpected source maximum options')
            return self.calculate(site, name, 'max_columns', receiver.columns,
                lambda *columns: pd.concat(columns, axis=1).max(axis=1))
        original = getattr(receiver, method)
        if self.force_mature:
            return self.fallback(site, name, method, original, args, kwargs)
        operation = method
        operands = [receiver, *args]
        if method in ('contains', 'strip'):
            from pandas.core.strings.accessor import StringMethods
            if type(receiver) is not StringMethods:
                return self.fallback(site, name, method, original, args, kwargs)
            operands[0] = receiver._orig
            allowed = {} if method == 'strip' else {'regex': True}
            if kwargs and kwargs != allowed:
                return self.fallback(site, name, method, original, args, kwargs)
        elif method == 'astype':
            if args != (str,) or kwargs:
                return self.fallback(site, name, method, original, args, kwargs)
            operation, operands = 'as_str', [receiver]
        elif method == 'map':
            if len(args) != 1 or type(args[0]) is not dict or kwargs != {'na_action': 'ignore'}:
                return self.fallback(site, name, method, original, args, kwargs)
            mapping = args[0]
            operands = [receiver, list(mapping), list(mapping.values())]
        elif method not in ('fillna', 'notna', 'isna', 'sum', 'any', 'unique') or kwargs:
            return self.fallback(site, name, method, original, args, kwargs)
        # A Boolean fill of floating storage is not a numeric 0/1 fill in
        # pandas. It can upcast/downcast with a warning, or nullable Float64
        # can reject it. Preserve that actual mature call/failure at this
        # frontier; this operation is explicitly mixed, not native coverage.
        if (method == 'fillna' and type(receiver) is pd.Series and len(args) == 1
                and type(args[0]) in (bool, np.bool_)):
            dtype = receiver.dtype
            if ((type(dtype) is type(np.dtype('float64')) and dtype == np.dtype('float64'))
                    or type(dtype) is pd.Float64Dtype):
                return self.fallback(site, name, method, original, args, kwargs)
        # Object fillna can downcast with an installed-pandas FutureWarning.
        # That mature operation stays at its actual source site, explicitly
        # mixed, until separately ported; no suppression/handwritten warning.
        if method == 'fillna' and type(receiver) is pd.Series and receiver.dtype == np.dtype('O'):
            raw = receiver._mgr.blocks[0].values
            missing = [value is None or value is pd.NA or type(value) in (float, np.float64) and np.isnan(value) for value in raw]
            replacement = args[0] if len(args) == 1 else _MISSING
            filled = [replacement if absent else value for value, absent in zip(raw, missing)]
            # String/object output has no downcast warning. Only the mature
            # object->numeric/bool conversion boundary remains explicit here.
            if any(missing) and all(type(value) in (bool, int, float, np.bool_, np.int64, np.float64) for value in filled):
                return self.fallback(site, name, method, original, args, kwargs)
        self.session.metrics['source_operations'][operation] += 1
        try: return self.session.request(operation, operands)
        except Unsupported:
            return self.fallback(site, name, method, original, args, kwargs)

    def numeric(self, site, name, function, value, *, errors):
        if errors != 'coerce': raise TransportError('unexpected source numeric options')
        return self.calculate(site, name, 'numeric', [value], lambda value: function(value, errors=errors))

    def concat(self, site, name, function, values, *, axis):
        if self.force_mature or axis != 1 or type(values) is not list or any(type(value) is not pd.Series for value in values):
            return self.fallback(site, name, 'concat', function, [values], {'axis': axis})
        # This is ownership of whole column operands; maximum/reduction remains
        # entirely native. It neither filters nor computes a scientific mask.
        return _Matrix(values)

    def zeros(self, site, name, function, count, *, dtype):
        if dtype is not int: raise TransportError('unexpected criterion output dtype')
        return self.calculate(site, name, 'zeros', [count], lambda count: function(count, dtype=dtype))

    def sorted(self, site, name, function, value):
        return self.calculate(site, name, 'sorted', [value], function)

    def assign(self, site, name, array, mask, value):
        def original(array, mask, value):
            array[mask] = value
            return array
        if not self.force_mature and type(mask) is pd.Series and len(mask._mgr.blocks) == 1:
            raw = mask._mgr.blocks[0].values
            if type(raw) is np.ndarray and raw.dtype == np.dtype('bool'):
                # NumPy's original array[Series] assignment performs this
                # ownership conversion, including Series/axis protocol calls.
                # Use its actual result as the native positional operand; do
                # not populate cache keys or discard a ceremonial conversion.
                mask = np.asarray(mask)
        return self.calculate(site, name, 'assign', [array, mask, value], original)


class _NativeCalls(ast.NodeTransformer):
    """Replace table expressions, preserving every original statement frontier."""
    def __init__(self): self.function_name = None; self.changes = Counter()

    def visit_FunctionDef(self, node):
        old = self.function_name
        self.function_name = node.name
        # 27 October 2026 CT: annotations/defaults are original CPython
        # definition work, not table expressions. In particular float | None
        # must stay a type union at the nested definition's original frontier.
        node.body = [self.visit(statement) for statement in node.body]
        self.function_name = old
        return node

    def call(self, node, procedure, args, keywords=()):
        self.changes[procedure] += 1
        # CPython attributes a chained method CALL to the method's attribute
        # line, which can differ from the receiver expression's first line.
        # The original LoF fillna CALL/co_lines is39, not AST receiver line37.
        site = node.func.end_lineno if procedure == 'method' else node.lineno
        result = ast.Call(func=ast.Attribute(value=ast.Name(id=_GLOBAL, ctx=ast.Load()), attr=procedure, ctx=ast.Load()),
            args=[ast.Constant(site), ast.Constant(self.function_name), *args], keywords=list(keywords))
        return ast.copy_location(result, node)

    def visit_BinOp(self, node):
        node = self.generic_visit(node)
        operations = {ast.BitAnd: 'and', ast.BitOr: 'or'}
        if type(node.op) in operations:
            return self.call(node, 'binary', [ast.Constant(operations[type(node.op)]), node.left, node.right])
        return node

    def visit_Compare(self, node):
        node = self.generic_visit(node)
        operations = {ast.Eq: 'eq', ast.NotEq: 'ne', ast.Lt: 'lt', ast.LtE: 'le', ast.Gt: 'gt', ast.GtE: 'ge'}
        if len(node.ops) == 1 and type(node.ops[0]) in operations:
            return self.call(node, 'binary', [ast.Constant(operations[type(node.ops[0])]), node.left, node.comparators[0]])
        return node

    def visit_UnaryOp(self, node):
        node = self.generic_visit(node)
        if type(node.op) is ast.Invert:
            return self.call(node, 'unary', [ast.Constant('invert'), node.operand])
        return node

    def visit_Call(self, node):
        node = self.generic_visit(node)
        function = node.func
        if type(function) is ast.Attribute and type(function.value) is ast.Name:
            if function.value.id == 'pd' and function.attr in ('to_numeric', 'concat'):
                return self.call(node, 'numeric' if function.attr == 'to_numeric' else 'concat', [function, *node.args], node.keywords)
            if function.value.id == 'np' and function.attr == 'logical_not':
                return self.call(node, 'unary_call', [ast.Constant('not'), function, *node.args], node.keywords)
            if function.value.id == 'np' and function.attr == 'zeros':
                return self.call(node, 'zeros', [function, *node.args], node.keywords)
            if function.value.id == 'logger':
                # Original logger call/locations remain inside this outer
                # post-callback check, rather than moving emission into us.
                result = ast.Call(func=ast.Attribute(value=ast.Name(id=_GLOBAL, ctx=ast.Load()), attr='after_log', ctx=ast.Load()), args=[node], keywords=[])
                return ast.copy_location(result, node)
        if type(function) is ast.Name and function.id == 'sorted':
            return self.call(node, 'sorted', [function, *node.args], node.keywords)
        if type(function) is ast.Attribute and function.attr in (
                'contains', 'fillna', 'astype', 'strip', 'notna', 'isna', 'sum', 'any', 'unique', 'map', 'max'):
            return self.call(node, 'method', [function.value, ast.Constant(function.attr), *node.args], node.keywords)
        return node

    def visit_Assign(self, node):
        node = self.generic_visit(node)
        if len(node.targets) == 1 and type(node.targets[0]) is ast.Subscript:
            target = node.targets[0]
            if type(target.value) is ast.Name and target.value.id in ('pp3_array', 'bp4_array', 'bp7_array'):
                call = self.call(node, 'assign', [ast.Name(id=target.value.id, ctx=ast.Load()), target.slice, node.value])
                return ast.copy_location(ast.Assign(targets=[ast.Name(id=target.value.id, ctx=ast.Store())], value=call), node)
        return node


def transformed(reference):
    raw = guards.reference_guard(reference)
    nodes = copy.deepcopy(guards.source_tree(raw))
    transformer = _NativeCalls()
    nodes = [transformer.visit(node) for node in nodes]
    tree = ast.fix_missing_locations(ast.Module(body=nodes, type_ignores=[]))
    compiled = compile(tree, str(reference.__file__), 'exec')
    result = {}
    for code in compiled.co_consts:
        if type(code) is CodeType and code.co_name in ('PP3_BP4_criteria', 'BP7_criteria'):
            original = getattr(reference, code.co_name)
            function = FunctionType(code, vars(reference), code.co_name, original.__defaults__)
            result[code.co_name] = function
    if len(result) != 2 or not transformer.changes:
        raise RuntimeError('complete original procedure transformation required')
    return result, dict(transformer.changes)


class Bridge:
    """Explicit optional bridge. It never installs itself into the caller."""
    def __init__(self, *, binary, binary_sha256, reference=_reference, allow_original_fallback=True):
        self.binary = Path(binary)
        self.binary_sha256 = binary_sha256
        self.reference = reference
        self.allow_original_fallback = allow_original_fallback
        self.last_metrics = None

    def call(self, procedure, df, *args, **kwargs):
        if procedure not in ('PP3_BP4_criteria', 'BP7_criteria'):
            raise ValueError('unknown complete insilico procedure')
        original = getattr(self.reference, procedure)
        with _LOCK:
            self.last_metrics = dict(mode='original', whole_fallback=None, processes=0)
            try:
                functions, changes = transformed(self.reference)
                bound = __import__('inspect').signature(original).bind(df, *args, **kwargs)
                bound.apply_defaults()
                pvs1 = bound.arguments.get('pvs1_criteria')
                guards.frame_guard(df, pvs1)
                confidence = bound.arguments.get('high_confidence_status', {})
                if type(confidence) is not dict:
                    raise ValueError('custom confidence mapping')
                for key, value in confidence.items(): _cell(key); _cell(value)
                if hashlib.sha256(guards.regular_bytes(self.binary)).hexdigest() != self.binary_sha256:
                    raise ValueError('optional native executable changed')
            except (OSError, ValueError, TypeError, AttributeError) as failure:
                self.last_metrics['whole_fallback'] = str(failure)
                if not self.allow_original_fallback:
                    raise Unsupported(str(failure)) from failure
                return original(df, *args, **kwargs)
            session = NativeSession(self.binary, self.binary_sha256)
            runtime = _Runtime(session, self.reference)
            old = vars(self.reference).get(_GLOBAL, _MISSING)
            if old is not _MISSING:
                raise RuntimeError('original namespace has a conflicting runtime slot')
            vars(self.reference)[_GLOBAL] = runtime
            failed = True
            try:
                result = functions[procedure](df, *args, **kwargs)
                session.close()
                failed = False
                return result
            finally:
                if failed: session.close(failed=True)
                vars(self.reference).pop(_GLOBAL, None)
                self.last_metrics = dict(session.metrics, changes=changes,
                    mode='mixed' if session.metrics['mature_fallbacks'] else 'native',
                    completed=not failed, whole_fallback=None,
                    callback_binding_failures=runtime.binding_failures,
                    hosted_boundaries=['raw ownership/admission', 'mature getters including boolean selection',
                                       'string-accessor construction', 'source logger/formatting',
                                       'np.asarray ownership', 'Series/array reconstruction'])


def PP3_BP4_criteria(df, pvs1_criteria=None,
        high_confidence_status=_reference.PP3_BP4_criteria.__defaults__[1], *,
        insilico_binary=None, insilico_binary_sha256=None):
    if insilico_binary is None:
        return _reference.PP3_BP4_criteria(df, pvs1_criteria, high_confidence_status)
    return Bridge(binary=insilico_binary, binary_sha256=insilico_binary_sha256).call(
        'PP3_BP4_criteria', df, pvs1_criteria, high_confidence_status)


def BP7_criteria(df, *, insilico_binary=None, insilico_binary_sha256=None):
    if insilico_binary is None: return _reference.BP7_criteria(df)
    return Bridge(binary=insilico_binary, binary_sha256=insilico_binary_sha256).call('BP7_criteria', df)
