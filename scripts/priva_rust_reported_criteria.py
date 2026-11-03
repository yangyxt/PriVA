"""Optional source-ordered reported-evidence procedures.

The Python default calls the historical function directly. Native opt-in uses
an ordinal session; unsupported admission returns to Python before effects.
After a public function starts, transport faults raise and NEVER replay it.

Pandas parsing, apply/row boxing, indexing, CPython set iteration/formatting and
math.log remain explicit compatibility effects. Native decisions cover all
thirteen interpreter families and their twenty-five subtypes. Rust selects
every next public-unit effect; no original function is lowered or rewritten
as a private Python control program.
"""
import math
import os
from pathlib import Path
import selectors
import stat
import struct
import subprocess
import tempfile
import time

import numpy as np
import pandas as pd

import acmg_ps3_bs3_pp5_bp6_reported as original
import mavedb_interpreter as interpreter
from priva_rust_reported_criteria_source import FAMILIES, PUBLIC, guard, compile_effects


PROTOCOL = 'reported-criteria-v1'
_DEFAULT = object()
_TIMEOUT = 30
_MAX_LINE = 16 * 1024 * 1024
_PRIMITIVES = (type(None), str, bool, int, float, np.bool_, np.int32, np.int64,
               np.uint32, np.uint64, np.float32, np.float64)


def _admit_text(value):
    if type(value) is str:
        try:
            value.encode('utf-8')
        except UnicodeError as error:
            raise ValueError('non-UTF-8 text requires original Python') from error


class NativeTransportError(RuntimeError):
    """No public-function replay is legal after this error."""


def _text(value):
    if type(value) is not str:
        raise NativeTransportError('nonstring native text payload')
    return value.encode('utf-8').hex()


def _untext(value):
    try:
        if type(value) is not str or value != value.lower() or len(value) % 2:
            raise ValueError('noncanonical hex')
        raw = bytes.fromhex(value)
        if raw.hex() != value:
            raise ValueError('noncanonical hex')
        return raw.decode('utf-8')
    except (ValueError, UnicodeError) as error:
        raise NativeTransportError('invalid native text payload') from error


def _bits(value):
    return '-' if value is None else struct.pack('>d', value).hex()


def _float(value):
    if value == '-':
        return None
    if type(value) is not str or len(value) != 16 or value != value.lower():
        raise NativeTransportError('invalid native f64 payload')
    try:
        raw = bytes.fromhex(value)
        if raw.hex() != value:
            raise ValueError('noncanonical hex')
        return struct.unpack('>d', raw)[0]
    except (ValueError, struct.error) as error:
        raise NativeTransportError('invalid native f64 payload') from error


def _bool(value):
    if value not in ('0', '1'):
        raise NativeTransportError('invalid native Boolean payload')
    return value == '1'


def _positions(value):
    if not value:
        return []
    try:
        result = [int(item) for item in value.split(',')]
        if any(number < 0 or str(number) != item for number, item in zip(result, value.split(','))):
            raise ValueError('noncanonical ordinal')
        return result
    except ValueError as error:
        raise NativeTransportError('invalid native ordinal list') from error


class Session:
    """One bounded native session, opened before any scientific host effect."""
    def __init__(self, binary, *, scratch=None):
        path = Path(binary)
        before = path.stat()
        if not stat.S_ISREG(before.st_mode) or not os.access(path, os.X_OK):
            raise ValueError('reported native executable unavailable')
        self._temporary = tempfile.TemporaryDirectory(prefix='priva-reported-', dir=scratch)
        self._error = open(Path(self._temporary.name) / 'stderr', 'w+b')
        self._process = None
        self.ordinal = 0
        self.operation_counts = {}
        self.program_counts = {}
        self.host_compatibility_counts = {}
        try:
            self._process = subprocess.Popen([str(path), 'reported-criteria', 'serve'],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self._error,
                bufsize=0, close_fds=True)
            os.set_blocking(self._process.stdin.fileno(), False)
            if self.exchange('hello', PROTOCOL) != [PROTOCOL]:
                raise NativeTransportError('reported native protocol changed')
            after = path.stat()
            fields = ('st_dev', 'st_ino', 'st_mode', 'st_size', 'st_mtime_ns', 'st_ctime_ns')
            if any(getattr(before, name) != getattr(after, name) for name in fields):
                raise NativeTransportError('native executable changed during startup')
        except BaseException:
            self.close()
            raise

    def exchange(self, operation, *fields):
        values = [str(self.ordinal), operation, *fields]
        if any(type(value) is not str or not value.isascii() or any(char in value for char in '\t\r\n') for value in values):
            raise NativeTransportError('invalid reported request fields')
        raw = ('\t'.join(values) + '\n').encode('ascii')
        if len(raw) > _MAX_LINE:
            raise NativeTransportError('reported request too large')
        process = self._process
        try:
            deadline = time.monotonic() + _TIMEOUT
            remaining = memoryview(raw)
            with selectors.DefaultSelector() as writable:
                writable.register(process.stdin, selectors.EVENT_WRITE)
                while remaining:
                    if not writable.select(max(0, deadline - time.monotonic())):
                        raise NativeTransportError('reported native request timeout')
                    try:
                        written = os.write(process.stdin.fileno(), remaining)
                    except BlockingIOError:
                        continue
                    if written <= 0:
                        raise NativeTransportError('reported native input closed')
                    remaining = remaining[written:]
            result = bytearray()
            with selectors.DefaultSelector() as selector:
                selector.register(process.stdout, selectors.EVENT_READ)
                while not result.endswith(b'\n'):
                    if not selector.select(max(0, deadline - time.monotonic())):
                        raise NativeTransportError('reported native response timeout')
                    chunk = os.read(process.stdout.fileno(), 65536)
                    if not chunk:
                        raise NativeTransportError('reported native session ended before response')
                    result.extend(chunk)
                    if len(result) > _MAX_LINE or b'\n' in result[:-1]:
                        raise NativeTransportError('invalid native response framing')
            row = result[:-1].decode('ascii').split('\t')
        except (OSError, UnicodeError) as error:
            raise NativeTransportError('reported native I/O failed') from error
        if len(row) < 2 or row[0] != str(self.ordinal):
            raise NativeTransportError('native reported source ordinal changed')
        self.ordinal += 1
        if row[1] != 'ok':
            message = _untext(row[2]) if len(row) == 3 and row[1] == 'error' else 'malformed native rejection'
            raise NativeTransportError(message)
        self.operation_counts[operation] = self.operation_counts.get(operation, 0) + 1
        return row[2:]

    def close(self):
        process = self._process
        if process is not None:
            if process.stdin is not None:
                process.stdin.close()
            try:
                process.wait(timeout=1)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            if process.stdout is not None:
                process.stdout.close()
            self._process = None
        if hasattr(self, '_error'):
            self._error.close()
        if hasattr(self, '_temporary'):
            self._temporary.cleanup()


def _raw_values(value):
    """Inspect ordinary storage without item-cache or custom dtype callbacks."""
    from pandas.core.internals.blocks import ExtensionBlock, NumpyBlock
    from pandas.core.internals.managers import BlockManager, SingleBlockManager
    if type(value) not in (pd.DataFrame, pd.Series):
        raise ValueError('custom frame/row requires original Python')
    manager = object.__getattribute__(value, '_mgr')
    expected = BlockManager if type(value) is pd.DataFrame else SingleBlockManager
    if type(manager) is not expected:
        raise ValueError('custom reported-evidence manager')
    for block in manager.blocks:
        if type(block) not in (NumpyBlock, ExtensionBlock):
            raise ValueError('custom reported-evidence block')
        raw = block.values
        if type(raw) is np.ndarray:
            arrays = (raw,)
        elif type(raw) in (pd.arrays.IntegerArray, pd.arrays.FloatingArray, pd.arrays.BooleanArray):
            if any(bool(item) for item in raw._mask):
                raise ValueError('missing extension values require original Python')
            arrays = (raw._data, raw._mask)
        elif type(raw) in (pd.arrays.StringArray, pd.arrays.NumpyExtensionArray):
            arrays = (raw._ndarray,)
        else:
            raise ValueError('custom reported-evidence extension storage')
        if any(type(array) is not np.ndarray for array in arrays):
            raise ValueError('custom reported-evidence backing array')
        for array in arrays:
            if array.dtype.kind == 'O':
                for item in array.flat:
                    if item is pd.NA or type(item) not in _PRIMITIVES:
                        raise ValueError('custom cell callback requires original Python')
                    _admit_text(item)
                    if type(item) is int and abs(item) > 2 ** 53:
                        raise ValueError('large Python integers require original conversion')


def _axis_values(axis):
    """Read labels without populating Index/RangeIndex lookup caches."""
    if type(axis) is pd.RangeIndex:
        return object.__getattribute__(axis, '_range')
    if type(axis) is pd.Index:
        raw = object.__getattribute__(axis, '_data')
        if type(raw) is np.ndarray and raw.ndim == 1:
            if all(type(item) in _PRIMITIVES for item in raw):
                return raw
    raise ValueError('custom reported-evidence axis')


def _metadata_families(frame):
    """Reject nonstring family equality before native family factorization.

    Python set membership has identity shortcuts for NaNs, while the source's
    later `== family` filtering does not. Such inputs retain that original
    program; converting them to a native string or integer group is unsafe.
    """
    manager = object.__getattribute__(frame, '_mgr')
    keys = list(_axis_values(manager.axes[0]))
    if 'assay_family' not in keys:
        return
    position = keys.index('assay_family')
    for block in manager.blocks:
        places = list(block._mgr_locs)
        if position not in places:
            continue
        raw = block.values
        if type(raw) is np.ndarray:
            values = raw[places.index(position)]
        elif type(raw) in (pd.arrays.StringArray, pd.arrays.NumpyExtensionArray):
            values = raw._ndarray
        else:
            raise ValueError('nonstring assay-family storage requires original Python')
        if any(type(item) is not str for item in values):
            raise ValueError('nonstring assay-family equality requires original Python')
        return
    raise ValueError('missing raw assay-family block')


def _column_inputs(frame, name):
    """The raw column and any stale cached Series, without a frame getter.

    Admission must cover the value the original getter will actually use, not
    silently substitute raw-manager values for an existing user cache.
    """
    manager = object.__getattribute__(frame, '_mgr')
    keys = list(_axis_values(manager.axes[0]))
    if name not in keys:
        raise ValueError('missing source column requires original failure')
    position = keys.index(name)
    raw = manager.iget_values(position)
    result = [raw]
    cache = object.__getattribute__(frame, '__dict__').get('_item_cache')
    if type(cache) is not dict:
        raise ValueError('custom source column cache')
    if name in cache:
        cached = cache[name]
        _raw_values(cached)
        if type(cached) is not pd.Series or cached.index is not frame.index:
            raise ValueError('custom cached source column')
        from priva_rust_pm4_bp3_packages import instance_guard
        instance_guard(cached, series=True)
        result.append(object.__getattribute__(cached, '_mgr').blocks[0].values)
    return result


def _cells(raw):
    if type(raw) is np.ndarray and raw.ndim == 1:
        return raw
    if type(raw) in (pd.arrays.StringArray, pd.arrays.NumpyExtensionArray):
        return raw._ndarray
    raise ValueError('unsupported native text column storage')


def _quiet_missing(value):
    # Integer bit inspection avoids invoking a NumPy comparison or changing its
    # error-state callbacks while deciding whether admission is supported.
    if value is None:
        return True
    if type(value) is float:
        bits = int.from_bytes(struct.pack('>d', value), 'big')
        return bits & 0x7ff0000000000000 == 0x7ff0000000000000 and bool(bits & 0x0008000000000000)
    return False


def _rank_admitted(value):
    if type(value) not in (type(None), bool, int, float):
        return False
    if type(value) is float:
        bits = int.from_bytes(struct.pack('>d', value), 'big')
        exponent = bits & 0x7ff0000000000000
        fraction = bits & 0x000fffffffffffff
        # A signaling NaN can trigger NumPy's original invalid callback. It
        # keeps that original path; ordinary NaNs and infinities are supported.
        return exponent != 0x7ff0000000000000 or fraction == 0 or bool(fraction & 0x0008000000000000)
    return True


def _text_column(frame, name):
    for raw in _column_inputs(frame, name):
        if any(type(value) is not str and not _quiet_missing(value) for value in _cells(raw)):
            raise ValueError('unsupported native text comparison requires original Python')


def _table_admission(name, values):
    if name in ('mavedb_interpretation_per_row', 'mavedb_score_interpretation') or (name == 'pick_interpreter_from_metadata' and len(values) == 3):
        metadata = values[1] if name != 'pick_interpreter_from_metadata' else values[2]
        if type(metadata) is pd.DataFrame:
            _text_column(metadata, 'urn')
    if name not in ('PS3_BS3_criteria', 'PP5_BP6_criteria'):
        return
    frame = values[0]
    if type(frame) is not pd.DataFrame:
        raise ValueError('native criteria require an ordinary frame')
    _text_column(frame, 'CLNSIG')
    if name == 'PS3_BS3_criteria':
        _column_inputs(frame, 'CLNREVSTAT')
        status = values[2] if len(values) == 3 else original.PS3_BS3_criteria.__defaults__[1]
        if type(status) is not dict or any(not _rank_admitted(value) for value in status.values()):
            raise ValueError('unsupported confidence conversion requires original Python')
    else:
        # Pandas alignment remains a mature effect, but unsupported Boolean
        # coercion/NA/scalar/index layouts keep the entire original function.
        # Duplicate/nonmonotonic labels on the SAME axis remain supported.
        for mask in values[1:]:
            if HostEffects._mask_values(mask) is None:
                raise ValueError('unsupported native exclusion mask')
            if len(mask) != len(frame) or (type(mask) is pd.Series and mask.index is not frame.index):
                raise ValueError('unsupported native exclusion alignment')


def _admit(name, values):
    if (pd.__version__, np.__version__) != ('2.2.3', '1.26.4'):
        raise ValueError('unqualified reported-evidence package context')
    # These existing source-derived guards protect actual package code/defaults
    # and recognized closures, rather than baselining mutable current methods.
    from priva_rust_pm4_bp3 import _package_profile
    from priva_rust_pm4_bp3_packages import guard as packages, instance_guard, logger_instance_guard
    from priva_rust_posterior import _package_guard as row_boxing
    binding, digest = _package_profile()
    packages(binding, digest)
    row_boxing()
    from priva_rust_reported_criteria_source import package_effect_guard
    package_effect_guard()
    logger_instance_guard(original.logger)
    current_logger = original.logger
    while current_logger is not None:
        state = object.__getattribute__(current_logger, '__dict__')
        if 'warning' in state or 'debug' in state:
            raise ValueError('per-instance reported logger callback requires original Python')
        current_logger = state.get('parent') if state.get('propagate') else None
    # PM4's canonical-source helper is reused for two additional logger methods;
    # these are named declarations, never identities captured after mutation.
    import logging
    from priva_rust_pm4_bp3_packages import plain_function, source_context
    context = source_context({'logging'}, {'logging': '43f86bbc08fdd5c7b6e697abffe1381f534a92bb32e1f1aee8360d6a142592a1'})
    for method in ('warning', 'debug'):
        plain_function(getattr(logging.Logger, method), 'logging:Logger.' + method, context)
    for position, value in enumerate(values):
        if type(value) in (pd.DataFrame, pd.Series):
            _raw_values(value)
            instance_guard(value, series=type(value) is pd.Series)
            state = object.__getattribute__(value, '__dict__')
            if any(key in state for key in ('apply', 'map', 'to_dict', 'to_string', 'head', 'drop_duplicates',
                    'sum', 'get', '__getitem__', '__setitem__', 'loc', 'iloc')):
                raise ValueError('per-instance mature API callback requires original Python')
            manager = object.__getattribute__(value, '_mgr')
            for axis in manager.axes:
                for label in _axis_values(axis):
                    _admit_text(label)
            if type(value) is pd.DataFrame or (name == 'mavedb_interpretation_per_row' and position == 0):
                # Mask Series are different: their ordinary duplicate or
                # nonmonotonic index remains a pandas alignment effect.
                keys = _axis_values(manager.axes[0])
                if any(type(key) is not str for key in keys) or len(set(keys)) != len(keys):
                    raise ValueError('custom/duplicate reported key axis')
        elif type(value) is dict:
            if any(type(key) is not str or type(item) not in _PRIMITIVES
                   or (type(item) is int and abs(item) > 2 ** 53) for key, item in value.items()):
                raise ValueError('custom reported mapping')
            for key, item in value.items():
                _admit_text(key)
                _admit_text(item)
        elif type(value) not in _PRIMITIVES:
            raise ValueError('custom reported argument')
        elif type(value) is int and abs(value) > 2 ** 53:
            raise ValueError('large integer conversion requires original Python')
        else:
            _admit_text(value)
    if name in ('mavedb_interpretation_per_row', 'mavedb_score_interpretation'):
        _metadata_families(values[1])
    _table_admission(name, values)


class NativeEffects:
    """Typed native decisions; mature effects stay at their source positions."""
    def __init__(self, session):
        self.session = session

    def _spec(self, spec):
        return [spec.assay_family, spec.subtype.value, spec.directionality.value, spec.scale.value,
                _bits(spec.neutral_band[0]), _bits(spec.neutral_band[1]),
                _text(str(spec.neutral_band)), _text(str(spec.units)), str(int(bool(spec.flexibility_requirement)))]

    def _result(self, row, spec):
        if len(row) != 7 or row[0] != 'return' or row[6] not in ('text', 'units'):
            raise NativeTransportError('invalid native interpreter result')
        # The same source AssaySpec is changed only after the original raw math
        # effect succeeded. Never replace it with a deserialized copy.
        if _bool(row[5]):
            spec.directionality = interpreter.Directionality.HIGHER_IS_FUNCTION
        return {'call': _untext(row[1]), 'standardized_score': _float(row[2]),
                'reason': _untext(row[3]), 'notes': spec.units if row[6] == 'units' else _untext(row[4])}

    def _to_log2(self, value):
        try:
            decision = self.session.exchange('log2-needed', _bits(value))
            if len(decision) != 1 or not _bool(decision[0]):
                return None
            return math.log(value, 2)
        except NativeTransportError:
            raise
        except Exception:
            # Exact original broad scientific conversion/libm catch. Native
            # transport errors have their own type and must never be hidden.
            return None

    def standardize(self, score, spec):
        row = self.session.exchange('standardize', str(int(score is None)), spec.scale.value)
        if row == ['none']:
            return None
        if row == ['float']:
            return float(score)
        if row == ['log2']:
            return self._to_log2(float(score))
        raise NativeTransportError('invalid standardization effect')

    def interpret_score(self, score, spec, gene_specific_rules=None):
        raw = self._spec(spec)
        rules = gene_specific_rules or {}
        row = self.session.exchange('begin', *raw, str(int(bool(rules))))
        if row[0] == 'return':
            return self._result(row, spec)
        if row != ['standardize']:
            raise NativeTransportError('invalid interpreter initial effect')
        value = self.standardize(score, spec)
        # Only Rust selects this Ephys getter, after standardization succeeded.
        request = self.session.exchange('gof-direction-needed', spec.assay_family, _bits(value))
        if request not in (['0'], ['1']):
            raise NativeTransportError('invalid native Ephys getter route')
        gof_dir = rules.get('gof_dir', '') if request == ['1'] else ''
        row = self.session.exchange('after-standardize', *raw, _bits(value), _text(str(gof_dir)))
        if row[0] == 'return':
            return self._result(row, spec)
        if len(row) != 2 or row[0] != 'raw':
            raise NativeTransportError('invalid interpreter raw effect')
        value = float(score)
        if row[1] == 'dg-negate':
            value = -value
        elif row[1] == 'binding-kd':
            next_row = self.session.exchange('binding-float', _bits(value))
            if next_row[0] == 'return':
                return self._result(next_row, spec)
            if next_row != ['binding-log10']:
                raise NativeTransportError('invalid Binding log effect')
            value = math.log(1.0 / value, 10)
        elif row[1] not in ('dg-identity', 'binding-identity'):
            raise NativeTransportError('unknown interpreter raw conversion')
        return self._result(self.session.exchange('after-raw', *raw, _bits(value), _text(str(gof_dir))), spec)

    def format_evidence(self, evidence):
        fields = []
        for item in evidence:
            fields.extend(map(_text, (f"{item['urn']}", f"{item['assay_family']}",
                f"{item['assay_subtype']}", f"{item['score']:.4f}", f"{item['call']}")))
        return self.render('render-evidence', *fields)

    def render(self, procedure, *fields):
        row = self.session.exchange(procedure, *fields)
        if len(row) != 1:
            raise NativeTransportError('invalid native formatted result')
        return _untext(row[0])


class HostEffects:
    """Finite mature-effect dispatcher; Rust owns the procedure state machine.

    No measured function is rewritten or executed as a private Python program.
    Each request is one named mature operation, primitive formatting/conversion,
    or one actual source publication. Source-owned conditions live in Rust.
    """
    def __init__(self, session, expressions):
        self.session = session
        self.native = NativeEffects(session)
        self.expressions = expressions
        self.failed_context = None

    def expression(self, name, context):
        code = self.expressions[name]
        namespace = interpreter.__dict__ if name.startswith(('pick.', 'band.')) else original.__dict__
        return eval(code, namespace, context)

    @staticmethod
    def done():
        return ['done']

    @staticmethod
    def boolean(value):
        return ['bool', str(int(bool(value)))]

    @staticmethod
    def textual(value):
        return ['text', _text(value)]

    @staticmethod
    def fields(arguments, count):
        if len(arguments) != count:
            raise NativeTransportError('native mature-effect arity changed')
        return arguments

    def run(self, name, arguments):
        names = {
            'mavedb_interpretation_per_row': ('row', 'metadata_df'),
            'mavedb_score_interpretation': ('df', 'mavedb_metadata'),
            'PS3_BS3_criteria': ('df', 'mavedb_metadata_tsv', 'high_confidence_status'),
            'PP5_BP6_criteria': ('df', 'clinvar_patho', 'clinvar_benign'),
            'pick_interpreter_from_metadata': ('urn', 'score', 'metadata_df'),
            'parse_neutral_band_from_metadata': ('meta',),
        }
        if name not in names or len(arguments) != len(names[name]):
            raise NativeTransportError('unknown native public procedure signature')
        context = dict(zip(names[name], arguments))
        context['_program'] = name
        context['_effects'] = []
        context['_row_publications'] = []
        context['_metadata_publications'] = []
        context['_metadata_final_publications'] = []
        context['_array_publications'] = []
        context['_result_publications'] = []
        reply = self.session.exchange('start', name)
        if len(reply) < 2 or not reply[0].isdigit() or str(int(reply[0])) != reply[0]:
            raise NativeTransportError('invalid native program identity')
        identity = reply[0]
        self.session.program_counts[name] = self.session.program_counts.get(name, 0) + 1
        try:
            while True:
                if reply[0] != identity or len(reply) < 2:
                    raise NativeTransportError('native context/source order changed')
                if reply[1] == 'return':
                    fields = reply[2:]
                    expected = {
                        'mavedb_interpretation_per_row': ['row'],
                        'mavedb_score_interpretation': ['df'],
                        'PS3_BS3_criteria': ['result_dict', 'df'],
                        'PP5_BP6_criteria': ['pp5_array', 'bp6_array'],
                        'pick_interpreter_from_metadata': ['result'],
                        'parse_neutral_band_from_metadata': ['_band_result'],
                    }[name]
                    if fields != expected or any(key not in context for key in fields):
                        raise NativeTransportError('invalid native return/publication plan')
                    if name == 'mavedb_interpretation_per_row' and context['_row_publications'] != ['MaveDB_PS3', 'MaveDB_BS3', 'MaveDB_score_interpretation']:
                        raise NativeTransportError('native return skipped row publications')
                    if name in ('PS3_BS3_criteria', 'PP5_BP6_criteria'):
                        expected_arrays = ['ps3_array', 'bs3_array'] if name == 'PS3_BS3_criteria' else ['pp5_array', 'bp6_array']
                        if context['_array_publications'] != expected_arrays:
                            raise NativeTransportError('native return skipped ordered strength writes')
                    if name == 'PS3_BS3_criteria' and context['_result_publications'] != ['clinvar_patho', 'clinvar_benign', 'PS3', 'BS3']:
                        raise NativeTransportError('native return skipped criteria dictionary writes')
                    if name == 'pick_interpreter_from_metadata':
                        if context.get('_early_na'):
                            if context['_metadata_publications'] or context['_metadata_final_publications']:
                                raise NativeTransportError('native NA route published later metadata')
                        else:
                            if context['_metadata_publications'] != ['metadata_flags', 'urn', 'n_syn_used', 'band_source', 'assay_family', 'subtype', 'neutral_band_parsed', 'MaveDB_BS3', 'MaveDB_PS3']:
                                raise NativeTransportError('native return skipped metadata writes')
                            call = context['result']['call']
                            expected = ['MaveDB_PS3'] if call in ('LoF', 'GoF') else ['MaveDB_BS3'] if call == 'Neutral' else []
                            if context['_metadata_final_publications'] != expected:
                                raise NativeTransportError('native return skipped final criteria flag')
                    values = tuple(context[key] for key in fields)
                    return values[0] if len(values) == 1 else values
                if reply[1] != 'effect' or len(reply) < 3:
                    raise NativeTransportError('invalid native source effect')
                if not context['_effects']:
                    first = {'mavedb_interpretation_per_row': 'row.missing', 'mavedb_score_interpretation': 'score.apply',
                        'PS3_BS3_criteria': 'criteria.result-init', 'PP5_BP6_criteria': 'reported.mask',
                        'pick_interpreter_from_metadata': 'pick.mode', 'parse_neutral_band_from_metadata': 'band.initial-text'}[name]
                    if reply[2] != first:
                        raise NativeTransportError('native program skipped its first source effect')
                result = self.effect(reply[2], reply[3:], context)
                context['_effects'].append(reply[2])
                reply = self.session.exchange('next', identity, *result)
        except BaseException:
            self.failed_context = context
            raise

    def score(self, df, metadata):
        # Parsing and the enclosing PS3 prefix may already have completed.
        # Unsupported parsed metadata falls back at this helper entry only.
        try:
            _raw_values(metadata)
            _metadata_families(metadata)
            _text_column(metadata, 'urn')
        except (ValueError, AttributeError, TypeError):
            counts = self.session.host_compatibility_counts
            counts['score-helper-admission'] = counts.get('score-helper-admission', 0) + 1
            return original.mavedb_score_interpretation(df, metadata)
        return self.run('mavedb_score_interpretation', (df, metadata))

    def _row(self, name, arguments, c):
        if name == 'row.missing':
            self.fields(arguments, 0)
            result = self.boolean(self.expression(name, c))
            c['_missing_result'] = result[1]
            return result
        if name == 'row.text':
            if len(arguments) != 1:
                if len(arguments) < 2 or arguments[1] not in ('high_confs', 'pvalues', 'scores'):
                    raise NativeTransportError('invalid native split publication')
                c[arguments[1]] = [_untext(value) for value in arguments[2:]]
            key = arguments[0]
            if key not in ('MaveDB_high_conf', 'MaveDB_pvalue', 'MaveDB_score', 'MaveDB_urn'):
                raise NativeTransportError('invalid source row text key')
            return self.textual(self.expression('row.text.' + key, c))
        if name == 'row.init-evidence':
            c['urn_sets'] = [_untext(value) for value in arguments]
            c['score_interpretation'] = []
            c['ps3_evidence'] = []
            c['bs3_evidence'] = []
        elif name == 'row.urn':
            urn, ordinal = self.fields(arguments, 2)
            if not ordinal.isdigit() or str(int(ordinal)) != ordinal:
                raise NativeTransportError('invalid native assay ordinal')
            c['i'] = int(ordinal)
            c['urn'] = _untext(urn)
            values = self.expression('row.urn-values', c)
            result, _mask = self._urn_mask(values, c['urn'])
            return ['bool', result[0]]
        elif name == 'row.metadata':
            self.fields(arguments, 0)
            column = self.expression('row.metadata-column', c)
            mask = self._selection_mask(column, c['urn'])
            c['meta_row'] = c['metadata_df'].loc[mask].iloc[0]
        elif name == 'row.rationale':
            self.fields(arguments, 0)
            c['interpretation'] = self.expression(name, c)
        elif name == 'row.family':
            self.fields(arguments, 0)
            c['assay_family'] = self.expression(name, c)
            if type(c['assay_family']) is not str:
                raise NativeTransportError('nonstring family bypassed metadata admission')
            return self.textual(c['assay_family'])
        elif name == 'row.subtype':
            self.fields(arguments, 0)
            c['assay_subtype'] = self.expression(name, c)
        elif name == 'row.append-rationale':
            self.fields(arguments, 0)
            c['score_interpretation'].append(c['interpretation'])
        elif name in ('row.index-score', 'row.index-high-conf', 'row.index-pvalue'):
            value, = self.fields(arguments, 1)
            target = {'row.index-score': 'score', 'row.index-high-conf': 'high_conf', 'row.index-pvalue': 'pvalue'}[name]
            c[target] = _untext(value)
        elif name in ('row.float-score', 'row.float-pvalue'):
            self.fields(arguments, 0)
            target = 'score_float' if name == 'row.float-score' else 'pvalue_float'
            try:
                c[target] = self.expression(name, c)
            except ValueError:
                if target == 'pvalue_float':
                    c[target] = np.nan
                return ['ValueError']
            return ['number', _bits(c[target])]
        elif name == 'row.interpret':
            self.fields(arguments, 0)
            c['result'] = self.run('pick_interpreter_from_metadata', (c['urn'], c['score_float'], c['metadata_df']))
            c['call'] = self.expression('row.call', c)
        elif name in ('row.interpreter-ps3', 'row.interpreter-bs3'):
            self.fields(arguments, 0)
            return self.boolean(self.expression(name, c))
        elif name == 'row.append-evidence':
            side, route = self.fields(arguments, 2)
            if side not in ('PS3', 'BS3') or route not in ('interpreted', 'pvalue'):
                raise NativeTransportError('invalid native evidence append')
            if route == 'interpreted':
                call, reason = c['call'], c['result'].get('reason', 'Unknown')
            else:
                rendered = self.session.exchange('render-pvalue-evidence', side, _text(f"{c['pvalue_float']:.4f}"))
                if len(rendered) != 2:
                    raise NativeTransportError('invalid native pvalue evidence format')
                call, reason = map(_untext, rendered)
            c['ps3_evidence' if side == 'PS3' else 'bs3_evidence'].append({
                'urn': c['urn'], 'assay_family': c['assay_family'], 'assay_subtype': c['assay_subtype'],
                'score': c['score_float'], 'call': call, 'reason': reason})
        elif name == 'row.missing-warning':
            self.fields(arguments, 0)
            original.logger.warning(self.native.render('render-missing-warning', _text(c['urn'])))
        elif name == 'row.variant':
            self.fields(arguments, 0)
            # Each original formatted value executes in order before Rust joins
            # the four strings; no vectorized/eager source getter pass.
            return ['texts', *(_text(f"{c['row'][key]}") for key in ('chrom', 'pos', 'ref', 'alt'))]
        elif name == 'row.family-sets':
            self.fields(arguments, 0)
            c['ps3_families'] = set(e['assay_family'] for e in c['ps3_evidence'])
            c['bs3_families'] = set(e['assay_family'] for e in c['bs3_evidence'])
            c['ps3_only_families'] = c['ps3_families'] - c['bs3_families']
            c['bs3_only_families'] = c['bs3_families'] - c['ps3_families']
            c['conflicted_families'] = c['ps3_families'] & c['bs3_families']
            return ['conflicts', *(_text(family) for family in c['conflicted_families'])]
        elif name == 'row.clean-evidence':
            if len(arguments) < 2:
                raise NativeTransportError('incomplete native clean-evidence plan')
            rows = c['ps3_evidence'] + c['bs3_evidence']
            c['_evidence_rows'] = rows
            c['_families'] = [_untext(value) for value in arguments[2:]]
            c['ps3_clean_evidence'] = self._selected(rows, arguments[0])
            c['bs3_clean_evidence'] = self._selected(rows, arguments[1])
            c['ps3_resolved_evidence'] = []
            c['bs3_resolved_evidence'] = []
        elif name == 'row.conflict':
            family, left, right, winner = self.fields(arguments, 4)
            family = int(family)
            if family < 0 or family >= len(c['_families']) or winner not in ('PS3', 'BS3', 'tie'):
                raise NativeTransportError('invalid native conflict family')
            c['family'] = c['_families'][family]
            c['family_ps3'] = self._selected(c['_evidence_rows'], left)
            c['family_bs3'] = self._selected(c['_evidence_rows'], right)
            if winner == 'PS3':
                c['ps3_resolved_evidence'].extend(c['family_ps3'])
            elif winner == 'BS3':
                c['bs3_resolved_evidence'].extend(c['family_bs3'])
        elif name == 'row.final-evidence':
            left, right = self.fields(arguments, 2)
            c['final_ps3_evidence'] = self._selected(c['_evidence_rows'], left)
            c['final_bs3_evidence'] = self._selected(c['_evidence_rows'], right)
        elif name == 'row.final-log-parts':
            mode, = self.fields(arguments, 1)
            if mode == 'both':
                ps3 = f"{c['ps3_only_families'] | set(e['assay_family'] for e in c['ps3_resolved_evidence'])}"
                bs3 = f"{c['bs3_only_families'] | set(e['assay_family'] for e in c['bs3_resolved_evidence'])}"
                conflicts = ''
            elif mode == 'conflict':
                ps3 = bs3 = ''
                conflicts = f"{c['conflicted_families']}"
            else:
                raise NativeTransportError('invalid native final log route')
            left = self.native.format_evidence(c['final_ps3_evidence'])
            right = self.native.format_evidence(c['final_bs3_evidence'])
            return ['log-parts', *map(_text, (ps3, bs3, conflicts, left, right))]
        elif name == 'row.publish':
            key, value = self.fields(arguments, 2)
            if key not in ('MaveDB_PS3', 'MaveDB_BS3', 'MaveDB_score_interpretation') or value not in ('true', 'false', 'nan'):
                raise NativeTransportError('invalid native row publication')
            if ((key in ('MaveDB_PS3', 'MaveDB_BS3') and value == 'nan')
                    or (key == 'MaveDB_score_interpretation' and value != 'nan')):
                raise NativeTransportError('native row publication scalar kind changed')
            if key == 'MaveDB_PS3' and c.get('_missing_result') == '0' and 'final_ps3_evidence' not in c:
                raise NativeTransportError('native row write preceded evidence resolution')
            if c.get('_missing_result') == '1' and value != ('nan' if key == 'MaveDB_score_interpretation' else 'false'):
                raise NativeTransportError('native missing-row value changed')
            if c.get('_missing_result') == '0' and key == 'MaveDB_score_interpretation':
                raise NativeTransportError('native nonmissing row lost rationale text')
            self._publication(c, '_row_publications', key, ('MaveDB_PS3', 'MaveDB_BS3', 'MaveDB_score_interpretation'))
            c['row'][key] = {'true': True, 'false': False, 'nan': np.nan}[value]
            c['_row_publications'].append(key)
        elif name == 'row.rationale-texts':
            self.fields(arguments, 0)
            if c['_row_publications'] != ['MaveDB_PS3', 'MaveDB_BS3']:
                raise NativeTransportError('native rationale conversion preceded criteria writes')
            return ['texts', *(_text(str(value)) for value in c['score_interpretation'])]
        elif name == 'row.publish-text':
            key, value = self.fields(arguments, 2)
            if key != 'MaveDB_score_interpretation':
                raise NativeTransportError('invalid native interpretation write')
            if c.get('_missing_result') != '0':
                raise NativeTransportError('native missing row published rationale text')
            self._publication(c, '_row_publications', key, ('MaveDB_PS3', 'MaveDB_BS3', 'MaveDB_score_interpretation'))
            c['row'][key] = _untext(value)
            c['_row_publications'].append(key)
        else:
            raise NativeTransportError('unknown mature row effect')
        return self.done()

    @staticmethod
    def _selected(rows, encoded):
        positions = _positions(encoded)
        if any(index >= len(rows) for index in positions):
            raise NativeTransportError('native evidence ordinal out of bounds')
        if len(set(positions)) != len(positions):
            raise NativeTransportError('native evidence ordinal duplicated')
        return [rows[index] for index in positions]

    @staticmethod
    def _publication(context, ledger, key, order):
        published = context[ledger]
        if len(published) >= len(order) or key != order[len(published)]:
            raise NativeTransportError('native publication was duplicated or reordered')
        # The caller records the key only AFTER its actual successful STORE.

    @staticmethod
    def _mask_values(value):
        if type(value) is pd.Series:
            manager = object.__getattribute__(value, '_mgr')
            if len(manager.blocks) != 1:
                return None
            raw = manager.blocks[0].values
            if type(raw) is pd.arrays.BooleanArray:
                if any(bool(item) for item in raw._mask):
                    return None
                raw = raw._data
        elif type(value) is np.ndarray:
            raw = value
        else:
            return None
        if type(raw) is not np.ndarray or raw.dtype != np.dtype('bool') or raw.ndim != 1:
            return None
        return raw

    @staticmethod
    def _boolean_data(fields, length, *, nullable=False):
        if len(fields) != length or any(field not in ('0', '1') for field in fields):
            raise NativeTransportError('native Boolean mask shape/type changed')
        data = np.array([field == '1' for field in fields], dtype=bool)
        return pd.arrays.BooleanArray(data, np.zeros(length, dtype=bool)) if nullable else data

    @staticmethod
    def _nullable(value):
        return type(value) is pd.Series and type(object.__getattribute__(value, '_mgr').blocks[0].values) is pd.arrays.BooleanArray

    def _urn_mask(self, raw, urn):
        fields = [_text(urn)]
        for value in _cells(raw):
            if type(value) is str:
                fields.extend(('str', _text(value)))
            elif _quiet_missing(value):
                fields.extend(('missing', ''))
            else:
                raise NativeTransportError('unsupported URN escaped pre-admission')
        result = self.session.exchange('urn-mask', *fields)
        if len(result) != len(raw) + 2 or result[0] not in ('0', '1'):
            raise NativeTransportError('invalid native URN selection')
        positions = [index for index, bit in enumerate(result[2:]) if bit == '1']
        if result[1] != (str(positions[0]) if positions else '-') or result[0] != str(int(not positions)):
            raise NativeTransportError('native first duplicate URN plan changed')
        return result, self._boolean_data(result[2:], len(raw), nullable=type(raw) is pd.arrays.StringArray)

    def _selection_mask(self, column, urn):
        raw = object.__getattribute__(column, '_mgr').blocks[0].values
        _selection, data = self._urn_mask(raw, urn)
        return column._construct_result(data, name=column.name)

    def _contains(self, expression, pattern, context):
        filled = self.expression(expression, context)
        # The source's fillna and actual str accessor construction still occur
        # in order. Rust replaces its fixed literal regexp decision, not the
        # mature dtype/index/attrs/flags materializer.
        accessor = filled.str
        raw = object.__getattribute__(filled, '_mgr').blocks[0].values
        texts = _cells(raw)
        if any(type(value) is not str for value in texts):
            raise NativeTransportError('unsupported text escaped pre-admission')
        bits = self.session.exchange('contains', pattern, *map(_text, texts))
        data = self._boolean_data(bits, len(texts), nullable=type(raw) is pd.arrays.StringArray)
        return accessor._wrap_result(data, fill_value=None, returns_string=False)

    def _confidence(self, context, pattern):
        mapped = self.expression('criteria.ranks.' + pattern, context)
        raw = object.__getattribute__(mapped, '_mgr').blocks[0].values
        if type(raw) is not np.ndarray or raw.ndim != 1:
            raise NativeTransportError('unsupported mapped rank storage')
        fields = []
        for value in raw:
            if value is None:
                rank = None
            elif type(value) in (bool, int, float, np.bool_, np.int64, np.float64):
                rank = float(value)
            else:
                raise NativeTransportError('unsupported mapped confidence rank')
            fields.extend(('1', _bits(rank)))
        bits = self.session.exchange('high-confidence', *fields)
        return mapped._construct_result(self._boolean_data(bits, len(raw)), name=mapped.name)

    def _logical(self, left, right, operation, *, exclude=False):
        from pandas.core.ops.common import get_op_result_name
        result_name = get_op_result_name(left, right)
        left, right = left._align_for_op(right, align_asobject=True)
        first, second = self._mask_values(left), self._mask_values(right)
        if first is None or second is None or len(first) != len(second):
            # A mature operation can create a previously unsupported coercion
            # after earlier rows were published. Execute ONLY this current
            # operator once; never replay an enclosing function or its effects.
            self.session.host_compatibility_counts['mask-coercion'] = self.session.host_compatibility_counts.get('mask-coercion', 0) + 1
            from pandas.core import ops
            from pandas.core.construction import extract_array
            import operator
            values = extract_array(right, extract_numpy=True, extract_range=True)
            if exclude:
                values = operator.invert(values)
            result = ops.logical_op(left._values, values, operator.and_ if operation == 'and' else operator.or_)
            return left._construct_result(result, name=result_name)
        fields = [] if exclude else [operation]
        if exclude:
            fields.append('1')
        for a, b in zip(first, second):
            fields.extend((str(int(a)), str(int(b))))
        bits = self.session.exchange('combine' if exclude else 'logical', *fields)
        data = self._boolean_data(bits, len(first), nullable=self._nullable(left) or self._nullable(right))
        return left._construct_result(data, name=result_name)

    def _tables(self, name, arguments, c):
        if name == 'score.apply':
            self.fields(arguments, 0)
            def mavedb_interpretation_per_row(row, metadata_df):
                return self.run('mavedb_interpretation_per_row', (row, metadata_df))
            c['df'] = c['df'].apply(mavedb_interpretation_per_row, axis=1, metadata_df=c['mavedb_metadata'])
        elif name == 'score.select':
            side, = self.fields(arguments, 1)
            if side not in ('PS3', 'BS3'):
                raise NativeTransportError('invalid source score selection')
            c['mavedb_' + side.lower() + '_recs'] = self.expression('score.select.' + side, c)
        elif name == 'score.info':
            side, = self.fields(arguments, 1)
            if side not in ('PS3', 'BS3'):
                raise NativeTransportError('invalid source score info')
            count = f"{c['df']['MaveDB_' + side].sum()}"
            table = f"{c['mavedb_' + side.lower() + '_recs'][:10].to_string()}"
            original.logger.info(self.native.render('render-score-info', side, _text(count), _text(table)))
        elif name == 'criteria.result-init':
            self.fields(arguments, 0)
            c['result_dict'] = {'PS3': [], 'BS3': [], 'clinvar_patho': [], 'clinvar_benign': []}
        elif name == 'criteria.clinvar':
            pattern, target = self.fields(arguments, 2)
            if (pattern, target) not in (('Pathogenic', 'clinvar_lof'), ('Benign', 'high_conf_benign')):
                raise NativeTransportError('invalid ClinVar source pattern')
            left = self._contains('criteria.text.' + pattern, pattern, c)
            right = self._confidence(c, pattern)
            c[target] = self._logical(left, right, 'and')
        elif name == 'criteria.result-mask':
            key, source = self.fields(arguments, 2)
            if (key, source) not in (('clinvar_patho', 'clinvar_lof'), ('clinvar_benign', 'high_conf_benign')):
                raise NativeTransportError('invalid source result mask')
            self._publication(c, '_result_publications', key, ('clinvar_patho', 'clinvar_benign', 'PS3', 'BS3'))
            c['result_dict'][key] = c[source]
            c['_result_publications'].append(key)
        elif name == 'criteria.zeros':
            target, = self.fields(arguments, 1)
            if target not in ('ps3_array', 'bs3_array', 'pp5_array', 'bp6_array'):
                raise NativeTransportError('invalid native zero-array target')
            c[target] = np.zeros(len(c['df']), dtype=int)
        elif name == 'criteria.metadata-truth':
            self.fields(arguments, 0)
            return self.boolean(c['mavedb_metadata_tsv'])
        elif name == 'criteria.read-table':
            self.fields(arguments, 0)
            c['mavedb_metadata'] = self.expression(name, c)
        elif name == 'criteria.drop-duplicates':
            self.fields(arguments, 0)
            c['mavedb_metadata'].drop_duplicates(subset=['urn'], inplace=True)
        elif name == 'criteria.metadata-info':
            self.fields(arguments, 0)
            count = f"{len(c['mavedb_metadata'])}"
            table = f"{c['mavedb_metadata'].head().to_string(index=False)}"
            original.logger.info(self.native.render('render-metadata-info', _text(count), _text(table)))
        elif name == 'criteria.score':
            self.fields(arguments, 0)
            c['df'] = self.score(c['df'], c['mavedb_metadata'])
        elif name == 'criteria.combine':
            side, target = self.fields(arguments, 2)
            if (side, target) not in (('PS3', 'ps3_criteria'), ('BS3', 'bs3_criteria')):
                raise NativeTransportError('invalid native criteria combination')
            c[target] = self._logical(c['clinvar_lof' if side == 'PS3' else 'high_conf_benign'],
                self.expression('criteria.score.' + side, c), 'or')
        elif name == 'criteria.alias':
            source, target = self.fields(arguments, 2)
            if (source, target) not in (('clinvar_lof', 'ps3_criteria'), ('high_conf_benign', 'bs3_criteria')):
                raise NativeTransportError('invalid source criteria alias')
            c[target] = c[source]
        elif name == 'criteria.assign':
            target, mask, strength = self.fields(arguments, 3)
            if (target, mask, strength) not in (('ps3_array', 'ps3_criteria', '3'), ('bs3_array', 'bs3_criteria', '3'),
                    ('pp5_array', 'pp5_criteria', '1'), ('bp6_array', 'bp6_criteria', '1')):
                raise NativeTransportError('invalid native ordered strength write')
            order = ('ps3_array', 'bs3_array') if c['_program'] == 'PS3_BS3_criteria' else ('pp5_array', 'bp6_array')
            self._publication(c, '_array_publications', target, order)
            c[target][c[mask]] = int(strength)
            c['_array_publications'].append(target)
        elif name == 'criteria.result-array':
            key, value = self.fields(arguments, 2)
            if (key, value) not in (('PS3', 'ps3_array'), ('BS3', 'bs3_array')):
                raise NativeTransportError('invalid native result publication')
            if c['_array_publications'] != ['ps3_array', 'bs3_array']:
                raise NativeTransportError('native result preceded strength writes')
            self._publication(c, '_result_publications', key, ('clinvar_patho', 'clinvar_benign', 'PS3', 'BS3'))
            c['result_dict'][key] = c[value]
            c['_result_publications'].append(key)
        elif name == 'reported.mask':
            pattern, exclude, target = self.fields(arguments, 3)
            if (pattern, exclude, target) not in (('athogenic', 'clinvar_patho', 'pp5_criteria'), ('enign', 'clinvar_benign', 'bp6_criteria')):
                raise NativeTransportError('invalid native reported mask')
            c[target] = self._logical(self._contains('reported.text.' + pattern, pattern, c), c[exclude], 'and', exclude=True)
        else:
            raise NativeTransportError('unknown mature table effect')
        return self.done()

    def _band(self, name, arguments, c):
        if name == 'band.initial-text':
            self.fields(arguments, 0)
            c['neutral_band_str'] = self.expression(name, c)
            return self.textual(c['neutral_band_str'])
        if name == 'band.strip-tokens':
            # CPython strip's whitespace set is retained explicitly; Rust owns
            # the comma split, NA classification and interval priority.
            c['parsed_text'] = [_untext(value).strip() for value in arguments]
            return ['texts', *map(_text, c['parsed_text'])]
        if name == 'band.float-tokens':
            if len(arguments) % 2:
                raise NativeTransportError('invalid native band conversion plan')
            try:
                values = []
                for kind, token in zip(arguments[::2], arguments[1::2]):
                    if kind not in ('nan', 'float'):
                        raise NativeTransportError('invalid native band conversion kind')
                    values.append(np.nan if kind == 'nan' else float(_untext(token)))
                c['parsed'] = values
            except (ValueError, SyntaxError, TypeError):
                return ['parse-error']
            return ['numbers', *map(_bits, values)]
        if name == 'band.flexibility':
            self.fields(arguments, 0)
            c['flexibility_req'] = c['meta'].get('flexibility_requirement', False)
        elif name == 'band.flex-truth':
            self.fields(arguments, 0)
            try:
                return self.boolean(c['flexibility_req'])
            except (ValueError, SyntaxError, TypeError):
                return ['parse-error']
        elif name == 'band.percentiles':
            self.fields(arguments, 0)
            c['q05'] = c['meta'].get('q05')
            c['q95'] = c['meta'].get('q95')
            return ['none-bits', str(int(c['q05'] is None)), str(int(c['q95'] is None))]
        elif name == 'band.float-percentiles':
            self.fields(arguments, 0)
            try:
                c['q05'] = float(c['q05'])
                c['q95'] = float(c['q95'])
            except (ValueError, TypeError):
                return ['parse-error']
            return ['numbers', _bits(c['q05']), _bits(c['q95'])]
        elif name == 'band.result':
            lo, hi = self.fields(arguments, 2)
            c['_band_result'] = (_float(lo), _float(hi))
        else:
            raise NativeTransportError('unknown mature band effect')
        return self.done()

    def _pick(self, name, arguments, c):
        if name == 'pick.mode':
            self.fields(arguments, 0)
            return self.boolean(hasattr(c['metadata_df'], 'loc'))
        if name == 'pick.frame-lookup':
            self.fields(arguments, 0)
            try:
                column = self.expression('pick.urn-column', c)
                mask = self._selection_mask(column, c['urn'])
                c['meta'] = c['metadata_df'][mask].iloc[0].to_dict()
            except (IndexError, KeyError):
                return ['not-found']
        elif name == 'pick.dict-lookup':
            self.fields(arguments, 0)
            c['meta'] = self.expression(name, c)
        elif name == 'pick.meta-truth':
            self.fields(arguments, 0)
            return self.boolean(c['meta'])
        elif name == 'pick.flags':
            self.fields(arguments, 0)
            meta = c['meta']
            c['flags'] = {'ambiguous': meta.get('ambiguous', False),
                'nonmonotonic': meta.get('nonmonotonic', False),
                'needs_gene_rules': meta.get('needs_gene_specific_rules', False),
                'flexibility_required': meta.get('flexibility_requirement', False)}
        elif name == 'pick.flag-truth':
            key, = self.fields(arguments, 1)
            if key not in ('ambiguous', 'nonmonotonic', 'needs_gene_rules'):
                raise NativeTransportError('invalid native quality flag')
            return self.boolean(c['flags'][key])
        elif name == 'pick.parse-band':
            self.fields(arguments, 0)
            c['lo'], c['hi'] = self.run('parse_neutral_band_from_metadata', (c['meta'],))
            return ['band', _bits(c['lo']), _bits(c['hi'])]
        elif name == 'pick.band-ready':
            self.fields(arguments, 0)
            c['neutral_band'] = (c['lo'], c['hi'])
        elif name == 'pick.family-subtype':
            self.fields(arguments, 0)
            try:
                c['family'] = c['meta'].get('assay_family')
                c['subtype_value'] = c['meta'].get('subtype_value')
            except (KeyError, AttributeError, TypeError) as error:
                c['_caught'] = error
                return ['construct-error']
            return ['family-subtype', *sum(([('str' if type(c[key]) is str else 'other'),
                (_text(c[key]) if type(c[key]) is str else '')] for key in ('family', 'subtype_value')), [])]
        elif name == 'pick.subtype':
            family, = self.fields(arguments, 1)
            if family not in FAMILIES or c['family'] != family:
                raise NativeTransportError('native family binding changed')
            try:
                c['subtype_enum_cls'] = interpreter.FAMILY_TO_SUBTYPE_ENUM.get(c['family'])
                c['subtype'] = c['subtype_enum_cls'][c['subtype_value']]
            except KeyError as error:
                c['_caught'] = error
                return ['KeyError']
            except (AttributeError, TypeError) as error:
                c['_caught'] = error
                return ['construct-error']
        elif name == 'pick.enum-strings':
            self.fields(arguments, 0)
            try:
                c['directionality_str'] = c['meta'].get('directionality', 'HIGHER_IS_FUNCTION')
                c['scale_str'] = c['meta'].get('scale', 'CENTERED_0')
            except (KeyError, AttributeError, TypeError) as error:
                c['_caught'] = error
                return ['construct-error']
        elif name == 'pick.enums':
            self.fields(arguments, 0)
            try:
                c['directionality'] = interpreter.Directionality[c['directionality_str']]
                c['scale'] = interpreter.Scale[c['scale_str']]
            except KeyError as error:
                c['_caught'] = error
                return ['KeyError']
            except (AttributeError, TypeError) as error:
                c['_caught'] = error
                return ['construct-error']
        elif name == 'pick.construct':
            self.fields(arguments, 0)
            try:
                c['flexibility_requirement'] = c['meta'].get('flexibility_requirement', False)
                c['spec'] = interpreter.AssaySpec(assay_family=c['family'], subtype=c['subtype'],
                    directionality=c['directionality'], scale=c['scale'], baseline=c['meta'].get('baseline', ''),
                    units=c['meta'].get('units', ''), neutral_band=c['neutral_band'],
                    flexibility_requirement=c['flexibility_requirement'])
            except (KeyError, AttributeError, TypeError) as error:
                c['_caught'] = error
                return ['construct-error']
        elif name == 'pick.interpret':
            self.fields(arguments, 0)
            c['result'] = self.native.interpret_score(c['score'], c['spec'])
        elif name == 'pick.publish':
            key, mode = self.fields(arguments, 2)
            if mode == 'true' and key in ('MaveDB_PS3', 'MaveDB_BS3'):
                if c['_metadata_publications'] != ['metadata_flags', 'urn', 'n_syn_used', 'band_source', 'assay_family', 'subtype', 'neutral_band_parsed', 'MaveDB_BS3', 'MaveDB_PS3']:
                    raise NativeTransportError('native criteria flag preceded metadata publications')
                expected = 'MaveDB_PS3' if c['result']['call'] in ('LoF', 'GoF') else 'MaveDB_BS3' if c['result']['call'] == 'Neutral' else None
                if key != expected or c['_metadata_final_publications']:
                    raise NativeTransportError('native criteria flag duplicated or changed')
                value = True
            elif mode == 'source':
                self._publication(c, '_metadata_publications', key, ('metadata_flags', 'urn', 'n_syn_used', 'band_source', 'assay_family', 'subtype', 'neutral_band_parsed', 'MaveDB_BS3', 'MaveDB_PS3'))
                values = {'metadata_flags': lambda: c['flags'], 'urn': lambda: c['urn'],
                    'n_syn_used': lambda: c['meta'].get('n_syn_used', 0),
                    'band_source': lambda: c['meta'].get('band_source', ''),
                    'assay_family': lambda: c['family'], 'subtype': lambda: c['subtype_value'],
                    'neutral_band_parsed': lambda: c['neutral_band'], 'MaveDB_BS3': lambda: False, 'MaveDB_PS3': lambda: False}
                if key not in values:
                    raise NativeTransportError('invalid native metadata publication')
                value = values[key]()
            else:
                raise NativeTransportError('invalid native metadata publication value')
            c['result'][key] = value
            if mode == 'source':
                c['_metadata_publications'].append(key)
            else:
                c['_metadata_final_publications'].append(key)
        elif name == 'pick.call':
            self.fields(arguments, 0)
            return self.textual(c['result']['call'])
        elif name == 'pick.failure-parts':
            kind, = self.fields(arguments, 1)
            formats = {
                'missing': lambda: (), 'ambiguous': lambda: (), 'rules': lambda: (),
                'band': lambda: (f"{c['meta'].get('neutral_bands')}", f"{c['meta'].get('q05')}", f"{c['meta'].get('q95')}"),
                'family': lambda: (f"{c['family']}",),
                'subtype': lambda: (f"{c['subtype_value']}", f"{c['family']}"),
                'enum': lambda: (f"{c['_caught']}", f"{c['directionality_str']}", f"{c['scale_str']}"),
                'construct': lambda: (f"{c['_caught']}", str(c['meta'])),
            }
            if kind not in formats:
                raise NativeTransportError('invalid native metadata failure route')
            return ['failure-parts', *map(_text, formats[kind]())]
        elif name == 'pick.na':
            kind, reason, notes = self.fields(arguments, 3)
            if kind == 'missing':
                note, flags = c['urn'], {}
            elif kind in ('ambiguous', 'rules'):
                note, flags = c['meta'].get('rationale', ''), c['flags']
            elif kind in ('band', 'enum', 'construct'):
                note, flags = _untext(notes), c['flags']
            elif kind in ('family', 'subtype'):
                note, flags = '', c['flags']
            else:
                raise NativeTransportError('unknown native NA construction')
            c['result'] = {'call': 'NA', 'standardized_score': None, 'reason': _untext(reason), 'notes': note,
                'metadata_flags': flags, 'urn': c['urn'], 'MaveDB_PS3': False, 'MaveDB_BS3': False}
            c['_early_na'] = True
        else:
            raise NativeTransportError('unknown mature metadata effect')
        return self.done()

    def effect(self, name, arguments, context):
        if name.startswith('row.'):
            return self._row(name, arguments, context)
        if name.startswith(('score.', 'criteria.', 'reported.')):
            return self._tables(name, arguments, context)
        if name.startswith('band.'):
            return self._band(name, arguments, context)
        if name.startswith('pick.'):
            return self._pick(name, arguments, context)
        if name == 'log':
            level, message = self.fields(arguments, 2)
            if level not in ('DEBUG', 'INFO'):
                raise NativeTransportError('invalid native source log level')
            {'DEBUG': original.logger.debug, 'INFO': original.logger.info}[level](_untext(message))
            return self.done()
        if name == 'raise.list-index':
            source, = self.fields(arguments, 1)
            if source not in ('scores', 'high_confs', 'pvalues'):
                raise NativeTransportError('invalid source list-index failure')
            context[source][context['i']]
            raise NativeTransportError('native missing-index premise disagrees with source list')
        raise NativeTransportError('unknown native mature effect')


def _call(name, args, binary, scratch):
    function = getattr(original, name)
    if binary is None:
        return function(*args)
    session = None
    try:
        raws = guard(original, interpreter)
        _admit(name, args)
        expressions = compile_effects(raws, {original.__name__: original.__file__, interpreter.__name__: interpreter.__file__})
        session = Session(binary, scratch=scratch)
    except (ValueError, OSError, ImportError, AttributeError, NativeTransportError):
        if session is not None:
            session.close()
        return function(*args)
    if name == 'PS3_BS3_criteria' and len(args) == 2:
        args = (*args, original.PS3_BS3_criteria.__defaults__[1])
    try:
        return HostEffects(session, expressions).run(name, args)
    finally:
        session.close()


def mavedb_interpretation_per_row(row, metadata_df, *, ranking_binary=None, scratch=None):
    return _call('mavedb_interpretation_per_row', (row, metadata_df), ranking_binary, scratch)


def mavedb_score_interpretation(df, mavedb_metadata, *, ranking_binary=None, scratch=None):
    return _call('mavedb_score_interpretation', (df, mavedb_metadata), ranking_binary, scratch)


def PS3_BS3_criteria(df, mavedb_metadata_tsv='', high_confidence_status=_DEFAULT, *, ranking_binary=None, scratch=None):
    args = (df, mavedb_metadata_tsv) if high_confidence_status is _DEFAULT else (df, mavedb_metadata_tsv, high_confidence_status)
    return _call('PS3_BS3_criteria', args, ranking_binary, scratch)


def PP5_BP6_criteria(df, clinvar_patho, clinvar_benign, *, ranking_binary=None, scratch=None):
    return _call('PP5_BP6_criteria', (df, clinvar_patho, clinvar_benign), ranking_binary, scratch)


def _admit_spec(spec):
    if type(spec) is not interpreter.AssaySpec:
        raise ValueError('custom assay spec requires original Python')
    fields = object.__getattribute__(spec, '__dict__')
    expected = {'assay_family', 'subtype', 'directionality', 'scale', 'baseline', 'units', 'neutral_band', 'flexibility_requirement'}
    if type(fields) is not dict or set(fields) != expected:
        raise ValueError('noncanonical assay state requires original Python')
    family = fields['assay_family']
    if type(family) is not str or family not in FAMILIES:
        raise ValueError('unknown assay family uses original registry behavior')
    enum_name = FAMILIES[family][0]
    if (type(fields['subtype']) is not getattr(interpreter, enum_name)
            or type(fields['directionality']) is not interpreter.Directionality
            or type(fields['scale']) is not interpreter.Scale):
        raise ValueError('noncanonical assay enums require original Python')
    band = fields['neutral_band']
    if (type(band) not in (tuple, list) or len(band) != 2
            or any(type(value) not in (bool, int, float) for value in band)):
        raise ValueError('unsupported neutral-band comparison requires original Python')
    if any(type(value) in (int, np.int32, np.int64, np.uint32, np.uint64) and abs(int(value)) > 2 ** 53 for value in band):
        raise ValueError('large exact-integer band requires original Python')
    if any(type(fields[key]) not in _PRIMITIVES for key in ('baseline', 'units', 'flexibility_requirement')):
        raise ValueError('custom assay field requires original Python')
    for key in ('assay_family', 'baseline', 'units', 'flexibility_requirement'):
        _admit_text(fields[key])


def _helper_call(name, arguments, binary, scratch):
    functions = {
        '_to_log2': interpreter._to_log2,
        'validate_spec': interpreter.validate_spec,
        'interpret_score': interpreter.interpret_score,
        'parse_neutral_band_from_metadata': interpreter.parse_neutral_band_from_metadata,
        'pick_interpreter_from_metadata': interpreter.pick_interpreter_from_metadata,
        'BaseInterpreter.standardize': interpreter.BaseInterpreter.standardize,
        'BaseInterpreter.call_from_band': interpreter.BaseInterpreter.call_from_band,
    }
    if name not in functions:
        raise ValueError('unknown reported interpreter procedure')
    original_function = functions[name]
    if binary is None:
        return original_function(*arguments)
    session = None
    try:
        raws = guard(original, interpreter)
        if name == '_to_log2' and type(arguments[0]) not in (type(None), bool, int, float):
            raise ValueError('NumPy scalar comparison requires original log helper')
        if name == 'BaseInterpreter.call_from_band' and type(arguments[1]) not in (type(None), bool, int, float):
            raise ValueError('NumPy scalar comparisons require original band helper')
        if name in ('validate_spec', 'interpret_score', 'BaseInterpreter.standardize', 'BaseInterpreter.call_from_band'):
            position = 0 if name == 'validate_spec' else 1 if name == 'interpret_score' else 2
            _admit_spec(arguments[position])
            primitive = [value for index, value in enumerate(arguments)
                         if index != position and not (name.startswith('BaseInterpreter.') and index == 0)]
            _admit(name, primitive)
        elif name == 'pick_interpreter_from_metadata':
            urn, score, metadata = arguments
            if type(urn) is not str:
                raise ValueError('custom URN requires original Python')
            if type(metadata) is dict:
                _admit(name, (urn, score))
                for key, value in metadata.items():
                    if type(key) is not str or type(value) is not dict:
                        raise ValueError('custom metadata lookup requires original Python')
                    _admit(name, (value,))
            else:
                _admit(name, arguments)
        else:
            _admit(name, arguments)
        expressions = compile_effects(raws, {original.__name__: original.__file__, interpreter.__name__: interpreter.__file__})
        session = Session(binary, scratch=scratch)
    except (ValueError, OSError, ImportError, AttributeError, NativeTransportError):
        if session is not None:
            session.close()
        return original_function(*arguments)
    native = NativeEffects(session)
    try:
        if name in ('parse_neutral_band_from_metadata', 'pick_interpreter_from_metadata'):
            return HostEffects(session, expressions).run(name, arguments)
        if name == '_to_log2':
            return native._to_log2(*arguments)
        if name == 'validate_spec':
            row = session.exchange('validate', *native._spec(arguments[0]))
            if len(row) != 1:
                raise NativeTransportError('invalid native spec validation')
            return _bool(row[0])
        if name == 'interpret_score':
            return native.interpret_score(*arguments)
        if name == 'BaseInterpreter.standardize':
            return native.standardize(*arguments[1:])
        if name == 'BaseInterpreter.call_from_band':
            _instance, value, spec = arguments
            row = session.exchange('band-call', *native._spec(spec), _bits(value))
            if len(row) != 1 or row[0] not in ('Unknown', 'Neutral', 'LoF', 'GoF'):
                raise NativeTransportError('invalid native band call')
            return row[0]
        raise NativeTransportError('native helper dispatch incomplete')
    finally:
        session.close()


def interpret_score(score, spec, gene_specific_rules=None, *, ranking_binary=None, scratch=None):
    return _helper_call('interpret_score', (score, spec, gene_specific_rules), ranking_binary, scratch)


def validate_spec(spec, *, ranking_binary=None, scratch=None):
    return _helper_call('validate_spec', (spec,), ranking_binary, scratch)


def _to_log2(value, *, ranking_binary=None, scratch=None):
    return _helper_call('_to_log2', (value,), ranking_binary, scratch)


def parse_neutral_band_from_metadata(meta, *, ranking_binary=None, scratch=None):
    return _helper_call('parse_neutral_band_from_metadata', (meta,), ranking_binary, scratch)


def pick_interpreter_from_metadata(urn, score, metadata_df, *, ranking_binary=None, scratch=None):
    return _helper_call('pick_interpreter_from_metadata', (urn, score, metadata_df), ranking_binary, scratch)
