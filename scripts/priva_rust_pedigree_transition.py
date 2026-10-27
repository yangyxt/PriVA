"""Private v2 session: complete stateless transitions, bounded pure-outcome reuse.

Entry point is priva_rust_pedigree_filtration.py filter --protocol
priva-pedigree-transition-v2. V1 and the Python default remain separate controls.
No cached item contains a record, parser observation/getter, log or output action.
The actual current getter/log/write program remains in the original host facade.
"""
from collections import OrderedDict
import hashlib
from pathlib import Path
import select
import sys

VERSION = 'priva-pedigree-transition-v2'
MAX_ENTRIES = 4096
MAX_BYTES = 16 * 1024 * 1024
_PHASES = frozenset(('Idle', 'Patient', 'Control', 'Chrom', 'Proband', 'Father', 'Mother', 'DebugAuto', 'DebugDrop', 'Done'))
_SOURCES = frozenset(('UNKNOWN', 'PATERNAL', 'MATERNAL', 'DENOVO', 'INHERITED'))


class _Profile:
    """Strong immutable source/profile token; equality is session-local identity.

    Contents include the complete admitted bindings. A new admitted namespace
    creates a new token, so a bare object-id/digest cannot alias a prior profile.
    """
    __slots__ = ('contents',)
    def __init__(self, contents):
        object.__setattr__(self, 'contents', contents)
    def __setattr__(self, name, value):
        raise AttributeError('immutable pedigree profile')


def _size(value, seen=None):
    """Conservative owned key/outcome graph size, including shared token contents."""
    seen = set() if seen is None else seen
    if id(value) in seen:
        return 0
    seen.add(id(value))
    result = sys.getsizeof(value)
    if type(value) is tuple:
        result += sum(_size(item, seen) for item in value)
    elif type(value) is _Profile:
        result += _size(value.contents, seen)
    return result


class _Memo:
    def __init__(self, *, entries=MAX_ENTRIES, budget=MAX_BYTES):
        self.entries, self.budget = entries, budget
        self.values = OrderedDict()
        self.bytes = 0
        self.evictions = self.oversize = 0

    def get(self, key):
        found = self.values.get(key)
        if found is None:
            return None
        self.values.move_to_end(key)
        return found[0]

    def put(self, key, value):
        size = _size((key, value))
        if size > self.budget or self.entries == 0:
            self.oversize += 1
            return
        previous = self.values.pop(key, None)
        if previous is not None:
            self.bytes -= previous[1]
        while self.values and (len(self.values) >= self.entries or self.bytes + size > self.budget):
            _, (_, removed) = self.values.popitem(last=False)
            self.bytes -= removed
            self.evictions += 1
        self.values[key] = (value, size)
        self.bytes += size


class TransitionSession:
    """One real child; hits reuse exact pure outcomes, never claim fresh execution."""
    def __init__(self, base, binary, *, timeout=300, statistics=None):
        self.base = base
        self._statistics = statistics
        self._host_sources = tuple(Path(path) for path in (base.__file__, __file__))
        source_rows = []
        for path in self._host_sources:
            before = base._identity(path)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if base._identity(path) != before:
                raise base.NativeProtocolError('pedigree transition source changed during admission')
            source_rows.append((str(path), before, digest))
        self._host_source_rows = tuple(source_rows)
        self.wire = base.NativeSession(binary, timeout=timeout, protocol=VERSION)
        self.path, self.identity, self.sha256 = self.wire.path, self.wire.identity, self.wire.sha256
        self.starts = self.wire.starts
        self.admitted_namespaces = {}
        self._namespace_profiles = {}
        self._source_profile = (VERSION, self.sha256, self.identity,
            self._host_source_rows,
            tuple(sorted(base._EXPECTED_BODIES.items())),
            tuple(sys.version_info[:3]), base.pd.__version__, base.np.__version__, base.pysam.__version__)
        self._profile = _Profile((self._source_profile, ()))
        self._memo = _Memo()
        self._context = self._state = ()
        self._ordinal = 0
        self._poisoned = False
        self.sequence = 1
        self.operations = {'hello': 1}
        self.completed_operations = {'hello': 1}
        self.cache_hits = {}
        self.native_transitions = {}

    @property
    def closed(self):
        return self.wire.closed

    def admit_namespace(self, namespace, bindings):
        """Called only at the facade's unchanged source-admission frontier."""
        key = id(namespace)
        if key in self._namespace_profiles:
            old_namespace, old_bindings = self._namespace_profiles[key]
            if old_namespace is not namespace or len(old_bindings) != len(bindings) or any(left is not right for left, right in zip(old_bindings, bindings)):
                raise self.base.NativeProtocolError('pedigree admitted transition profile changed')
            return
        self._namespace_profiles[key] = (namespace, bindings)
        # Strong namespace references prevent id reuse; bindings remain complete
        # identity values, not a digest-only substitute for source equality.
        content = tuple((name, value[1]) for name, value in self._namespace_profiles.items())
        self._profile = _Profile((self._source_profile, content))

    def _fail(self, message):
        self._poisoned = True
        raise self.base.NativeProtocolError(message)

    def _available(self):
        if self._poisoned or self.wire.closed or self.wire.process.poll() is not None:
            self._fail('native pedigree transition process is not available')

    def _fields(self, fields):
        if type(fields) is not dict or any(type(key) is not str or not key or type(value) is not str for key, value in fields.items()):
            self._fail('invalid typed transition fields')
        if 'op' in fields or 'sequence' in fields:
            self._fail('reserved transition message fields')
        return tuple(sorted(fields.items()))

    def _number(self, value):
        if type(value) is not str or not value or len(value) > 20 or not value.isascii() or not value.isdecimal() or (len(value) > 1 and value[0] == '0'):
            self._fail('noncanonical decision integer')
        number = int(value)
        if number > (1 << 64) - 1:
            self._fail('decision integer exceeds native bound')
        return number

    def _validate_gt(self, fields, role):
        kind = fields.pop(role + '.kind', None)
        if kind == 'none':
            return
        if kind in ('int', 'bool', 'str'):
            value = fields.pop(role + '.value', None)
            if value is None or (kind == 'bool' and value not in ('0', '1')):
                self._fail('invalid typed scalar decision GT')
            if kind == 'int':
                body = value[1:] if value.startswith('-') else value
                if value != '0' and (not body or body.startswith('0') or not body.isascii() or not body.isdecimal()):
                    self._fail('noncanonical decision integer allele')
            return
        if kind not in ('tuple', 'list'):
            self._fail('unqualified typed decision GT')
        count = self._number(fields.pop(role + '.count', None))
        if count > len(fields):
            self._fail('decision GT count exceeds fields')
        for index in range(count):
            prefix = role + '.' + str(index)
            allele = fields.get(prefix + '.kind')
            if allele == 'none':
                fields.pop(prefix + '.kind')
                continue
            if allele not in ('int', 'bool', 'str'):
                self._fail('unqualified typed decision allele')
            self._validate_gt(fields, prefix)

    def _validate_state(self, state):
        if not self._context:
            if state:
                self._fail('decision state without initialized context')
            return ()
        fields = dict(state)
        phase = fields.pop('phase', None)
        index = self._number(fields.pop('index', None))
        chrom = fields.pop('chrom_type', None)
        if phase not in _PHASES or chrom not in ('', 'AUTOSOMAL', 'X', 'Y') or (phase not in ('Patient', 'Control') and index):
            self._fail('invalid complete decision state')
        context = dict(self._context)
        if phase in ('Patient', 'Control') and index >= self._number(context[('patients' if phase == 'Patient' else 'controls') + '.count']):
            self._fail('decision participant index exceeds context')
        for role in ('proband', 'father', 'mother'):
            self._validate_gt(fields, role)
        if fields:
            self._fail('unknown complete decision state fields')
        return tuple(sorted(state.items()))

    def _decode(self, response, operation):
        state, outcome = {}, {}
        for key, value in response.items():
            if key.startswith('state.'):
                state[key[6:]] = value
            elif key.startswith('outcome.'):
                outcome[key[8:]] = value
            else:
                self._fail('foreign transition response field')
        next_state = self._validate_state(state)
        status = outcome.get('status')
        if status == 'source-error':
            if set(outcome) != {'status', 'class', 'message'} or outcome['class'] not in ('TypeError', 'IndexError'):
                self._fail('unqualified transition source exception')
        elif status == 'ok':
            value = {key: value for key, value in outcome.items() if key != 'status'}
            if operation in ('genotype', 'chromosome', 'route', 'source'):
                allowed = {'genotype': ('true', 'false'), 'chromosome': ('AUTOSOMAL', 'X', 'Y'),
                           'route': ('auto', 'X', 'Y', 'unknown'), 'source': _SOURCES}[operation]
                if set(value) != {'value'} or value['value'] not in allowed:
                    self._fail('invalid direct transition result')
            else:
                action = value.get('action')
                fields = {'patient-gt': {'action', 'sample'}, 'control-gt': {'action', 'sample'},
                    'source-gt': {'action', 'sample'}, 'constructor-chrom': {'action'}, 'debug-autosomal': {'action'},
                    'debug-no-patient': {'action'}, 'debug-control-reject': {'action', 'sample'},
                    'drop': {'action'}, 'write': {'action', 'value'}}
                if action not in fields or set(value) != fields[action] or (action == 'write' and value['value'] not in _SOURCES):
                    self._fail('invalid filter transition action')
        else:
            self._fail('missing transition outcome status')
        if operation in ('genotype', 'route', 'source') and next_state != self._state:
            self._fail('pure direct operation changed decision state')
        return next_state, tuple(sorted(outcome.items()))

    def _validate_context(self, fields):
        remaining = dict(fields)
        for name in ('proband', 'father', 'mother', 'sex'):
            if name not in remaining:
                self._fail('missing initialized context field')
            remaining.pop(name)
        for role in ('father', 'mother'):
            kind = remaining.pop(role + '.original_kind', None)
            if kind not in ('none', 'str') or (kind == 'none' and fields[role] != ''):
                self._fail('invalid initialized parent kind')
        if remaining.pop('sex.original_kind', None) != 'str':
            self._fail('invalid initialized sex kind')
        for role in ('patients', 'controls'):
            count = self._number(remaining.pop(role + '.count', None))
            if count > len(remaining):
                self._fail('initialized participant count exceeds fields')
            for index in range(count):
                if remaining.pop(role + '.' + str(index), None) is None:
                    self._fail('missing initialized participant')
        if remaining:
            self._fail('foreign initialized context fields')

    def _initialize(self, fields):
        if self._context or self._state:
            self._fail('filter context already initialized')
        context = self._fields(fields)
        self._validate_context(fields)
        response = self.wire.request('kernel-init', **{'context.' + name: value for name, value in context})
        echoed = tuple(sorted((key[8:], value) for key, value in response.items() if key.startswith('context.')))
        state = {key[6:]: value for key, value in response.items() if key.startswith('state.')}
        if echoed != context or response.get('value') != 'ready' or set(response) != {'value'} | {'context.' + key for key, _ in context} | {'state.' + key for key in state}:
            self._fail('kernel initialization context/coverage changed')
        self._context = context
        self._state = self._validate_state(state)
        return {'value': 'ready'}

    def request(self, operation, **fields):
        self._available()
        self.sequence += 1
        self.operations[operation] = self.operations.get(operation, 0) + 1
        try:
            if operation in ('family', 'parents', 'quit'):
                result = self.wire.request(operation, **fields)
            elif operation == 'filter-init':
                result = self._initialize(fields)
            else:
                message = self._fields(fields)
                ordinal = None
                if operation == 'filter-begin':
                    if not self._context or dict(self._state).get('phase') not in ('Idle', 'Done') or set(fields) != {'ordinal'}:
                        self._fail('previous record action incomplete or filter uninitialized')
                    ordinal = self._number(fields['ordinal'])
                    if ordinal != self._ordinal + 1:
                        self._fail('record source ordinal changed')
                    message = ()
                elif operation not in ('genotype', 'chromosome', 'route', 'source', 'filter-reply'):
                    self._fail('unqualified logical transition operation')
                key = (self._profile, self._context, self._state, operation, message)
                cached = self._memo.get(key)
                if cached is None:
                    query = {'procedure': operation}
                    for prefix, values in (('context.', self._context), ('state.', self._state), ('message.', message)):
                        query.update((prefix + name, value) for name, value in values)
                    response = self.wire.request('transition', **query)
                    cached = self._decode(response, operation)
                    self.native_transitions[operation] = self.native_transitions.get(operation, 0) + 1
                    self._memo.put(key, cached)
                else:
                    # No request is pending: readable output here is unsolicited
                    # corruption/EOF, not a response that a cache hit may hide.
                    if select.select([self.wire.process.stdout], [], [], 0)[0]:
                        self._fail('unexpected native output before cached transition')
                    self.cache_hits[operation] = self.cache_hits.get(operation, 0) + 1
                self._state, outcome = cached
                if ordinal is not None:
                    self._ordinal = ordinal
                result = dict(outcome)
                status = result.pop('status')
                if status == 'source-error':
                    self.completed_operations[operation] = self.completed_operations.get(operation, 0) + 1
                    cls = {'TypeError': TypeError, 'IndexError': IndexError}[result['class']]
                    raise cls(result['message'])
            self.completed_operations[operation] = self.completed_operations.get(operation, 0) + 1
            return result
        except self.base.NativeProtocolError:
            self._poisoned = True
            raise

    def close(self, *, suppress=False):
        try:
            self.wire.close(suppress=suppress, _request=self.request)
        finally:
            if self._statistics is not None:
                self._statistics.update(self.metrics())
            if not suppress:
                for path, before, _digest in self._host_source_rows:
                    if self.base._identity(Path(path)) != before:
                        raise self.base.NativeProtocolError('pedigree transition host source changed during session')

    def __enter__(self):
        return self

    def __exit__(self, cls, value, traceback):
        self.close(suppress=cls is not None)

    def metrics(self):
        return dict(transport_profile=VERSION, logical_sequence_including_close=self.sequence,
            logical_operations=dict(self.operations), wire_sequence=self.wire.sequence,
            actual_wire_operations=dict(self.wire.operations), native_transition_executions=dict(self.native_transitions),
            reused_pure_outcomes=dict(self.cache_hits), cache_entries=len(self._memo.values), cache_key_outcome_bytes=self._memo.bytes,
            cache_entry_limit=self._memo.entries, cache_byte_limit=self._memo.budget, cache_evictions=self._memo.evictions,
            cache_oversize_uncached=self._memo.oversize,
            transport_scope='Explicit v2: cache hits reuse exact pure outcomes, not fresh Rust execution. Timeout/corruption detection timing differs from v1. No scientific fallback/replay.')
