"""Complete optional native pedigree filter with mature pandas/pysam I/O.

The original Python module remains the default and unsupported-source route.
Use ``filter --vcf VCF --ped-table PED --target-family FAMILY --output VCF
--binary PRIVA`` for an explicit native whole-file run. ``genotype``, ``chromosome``
and ``source`` expose individual procedures. A whole filter owns ONE persistent
native session; no subprocess, temporary file or VCF prefetch occurs per record.

2026-10-27 source contract: carrier accepts 1 OR 2; source requires complete 1.
Unknown actual pysam sample names raise when reached. Do not eagerly read GTs,
normalize ALT/ploidy/sex, sort records, or replay a failed written prefix.
"""
import argparse
import ast
import atexit
import builtins
import csv
import hashlib
import importlib
import io
import logging
import math
import os
from pathlib import Path
import select
import struct
import subprocess
import sys
import time
import types

import numpy as np
import pandas as pd
import pysam

from priva_rust_pedigree_transition import TransitionSession, VERSION as TRANSITION_PROTOCOL

PROTOCOL = 'priva-pedigree-v1'
MAX_FRAME = 64 * 1024 * 1024
_EXPECTED_BODIES = {'VariantSource.__init__': '568c725abd16882ca406567ee5e499b54a8e452b29175a5d79d76e66000d876a',
 'VariantSource._determine_autosomal_source': '38b779dbdb4e4261dfad60e1ce20e0b407c9ffde82d120849889707403e7d39d',
 'VariantSource._determine_chrX_source': 'c4963963d5422351186abb6f05af02cc412489853241ab262a0b1cfdacf08ab3',
 'VariantSource._determine_chrY_source': 'f043be81ab782973743870c64eebe80d22cd2ca07ce61ea59f377532fc11e6f1',
 'VariantSource._genotype_exists': 'bbd48f0d0503b34716dffdea41cafc3fdcde7a9fef645ede7d2a89fed6ffba4f',
 'VariantSource._get_chrom_type': '8dbc7bf554de7225f9bee8cdd928d2233f532a3d143a3d4fe8f60a22a40d4079',
 'VariantSource._get_genotype': '7d89d23826b832f770343419e39893ffa6270aba73406691355d192bca63a02a',
 'VariantSource._has_variant': 'fabe969baa26aaebdf8f318a12a9d89ff2c61b86e8a8d950f1dbd37627a35796',
 'VariantSource.determine_source': 'f9b13c12c58702a6a14d40080cdf77775e180231fa37cab8fe842bf280f1f29d',
 'is_homozygous_or_hemizygous': 'fd59b601c9890edca83fe7b6b992f932ba47a20f5e7a7d53a24c6f94f6a7e4ca',
 'is_variant_carrier': '2cc0c8a8349842bdcf0106ac801bbfc320de379a3de362d03ca4a9a4ec00fd82',
 'pedigree_filter': '2469c2be80bafbc66c912b53057990ac95dc694a4bcc43347547bf3fc902e075'}
_SESSIONS = {}
_EXECUTABLES = {}
_COLUMNS = ['#FamilyID', 'IndividualID', 'PaternalID', 'MaternalID', 'Sex', 'Phenotype']


class Unsupported(ValueError):
    """Admission failed before an observable operation; reference remains valid."""


class NativeProtocolError(RuntimeError):
    """A native failure is never converted into whole-file replay after effects."""


def _reference(namespace=None):
    return namespace if namespace is not None else vars(importlib.import_module('pedigree_filtration_per_fam'))


def _code_shape(code):
    return (code.co_code, code.co_names, code.co_varnames, code.co_argcount,
            code.co_kwonlyargcount, code.co_posonlyargcount, code.co_flags,
            code.co_freevars, code.co_cellvars, code.co_stacksize, code.co_exceptiontable,
            tuple(_code_shape(value) if isinstance(value, types.CodeType) else value
                  for value in code.co_consts))


def _admit(namespace):
    """Inspect bodies/callables without executing a second copy of the module."""
    if sys.version_info[:3] != (3, 12, 3):
        raise Unsupported('unqualified Python exception/parser profile')
    if (pd.__version__, np.__version__, pysam.__version__) != ('2.2.3', '1.26.4', '0.22.1'):
        raise Unsupported('unqualified parser packages')
    source_path = Path(namespace['__file__'])
    source = source_path.read_text(encoding='utf-8')
    tree = ast.parse(source)
    expected_code = {}
    def visit(code):
        if code.co_name != '<module>':
            expected_code[code.co_qualname] = _code_shape(code)
        for value in code.co_consts:
            if isinstance(value, types.CodeType):
                visit(value)
    visit(compile(source, str(source_path), 'exec', optimize=sys.flags.optimize))
    functions = {}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef):
            functions[node.name] = (node, namespace.get(node.name))
        if isinstance(node, ast.ClassDef) and node.name == 'VariantSource':
            cls = namespace.get(node.name)
            if type(cls) is not type:
                raise Unsupported('custom source class')
            for child in node.body:
                if isinstance(child, ast.FunctionDef):
                    functions['VariantSource.' + child.name] = (child, vars(cls).get(child.name))
    if set(functions) != set(_EXPECTED_BODIES):
        raise Unsupported('source definition membership changed')
    expected_class = {name.split('.')[1] for name in _EXPECTED_BODIES if name.startswith('VariantSource.')}
    expected_class.update(('__module__', '__dict__', '__weakref__', '__doc__'))
    if set(vars(namespace['VariantSource'])) != expected_class or namespace['VariantSource'].__bases__ != (object,):
        raise Unsupported('source class hooks changed')
    for name, (node, function) in functions.items():
        if hashlib.sha256(ast.get_source_segment(source, node).encode()).hexdigest() != _EXPECTED_BODIES[name]:
            raise Unsupported('source body changed: ' + name)
        if type(function) is not types.FunctionType or _code_shape(function.__code__) != expected_code[name]:
            raise Unsupported('source callable changed: ' + name)
        if function.__defaults__ is not None or function.__kwdefaults__ is not None:
            raise Unsupported('source defaults changed')
        if function.__globals__ is not namespace:
            raise Unsupported('source global namespace changed')
    if any(namespace.get(name) is not value for name, value in (('pd', pd), ('np', np), ('pysam', pysam))):
        raise Unsupported('source parser binding changed')
    if pd.read_csv.__module__ != 'pandas.io.parsers.readers' or type(namespace.get('logger')) is not logging.Logger:
        raise Unsupported('source parser/logger callable changed')
    if any(name in vars(namespace['logger']) for name in
           ('debug', 'info', 'error', 'isEnabledFor', '_log', 'makeRecord', 'handle')):
        raise Unsupported('source logger instance methods changed')
    original_builtins = namespace.get('__builtins__', vars(builtins))
    original_builtins = vars(original_builtins) if isinstance(original_builtins, types.ModuleType) else original_builtins
    if any(original_builtins.get(name) is not getattr(builtins, name) for name in ('set', 'len', 'all', 'str')):
        raise Unsupported('source builtin changed')


def _bindings(namespace):
    cls = namespace['VariantSource']
    functions = (namespace['is_homozygous_or_hemizygous'], namespace['is_variant_carrier'],
                 namespace['pedigree_filter'], *(vars(cls).get(name.split('.')[-1])
                     for name in _EXPECTED_BODIES if name.startswith('VariantSource.')))
    return (namespace['is_homozygous_or_hemizygous'], namespace['is_variant_carrier'],
            namespace['pedigree_filter'], cls, namespace['pd'], namespace['np'], namespace['pysam'],
            namespace['logger'], pd.read_csv,
            cls.__new__, cls.__getattribute__, cls.__setattr__,
            *(vars(namespace['logger']).get(name) for name in
                ('debug', 'info', 'error', 'isEnabledFor', '_log', 'makeRecord', 'handle')),
            *(item for function in functions for item in (function, getattr(function, '__code__', None),
                getattr(function, '__defaults__', None), getattr(function, '__kwdefaults__', None))))


def _admit_session(namespace, session, *, already_checked=False):
    key = id(namespace)
    current = _bindings(namespace)
    if key not in session.admitted_namespaces:
        if not already_checked:
            _admit(namespace)
        session.admitted_namespaces[key] = current
    elif any(left is not right for left, right in zip(session.admitted_namespaces[key], current)):
        raise NativeProtocolError('pedigree source bindings changed during session')
    if isinstance(session, TransitionSession):
        session.admit_namespace(namespace, current)


def _identity(path):
    def row(value):
        return (value.st_ino, value.st_mode, value.st_size, value.st_mtime_ns,
                value.st_ctime_ns, value.st_nlink, value.st_uid, value.st_gid)
    return row(path.lstat()), row(path.stat()), str(path.resolve(strict=True))


def _executable(path):
    path = Path(os.path.abspath(os.fspath(path)))
    before = _identity(path)
    if not path.is_file() or path.stat().st_size == 0 or not os.access(path, os.X_OK):
        raise NativeProtocolError('pedigree binary must be a nonempty executable')
    cached = _EXECUTABLES.get(str(path))
    if cached is not None and cached[0] == before:
        return path, before, cached[1]
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    if _identity(path) != before:
        raise NativeProtocolError('pedigree binary changed during admission')
    _EXECUTABLES[str(path)] = (before, digest.hexdigest())
    return path, before, digest.hexdigest()


def _frame(fields):
    if type(fields) is not dict or any(type(k) is not str or not k or type(v) is not str for k, v in fields.items()):
        raise NativeProtocolError('invalid pedigree request fields')
    stream = io.StringIO(newline='')
    writer = csv.writer(stream, lineterminator='\n')
    writer.writerow(('name', 'value'))
    writer.writerows(fields.items())
    value = stream.getvalue().encode('utf-8')
    if not 0 < len(value) <= MAX_FRAME:
        raise NativeProtocolError('pedigree frame exceeds qualified bound')
    return struct.pack('>Q', len(value)) + value


def _parse_frame(value):
    try:
        reader = csv.reader(io.StringIO(value.decode('utf-8'), newline=''), strict=True)
        if next(reader) != ['name', 'value']:
            raise ValueError('columns')
        result = {}
        for row in reader:
            if len(row) != 2 or not row[0] or row[0] in result:
                raise ValueError('duplicate/empty field')
            result[row[0]] = row[1]
        return result
    except (UnicodeError, csv.Error, ValueError, StopIteration) as error:
        raise NativeProtocolError('invalid native pedigree response') from error


class NativeSession:
    """One process, strict sequence and no retries; stderr stays off the wire."""
    def __init__(self, binary, *, timeout=300, protocol=PROTOCOL):
        self.path, self.identity, self.sha256 = _executable(binary)
        self.timeout = timeout
        self.sequence = 0
        self.closed = False
        self.admitted_namespaces = {}
        self.operations = {}
        self.completed_operations = {}
        env = os.environ.copy()
        env.setdefault('POLARS_MAX_THREADS', '1')
        self.process = subprocess.Popen([str(self.path), 'pedigree-filtration', 'session', '--protocol', protocol],
                                        stdin=subprocess.PIPE, stdout=subprocess.PIPE, env=env, bufsize=0)
        self.starts = 1
        os.set_blocking(self.process.stdin.fileno(), False)
        try:
            if self.request('hello', protocol=protocol) != {'value': protocol}:
                raise NativeProtocolError('unqualified pedigree handshake')
        except BaseException:
            self.close(suppress=True)
            raise

    def _read(self, count, deadline):
        result = bytearray()
        while len(result) < count:
            left = deadline - time.monotonic()
            if left <= 0 or not select.select([self.process.stdout], [], [], max(0, left))[0]:
                raise NativeProtocolError('native pedigree response timed out')
            block = os.read(self.process.stdout.fileno(), count - len(result))
            if not block:
                raise NativeProtocolError('native pedigree response ended early')
            result.extend(block)
        return bytes(result)

    def request(self, operation, **fields):
        if self.closed or self.process.poll() is not None:
            raise NativeProtocolError('native pedigree process is not available')
        self.sequence += 1
        self.operations[operation] = self.operations.get(operation, 0) + 1
        request = dict(fields, op=operation, sequence=str(self.sequence))
        value = _frame(request)
        try:
            view = memoryview(value)
            deadline = time.monotonic() + self.timeout
            while view:
                left = deadline - time.monotonic()
                if left <= 0 or not select.select([], [self.process.stdin], [], max(0, left))[1]:
                    raise NativeProtocolError('native pedigree request timed out')
                try:
                    written = os.write(self.process.stdin.fileno(), view)
                except BlockingIOError:
                    continue
                if written <= 0:
                    raise NativeProtocolError('native pedigree short request write')
                view = view[written:]
            deadline = time.monotonic() + self.timeout
            count, = struct.unpack('>Q', self._read(8, deadline))
            if not 0 < count <= MAX_FRAME:
                raise NativeProtocolError('invalid native pedigree response length')
            response = _parse_frame(self._read(count, deadline))
        except (OSError, BrokenPipeError) as error:
            raise NativeProtocolError('native pedigree transport failed') from error
        if response.pop('sequence', None) != str(self.sequence):
            raise NativeProtocolError('native pedigree response ordinal changed')
        status = response.pop('status', None)
        if status == 'source-error' and set(response) == {'class', 'message'}:
            cls = {'TypeError': TypeError, 'IndexError': IndexError}.get(response['class'])
            if cls is None:
                raise NativeProtocolError('unqualified source exception')
            self.completed_operations[operation] = self.completed_operations.get(operation, 0) + 1
            raise cls(response['message'])
        if status != 'ok':
            raise NativeProtocolError('native pedigree response status changed')
        self.completed_operations[operation] = self.completed_operations.get(operation, 0) + 1
        return response

    def close(self, *, suppress=False, _request=None):
        if self.closed:
            return
        error = None
        try:
            if self.process.poll() is None:
                if (_request or self.request)('quit') != {'value': 'closed'}:
                    raise NativeProtocolError('native pedigree close response changed')
            self.process.stdin.close()
            if self.process.wait(timeout=5) != 0:
                raise NativeProtocolError('native pedigree process failed')
            if self.process.stdout.read(1):
                raise NativeProtocolError('extra native pedigree output')
            if _identity(self.path) != self.identity:
                raise NativeProtocolError('native pedigree executable changed during session')
        except BaseException as observed:
            error = observed
            if self.process.poll() is None:
                self.process.terminate()
                try:
                    self.process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    self.process.kill()
                    self.process.wait(timeout=5)
        finally:
            self.closed = True
            self.process.stdin.close()
            self.process.stdout.close()
        if error is not None and not suppress:
            raise error

    def __enter__(self):
        return self

    def __exit__(self, cls, value, traceback):
        self.close(suppress=cls is not None)


def _session(binary=None, session=None):
    if session is not None:
        return session
    if not binary:
        raise Unsupported('native binary not selected')
    key = os.path.abspath(os.fspath(binary))
    if key not in _SESSIONS or _SESSIONS[key].closed:
        _SESSIONS[key] = NativeSession(key)
    return _SESSIONS[key]


def close_sessions():
    for session in _SESSIONS.values():
        session.close(suppress=True)
    _SESSIONS.clear()


atexit.register(close_sessions)


def _encode_gt(value, prefix='gt'):
    if value is None:
        return {prefix + '.kind': 'none'}
    if type(value) in (int, bool, str):
        kind = {int: 'int', bool: 'bool', str: 'str'}[type(value)]
        return {prefix + '.kind': kind, prefix + '.value': str(int(value)) if type(value) is bool else str(value)}
    if type(value) not in (tuple, list):
        raise Unsupported('custom genotype container')
    result = {prefix + '.kind': type(value).__name__, prefix + '.count': str(len(value))}
    for index, allele in enumerate(value):
        if allele is not None and type(allele) not in (int, bool, str):
            raise Unsupported('custom/unhashable genotype allele')
        stem = prefix + '.' + str(index)
        result[stem + '.kind'] = 'none' if allele is None else {int: 'int', bool: 'bool', str: 'str'}[type(allele)]
        if allele is not None:
            result[stem + '.value'] = str(int(allele)) if type(allele) is bool else str(allele)
    return result


def _value(response):
    if set(response) != {'value'}:
        raise NativeProtocolError('native scalar fields changed')
    return response['value']


def _predicate(name, value, *, binary=None, session=None, namespace=None):
    namespace = _reference(namespace)
    try:
        fields = _encode_gt(value)
        if session is None:
            _admit(namespace)
        client = _session(binary, session)
        _admit_session(namespace, client)
    except Unsupported:
        if name == 'homo':
            return namespace['is_homozygous_or_hemizygous'](value)
        if name == 'carrier':
            return namespace['is_variant_carrier'](value)
        method = namespace['VariantSource']._genotype_exists if name == 'exists' else namespace['VariantSource']._has_variant
        return method(_ReferencePredicates(namespace), value)
    result = _value(client.request('genotype', predicate=name, **fields))
    if result not in ('true', 'false'):
        raise NativeProtocolError('native Boolean value changed')
    return result == 'true'


class _ReferencePredicates:
    def __init__(self, namespace):
        self.namespace = namespace

    def _genotype_exists(self, genotype):
        return self.namespace['VariantSource']._genotype_exists(self, genotype)


def is_homozygous_or_hemizygous(genotype, *, binary=None, session=None, namespace=None):
    return _predicate('homo', genotype, binary=binary, session=session, namespace=namespace)


def is_variant_carrier(genotype, *, binary=None, session=None, namespace=None):
    return _predicate('carrier', genotype, binary=binary, session=session, namespace=namespace)


def _source_log(namespace, level, message, function, line):
    logger = namespace['logger']
    if logger.isEnabledFor(level):
        # This is an explicit source-map label, not a fabricated adapter stack.
        record = logger.makeRecord(logger.name, level, namespace['__file__'], line,
                                   message, (), None, function)
        logger.handle(record)


class VariantSource:
    """All original source procedures, with shared-session native predicates."""
    def __init__(self, record, proband, father, mother, proband_sex, *, binary=None, session=None, namespace=None):
        self._namespace = _reference(namespace)
        if binary and session is None:
            _admit(self._namespace)
        self._session = _session(binary, session) if binary or session else None
        if self._session is not None:
            _admit_session(self._namespace, self._session)
        self.record, self.proband, self.father, self.mother = record, proband, father, mother
        self.proband_sex = proband_sex
        self.chrom = self.record.chrom
        self.chrom_type = self._get_chrom_type()

    def _get_chrom_type(self):
        if self._session is None or type(self.chrom) is not str:
            return self._namespace['VariantSource']._get_chrom_type(self)
        return _value(self._session.request('chromosome', chrom=self.chrom))

    def _get_genotype(self, sample):
        # Actual __contains__ and .get behavior is the mature-parser boundary.
        return self._namespace['VariantSource']._get_genotype(self, sample)

    def _genotype_exists(self, genotype):
        return _predicate('exists', genotype, session=self._session, namespace=self._namespace)

    def _has_variant(self, genotype):
        return _predicate('has', genotype, session=self._session, namespace=self._namespace)

    def _route(self, chromosome_type):
        if self._session is None or type(chromosome_type) is not str or type(self.proband_sex) is not str:
            raise Unsupported('reference source dispatcher required')
        return _value(self._session.request('route', chrom_type=chromosome_type, sex=self.proband_sex))

    def determine_source(self):
        if self._session is None:
            return self._namespace['VariantSource'].determine_source(self)
        proband_gt = self._get_genotype(self.proband)
        if not self._has_variant(proband_gt):
            return 'UNKNOWN'
        father_gt = self._get_genotype(self.father)
        mother_gt = self._get_genotype(self.mother)
        try:
            route = self._route(self.chrom_type)
        except Unsupported:
            # No extra record access is repeated: the same supplied GTs remain.
            if self.chrom_type == 'AUTOSOMAL':
                return self._determine_autosomal_source(proband_gt, father_gt, mother_gt)
            if self.chrom_type == 'X':
                return self._determine_chrX_source(proband_gt, father_gt, mother_gt)
            if self.chrom_type == 'Y':
                return self._determine_chrY_source(proband_gt, father_gt)
            return 'UNKNOWN'
        if route == 'auto':
            return self._determine_autosomal_source(proband_gt, father_gt, mother_gt)
        if route == 'X':
            return self._determine_chrX_source(proband_gt, father_gt, mother_gt)
        if route == 'Y':
            return self._determine_chrY_source(proband_gt, father_gt)
        if route != 'unknown':
            raise NativeProtocolError('source dispatch route changed')
        return 'UNKNOWN'

    def _debug_autosomal(self, proband_gt, father_gt, mother_gt):
        message = (f'For variant {self.record.chrom}:{self.record.pos}:{self.record.id}:{self.record.ref}:{self.record.alts[0]}, '
                   f'proband {self.proband} GT: {proband_gt}, has_variant: {self._has_variant(proband_gt)}, '
                   f'father {self.father} GT: {father_gt}, has_variant: {self._has_variant(father_gt)}, '
                   f'mother {self.mother} GT: {mother_gt}, has_variant: {self._has_variant(mother_gt)}')
        _source_log(self._namespace, logging.DEBUG, message, '_determine_autosomal_source', 123)

    def _native_source(self, method, proband_gt, father_gt, mother_gt):
        if type(self.proband_sex) is not str:
            raise Unsupported('raw sex equality requires provided-GT reference method')
        fields = dict(_encode_gt(proband_gt, 'proband'), **_encode_gt(father_gt, 'father'), **_encode_gt(mother_gt, 'mother'))
        return _value(self._session.request('source', method=method, sex=self.proband_sex, **fields))

    def _determine_autosomal_source(self, proband_gt, father_gt, mother_gt):
        if self._session is None:
            return self._namespace['VariantSource']._determine_autosomal_source(self, proband_gt, father_gt, mother_gt)
        self._debug_autosomal(proband_gt, father_gt, mother_gt)
        try:
            fields = dict(_encode_gt(proband_gt, 'proband'), **_encode_gt(father_gt, 'father'), **_encode_gt(mother_gt, 'mother'))
        except Unsupported:
            # Custom objects were not admitted; do not duplicate the debug call.
            return self._autosomal_reference_after_debug(proband_gt, father_gt, mother_gt)
        return _value(self._session.request('source', method='auto', sex='', **fields))

    def _autosomal_reference_after_debug(self, proband_gt, father_gt, mother_gt):
        if not self._genotype_exists(father_gt) or not self._genotype_exists(mother_gt):
            return 'INHERITED' if self._has_variant(father_gt) or self._has_variant(mother_gt) else 'UNKNOWN'
        self._has_variant(proband_gt)
        father, mother = self._has_variant(father_gt), self._has_variant(mother_gt)
        if not father and not mother:
            return 'DENOVO'
        if father and not mother:
            return 'PATERNAL'
        if not father and mother:
            return 'MATERNAL'
        if father and mother:
            return 'INHERITED'
        return 'UNKNOWN'

    def _determine_chrX_source(self, proband_gt, father_gt, mother_gt):
        if self._session is None or type(self.proband_sex) is not str:
            return self._namespace['VariantSource']._determine_chrX_source(self, proband_gt, father_gt, mother_gt)
        if self._route('X') == 'auto':
            return self._determine_autosomal_source(proband_gt, father_gt, mother_gt)
        try:
            return self._native_source('X', proband_gt, father_gt, mother_gt)
        except Unsupported:
            return self._namespace['VariantSource']._determine_chrX_source(self, proband_gt, father_gt, mother_gt)

    def _determine_chrY_source(self, proband_gt, father_gt):
        if self._session is None or type(self.proband_sex) is not str:
            return self._namespace['VariantSource']._determine_chrY_source(self, proband_gt, father_gt)
        try:
            return self._native_source('Y', proband_gt, father_gt, None)
        except Unsupported:
            return self._namespace['VariantSource']._determine_chrY_source(self, proband_gt, father_gt)


def _vector(response, prefix):
    try:
        raw = response.pop(prefix + '.count')
        count = int(raw)
        if count < 0 or str(count) != raw:
            raise ValueError('noncanonical count')
        return [response.pop(prefix + '.' + str(index)) for index in range(count)]
    except (KeyError, ValueError) as error:
        raise NativeProtocolError('native vector coverage changed') from error


def _put_vector(fields, prefix, values):
    fields[prefix + '.count'] = str(len(values))
    for index, value in enumerate(values):
        if type(value) is not str:
            raise Unsupported('custom family sample identifier')
        fields[prefix + '.' + str(index)] = value


def _family(client, ped_df, target):
    if type(target) not in (str, int, bool, float, type(None)):
        raise Unsupported('target equality requires original pandas comparison')
    fields = dict(rows=str(len(ped_df)), **{'target.kind': 'str' if type(target) is str else 'other', 'target.value': target if type(target) is str else ''})
    for row, values in enumerate(ped_df.itertuples(index=False, name=None)):
        for column, value in enumerate(values):
            if type(value) is not str:
                raise Unsupported('custom parsed pedigree cell')
            fields[f'row.{row}.{column}'] = value
    response = client.request('family', **fields)
    if response.pop('action', None) != 'family':
        raise NativeProtocolError('native family action changed')
    positions = _vector(response, 'positions')
    patients, controls = _vector(response, 'patients'), _vector(response, 'controls')
    if set(response) != {'proband', 'sex'}:
        raise NativeProtocolError('native family fields changed')
    try:
        ordinals = [int(value) for value in positions]
        if any(str(value) != raw for value, raw in zip(ordinals, positions)) or ordinals != sorted(set(ordinals)) or any(not 0 <= value < len(ped_df) for value in ordinals):
            raise ValueError('family positions')
    except ValueError as error:
        raise NativeProtocolError('native family row order changed') from error
    mask = np.zeros(len(ped_df), dtype=bool)
    mask[ordinals] = True
    family = ped_df[pd.Series(mask, index=ped_df.index)]
    return family, set(patients), set(controls), response['proband'], response['sex'], bool(patients)


def _parents(client, family):
    fields = {}
    for role, column in (('father', 'PaternalID'), ('mother', 'MaternalID')):
        fields[role + '.count'] = str(len(family))
        for index, value in enumerate(family[column]):
            prefix = f'{role}.{index}'
            if type(value) is str:
                fields[prefix + '.kind'], fields[prefix + '.value'] = 'str', value
            elif type(value) is float and math.isnan(value):
                fields[prefix + '.kind'] = 'none'
            else:
                raise NativeProtocolError('unqualified normalized parent cell')
    response = client.request('parents', **fields)
    if response.pop('action', None) != 'parents':
        raise NativeProtocolError('native parent action changed')
    selected = []
    for role, column in (('father', 'PaternalID'), ('mother', 'MaternalID')):
        present = response.pop(role + '.present', None)
        if present == '0':
            selected.append(None)
        elif present == '1':
            raw = response.pop(role + '.position')
            position = int(raw)
            value = response.pop(role + '.value')
            if str(position) != raw or not 0 <= position < len(family) or family[column].iloc[position] != value:
                raise NativeProtocolError('native parent source position changed')
            selected.append(family[column].iloc[position])
        else:
            raise NativeProtocolError('native parent presence changed')
    if response:
        raise NativeProtocolError('extra native parent metadata')
    return selected


def _filter_record(client, record, ordinal, context, namespace):
    response = client.request('filter-begin', ordinal=str(ordinal))
    instance, source_gts = None, []
    while True:
        action = response.pop('action', None)
        if action in ('patient-gt', 'control-gt'):
            if set(response) != {'sample'}:
                raise NativeProtocolError('native sample action fields changed')
            sample = response['sample']
            # The original actual C-extension membership may itself raise.
            if sample in record.samples:
                fields = dict(present='1', **_encode_gt(record.samples[sample]['GT']))
            else:
                fields = dict(present='0')
        elif action == 'constructor-chrom':
            if response:
                raise NativeProtocolError('native constructor fields changed')
            instance = VariantSource(record, *context, session=client, namespace=namespace)
            fields = dict(ack='1')
        elif action == 'source-gt':
            if set(response) != {'sample'} or instance is None:
                raise NativeProtocolError('native source getter fields changed')
            # Empty string and None both skip; the context preserves the actual
            # original value for logs and direct facades.
            sample = context[len(source_gts)]
            if response['sample'] != (sample or ''):
                raise NativeProtocolError('native source getter identity changed')
            value = instance._get_genotype(sample)
            source_gts.append(value)
            fields = _encode_gt(value)
        elif action == 'debug-autosomal':
            if response or instance is None or len(source_gts) != 3:
                raise NativeProtocolError('native autosomal debug state changed')
            instance._debug_autosomal(*source_gts)
            fields = dict(ack='1')
        elif action == 'debug-no-patient':
            if response:
                raise NativeProtocolError('native patient rejection fields changed')
            message = f'No patient ({context.patients}) carries the variant {record}'
            _source_log(namespace, logging.DEBUG, message, 'pedigree_filter', 251)
            fields = dict(ack='1')
        elif action == 'debug-control-reject':
            if set(response) != {'sample'}:
                raise NativeProtocolError('native control rejection fields changed')
            message = f'Control sample {response["sample"]} is homozygous or hemizygous for the variant {record}'
            _source_log(namespace, logging.DEBUG, message, 'pedigree_filter', 260)
            fields = dict(ack='1')
        elif action == 'drop':
            if response:
                raise NativeProtocolError('native drop fields changed')
            return None
        elif action == 'write':
            if set(response) != {'value'} or response['value'] not in ('UNKNOWN', 'PATERNAL', 'MATERNAL', 'DENOVO', 'INHERITED'):
                raise NativeProtocolError('native source value changed')
            return response['value']
        else:
            raise NativeProtocolError('unknown native current-record action')
        response = client.request('filter-reply', **fields)


class _Context(tuple):
    def __new__(cls, proband, father, mother, sex, patients):
        result = super().__new__(cls, (proband, father, mother, sex))
        result.patients = patients
        return result


def pedigree_filter(vcf_path, ped_path, target_fam, output, *, binary=None, namespace=None, statistics=None, protocol=PROTOCOL):
    """Run the entire ordered native procedure, or admit the reference first."""
    namespace = _reference(namespace)
    try:
        if not binary or type(target_fam) not in (str, int, bool, float, type(None)):
            raise Unsupported('explicit binary/plain family target required')
        _admit(namespace)
    except Unsupported as error:
        if statistics is not None:
            statistics.update(backend='python', native_processes=0, fallback_reason=str(error))
        return namespace['pedigree_filter'](vcf_path, ped_path, target_fam, output)
    # Startup/source admission occurs before PED/log/output effects. A native
    # startup/transport error is explicit; no blanket catch hides it as Python.
    if protocol not in (PROTOCOL, TRANSITION_PROTOCOL):
        raise Unsupported('unqualified pedigree transport profile')
    client = NativeSession(binary) if protocol == PROTOCOL else TransitionSession(sys.modules[__name__], binary, statistics=statistics)
    with client:
        _admit_session(namespace, client, already_checked=True)
        if statistics is not None:
            statistics.update(backend='rust', native_processes=client.starts,
                              observed_native_operations=client.operations,
                              completed_native_operations=client.completed_operations,
                              binary_sha256=client.sha256, reference_fallback=False,
                              parsed_records=0, written_records=0)
        ped_df = pd.read_csv(ped_path, sep='\t', header=None, names=_COLUMNS,
                             usecols=range(6), dtype=str, keep_default_na=False)
        family, patients, controls, proband, sex, has_patient = _family(client, ped_df, target_fam)
        if family.empty:
            _source_log(namespace, logging.ERROR, f'Family {target_fam} not found in pedigree file.', 'pedigree_filter', 211)
            return None
        if not has_patient:
            raise IndexError('list index out of range')
        na_values = ['0', '.', '', 'NA', 'nan', 'NaN', '-', -9, 0, 0.0]
        family.loc[:, 'PaternalID'] = family['PaternalID'].replace(na_values, np.nan)
        family.loc[:, 'MaternalID'] = family['MaternalID'].replace(na_values, np.nan)
        _source_log(namespace, logging.INFO, f'The pedigree table looks like: \n{family.to_string(index=False)}', 'pedigree_filter', 225)
        _source_log(namespace, logging.INFO, f'The proband is {proband} and the proband sex is {sex}', 'pedigree_filter', 226)
        _source_log(namespace, logging.INFO, f'The patients are {patients}', 'pedigree_filter', 227)
        father, mother = _parents(client, family)
        fields = dict(proband=proband, father=father or '', mother=mother or '', sex=str(sex))
        if isinstance(client, TransitionSession):
            fields.update({'father.original_kind': 'none' if father is None else 'str',
                           'mother.original_kind': 'none' if mother is None else 'str',
                           'sex.original_kind': 'str'})
        _put_vector(fields, 'patients', list(patients))
        _put_vector(fields, 'controls', list(controls))
        if client.request('filter-init', **fields) != {'value': 'ready'}:
            raise NativeProtocolError('native complete-filter initialization changed')
        context = _Context(proband, father, mother, str(sex), patients)
        rows = kept = 0
        with pysam.VariantFile(vcf_path) as vcf_in:
            vcf_in.header.info.add('VARIANT_SOURCE', '1', 'String', 'Source of the variant')
            with pysam.VariantFile(output, 'w', header=vcf_in.header) as vcf_out:
                for rows, record in enumerate(vcf_in, 1):
                    if statistics is not None:
                        statistics['parsed_records'] = rows
                    source = _filter_record(client, record, rows, context, namespace)
                    if source is not None:
                        record.info['VARIANT_SOURCE'] = source
                        vcf_out.write(record)
                        kept += 1
                        if statistics is not None:
                            statistics['written_records'] = kept
        if statistics is not None:
            statistics.update(backend='rust', native_processes=client.starts,
                              protocol_requests=client.sequence, rows=rows, kept=kept,
                              binary_sha256=client.sha256, source_map='original log labels; actual adapter/native stack retained')
    _source_log(namespace, logging.INFO, f'Filtered VCF written to {output}', 'pedigree_filter', 278)
    if statistics is not None and isinstance(client, TransitionSession):
        statistics.update(client.metrics())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='procedure', required=True)
    whole = commands.add_parser('filter')
    for name in ('vcf', 'ped-table', 'target-family', 'output', 'binary'):
        whole.add_argument('--' + name, required=True)
    whole.add_argument('--protocol', choices=(PROTOCOL, TRANSITION_PROTOCOL), default=PROTOCOL)
    genotype = commands.add_parser('genotype')
    genotype.add_argument('--predicate', choices=('homo', 'carrier', 'exists', 'has'), required=True)
    genotype.add_argument('--allele', action='append', default=[], help='integer allele or . for missing; repeat for ploidy')
    genotype.add_argument('--binary', required=True)
    chrom = commands.add_parser('chromosome')
    chrom.add_argument('--chromosome', required=True)
    chrom.add_argument('--binary', required=True)
    source = commands.add_parser('source')
    source.add_argument('--method', choices=('auto', 'X', 'Y'), required=True)
    source.add_argument('--sex', required=True)
    source.add_argument('--proband-allele', action='append', default=[])
    source.add_argument('--father-allele', action='append', default=[])
    source.add_argument('--mother-allele', action='append', default=[])
    source.add_argument('--binary', required=True)
    args = parser.parse_args()
    if args.procedure == 'filter':
        return pedigree_filter(args.vcf, args.ped_table, args.target_family, args.output, binary=args.binary, protocol=args.protocol)
    with NativeSession(args.binary) as client:
        if args.procedure == 'genotype':
            value = tuple(None if value == '.' else int(value) for value in args.allele)
            print(_predicate(args.predicate, value, session=client))
        elif args.procedure == 'chromosome':
            print(_value(client.request('chromosome', chrom=args.chromosome)))
        else:
            fields = {}
            for name in ('proband', 'father', 'mother'):
                raw = getattr(args, name + '_allele')
                fields.update(_encode_gt(tuple(None if value == '.' else int(value) for value in raw) if raw else None, name))
            print(_value(client.request('source', method=args.method, sex=args.sex, **fields)))


if __name__ == '__main__':
    main()
