#!/usr/bin/env python3
"""Reuse a completed DAS build after proven ClinVar consumer-input equivalence.

The original completion receipt and scientific outputs remain untouched. A
separate, explicitly prepared proof must bind exact ordered consumer inputs,
reviewed publication and the unchanged analysis source. This checker verifies
that proof's integrity and bindings; it never infers equivalence from mtimes or
rescans a VCF. Without this optional receipt, normal DAS completion rules apply.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import sys
import tempfile

import das_install_state as state


KIND = 'das-clinvar-input-continuity-v1'
PROOF_KIND = 'das-clinvar-consumer-equivalence-v1'
MAX_METADATA_BYTES = 16 * 1024 * 1024


def require(ok, message):
    if not ok:
        raise ValueError(message)


def metadata(path):
    """Read bounded, stable regular metadata without following a link."""
    path = Path(path)
    before = path.lstat()
    require(stat.S_ISREG(before.st_mode) and 0 < before.st_size <= MAX_METADATA_BYTES,
            'Expected nonempty bounded regular metadata: ' + str(path))
    identity = lambda value: (value.st_dev, value.st_ino, value.st_size,
                               value.st_mtime_ns, value.st_ctime_ns, value.st_mode)
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
    with os.fdopen(descriptor, 'rb') as stream:
        require(identity(os.fstat(stream.fileno())) == identity(before),
                'Metadata changed before reading: ' + str(path))
        raw = stream.read(MAX_METADATA_BYTES + 1)
        require(len(raw) == before.st_size
                and identity(os.fstat(stream.fileno())) == identity(before),
                'Metadata changed during reading: ' + str(path))
    require(identity(path.lstat()) == identity(before),
            'Metadata path changed during reading: ' + str(path))
    return raw


def pin(path):
    path = Path(path).absolute()
    before = state.identities([path])[0]
    raw = metadata(path)
    require(state.identities([path])[0] == before, 'Metadata identity changed: ' + str(path))
    return {'path': str(path), 'sha256': hashlib.sha256(raw).hexdigest(), 'identity': before}


def pinned_bytes(value):
    require(isinstance(value, dict) and set(value) == {'path', 'sha256', 'identity'}
            and isinstance(value['path'], str) and Path(value['path']).is_absolute()
            and isinstance(value['sha256'], str)
            and re.fullmatch('[0-9a-f]{64}', value['sha256']), 'Invalid metadata proof pin')
    before = state.identities([value['path']])[0]
    raw = metadata(value['path'])
    require(before == value['identity'] == state.identities([value['path']])[0]
            and hashlib.sha256(raw).hexdigest() == value['sha256'],
            'Pinned metadata changed: ' + value['path'])
    return raw


def read_pin(value):
    return json.loads(pinned_bytes(value))


def proof(path, checksum, profile, old, new, sources):
    binding = pin(path)
    require(binding['sha256'] == checksum, 'Explicit consumer-equivalence proof changed')
    value = read_pin(binding)
    require(value.get('kind') == PROOF_KIND and type(value.get('schema_version')) is int
            and value['schema_version'] == 1
            and value.get('status') == 'completed' and value.get('profile') == profile,
            'Consumer-equivalence proof is incomplete or belongs to another profile')
    require(value.get('old_input') == old and value.get('new_input') == new,
            'Consumer-equivalence proof input identities differ')
    gates = ('full_record_census_equal', 'full_global_order_equal',
             'protected_annotation_fields_equal', 'headers_reviewed',
             'publication_identity_exact')
    require(isinstance(value.get('gates'), dict)
            and all(value['gates'].get(name) is True for name in gates),
            'Consumer-equivalence publication/header/row gates are incomplete')
    consumer = value.get('consumer', {})
    count = consumer.get('ordered_rows')
    digest = consumer.get('old_sha256')
    analysis = [item['sha256'] for item in sources if item['role'] == 'analysis']
    require(consumer.get('role') == 'das_and_dependent_pm1'
            and consumer.get('exact') is True and type(consumer.get('differing_rows')) is int
            and consumer['differing_rows'] == 0 and type(count) is int and count > 0
            and isinstance(digest, str) and re.fullmatch('[0-9a-f]{64}', digest)
            and consumer.get('new_sha256') == digest and len(analysis) == 1
            and consumer.get('analysis_source_sha256') == analysis[0],
            'Exact ordered DAS consumer/source binding is missing')
    evidence = value.get('evidence')
    require(isinstance(evidence, dict) and {'plan', 'comparison', 'publication'} <= set(evidence),
            'Underlying completed comparison/publication evidence is missing')
    for item in evidence.values():
        pinned_bytes(item)
    return binding, list(evidence.values())


def candidate(args, old_path, new_path, proof_path, proof_sha256):
    current = {'version': state.VERSION, **state.input_state(args),
               'outputs': state.identities(args.output)}
    original_pin = pin(args.receipt)
    original = read_pin(original_pin)
    old, new = state.identities([old_path, new_path])
    require(old != new, 'Continuity requires a distinct selected ClinVar input')
    positions = [i for i, item in enumerate(current['inputs']) if item == new]
    require(len(positions) == 1, 'Exactly one selected new ClinVar input is required')
    ordinal = positions[0]
    prior = dict(current, inputs=list(current['inputs']))
    prior['inputs'][ordinal] = old
    require(original == prior,
            'Original DAS source, parameters, bundled inputs, other data or outputs differ')
    protected = [args.receipt, *args.input, *args.output, proof_path]
    destination = Path(args.continuity).absolute()
    require(all(destination != Path(path).absolute() and
                (not destination.exists() or not os.path.samefile(destination, path))
                for path in protected), 'Continuity receipt aliases a protected file')
    proof_pin, evidence_pins = proof(proof_path, proof_sha256, args.profile,
                                     old, new, current['sources'])
    require({'version': state.VERSION, **state.input_state(args),
             'outputs': state.identities(args.output)} == current
            and state.identities([old_path, new_path]) == [old, new]
            and pin(args.receipt) == original_pin,
            'DAS inputs, outputs or original receipt changed during proof verification')
    # Treat the certificate and all underlying proof/source files as one pinned
    # snapshot. Rechecking only scientific dependencies would miss an evidence
    # change that occurred after its first scan, including during publication.
    for item in [proof_pin, *evidence_pins]:
        pinned_bytes(item)
    return {'kind': KIND, 'schema_version': 1, 'producer_execution': False,
            'checker_sha256': pin(__file__)['sha256'], 'original_receipt': original_pin,
            'old_input': old, 'new_input': new, 'input_ordinal': ordinal,
            'current_state': current, 'equivalence_proof': proof_pin}


def record(args):
    value = candidate(args, args.old_input, args.new_input,
                      args.proof, args.proof_sha256)
    raw = (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()
    path = Path(args.continuity)
    if os.path.lexists(path):
        require(metadata(path) == raw, 'Existing continuity differs; preserved')
        return
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + '.', dir=path.parent)
    try:
        with os.fdopen(descriptor, 'wb') as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        require(candidate(args, args.old_input, args.new_input, args.proof,
                          args.proof_sha256) == value,
                'DAS continuity changed before publication')
        try:
            os.link(temporary, path)
        except FileExistsError:
            require(metadata(path) == raw, 'Concurrent continuity differs; preserved')
    finally:
        os.unlink(temporary)


def check(args):
    value = json.loads(metadata(args.continuity))
    require(value.get('kind') == KIND, 'Unknown DAS continuity receipt')
    expected = candidate(args, value['old_input']['path'], value['new_input']['path'],
                         value['equivalence_proof']['path'], value['equivalence_proof']['sha256'])
    require(json.dumps(value, sort_keys=True, allow_nan=False)
            == json.dumps(expected, sort_keys=True, allow_nan=False),
            'DAS consumer-input continuity no longer matches')
    print('Reusing original DAS outputs through exact ClinVar consumer-input continuity')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for name, function in (('record', record), ('check', check)):
        command = commands.add_parser(name)
        for option in ('receipt', 'continuity', 'profile'):
            command.add_argument('--' + option, required=True)
        command.add_argument('--input', action='append', required=True)
        command.add_argument('--output', action='append', required=True)
        for option in ('source', 'bundled', 'parameter'):
            command.add_argument('--' + option, action='append', default=[], metavar='NAME=VALUE')
        if name == 'record':
            for option in ('old-input', 'new-input', 'proof', 'proof-sha256'):
                command.add_argument('--' + option, required=True)
        command.set_defaults(function=function)
    args = parser.parse_args()
    try:
        args.function(args)
    except (OSError, ValueError, TypeError, KeyError) as error:
        print('DAS input continuity: ' + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
