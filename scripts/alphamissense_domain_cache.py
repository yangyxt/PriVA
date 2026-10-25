#!/usr/bin/env python3
"""Checkpoint the unchanged AlphaMissense domain collector for installation.

2026-10-24 CT requirement: preserve the producer's first-wins/max/min rules,
CSQ/domain multiplicity and chromosome order while bounding collector count.
Raw chromosome dictionaries stay on disk instead of crossing Pool pipes. A
gene is finalized only after its last contributing chromosome, including PAR
or other nonadjacent contributions. This is software cache construction, not
scientific/model/family qualification. The scientific producer is not edited.

All destinations are beneath --cache-root. Source paths are provenance only;
generation identity uses source CONTENT hashes, so matching O/N snapshots reuse
the same completed generation. Input identity includes both available indexes
and exact nanoseconds as strings; no large input content hash is introduced.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.util
import json
from multiprocessing import Pool
import os
from pathlib import Path
import pickle
import resource
import stat
import sys
import time
import uuid

sys.dont_write_bytecode = True
HERE = Path(__file__).absolute().parent
PRODUCER = HERE / 'stat_protein_domain_amscores.py'
PRODUCER_SHA256 = 'cc7b2d60665e8052e6362107e359cc133e31c398f243589ee63ad7e8f90341c4'
SCHEMA = 'priva-alphamissense-domain-checkpoints-v1'


class InvalidCache(ValueError):
    """An owned completed checkpoint failed verification; preserve its history."""


def require(value, message):
    if not value:
        raise ValueError(message)


def utc():
    return datetime.now(timezone.utc).isoformat()


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n').encode()


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def identity(path, *, device=False):
    path = Path(path)
    value = path.stat()
    require(stat.S_ISREG(value.st_mode), 'Expected a regular file: ' + str(path))
    fields = ('st_ino', 'st_mode', 'st_size', 'st_mtime_ns', 'st_ctime_ns', 'st_nlink')
    if device:
        fields = ('st_dev', *fields)
    return {field: str(getattr(value, field)) for field in fields}


def inputs(path, *, device=False):
    path = Path(path).absolute()
    members = [path, *(Path(str(path) + suffix) for suffix in ('.tbi', '.csi')
                        if os.path.lexists(str(path) + suffix))]
    require(len(members) > 1, 'An indexed input VCF is required')
    return [dict(path=str(member), resolved=str(member.resolve(strict=True)),
                 identity=identity(member, device=device)) for member in members]


def source_state():
    result = {}
    for name, path in (('producer', PRODUCER), ('helper', Path(__file__).absolute())):
        before = identity(path, device=True)
        sha = digest(path)
        require(identity(path, device=True) == before, 'Source changed while reading: ' + str(path))
        result[name] = dict(path=str(path), identity=before, sha256=sha)
    require(result['producer']['sha256'] == PRODUCER_SHA256, 'Scientific producer differs from its pinned bytes')
    return result


def producer():
    require(digest(PRODUCER) == PRODUCER_SHA256, 'Scientific producer differs from its pinned bytes')
    spec = importlib.util.spec_from_file_location('_priva_domain_scientific_producer', PRODUCER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def managed(path, cache_root):
    path, root = Path(path).absolute(), Path(cache_root).absolute()
    require('..' not in path.parts and '..' not in root.parts and path.is_relative_to(root),
            'Destination leaves configured cache_root: ' + str(path))
    require(path.resolve().is_relative_to(root.resolve()), 'Destination aliases outside cache_root: ' + str(path))
    for item in (path, *path.parents):
        if item == root:
            break
        require(not item.is_symlink(), 'Preserve linked destination: ' + str(item))
    return path


def no_input_alias(path, input_members):
    path = Path(path)
    for member in input_members:
        source = Path(member['path'])
        require(path.resolve() != source.resolve() and
                not (path.exists() and os.path.samefile(path, source)),
                'Destination aliases selected input/index: ' + str(path))


@contextmanager
def lock(folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / '.build.lock').open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def sync_directory(folder):
    descriptor = os.open(folder, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_json(path, value):
    path = Path(path)
    partial = path.with_name(path.name + '.partial-' + uuid.uuid4().hex)
    with partial.open('xb') as stream:
        stream.write(encoded(value)); stream.flush(); os.fsync(stream.fileno())
    os.replace(partial, path)
    sync_directory(path.parent)


class HashWriter:
    def __init__(self, stream):
        self.stream, self.hash = stream, hashlib.sha256()

    def write(self, value):
        written = self.stream.write(value)
        # 2026-10-24 CT: protocol 5 may pass a PickleBuffer for a large NumPy
        # array. It has no len(); typed buffers also count elements, not bytes.
        with memoryview(value) as view:
            require(written == view.nbytes, 'Short checkpoint write')
        self.hash.update(value)
        return written


def write_pickle(path, value):
    """Retain an interrupted partial; expose its fresh name only after close."""
    path = Path(path)
    require(not os.path.lexists(path), 'Never overwrite a retained attempt: ' + str(path))
    partial = path.with_name(path.name + '.partial-' + uuid.uuid4().hex)
    with partial.open('xb') as stream:
        writer = HashWriter(stream)
        pickle.dump(value, writer, protocol=pickle.HIGHEST_PROTOCOL)
        stream.flush(); os.fsync(stream.fileno())
    os.replace(partial, path); sync_directory(path.parent)
    return dict(path=str(path), sha256=writer.hash.hexdigest(), identity=identity(path))


def file_entry(path):
    path = Path(path)
    before = identity(path)
    sha = digest(path)
    require(identity(path) == before, 'File changed while checksumming: ' + str(path))
    return dict(path=str(path), sha256=sha, identity=before)


def checked_file(entry, plan):
    path = managed(entry['path'], plan['cache_root'])
    try:
        require(identity(path) == entry['identity'] and digest(path) == entry['sha256'] and
                identity(path) == entry['identity'], 'Completed bytes/identity differ: ' + str(path))
    except (OSError, ValueError) as error:
        raise InvalidCache(str(error)) from error
    return path


def read_json(path):
    require(Path(path).stat().st_size <= 16 * 1024 * 1024, 'Oversized cache metadata: ' + str(path))
    with Path(path).open('rb') as stream:
        return json.load(stream)


def paths(plan):
    output = managed(plan['output_pickle'], plan['cache_root'])
    base = managed(str(output) + '.chromosomes', plan['cache_root'])
    root = managed(base / plan['generation'], plan['cache_root'])
    return output, base, root


def generation(recipe):
    return hashlib.sha256(encoded(recipe)).hexdigest()


def create_plan(*, input_vcf, output_pickle, cache_root):
    before_sources = source_state()
    before_inputs = inputs(input_vcf, device=True)
    current = inputs(input_vcf)
    output = managed(output_pickle, cache_root)
    base = managed(str(output) + '.chromosomes', cache_root)
    for path in (output, Path(str(output) + '.complete.json'), Path(str(output) + '.publish.pending.json'), base):
        no_input_alias(path, current)
    require(not any(Path(item['resolved']).is_relative_to(base.resolve()) for item in current),
            'Selected input/index lies in the checkpoint namespace')
    module = producer()
    with module.pysam.VariantFile(str(input_vcf)) as vcf:
        contigs = list(vcf.header.contigs)  # literal source order; no canonical-contig filter
    require(inputs(input_vcf, device=True) == before_inputs and source_state() == before_sources,
            'Input/index/source changed while planning')
    recipe = dict(schema=SCHEMA, inputs=current, contigs=contigs,
                  source_sha256={key: value['sha256'] for key, value in before_sources.items()})
    plan = dict(recipe=recipe, generation=generation(recipe), cache_root=str(Path(cache_root).absolute()),
                output_pickle=str(output), input_vcf=str(Path(input_vcf).absolute()),
                source_provenance=before_sources, created_utc=utc(), installation_accepted=False)
    _, _, root = paths(plan)
    publication_owned(plan)
    with lock(root):
        plan_path = root / 'plan.json'
        if plan_path.exists():
            old = read_json(plan_path)
            require(old['recipe'] == recipe and old['generation'] == plan['generation'] and
                    old['output_pickle'] == str(output), 'Existing plan differs; preserve it')
        else:
            atomic_json(plan_path, plan)
    return plan_path


def load_plan(path):
    path = Path(path).absolute()
    plan = read_json(path)
    recipe = plan['recipe']
    require(recipe['schema'] == SCHEMA and generation(recipe) == plan['generation'] and
            recipe['source_sha256'] == {key: value['sha256'] for key, value in source_state().items()} and
            recipe['inputs'] == inputs(plan['input_vcf']), 'Plan input/index/source generation changed')
    require(isinstance(recipe['contigs'], list) and all(isinstance(value, str) for value in recipe['contigs']) and
            len(recipe['contigs']) == len(set(recipe['contigs'])), 'Malformed literal contig inventory')
    output, _, root = paths(plan)
    require(path == root / 'plan.json', 'Plan is outside its generation')
    for member in (output, Path(str(output) + '.complete.json'), Path(str(output) + '.publish.pending.json')):
        no_input_alias(member, recipe['inputs'])
    return plan


def chromosome_paths(plan, ordinal):
    require(type(ordinal) is int and 0 <= ordinal < len(plan['recipe']['contigs']), 'Unknown chromosome ordinal')
    folder = managed(paths(plan)[2] / ('chromosome-%05d' % ordinal), plan['cache_root'])
    managed(folder / 'complete.json', plan['cache_root'])
    no_input_alias(folder / 'complete.json', plan['recipe']['inputs'])
    return folder, folder / 'complete.json'


def checkpoint(plan, ordinal, *, verify_raw=True):
    folder, ready = chromosome_paths(plan, ordinal)
    try:
        value = read_json(ready)
        require(value['generation'] == plan['generation'] and value['ordinal'] == ordinal and
                value['chromosome'] == plan['recipe']['contigs'][ordinal] and value['complete'] is True,
                'Checkpoint generation/chromosome is incomplete or different')
        raw = managed(value['raw']['path'], plan['cache_root'])
        require(identity(raw) == value['raw']['identity'], 'Raw checkpoint identity differs')
        if verify_raw:
            checked_file(value['raw'], plan)
        genes = read_json(checked_file(value['genes'], plan))
        require(raw.parent == folder and Path(value['genes']['path']).parent == folder and
                isinstance(genes, list) and all(isinstance(gene, str) for gene in genes) and len(genes) == len(set(genes)),
                'Malformed checkpoint gene inventory/path')
        return value, genes
    except (OSError, KeyError, TypeError, json.JSONDecodeError, InvalidCache, ValueError) as error:
        raise InvalidCache(str(error)) from error


def collect_chromosome(*, plan_path, ordinal):
    plan = load_plan(plan_path)
    publication_owned(plan)
    folder, ready = chromosome_paths(plan, ordinal)
    with lock(folder):
        if ready.exists():
            try:
                observed_inputs, observed_sources = inputs(plan['input_vcf'], device=True), source_state()
                value, genes = checkpoint(plan, ordinal)
                require(inputs(plan['input_vcf'], device=True) == observed_inputs and source_state() == observed_sources,
                        'Input/index/source changed during checkpoint reuse verification')
                return dict(receipt=str(ready), sha256=digest(ready), genes=len(genes), reused=True)
            except InvalidCache:
                # Retain the invalid receipt AND every raw/partial member it names.
                os.replace(ready, ready.with_name('complete.invalid-' + uuid.uuid4().hex + '.json'))
        before_inputs = inputs(plan['input_vcf'], device=True)
        before_sources = source_state()
        module = producer()
        collector = module.DomainAMScoreCollector(plan['input_vcf'])
        try:
            collector.collect_scores(region=plan['recipe']['contigs'][ordinal])
        finally:
            require(inputs(plan['input_vcf'], device=True) == before_inputs and source_state() == before_sources,
                    'Input/index/source changed during original collection')
        genes = list(collector.scores_dict)
        # Convert/release only one raw gene tree at once, never a whole-chromosome clone.
        for gene in genes:
            collector.scores_dict[gene] = module.to_plain_dict(collector.scores_dict[gene])
        attempt = uuid.uuid4().hex
        raw = write_pickle(folder / (attempt + '.raw.pkl'), dict(collector.scores_dict))
        del collector
        gene_path = folder / (attempt + '.genes.json')
        atomic_json(gene_path, genes)
        require(inputs(plan['input_vcf'], device=True) == before_inputs and source_state() == before_sources,
                'Input/index/source changed before checkpoint completion')
        value = dict(generation=plan['generation'], ordinal=ordinal, chromosome=plan['recipe']['contigs'][ordinal],
                     raw=raw, genes=file_entry(gene_path), source_provenance=before_sources,
                     complete=True, finished_utc=utc(), installation_accepted=False)
        atomic_json(ready, value)
        return dict(receipt=str(ready), sha256=digest(ready), genes=len(genes), reused=False)


def checked_receipt_refs(plan, refs):
    """Bind every small completion receipt; never reread raw chromosome models."""
    contigs = plan['recipe']['contigs']
    require(isinstance(refs, list) and len(refs) == len(contigs), 'Missing chromosome completion references')
    for ordinal, (chromosome, entry) in enumerate(zip(contigs, refs)):
        expected = chromosome_paths(plan, ordinal)[1]
        require(isinstance(entry, dict) and set(entry) == {'path', 'sha256'} and entry['path'] == str(expected),
                'Chromosome receipt references are not in exact header order')
        before = identity(expected, device=True)
        require(0 < int(before['st_size']) <= 16 * 1024 * 1024 and digest(expected) == entry['sha256'],
                'Chromosome completion receipt checksum differs')
        value = read_json(expected)
        require(identity(expected, device=True) == before and value.get('generation') == plan['generation'] and
                type(value.get('ordinal')) is int and value['ordinal'] == ordinal and
                value.get('chromosome') == chromosome and value.get('complete') is True,
                'Chromosome receipt generation/ordinal/contig/completion differs')


def output_ready(plan):
    output = paths(plan)[0]
    pending, ready = Path(str(output) + '.publish.pending.json'), Path(str(output) + '.complete.json')
    if pending.exists() or not ready.exists():
        return None
    try:
        value = read_json(ready)
        require(value['schema'] == SCHEMA and value['generation'] == plan['generation'] and value['complete'] is True and
                value['output']['path'] == str(output) and value['source_sha256'] == plan['recipe']['source_sha256'] and
                value['plan']['path'] == str(paths(plan)[2] / 'plan.json'), 'Final generation differs')
        checked_file(value['plan'], plan)
        checked_receipt_refs(plan, value.get('checkpoint_receipts'))
        checked_file(value['output'], plan)
        return dict(output=str(output), receipt=str(ready), sha256=value['output']['sha256'], reused=True)
    except (OSError, KeyError, TypeError, ValueError):
        return None


def publication_owned(plan):
    """Containment is not ownership of an existing user-supplied payload."""
    output = paths(plan)[0]
    records = []
    for suffix in ('.complete.json', '.publish.pending.json'):
        path = managed(str(output) + suffix, plan['cache_root'])
        no_input_alias(path, plan['recipe']['inputs'])
        if path.exists():
            value = read_json(path)
            require(value.get('schema') == SCHEMA and value.get('output_path') == str(output),
                    'Preserve unowned/unknown publication sidecar: ' + str(path))
            records.append(value)
    if output.exists():
        candidates = [record.get('output', record.get('new', {})) for record in records]
        require(candidates, 'Preserve unowned canonical output with no helper receipt: ' + str(output))
        before = identity(output)
        sha = digest(output)
        require(identity(output) == before and any(entry.get('sha256') == sha and
                entry.get('identity', {}).get('st_size') == before['st_size'] for entry in candidates),
                'Preserve unowned or changed canonical output; no helper completion proof: ' + str(output))


def publish(plan, built, checkpoint_receipts):
    """One atomic payload rename; pending gates reuse until its receipt closes.

    Existing payload/receipt/pending files are retained by exact name, never
    swept. A failed publication remains pending; a retry rebuilds from verified
    raw checkpoints and retains that interrupted publication before replacing it.
    """
    output, base, _ = paths(plan)
    ready, pending = Path(str(output) + '.complete.json'), Path(str(output) + '.publish.pending.json')
    for path in (output, ready, pending):
        managed(path, plan['cache_root']); no_input_alias(path, plan['recipe']['inputs'])
    publication_owned(plan)
    checked_receipt_refs(plan, checkpoint_receipts)
    require(inputs(plan['input_vcf']) == plan['recipe']['inputs'] and
            {key: value['sha256'] for key, value in source_state().items()} == plan['recipe']['source_sha256'],
            'Input/index/source changed before publication')
    retained = base / ('retained-' + uuid.uuid4().hex)
    retained.mkdir()
    history = []
    for path in (output, ready, pending):
        if os.path.lexists(path):
            require(path.is_file() and not path.is_symlink(), 'Preserve nonregular publication member: ' + str(path))
            history.append(dict(path=str(path), retained=str(retained / path.name), identity=identity(path)))
    # Record old identities before any rename, including interrupted history.
    atomic_json(retained / 'publication.json', dict(generation=plan['generation'], prior=history, new=built))
    for entry in history:
        os.replace(entry['path'], entry['retained'])
    atomic_json(pending, dict(schema=SCHEMA, output_path=str(output), generation=plan['generation'], new=built, retained=str(retained)))
    require(inputs(plan['input_vcf']) == plan['recipe']['inputs'], 'Input/index changed before publication')
    os.replace(built['path'], output)
    value = dict(schema=SCHEMA, output_path=str(output), generation=plan['generation'],
                 output=dict(path=str(output), sha256=built['sha256'], identity=identity(output)),
                 plan=file_entry(paths(plan)[2] / 'plan.json'), source_sha256=plan['recipe']['source_sha256'],
                 checkpoint_receipts=checkpoint_receipts, complete=True, finished_utc=utc(), installation_accepted=False)
    atomic_json(ready, value)
    require(inputs(plan['input_vcf']) == plan['recipe']['inputs'] and
            {key: value['sha256'] for key, value in source_state().items()} == plan['recipe']['source_sha256'],
            'Input/index/source changed while publishing; retain pending marker')
    pending.unlink(); sync_directory(output.parent)
    return dict(output=str(output), receipt=str(ready), sha256=built['sha256'], reused=False)


def merge(*, plan_path):
    plan = load_plan(plan_path)
    output, base, root = paths(plan)
    with lock(base):
        publication_owned(plan)  # fail before collecting/moving any unowned canonical
        reused = output_ready(plan)
        if reused:
            return reused
        before_inputs, before_sources = inputs(plan['input_vcf'], device=True), source_state()
        inventories, receipts, first, last = [], [], {}, {}
        for ordinal in range(len(plan['recipe']['contigs'])):
            value, genes = checkpoint(plan, ordinal, verify_raw=False)  # no missing-checkpoint fallback
            inventories.append((value, genes))
            receipts.append(dict(path=str(chromosome_paths(plan, ordinal)[1]), sha256=digest(chromosome_paths(plan, ordinal)[1])))
            for gene in genes:
                first.setdefault(gene, None); last[gene] = ordinal
        module = producer()
        final, pending = dict.fromkeys(first), {}
        for ordinal, (value, genes) in enumerate(inventories):
            raw_path = checked_file(value['raw'], plan)
            before = identity(raw_path)
            with raw_path.open('rb') as stream:
                part = pickle.load(stream)
            require(identity(raw_path) == before and type(part) is dict and list(part) == genes,
                    'Raw checkpoint changed or disagrees with its ordered gene inventory')
            for gene in list(part):
                module._merge_scores(pending, {gene: part.pop(gene)})
                if last[gene] == ordinal:
                    completed = {gene: pending.pop(gene)}
                    module._finalize_scores_dict(completed)
                    final[gene] = module.to_plain_dict(completed[gene])
                    del completed
            del part
        require(not pending and all(value is not None for value in final.values()), 'A gene was not completely reduced')
        require(inputs(plan['input_vcf'], device=True) == before_inputs and source_state() == before_sources,
                'Input/index/source changed during ordered merge')
        built = write_pickle(root / ('final-' + uuid.uuid4().hex + '.pkl'), final)
        del final
        require(inputs(plan['input_vcf'], device=True) == before_inputs and source_state() == before_sources,
                'Input/index/source changed before final publication')
        return publish(plan, built, receipts)


def collect_task(arguments):
    return collect_chromosome(plan_path=arguments[0], ordinal=arguments[1])


def build(*, input_vcf, output_pickle, cache_root, workers=1):
    require(type(workers) is int and workers > 0, '--workers must be a positive explicit process cap')
    cpu_limit = int(os.environ.get('SLURM_CPUS_PER_TASK') or (len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else os.cpu_count() or 1))
    require(workers <= cpu_limit, 'Collector count exceeds the actual CPU allocation')
    plan_path = create_plan(input_vcf=input_vcf, output_pickle=output_pickle, cache_root=cache_root)
    plan = load_plan(plan_path)
    reused = output_ready(plan)
    if reused:
        return reused
    tasks = [(str(plan_path), ordinal) for ordinal in range(len(plan['recipe']['contigs']))]
    if workers == 1:
        for task in tasks:
            collect_task(task)
    elif tasks:
        # Only small receipt metadata crosses IPC. Finish all collectors before
        # constructing the final model; recycled children never inherit it.
        with Pool(min(workers, len(tasks)), maxtasksperchild=1) as pool:
            list(pool.imap(collect_task, tasks, chunksize=1))
    return merge(plan_path=plan_path)


def assert_equal(expected, actual, path='model'):
    """Exact logical comparator, including insertion order and floating bits."""
    import numpy as np
    import struct
    require(type(expected) is type(actual), 'Type differs at ' + path)
    if isinstance(expected, np.ndarray):
        require(expected.dtype == actual.dtype and expected.shape == actual.shape and
                expected.tobytes() == actual.tobytes(), 'Array differs at ' + path)
    elif isinstance(expected, dict):
        require(list(expected) == list(actual), 'Ordered keys differ at ' + path)
        for key in expected:
            assert_equal(expected[key], actual[key], path + '/' + str(key))
    elif isinstance(expected, (list, tuple)):
        require(len(expected) == len(actual), 'Length differs at ' + path)
        for ordinal, (left, right) in enumerate(zip(expected, actual)):
            assert_equal(left, right, path + '/' + str(ordinal))
    elif isinstance(expected, float):
        require(struct.pack('>d', expected) == struct.pack('>d', actual), 'Floating bits differ at ' + path)
    else:
        require(expected == actual, 'Value differs at ' + path)


def compare(*, expected_pickle, actual_pickle):
    """Later full-model benchmark; explicit invocation reads both named models."""
    started = time.monotonic()
    with Path(expected_pickle).open('rb') as stream:
        expected = pickle.load(stream)
    with Path(actual_pickle).open('rb') as stream:
        actual = pickle.load(stream)
    assert_equal(expected, actual)
    return dict(equal=True, wall_seconds=time.monotonic() - started,
                maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                installation_accepted=False, scientific_qualification=False)


def compare_chromosome(*, plan_path, chromosome):
    """Independent legacy collection versus one completed raw chromosome cache.

    Explicit later benchmark only. It reads the entire named chromosome, keeps
    both small-chromosome models for exact comparison, and reports this process's
    high-water RSS, not a bound on a concurrent build or a full-model merger.
    """
    plan = load_plan(plan_path)
    ordinal = plan['recipe']['contigs'].index(chromosome)
    before_inputs, before_sources = inputs(plan['input_vcf'], device=True), source_state()
    module = producer()
    started = time.monotonic()
    legacy = {}
    module._merge_scores(legacy, module._collect_worker((plan['input_vcf'], chromosome)))
    module._finalize_scores_dict(legacy)
    legacy = module.to_plain_dict(legacy)
    legacy_seconds = time.monotonic() - started
    value, genes = checkpoint(plan, ordinal)
    started = time.monotonic()
    with Path(value['raw']['path']).open('rb') as stream:
        raw = pickle.load(stream)
    require(list(raw) == genes, 'Checkpoint gene inventory differs')
    module._finalize_scores_dict(raw)
    actual = module.to_plain_dict(raw)
    assert_equal(legacy, actual)
    require(inputs(plan['input_vcf'], device=True) == before_inputs and source_state() == before_sources,
            'Input/index/source changed during benchmark')
    return dict(chromosome=chromosome, equal=True, genes=len(genes), legacy_seconds=legacy_seconds,
                checkpoint_reduce_compare_seconds=time.monotonic() - started,
                maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                installation_accepted=False, scientific_qualification=False)


def main(argv=None):
    global PRODUCER
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='procedure', required=True)
    for name in ('plan', 'build'):
        child = commands.add_parser(name)
        for key in ('input-vcf', 'output-pickle', 'cache-root'):
            child.add_argument('--' + key, type=Path, required=True)
        if name == 'build':
            child.add_argument('--workers', type=int, default=1)
    child = commands.add_parser('collect-chromosome'); child.add_argument('--plan', type=Path, required=True)
    selection = child.add_mutually_exclusive_group(required=True)
    selection.add_argument('--ordinal', type=int); selection.add_argument('--chromosome')
    for name in ('merge', 'status', 'check'):
        child = commands.add_parser(name); child.add_argument('--plan', type=Path, required=True)
    child = commands.add_parser('compare')
    child.add_argument('--expected-pickle', type=Path, required=True); child.add_argument('--actual-pickle', type=Path, required=True)
    child = commands.add_parser('compare-chromosome')
    child.add_argument('--plan', type=Path, required=True); child.add_argument('--chromosome', required=True)
    for child in commands.choices.values():
        child.add_argument('--source-script', type=Path, default=PRODUCER,
                           help='Exact unchanged scientific producer; only its content SHA enters generation identity')
    args = parser.parse_args(argv)
    PRODUCER = args.source_script.absolute()
    if args.procedure in {'plan', 'build'}:
        values = dict(input_vcf=args.input_vcf, output_pickle=args.output_pickle, cache_root=args.cache_root)
        value = build(**values, workers=args.workers) if args.procedure == 'build' else dict(plan=str(create_plan(**values)))
    elif args.procedure == 'collect-chromosome':
        plan = load_plan(args.plan)
        ordinal = args.ordinal if args.ordinal is not None else plan['recipe']['contigs'].index(args.chromosome)
        value = collect_chromosome(plan_path=args.plan, ordinal=ordinal)
    elif args.procedure == 'merge':
        value = merge(plan_path=args.plan)
    elif args.procedure == 'status':
        plan = load_plan(args.plan)
        value = dict(generation=plan['generation'], contigs=plan['recipe']['contigs'],
                     receipt_present=[chromosome_paths(plan, i)[1].exists() for i in range(len(plan['recipe']['contigs']))],
                     content_verified=False, installation_accepted=False)
    elif args.procedure == 'check':
        value = output_ready(load_plan(args.plan))
        require(value is not None, 'No verified completed domain model for this plan')
        value.update(content_verified=True, installation_accepted=False, scientific_qualification=False)
    elif args.procedure == 'compare':
        value = compare(expected_pickle=args.expected_pickle, actual_pickle=args.actual_pickle)
    else:
        value = compare_chromosome(plan_path=args.plan, chromosome=args.chromosome)
    print(json.dumps(value, sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
