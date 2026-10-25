#!/usr/bin/env python3
"""Guarded shell-workflow bridge to the installed original-Python CADD stages.

Installation never qualifies this route. A selected qualified profile must bind
both native single-VEP selection and the original publisher-prescore policy.
The original CADD prepare/extract/sort commands remain the scientific authority;
this module proves their alleles can reuse the already emitted PriVA CSQs.
Changed/ambiguous allele representations need qualification, never another VEP.
"""
import argparse
import gzip
import hashlib
import importlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
import yaml

ROOT = Path(__file__).resolve().parents[1]
POLICY = 'upstream-v1.7.2-inclAnno-v1'
CHECKS = {'prepared_allele_mapping': 'exact', 'publisher_import': 'exact',
          'novel_feature_joins': 'exact', 'source_order': 'preserved', 'missingness': 'preserved'}


class WorkflowError(ValueError):
    pass


def require(ok, message):
    if not ok:
        raise WorkflowError(message)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def adapter_identity():
    names = ('scripts/cadd_shared_workflow.py', 'scripts/cadd_native.py', 'scripts/annotation_vcf.sh',
             'scripts/anno_main_args', 'scripts/anno_vep_args', 'scripts/common_bash_utils.sh',
             'bin/cadd_task.py', 'bin/cadd_shared_vep.py', 'bin/cadd_select_view.py',
             'bin/cadd/CaddSelection.pm', 'bin/cadd/upstream_annotate.sh', 'bin/cadd/upstream_calculate.sh',
             'conf/cadd/GRCh37.reference.json', 'conf/cadd/GRCh38.reference.json')
    return {name: digest(ROOT / name) for name in names}


def read_json(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write('\n')


def text_file(path):
    with Path(path).open('rb') as handle:
        compressed = handle.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if compressed else Path(path).open()


def vcf(path):
    headers = []
    with text_file(path) as handle:
        for line in handle:
            if line.startswith('#'):
                headers.append(line.rstrip('\r\n'))
            else:
                break
    require(any(line.startswith('#CHROM\t') for line in headers), 'Missing VCF column header')
    return headers, vcf_rows(path)


def vcf_rows(path):
    with text_file(path) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            require(len(fields) >= 5, 'Incomplete VCF record')
            yield fields


def key(row):
    return tuple(row[i] for i in (0, 1, 3, 4))


def csq(row):
    require(len(row) >= 8, 'Shared VEP source needs INFO and complete CSQs')
    values = [item[4:] for item in row[7].split(';') if item.startswith('CSQ=')]
    require(len(values) == 1 and values[0], 'Every shared source occurrence needs one complete CSQ payload')
    return values[0]


def native_modules():
    sys.path.insert(0, str(ROOT / 'bin'))
    task = importlib.import_module('cadd_task')
    shared = importlib.import_module('cadd_shared_vep')
    view = importlib.import_module('cadd_select_view')
    spec = importlib.util.spec_from_file_location('priva_cadd_native_bridge', ROOT / 'scripts/cadd_native.py')
    native = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(native)
    return task, shared, view, native


def selection(args, config):
    selected = args.profile or config.get('cadd_profile')
    if not selected:
        candidates = Path(config.get('cache_root', '<cache_root>')) / 'runtime/cadd-1.7.2/profiles'
        raise WorkflowError('Select an explicitly qualified CADD profile with --cadd-profile or config cadd_profile; '
                            'installed candidates: ' + str(candidates))
    path = Path(selected)
    require(path.is_absolute() and path.is_file(), 'CADD profile is missing; install resources, then explicitly qualify the runtime profile')
    wanted = config.get('cadd_profile_sha256')
    if wanted:
        require(digest(path) == wanted, 'Configured CADD profile checksum differs')
    return path


def file_record(task, path, expected):
    require(task.record(path) == {name: value for name, value in expected.items() if name != 'sha256'},
            'Publisher prescore resource identity changed: ' + str(path))
    if 'sha256' in expected:
        require(digest(path) == expected['sha256'], 'Publisher prescore resource checksum changed: ' + str(path))


def check_prescoring(profile, task):
    block = profile.get('prescoring', {})
    require(block.get('policy_id') == POLICY, 'Profile lacks the original incl_anno prescore policy')
    for name in ('configuration', 'prepare', 'extract', 'snakefile'):
        item = block.get(name, {})
        require(item.get('path') and Path(item['path']).is_absolute() and digest(item['path']) == item.get('sha256'),
                'Original prescore source pin differs: ' + name)
    config = yaml.safe_load(Path(block['configuration']['path']).read_text())
    require(config.get('Annotation') is True and isinstance(config.get('Header'), str),
            'Prescore route requires the original annotation-including scientific configuration')
    source = Path(profile['annotation']['env']['CADD_SOURCE_ROOT'])
    require((source / config['PrescoredFolder']).resolve() == Path(block['folder']).resolve(),
            'Prescore folder differs from the original selected configuration')
    entries = block.get('files')
    require(isinstance(entries, list), 'Profile must declare the ordered publisher prescore inputs')
    actual = sorted(Path(block['folder']).glob('*.tsv.gz'), key=lambda path: os.fsencode(str(path)))
    require([str(path) for path in actual] == [entry['path'] for entry in entries],
            'Ordered publisher prescore membership changed')
    for entry in entries:
        file_record(task, entry['path'], entry['data'])
        file_record(task, str(entry['path']) + '.tbi', entry['index'])
    qualification = block.get('qualification')
    require(qualification and digest(qualification['path']) == qualification.get('sha256'),
            'Publisher/imported and novel CADD routes are not jointly qualified')
    evidence = read_json(qualification['path'])
    expected = {'schema_version': 1, 'status': 'qualified', 'policy_id': POLICY,
                'profile_id': profile['profile_id'], 'resource_manifest_sha256': profile['contract']['resource_manifest_sha256'],
                'stub': False, 'checks': CHECKS}
    require(all(evidence.get(name) == value for name, value in expected.items()),
            'Prescore-route qualification does not cover this exact profile and policy')
    proof = evidence.get('evidence', {})
    require(proof.get('path') and digest(proof['path']) == proof.get('sha256'),
            'Prescore-route qualification evidence is missing or changed')
    return block, config


def preflight(args):
    config = yaml.safe_load(Path(args.config).read_text())
    require(isinstance(config, dict), 'CONFIG must be a YAML mapping')
    profile_path = selection(args, config)
    task, shared, view, native = native_modules()
    profile = native.preflight(config, profile_path)
    block, original = check_prescoring(profile, task)
    require(Path(config.get('tmp_dir', '')).is_absolute() and Path(config['tmp_dir']).is_dir(),
            'Configured tmp_dir must be an existing absolute directory')
    return config, profile_path, profile, block, original, (task, shared, view, native)


def native_rows(source, profile, modules, temporary):
    """Validate every original CSQ ordinal without altering the source artifact."""
    _, shared, view, _ = modules
    context = shared.native_context(profile)
    view.select_rows(source, temporary, context)
    return context


class OccurrenceMap:
    """Bounded SQLite index; payloads and occurrence proofs stay on disk.

    Original chromosome renaming is proved with complete-record SHA256 and
    occurrence counts. Prepared sort/uniq may reorder/collapse only exact allele
    groups with identical complete CSQs. No prepared allele trim is guessed.
    """
    def __init__(self, original, renamed, prepared, contig_mapping, work, proof):
        self.connection = None
        self.previous_sqlite_tmp = os.environ.get('SQLITE_TMPDIR')
        os.environ['SQLITE_TMPDIR'] = str(work)
        self.connection = sqlite3.connect(str(Path(work) / 'occurrences.sqlite'))
        db = self.connection
        db.execute('PRAGMA cache_size=-2048')
        db.execute('PRAGMA temp_store=FILE')
        # Explicitly direct SQLite sort/index spills too; this is a dedicated
        # short-lived CLI connection, never a shared application database.
        db.execute("PRAGMA temp_store_directory='" + str(work).replace("'", "''") + "'")
        db.execute('CREATE TABLE source (ordinal INTEGER PRIMARY KEY, allele TEXT, csq TEXT, record_sha TEXT, used INTEGER DEFAULT 0, renamed_ordinal INTEGER)')
        db.execute('CREATE INDEX source_records ON source(record_sha,used,ordinal)')
        db.execute('CREATE INDEX source_alleles ON source(allele,renamed_ordinal)')
        db.execute('CREATE TABLE prepared (ordinal INTEGER PRIMARY KEY, allele TEXT UNIQUE, source_ordinal INTEGER)')
        try:
            with db, Path(proof).open('x') as mapping_output:
                for ordinal, row in enumerate(vcf_rows(original)):
                    converted = row[:]
                    converted[0] = contig_mapping.get(row[0], row[0])
                    db.execute('INSERT INTO source(ordinal,allele,csq,record_sha) VALUES(?,?,?,?)',
                        (ordinal, '\t'.join(key(converted)), csq(row), hashlib.sha256('\t'.join(converted).encode()).hexdigest()))
                for ordinal, row in enumerate(vcf_rows(renamed)):
                    record_sha = hashlib.sha256('\t'.join(row).encode()).hexdigest()
                    found = db.execute('SELECT ordinal FROM source WHERE record_sha=? AND used=0 ORDER BY ordinal LIMIT 1', (record_sha,)).fetchone()
                    require(found is not None, 'Original chromosome renaming changed records, source CSQs, genotypes or occurrence counts')
                    db.execute('UPDATE source SET used=1,renamed_ordinal=? WHERE ordinal=?', (ordinal, found[0]))
                    mapping_output.write(json.dumps({'stage':'chromosome_rename', 'renamed_occurrence':ordinal,
                                                     'source_occurrence':found[0]}) + '\n')
                require(db.execute('SELECT 1 FROM source WHERE used=0 LIMIT 1').fetchone() is None,
                        'Original chromosome renaming lost source occurrences')
                for ordinal, row in enumerate(vcf_rows(prepared)):
                    allele = '\t'.join(key(row))
                    found = db.execute('SELECT ordinal,csq,renamed_ordinal FROM source WHERE allele=? ORDER BY renamed_ordinal', (allele,))
                    first = found.fetchone()
                    require(first is not None, 'Original CADD preparation changed/dropped allele representation; shared-CSQ remapping is unqualified')
                    try:
                        db.execute('INSERT INTO prepared VALUES(?,?,?)', (ordinal, allele, first[0]))
                    except sqlite3.IntegrityError:
                        raise WorkflowError('Prepared allele set is not unique')
                    for occurrence in itertools.chain((first,), found):
                        require(occurrence[1] == first[1], 'Repeated allele has differing complete CSQs; shared source mapping is ambiguous')
                        mapping_output.write(json.dumps({'stage':'prepared_occurrence', 'prepared_occurrence':ordinal,
                            'allele':list(key(row)), 'source_occurrence':occurrence[0], 'renamed_occurrence':occurrence[2]}) + '\n')
                require(db.execute('SELECT 1 FROM source s WHERE NOT EXISTS (SELECT 1 FROM prepared p WHERE p.allele=s.allele) LIMIT 1').fetchone() is None,
                        'Original CADD preparation changed/dropped allele representation; shared-CSQ remapping is unqualified')
        except BaseException:
            self.close()
            raise

    def shared_records(self, prepared):
        for row in vcf_rows(prepared):
            found = self.connection.execute('SELECT s.csq FROM prepared p JOIN source s ON s.ordinal=p.source_ordinal WHERE p.allele=?',
                                            ('\t'.join(key(row)),)).fetchone()
            require(found is not None, 'Novel allele is outside the proven prepared/shared source map')
            yield [row[0], row[1], '.', row[3], row[4], '.', '.', 'CSQ=' + found[0]]

    def close(self):
        if self.connection is not None:
            self.connection.close()
            self.connection = None
        if self.previous_sqlite_tmp is None:
            os.environ.pop('SQLITE_TMPDIR', None)
        else:
            os.environ['SQLITE_TMPDIR'] = self.previous_sqlite_tmp

    def __enter__(self):
        return self

    def __exit__(self, *unused):
        self.close()


def write_shared(path, headers, records):
    # Upstream ESM indel readers assume the first INFO declaration is CSQ and
    # split row INFO by commas. Inherited AF/model tags and samples must not leak
    # into CADD. All source CSQs are retained unchanged in this private artifact.
    schemas = [line for line in headers if line.startswith('##INFO=<ID=CSQ,')]
    require(len(schemas) == 1, 'Need exactly one original shared CSQ schema')
    retained = ['##fileformat=VCFv4.2'] + [line for line in headers if line.startswith('##contig=')]
    with gzip.open(path, 'xt') as handle:
        handle.write('\n'.join(retained + schemas + ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO']) + '\n')
        for row in records:
            handle.write('\t'.join(row) + '\n')


def table(path):
    header = None
    with text_file(path) as handle:
        for line in handle:
            line = line.rstrip('\r\n')
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                header = line
    return header, score_rows(path)


def score_rows(path):
    seen = False
    with text_file(path) as handle:
        for line in handle:
            text = line.rstrip('\r\n')
            if text.startswith('##'):
                continue
            if text.startswith('#'):
                seen = True
            else:
                require(seen, 'Score rows precede their header')
                yield text


def score_key(line):
    fields = line.split('\t')
    require(len(fields) >= 4, 'Incomplete publisher score row')
    return tuple(fields[:4])


def prove_split(previous, remaining, found):
    remaining, found = iter(remaining), iter(found)
    next_remaining, next_found = next(remaining, None), next(found, None)
    for row in previous:
        allele = key(row)
        kept = next_remaining is not None and key(next_remaining) == allele
        imported = next_found is not None and score_key(next_found) == allele
        require(kept != imported, 'Publisher lookup changed allele representation, split membership or order, or returned an unrelated allele')
        if kept:
            next_remaining = next(remaining, None)
        else:
            while next_found is not None and score_key(next_found) == allele:
                next_found = next(found, None)
    require(next_remaining is None and next_found is None, 'Publisher lookup returned an unrelated allele or leftover split rows')


def join_scores(sources, output, environment):
    # Snakefile join selects the last non-## header, sorts all data with the
    # default whole-line tie fallback, and does NOT uniq publisher score rows.
    header = None
    with tempfile.TemporaryDirectory(prefix='cadd-join-', dir=environment['TMPDIR']) as temporary:
        data = Path(temporary) / 'rows.tsv'
        with data.open('w') as raw:
            for source in sources:
                with text_file(source) as handle:
                    for line in handle:
                        if line.startswith('##'):
                            continue
                        if line.startswith('#'):
                            header = line.rstrip('\r\n')
                        else:
                            raw.write(line)
        require(header is not None, 'Original CADD join has no score header')
        with Path(output).open('xb') as handle, data.open('rb') as incoming:
            handle.write((header + '\n').encode()); handle.flush()
            subprocess.run(['sort', '-k1,1', '-k2,2n', '-k3,3', '-k4,4'],
                           stdin=incoming, stdout=handle, check=True, env=environment)


def verify_score_schema(path, profile, prepared):
    header, rows = table(path)
    expected = profile['contract']['columns'] + ['RawScore', 'PHRED']
    require(header is not None and header.lstrip('#').split('\t') == expected,
            'Publisher/native final score schema differs from the complete ordered incl_anno contract')
    def checked_rows():
        for line in rows:
            require(len(line.split('\t')) == len(expected), 'Final CADD row width differs')
            yield line
    groups = (allele for allele, group in itertools.groupby(checked_rows(), key=score_key))
    for actual, wanted in itertools.zip_longest(groups, (key(row) for row in prepared)):
        require(actual == wanted, 'Final CADD scores lost an allele or contain an unrelated allele')


def calculate(args):
    config, profile_path, profile, block, _, modules = preflight(args)
    task, shared, view, native = modules
    source, output = Path(args.input).resolve(strict=True), Path(args.output).absolute()
    fingerprint = {'input_sha256': digest(source), 'profile_sha256': digest(profile_path),
                   'adapter_sources': adapter_identity(), 'resource_manifest_sha256': profile['contract']['resource_manifest_sha256'],
                   'prescore_inputs_sha256': json_hash(block)}
    with tempfile.TemporaryDirectory(prefix='cadd-source-check-', dir=config['tmp_dir']) as temporary:
        native_rows(source, profile, modules, Path(temporary) / 'selection.gz')
    receipt_path = Path(str(output) + '.cadd-native.json')
    if output.exists() or receipt_path.exists():
        require(output.is_file() and receipt_path.is_file(),
                'Existing CADD output has no matching closed native receipt; retained without overwrite')
        receipt = read_json(receipt_path)
        require(receipt.get('inputs') == fingerprint and receipt.get('output_sha256') == digest(output),
                'Existing CADD output belongs to different input/profile/source/resources; retained')
        return output
    require(output.parent.is_dir(), 'Output parent must exist')
    root = Path(str(output) + '.native')
    require(not root.exists(), 'Native CADD evidence directory exists; retain it and select a fresh output')
    root.mkdir()
    work = Path(tempfile.mkdtemp(prefix='priva-cadd-', dir=config['tmp_dir']))
    environment = {name: value for name, value in os.environ.items() if not name.startswith('CADD_')}
    environment.update(LC_ALL='C', TMPDIR=str(work), PYTHONDONTWRITEBYTECODE='1', PYTHONNOUSERSITE='1')
    with source.open('rb') as incoming:
        compressed_source = incoming.read(2) == b'\x1f\x8b'
    frozen_source = root / ('source.original.vcf.gz' if compressed_source else 'source.original.vcf')
    shutil.copyfile(source, frozen_source)
    require(digest(source) == digest(frozen_source) == fingerprint['input_sha256'], 'Shared source changed while being frozen')
    renamed = root / 'source.grch.vcf.gz'
    contig_map = ROOT / 'data/liftover/ucsc_to_GRC.contig.map.tsv'
    require(contig_map.is_file() and contig_map.stat().st_size > 0, 'Bundled chromosome map missing or empty')
    subprocess.run(['bash', str(ROOT / 'scripts/common_bash_utils.sh'), 'liftover_from_ucsc_to_GRCh',
                    str(frozen_source), str(contig_map), str(renamed)], check=True, env=environment)
    mapping = dict(line.split() for line in contig_map.read_text().splitlines() if line.strip())
    prepared = root / 'prepared.vcf'
    subprocess.run(['bash', '-c', 'set -euo pipefail; gzip -cd "$1" | "$2" "$3" | sort -k1,1 -k2,2n -k4,4 -k5,5 | uniq > "$4"', '_',
                    str(renamed), profile['annotation']['env']['CADD_PYTHON2'], block['prepare']['path'], str(prepared)],
                   check=True, env=environment)
    mapping_path = root / 'source-occurrences.jsonl'
    with OccurrenceMap(frozen_source, renamed, prepared, mapping, work, mapping_path) as occurrences:
        novel_count = score_prepared(config, profile_path, profile, block, modules, fingerprint, root, work,
                                     environment, frozen_source, renamed, prepared, mapping_path, occurrences)
    preflight(args)
    require(digest(source) == fingerprint['input_sha256'] and digest(profile_path) == fingerprint['profile_sha256'] and
            adapter_identity() == fingerprint['adapter_sources'], 'CADD input/profile/adapter changed during execution')
    joined = root / 'scores.tsv'
    with output.open('xb') as handle, joined.open('rb') as result:
        shutil.copyfileobj(result, handle)
    save(receipt_path, {'schema_version': 1, 'stage': 'priva-cadd-native-complete', 'inputs': fingerprint,
        'output_sha256': digest(output), 'evidence_root': str(root), 'mapping_sha256': digest(mapping_path),
        'publisher_files_processed': len(block['files']), 'novel_alleles': novel_count,
        'production_qualification_created': False})
    return output


def score_prepared(config, profile_path, profile, block, modules, fingerprint, root, work,
                   environment, frozen_source, renamed, prepared, mapping_path, occurrences):
    _, _, view, native = modules
    pre = root / 'publisher.pre.tsv'
    pre.write_text('## Prescored variant file\n')
    current = prepared
    split_steps = []
    for ordinal, entry in enumerate(block['files']):
        found, remaining = root / ('publisher.%d.tsv' % ordinal), root / ('remaining.%d.vcf' % ordinal)
        with current.open() as incoming, remaining.open('x') as outgoing:
            subprocess.run([profile['annotation']['env']['CADD_PYTHON2'], block['extract']['path'], '--header',
                            '-p', entry['path'], '--found_out=' + str(found)], stdin=incoming, stdout=outgoing,
                           check=True, env=environment)
        prove_split(vcf_rows(current), vcf_rows(remaining), score_rows(found))
        split_steps.append({'publisher_path': entry['path'],
            'input': {'path': str(current), 'sha256': digest(current)},
            'found': {'path': str(found), 'sha256': digest(found)},
            'remaining': {'path': str(remaining), 'sha256': digest(remaining)}})
        with pre.open('ab') as handle, found.open('rb') as acquired:
            shutil.copyfileobj(acquired, handle)
        current = remaining
        # Deliberately continue after all hits: later indexes/errors and tabix
        # headers are part of the unchanged original prescore route.
    sources = [pre]
    keys_path = root / 'novel-keys.tsv'
    novel_count = 0
    with keys_path.open('x') as keys_output:
        for row in vcf_rows(current):
            keys_output.write('\t'.join(key(row)) + '\n')
            novel_count += 1
    split_path = root / 'publisher-split.json'
    save(split_path, {'schema_version': 1, 'stage': 'upstream-prescore-novel-subset-v1', 'status': 'completed', 'policy_id': POLICY,
        'profile_sha256': fingerprint['profile_sha256'], 'prescore_inputs_sha256': fingerprint['prescore_inputs_sha256'],
        'full_source': {'path': str(frozen_source), 'sha256': digest(frozen_source)},
        'mapping': {'path': str(mapping_path), 'sha256': digest(mapping_path)},
        'prepared': {'path': str(prepared), 'sha256': digest(prepared)}, 'steps': split_steps,
        'novel': {'path': str(current), 'sha256': digest(current)},
        'novel_keys': {'path': str(keys_path), 'sha256': digest(keys_path), 'count': novel_count}})
    if novel_count:
        headers, _ = vcf(renamed)
        novel_source = root / 'novel.shared.vcf.gz'
        write_shared(novel_source, headers, occurrences.shared_records(current))
        bundle = root / 'shared_vep_bundle'
        view.build(novel_source, profile, bundle)
        receipt = read_json(bundle / 'manifest.json')
        receipt['pipeline_provenance'] = {'stage': 'upstream-prescore-novel-subset-v1', **fingerprint,
            'full_source': {'path': str(frozen_source), 'sha256': digest(frozen_source)},
            'mapping': {'path': str(mapping_path), 'sha256': digest(mapping_path)},
            'prescore_split': {'path': str(split_path), 'sha256': digest(split_path)},
            'novel_keys_sha256': digest(keys_path)}
        (bundle / 'manifest.json').write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
        native.run(config, profile_path, bundle, root / 'native_stages')
        sources.append(root / 'native_stages/calculation_bundle/scores.tsv')
    joined = root / 'scores.tsv'
    join_scores(sources, joined, environment)
    verify_score_schema(joined, profile, vcf_rows(prepared))
    return novel_count


def verify_source(args):
    config, _, profile, _, _, modules = preflight(args)
    with tempfile.TemporaryDirectory(prefix='cadd-source-check-', dir=config['tmp_dir']) as temporary:
        native_rows(args.input, profile, modules, Path(temporary) / 'selection.gz')


def vep_options(args):
    _, profile_path, profile, _, _, (_, shared, _, _) = preflight(args)
    settings, context_sha = shared.verify_profile(profile)
    context = settings['selection_context']
    values = [str(ROOT / 'bin/cadd'), '+'.join(context['cadd_transcript_sources']), settings['vep_release'], context['perl_hash_seed']]
    require(all(value and '\n' not in value for value in values), 'Invalid qualified VEP settings')
    print(json.dumps({'event': 'qualified_shared_cadd_vep', 'profile': str(profile_path),
        'profile_sha256': digest(profile_path), 'selection_context_sha256': context_sha,
        'annotation_options': context['annotation_options']}, sort_keys=True), file=sys.stderr, flush=True)
    print('\n'.join(values))


def hub(args, record=False):
    _, profile_path, profile, _, _, _ = preflight(args)
    path = Path(args.hub)
    marker = Path(str(path) + '.cadd-profile.json')
    if not path.exists():
        require(not record, 'Cannot record a missing hub')
        return
    require(path.is_file() and path.stat().st_size > 0, 'Configured hub is not a nonempty file')
    identity = {'profile_sha256': digest(profile_path), 'resource_manifest_sha256': profile['contract']['resource_manifest_sha256'],
                'hub_sha256': digest(path)}
    if record:
        temporary = marker.with_name(marker.name + '.partial')
        save(temporary, {'schema_version': 1, **identity})
        os.replace(temporary, marker)
    else:
        require(marker.is_file() and read_json(marker) == {'schema_version': 1, **identity},
                'Existing hub lacks matching qualified CADD provenance; retained without adoption: ' + str(path))


def extract_misses(args):
    preflight(args)
    header, covered = table(args.covered)
    require(header and all(name in header.lstrip('#').split('\t') for name in ('Chrom', 'Pos', 'Ref', 'Alt')),
            'Covered CADD hub schema is incomplete')
    columns = header.lstrip('#').split('\t')
    indices = [columns.index(name) for name in ('Chrom', 'Pos', 'Ref', 'Alt')]
    headers, rows = vcf(args.input)
    require(not Path(args.output).exists(), 'Miss-subset output already exists')
    config = yaml.safe_load(Path(args.config).read_text())
    with tempfile.TemporaryDirectory(prefix='cadd-exact-misses-', dir=config['tmp_dir']) as temporary:
        database = sqlite3.connect(str(Path(temporary) / 'covered.sqlite'))
        database.execute('PRAGMA cache_size=-2048')
        database.execute('PRAGMA temp_store=FILE')
        database.execute("PRAGMA temp_store_directory='" + temporary.replace("'", "''") + "'")
        database.execute('CREATE TABLE covered (allele TEXT PRIMARY KEY,seen INTEGER DEFAULT 0)')
        plain = Path(temporary) / 'misses.vcf'
        try:
            with database, plain.open('x') as handle:
                for line in covered:
                    fields = line.split('\t')
                    allele = '\t'.join(fields[index] for index in indices)
                    database.execute('INSERT OR IGNORE INTO covered(allele) VALUES(?)', (allele,))
                handle.write('\n'.join(headers) + '\n')
                for row in rows:
                    allele = '\t'.join(key(row))
                    if database.execute('SELECT 1 FROM covered WHERE allele=?', (allele,)).fetchone():
                        database.execute('UPDATE covered SET seen=1 WHERE allele=?', (allele,))
                    else:
                        handle.write('\t'.join(row) + '\n')
                require(database.execute('SELECT 1 FROM covered WHERE seen=0 LIMIT 1').fetchone() is None,
                        'Covered hub returned unrelated alleles')
        finally:
            database.close()
        with Path(args.output).open('xb') as output:
            subprocess.run(['bgzip', '-c', str(plain)], stdout=output, check=True)
        subprocess.run(['tabix', '-p', 'vcf', args.output], check=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('preflight', 'vep-options', 'verify-source', 'calculate', 'hub-check', 'hub-record', 'extract-misses'):
        command = commands.add_parser(name)
        command.add_argument('--config', required=True)
        command.add_argument('--profile')
        if name in ('verify-source', 'calculate', 'extract-misses'):
            command.add_argument('--input', required=True)
        if name in ('calculate', 'extract-misses'):
            command.add_argument('--output', required=True)
        if name in ('hub-check', 'hub-record'):
            command.add_argument('--hub', required=True)
        if name == 'extract-misses':
            command.add_argument('--covered', required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == 'preflight':
            print(preflight(args)[1])
        elif args.command == 'vep-options':
            vep_options(args)
        elif args.command == 'verify-source':
            verify_source(args)
        elif args.command == 'calculate':
            calculate(args)
        elif args.command == 'extract-misses':
            extract_misses(args)
        else:
            hub(args, record=args.command == 'hub-record')
    except (OSError, ValueError, RuntimeError, KeyError, subprocess.CalledProcessError) as error:
        print('PriVA CADD shared workflow: ' + str(error), file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
