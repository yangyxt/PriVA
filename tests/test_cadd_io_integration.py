"""Portable CADD I/O integration tests with synthetic upstream tools only.

Run from a source checkout with ``python -B -m unittest discover -s tests
-p test_cadd_io_integration.py -v``. PRIVA_TEST_TMPDIR selects disposable test
storage. These tests execute the adapter and transport utility, never CADD
predictors, model libraries, VEP, downloads or an installed scientific profile.
The gzip fixture checks adapter wiring, not bgzip or container equivalence.
"""
import argparse
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import sys
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / 'bin/cadd/io_transport.py'


def load_module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H = load_module('cadd_io_hub_fixture', 'tests/test_cadd_task_model_hub.py')
F = load_module('cadd_io_native_fixture', 'tests/test_cadd_native.py')
W = load_module('cadd_io_shared_workflow', 'scripts/cadd_shared_workflow.py')
N, T, S = F.N, F.T, F.S


TRACE = r'''
def trace(kind, **details):
    record = dict(kind=kind, argv=sys.argv[:], cwd=os.getcwd(), path0=sys.path[0],
        LC_ALL=os.environ.get('LC_ALL'), PYTHONUTF8=os.environ.get('PYTHONUTF8'),
        PYTHONIOENCODING=os.environ.get('PYTHONIOENCODING'),
        KERAS_BACKEND=os.environ.get('KERAS_BACKEND'), **details)
    with open(os.environ['FIXTURE_IO_CALLS'], 'a', encoding='utf-8') as target:
        target.write(json.dumps(record, sort_keys=True) + '\n')
'''

REGSEQ = r'''
import gzip, json, os, sys
from pathlib import Path
from fixture_sibling import SENTINEL
''' + TRACE + r'''
trace('regseq_predictor', sibling=SENTINEL)
options = ['--variants', '--model', '--weights', '--reference', '--genome', '--output']
assert len(sys.argv) == 13 and sys.argv[1::2] == options
with gzip.open(sys.argv[2], 'rt', encoding='utf-8') as incoming:
    payload = incoming.read()
# Deliberately use the original implicit plain-writer encoding. The command
# must start with its scoped UTF-8 mode before this synthetic writer runs.
with open(sys.argv[-1], 'x') as target:
    target.write(payload)
if os.environ.get('FIXTURE_FAIL_IO') == 'regseq':
    raise SystemExit(73)
'''

MMSPLICE = r'''
import gzip, json, os, sys
from fixture_sibling import SENTINEL
''' + TRACE + r'''
trace('mmsplice_predictor', sibling=SENTINEL)
assert len(sys.argv) == 7 and sys.argv[1::2] == ['-i', '-g', '-f']
with gzip.open(sys.argv[2], 'rt') as incoming:
    payload = incoming.read()
sys.stdout.write('Variant(CHROM=FIXTURE_DEBUG_TO_REMOVE\n')
sys.stdout.write(payload)
if os.environ.get('FIXTURE_FAIL_IO') == 'mmsplice':
    raise SystemExit(74)
'''


class CaddIOIntegration(unittest.TestCase):
    def hub_fixture(self, assembly):
        fixture = H.TaskModelHubTests(methodName='runTest')
        fixture.setUp()
        self.addCleanup(fixture.tearDown)
        self.assertTrue(HELPER.is_file(), 'The reviewed transport integration must be present')
        fixture.env.update(CADD_ASSEMBLY=assembly, FIXTURE_IO_CALLS=str(fixture.task / 'io-calls.jsonl'))
        for name in ('PYTHONUTF8', 'PYTHONIOENCODING', 'KERAS_BACKEND'):
            fixture.env.pop(name, None)
        # Retain the existing ESM fixture and all its ordered-model/cache
        # assertions. A new strict route actually executes the real helper;
        # it cannot mistake the helper for the old RegSeq command.
        delegate = fixture.tools / 'esm-fixture-python'
        delegate.write_bytes(fixture.interpreter.read_bytes())
        fixture.env['FIXTURE_ESM_DELEGATE'] = str(delegate)
        fixture.interpreter.write_text('#!' + sys.executable + '\n' + r'''
import hashlib, json, os, sys
from pathlib import Path
''' + TRACE + r'''
name = Path(sys.argv[1]).name
trace('route', script=name)
if name == 'io_transport.py':
    os.execv(sys.executable, [sys.executable, '-B'] + sys.argv[1:])
elif name == 'annotateVEPvcf.py':
    payload = sys.stdin.buffer.read()
    trace('annotation_consumer', bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
    sys.stdout.buffer.write(b'#Chrom\tPos\tFixtureOnly\n1\t1\tNA\n')
elif name.startswith('esmScore_'):
    os.execv(sys.executable, [sys.executable, '-B', os.environ['FIXTURE_ESM_DELEGATE']] + sys.argv[1:])
else:
    raise SystemExit('Unexpected scientific tool or bypassed transport: ' + name)
''', encoding='utf-8')
        (fixture.tools / 'bgzip').write_text('#!/bin/sh\n# Adapter fixture only; not a BGZF implementation.\nexec gzip "$@"\n')
        scripts = fixture.source / 'src/scripts'
        for relative, body in (('lib/tools/regulatorySequence/predictVariants.py', REGSEQ),
                               ('lib/tools/MMSplice.py', MMSPLICE)):
            path = scripts / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(body, encoding='utf-8')
            (path.parent / 'fixture_sibling.py').write_text('SENTINEL = "original-script-directory"\n')
        gtf = fixture.data / 'transcripts.gtf'
        gtf.write_text('FIXTURE_ONLY\n')
        tabix = fixture.tools / 'tabix'
        tabix.write_text('#!' + sys.executable + '\n' + r'''
import json, os, sys
from pathlib import Path
''' + TRACE + r'''
trace('tabix')
assert sys.argv[1:3] == ['-p', 'vcf']
Path(sys.argv[3] + '.tbi').write_bytes(b'FIXTURE_INDEX_ONLY\n')
''', encoding='utf-8')
        tabix.chmod(0o755)
        fixture.env.update(CADD_MMSPLICE_PYTHON=str(fixture.interpreter),
            CADD_MMSPLICE_GTF=str(gtf), CADD_REFERENCE_FASTA=str(fixture.data / 'reference.fa'),
            CADD_TABIX=str(tabix))
        # Unicode and duplicate rows are deliberately synthetic. Scientific
        # parsing/model behavior and the actual Python 3.6 codec are separate.
        fixture.payload = ('##fileformat=VCFv4.2\n##fixture=\u2212\u2013\u2265\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
            '1\t1\t.\tA\tC\t.\t.\t.\n').encode('utf-8')
        fixture.payload += b'1\t1\t.\tA\tC\t.\t.\t.\n'
        fixture.input.write_bytes(gzip.compress(fixture.payload))
        return fixture

    def native_fixture(self):
        fixture = F.NativeCadd(methodName='runTest')
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        return fixture

    def snapshot(self, fixture):
        paths = [fixture.input]
        for directory in (fixture.source, fixture.data, fixture.hub):
            paths += [path for path in directory.rglob('*') if path.is_file()]
        result = {}
        for path in paths:
            value = path.lstat()
            self.assertTrue(stat.S_ISREG(value.st_mode))
            result[str(path)] = tuple(getattr(value, name) for name in
                ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns', 'st_mode', 'st_nlink')) + (
                hashlib.sha256(path.read_bytes()).hexdigest(),)
        return result

    def records(self, fixture):
        return [json.loads(line) for line in Path(fixture.env['FIXTURE_IO_CALLS']).read_text().splitlines()]

    def stage(self, fixture):
        paths = list(fixture.task.glob('annotation.*'))
        self.assertEqual(len(paths), 1)
        return paths[0]

    def check_regseq(self, fixture, records):
        stage = self.stage(fixture)
        routes = [item for item in records if item['kind'] == 'route' and item['script'] == 'io_transport.py']
        original = str(fixture.source / 'src/scripts/lib/tools/regulatorySequence/predictVariants.py')
        arguments = [original, '--variants', str(stage / 'esm.vcf.gz'),
            '--model', fixture.env['CADD_REGSEQ_MODEL'], '--weights', fixture.env['CADD_REGSEQ_WEIGHTS'],
            '--reference', fixture.env['CADD_REGSEQ_REFERENCE'], '--genome', fixture.env['CADD_REGSEQ_GENOME'],
            '--output', str(stage / 'regseq.vcf.gz')]
        self.assertEqual(routes[0]['argv'][1:], [str(HELPER), 'regseq', '--bgzip', fixture.env['CADD_BGZIP'],
            '--work-dir', str(stage), '--'] + arguments)
        record, = [item for item in records if item['kind'] == 'regseq_predictor']
        self.assertEqual(record['argv'][:-1], arguments[:-1])
        plain = Path(record['argv'][-1])
        self.assertEqual(plain.parent.parent, stage)
        self.assertTrue(plain.parent.name.startswith('regseq-utf8-'))
        self.assertEqual(plain.name, 'regseq.vcf')
        self.assertEqual(plain.read_bytes(), fixture.payload)
        self.assertEqual(record['path0'], str(Path(original).parent))
        self.assertEqual(record['sibling'], 'original-script-directory')
        self.assertEqual((record['LC_ALL'], record['PYTHONUTF8'], record['PYTHONIOENCODING']), ('C', '1', None))
        self.assertEqual(gzip.decompress((stage / 'regseq.vcf.gz').read_bytes()), fixture.payload)
        return stage

    def check_annotation(self, fixture, records):
        record, = [item for item in records if item['kind'] == 'annotation_consumer']
        self.assertEqual((record['bytes'], record['sha256']),
            (len(fixture.payload), hashlib.sha256(fixture.payload).hexdigest()))
        self.assertEqual((record['LC_ALL'], record['PYTHONUTF8'], record['PYTHONIOENCODING'], record['KERAS_BACKEND']),
                         ('C', None, None, None))
        self.assertEqual(gzip.decompress(fixture.output.read_bytes()), b'#Chrom\tPos\tFixtureOnly\n1\t1\tNA\n')

    def test_shared_grch37_executes_regseq_helper_with_original_operands(self):
        fixture = self.hub_fixture('GRCh37')
        before = self.snapshot(fixture)
        done = fixture.run_adapter()
        self.assertEqual(done.returncode, 0, done.stderr.decode('utf-8', errors='replace'))
        records = self.records(fixture)
        self.check_regseq(fixture, records)
        self.check_annotation(fixture, records)
        fixture.check_calls(H.MODELS)
        self.assertFalse(any(item['kind'] in ('tabix', 'mmsplice_predictor') for item in records))
        self.assertEqual(self.snapshot(fixture), before)

    def test_shared_grch38_executes_mmsplice_helper_and_original_pipeline(self):
        fixture = self.hub_fixture('GRCh38')
        before = self.snapshot(fixture)
        done = fixture.run_adapter()
        self.assertEqual(done.returncode, 0, done.stderr.decode('utf-8', errors='replace'))
        records = self.records(fixture)
        stage = self.check_regseq(fixture, records)
        self.check_annotation(fixture, records)
        fixture.check_calls(H.MODELS)
        original = str(fixture.source / 'src/scripts/lib/tools/MMSplice.py')
        arguments = [original, '-i', str(stage / 'regseq.vcf.gz'), '-g', fixture.env['CADD_MMSPLICE_GTF'],
                     '-f', fixture.env['CADD_REFERENCE_FASTA']]
        route, = [item for item in records if item['kind'] == 'route'
                  and item['script'] == 'io_transport.py' and item['argv'][2] == 'mmsplice']
        self.assertEqual(route['argv'][1:], [str(HELPER), 'mmsplice', '--input', str(stage / 'regseq.vcf.gz'), '--'] + arguments)
        record, = [item for item in records if item['kind'] == 'mmsplice_predictor']
        self.assertEqual(record['argv'], arguments)
        self.assertEqual(record['cwd'], str(fixture.task))
        self.assertEqual(record['path0'], str(Path(original).parent))
        self.assertEqual(record['sibling'], 'original-script-directory')
        self.assertEqual((record['LC_ALL'], record['PYTHONUTF8'], record['PYTHONIOENCODING'], record['KERAS_BACKEND']),
                         ('C', None, 'utf-8', 'tensorflow'))
        self.assertEqual(gzip.decompress((stage / 'mmsplice.vcf.gz').read_bytes()), fixture.payload)
        self.assertTrue((stage / 'regseq.vcf.gz.tbi').is_file())
        kinds = [item['kind'] for item in records]
        self.assertLess(kinds.index('regseq_predictor'), kinds.index('tabix'))
        self.assertLess(kinds.index('tabix'), kinds.index('mmsplice_predictor'))
        self.assertLess(kinds.index('mmsplice_predictor'), kinds.index('annotation_consumer'))
        self.assertEqual(self.snapshot(fixture), before)

    def test_regseq_failure_keeps_partial_and_stops_later_stages(self):
        fixture = self.hub_fixture('GRCh38')
        fixture.env['FIXTURE_FAIL_IO'] = 'regseq'
        before = self.snapshot(fixture)
        done = fixture.run_adapter()
        self.assertEqual(done.returncode, 73)
        records = self.records(fixture)
        stage = self.stage(fixture)
        record, = [item for item in records if item['kind'] == 'regseq_predictor']
        self.assertEqual(Path(record['argv'][-1]).read_bytes(), fixture.payload)
        self.assertFalse((stage / 'regseq.vcf.gz').exists())
        self.assertFalse(fixture.output.exists())
        self.assertFalse(any(item['kind'] in ('tabix', 'mmsplice_predictor', 'annotation_consumer') for item in records))
        self.assertEqual(self.snapshot(fixture), before)

    def test_mmsplice_failure_keeps_pipeline_partial_and_stops_annotation(self):
        fixture = self.hub_fixture('GRCh38')
        fixture.env['FIXTURE_FAIL_IO'] = 'mmsplice'
        before = self.snapshot(fixture)
        done = fixture.run_adapter()
        self.assertEqual(done.returncode, 74)
        records = self.records(fixture)
        stage = self.check_regseq(fixture, records)
        self.assertTrue((stage / 'mmsplice.vcf.gz').is_file())
        self.assertTrue(any(item['kind'] == 'mmsplice_predictor' for item in records))
        self.assertFalse(any(item['kind'] == 'annotation_consumer' for item in records))
        self.assertFalse(fixture.output.exists())
        self.assertEqual(self.snapshot(fixture), before)

    def test_helper_is_bound_by_both_task_scopes_and_shared_workflow(self):
        key = 'bin/cadd/io_transport.py'
        self.assertEqual(T.source_identity('native')[key], T.digest(HELPER))
        # Check the Nextflow identity membership in either repository, even
        # where the original PriVA tree deliberately has no Nextflow files.
        with mock.patch.object(T, 'digest', return_value='a' * 64):
            native, nextflow = T.source_identity('native'), T.source_identity('nextflow')
        self.assertIn(key, native)
        self.assertIn(key, nextflow)
        self.assertNotIn('cadd.nf', native)
        self.assertIn('cadd.nf', nextflow)
        self.assertEqual(W.adapter_identity()[key], W.digest(HELPER))

    def test_helper_drift_rejects_native_preflight_and_task_completion(self):
        fixture = self.native_fixture()
        path, _ = fixture.install()
        before, digest = T.source_identity('native'), T.digest
        def changed(value):
            return '0' * 64 if Path(value) == HELPER else digest(value)
        # Qualification is mocked solely to reach the engineering source
        # guard. No fake profile becomes scientifically qualified or runs.
        with mock.patch.object(T, 'digest', side_effect=changed), mock.patch.object(S, 'verify_profile'), \
                mock.patch.object(N, 'verify_prescoring_profile'), mock.patch.object(T, 'stage_preflight') as stages:
            with self.assertRaisesRegex(T.ContractError, 'adapter source changed'):
                N.preflight(fixture.config, path)
            with self.assertRaisesRegex(T.ContractError, 'adapter source changed'):
                T.receipt({}, 'unused', 'annotation', 'unused', 'unused', 0, False, before, source_scope='native')
        stages.assert_not_called()
        self.assertEqual(list(fixture.work.iterdir()), [])

    def test_helper_drift_rejects_shared_workflow_completed_output_reuse(self):
        fixture = self.native_fixture()
        source = fixture.write(fixture.root / 'synthetic.input', b'FIXTURE_ONLY\n')
        profile_path = fixture.write(fixture.root / 'synthetic.profile', b'FIXTURE_PROFILE_ONLY\n')
        output = fixture.write(fixture.root / 'synthetic.output', b'FIXTURE_OUTPUT_ONLY\n')
        profile = {'contract': {'resource_manifest_sha256': 'a' * 64}}
        fingerprint = dict(input_sha256=W.digest(source), profile_sha256=W.digest(profile_path),
            adapter_sources=W.adapter_identity(), resource_manifest_sha256='a' * 64,
            prescore_inputs_sha256=W.json_hash({}))
        receipt = fixture.write(Path(str(output) + '.cadd-native.json'),
            N.encoded(dict(inputs=fingerprint, output_sha256=W.digest(output))))
        saved = {str(p): (p.stat().st_ino, p.stat().st_mtime_ns, p.read_bytes()) for p in (output, receipt)}
        args = argparse.Namespace(input=str(source), output=str(output))
        setup = (dict(tmp_dir=str(fixture.work)), profile_path, profile, {}, None, (None, None, None, None))
        digest = W.digest
        def changed(value):
            return '0' * 64 if Path(value) == HELPER else digest(value)
        with mock.patch.object(W, 'preflight', return_value=setup), mock.patch.object(W, 'native_rows'), \
                mock.patch('subprocess.run', side_effect=AssertionError('No models or other subprocesses')):
            self.assertEqual(W.calculate(args), output)
            with mock.patch.object(W, 'digest', side_effect=changed), \
                    self.assertRaisesRegex(W.WorkflowError, 'different input/profile/source/resources'):
                W.calculate(args)
        self.assertEqual(saved, {str(p): (p.stat().st_ino, p.stat().st_mtime_ns, p.read_bytes()) for p in (output, receipt)})
        self.assertFalse(Path(str(output) + '.native').exists())

    def test_source_refresh_creates_distinct_unqualified_profile_without_reinstall(self):
        fixture = self.native_fixture()
        original_path, original = fixture.install()
        old_bytes = original_path.read_bytes()
        old_stat = original_path.stat()
        component_bytes = fixture.receipt.read_bytes()
        manifest_path = Path(original['resource_manifest_path'])
        manifest_bytes = manifest_path.read_bytes()
        digest = T.digest
        def changed(value):
            return '0' * 64 if Path(value) == HELPER else digest(value)
        with mock.patch.object(T, 'digest', side_effect=changed), \
                mock.patch('subprocess.run', side_effect=AssertionError('No solve, install, probe or models')), \
                mock.patch('urllib.request.urlopen', side_effect=AssertionError('No download')):
            with self.assertRaisesRegex(T.ContractError, 'adapter source changed'):
                N.check(fixture.config, original_path)
            current_path, current = fixture.install()
            self.assertEqual(N.check(fixture.config, current_path), current)
            self.assertEqual(N.install_profile(fixture.config, fixture.receipt, verify_only=True)['path'], str(current_path))
        self.assertNotEqual(current_path, original_path)
        self.assertNotEqual(T.digest(current_path), T.digest(original_path))
        self.assertEqual(current['profile_id'], original['profile_id'])
        self.assertEqual(current['software_id'], original['software_id'])
        self.assertEqual(current['resource_manifest_path'], original['resource_manifest_path'])
        self.assertEqual(current['calculation'], original['calculation'])
        self.assertFalse(current['production_activation'])
        self.assertFalse(current['scientific_qualification'])
        self.assertIsNone(current['annotation']['shared_vep']['qualification'])
        self.assertIsNone(current['prescoring']['qualification'])
        self.assertEqual(original_path.read_bytes(), old_bytes)
        self.assertEqual((original_path.stat().st_ino, original_path.stat().st_mtime_ns), (old_stat.st_ino, old_stat.st_mtime_ns))
        self.assertEqual(fixture.receipt.read_bytes(), component_bytes)
        self.assertEqual(manifest_path.read_bytes(), manifest_bytes)


if __name__ == '__main__':
    unittest.main()
