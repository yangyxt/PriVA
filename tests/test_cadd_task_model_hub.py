"""Task-local checkpoint transport controls; fake tools, no scientific models."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
WRAPPER = ROOT / 'bin/cadd/upstream_annotate.sh'
MODELS = ['esm1v_t33_650M_UR90S_' + str(n) for n in range(1, 6)]


class TaskModelHubTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='cadd-task-hub-',
            dir=os.environ.get('PRIVA_TEST_TMPDIR', os.environ.get('TMPDIR')))
        self.root = Path(self.temporary.name)
        self.task = self.root / 'task'
        self.task.mkdir()
        self.source = self.root / 'unchanged source'
        scripts = self.source / 'src/scripts'
        for name in ('VCF2vepVCF.py', 'annotateVEPvcf.py',
                     'lib/tools/esmScore/esmScore_missense_av_fast.py',
                     'lib/tools/esmScore/esmScore_frameshift_av.py',
                     'lib/tools/esmScore/esmScore_inFrame_av.py',
                     'lib/tools/regulatorySequence/predictVariants.py'):
            path = scripts / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('# Fixture tool identity only; fake interpreter below.\n')
        self.hub = self.root / 'external model hub'
        (self.hub / 'checkpoints').mkdir(parents=True)
        for model in MODELS:
            (self.hub / 'checkpoints' / (model + '.pt')).write_bytes(('FIXTURE_ONLY:' + model).encode())
        self.tools = self.root / 'fixture tools'
        self.tools.mkdir()
        git = self.tools / 'git'
        git.write_text('#!/bin/sh\ncase "$*" in *rev-parse*) printf "%s\\n" 1fead0c1e02d3f0a8599bdb04e6b3a2dd70b4866;; *diff*) exit 0;; *) exit 99;; esac\n')
        git.chmod(0o755)
        bgzip = self.tools / 'bgzip'
        bgzip.write_text('#!/bin/sh\nexec gzip -c\n')
        bgzip.chmod(0o755)
        self.interpreter = self.tools / 'fixture-python'
        self.interpreter.write_text('#!' + sys.executable + '\n' + r'''
import gzip, json, os, sys
from pathlib import Path
name=Path(sys.argv[1]).name
args=sys.argv[2:]
def value(option):return args[args.index(option)+1]
if name=='annotateVEPvcf.py':
    sys.stdin.buffer.read();sys.stdout.write('#Chrom\tPos\tFixtureOnly\n1\t1\tNA\n')
else:
    if name.startswith('esmScore_'):
        hub=Path(value('--model-directory'));checkpoints=hub/'checkpoints'
        models=[args[i+1] for i,arg in enumerate(args) if arg=='--model']
        try:os.makedirs(checkpoints)
        except FileExistsError:pass
        # A cache bookkeeping write must be task-local, even if scientific
        # payloads themselves are already complete and immutable.
        (hub/(name+'.bookkeeping')).write_text('FIXTURE_ONLY\n')
        record={'tool':name,'models':models,'hub':str(hub),'args':args,
                'targets':[str((checkpoints/(model+'.pt')).resolve(strict=True)) for model in models]}
        with open(os.environ['FIXTURE_CALLS'],'a')as out:out.write(json.dumps(record)+'\n')
        if os.environ.get('FIXTURE_FAIL_TOOL')==name:raise SystemExit(71)
        source,target=value('--input'),value('--output')
    else:source,target=value('--variants'),value('--output')
    Path(target).write_bytes(Path(source).read_bytes())
''')
        self.interpreter.chmod(0o755)
        self.input = self.root / 'input.vcf.gz'
        self.input.write_bytes(gzip.compress(b'##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'))
        self.data = self.root / 'fixture data'
        self.data.mkdir()
        files = {}
        for name in ('peptide.fa', 'model.json', 'weights.h5', 'reference.fa',
                     'reference.fa.fai', 'reference.fa.genome', 'references.cfg'):
            path = self.data / name
            path.write_text('FIXTURE_ONLY\n')
            files[name] = str(path)
        self.env = dict(os.environ, PATH=str(self.tools) + ':' + os.environ['PATH'],
            CADD_SOURCE_ROOT=str(self.source), CADD_PYTHON2=str(self.interpreter), CADD_ASSEMBLY='GRCh37',
            CADD_ESM_PYTHON=str(self.interpreter), CADD_ESM_MODEL_DIR=str(self.hub),
            CADD_ESM_MODELS=','.join(MODELS), CADD_ESM_BATCH_SIZE='1', CADD_ANNOTATION_INPUT='shared-vep',
            CADD_ESM_TRANSCRIPTS=files['peptide.fa'],
            CADD_REGSEQ_PYTHON=str(self.interpreter), CADD_REGSEQ_MODEL=files['model.json'],
            CADD_REGSEQ_WEIGHTS=files['weights.h5'], CADD_REGSEQ_REFERENCE=files['reference.fa'],
            CADD_REGSEQ_GENOME=files['reference.fa.genome'], CADD_REFERENCE_CONFIG=files['references.cfg'],
            CADD_BGZIP=str(bgzip), PYTHONDONTWRITEBYTECODE='1', PYTHONNOUSERSITE='1',
            FIXTURE_CALLS=str(self.task / 'calls.jsonl'))
        self.output = self.task / 'annotations.tsv.gz'

    def tearDown(self):
        self.temporary.cleanup()

    def snapshot(self):
        return {p.name: (p.stat().st_ino, p.stat().st_mtime_ns, hashlib.sha256(p.read_bytes()).hexdigest())
                for p in (self.hub / 'checkpoints').iterdir()}

    def run_adapter(self, *, readonly=False):
        command = ['bash', str(WRAPPER), str(self.input), str(self.output), str(self.task)]
        if readonly:
            command = ['bwrap', '--die-with-parent', '--unshare-net', '--ro-bind', '/', '/',
                '--bind', str(self.task), str(self.task), '--tmpfs', '/tmp', '--proc', '/proc',
                '--dev', '/dev', '--chdir', str(self.task), '--', *command]
        return subprocess.run(command, cwd=self.task, env=self.env, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, timeout=30)

    def check_calls(self, expected_models):
        calls = [json.loads(line) for line in Path(self.env['FIXTURE_CALLS']).read_text().splitlines()]
        self.assertEqual([x['tool'] for x in calls], ['esmScore_missense_av_fast.py',
            'esmScore_frameshift_av.py', 'esmScore_inFrame_av.py'])
        for call in calls:
            self.assertEqual(call['models'], expected_models)
            self.assertTrue(Path(call['hub']).is_relative_to(self.task))
            self.assertNotEqual(Path(call['hub']), self.hub)
            self.assertEqual(call['args'][call['args'].index('--transcripts')+1], self.env['CADD_ESM_TRANSCRIPTS'])
            self.assertEqual(call['args'][call['args'].index('--batch-size')+1], '1')
            for model, target in zip(expected_models, call['targets']):
                local = Path(call['hub']) / 'checkpoints' / (model + '.pt')
                self.assertTrue(local.is_symlink())
                self.assertTrue(os.path.samefile(local, self.hub / 'checkpoints' / (model + '.pt')))
                self.assertEqual(target, str((self.hub / 'checkpoints' / (model + '.pt')).resolve()))

    def test_every_original_model_occurrence_uses_same_payload_in_task_cache(self):
        models = [MODELS[1], MODELS[0], MODELS[1], *MODELS[2:]]
        self.env['CADD_ESM_MODELS'] = ','.join(models)
        before = self.snapshot()
        done = self.run_adapter()
        self.assertEqual(done.returncode, 0, done.stderr.decode())
        self.check_calls(models)
        self.assertEqual(self.snapshot(), before)
        self.assertEqual(gzip.decompress(self.output.read_bytes()), b'#Chrom\tPos\tFixtureOnly\n1\t1\tNA\n')

    @unittest.skipUnless(shutil.which('bwrap'), 'Linux Bubblewrap is needed for a real read-only mount')
    def test_all_three_tools_create_only_task_local_cache_under_readonly_boundary(self):
        before = self.snapshot()
        done = self.run_adapter(readonly=True)
        self.assertEqual(done.returncode, 0, done.stderr.decode())
        self.check_calls(MODELS)
        self.assertEqual(self.snapshot(), before)
        self.assertEqual(set(p.name for p in self.hub.iterdir()), {'checkpoints'})

    def test_missing_required_payload_fails_before_tool_or_staged_cache(self):
        (self.hub / 'checkpoints' / (MODELS[-1] + '.pt')).unlink()
        done = self.run_adapter()
        self.assertNotEqual(done.returncode, 0)
        self.assertIn(b'required file is not readable', done.stderr)
        self.assertFalse(Path(self.env['FIXTURE_CALLS']).exists())
        self.assertFalse(self.output.exists())
        self.assertEqual(list(self.task.glob('annotation.*')), [])

    def test_upstream_tool_failure_exit_and_unpublished_output_are_preserved(self):
        self.env['FIXTURE_FAIL_TOOL'] = 'esmScore_inFrame_av.py'
        before = self.snapshot()
        done = self.run_adapter()
        self.assertEqual(done.returncode, 71)
        self.assertFalse(self.output.exists())
        self.assertEqual(self.snapshot(), before)
        self.check_calls(MODELS)


if __name__ == '__main__':
    unittest.main()
