"""Lossless partition and worker-failure fixtures; these do not run VEP/models."""
import argparse
from contextlib import redirect_stdout
import gzip
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

SCRIPTS = Path(__file__).parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location("vep_chromosome_jobs", SCRIPTS / "vep_chromosome_jobs.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
HEADER = (b'##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n'
          b'##contig=<ID=chr2,length=100>\n##contig=<ID=NT_1.1,length=100>\n'
          b'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
ROWS = [b'chr1\t10\tx\tGG\tAA,C\t.\tPASS\t.\tGT\t./.\n',
        b'chr1\t10\tx\tGG\tAA,C\t.\tPASS\t.\tGT\t./.\n',
        b'chr2\t20\tnonvariant\tA\t.\t.\t.\t.\tGT\t0/0\n',
        b'NT_1.1\t2\tbnd\tC\tC]chr1:10]\t.\tPASS\t.\tGT\t1\n']
ANNOTATED_HEADER = HEADER.replace(b'#CHROM',b'##INFO=<ID=CSQ,Number=.,Type=String,Description="Allele|Feature|HGNC_ID">\n#CHROM')


def annotate(row):
    fields=row.rstrip(b'\n').split(b'\t')
    fields[7]=b'CSQ=A|ENST1|' if fields[7]==b'.' else fields[7]+b';CSQ=A|ENST1|'
    return b'\t'.join(fields)+b'\n'


class JobTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source_root = self.root / 'repo'
        (self.source_root / 'scripts').mkdir(parents=True)
        for name in ('annotation_vcf.sh', 'vep_chromosome_cache.py', 'vep_chromosome_jobs.py',
                     'clinvar_gene_refresh.py', 'argparse.bash', 'anno_vep_args'):
            (self.source_root / 'scripts' / name).write_text('fixture code\n')
        self.config = self.root / 'config.yaml'
        self.config.write_text('cache_root: ' + str(self.root) + '\nassembly: GRCh38\n')
        self.source = self.root / 'source.vcf.gz'
        self.set_rows(ROWS)
        self.work = self.root / 'work'
        self.args = argparse.Namespace(config=self.config, source_root=self.source_root,
                    input_vcf=self.source, work_dir=self.work, bcftools='unused',
                    forks=8, cpus_per_chromosome=10, max_parallel=12)

    def set_rows(self, rows):
        with gzip.open(self.source, 'wb') as out:
            out.write(HEADER + b''.join(rows))
        self.indexed = dict(input=M.identity(self.source), indexes=[],
                            header_sha256=hashlib.sha256(HEADER).hexdigest(),
                            contig_order=['chr1','chr2','NT_1.1'],
                            records_by_contig={'chr1':2,'chr2':1,'NT_1.1':1},total_records=4)

    def split(self):
        with mock.patch.object(M, 'index_inventory', return_value=self.indexed), \
             mock.patch.object(M, 'profile', return_value={'files':{}, 'fixture':True}), \
             redirect_stdout(io.StringIO()):
            M.split(self.args)
        self.manifest = self.work / 'chromosomes.json'
        return json.loads(self.manifest.read_text())

    def worker_args(self, index=0):
        return argparse.Namespace(manifest=self.manifest, manifest_sha256=M.digest(self.manifest),contig_index=index)

    def run_worker(self, result=0, *, profiles=None):
        args = self.worker_args()
        def execute(argv, **kwargs):
            output = Path(argv[argv.index('--output-vcf')+1])
            stats = Path(argv[argv.index('--stats-file')+1])
            output.write_bytes(b'fixture annotation or retained partial\n')
            stats.write_bytes(b'fixture stats\n')
            return subprocess.CompletedProcess(argv,result)
        profile_options = dict(side_effect=profiles) if profiles else dict(return_value={'files':{},'fixture':True})
        with mock.patch.dict(os.environ, {'SLURM_CPUS_PER_TASK':'10'}), \
             mock.patch.object(M, 'profile', **profile_options), \
             mock.patch.object(M.subprocess, 'run', side_effect=execute), redirect_stdout(io.StringIO()):
            M.worker(args)

    def test_split_preserves_every_byte_duplicate_nonvariant_and_noncanonical(self):
        plan = self.split()
        self.assertEqual(plan['source_records'],4)
        bodies = [Path(row['input_path']).read_bytes() for row in plan['shards']]
        self.assertTrue(all(body.startswith(HEADER) for body in bodies))
        self.assertEqual(b''.join(body[len(HEADER):] for body in bodies),b''.join(ROWS))
        self.assertEqual([(r['source_ordinal_start'],r['source_ordinal_end']) for r in plan['shards']],[(1,2),(3,3),(4,4)])
        self.assertEqual(plan['forks'],8)
        self.assertEqual(plan['cpus_per_chromosome'],10)
        self.assertEqual(plan['max_parallel'],12)
        self.assertFalse(plan['upstream_record_filtering_applied_during_split'])
        self.assertFalse(plan['canonical_publication'])
        self.assertEqual(Path(plan['worker_config']['identity']['path']).read_bytes(),self.config.read_bytes())

    def test_nine_forks_do_not_fit_ten_cpus_for_pinned_spawn_loop(self):
        self.args.forks = 9
        with self.assertRaisesRegex(ValueError,'forks \+ two'):
            self.split()
        self.assertFalse(self.work.exists())

    def test_wrong_index_count_retains_shards_without_ready_manifest(self):
        self.indexed['records_by_contig']['chr1'] = 3
        with self.assertRaisesRegex(ValueError,'record counts differ'):
            self.split()
        self.assertFalse((self.work/'chromosomes.json').exists())
        self.assertTrue((self.work/'contig-000/input.vcf').exists())

    def test_corrupt_closed_copy_cannot_be_blessed_as_the_source(self):
        original_pinned=M.pinned
        def corrupt_copy(path):
            path=Path(path)
            if path.name == 'input.vcf':
                with path.open('ab') as out:out.write(b'CORRUPT')
            return original_pinned(path)
        with mock.patch.object(M,'pinned',side_effect=corrupt_copy):
            with self.assertRaisesRegex(ValueError,'copy bytes differ'):
                self.split()
        self.assertFalse((self.work/'chromosomes.json').exists())
        self.assertTrue((self.work/'contig-000/input.vcf').exists())

    def test_reordered_source_fails_without_a_plan(self):
        self.set_rows([ROWS[0],ROWS[2],ROWS[1],ROWS[3]])
        with self.assertRaisesRegex(ValueError,'order regressed'):
            self.split()
        self.assertFalse((self.work/'chromosomes.json').exists())

    def test_existing_directory_is_preserved(self):
        self.work.mkdir()
        sentinel = self.work/'existing';sentinel.write_text('keep')
        with self.assertRaisesRegex(ValueError,'Fresh work directory'):
            self.split()
        self.assertEqual(sentinel.read_text(),'keep')

    def test_success_is_explicitly_pending_validation(self):
        self.split();self.run_worker()
        value=json.loads((self.work/'contig-000/worker-result.json').read_text())
        self.assertEqual(value['state'],'annotation_finished_pending_validation')
        self.assertEqual(value['returncode'],0)
        self.assertTrue(value['scientific_resources_unchanged'])
        self.assertFalse(value['output_equivalence_validated'])
        self.assertFalse(value['canonical_publication'])
        self.assertEqual(value['argv'][-1],'8')

    def test_failed_worker_retains_partial_and_nonzero_result(self):
        self.split()
        with self.assertRaises(SystemExit) as error:
            self.run_worker(42)
        self.assertEqual(error.exception.code,42)
        value=json.loads((self.work/'contig-000/worker-result.json').read_text())
        self.assertEqual(value['state'],'failed')
        self.assertEqual(value['returncode'],42)
        self.assertTrue((self.work/'contig-000/output.vcf').exists())

    def test_changed_resources_fail_even_when_vep_exit_is_zero(self):
        self.split()
        with self.assertRaisesRegex(ValueError,'resources changed while'):
            self.run_worker(profiles=[{'files':{},'fixture':True},{'files':{},'fixture':False}])
        value=json.loads((self.work/'contig-000/worker-result.json').read_text())
        self.assertEqual(value['returncode'],0)
        self.assertEqual(value['state'],'failed')
        self.assertFalse(value['scientific_resources_unchanged'])

    def test_postcondition_exception_records_closed_child_failure(self):
        self.split()
        with self.assertRaisesRegex(ValueError,'resources changed while'):
            self.run_worker(profiles=[{'files':{},'fixture':True},ValueError('fixture policy parse failure')])
        value=json.loads((self.work/'contig-000/worker-result.json').read_text())
        self.assertEqual(value['returncode'],0)
        self.assertEqual(value['state'],'failed')
        self.assertEqual(value['postcondition_error'],'ValueError: fixture policy parse failure')

    def test_changed_argparse_source_is_bound_to_the_worker_plan(self):
        self.split()
        (self.source_root/'scripts/argparse.bash').write_text('changed')
        with self.assertRaisesRegex(ValueError,'Pinned file identity changed'):
            self.run_worker()
        self.assertFalse((self.work/'contig-000/started.json').exists())

    def test_changed_source_code_fails_before_worker_starts(self):
        self.split()
        (self.source_root/'scripts/annotation_vcf.sh').write_text('changed')
        with self.assertRaisesRegex(ValueError,'Pinned file identity changed'):
            self.run_worker()
        self.assertFalse((self.work/'contig-000/started.json').exists())

    def test_insufficient_slurm_allocation_fails_before_start(self):
        self.split()
        with mock.patch.dict(os.environ, {'SLURM_CPUS_PER_TASK':'9'}), \
             mock.patch.object(M,'profile', return_value={'files':{},'fixture':True}):
            with self.assertRaisesRegex(ValueError,'declared Slurm CPU allocation'):
                M.worker(self.worker_args())
        self.assertFalse((self.work/'contig-000/started.json').exists())

    def test_retry_never_overwrites_previous_attempt(self):
        self.split();self.run_worker()
        before=(self.work/'contig-000/output.vcf').read_bytes()
        with self.assertRaisesRegex(ValueError,'attempt already exists'):
            self.run_worker()
        self.assertEqual((self.work/'contig-000/output.vcf').read_bytes(),before)

    def test_manifest_checksum_mismatch_is_rejected(self):
        self.split();args=self.worker_args();args.manifest_sha256='0'*64
        with self.assertRaisesRegex(ValueError,'Manifest checksum'):
            M.read_plan(args)

    def test_status_distinguishes_finished_from_ready(self):
        self.split();self.run_worker()
        with redirect_stdout(io.StringIO()) as out:
            M.status(self.worker_args())
        value=json.loads(out.getvalue())
        self.assertEqual([r['state'] for r in value['chromosomes']],
                         ['annotation_finished_pending_validation','pending','pending'])
        self.assertFalse(value['installation_accepted'])

    def fixture_annotation(self, index=0, rows=None):
        plan=self.split()
        shard=plan['shards'][index]
        if rows is None:
            rows=ROWS[:2] if index==0 else [ROWS[index+1]]
        folder=Path(shard['folder'])
        output=folder/'output.vcf'
        output.write_bytes(ANNOTATED_HEADER+b''.join(annotate(row) for row in rows))
        result=dict(kind='priva-vep-chromosome-worker-result-v1',
                    manifest_sha256=M.digest(self.manifest),contig=shard['contig'],
                    state='annotation_finished_pending_validation',returncode=0,
                    scientific_resources_unchanged=True,output=M.identity(output))
        (folder/'worker-result.json').write_text(json.dumps(result))
        args=self.worker_args(index);args.vep_source_root=None
        return folder,args

    def validate_fixture(self,args):
        with mock.patch.object(M,'profile',return_value={'files':{},'fixture':True}), \
             redirect_stdout(io.StringIO()):
            M.validate(args)

    def test_validator_preserves_duplicate_occurrences_and_missing_feature_values(self):
        folder,args=self.fixture_annotation()
        self.validate_fixture(args)
        value=json.loads((folder/'source-occurrences.validated.json').read_text())
        self.assertEqual((value['source_records'],value['emitted_records'],value['upstream_omitted_records']),(2,2,0))
        self.assertEqual(value['output']['sha256'],M.digest(folder/'output.vcf'))
        self.assertTrue(value['source_occurrences_validated'])
        self.assertFalse(value['annotation_equivalence_validated'])
        self.assertFalse(value['canonical_publication'])
        with redirect_stdout(io.StringIO()) as out:M.status(args)
        self.assertEqual(json.loads(out.getvalue())['chromosomes'][0]['state'],
                         'source_occurrences_validated_annotations_pending_equivalence')

    def test_validator_reports_changed_allele_and_retains_failure_ledger(self):
        folder,args=self.fixture_annotation(rows=[ROWS[0].replace(b'AA,C',b'AT,C'),ROWS[1]])
        with self.assertRaisesRegex(ValueError,'changed output occurrence'):
            self.validate_fixture(args)
        self.assertTrue((folder/'validation-omitted-source-occurrences.jsonl').exists())
        self.assertFalse((folder/'source-occurrences.validated.json').exists())

    def test_validator_rejects_extra_emitted_occurrence(self):
        folder,args=self.fixture_annotation(rows=ROWS[:2]+[ROWS[0]])
        with self.assertRaisesRegex(ValueError,'Extra output occurrence'):
            self.validate_fixture(args)
        self.assertFalse((folder/'source-occurrences.validated.json').exists())

    def test_nonvariant_omissions_require_explicit_source_policy_proof(self):
        folder,args=self.fixture_annotation(index=1,rows=[])
        with self.assertRaisesRegex(ValueError,'Incomplete VCF occurrence'):
            self.validate_fixture(args)
        self.assertFalse((folder/'source-occurrences.validated.json').exists())

    def test_authorized_upstream_nonvariant_omission_has_exact_ordinal_ledger(self):
        folder,args=self.fixture_annotation(index=1,rows=[])
        policy=dict(kind='fixture-policy',nonvariant_omissions_authorized=True,allow_non_variant=False)
        with mock.patch.object(M,'vep116_omission_policy',return_value=policy):
            self.validate_fixture(args)
        value=json.loads((folder/'source-occurrences.validated.json').read_text())
        self.assertEqual((value['source_records'],value['emitted_records'],value['upstream_omitted_records']),(1,0,1))
        ledger=json.loads((folder/'validation-omitted-source-occurrences.jsonl').read_text())
        self.assertEqual(ledger['source_ordinal'],3)
        self.assertEqual(ledger['source_row_sha256'],hashlib.sha256(ROWS[2]).hexdigest())

    def test_validator_refuses_output_changed_since_child_closed(self):
        folder,args=self.fixture_annotation()
        with (folder/'output.vcf').open('ab') as out:out.write(b'changed')
        with self.assertRaisesRegex(ValueError,'output changed since its child closed'):
            self.validate_fixture(args)
        self.assertFalse((folder/'source-occurrences.validated.json').exists())


if __name__ == '__main__':
    unittest.main()
