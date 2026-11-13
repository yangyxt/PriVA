"""Alpha deployment contracts; no annotation, motif calculation or licensed data."""
import hashlib
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
WRAPPER = ROOT / 'scripts/build_cluster_alpha.sbatch'
COMMAND_SHA = '569619479a4267693bc2f9caca56cc12761d8abb87ee3a232138b0624ae57b99'


def original_motif_command():
    # Exact invocation from install_utils.sh at scientific source
    # 09223ef5a8b0ceb4c6d61055c6eefe62b51ba903. Today's installer has a separate
    # publication route; it cannot supply this historical wrapper's oracle.
    block = (ROOT / 'tests/fixtures/alpha_motif_command_09223ef.sh').read_text()
    if hashlib.sha256(block.encode()).hexdigest() != COMMAND_SHA:
        raise AssertionError('Frozen motif command block changed; review source provenance')
    return shlex.split(block.replace('\\\n', ''))


class AlphaDeployment(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='alpha_wrapper_', dir=ROOT / 'tests')
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.reference = self.folder / 'reference'
        (self.reference / 'scripts').mkdir(parents=True)
        self.motif_py = self.reference / 'scripts/am_pick_missense_constrained_segments.py'
        self.motif_py.write_text('# fixture timestamp only; never executed\n')
        self.mapping = self.folder / 'mapping.tsv'
        self.mapping.write_text('Gene stable ID\tTranscript stable ID\tHGNC symbol\tUniProt ID\tRefSeq mRNA ID\n'
                                'ENSG1\tENST1\tGENE1\tP1\tNM_1\n')
        self.mapping_sha = hashlib.sha256(self.mapping.read_bytes()).hexdigest()
        self.output = self.folder / 'motifs.pkl'
        self.vcf = self.folder / 'annotated.vcf.gz'
        self.vcf.write_text('placeholder; not used as genomic data\n')
        os.utime(self.vcf, (100, 100))
        os.utime(self.motif_py, (50, 50))
        self.config = self.folder / 'separate.yaml'
        self.config.write_text(f'assembly: hg19\nalphamissense_vep_vcf: before_annotation.vcf.gz\n'
                               f'alphamissense_intolerant_motifs: {self.output}\nthreads: 4\n')
        self.trace = self.folder / 'calls.txt'
        self.env = dict(os.environ, TEST_CONFIG=str(self.config), TEST_REFERENCE=str(self.reference),
                        TEST_MAPPING=str(self.mapping), TEST_SHA=self.mapping_sha, TEST_TRACE=str(self.trace),
                        TEST_VCF=str(self.vcf), SLURM_JOB_ID='fixture', SLURM_CPUS_PER_TASK='4')

    def run_wrapper(self, body='alpha_run_stages "$TEST_CONFIG" "$TEST_REFERENCE" "$TEST_MAPPING" "$TEST_SHA"'):
        stubs = r'''
            record_call() { printf '%q ' "$@" >> "$TEST_TRACE"; printf '\n' >> "$TEST_TRACE"; }
            bash() {
                record_call bash "$@"
                if [[ ${2:-} == AlphaMissense_anno ]]; then
                    [[ ${ANNOTATION_STATUS:-0} == 0 ]] || return "$ANNOTATION_STATUS"
                    sed -i "s|^alphamissense_vep_vcf:.*|alphamissense_vep_vcf: $TEST_VCF|" "$TEST_CONFIG"
                    if [[ ${MUTATE_DURING_ANNOTATION:-0} == 1 ]]; then printf '\nchanged' >> "$TEST_MAPPING"; fi
                fi
                return 0
            }
            python() {
                record_call python "$@"
                if [[ ${MUTATE_DURING_MOTIF:-0} == 1 ]]; then printf '\nchanged' >> "$TEST_MAPPING"; fi
                if [[ ${MUTATE_DURING_MOTIF:-0} == 2 ]]; then touch "$TEST_MAPPING"; fi
                return "${MOTIF_STATUS:-0}"
            }
        '''
        return subprocess.run(['bash', '-c', 'source "$1"\n' + stubs + '\n' + body,
                               'fixture', str(WRAPPER)], env=self.env, text=True, capture_output=True)

    def calls(self):
        return [shlex.split(line) for line in self.trace.read_text().splitlines()] if self.trace.exists() else []

    def test_missing_mapping_stops_before_source_or_environment_access(self):
        self.env['TEST_MAPPING'] = str(self.folder / 'absent.tsv')
        result = self.run_wrapper('main "$TEST_CONFIG" "$TEST_REFERENCE" "$TEST_MAPPING" "$TEST_SHA"')
        self.assertEqual(result.returncode, 67, result.stderr)
        self.assertIn('explicit absolute', result.stderr)
        self.assertNotIn('REFERENCE_COMMIT', result.stderr)
        self.assertEqual(self.calls(), [])

    def test_missing_hash_stops_before_annotation(self):
        self.env['TEST_SHA'] = ''
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 67)
        self.assertIn('approved mapping SHA256', result.stderr)
        self.assertEqual(self.calls(), [])

    def test_wrong_hash_stops_before_annotation(self):
        self.env['TEST_SHA'] = '0' * 64
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 67)
        self.assertIn('checksum mismatch', result.stderr)
        self.assertEqual(self.calls(), [])

    def test_empty_mapping_stops_before_annotation(self):
        self.mapping.write_bytes(b'')
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 67)
        self.assertEqual(self.calls(), [])

    def test_relative_mapping_is_rejected(self):
        self.env['TEST_MAPPING'] = 'mapping.tsv'
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 67)
        self.assertEqual(self.calls(), [])

    def test_original_command_order_defaults_and_updated_vcf_are_preserved(self):
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.calls()
        self.assertEqual(calls[0], ['bash', str(self.reference / 'scripts/install_utils.sh'),
                                   'AlphaMissense_anno', str(self.config)])
        substitutions = {'${pick_intolerant_motifs_py}': str(self.motif_py),
                         '${alphamissense_vep_vcf}': str(self.vcf), '${SCRIPT_DIR}': str(self.reference / 'scripts'),
                         '${assembly}': 'hg19', '${alphamissense_intolerant_motifs}': str(self.output),
                         '${threads}': '4'}
        original = original_motif_command()
        for index, argument in enumerate(original):
            for key, value in substitutions.items():
                argument = argument.replace(key, value)
            original[index] = argument
        self.assertEqual(calls[1], original + ['--nm_to_enst_tsv', str(self.mapping)])
        self.assertEqual(calls[2], ['bash', str(self.reference / 'scripts/common_bash_utils.sh'),
                                   'update_yaml', str(self.config), 'alphamissense_intolerant_motifs', str(self.output)])
        self.assertIn('Fresh motif generation completed', result.stdout)

    def test_original_timestamp_reuse_preserves_output_and_reports_unverified_provenance(self):
        self.output.write_bytes(b'existing motif output')
        os.utime(self.output, (200, 200))
        before = self.output.stat().st_mtime_ns
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(self.calls()), 1)  # Annotation may itself reuse its own output.
        self.assertEqual(self.output.read_bytes(), b'existing motif output')
        self.assertEqual(self.output.stat().st_mtime_ns, before)
        self.assertIn('skip this step', result.stderr)
        self.assertIn('mapping provenance is unverified', result.stdout)
        self.assertNotIn('Fresh motif generation', result.stdout)

    def test_equal_timestamp_does_not_reuse(self):
        self.output.write_bytes(b'existing')
        os.utime(self.output, (100, 100))
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(self.calls()), 3)

    def test_annotation_failure_does_not_run_motif(self):
        self.env['ANNOTATION_STATUS'] = '41'
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 41)
        self.assertEqual(len(self.calls()), 1)

    def test_motif_failure_does_not_update_config(self):
        self.env['MOTIF_STATUS'] = '42'
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 42)
        self.assertEqual(len(self.calls()), 2)
        self.assertNotIn('Fresh motif generation completed', result.stdout)

    def test_mapping_change_during_annotation_blocks_motif(self):
        self.env['MUTATE_DURING_ANNOTATION'] = '1'
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 67)
        self.assertEqual(len(self.calls()), 1)

    def test_mapping_change_during_motif_blocks_completion(self):
        self.env['MUTATE_DURING_MOTIF'] = '1'
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 67)
        self.assertEqual(len(self.calls()), 2)
        self.assertNotIn('Fresh motif generation completed', result.stdout)

    def test_same_bytes_with_changed_mapping_metadata_block_completion(self):
        self.env['MUTATE_DURING_MOTIF'] = '2'
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 67)
        self.assertIn('changed during motif generation', result.stderr)
        self.assertEqual(hashlib.sha256(self.mapping.read_bytes()).hexdigest(), self.mapping_sha)
        self.assertEqual(len(self.calls()), 2)
        self.assertNotIn('Fresh motif generation completed', result.stdout)

    def test_changed_frozen_source_is_rejected(self):
        (self.reference / 'REFERENCE_COMMIT').write_text('09223ef5a8b0ceb4c6d61055c6eefe62b51ba903\n')
        (self.reference / 'scripts/install_utils.sh').write_text('# different source\n')
        result = self.run_wrapper('alpha_validate_source "$TEST_REFERENCE"')
        self.assertEqual(result.returncode, 65)
        self.assertIn('source checksum mismatch', result.stderr)
        self.assertEqual(self.calls(), [])


if __name__ == '__main__':
    unittest.main()
