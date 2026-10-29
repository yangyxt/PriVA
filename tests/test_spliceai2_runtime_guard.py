"""Runtime check of the SpliceAI score files: one version per run, each with its own cutoffs."""
import gzip
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/annotation_vcf.sh'
REFERENCE = ROOT / 'tests/fixtures/installer-054b494/scripts/annotation_vcf.sh'
HEADER = b'##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
V2_MARKERS = (b'##SpliceAI2Version=2', b'##source=PriVA_SpliceAI2_projection',
              b'##SpliceAI2Policy=annotation_only_no_ACMG_or_CADD_activation')


class SpliceAI2RuntimeGuard(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='spliceai2_guard_',
            dir=os.environ.get('PRIVA_TEST_TMPDIR', os.environ.get('TMPDIR')))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.scripts = self.folder / 'source/scripts'
        self.scripts.mkdir(parents=True)
        for name in ('annotation_vcf.sh', 'common_bash_utils.sh', 'anno_vep_args', 'anno_main_args',
                     'spliceai_cache_input.py', 'spliceai_cutoffs.py'):
            shutil.copyfile(ROOT / 'scripts' / name, self.scripts / name)
        self.script = self.scripts / 'annotation_vcf.sh'
        # Parser fixture retains the argument-to-local-variable behavior needed
        # at both real entry points; scientific tools below are trace-only.
        (self.scripts / 'argparse.bash').write_text('''argparse() {
            while [[ $# -gt 0 ]]; do
                case "$1" in
                    --config) config=$2 ;;
                    --input_vcf) input_vcf=$2 ;;
                    --spliceai_snv_prescore) spliceai_snv_prescore=$2 ;;
                    --spliceai_indel_prescore) spliceai_indel_prescore=$2 ;;
                    *) return 87 ;;
                esac
                shift 2
            done
        }
        ''')
        self.config = self.folder / 'config.yaml'
        self.snv = self.folder / 'snv.vcf.gz'
        self.indel = self.folder / 'indel.vcf.gz'
        self.snv.write_bytes(gzip.compress(HEADER))
        self.indel.write_bytes(gzip.compress(HEADER))
        self.trace = self.folder / 'tool-trace.txt'
        self.input = self.folder / 'input.vcf.gz'
        self.input.write_bytes(gzip.compress(HEADER))
        self.settings = dict(assembly='hg19', spliceai_snv_prescore=str(self.snv),
                             spliceai_indel_prescore=str(self.indel),
                             hub_vcf_file=str(self.folder / 'never-create-hub/hub.vcf.gz'),
                             hub_cadd_file=str(self.folder / 'never-create-cadd/hub.tsv.gz'),
                             ref_genome='/reference.fa', vep_cache_dir='/cache', vep_plugins_dir='/plugins',
                             vep_plugins_cachedir='/plugin-cache', loftee_repo='/loftee', threads=4, loeuf_prescore='loeuf.tsv.gz',
                             pext_enabled=False)
        self.env = dict(os.environ, TEST_CONFIG=str(self.config), TEST_TRACE=str(self.trace),
                        TEST_INPUT=str(self.input), TEST_SNV=str(self.snv), TEST_INDEL=str(self.indel))
        self.save_config()

    def save_config(self):
        self.config.write_text(yaml.safe_dump(self.settings, sort_keys=False))

    def run_shell(self, body):
        return subprocess.run(['bash', '-c', 'load_script() { source "$1"; }; load_script "$1"; ' + body, 'test', str(self.script)],
                              env=self.env, text=True, capture_output=True)

    def guard(self, paths=None):
        self.save_config()
        body = 'priva_check_spliceai_versions "$TEST_CONFIG"'
        if paths is not None:
            self.env['TEST_OVERRIDE_SNV'], self.env['TEST_OVERRIDE_INDEL'] = map(str, paths)
            body += ' "$TEST_OVERRIDE_SNV" "$TEST_OVERRIDE_INDEL"'
        return self.run_shell(body)

    def test_original_annotation_oracle_is_unchanged(self):
        self.assertEqual(hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),
                         '6b7ad6fd8eb118841cf36819193f7ce2ee14363432c66f0da75d7788c2d8d280')

    def v2_header(self, path, marker=V2_MARKERS[0], compressed=True):
        text = HEADER.replace(b'#CHROM', marker + b'\n#CHROM')
        path.write_bytes(gzip.compress(text) if compressed else text)

    def test_native_tsv_must_be_converted_before_use(self):
        columns = ['chrom', 'pos', 'ref', 'alt', 'strand']
        for effect in ('donor_gain', 'donor_loss', 'acceptor_gain', 'acceptor_loss'):
            columns += [f'{effect}_{field}_{rank}' for field in ('delta_score', 'ref_score', 'alt_score', 'dist') for rank in range(3)]
        for effect in ('jxn_gain', 'jxn_loss'):
            columns += [f'{effect}_{field}_{rank}' for field in ('delta_score', 'ref_score', 'alt_score', 'donor_dist', 'acceptor_dist') for rank in range(3)]
        table = gzip.compress(('\t'.join(columns + ['spliceai2_summary_score']) + '\n').encode())
        self.snv.write_bytes(table)
        self.indel.write_bytes(table)
        self.settings['assembly'] = 'hg38'
        result = self.guard()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('must be converted to indexed VCF', result.stderr)

    def test_v1_uses_only_original_two_keys(self):
        result = self.guard()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('SpliceAI v1 cutoffs: splice_altering 0.8', result.stderr)

    def test_v2_pair_is_accepted_with_spliceai2_cutoffs(self):
        self.settings['assembly'] = 'hg38'
        for marker in V2_MARKERS:
            for compressed in (True, False):
                with self.subTest(marker=marker, compressed=compressed):
                    self.v2_header(self.snv, marker, compressed)
                    self.v2_header(self.indel, marker, compressed)
                    result = self.guard()
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn('SpliceAI2 cutoffs: splice_altering', result.stderr)

    def test_each_v2_header_marker_beside_v1_is_refused_as_mixed_versions(self):
        self.settings['assembly'] = 'hg38'
        for path in (self.snv, self.indel):
            for marker in V2_MARKERS:
                for compressed in (True, False):
                    with self.subTest(role=path.name, marker=marker, compressed=compressed):
                        self.v2_header(path, marker, compressed)
                        result = self.guard()
                        self.assertNotEqual(result.returncode, 0)
                        self.assertIn('different versions', result.stderr)
                        path.write_bytes(gzip.compress(HEADER))

    def test_v2_pair_is_refused_for_a_grch37_config(self):
        self.v2_header(self.snv)
        self.v2_header(self.indel)
        for assembly in ('hg19', 'GRCh37'):
            with self.subTest(assembly=assembly):
                self.settings['assembly'] = assembly
                result = self.guard()
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('configured assembly is GRCh37', result.stderr)

    def test_effective_cli_paths_take_precedence_over_config_paths(self):
        v2 = self.folder / 'override.vcf.gz'
        self.v2_header(v2)
        result = self.guard((v2, self.indel))
        self.assertNotEqual(result.returncode, 0)
        self.settings['spliceai_snv_prescore'] = str(v2)
        result = self.guard((self.snv, self.indel))
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_valid_v1_headers_and_unrelated_mentions_remain_allowed(self):
        self.snv.write_bytes(HEADER.replace(b'#CHROM', b'##note=SpliceAI2 is a different resource\n#CHROM'))
        result = self.guard()
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_missing_resources_keep_existing_legacy_validation(self):
        self.snv.unlink()
        self.indel.unlink()
        result = self.guard()
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_header_scan_does_not_decompress_or_qualify_the_body(self):
        # The deliberately corrupt far-end CRC must not be read. Existing cache
        # installation does full integrity validation; this guard only detects
        # a version marker before the first VCF record.
        compressed = bytearray(gzip.compress(HEADER + b'chr1\t10\t.\tA\tC\t.\t.\t.\n' * 100000))
        compressed[-8] ^= 1
        self.snv.write_bytes(compressed)
        result = self.guard()
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_header_length_and_line_count_have_hard_bounds(self):
        for text in (b'##' + b'x' * 65536 + b'\n' + HEADER,
                     b'##x\n' * 4096 + HEADER,
                     (b'##' + b'x' * 65532 + b'\n') * 17 + HEADER):
            with self.subTest(length=len(text)):
                self.snv.write_bytes(gzip.compress(text))
                result = self.guard()
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('Cannot use the SpliceAI score files', result.stderr)

    def test_unreadable_truncated_and_incomplete_headers_fail_closed(self):
        for data in (b'not a VCF\n', b'##fileformat=VCFv4.2\n', b'\x1f\x8b\x08'):
            with self.subTest(data=data):
                self.snv.write_bytes(data)
                result = self.guard()
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('SpliceAI', result.stderr)

    def test_main_mixed_versions_stop_before_tools_or_hub_creation(self):
        self.settings['assembly'] = 'hg38'
        self.v2_header(self.snv)
        self.save_config()
        result = self.run_shell('''
            bcftools() { printf 'bcftools\\n' >> "$TEST_TRACE"; }
            check_path() { printf 'check_path\\n' >> "$TEST_TRACE"; }
            main_workflow --config "$TEST_CONFIG" --input_vcf "$TEST_INPUT"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('different versions', result.stderr)
        self.assertFalse(self.trace.exists())
        self.assertFalse((self.folder / 'never-create-hub').exists())
        self.assertFalse((self.folder / 'never-create-cadd').exists())

    def test_direct_vep_mixed_versions_stop_before_any_path_checks_or_tools(self):
        self.settings['assembly'] = 'hg38'
        self.save_config()
        self.v2_header(self.snv)
        result = self.run_shell('''
            check_path() { printf 'check_path\\n' >> "$TEST_TRACE"; }
            vep() { printf 'vep\\n' >> "$TEST_TRACE"; }
            anno_VEP_data --config "$TEST_CONFIG" --input_vcf "$TEST_INPUT"
        ''')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('different versions', result.stderr)
        self.assertFalse(self.trace.exists())

    def test_direct_vep_v2_pair_runs_with_spliceai2_files(self):
        self.settings['assembly'] = 'hg38'
        self.v2_header(self.snv)
        self.v2_header(self.indel)
        self.save_config()
        result = self.run_shell('''
            check_path() { return 0; }
            check_vcf_infotags() { return 1; }
            randomID() { printf 'FIXED\\n'; }
            capture() { printf '%s\\n' "$@" >> "$TEST_TRACE"; }
            vep() { capture vep "$@"; }
            bcftools() { capture bcftools "$@"; }
            tabix() { capture tabix "$@"; }
            announce_remove_tmps() { capture announce_remove_tmps "$@"; }
            display_vcf() { capture display_vcf "$@"; }
            anno_VEP_data --config "$TEST_CONFIG" --input_vcf "$TEST_INPUT"
        ''')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('SpliceAI2 cutoffs', result.stderr)
        self.assertIn(f'SpliceAI,snv={self.snv},indel={self.indel},cutoff=0.5'.encode(), self.trace.read_bytes())

    def test_v1_arguments_match_frozen_source_with_reviewed_plugin_directory_fix(self):
        body = '''
            check_path() { return 0; }
            check_vcf_infotags() { return 1; }
            randomID() { printf 'FIXED\\n'; }
            capture() { printf '%s\\n' "$@" >> "$TEST_TRACE"; }
            vep() { capture vep "$@"; }
            bcftools() { capture bcftools "$@"; }
            tabix() { capture tabix "$@"; }
            announce_remove_tmps() { capture announce_remove_tmps "$@"; }
            display_vcf() { capture display_vcf "$@"; }
            anno_VEP_data --config "$TEST_CONFIG" --input_vcf "$TEST_INPUT"
        '''
        current = self.run_shell(body)
        self.assertEqual(current.returncode, 0, current.stderr)
        current_trace = self.trace.read_bytes()
        self.assertIn(b'SpliceAI,snv=', current_trace)
        self.assertIn(b'cutoff=0.5', current_trace)
        self.trace.unlink()
        self.script.write_bytes(REFERENCE.read_bytes())
        original = self.run_shell(body)
        self.assertEqual(original.returncode, 0, original.stderr)
        original_trace = self.trace.read_bytes()
        # The reviewed VEP116 loader repair changes this single argument from
        # one nonexistent comma-joined path to VEP's whitespace-separated list.
        # Every other argument and following command remains byte-exact; this
        # trace fixture is not a scientific equality claim across VEP profiles.
        old_directory = b'--dir_plugins\n/plugins,/loftee\n'
        corrected_directory = b'--dir_plugins\n/plugins /loftee\n'
        self.assertEqual(original_trace.count(old_directory), 1)
        self.assertEqual(current_trace.count(corrected_directory), 1)
        self.assertEqual(original_trace.replace(old_directory, corrected_directory, 1), current_trace)


if __name__ == '__main__':
    unittest.main()
