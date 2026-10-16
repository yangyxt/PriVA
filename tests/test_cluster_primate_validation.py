"""Exact PrimateAI validation against the literal original formatting pipeline."""
import contextlib
import gzip
import io
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/import_cluster_licensed.sbatch'
BODY = SCRIPT.read_text().split("<<'PY'\n", 1)[1].split('\nPY\n', 1)[0]
ORIGINAL = r'''gunzip -cf "$1" | sed '12s/.*/#&/' | sed '/^$/d' | awk 'NR<12{print $0;next}{print $0 | "sort -k1,1 -k 2,2n -V"}' | bgzip'''


class ExactPrimateValidation(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='primate_exact_fixture_', dir=ROOT / 'tests')
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.raw = self.folder / 'original.tsv.gz'
        self.formatted = self.folder / 'formatted.tsv.bgz'
        self.header = b''.join(b'# metadata ' + str(i).encode() + b'\n' for i in range(1, 11))
        self.header += b'\nchr\tpos\tref\talt\tscore\n'
        self.rows = [b'chr10\t10\tA\tG\t0.5', b'chr2\t2\tG\tT\t0.4',
                     b'chr2\t2\tA\tC\t0.1', b'chr2\t2\tA\tC\t0.1']
        self.namespace = {'__name__': 'primate_validation_fixture'}
        exec(compile(BODY, str(SCRIPT), 'exec'), self.namespace)
        self.validator = self.namespace['validate_primateai_format']
        self.addCleanup(mock.patch.stopall)
        mock.patch.dict(os.environ, {'TMPDIR': str(self.folder)}).start()
        self.prepare(self.header + b'\n'.join(self.rows) + b'\n')

    def prepare(self, content):
        self.raw.write_bytes(gzip.compress(content))
        result = subprocess.run(['bash', '-o', 'pipefail', '-c', ORIGINAL, 'original', str(self.raw)],
                                capture_output=True, check=True)
        self.formatted.write_bytes(result.stdout)

    def validate(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return self.validator(self.raw, self.formatted)

    def assert_clean(self):
        self.assertEqual(list(self.folder.glob('primateai-exact-*')), [])

    def test_normal_duplicate_tie_and_header_bytes_match_literal_original(self):
        original_bytes = self.raw.read_bytes()
        derived_bytes = self.formatted.read_bytes()
        receipt = self.validate()
        self.assertEqual(receipt['status'], 'exact_uncompressed_match')
        self.assertTrue(receipt['raw_parity_verified'])
        self.assertEqual(receipt['matched_uncompressed_bytes'], len(gzip.decompress(derived_bytes)))
        self.assertEqual(receipt['header_records'], 11)
        self.assertEqual(set(receipt['exit_statuses'].values()), {0})
        self.assertEqual(receipt['commands'][3]['argv'], ['sort', '-k1,1', '-k', '2,2n', '-V'])
        self.assertIn('GNU coreutils', receipt['tool_versions']['sort'])
        self.assertIn('LC_ALL', receipt['locale'])
        self.assertEqual(self.raw.read_bytes(), original_bytes)
        self.assertEqual(self.formatted.read_bytes(), derived_bytes)
        self.assert_clean()

    def test_unterminated_final_data_record_matches_awk_newline(self):
        self.prepare(self.header + b'\n'.join(self.rows))
        self.assertTrue(gzip.decompress(self.formatted.read_bytes()).endswith(b'\n'))
        self.validate()
        self.assert_clean()

    def test_unterminated_header_only_input_matches_awk_newline(self):
        self.prepare(b'# first\n# last')
        self.assertEqual(gzip.decompress(self.formatted.read_bytes()), b'# first\n# last\n')
        self.assertEqual(self.validate()['header_records'], 2)
        self.assert_clean()

    def test_cr_and_whitespace_only_lines_are_not_stripped(self):
        self.prepare((self.header + b'\n'.join(self.rows) + b'\n').replace(b'\n', b'\r\n'))
        self.validate()
        self.prepare(self.header.replace(b'\nchr\t', b' \nchr\t') + b'\n'.join(self.rows) + b'\n')
        self.validate()
        self.assert_clean()

    def test_dropped_duplicate_row_is_rejected(self):
        lines = gzip.decompress(self.formatted.read_bytes()).splitlines(keepends=True)
        del lines[11]  # Remove one occurrence of the duplicated first data row.
        self.formatted.write_bytes(gzip.compress(b''.join(lines)))
        with self.assertRaisesRegex(RuntimeError, 'byte_compare exited 1'):
            self.validate()
        self.assert_clean()

    def test_replaced_row_with_same_count_is_rejected(self):
        data = gzip.decompress(self.formatted.read_bytes()).replace(b'0.5', b'0.6')
        self.formatted.write_bytes(gzip.compress(data))
        with self.assertRaisesRegex(RuntimeError, 'byte_compare exited 1'):
            self.validate()
        self.assert_clean()

    def test_changed_tie_order_is_rejected(self):
        lines = gzip.decompress(self.formatted.read_bytes()).splitlines(keepends=True)
        lines[11], lines[13] = lines[13], lines[11]
        self.formatted.write_bytes(gzip.compress(b''.join(lines)))
        with self.assertRaisesRegex(RuntimeError, 'byte_compare exited 1'):
            self.validate()
        self.assert_clean()

    def test_first_byte_and_trailing_byte_mismatches_are_rejected(self):
        original = gzip.decompress(self.formatted.read_bytes())
        for data in (b'!' + original[1:], original + b'X'):
            with self.subTest(data_length=len(data)):
                self.formatted.write_bytes(gzip.compress(data))
                with self.assertRaisesRegex(RuntimeError, 'byte_compare exited 1'):
                    self.validate()
                self.assert_clean()

    @staticmethod
    def corrupt_crc(path):
        data = bytearray(path.read_bytes())
        data[-8] ^= 1
        path.write_bytes(data)

    def test_source_decompression_error_is_rejected(self):
        self.corrupt_crc(self.raw)
        with self.assertRaisesRegex(RuntimeError, 'source_gunzip exited'):
            self.validate()
        self.assert_clean()

    def test_derived_error_after_complete_matching_output_is_rejected(self):
        # Use one gzip member so its CRC can fail after emitting all expected bytes.
        self.formatted.write_bytes(gzip.compress(gzip.decompress(self.formatted.read_bytes())))
        self.corrupt_crc(self.formatted)
        with self.assertRaisesRegex(RuntimeError, 'derived_gunzip exited'):
            self.validate()
        self.assert_clean()

    def test_partial_sort_output_with_nonzero_exit_is_rejected(self):
        fake_bin = self.folder / 'bin'
        fake_bin.mkdir()
        sorter = fake_bin / 'sort'
        real_sort = shutil.which('sort')
        sorter.write_text('#!/usr/bin/env bash\n'
                          f'if [[ $1 == --version ]]; then exec {real_sort} "$@"; fi\n'
                          'IFS= read -r first || exit 23\ncat >/dev/null\nprintf "%s\\n" "$first"\nexit 23\n')
        sorter.chmod(0o755)
        with mock.patch.dict(os.environ, {'PATH': str(fake_bin) + ':' + os.environ['PATH']}):
            with self.assertRaisesRegex(RuntimeError, 'sort exited 23'):
                self.validate()
        self.assert_clean()

    def test_early_mismatch_reaps_every_owned_child(self):
        self.prepare(self.header + (b'chr1\t1\tA\tC\t0.1\n' * 100000))
        data = gzip.decompress(self.formatted.read_bytes())
        self.formatted.write_bytes(gzip.compress(b'!' + data[1:]))
        processes = []
        actual_popen = subprocess.Popen
        def record(*args, **kwargs):
            process = actual_popen(*args, **kwargs)
            processes.append(process)
            return process
        with mock.patch('subprocess.Popen', side_effect=record):
            with self.assertRaisesRegex(RuntimeError, 'byte_compare exited 1'):
                self.validate()
        self.assertTrue(processes)
        self.assertTrue(all(process.poll() is not None for process in processes))
        self.assert_clean()

    def assert_changed_identity_rejected(self, target):
        actual = self.namespace['file_identity']
        calls = 0
        def change_before_final_check(path):
            nonlocal calls
            if path == target:
                calls += 1
                if calls == 2:
                    stat = path.stat()
                    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
            return actual(path)
        with mock.patch.dict(self.namespace, {'file_identity': change_before_final_check}):
            with self.assertRaisesRegex(RuntimeError, 'changed during exact validation'):
                self.validate()
        self.assert_clean()

    def test_source_identity_change_is_rejected(self):
        self.assert_changed_identity_rejected(self.raw)

    def test_derived_identity_change_is_rejected(self):
        self.assert_changed_identity_rejected(self.formatted)


if __name__ == '__main__':
    unittest.main()
