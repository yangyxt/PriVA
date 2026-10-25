"""Official gnomAD BGZF suffix support with real htslib format/index checks."""
import gzip
import hashlib
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'scripts/common_bash_utils.sh'
REFERENCE = ROOT / 'tests/fixtures/installer-054b494/scripts/common_bash_utils.sh'
VCF = (b'##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
       b'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
       b'chr1\t10\t.\tA\tC\t.\tPASS\t.\nchr1\t20\t.\tG\tT\t.\tPASS\t.\n')


class VcfBgzValidation(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='vcf_bgz_', dir=ROOT / 'tests')
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)

    def bgzf(self, suffix, data=VCF):
        path = self.folder / ('input' + suffix)
        with path.open('wb') as stream:
            subprocess.run(['bgzip', '-c'], input=data, stdout=stream, check=True)
        return path

    def validate(self, path, function='check_vcf_validity', reference=False, expected=1):
        return subprocess.run(['bash', '-c', 'load() { source "$1"; }; load "$1"; "$2" "$3" "$4"',
                               'test', str(REFERENCE if reference else SCRIPT), function, str(path), str(expected)],
                              text=True, capture_output=True)

    def test_historical_reference_rejects_official_suffix_and_keeps_its_hash(self):
        self.assertEqual(hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),
                         '4cd2cffed0d6a46d1280eac97b918c1dd6c717c881dc41b0a819fbecaed439b9')
        path = self.bgzf('.vcf.bgz')
        result = self.validate(path, reference=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Unsupported VCF file format', result.stderr)
        result = self.validate(path)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_both_suffixes_create_queryable_indexes_without_changing_vcf_bytes(self):
        for suffix in ('.vcf.gz', '.vcf.bgz'):
            with self.subTest(suffix=suffix):
                path = self.bgzf(suffix)
                original = path.read_bytes()
                result = self.validate(path)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(path.read_bytes(), original)
                self.assertTrue(Path(str(path) + '.tbi').is_file())
                query = subprocess.check_output(['tabix', str(path), 'chr1:10-10'])
                self.assertEqual(query, b'chr1\t10\t.\tA\tC\t.\tPASS\t.\n')

    def test_record_count_check_is_not_bypassed(self):
        for suffix in ('.vcf.gz', '.vcf.bgz'):
            with self.subTest(suffix=suffix):
                result = self.validate(self.bgzf(suffix), expected=3)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('expected to have 3 lines', result.stderr)

    def test_existing_current_tbi_is_reused_for_both_suffixes(self):
        for suffix in ('.vcf.gz', '.vcf.bgz'):
            with self.subTest(suffix=suffix):
                path = self.bgzf(suffix)
                subprocess.run(['tabix', '-p', 'vcf', str(path)], check=True)
                index = Path(str(path) + '.tbi')
                os.utime(path, (time.time() - 20, time.time() - 20))
                before = (index.read_bytes(), index.stat().st_mtime_ns)
                result = self.validate(path, 'check_vcf_index')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual((index.read_bytes(), index.stat().st_mtime_ns), before)

    def test_stale_csi_is_rebuilt_using_existing_preference(self):
        for suffix in ('.vcf.gz', '.vcf.bgz'):
            with self.subTest(suffix=suffix):
                path = self.bgzf(suffix)
                subprocess.run(['bcftools', 'index', '-f', str(path)], check=True)
                index = Path(str(path) + '.csi')
                os.utime(path, (time.time() - 20, time.time() - 20))
                os.utime(index, (time.time() - 40, time.time() - 40))
                result = self.validate(path, 'check_vcf_index')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertGreater(index.stat().st_mtime_ns, path.stat().st_mtime_ns)
                self.assertFalse(Path(str(path) + '.tbi').exists())

    def test_nonvcf_bgzf_and_truncated_input_still_fail(self):
        for suffix in ('.vcf.gz', '.vcf.bgz'):
            for data in (b'not\ta\tVCF\n', b''):
                with self.subTest(suffix=suffix, data=data):
                    path = self.bgzf(suffix, data)
                    result = self.validate(path)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertFalse(Path(str(path) + '.tbi').exists())
            with self.subTest(suffix=suffix, truncated=True):
                path = self.bgzf(suffix)
                path.write_bytes(path.read_bytes()[:20])
                result = self.validate(path)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(Path(str(path) + '.tbi').exists())

    def test_regular_gzip_cannot_masquerade_as_bgzf(self):
        for suffix in ('.vcf.gz', '.vcf.bgz'):
            with self.subTest(suffix=suffix):
                path = self.folder / ('ordinary' + suffix)
                path.write_bytes(gzip.compress(VCF))
                result = self.validate(path)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(Path(str(path) + '.tbi').exists())

    def test_other_extensions_remain_rejected(self):
        for suffix in ('.tsv.bgz', '.vcf.bgzip', '.txt.gz'):
            with self.subTest(suffix=suffix):
                path = self.bgzf(suffix)
                result = self.validate(path, 'check_vcf_index')
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('Unsupported VCF file format', result.stderr)


if __name__ == '__main__':
    unittest.main()
