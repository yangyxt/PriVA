"""Bounded derived-header compatibility; no raw catalogue or VEP execution."""
import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "scripts/alphamissense_cache.py"
SPEC = importlib.util.spec_from_file_location("alphamissense_annotation_header", HELPER)
CACHE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CACHE)
COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


class AlphaAnnotationHeader(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="alpha_header_", dir=os.environ.get("PRIVA_TEST_TMPDIR"))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.path = self.folder / "selected.annotation"

    def write(self, headers=None, *, compressed=True, columns=True):
        headers = headers if headers is not None else ['##VEP="v116.3" API="v116" assembly="GRCh38.p14"']
        content = "##fileformat=VCFv4.2\n" + "\n".join(headers) + "\n"
        if columns:
            content += COLUMNS + "chr1\t1\t.\tA\tC\t.\tPASS\tCSQ=C|missense_variant\n"
        self.path.write_bytes(gzip.compress(content.encode()) if compressed else content.encode())
        return self.path

    def check(self, *, assembly="hg38", major="116"):
        return CACHE.check_annotation(input_file=self.path, assembly=assembly, vep_major=major)

    def test_real_writer_style_minor_release_and_patch_assembly_preserve_source(self):
        self.write(['##VEP="v116.3" API="v116" time="2026-10-17 17:54:48" '
                    'cache="/cache/homo_sapiens_merged/116_GRCh38" ensembl=116.0d85231 '
                    'assembly="GRCh38.p14" refseq="GCF_000001405.40 - GCF_GRCh38.p14"',
                    "##VEP-command-line='vep --assembly GRCh38 --offline'"])
        before = CACHE.file_identity(self.path), self.path.read_bytes()
        self.assertEqual(self.check()["vep_version"], "116.3")
        self.assertEqual((CACHE.file_identity(self.path), self.path.read_bytes()), before)

    def test_both_configuration_assembly_aliases_match_the_header_family(self):
        for header, aliases in (("GRCh37.p13", ("hg19", "GRCh37")), ("GRCh38.p14", ("hg38", "GRCh38"))):
            self.write([f'##VEP="v116" assembly="{header}"'])
            for alias in aliases:
                with self.subTest(header=header, alias=alias):
                    self.assertEqual(self.check(assembly=alias)["assembly"], aliases[0])

    def test_new_requested_major_refuses_old_annotation_without_modifying_it(self):
        self.write()
        before = CACHE.file_identity(self.path), self.path.read_bytes()
        with self.assertRaisesRegex(ValueError, "differs from requested major 117"):
            self.check(major="117")
        self.assertEqual((CACHE.file_identity(self.path), self.path.read_bytes()), before)

    def test_recent_numeric_major_is_supported_without_a_fixed_release_list(self):
        self.write(['##VEP="v117.1" assembly="GRCh38"'])
        self.assertEqual(self.check(major="117")["vep_major"], "117")

    def test_wrong_selected_assembly_is_refused(self):
        self.write(['##VEP="v116.3" assembly="GRCh37.p13"'])
        with self.assertRaisesRegex(ValueError, "differs from requested"):
            self.check()

    def test_missing_release_or_assembly_cannot_establish_reuse(self):
        for headers in ([], ['##VEP="v116"'], ['##VEP="" assembly="GRCh38"'], ['##VEP="unknown" assembly="GRCh38"']):
            with self.subTest(headers=headers):
                self.write(headers)
                with self.assertRaises(ValueError):
                    self.check()

    def test_duplicate_identical_or_conflicting_release_headers_are_refused(self):
        for second in ("116.3", "117"):
            with self.subTest(second=second):
                self.write(['##VEP="v116.3" assembly="GRCh38"', f'##VEP="v{second}" assembly="GRCh38"'])
                with self.assertRaisesRegex(ValueError, "exactly one"):
                    self.check()

    def test_conflicting_or_duplicate_assembly_metadata_is_refused(self):
        for second in ("GRCh38", "GRCh37"):
            with self.subTest(second=second):
                self.write([f'##VEP="v116.3" assembly="GRCh38" assembly="{second}"'])
                with self.assertRaisesRegex(ValueError, "unambiguous"):
                    self.check()

    def test_plain_vcf_and_compression_magic_are_supported_independent_of_name(self):
        self.write(compressed=False)
        self.assertEqual(self.check()["vep_major"], "116")
        self.write(compressed=True)
        self.assertEqual(self.check()["vep_major"], "116")

    def test_header_limit_and_missing_column_header_fail_closed(self):
        self.write(['##padding=' + "x" * 256, '##VEP="v116" assembly="GRCh38"'])
        with mock.patch.object(CACHE, "ANNOTATION_HEADER_BYTES", 128):
            with self.assertRaisesRegex(ValueError, "bounded"):
                self.check()
        self.write(columns=False)
        with self.assertRaisesRegex(ValueError, "complete bounded"):
            self.check()

    def test_cli_reports_compatible_identity_and_nonzero_on_changed_major(self):
        self.write()
        command = [sys.executable, str(HELPER), "check-annotation", "--input", str(self.path),
                   "--assembly", "GRCh38", "--vep-major", "116"]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)["assembly"], "hg38")
        result = subprocess.run(command[:-1] + ["117"], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("differs from requested major", result.stderr)

    def test_changed_header_aborts_success_mismatch_and_parse_error(self):
        original = CACHE._check_vep_annotation_header
        for headers in (None, ['##VEP="v113" assembly="GRCh38"'], ['##VEP="unterminated']):
            with self.subTest(headers=headers):
                self.write(headers)
                def changed(*args):
                    try:
                        return original(*args)
                    finally:
                        self.path.write_bytes(self.path.read_bytes() + b'changed')
                with mock.patch.object(CACHE, '_check_vep_annotation_header', side_effect=changed):
                    with self.assertRaises(CACHE.AnnotationHeaderChangedError):
                        self.check()

    def test_disappearing_header_is_a_change_even_when_parsing_raises(self):
        self.write([])
        original = CACHE._check_vep_annotation_header
        def removed(*args):
            try:
                return original(*args)
            finally:
                self.path.unlink()
        with mock.patch.object(CACHE, '_check_vep_annotation_header', side_effect=removed):
            with self.assertRaises(CACHE.AnnotationHeaderChangedError):
                self.check()


if __name__ == "__main__":
    unittest.main()
