"""Actual CDS reader I/O: identical sequence results across safe cache-name aliases."""
import gzip
import hashlib
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from alternative_start_codon import AlternativeStartCodonDetector


FASTA = b">ENST00000000233.1 cds chromosome:GRCh37:1:1:18:1\nATGAAAATGCCCATGTAA\n"


class CDSInputIO(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_cds_io_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)

    @staticmethod
    def write_fasta(path):
        path.parent.mkdir(parents=True, exist_ok=True)
        if str(path).endswith(".gz"):
            with gzip.open(path, "wb") as handle:
                handle.write(FASTA)
        else:
            path.write_bytes(FASTA)

    def test_standard_basename_keeps_exact_existing_cache_path_and_reuses_without_parsing(self):
        selected = self.folder / "standard.fa.gz"
        self.write_fasta(selected)
        before = hashlib.sha256(selected.read_bytes()).hexdigest()
        detector = AlternativeStartCodonDetector(str(selected))
        self.assertEqual(detector.cache_pkl_path, str(self.folder / "standard.cds_cache.pkl.gz"))
        expected = detector.find_alternative_start("ENST00000000233")
        with patch.object(AlternativeStartCodonDetector, "_load_fasta", side_effect=AssertionError("must reuse completed cache")):
            reused = AlternativeStartCodonDetector(str(selected))
        self.assertEqual(reused.find_alternative_start("ENST00000000233"), expected)
        self.assertEqual(hashlib.sha256(selected.read_bytes()).hexdigest(), before)

    def test_plain_generic_gzip_parent_token_and_spaces_preserve_bytes_and_actual_reader_results(self):
        standard = self.folder / "standard.fa.gz"
        self.write_fasta(standard)
        baseline = AlternativeStartCodonDetector(str(standard))
        expected = baseline.find_alternative_start("ENST00000000233")
        cases = ("plain.fa", "generic.other.gz", "selected input's.fa", "selected input's.fa.gz", "selected input's.other.gz",
                 "refs.fa.gz/plain.fa", "refs.fa.gz/generic.other.gz", "refs.fa.gz/standard.fa.gz")
        for relative in cases:
            with self.subTest(relative=relative):
                selected = self.folder / relative
                self.write_fasta(selected)
                before = selected.read_bytes()
                detector = AlternativeStartCodonDetector(str(selected))
                expected_cache = (str(selected)[:-6] + ".cds_cache.pkl.gz"
                                  if selected.name.endswith(".fa.gz") else str(selected) + ".cds_cache.pkl.gz")
                self.assertEqual(detector.cache_pkl_path, expected_cache)
                self.assertNotEqual(detector.cache_pkl_path, str(selected))
                self.assertEqual(detector.cds_sequences, baseline.cds_sequences)
                self.assertEqual(detector.transcript_to_key, baseline.transcript_to_key)
                self.assertEqual(detector.find_alternative_start("ENST00000000233"), expected)
                self.assertEqual(selected.read_bytes(), before)
                with patch.object(AlternativeStartCodonDetector, "_load_fasta", side_effect=AssertionError("must reuse completed cache")):
                    reused = AlternativeStartCodonDetector(str(selected))
                self.assertEqual(reused.find_alternative_start("ENST00000000233"), expected)
                self.assertEqual(selected.read_bytes(), before)

    def test_explicit_cache_path_remains_exact(self):
        selected = self.folder / "selected.fa.gz"
        explicit_cache = self.folder / "explicit_cache.pkl.gz"
        self.write_fasta(selected)
        detector = AlternativeStartCodonDetector(str(selected), str(explicit_cache))
        self.assertEqual(detector.cache_pkl_path, str(explicit_cache))
        self.assertTrue(explicit_cache.is_file())


if __name__ == "__main__":
    unittest.main()
