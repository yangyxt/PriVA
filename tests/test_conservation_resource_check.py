"""Actual small BigWig integrity/publication tests; no remote resource transfer."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


HERE = Path(__file__).resolve().parent
HELPER = HERE / "conservation_resource_check.py"
if not HELPER.is_file():
    HELPER = HERE.parent / "scripts/conservation_resource_check.py"
SPEC = importlib.util.spec_from_file_location("conservation_resource_check", HELPER)
guard = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(guard)
PERL = Path(sys.executable).parent / "perl"


def make_bigwig(folder, name, chromosomes):
    sizes = folder / (name + ".sizes")
    wig = folder / (name + ".wig")
    output = folder / (name + ".bw")
    sizes.write_text("".join(f"{chrom}\t{length}\n" for chrom, length in chromosomes.items()))
    wig.write_text("".join(f"variableStep chrom={chrom}\n1 1.5\n10 -0.5\n" for chrom in chromosomes))
    subprocess.run([str(PERL), "-MBio::DB::BigFile", "-e",
                    "Bio::DB::BigFile->createBigWig($ARGV[0], $ARGV[1], $ARGV[2], {});",
                    str(wig), str(sizes), str(output)], check=True, capture_output=True, text=True)
    return output


def make_reference(folder, assembly, chromosomes=None):
    chromosomes = chromosomes or {"chr1": guard.CHR1_LENGTH[guard.ASSEMBLIES[assembly]]}
    path = folder / (assembly + ".reference.fa")
    path.write_text(">chr1\nA\n")
    Path(str(path) + ".fai").write_text("".join(
        f"{name}\t{length}\t6\t1\t2\n" for name, length in chromosomes.items()))
    return path


class ConservationIntegrity(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_conservation_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.reference = make_reference(self.folder, "hg38", {"chr1": 248956422, "chr2": 242193529})
        self.track = make_bigwig(self.folder, "valid", {"1": 248956422, "2": 242193529})

    def test_actual_bigwig_and_reference_chromosome_aliases_validate(self):
        record = guard.validate(self.track, "GRCh38", self.reference)
        self.assertEqual(record["shared_reference_chromosomes"], ["1", "2"])
        self.assertEqual(record["sha256"], hashlib.sha256(self.track.read_bytes()).hexdigest())
        self.assertEqual(record["reader"], "Bio::DB::BigFile")
        self.assertIsNone(record["release"])
        self.assertIsNone(record["source_url"])
        self.assertFalse(record["scientific_qualified"])

    def test_actual_hg19_bigwig_is_accepted_without_changing_its_source_contract(self):
        reference = make_reference(self.folder, "hg19")
        track = make_bigwig(self.folder, "hg19", {"chr1": 249250621})
        self.assertEqual(guard.validate(track, "GRCh37", reference)["assembly"], "GRCh37")

    def test_corrupt_text_is_rejected_and_retained(self):
        path = self.folder / "corrupt.bw"
        original = b"nonempty successful HTTP body is not a BigWig\n"
        path.write_bytes(original)
        with self.assertRaisesRegex(ValueError, "Cannot validate Conservation BigWig"):
            guard.validate(path, "hg38", self.reference)
        self.assertEqual(path.read_bytes(), original)

    def test_truncated_actual_bigwig_is_rejected_and_retained(self):
        path = self.folder / "truncated.bw.partial"
        original = self.track.read_bytes()[:len(self.track.read_bytes()) // 2]
        path.write_bytes(original)
        with self.assertRaisesRegex(ValueError, "Cannot validate Conservation BigWig"):
            guard.validate(path, "hg38", self.reference)
        self.assertEqual(path.read_bytes(), original)

    def test_wrong_track_assembly_is_rejected(self):
        track = make_bigwig(self.folder, "wrong_assembly", {"chr1": 249250621})
        with self.assertRaisesRegex(ValueError, "BigWig chr1 length does not identify GRCh38"):
            guard.validate(track, "hg38", self.reference)

    def test_wrong_reference_assembly_is_rejected(self):
        reference = make_reference(self.folder, "hg19")
        with self.assertRaisesRegex(ValueError, "Configured reference chr1 length does not identify GRCh38"):
            guard.validate(self.track, "hg38", reference)

    def test_every_shared_chromosome_length_must_match_reference(self):
        track = make_bigwig(self.folder, "wrong_chr2", {"chr1": 248956422, "chr2": 242193528})
        with self.assertRaisesRegex(ValueError, "chromosome length mismatch for 2"):
            guard.validate(track, "hg38", self.reference)

    def test_extra_track_contigs_are_reported_without_changing_reference_scope(self):
        track = make_bigwig(self.folder, "extra", {"chr1": 248956422, "oldTrackScaffold": 1000})
        record = guard.validate(track, "hg38", self.reference)
        self.assertEqual(record["extra_non_reference_chromosomes"], ["oldTrackScaffold"])
        self.assertEqual(record["shared_reference_chromosomes"], ["1"])

    def test_missing_reference_or_index_is_rejected_without_substitution(self):
        for path in (self.folder / "missing.fa", self.reference):
            with self.subTest(path=path):
                if path == self.reference:
                    Path(str(path) + ".fai").unlink()
                with self.assertRaisesRegex(ValueError, "configured reference"):
                    guard.validate(self.track, "hg38", path)

    def test_known_default_publisher_pin_is_exact_and_size_is_enforced(self):
        self.assertEqual(guard.GRCH38_DEFAULT_PIN, {"bytes": 9611134022, "bsd_sum": [30633, 9385874]})
        with self.assertRaisesRegex(ValueError, "publisher byte count mismatch"):
            guard.validate(self.track, "hg38", self.reference, guard.GRCH38_DEFAULT_URL)

    def test_known_default_bsd_pin_is_enforced_after_real_bigwig_validation(self):
        measured = guard.bsd_sum(self.track)
        pin = {"bytes": self.track.stat().st_size, "bsd_sum": [measured[0] ^ 1, measured[1]]}
        with patch.object(guard, "GRCH38_DEFAULT_PIN", pin), self.assertRaisesRegex(ValueError, "BSD checksum mismatch"):
            guard.validate(self.track, "hg38", self.reference, guard.GRCH38_DEFAULT_URL)
        pin["bsd_sum"] = measured
        with patch.object(guard, "GRCH38_DEFAULT_PIN", pin):
            record = guard.validate(self.track, "hg38", self.reference, guard.GRCH38_DEFAULT_URL)
        self.assertEqual(record["release"], 114)
        self.assertEqual(record["publisher_bsd_sum"], measured)
        self.assertFalse(record["publisher_cryptographic_checksum_claimed"])

    def test_configured_custom_identity_is_not_inferred_from_default_source_setting(self):
        # Existing-file validation deliberately receives no download URL.
        record = guard.validate(self.track, "hg38", self.reference)
        self.assertEqual(record["source_kind"], "configured_existing_identity")
        self.assertIsNone(record["release"])
        self.assertIsNone(record["publisher_bsd_sum"])

    def test_publication_failure_retains_verified_partial(self):
        destination = self.folder / "publish_failure.bw"
        before = self.track.read_bytes()
        with patch.object(guard.os, "link", side_effect=OSError("fixture publication failure")), \
                self.assertRaisesRegex(OSError, "publication failure"):
            guard.publish(self.track, destination, "hg38", self.reference, "https://example.invalid/custom.bw")
        self.assertEqual(self.track.read_bytes(), before)
        self.assertFalse(destination.exists())

    def test_existing_publication_target_is_never_overwritten(self):
        destination = self.folder / "existing.bw"
        destination.write_bytes(b"retained prior target")
        with self.assertRaises(FileExistsError):
            guard.publish(self.track, destination, "hg38", self.reference, "https://example.invalid/custom.bw")
        self.assertTrue(self.track.is_file())
        self.assertEqual(destination.read_bytes(), b"retained prior target")

    def test_successful_publication_preserves_verified_bytes(self):
        destination = self.folder / "published.bw"
        before = self.track.read_bytes()
        record = guard.publish(self.track, destination, "hg38", self.reference, "https://example.invalid/custom.bw")
        self.assertFalse(self.track.exists())
        self.assertEqual(destination.read_bytes(), before)
        self.assertEqual(record["path"], str(destination))
        self.assertEqual(record["sha256"], hashlib.sha256(before).hexdigest())
        self.assertIsNone(record["release"])

    def test_cli_reports_corrupt_input_without_publication(self):
        path = self.folder / "invalid.partial"
        path.write_text("not a BigWig\n")
        destination = self.folder / "invalid.final.bw"
        result = subprocess.run([sys.executable, str(HELPER), "--input", str(path), "--assembly", "hg38",
                                 "--reference", str(self.reference), "--source-url", "https://example.invalid/custom.bw",
                                 "--publish-to", str(destination)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Conservation validation failed", result.stderr)
        self.assertTrue(path.exists())
        self.assertFalse(destination.exists())


if __name__ == "__main__":
    unittest.main()
