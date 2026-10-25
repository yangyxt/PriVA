"""Real BGZF/tabix installation fixtures; no full catalogue or scientific run."""
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
SPEC = importlib.util.spec_from_file_location("alphamissense_cache", ROOT / "scripts/alphamissense_cache.py")
CACHE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CACHE)


class AlphaMissenseCacheTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="alphamissense_cache_",
                                                dir=os.environ.get("PRIVA_TEST_TMPDIR"))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.cache_dir = self.folder / "managed plugins"
        self.raw = self.folder / "AlphaMissense_hg38.tsv.gz"
        content = ("# Copyright fixture\n" + CACHE.HEADER + "\n"
                   "chr1\t11\tA\tC\thg38\tP1\tENST1\tA1C\t0.2500\tambiguous\n"
                   "chr1\t12\tG\tT\thg38\tP1\tENST1\tG2T\t0.001\tbenign\n")
        self.content = content.encode()
        with self.raw.open("wb") as output:
            subprocess.run(["bgzip", "-c"], input=self.content, stdout=output, check=True)
        self.original = self.raw.read_bytes()
        self.pin = dict(bytes=len(self.original), md5=CACHE.checksum(self.raw), generation="fixture-generation")
        patcher = mock.patch.dict(CACHE.SOURCES, {"hg38": self.pin})
        patcher.start()
        self.addCleanup(patcher.stop)
        self.destination = self.cache_dir / "AlphaMissense/AlphaMissense_hg38.tsv.gz"

    def install(self, source=None):
        return CACHE.install(cache_dir=self.cache_dir, assembly="GRCh38", source=source)

    def download(self, path, assembly):
        self.assertEqual(assembly, "hg38")
        path.write_bytes(self.original)
        return "https://fixture.invalid/original"

    def test_supplied_original_is_adopted_as_managed_bytes_and_reindexed(self):
        before = CACHE.file_identity(self.raw)
        with mock.patch.object(CACHE, "download", side_effect=AssertionError("unexpected network")):
            result = self.install(self.raw)
        self.assertEqual(result, self.destination)
        self.assertEqual(result.read_bytes(), self.original)
        self.assertFalse(result.is_symlink())
        self.assertEqual(CACHE.file_identity(self.raw), before)
        self.assertFalse(Path(str(self.raw) + ".tbi").exists())
        self.assertEqual(CACHE.check_index(result, "hg38"), ["chr1"])
        receipt = json.loads((result.parent / "AlphaMissense_hg38.ready.json").read_text())
        self.assertEqual(receipt["publisher"], self.pin)
        self.assertEqual(receipt["supplied_source"], str(self.raw))
        self.assertFalse(receipt["scientific_acceptance"])

    def test_selected_assembly_download_is_automatic_and_no_other_assembly_is_created(self):
        with mock.patch.object(CACHE, "download", side_effect=self.download) as download:
            self.install()
        download.assert_called_once()
        self.assertFalse((self.destination.parent / "AlphaMissense_hg19.tsv.gz").exists())
        self.assertEqual(self.destination.read_bytes(), self.original)

    def test_invalid_supplied_source_is_never_replaced_or_silently_downloaded(self):
        self.raw.write_bytes(b"invalid supplied bytes")
        with mock.patch.object(CACHE, "download", side_effect=AssertionError("unexpected fallback")):
            with self.assertRaisesRegex(ValueError, "size mismatch"):
                self.install(self.raw)
        self.assertEqual(self.raw.read_bytes(), b"invalid supplied bytes")
        self.assertFalse(self.destination.exists())

    def test_wrong_assembly_source_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "size mismatch"):
            CACHE.install(cache_dir=self.cache_dir, assembly="hg19", source=self.raw)
        self.assertEqual(self.raw.read_bytes(), self.original)

    def test_complete_receipt_reuses_index_without_revalidating_large_bytes(self):
        self.install(self.raw)
        before = CACHE.file_identity(self.destination)
        with mock.patch.object(CACHE, "validate", side_effect=AssertionError("completed source reread")), \
                mock.patch.object(CACHE, "download", side_effect=AssertionError("unexpected network")):
            self.install()
        self.assertEqual(CACHE.file_identity(self.destination), before)

    def test_bad_index_is_rebuilt_and_raw_identity_is_unchanged(self):
        self.install(self.raw)
        before = CACHE.file_identity(self.destination)
        index = Path(str(self.destination) + ".tbi")
        index.write_bytes(b"bad index")
        with mock.patch.object(CACHE, "download", side_effect=AssertionError("unexpected network")):
            self.install()
        self.assertEqual(CACHE.file_identity(self.destination), before)
        self.assertEqual(CACHE.check_index(self.destination, "hg38"), ["chr1"])
        self.assertEqual(len(list(index.parent.glob(index.name + ".invalid-*"))), 1)

    def test_invalid_managed_payload_is_retained_while_original_is_recovered(self):
        self.destination.parent.mkdir(parents=True)
        self.destination.write_bytes(b"interrupted old final")
        with mock.patch.object(CACHE, "download", side_effect=self.download):
            self.install()
        retained = list(self.destination.parent.glob(self.destination.name + ".invalid-*"))
        self.assertEqual([path.read_bytes() for path in retained], [b"interrupted old final"])
        self.assertEqual(self.destination.read_bytes(), self.original)

    def test_old_symlink_view_is_adopted_without_writing_its_external_target(self):
        self.destination.parent.mkdir(parents=True)
        self.destination.symlink_to(self.raw)
        before = CACHE.file_identity(self.raw)
        self.install()
        self.assertFalse(self.destination.is_symlink())
        self.assertEqual(CACHE.file_identity(self.raw), before)
        self.assertFalse(Path(str(self.raw) + ".tbi").exists())

    def test_interrupted_download_keeps_partial_and_no_ready_receipt(self):
        def fail(path, assembly):
            path.write_bytes(self.original[:20])
            raise subprocess.CalledProcessError(18, ["curl"])
        with mock.patch.object(CACHE, "download", side_effect=fail):
            with self.assertRaises(subprocess.CalledProcessError):
                self.install()
        self.assertFalse(self.destination.exists())
        self.assertFalse((self.destination.parent / "AlphaMissense_hg38.ready.json").exists())
        self.assertEqual(Path(str(self.destination) + ".partial").read_bytes(), self.original[:20])

    def test_download_resumes_the_pinned_generation_and_retains_corrupt_attempt(self):
        partial = self.folder / "transfer.partial"
        partial.write_bytes(self.original[:20])
        def resume(command, **kwargs):
            self.assertIn("--continue-at", command)
            self.assertTrue(command[-1].endswith("?generation=fixture-generation"))
            self.assertEqual(Path(command[command.index("--output") + 1]).read_bytes(), self.original[:20])
            partial.write_bytes(self.original)
        with mock.patch.object(CACHE, "run", side_effect=resume) as run:
            CACHE.download(partial, "hg38")
        run.assert_called_once()
        partial.write_bytes(b"x" * len(self.original))
        def clean(command, **kwargs):
            self.assertFalse(partial.exists())
            partial.write_bytes(self.original)
        with mock.patch.object(CACHE, "run", side_effect=clean):
            CACHE.download(partial, "hg38")
        self.assertEqual(len(list(self.folder.glob("transfer.partial.invalid-*"))), 1)

    def test_index_failure_never_promotes_raw_or_receipt(self):
        original_run = CACHE.run
        def fail(command, **kwargs):
            if command[0] == "tabix":
                raise subprocess.CalledProcessError(2, command)
            return original_run(command, **kwargs)
        with mock.patch.object(CACHE, "run", side_effect=fail):
            with self.assertRaises(subprocess.CalledProcessError):
                self.install(self.raw)
        self.assertFalse(self.destination.exists())
        self.assertFalse((self.destination.parent / "AlphaMissense_hg38.ready.json").exists())
        self.assertTrue(Path(str(self.destination) + ".partial").is_file())


if __name__ == "__main__":
    unittest.main()
