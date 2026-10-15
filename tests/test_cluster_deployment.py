"""Installation integrity tests; no network, solver or production caches."""
import hashlib
import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location(
    "resources", Path(__file__).resolve().parents[1] / "scripts/deploy_cluster_resources.py")
resources = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(resources)


class DownloadIntegrity(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.destination = Path(self.temporary.name) / "cache.bin"
        self.payload = b"expected annotation cache\n"
        self.md5 = hashlib.md5(self.payload).hexdigest()

    def fake_transfer(self, command):
        Path(command[command.index("--output") + 1]).write_bytes(self.payload)

    def test_verified_partial_is_published(self):
        with patch.object(resources, "run", self.fake_transfer):
            resources.download("https://example.invalid/data", self.destination, self.md5)
        self.assertEqual(self.destination.read_bytes(), self.payload)
        self.assertFalse(Path(str(self.destination) + ".partial").exists())

    def test_corrupt_transfer_is_not_published(self):
        with patch.object(resources, "run", self.fake_transfer), self.assertRaises(RuntimeError):
            resources.download("https://example.invalid/data", self.destination, "0" * 32)
        self.assertFalse(self.destination.exists())
        self.assertTrue(Path(str(self.destination) + ".partial").exists())

    def test_corrupt_existing_file_is_not_overwritten(self):
        self.destination.write_bytes(b"different cache identity")
        with patch.object(resources, "run") as transfer, self.assertRaises(RuntimeError):
            resources.download("https://example.invalid/data", self.destination, self.md5)
        transfer.assert_not_called()
        self.assertEqual(self.destination.read_bytes(), b"different cache identity")

    def test_valid_existing_file_never_downloads(self):
        self.destination.write_bytes(self.payload)
        with patch.object(resources, "run") as transfer:
            resources.download("https://example.invalid/data", self.destination, self.md5)
        transfer.assert_not_called()

    def test_unexpected_size_is_not_published(self):
        with patch.object(resources, "run", self.fake_transfer), self.assertRaises(RuntimeError):
            resources.download("https://example.invalid/data", self.destination, expected_size=42)
        self.assertFalse(self.destination.exists())


class ReuseIntegrity(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.args = SimpleNamespace(cache_root=self.root / "cache", state_dir=self.root / "state",
                                    software_root=self.root / "software", tmp_root=self.root / "tmp",
                                    samtools="samtools")
        self.args.state_dir.mkdir()

    def test_changed_existing_fasta_is_rejected_against_archive(self):
        folder = self.args.cache_root / "references"
        folder.mkdir(parents=True)
        original = b">chr1\nACGT\n"
        (folder / "ucsc.hg19.fasta").write_bytes(b">chr1\nAAAA\n")
        archive = folder / "ucsc.hg19.fasta.gz"
        with gzip.open(archive, "wb") as stream:
            stream.write(original)
        manifest = {"reference": {"url": "https://example.invalid/ref", "md5": resources.checksum(archive, "md5")}}
        with self.assertRaisesRegex(RuntimeError, "does not match"):
            resources.reference(self.args, manifest)

    def test_unidentified_existing_image_is_rejected(self):
        folder = self.args.software_root / "cadd-1.7.2"
        folder.mkdir(parents=True)
        (folder / "cadd-upstream-0.1.1.sif").write_bytes(b"unverified image")
        with patch.object(resources, "source", return_value=self.root), self.assertRaisesRegex(RuntimeError, "no verified"):
            resources.cadd_runtime(self.args, {"cadd": {"image": "docker://example@sha256:abc"}})

    def test_mutated_existing_image_is_rejected(self):
        folder = self.args.software_root / "cadd-1.7.2"
        folder.mkdir(parents=True)
        image = folder / "cadd-upstream-0.1.1.sif"
        image.write_bytes(b"original image")
        origin = {"image_source": "docker://example@sha256:abc", "image_sha256": resources.checksum(image)}
        resources.write_json(folder / "image.origin.json", origin)
        image.write_bytes(b"changed image")
        with patch.object(resources, "source", return_value=self.root), self.assertRaisesRegex(RuntimeError, "differs"):
            resources.cadd_runtime(self.args, {"cadd": {"image": origin["image_source"]}})


if __name__ == "__main__":
    unittest.main()
