"""Installation integrity tests; no network, solver or production caches."""
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest
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


if __name__ == "__main__":
    unittest.main()
