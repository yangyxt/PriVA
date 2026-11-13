"""Installation integrity tests; no network, solver or production caches."""
import hashlib
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch
from urllib.parse import urlparse

SPEC = importlib.util.spec_from_file_location(
    "resources", Path(__file__).resolve().parents[1] / "scripts/deploy_cluster_resources.py")
resources = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(resources)
READINESS_SPEC = importlib.util.spec_from_file_location(
    "readiness", Path(__file__).resolve().parents[1] / "scripts/check_cluster_readiness.py")
readiness = importlib.util.module_from_spec(READINESS_SPEC)
READINESS_SPEC.loader.exec_module(readiness)


class CaddPrescoreIdentity(unittest.TestCase):
    def setUp(self):
        self.manifest = json.loads((resources.REPO / "envs/cluster-resources.hg19.json").read_text())
        self.cadd = self.manifest["cadd"]

    def test_washington_grch37_v17_data_and_index_identities_are_pinned(self):
        # Data hashes/URLs agree with historical PriVA 09223ef config.yaml;
        # index hashes identify the matching Washington publisher pair.
        self.assertEqual(self.manifest["assembly"], "GRCh37")
        self.assertEqual(self.cadd["source_commit"], "1fead0c1e02d3f0a8599bdb04e6b3a2dd70b4866")
        self.assertEqual(self.cadd["prescore_url_root"],
                         "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/")
        self.assertEqual(self.cadd["prescores"], {
            "whole_genome_SNVs_inclAnno.tsv.gz": "1b89cfd7da4e2a8f336c1ec33a6f44eb",
            "whole_genome_SNVs_inclAnno.tsv.gz.tbi": "b3f770e0ba669e23c65a109b8221a1b5",
            "gnomad.genomes-exomes.r4.0.indel_inclAnno.tsv.gz": "3d88f8fde9df07060b4f8fe07a0d5d5d",
            "gnomad.genomes-exomes.r4.0.indel_inclAnno.tsv.gz.tbi": "c583ed9f852c9896081205c3523f9f23",
        })

    def test_prescore_urls_and_checksum_source_share_publisher_and_release(self):
        self.assertEqual(self.cadd["checksum_source"], self.cadd["prescore_url_root"] + "MD5SUMs")
        checksum_source = urlparse(self.cadd["checksum_source"])
        for filename in self.cadd["prescores"]:
            with self.subTest(filename=filename):
                download_source = urlparse(self.cadd["prescore_url_root"] + filename)
                self.assertEqual(download_source.netloc, checksum_source.netloc)
                self.assertEqual(download_source.path.rsplit("/", 1)[0],
                                 checksum_source.path.rsplit("/", 1)[0])


class IndependentLicensedReadiness(unittest.TestCase):
    def test_imported_spliceai_does_not_clear_missing_primateai(self):
        records = [{"key": key, "status": "present_not_scientifically_qualified"}
                   for key in ("spliceai_snv_prescore", "spliceai_indel_prescore")]
        records.append({"key": "primateai_prescore", "status": "missing"})
        blockers = readiness.licensed_blockers(records)
        self.assertEqual([item["resource"] for item in blockers], ["PrimateAI v0.2"])
        self.assertEqual(blockers[0]["unavailable_config_keys"], ["primateai_prescore"])

    def test_missing_indel_keeps_spliceai_blocked(self):
        records = [{"key": key, "status": "present_not_scientifically_qualified"}
                   for key in ("spliceai_snv_prescore", "primateai_prescore")]
        blockers = readiness.licensed_blockers(records)
        self.assertEqual([item["resource"] for item in blockers], ["SpliceAI v1"])
        self.assertEqual(blockers[0]["unavailable_config_keys"], ["spliceai_indel_prescore"])

    def test_all_present_only_clears_installation_blockers(self):
        records = [{"key": key, "status": "present_not_scientifically_qualified"}
                   for key in ("spliceai_snv_prescore", "spliceai_indel_prescore", "primateai_prescore")]
        self.assertEqual(readiness.licensed_blockers(records), [])
        self.assertTrue(all(item["status"] == "present_not_scientifically_qualified" for item in records))


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

    def test_completed_verified_partial_is_published_without_http(self):
        partial = Path(str(self.destination) + ".partial")
        partial.write_bytes(self.payload)
        with patch.object(resources, "run") as transfer:
            result = resources.download("https://example.invalid/data", self.destination, self.md5,
                                        expected_size=len(self.payload))
        transfer.assert_not_called()
        self.assertEqual(result, self.destination)
        self.assertEqual(self.destination.read_bytes(), self.payload)
        self.assertFalse(partial.exists())

    def test_incomplete_partial_resumes_and_verifies_before_publication(self):
        partial = Path(str(self.destination) + ".partial")
        partial.write_bytes(self.payload[:5])
        with patch.object(resources, "run", side_effect=self.fake_transfer) as transfer:
            resources.download("https://example.invalid/data", self.destination, self.md5)
        transfer.assert_called_once()
        self.assertIn("--continue-at", transfer.call_args.args[0])
        self.assertEqual(self.destination.read_bytes(), self.payload)
        self.assertFalse(partial.exists())

    def test_corrupt_existing_partial_is_not_promoted(self):
        partial = Path(str(self.destination) + ".partial")
        corrupt = b"x" * len(self.payload)
        partial.write_bytes(corrupt)
        with patch.object(resources, "run") as transfer, self.assertRaisesRegex(RuntimeError, "Checksum mismatch"):
            resources.download("https://example.invalid/data", self.destination, self.md5)
        transfer.assert_called_once()
        self.assertFalse(self.destination.exists())
        self.assertEqual(partial.read_bytes(), corrupt)

    def test_verified_partial_with_wrong_size_is_not_promoted(self):
        partial = Path(str(self.destination) + ".partial")
        partial.write_bytes(self.payload)
        with patch.object(resources, "run") as transfer, self.assertRaisesRegex(RuntimeError, "Unexpected size"):
            resources.download("https://example.invalid/data", self.destination, self.md5,
                               expected_size=len(self.payload) + 1)
        transfer.assert_not_called()
        self.assertFalse(self.destination.exists())
        self.assertEqual(partial.read_bytes(), self.payload)

    def test_size_alone_does_not_skip_partial_transfer(self):
        partial = Path(str(self.destination) + ".partial")
        partial.write_bytes(b"x" * len(self.payload))
        with patch.object(resources, "run", side_effect=self.fake_transfer) as transfer:
            resources.download("https://example.invalid/data", self.destination,
                               expected_size=len(self.payload))
        transfer.assert_called_once()
        self.assertEqual(self.destination.read_bytes(), self.payload)

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

    def test_existing_final_is_rechecked_on_every_reuse(self):
        self.destination.write_bytes(self.payload)
        with patch.object(resources, "run") as transfer:
            resources.download("https://example.invalid/data", self.destination, self.md5)
            corrupt = b"x" * len(self.payload)
            self.destination.write_bytes(corrupt)
            with self.assertRaisesRegex(RuntimeError, "Checksum mismatch"):
                resources.download("https://example.invalid/data", self.destination, self.md5)
        transfer.assert_not_called()
        self.assertEqual(self.destination.read_bytes(), corrupt)

    def test_existing_final_with_wrong_size_is_rejected(self):
        self.destination.write_bytes(self.payload)
        with patch.object(resources, "run") as transfer, self.assertRaisesRegex(RuntimeError, "Unexpected size"):
            resources.download("https://example.invalid/data", self.destination, self.md5,
                               expected_size=len(self.payload) + 1)
        transfer.assert_not_called()
        self.assertEqual(self.destination.read_bytes(), self.payload)

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


if __name__ == "__main__":
    unittest.main()
