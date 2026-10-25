"""Physical prepared-v2 adoption: independent bytes and unchanged scientific gates."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import pysam

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import spliceai_cache_input as cache


class PreparedSpliceAI2Relocation(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="spliceai2_relocation_",
            dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.source = self.folder / "original"
        self.source.mkdir()
        self.root = self.folder / "managed"
        self.target = self.root / "VEP_plugin_cache/SpliceAI2/GRCh38"
        self.outputs = {}
        for role in ("snv", "indel"):
            path = self.source / f"SpliceAI2_GRCh38.{role}.vep.vcf.gz"
            with pysam.BGZFile(str(path), "wb") as handle:
                handle.write(("##fileformat=VCFv4.2\n##SpliceAI2Version=2\n##assembly=GRCh38\n"
                              "##contig=<ID=chr1,length=1000>\n"
                              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                              "chr1\t10\t.\tA\tC\t.\t.\t.\n").encode())
            pysam.tabix_index(str(path), preset="vcf", force=True)
            index = Path(str(path) + ".tbi")
            self.outputs[role] = {"path": str(path), "bytes": path.stat().st_size,
                "sha256": self.sha(path), "index_path": str(index), "index_sha256": self.sha(index)}
        context = self.source / "SpliceAI2_GRCh38.gene_context.tsv"
        context.write_text("chrom\tstart\tend\tstrand\tgene\nchr1\t1\t100\t+\tGENE\n")
        self.data = {"schema_version": 1, "resource": "SpliceAI2", "version": "2.0", "assembly": "GRCh38",
            "status": "prepared_annotation_only", "coverage": "complete_published", "scientific_qualified": False,
            "source_revision": cache.REVISION, "source_manifest_sha256": "a" * 64, "outputs": self.outputs,
            "gene_context": {"path": str(context), "sha256": self.sha(context), "scope": "unchanged gene association"},
            "projection": "unchanged score strings, missingness and row order", "thresholds": {"balanced": 0.25},
            "converter": {"path": "/historical/production/priva", "sha256": "b" * 64},
            "source_files": [{"path": "/historical/publisher/shard.tsv.gz", "sha256": "c" * 64}]}
        self.receipt = self.source / "adapter.ready.json"
        self.save_receipt()
        (self.source / ".install.lock").touch()
        self.originals = {path.name: (path.read_bytes(), cache._identity(path))
                          for path in self.source.iterdir()}

    @staticmethod
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def save_receipt(self):
        self.receipt.write_text(json.dumps(self.data, indent=2) + "\n")

    def adopt(self):
        return cache.adopt_completed_projection(source_dir=self.source, cache_dir=self.target, cache_root=self.root)

    def paths(self, folder):
        return [folder / Path(self.outputs[role]["path"]).name for role in ("snv", "indel")]

    def assert_originals_unchanged(self):
        for name, (data, before) in self.originals.items():
            self.assertEqual((self.source / name).read_bytes(), data)
            self.assertEqual(cache._identity(self.source / name), before)

    def test_copy_is_independent_and_receipt_changes_only_locations_and_provenance(self):
        result = self.adopt()
        self.assertEqual(result["status"], "complete")
        self.assertFalse(result["converter_executed"])
        self.assert_originals_unchanged()
        self.assertEqual((self.target / "adapter.original.ready.json").read_bytes(), self.receipt.read_bytes())
        copied = json.loads((self.target / "adapter.ready.json").read_text())
        expected = json.loads(json.dumps(self.data))
        for role in ("snv", "indel"):
            for field in ("path", "index_path"):
                expected["outputs"][role][field] = str(self.target / Path(expected["outputs"][role][field]).name)
        expected["gene_context"]["path"] = str(self.target / Path(expected["gene_context"]["path"]).name)
        expected["relocation"] = copied["relocation"]
        self.assertEqual(copied, expected)
        self.assertIs(copied["scientific_qualified"], False)
        for item in result["members"]:
            source, target = Path(item["source"]), Path(item["destination"])
            self.assertTrue(target.is_relative_to(self.root))
            self.assertFalse(target.is_symlink())
            self.assertNotEqual(source.stat().st_ino, target.stat().st_ino)
            self.assertEqual(self.sha(source), self.sha(target))
        with cache._completed_projection(*self.paths(self.target)):
            pass

    def test_relocated_publication_validates_without_original_source_paths(self):
        self.adopt()
        self.source.rename(self.folder / "retained_original_elsewhere")
        with cache._completed_projection(*self.paths(self.target)):
            pass

    def test_existing_matching_destination_reuses_without_mutation(self):
        self.adopt()
        before = {path.name: cache._identity(path) for path in self.target.iterdir()}
        result = self.adopt()
        self.assertEqual(result["status"], "reused")
        self.assertEqual({path.name: cache._identity(path) for path in self.target.iterdir()}, before)
        self.assert_originals_unchanged()

    def test_explicit_root_directory_alias_has_canonical_receipts_and_valid_alias_paths(self):
        self.root.mkdir()
        alias = self.folder / "selected-cache-alias"
        alias.symlink_to(self.root, target_is_directory=True)
        planned = alias / self.target.relative_to(self.root)
        result = cache.adopt_completed_projection(source_dir=self.source, cache_dir=planned, cache_root=alias)
        self.assertEqual(result["destination"], str(self.target))
        with cache._completed_projection(*self.paths(planned)):
            pass
        self.assert_originals_unchanged()

    def test_missing_source_lock_fails_without_creating_files_in_original_folder(self):
        lock = self.source / ".install.lock"
        lock.unlink()
        before = cache._identity(self.source)
        with self.assertRaisesRegex(ValueError, "must already contain its installation lock"):
            self.adopt()
        self.assertFalse(lock.exists())
        self.assertEqual(cache._identity(self.source), before)
        self.assertFalse(self.root.exists())

    def interrupt(self):
        original = cache._copy_verified_member
        def copy(source, target, checksum):
            if ".indel.vep.vcf.gz" in source.name and not source.name.endswith(".tbi"):
                target.with_name(target.name + ".copying").write_bytes(source.read_bytes()[:19])
                raise OSError("fixture copy interrupted")
            return original(source, target, checksum)
        with mock.patch.object(cache, "_copy_verified_member", side_effect=copy):
            with self.assertRaisesRegex(OSError, "fixture copy interrupted"):
                self.adopt()

    def test_owned_partial_copy_resumes_without_publishing_incomplete_cache(self):
        self.interrupt()
        self.assertFalse(self.target.exists())
        stage, = self.target.parent.glob(".GRCh38.adoption-*")
        self.assertFalse((stage / "adapter.ready.json").exists())
        self.assertTrue(any(stage.glob("*.copying")))
        self.adopt()
        self.assertTrue((self.target / "relocation.ready.json").exists())
        self.assert_originals_unchanged()

    def test_corrupt_partial_prefix_fails_without_discarding_evidence(self):
        self.interrupt()
        stage, = self.target.parent.glob(".GRCh38.adoption-*")
        partial, = stage.glob("*.copying")
        partial.write_bytes(b"incorrect prefix")
        with self.assertRaisesRegex(ValueError, "partial prefix differs"):
            self.adopt()
        self.assertEqual(partial.read_bytes(), b"incorrect prefix")
        self.assertFalse(self.target.exists())
        self.assert_originals_unchanged()

    def test_linked_staging_metadata_cannot_write_through_to_original_receipt(self):
        self.interrupt()
        stage, = self.target.parent.glob(".GRCh38.adoption-*")
        (stage / "adapter.original.ready.json").symlink_to(self.receipt)
        with self.assertRaisesRegex(ValueError, "linked adoption metadata"):
            self.adopt()
        self.assertFalse(self.target.exists())
        self.assert_originals_unchanged()

    def test_hardlinked_partial_cannot_overwrite_original_input(self):
        self.interrupt()
        stage, = self.target.parent.glob(".GRCh38.adoption-*")
        partial, = stage.glob("*.copying")
        partial.unlink()
        original = self.source / partial.name.removesuffix(".copying")
        os.link(original, partial)  # Deliberately hostile disposable fixture.
        before, content = cache._identity(original), original.read_bytes()
        with self.assertRaisesRegex(ValueError, "hardlinked adoption partial"):
            self.adopt()
        self.assertEqual(cache._identity(original), before)
        self.assertEqual(original.read_bytes(), content)
        self.assertFalse(self.target.exists())

    def test_a_relocated_cache_can_be_imported_again_without_changing_either_source(self):
        self.adopt()
        next_root = self.folder / "second_managed"
        next_target = next_root / "SpliceAI2/GRCh38"
        before = (self.target / "adapter.ready.json").read_bytes()
        cache.adopt_completed_projection(source_dir=self.target, cache_dir=next_target, cache_root=next_root)
        self.assertEqual((self.target / "adapter.ready.json").read_bytes(), before)
        with cache._completed_projection(*self.paths(next_target)):
            pass
        self.assert_originals_unchanged()

    def test_changed_original_checksum_fails_before_publication(self):
        self.data["outputs"]["snv"]["sha256"] = "d" * 64
        self.save_receipt()
        with self.assertRaisesRegex(ValueError, "Original prepared member SHA256"):
            self.adopt()
        self.assertFalse(self.target.exists())

    def test_changed_copy_is_not_repaired_or_reused_silently(self):
        self.adopt()
        path = Path(str(self.paths(self.target)[0]) + ".tbi")
        data = bytearray(path.read_bytes()); data[-1] ^= 1
        path.write_bytes(data)
        with self.assertRaisesRegex(ValueError, "staging bytes differ"):
            self.adopt()
        self.assertEqual(path.read_bytes(), bytes(data))
        self.assert_originals_unchanged()

    def test_missing_original_receipt_copy_rejects_relocated_registration(self):
        self.adopt()
        (self.target / "adapter.original.ready.json").unlink()
        with self.assertRaises((ValueError, OSError)):
            with cache._completed_projection(*self.paths(self.target)):
                pass

    def test_scientific_metadata_cannot_be_changed_by_a_relocation(self):
        self.adopt()
        receipt = self.target / "adapter.ready.json"
        value = json.loads(receipt.read_text()); value["thresholds"]["balanced"] = 0.125
        receipt.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, "more than cache locations"):
            with cache._completed_projection(*self.paths(self.target)):
                pass

    def test_subset_receipt_keeps_existing_rejection_gate(self):
        self.data["coverage"] = "validation_subset"; self.save_receipt()
        with self.assertRaisesRegex(ValueError, "completed full GRCh38"):
            self.adopt()
        self.assertFalse(self.target.exists())

    def test_external_symlinked_or_occupied_destination_is_never_replaced(self):
        with self.assertRaisesRegex(ValueError, "inside cache_root"):
            cache.adopt_completed_projection(source_dir=self.source, cache_dir=self.folder / "external", cache_root=self.root)
        self.root.mkdir()
        (self.root / "VEP_plugin_cache").symlink_to(self.source, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "destination symlinks"):
            self.adopt()
        (self.root / "VEP_plugin_cache").unlink()
        self.target.mkdir(parents=True)
        (self.target / "retain").write_bytes(b"unknown occupied cache")
        with self.assertRaises((OSError, ValueError)):
            self.adopt()
        self.assertEqual((self.target / "retain").read_bytes(), b"unknown occupied cache")

    def test_named_cli_exposes_physical_adoption(self):
        result = subprocess.run([sys.executable, str(ROOT / "scripts/spliceai_cache_input.py"), "adopt-prepared",
            "--source-dir", str(self.source), "--output-dir", str(self.target), "--cache-root", str(self.root)],
            text=True, capture_output=True, timeout=20)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)["status"], "complete")
        self.assert_originals_unchanged()


if __name__ == "__main__":
    unittest.main()
