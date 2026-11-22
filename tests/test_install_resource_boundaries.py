"""Selected reverse-chain identity and the documented repeat filename boundary."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import fake_priva_binary


ROOT = Path(__file__).resolve().parents[1]


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / filename)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


PATHS = module("resource_boundary_paths", "installer_config_paths.py")
PREREQUISITES = module("resource_boundary_prerequisites", "installer_prerequisites.py")


class InstallResourceBoundaries(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="resource_boundaries_", dir=os.environ.get("PRIVA_TEST_TMPDIR"))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.source_root = self.folder / "source"
        self.bundle = self.source_root / "data/liftover/hg38ToHg19.over.chain.gz"
        self.bundle.parent.mkdir(parents=True)
        self.bundle.write_bytes(b"opaque genuine bundle fixture")
        self.cache = self.folder / "managed"
        self.config = self.folder / "config.json"
        self.events = self.folder / "events"
        self.events.touch()
        self.reference = self.folder / "reference.fa"
        self.reference.write_bytes(b">chr1\nA\n")
        self.settings = dict(cache_root=str(self.cache), tmp_dir=str(self.folder / "tmp"),
                             assembly="hg19", vep_installed_version="116", ref_genome=str(self.reference),
                             gnomad_vcf_dir=str(self.cache / "gnomAD"))
        self.env = dict(os.environ, TEST_PYTHON=sys.executable, TEST_CONFIG=str(self.config),
                        TEST_SOURCE_ROOT=str(self.source_root), TEST_EVENTS=str(self.events),
                        TEST_CACHE=str(self.cache / "VEP_plugin_cache"),
                        **fake_priva_binary.environment(self.folder))

    def write(self):
        self.config.write_text(json.dumps(self.settings))

    def resolve(self):
        PATHS.resolve_config(argparse.Namespace(config=self.config, source_root=self.source_root,
                                               managed_updates_json=None, check=False, output=None))


    def supplied_chain(self, name="hg38ToHg19.over.chain.gz"):
        source = self.folder / "manual" / name
        source.parent.mkdir(exist_ok=True)
        source.write_bytes(b"exact selected external chain; never substitute bundle")
        return source

    def test_absolute_chain_is_adopted_with_source_provenance_and_independent_bytes(self):
        source = self.supplied_chain()
        companion = Path(str(source) + ".provenance.json")
        companion.write_bytes(b'{"origin":"selected original"}\n')
        self.settings["hg38_hg19_chain"] = str(source)
        self.write()
        before = PATHS.cache_identity(source), source.read_bytes(), self.bundle.read_bytes()
        self.resolve()
        selected = Path(json.loads(self.config.read_text())["hg38_hg19_chain"])
        self.assertEqual(selected, self.cache / "liftOver/supplied/hg38_hg19_chain" / source.name)
        self.assertEqual(selected.read_bytes(), source.read_bytes())
        self.assertNotEqual(selected.stat().st_ino, source.stat().st_ino)
        self.assertEqual(Path(str(selected) + ".provenance.json").read_bytes(), companion.read_bytes())
        self.assertEqual((PATHS.cache_identity(source), source.read_bytes(), self.bundle.read_bytes()), before)

    def test_legacy_canonical_chain_basename_reuses_opaque_bundle_without_hashing(self):
        self.settings.pop("ref_genome")
        self.settings["hg38_hg19_chain"] = self.bundle.name
        self.write()
        before = PATHS.cache_identity(self.bundle)
        with mock.patch.object(PATHS, "sha256_file", side_effect=AssertionError("bundle content must not be hashed")):
            self.resolve()
        self.assertEqual(json.loads(self.config.read_text())["hg38_hg19_chain"], str(self.bundle))
        self.assertEqual(PATHS.cache_identity(self.bundle), before)
        self.assertFalse((self.cache / "installer-source-adoptions.json").exists())

    def test_existing_custom_bundled_directory_basename_keeps_its_source_identity(self):
        source = self.bundle.with_name("custom.reverse.chain.gz")
        source.write_bytes(b"selected custom basename")
        self.settings["hg38_hg19_chain"] = source.name
        self.write()
        self.resolve()
        selected = Path(json.loads(self.config.read_text())["hg38_hg19_chain"])
        self.assertTrue(selected.is_relative_to(self.cache))
        self.assertEqual(selected.read_bytes(), source.read_bytes())
        self.assertEqual(source.read_bytes(), b"selected custom basename")

    def test_missing_custom_chain_including_under_cache_root_fails_before_mutation(self):
        for source in (self.folder / "missing.chain.gz", self.cache / "missing.chain.gz"):
            with self.subTest(source=source):
                self.settings["hg38_hg19_chain"] = str(source)
                self.write()
                before = self.config.read_bytes()
                with self.assertRaisesRegex(ValueError, "hg38_hg19_chain source is missing"):
                    self.resolve()
                self.assertEqual(self.config.read_bytes(), before)
                self.assertFalse(self.cache.exists())

    def test_relative_chain_path_syntax_is_rejected_instead_of_rebased(self):
        self.settings["hg38_hg19_chain"] = "../manual/custom.chain.gz"
        self.write()
        before = self.config.read_bytes()
        with self.assertRaisesRegex(ValueError, "bundled filename or an absolute supplied path"):
            self.resolve()
        self.assertEqual(self.config.read_bytes(), before)
        self.assertFalse(self.cache.exists())


    def test_repeat_path_syntax_is_rejected_before_config_or_cache_mutation(self):
        for value in (str(self.folder / "selected.bed.gz"), "../selected.bed.gz", "subdir/selected.bed.gz", "subdir\\selected.bed.gz", ".", ".."):
            with self.subTest(value=value):
                self.settings["repeat_region_file_name"] = value
                self.write()
                before = self.config.read_bytes()
                with self.assertRaisesRegex(ValueError, "repeat_region_file_name must be one filename"):
                    self.resolve()
                with self.assertRaisesRegex(ValueError, "repeat_region_file_name must be one filename"):
                    PREREQUISITES.validate_bundled_snapshots(self.settings, self.source_root)
                self.assertEqual(self.config.read_bytes(), before)
                self.assertFalse(self.cache.exists())

    def test_custom_repeat_basename_remains_supported_without_content_reads(self):
        self.settings["repeat_region_file_name"] = "custom.reviewed.repeat.bed.gz"
        paths = ["HCSeeker/HC_spots.hg19.tsv", "gnomAD_RMC/gnomad_v2.1.1_rmc_flat.hg19.tsv",
                 "alphamissense/alphamissense_mean_score.tsv", "vep/csq_interpretations.tsv",
                 "repeats/custom.reviewed.repeat.bed.gz", "hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv",
                 "dispensable_genes/dispensable_gene_list.txt", "gofcards/gofcards_mechanism_reviews.tsv",
                 "mondo/disease_scope_overrides.tsv", "patho_mechanism/gene_nonlof_mechanism_curated_assertions.schema.json"]
        for relative in paths:
            path = self.source_root / "data" / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"opaque nonempty bundled fixture")
        self.write()
        self.resolve()
        self.assertEqual(json.loads(self.config.read_text())["repeat_region_file_name"], self.settings["repeat_region_file_name"])
        with mock.patch.object(Path, "open", side_effect=AssertionError("bundle contents must not be read")):
            PREREQUISITES.validate_bundled_snapshots(self.settings, self.source_root)


if __name__ == "__main__":
    unittest.main()
