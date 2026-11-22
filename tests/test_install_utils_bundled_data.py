"""Exact canonical bundles reuse opaque nonzero bytes; external builds retain their routes."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import fake_priva_binary


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("bundled_portable_paths", ROOT / "scripts/installer_config_paths.py")
PATHS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PATHS)


class BundledDataPolicy(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix="priva_bundled_policy_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temp.cleanup)
        self.folder = Path(temp.name)
        self.source = self.folder / "checkout"
        self.cache = self.folder / "managed_cache"
        self.config_path = self.folder / "config.json"
        self.config = {"assembly": "hg19", "vep_installed_version": "116",
                       "cache_root": str(self.cache), "tmp_dir": str(self.folder),
                       "gene_mechanism_cache_dir": str(self.cache / "patho_mechanism"),
                       "gene_mechanism_raw_dir": str(self.cache / "patho_mechanism/raw"),
                       "gofcards_workdir": str(self.cache / "gofcards/build_work"),
                       "clinvar_vcv_xml": str(self.cache / "patho_mechanism/raw/clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz")}
        self.bundles = {}
        for key, relative in PATHS.BUNDLED_CACHES.items():
            path = self.source / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"opaque nonempty original bundle, not decoded by installer\n")
            self.bundles[key] = path
            self.config[key] = str(path)
        for assembly in ("hg19", "hg38"):
            for key, name in (("clingen_map", f"clingen_map.{assembly}.pkl.gz"),
                              ("gene_dosage_sensitivity", f"gene_dosage_sensitivity.{assembly}.tsv")):
                path = self.source / "data/clingen" / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"opaque original ClinGen bundle\n")
                self.bundles[key + "." + assembly] = path
                if assembly == "hg19":
                    self.config[key] = str(path)
        erepo = self.source / "data/clingen/clingen_erepo.tsv"
        erepo.write_bytes(b"original bundled eRepo source\n")
        self.bundles["clingen_erepo"] = erepo
        self.env = {key: value for key, value in os.environ.items() if not key.startswith("PRIVA_FORCE_")}
        self.env.update(TEST_SOURCE=str(self.source), TEST_CONFIG=str(self.config_path),
                        TEST_PYTHON=sys.executable, TEST_FOLDER=str(self.folder), PYTHONDONTWRITEBYTECODE="1",
                        **fake_priva_binary.environment(self.folder))


    def test_unset_outputs_select_bundles_and_force_selects_managed_outputs(self):
        config = {key: self.config[key] for key in ("cache_root", "tmp_dir", "assembly", "vep_installed_version")}
        force_env = {key: "0" for key in os.environ if key.startswith("PRIVA_FORCE_")}
        with mock.patch.dict(os.environ, force_env):
            updates = PATHS.portable_updates(config, self.source)
        for key in PATHS.BUNDLED_CACHES:
            self.assertEqual(updates[key], str(self.bundles[key]))
        self.assertEqual(updates["clingen_map"], str(self.bundles["clingen_map.hg19"]))
        with mock.patch.dict(os.environ, PRIVA_FORCE_ALL_CACHES="1"):
            updates = PATHS.portable_updates(config, self.source)
        for key in (*PATHS.BUNDLED_CACHES, "clingen_map", "gene_dosage_sensitivity"):
            self.assertTrue(Path(updates[key]).is_relative_to(self.cache), key)

    def test_missing_empty_or_external_link_bundle_defaults_to_managed_output(self):
        config = {key: self.config[key] for key in ("cache_root", "tmp_dir", "assembly", "vep_installed_version")}
        for state in ("missing", "empty", "external_link"):
            with self.subTest(state=state):
                selected = self.bundles["gene_nonlof_mechanism_json"]
                if selected.exists() or selected.is_symlink():
                    selected.unlink()
                if state == "empty":
                    selected.touch()
                if state == "external_link":
                    external = self.folder / "external_nonlof.json.gz"
                    external.write_bytes(b"external bytes retain their policy")
                    selected.symlink_to(external)
                updates = PATHS.portable_updates(config, self.source)
                self.assertEqual(updates["gene_nonlof_mechanism_json"], str(self.cache / "patho_mechanism/gene_nonlof_mechanism_curated_assertions.json.gz"))


if __name__ == "__main__":
    unittest.main()
