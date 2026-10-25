"""Fixture-only receipt checks; no downloads, annotations or real cache edits."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("local_plugin_test", ROOT / "scripts/vep_plugins_release.py")
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)


class LocalPatchTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.plugins, self.stage = self.root / "plugins", self.root / "stage"
        self.plugins.mkdir()
        self.stage.mkdir()
        (self.stage / "plugin_config.txt").write_text("fixture config\n")
        (self.stage / "UTRAnnotator.pm").write_text("package UTRAnnotator;\n# upstream fixture\n1;\n")
        helper.install(self.plugins, "116", ["UTRAnnotator"], self.stage)
        self.original = self.root / "original.pm"
        self.original.write_bytes((self.plugins / "UTRAnnotator.pm").read_bytes())
        self.patch = self.root / "reporting.patch"
        self.patch.write_text("fixture patch provenance\n")
        (self.plugins / "UTRAnnotator.pm").write_text("package UTRAnnotator;\n# local fixture\n1;\n")
        self.receipt = helper.load_receipt(self.plugins)
        self.receipt["local_patches"] = {"UTRAnnotator.pm": {
            "profile": "fixture-all-effects-v1",
            "upstream_sha256": helper.digest(self.original),
            "installed_sha256": helper.digest(self.plugins / "UTRAnnotator.pm"),
            "original_file": str(self.original), "patch_file": str(self.patch),
            "patch_sha256": helper.digest(self.patch),
        }}
        helper.write_receipt(self.plugins / helper.RECEIPT, self.receipt)

    def test_recorded_patch_passes_and_same_upstream_install_preserves_it(self):
        self.assertTrue(helper.matches(self.plugins, self.receipt, "116", ["UTRAnnotator"]))
        before = (self.plugins / "UTRAnnotator.pm").stat().st_ino
        helper.install(self.plugins, "116", ["UTRAnnotator"], self.stage)
        self.assertEqual((self.plugins / "UTRAnnotator.pm").stat().st_ino, before)
        self.assertEqual(helper.load_receipt(self.plugins), self.receipt)

    def test_new_release_and_changed_same_release_source_preserve_patch(self):
        for release, changed in (("117", False), ("116", True)):
            with self.subTest(release=release):
                if changed:
                    (self.stage / "UTRAnnotator.pm").write_text("package UTRAnnotator;\n# changed upstream\n1;\n")
                before = (self.plugins / "UTRAnnotator.pm").read_bytes()
                with self.assertRaisesRegex(ValueError, "needs review"):
                    helper.install(self.plugins, release, ["UTRAnnotator"], self.stage)
                self.assertEqual((self.plugins / "UTRAnnotator.pm").read_bytes(), before)

    def test_tampered_module_original_or_patch_invalidates_receipt(self):
        for path in (self.plugins / "UTRAnnotator.pm", self.original, self.patch):
            with self.subTest(path=path.name):
                before = path.read_bytes()
                path.write_bytes(before + b"changed\n")
                self.assertFalse(helper.matches(self.plugins, self.receipt, "116", ["UTRAnnotator"]))
                with self.assertRaisesRegex(ValueError, "needs review"):
                    helper.install(self.plugins, "116", ["UTRAnnotator"], self.stage)
                path.write_bytes(before)

    def test_unrecorded_customization_still_fails(self):
        receipt = dict(self.receipt)
        receipt.pop("local_patches")
        helper.write_receipt(self.plugins / helper.RECEIPT, receipt)
        self.assertFalse(helper.matches(self.plugins, receipt, "116", ["UTRAnnotator"]))
        with self.assertRaisesRegex(ValueError, "differs from official"):
            helper.install(self.plugins, "116", ["UTRAnnotator"], self.stage)


if __name__ == "__main__":
    unittest.main()
