"""Driver identity gates; the Perl boundary test executes the official plugin."""
import contextlib
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("spliceai2_identity", ROOT / "scripts/spliceai2_to_vep_vcf.py")
DEPLOY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DEPLOY)
PINS = {
    "113": "500a534c403bcbcc15753718768f0a1f352cc6b93f9119e5eba25745f7205b9a",
    "116": "d08ce7d8c2e6b3229cef638beff145e465e224f0e3fcd43fd49d501ebc2740b1",
}


class PluginIdentityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="spliceai2_identity_", dir=ROOT / "tests")
        self.addCleanup(self.temp.cleanup)
        self.plugin = Path(self.temp.name) / "SpliceAI.pm"
        self.plugin.write_text("identity fixture; official plugin execution is tested separately\n")
        self.item = {"path": str(self.plugin), "bytes": self.plugin.stat().st_size,
                     "sha256": DEPLOY.digest(self.plugin)}

    def verify(self, item=None):
        with contextlib.redirect_stdout(io.StringIO()):
            return DEPLOY.verified_plugin(item or self.item, {})

    def test_only_explicit_113_and_116_pins_are_qualified(self):
        self.assertEqual(set(DEPLOY.QUALIFIED_PLUGINS), set(PINS.values()))
        self.assertEqual({info["release"] for info in DEPLOY.QUALIFIED_PLUGINS.values()}, set(PINS))

    def test_each_qualified_identity_returns_its_actual_selection(self):
        # Stand in for the known payload's digest only. The separate real
        # Perl test pins and executes release116 bytes; checksum failure below
        # uses the unmodified hash implementation.
        for release, sha256 in PINS.items():
            with self.subTest(release=release), mock.patch.object(DEPLOY, "digest", return_value=sha256):
                result = self.verify({**self.item, "sha256": sha256})
            self.assertEqual(result["sha256"], sha256)
            self.assertEqual(result["release"], release)
            self.assertEqual(result["path"], str(self.plugin))
            self.assertIn("normalization" if release == "113" else "Exact position/REF/ALT",
                          result["allele_matching"])

    def test_unknown_plugin_is_rejected_even_with_correct_checksum(self):
        with self.assertRaisesRegex(ValueError, "Unqualified official plugin revision"):
            self.verify()

    def test_claiming_a_qualified_digest_does_not_bypass_payload_verification(self):
        for sha256 in PINS.values():
            with self.subTest(sha256=sha256), self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
                self.verify({**self.item, "sha256": sha256})

    def test_claiming_a_qualified_digest_does_not_bypass_size_verification(self):
        with self.assertRaisesRegex(ValueError, "Input size mismatch"):
            self.verify({**self.item, "sha256": PINS["116"], "bytes": self.item["bytes"] + 1})


if __name__ == "__main__":
    unittest.main()
