"""Portable deployment of the approved reporting-only UTR116 implementation.

These are software/transaction tests, not full-cache or family acceptance.
The small original module fixture retains Ensembl's Apache-2.0 license header.
"""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("utr_install_test", ROOT / "scripts/vep_plugins_release.py")
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
ORIGINAL = ROOT / "tests/fixtures/UTRAnnotator116.original.pm"


class UtrReportingInstallTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.directory, self.stage = self.root / "plugins", self.root / "stage"
        self.directory.mkdir()
        self.stage.mkdir()
        (self.stage / "plugin_config.txt").write_text("fixture configuration\n")
        (self.stage / "UTRAnnotator.pm").write_bytes(ORIGINAL.read_bytes())
        (self.stage / "SpliceAI.pm").write_text("package SpliceAI;\n1;\n")
        helper.install(self.directory, "116", ["UTRAnnotator", "SpliceAI"], self.stage)
        self.target = self.directory / "UTRAnnotator.pm"
        self.receipt = self.directory / helper.RECEIPT
        self.pending = self.directory / ".priva-utr-reporting.pending.json"

    def ensure(self):
        return helper.ensure_utr_reporting(self.directory, "116")

    def test_installs_exact_approved_bytes_preserving_other_plugins_and_receipt_fields(self):
        previous = helper.load_receipt(self.directory)
        other = self.directory / "SpliceAI.pm"
        other_before = other.read_bytes(), other.stat().st_ino
        self.target.chmod(0o640)
        self.assertEqual(self.ensure()["status"], "installed")
        current = helper.load_receipt(self.directory)
        metadata = current.pop("local_patches")["UTRAnnotator.pm"]
        self.assertEqual(current, previous)
        self.assertEqual(other_before, (other.read_bytes(), other.stat().st_ino))
        self.assertEqual(helper.digest(self.target), helper.UTR_INSTALLED_SHA)
        self.assertEqual(self.target.stat().st_mode & 0o777, 0o640)
        self.assertEqual(Path(metadata["original_file"]).read_bytes(), ORIGINAL.read_bytes())
        self.assertEqual(helper.digest(Path(metadata["patch_file"])), metadata["patch_sha256"])
        self.assertTrue(helper.matches(self.directory, helper.load_receipt(self.directory), "116", ["UTRAnnotator", "SpliceAI"]))
        self.assertFalse(self.pending.exists())

    def test_repeated_install_is_a_live_file_identity_noop(self):
        self.ensure()
        before = [(p.read_bytes(), p.stat().st_ino, p.stat().st_mtime_ns) for p in (self.target, self.receipt)]
        self.assertEqual(self.ensure()["status"], "reused")
        self.assertEqual(before, [(p.read_bytes(), p.stat().st_ino, p.stat().st_mtime_ns) for p in (self.target, self.receipt)])

    def test_existing_approved_external_backup_paths_and_receipt_bytes_remain(self):
        self.ensure()
        receipt = helper.load_receipt(self.directory)
        patch = receipt["local_patches"]["UTRAnnotator.pm"]
        for field in ("original_file", "patch_file"):
            old = Path(patch[field])
            new = self.root / old.name
            old.rename(new)
            patch[field] = str(new)
        helper.write_receipt(self.receipt, receipt)
        (self.directory / ".priva-utr-reporting.lock").unlink()
        files_before = sorted(self.directory.iterdir())
        before = self.receipt.read_bytes()
        self.assertEqual(self.ensure()["status"], "reused")
        self.assertEqual(self.receipt.read_bytes(), before)
        self.assertEqual(sorted(self.directory.iterdir()), files_before)

    def test_other_release_does_not_change_or_create_files(self):
        before = sorted(self.directory.iterdir())
        self.assertEqual(helper.ensure_utr_reporting(self.directory, "117")["status"], "not_applicable")
        self.assertEqual(sorted(self.directory.iterdir()), before)

    def test_unknown_source_is_preserved(self):
        self.target.write_bytes(self.target.read_bytes() + b"# customization\n")
        before = self.target.read_bytes(), self.receipt.read_bytes()
        with self.assertRaisesRegex(ValueError, "pinned original"):
            self.ensure()
        self.assertEqual(before, (self.target.read_bytes(), self.receipt.read_bytes()))

    def test_unrecorded_known_output_is_not_silently_adopted(self):
        self.target.write_bytes(helper.UTR_REPORTING_MODULE.read_bytes())
        with self.assertRaisesRegex(ValueError, "Unrecorded"):
            self.ensure()
        self.assertFalse(self.pending.exists())

    def test_changed_bundle_fails_before_live_publication(self):
        source = self.root / "changed.pm"
        source.write_bytes(helper.UTR_REPORTING_MODULE.read_bytes() + b"# changed\n")
        before = self.target.read_bytes(), self.receipt.read_bytes()
        with mock.patch.object(helper, "UTR_REPORTING_MODULE", source):
            with self.assertRaisesRegex(ValueError, "approved profile"):
                self.ensure()
        self.assertEqual(before, (self.target.read_bytes(), self.receipt.read_bytes()))

    def test_receipt_failure_is_recoverable_without_replacing_the_original_backup(self):
        write = helper.write_receipt
        def fail_live(path, value):
            if path == self.receipt:
                raise OSError("simulated receipt publication failure")
            return write(path, value)
        with mock.patch.object(helper, "write_receipt", side_effect=fail_live):
            with self.assertRaisesRegex(OSError, "publication failure"):
                self.ensure()
        self.assertEqual(helper.digest(self.target), helper.UTR_INSTALLED_SHA)
        self.assertTrue(self.pending.exists())
        metadata = json.loads(self.pending.read_text())["receipt_after"]["local_patches"]["UTRAnnotator.pm"]
        original = Path(metadata["original_file"])
        original_inode = original.stat().st_ino
        self.ensure()
        self.assertEqual(original.stat().st_ino, original_inode)
        self.assertEqual(helper.digest(original), helper.UTR_UPSTREAM_SHA)
        self.assertFalse(self.pending.exists())

    def test_module_failure_is_recoverable_before_receipt_publication(self):
        replace = helper.os.replace
        def fail_module(source, destination):
            if destination == self.target:
                raise OSError("simulated module publication failure")
            return replace(source, destination)
        with mock.patch.object(helper.os, "replace", side_effect=fail_module):
            with self.assertRaisesRegex(OSError, "module publication failure"):
                self.ensure()
        self.assertEqual(helper.digest(self.target), helper.UTR_UPSTREAM_SHA)
        self.assertTrue(self.pending.exists())
        self.ensure()
        self.assertEqual(helper.digest(self.target), helper.UTR_INSTALLED_SHA)

    def test_conflicting_recovery_and_backup_are_preserved(self):
        self.pending.write_text(json.dumps({"receipt_before": "unrelated", "receipt_after": {}}))
        before = self.pending.read_bytes(), self.target.read_bytes(), self.receipt.read_bytes()
        with self.assertRaisesRegex(ValueError, "context changed"):
            self.ensure()
        self.assertEqual(before, (self.pending.read_bytes(), self.target.read_bytes(), self.receipt.read_bytes()))
        self.pending.unlink()
        original = self.directory / ".priva-local-patches/UTRAnnotator116_all_effects/UTRAnnotator.original.pm"
        original.parent.mkdir(parents=True)
        original.write_bytes(b"unrelated original backup")
        with self.assertRaisesRegex(ValueError, "Conflicting.*backup"):
            self.ensure()
        self.assertEqual(original.read_bytes(), b"unrelated original backup")

    def test_tampered_recorded_backups_and_symlink_targets_fail(self):
        self.ensure()
        receipt = helper.load_receipt(self.directory)
        original = Path(receipt["local_patches"]["UTRAnnotator.pm"]["original_file"])
        original.write_bytes(original.read_bytes() + b"changed")
        with self.assertRaisesRegex(ValueError, "approved profile"):
            self.ensure()
        self.target.unlink()
        self.target.symlink_to(helper.UTR_REPORTING_MODULE)
        with self.assertRaisesRegex(ValueError, "symlinked"):
            self.ensure()

    def test_malformed_receipts_fail_before_changing_live_files(self):
        initial = helper.load_receipt(self.directory)
        for patches in ([], {"UTRAnnotator.pm": "not a mapping"}):
            with self.subTest(patches=patches):
                receipt = dict(initial, local_patches=patches)
                helper.write_receipt(self.receipt, receipt)
                before = self.target.read_bytes(), self.receipt.read_bytes()
                with self.assertRaisesRegex(ValueError, "Malformed"):
                    self.ensure()
                self.assertEqual(before, (self.target.read_bytes(), self.receipt.read_bytes()))

    def test_malformed_recovery_is_never_overwritten(self):
        for value in (None, [], {"receipt_before": {}}):
            with self.subTest(value=value):
                self.pending.write_text(json.dumps(value))
                before = self.pending.read_bytes(), self.target.read_bytes(), self.receipt.read_bytes()
                with self.assertRaisesRegex(ValueError, "Malformed.*recovery"):
                    self.ensure()
                self.assertEqual(before, (self.pending.read_bytes(), self.target.read_bytes(), self.receipt.read_bytes()))

    def test_only_reporting_region_differs_from_the_pinned_original(self):
        old, new = ORIGINAL.read_text(), helper.UTR_REPORTING_MODULE.read_text()
        self.assertEqual(hashlib.sha256(ORIGINAL.read_bytes()).hexdigest(), helper.UTR_UPSTREAM_SHA)
        self.assertEqual(old[old.index("sub count_number_ATG"):], new[new.index("sub count_number_ATG"):])
        self.assertEqual(old.split("sub run {", 1)[1].split("  my $output_five_prime_flag;", 1)[0],
                         new.split("sub run {", 1)[1].split("  my %utr_effect = %{_format_all_effects(", 1)[0])

    def test_completed_transaction_cleanup_requires_exact_owned_context(self):
        before = helper.load_receipt(self.directory)
        self.ensure()
        after = helper.load_receipt(self.directory)
        invalid = ({"receipt_after": after},
                   {"receipt_before": before, "receipt_after": after, "unrelated": True},
                   {"receipt_before": after, "receipt_after": after},
                   {"receipt_before": dict(before, unrelated=True), "receipt_after": after})
        for value in invalid:
            with self.subTest(value=value):
                self.pending.write_text(json.dumps(value))
                saved = self.pending.read_bytes(), self.target.read_bytes(), self.receipt.read_bytes()
                with self.assertRaisesRegex(ValueError, "Conflicting.*recovery"):
                    self.ensure()
                self.assertEqual(saved, (self.pending.read_bytes(), self.target.read_bytes(), self.receipt.read_bytes()))
        self.pending.write_text(json.dumps({"receipt_before": before, "receipt_after": after}))
        self.assertEqual(self.ensure()["status"], "reused")
        self.assertFalse(self.pending.exists())


if __name__ == "__main__":
    unittest.main()
