"""Read-only receipt routing and bounded motif/PM1 installer upgrade fixtures."""
import json
import os
from pathlib import Path
import tempfile
import time
import unittest

from test_biomart_nm_to_enst import BODY, MAPPING
import test_installer_prerequisites as prerequisites
from test_installer_prerequisites import ROOT


def write_receipt(artifact, mapping, inputs=()):
    MAPPING.atomic_json(Path(str(artifact) + ".nm_to_enst.json"), {
        "mapping_identity": MAPPING.mapping_identity(mapping, inputs),
        "artifact_bytes": artifact.stat().st_size, "artifact_sha256": MAPPING.sha256(artifact),
    })


class MappingArtifactRoutes(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(
            prefix="priva_mapping_route_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.root = self.folder / "managed cache"
        self.root.mkdir()
        self.mapping = self.folder / "selected.tsv"
        self.mapping.write_bytes(BODY)
        self.old = self.folder / "old.pkl"
        self.old.write_bytes(b"old completed scientific output")
        self.anchor = self.folder / "anchor.tsv"
        self.anchor.write_bytes(b"original anchors\n")
        write_receipt(self.old, self.mapping, [self.anchor])
        self.receipt = Path(str(self.old) + ".nm_to_enst.json")

    def route(self, **overrides):
        args = dict(artifact=self.old, mapping=self.mapping, cache_root=self.root,
                    output_dir=self.root / "VEP_plugin_cache/AlphaMissense", assembly="hg38",
                    kind="motifs", inputs=[self.anchor])
        args.update(overrides)
        return MAPPING.artifact_route(**args)

    def change_mapping(self):
        self.mapping.write_bytes(BODY.replace(b"ENST000001", b"ENST000009"))

    def test_same_identity_reuses_without_creating_directories_or_changing_old_bytes(self):
        before = self.old.read_bytes(), self.receipt.read_bytes(), self.old.stat().st_mtime_ns
        self.assertEqual(self.route(), self.old)
        os.utime(self.mapping, (1, 1))
        self.assertEqual(self.route(), self.old)
        self.assertEqual((self.old.read_bytes(), self.receipt.read_bytes(), self.old.stat().st_mtime_ns), before)
        self.assertEqual(list(self.root.iterdir()), [])

    def test_changed_mapping_chooses_stable_new_managed_namespace_and_preserves_old(self):
        before = self.old.read_bytes(), self.receipt.read_bytes(), self.old.stat().st_ctime_ns
        self.change_mapping()
        target = self.route()
        self.assertNotEqual(target, self.old)
        self.assertTrue(target.is_relative_to(self.root))
        self.assertEqual(target, self.route())
        self.assertFalse(target.parent.exists())  # Routing itself is read-only.
        self.assertEqual((self.old.read_bytes(), self.receipt.read_bytes(), self.old.stat().st_ctime_ns), before)
        other = self.route(kind="pm1")
        self.assertNotEqual(target, other)
        self.assertEqual(other.name, "pm1_regions.hg38.pkl")

    def test_changed_anchor_selects_distinct_namespace_without_mtime_identity(self):
        self.change_mapping()
        target = self.route()
        os.utime(self.anchor, (1, 1))
        self.assertEqual(self.route(), target)
        self.anchor.write_bytes(b"changed selected anchors\n")
        self.assertNotEqual(self.route(), target)

    def test_completed_new_target_reuses_with_matching_receipt(self):
        self.change_mapping()
        target = self.route()
        target.parent.mkdir(parents=True)
        target.write_bytes(b"new completed scientific output")
        write_receipt(target, self.mapping, [self.anchor])
        before = target.read_bytes(), Path(str(target) + ".nm_to_enst.json").read_bytes(), target.stat().st_mtime_ns
        self.assertEqual(self.route(), target)
        self.assertEqual(self.route(artifact=target), target)
        self.assertEqual((target.read_bytes(), Path(str(target) + ".nm_to_enst.json").read_bytes(), target.stat().st_mtime_ns), before)

    def test_missing_corrupt_or_invalid_old_receipt_never_rotates(self):
        original = self.receipt.read_bytes()
        self.change_mapping()
        invalid_identity = json.loads(original)
        invalid_identity["mapping_identity"] = {"mapping_sha256": "not a complete identity"}
        for content in (None, b"invalid JSON", json.dumps(invalid_identity).encode()):
            with self.subTest(content=content):
                if content is None:
                    self.receipt.unlink()
                else:
                    self.receipt.write_bytes(content)
                with self.assertRaisesRegex(ValueError, "missing/stale NM-to-ENST receipt"):
                    self.route()
                self.assertEqual(list(self.root.iterdir()), [])
                self.receipt.write_bytes(original)

    def test_changed_old_artifact_is_not_mistaken_for_changed_input(self):
        self.change_mapping()
        self.old.write_bytes(b"corrupt or partial old artifact")
        with self.assertRaisesRegex(ValueError, "output content changed"):
            self.route()
        self.old.unlink()
        with self.assertRaisesRegex(ValueError, "output content changed"):
            self.route()

    def test_occupied_new_artifact_or_orphan_receipt_stays_error(self):
        self.change_mapping()
        target = self.route()
        target.parent.mkdir(parents=True)
        target.write_bytes(b"unknown occupied output")
        with self.assertRaisesRegex(ValueError, "missing/stale"):
            self.route()
        write_receipt(target, self.mapping, [self.anchor])
        target.write_bytes(b"partial replacement")
        with self.assertRaisesRegex(ValueError, "output content changed"):
            self.route()
        target.unlink()
        with self.assertRaisesRegex(ValueError, "missing/stale"):
            self.route()

    def test_destination_outside_root_and_link_escapes_are_rejected(self):
        self.change_mapping()
        outside = self.folder / "outside"
        outside.mkdir()
        with self.assertRaisesRegex(ValueError, "escapes cache_root"):
            self.route(output_dir=outside)
        plugin = self.root / "plugins"
        plugin.symlink_to(outside, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "Linked managed"):
            self.route(output_dir=plugin / "AlphaMissense")
        self.assertEqual(list(outside.iterdir()), [])

    def test_valid_selected_root_alias_is_canonicalized_once(self):
        self.change_mapping()
        alias = self.folder / "root alias"
        alias.symlink_to(self.root, target_is_directory=True)
        self.assertEqual(self.route(), self.route(cache_root=alias, output_dir=alias / "VEP_plugin_cache/AlphaMissense"))

    def test_linked_or_shared_new_target_and_linked_receipt_are_rejected(self):
        self.change_mapping()
        target = self.route()
        target.parent.mkdir(parents=True)
        target.symlink_to(self.old)
        with self.assertRaisesRegex(ValueError, "Linked managed"):
            self.route()
        target.unlink()
        os.link(self.old, target)
        with self.assertRaisesRegex(ValueError, "Shared managed"):
            self.route()
        target.unlink()
        Path(str(target) + ".nm_to_enst.json").symlink_to(self.receipt)
        with self.assertRaisesRegex(ValueError, "Linked managed"):
            self.route()


class MappingArtifactInstaller(unittest.TestCase):
    def setUp(self):
        self.fixture = prerequisites.InstallerPrerequisites("test_mapping_is_passed_verbatim_and_pm1_options_are_preserved")
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.folder / "managed cache"
        self.root.mkdir()
        self.fixture.config.update(cache_root=str(self.root), vep_plugins_cachedir=str(self.root / "plugins"))
        self.anchors = [ROOT / "data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.hg38.tsv",
                        ROOT / "data/HCSeeker/HC_spots.hg38.tsv"]

    def test_changed_mapping_build_and_completed_target_fast_reuse_register_correct_key(self):
        for function, key, name in (
            ("AlphaMissense_pick_missense_constrained_segments", "alphamissense_intolerant_motifs", "motifs.pkl"),
            ("prepare_pm1_regions", "pm1_regions_pkl", "pm1.pkl"),
        ):
            with self.subTest(function=function):
                old = self.fixture.folder / name
                old.write_bytes(b"intact old profile output")
                self.fixture.config[key] = str(old)
                write_receipt(old, self.fixture.mapping, self.anchors)
                receipt = Path(str(old) + ".nm_to_enst.json")
                before = old.read_bytes(), receipt.read_bytes(), old.stat().st_ctime_ns
                self.fixture.mapping.write_bytes(self.fixture.mapping.read_bytes() +
                    b"ENSG000009\tENST000009.4\tGENE\tP00009\tNM_000009.2\n")
                self.fixture.config["vep_installed_version"] = "116"
                setup = r'''
                update_yaml() {
                    "$TEST_PYTHON" -c 'import json,sys; p,k,v=sys.argv[1:]; d=json.load(open(p)); d[k]=v; open(p,"w").write(json.dumps(d))' "$@"
                    printf 'REGISTER:%s\n' "$2"
                }
                python() {
                    if [[ $1 == */installer_prerequisites.py || $1 == */biomart_nm_to_enst.py ]]; then
                        "$TEST_PYTHON" "$@"
                    else
                        "$TEST_PYTHON" -c 'import json,sys; from pathlib import Path; a=sys.argv[1:]; print("SCI:"+json.dumps(a)); Path(a[a.index("--output")+1]).write_bytes(b"new completed profile output")' "$@"
                    fi
                }
                '''
                result = self.fixture.run_shell(setup + f'\n{function} "$TEST_CONFIG"')
                self.assertEqual(result.returncode, 0, result.stderr)
                call = next(json.loads(line[4:]) for line in result.stdout.splitlines() if line.startswith("SCI:"))
                self.assertEqual(call[call.index("--nm_to_enst_tsv") + 1], str(self.fixture.mapping))
                if key == "pm1_regions_pkl":
                    self.assertIn("--skip_das", call)
                current = json.loads(self.fixture.config_file.read_text())
                target = Path(current[key])
                self.assertTrue(target.is_relative_to(self.root))
                self.assertNotEqual(target, old)
                MAPPING.check_receipt(target, self.fixture.mapping, self.anchors)
                self.assertEqual((old.read_bytes(), receipt.read_bytes(), old.stat().st_ctime_ns), before)
                os.utime(target, (time.time() + 100, time.time() + 100))
                complete_before = target.read_bytes(), Path(str(target) + ".nm_to_enst.json").read_bytes()
                # Starting from the old config again must find and register the
                # completed new namespace, without calling a scientific builder.
                setup = setup[:setup.index("                python() {")] + r'''
                python() {
                    if [[ $1 == */installer_prerequisites.py || $1 == */biomart_nm_to_enst.py ]]; then "$TEST_PYTHON" "$@"
                    else printf 'UNEXPECTED_SCIENTIFIC_COMMAND\n'; return 99; fi
                }
                '''
                result = self.fixture.run_shell(setup + f'\n{function} "$TEST_CONFIG"')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertNotIn("UNEXPECTED", result.stdout)
                self.assertEqual(json.loads(self.fixture.config_file.read_text())[key], str(target))
                self.assertIn("REGISTER:" + key, result.stdout)
                self.assertEqual((target.read_bytes(), Path(str(target) + ".nm_to_enst.json").read_bytes()), complete_before)
                self.assertEqual((old.read_bytes(), receipt.read_bytes(), old.stat().st_ctime_ns), before)


if __name__ == "__main__":
    unittest.main()
