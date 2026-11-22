"""Installer startup tests with isolated fake managers; no packages/downloads.

Run with Python's standard library: python -m unittest discover -s tests
-p test_install_utils_bootstrap.py. Also discoverable by pytest.
"""

import json
import contextlib
import importlib.util
import io
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
from unittest import mock

import fake_priva_binary


ROOT = Path(__file__).resolve().parents[1]


class PluginReleaseReceiptTests(unittest.TestCase):
    """Small text fixtures exercise publication; no upstream files or downloads."""

    @classmethod
    def setUpClass(cls):
        import importlib.util
        spec = importlib.util.spec_from_file_location("vep_plugins_release", ROOT / "scripts/vep_plugins_release.py")
        cls.helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.helper)

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="priva_plugin_release_", dir=ROOT / "tests")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.plugins = self.directory / "plugins"
        self.plugins.mkdir()
        self.names = ["SpliceAI", "SpliceVault"]

    def stage(self, release):
        directory = self.directory / ("stage-" + release)
        directory.mkdir(exist_ok=True)
        (directory / "plugin_config.txt").write_text("official config fixture " + release)
        for name in self.names:
            (directory / (name + ".pm")).write_text("package " + name + ";\n# fixture " + release + "\n1;\n")
        return directory

    def test_transition_uses_the_same_lexical_cache_root_normalization(self):
        actual = self.directory / "actual/deep"
        actual.mkdir(parents=True)
        link = self.directory / "link"
        link.symlink_to(actual, target_is_directory=True)
        plugins = self.directory / "cache/CustomPlugins"
        plugins.mkdir(parents=True)
        (plugins / "SpliceAI.pm").write_text("package SpliceAI;\n1;\n")
        prefix = self.directory / "environment"
        prefix.mkdir()
        transition = prefix / ".priva-vep-release-transition.json"
        self.helper.capture_transition(transition, {
            "cache_root": str(link) + "/../cache",
            "cache_folders": {"vep_plugin": "CustomPlugins"},
        }, "113", "116")
        record = json.loads(transition.read_text())
        self.assertEqual(record["plugins_directory"], str(plugins.resolve()))
        self.assertEqual(set(record["plugins"]), {"SpliceAI.pm"})

    def test_identical_unrecorded_files_are_verified_without_replacement(self):
        stage = self.stage("113")
        for name in self.names:
            shutil.copyfile(stage / (name + ".pm"), self.plugins / (name + ".pm"))
        before = {name: (self.plugins / (name + ".pm")).stat().st_ino for name in self.names}
        self.helper.install(self.plugins, "113", self.names, stage)
        self.assertTrue(self.helper.matches(self.plugins, self.helper.load_receipt(self.plugins), "113", self.names))
        self.assertEqual(before, {name: (self.plugins / (name + ".pm")).stat().st_ino for name in self.names})

    def test_html_or_wrong_package_cannot_publish_receipt_or_files(self):
        stage = self.stage("116")
        for text in ("<html>download temporarily unavailable</html>", "package WrongName;\n1;\n"):
            with self.subTest(text=text):
                (stage / "SpliceAI.pm").write_text(text)
                with self.assertRaisesRegex(ValueError, "does not declare package SpliceAI"):
                    self.helper.install(self.plugins, "116", self.names, stage)
                self.assertEqual(list(self.plugins.iterdir()), [])

    def test_partial_release_publication_can_retry_without_claiming_success(self):
        self.helper.install(self.plugins, "113", self.names, self.stage("113"))
        stage = self.stage("116")
        original_replace = self.helper.os.replace

        def fail_second(source, destination):
            if Path(destination).name == "SpliceVault.pm":
                raise OSError("fixture write failure")
            return original_replace(source, destination)

        with mock.patch.object(self.helper.os, "replace", side_effect=fail_second):
            with self.assertRaisesRegex(OSError, "fixture write failure"):
                self.helper.install(self.plugins, "116", self.names, stage)
        receipt = self.helper.load_receipt(self.plugins)
        self.assertEqual(receipt["release"], "113")
        self.assertFalse(self.helper.matches(self.plugins, receipt, "116", self.names))
        self.helper.install(self.plugins, "116", self.names, stage)
        self.assertTrue(self.helper.matches(self.plugins, self.helper.load_receipt(self.plugins), "116", self.names))
        self.assertFalse(list(self.plugins.glob(".priva-plugin-*")))

    def test_partial_retry_still_preserves_unknown_customized_bytes(self):
        self.helper.install(self.plugins, "113", self.names, self.stage("113"))
        stage = self.stage("116")
        shutil.copyfile(stage / "SpliceAI.pm", self.plugins / "SpliceAI.pm")
        (self.plugins / "SpliceVault.pm").write_text("package SpliceVault;\n# customized\n1;\n")
        with self.assertRaisesRegex(ValueError, "differs from official"):
            self.helper.install(self.plugins, "116", self.names, stage)
        self.assertIn("customized", (self.plugins / "SpliceVault.pm").read_text())
        self.assertEqual(self.helper.load_receipt(self.plugins)["release"], "113")


class CachePathTests(unittest.TestCase):
    """Exercise cache configuration entirely in memory; no test folders/files."""

    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("priva_cache_tests", ROOT / "scripts/installer_config_paths.py")
        helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(helper)
        cls.namespace = {
            "resolve_cache_paths": lambda config: helper.resolve_cache_paths(config, ROOT),
            "render_cache_config": helper.render_config,
            "configure_cache_paths": lambda path, config, mode: helper.configure_cache_paths(path, config, mode, ROOT),
        }

    def resolve(self, config):
        return self.namespace["resolve_cache_paths"](config)

    def test_defaults_need_only_a_root(self):
        directories, updates = self.resolve({"cache_root": "/scratch/caches/"})
        self.assertEqual(directories, {
            "vep_cache_dir": "/scratch/caches/VEP_cache",
            "vep_plugins_dir": "/scratch/caches/VEP_plugin",
            "vep_plugins_cachedir": "/scratch/caches/VEP_plugin_cache",
            "gnomad_vcf_dir": "/scratch/caches/gnomAD",
            "clinvar_vcf_dir": "/scratch/caches/ClinVar",
            "cadd_parent_dir": "/scratch/caches/CADD",
        })
        self.assertEqual(updates, directories)

    def test_partial_overrides_and_null_mapping(self):
        directories, _ = self.resolve({
            "cache_root": "/scratch/caches",
            "cache_folders": {"vep_cache": "VEP_113", "cadd": "CADD_17"},
        })
        self.assertEqual(directories["vep_cache_dir"], "/scratch/caches/VEP_113")
        self.assertEqual(directories["cadd_parent_dir"], "/scratch/caches/CADD_17")
        self.assertEqual(directories["gnomad_vcf_dir"], "/scratch/caches/gnomAD")
        self.assertEqual(self.resolve({"cache_root": "/scratch/caches", "cache_folders": None}),
                         self.resolve({"cache_root": "/scratch/caches"}))

    def test_legacy_configurations_keep_their_full_paths(self):
        self.assertEqual(self.resolve({"vep_cache_dir": "/shared/vep"}), ({}, {}))

    def test_root_change_plans_managed_children_and_external_cache_destinations(self):
        config = {
            "cache_root": "/new/caches", "vep_plugins_cachedir": "/old/plugins",
            "spliceai_snv_prescore": "/old/plugins/SpliceAI/snv.vcf.gz",
            "spliceai_indel_prescore": "/shared/SpliceAI/indel.vcf.gz",
            "primateai_prescore": "/old/plugins_extra/scores.tsv.bgz",
            "input_vcf": "/old/plugins/cohort.vcf.gz",
        }
        _, updates = self.resolve(config)
        self.assertEqual(updates["spliceai_snv_prescore"],
                         "/new/caches/VEP_plugin_cache/SpliceAI/snv.vcf.gz")
        self.assertTrue(updates["spliceai_indel_prescore"].startswith("/new/caches/VEP_plugin_cache/SpliceAI/"))
        self.assertTrue(updates["primateai_prescore"].startswith("/new/caches/VEP_plugin_cache/PrimateAI/"))
        self.assertNotIn("input_vcf", updates)

    def test_folder_rename_and_nested_old_folders(self):
        _, updates = self.resolve({
            "cache_root": "/cache", "cache_folders": {"vep_plugin": "Plugins"},
            "vep_cache_dir": "/cache/VEP", "vep_plugins_dir": "/cache/VEP/plugins",
            "loftee_repo": "/cache/VEP/plugins/LoFtee_repos/loftee-hg38",
        })
        self.assertEqual(updates["loftee_repo"], "/cache/Plugins/LoFtee_repos/loftee-hg38")

    def test_invalid_roots_are_rejected(self):
        for root in ("relative", "/", "/cache root", "/cache;$USER", 12, False, []):
            with self.subTest(root=root), self.assertRaises(ValueError):
                self.resolve({"cache_root": root})

    def test_folder_overrides_require_known_names_and_single_components(self):
        for folders in ([], False, "VEP", {"typo": "VEP"}, {"vep_cache": "../VEP"},
                        {"vep_cache": "/VEP"}, {"vep_cache": "one/two"},
                        {"vep_cache": ""}, {"vep_cache": None}, {"vep_cache": 5}):
            with self.subTest(folders=folders), self.assertRaises(ValueError):
                self.resolve({"cache_root": "/cache", "cache_folders": folders})
        with self.assertRaisesRegex(ValueError, "Indent"):
            self.resolve({"cache_root": "/cache", "vep_cache": "VEP"})

    def test_render_preserves_comments_and_unrelated_values(self):
        original = '# Header\nvep_cache_dir: "/old/#vep"  # keep comment\ninput_vcf: "keep.vcf.gz"\n'
        render = self.namespace["render_cache_config"]
        updated = render(original, {"vep_cache_dir": "/cache/VEP", "cadd_parent_dir": "/cache/CADD"})
        self.assertIn('vep_cache_dir: "/cache/VEP"  # keep comment\n', updated)
        self.assertIn('input_vcf: "keep.vcf.gz"\n', updated)
        self.assertTrue(updated.startswith('# Header\n'))
        self.assertEqual(render(updated, {"vep_cache_dir": "/cache/VEP", "cadd_parent_dir": "/cache/CADD"}), updated)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            render('vep_cache_dir: a\nvep_cache_dir: b\n', {"vep_cache_dir": "/cache/VEP"})

    def test_check_mode_never_writes_or_creates_directories(self):
        output = io.StringIO()
        with mock.patch.object(Path, "read_text", side_effect=AssertionError("Unexpected read")), \
             mock.patch.object(Path, "write_text", side_effect=AssertionError("Unexpected write")), \
             mock.patch.object(Path, "mkdir", side_effect=AssertionError("Unexpected directory")), \
             contextlib.redirect_stdout(output):
            self.namespace["configure_cache_paths"]("unused", {"cache_root": "/cache"}, "--check")
        self.assertEqual(json.loads(output.getvalue())["gnomad_vcf_dir"], "/cache/gnomAD")

    def test_write_is_idempotent_and_only_creates_selected_directories(self):
        config = {"cache_root": "/cache"}
        directories, updates = self.resolve(config)
        original = self.namespace["render_cache_config"]('cache_root: /cache\n', updates)
        with mock.patch.object(Path, "read_text", return_value=original), \
             mock.patch.object(Path, "write_text") as write, \
             mock.patch.object(Path, "mkdir", autospec=True) as mkdir, \
             contextlib.redirect_stdout(io.StringIO()):
            self.namespace["configure_cache_paths"]("unused", config, "write")
        write.assert_not_called()
        self.assertEqual({str(call.args[0]) for call in mkdir.call_args_list}, set(directories.values()))


if __name__ == "__main__":
    unittest.main()
