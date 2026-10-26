"""Stage boundaries and immutable CADD installation evidence; no real downloads."""
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("cadd_install", ROOT / "scripts/cadd_runtime_install.py")
installer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(installer)


class RuntimeInstaller(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="cadd_install_",
                                                     dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.mamba = self.root / "mamba"
        selected = patch.dict(os.environ, CONDA_PREFIX=str(self.mamba / "envs/priva_acmg"),
                              MAMBA_ROOT_PREFIX=str(self.mamba))
        selected.start()
        self.addCleanup(selected.stop)
        self.manifest = ROOT / "scripts/cadd_runtime_pins.json"
        self.pins = installer.load_pins(self.manifest)
        self.args = installer.arguments(["--config-json", "{}", "--stage", "calculate",
                                         "--backend", "container",
                                         "--state-dir", str(self.root / "state"),
                                         "--runtime-dir", str(self.root / "runtime")])
        self.config = {"assembly": "hg38", "cache_root": str(self.root / "cache"),
                       "vep_installed_version": "116"}

    def seeded_runtime(self):
        folder = self.args.runtime_dir
        folder.mkdir(parents=True)
        image = folder / "cadd-upstream-0.1.1.sif"
        image.write_bytes(b"fixed fixture image\n")
        (folder / "image.origin.json").write_text(json.dumps(
            {"image_source": self.pins["image"], "image_sha256": installer.sha256(image)}))
        installer.runtime.install_cadd_entry_points(folder)
        return folder, image

    def test_cli_defaults_to_both_and_rejects_unknown_stage(self):
        self.assertEqual(installer.arguments(["--config-json", "{}"]).stage, "both")
        with patch("sys.stderr", io.StringIO()), self.assertRaises(SystemExit):
            installer.arguments(["--config-json", "{}", "--stage", "score"])

    def test_legacy_source_without_cache_root_keeps_software_beside_priva(self):
        args = installer.arguments(["--config-json", "{}"])
        source = self.root / "legacy/CADD-scripts-1.7.2"
        paths = installer.installation_paths({"cadd_base_dir": str(source)}, args)
        runtime = self.mamba / "share/priva/cadd-1.7.2"
        self.assertEqual(paths, (None, source, runtime, runtime / "receipts", runtime / "receipts/tmp",
                                runtime / "bin/cadd-python-calculate"))
        self.assertFalse(source.exists())

    def test_missing_empty_or_null_root_never_uses_cwd_without_selected_source(self):
        for value in (None, "", "   "):
            with self.subTest(value=value), self.assertRaisesRegex(RuntimeError, "Set cache_root, cadd_base_dir or --source"):
                installer.installation_paths({"cache_root": value, "cadd_base_dir": ""}, self.args)
        self.assertFalse(self.args.runtime_dir.exists())

    def test_empty_root_with_explicit_paths_keeps_each_selected_path(self):
        self.args.source = self.root / "selected-source"
        self.args.calculation_program = self.root / "selected-python-route"
        paths = installer.installation_paths({"cache_root": ""}, self.args)
        self.assertEqual(paths[0], None)
        self.assertEqual(paths[1], self.args.source)
        self.assertEqual(paths[2], self.args.runtime_dir)
        self.assertEqual(paths[3], self.args.state_dir)
        self.assertEqual(paths[5], self.args.calculation_program)

    def test_fresh_configured_root_keeps_data_in_cache_and_runtime_beside_priva(self):
        args = installer.arguments(["--config-json", "{}"])
        cache_root = self.root / "caches"
        paths = installer.installation_paths({"cache_root": str(cache_root)}, args)
        self.assertEqual(paths[0], cache_root)
        self.assertEqual(paths[1], cache_root / "CADD/CADD-scripts-1.7.2")
        runtime = self.mamba / "share/priva/cadd-1.7.2"
        self.assertEqual(paths[2], runtime)
        self.assertEqual(paths[4], runtime / "receipts/tmp")
        self.assertEqual(paths[5], runtime / "bin/cadd-python-calculate")
        self.assertFalse(cache_root.exists())

    def test_configured_tmp_dir_selects_scratch_without_creating_it(self):
        selected = self.root / "configured-tmp"
        self.args.tmp_root = None
        paths = installer.installation_paths(dict(self.config, tmp_dir=str(selected)), self.args)
        self.assertEqual(paths[4], selected)
        self.assertFalse(selected.exists())

    def test_explicit_tmp_root_overrides_configured_scratch(self):
        selected = self.root / "explicit-tmp"
        args = installer.arguments(["--config-json", "{}", "--tmp-root", str(selected)])
        paths = installer.installation_paths(dict(self.config, tmp_dir=str(self.root / "configured-tmp")), args)
        self.assertEqual(paths[4], selected)
        self.assertFalse(selected.exists())

    def test_environment_tmp_root_overrides_configured_scratch(self):
        selected = self.root / "environment-tmp"
        with patch.dict(os.environ, {"PRIVA_CADD_TMP_ROOT": str(selected)}):
            args = installer.arguments(["--config-json", "{}"])
        paths = installer.installation_paths(dict(self.config, tmp_dir=str(self.root / "configured-tmp")), args)
        self.assertEqual(paths[4], selected)
        self.assertFalse(selected.exists())

    def test_sibling_runtime_is_not_selected_outside_configured_cache_root(self):
        cache_root = self.root / "caches"
        legacy = self.root / "softwares/cadd-1.7.2"
        self.args.runtime_dir = legacy
        _, image = self.seeded_runtime()
        before = installer.cache.identity(image)
        self.args.runtime_dir = None
        paths = installer.installation_paths({"cache_root": str(cache_root)}, self.args, self.pins)
        self.assertEqual(paths[2], cache_root / "runtime/cadd-1.7.2")
        self.assertEqual(installer.cache.identity(image), before)
        self.assertFalse(cache_root.exists())

    def test_unrelated_sibling_image_is_retained_and_does_not_block_installation(self):
        cache_root = self.root / "caches"
        legacy = self.root / "softwares/cadd-1.7.2"
        legacy.mkdir(parents=True)
        image = legacy / "cadd-upstream-0.1.1.sif"
        image.write_bytes(b"manual image")
        args = installer.arguments(["--config-json", "{}"])
        paths = installer.installation_paths({"cache_root": str(cache_root)}, args, self.pins)
        self.assertEqual(paths[2], self.mamba / "share/priva/cadd-1.7.2")
        self.assertEqual(image.read_bytes(), b"manual image")
        self.assertFalse(cache_root.exists())

    def test_sibling_reference_and_external_receipt_do_not_override_cache_root(self):
        cache_root = self.root / "caches"
        legacy = self.root / "softwares/cadd-1.7.2"
        legacy.mkdir(parents=True)
        (legacy / "cadd-upstream-0.1.1.sif").write_bytes(b"manual image")
        self.args.runtime_dir = None
        self.args.state_dir.mkdir()
        (self.args.state_dir / "cadd-runtime.ready.json").write_text(json.dumps(
            {"image_source": self.pins["image"], "image_sha256": "a" * 64}))
        paths = installer.installation_paths({"cache_root": str(cache_root)}, self.args, self.pins)
        self.assertEqual(paths[2], cache_root / "runtime/cadd-1.7.2")
        self.assertEqual(paths[3], self.args.state_dir)

    def test_legacy_calculation_configuration_reaches_complete_stage_receipt(self):
        self.config.pop("cache_root")
        self.config["cadd_base_dir"] = str(self.root / "legacy-source")
        receipt, annotation, _, _, _ = self.invoke_stage("calculate")
        self.assertEqual(json.loads(receipt.read_text())["source"], self.config["cadd_base_dir"])
        annotation.assert_not_called()

    def test_manifest_cannot_silently_change_scientific_source_or_python_commands(self):
        for field in ("source_commit", "wrapper_sha256", "implementation"):
            manifest = json.loads(self.manifest.read_text())
            if field != "source_commit":
                manifest["calculation"][field] = "changed"
            else:
                manifest[field] = "0" * 40
            path = self.root / "changed.json"
            path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(RuntimeError, "pin"):
                installer.load_pins(path)

    def test_suspended_standalone_manifest_and_cli_are_rejected_without_acquisition(self):
        manifest = dict(self.pins, standalone={"source_url": "https://example.invalid/cadd-rs"})
        path = self.root / "suspended.json"
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(RuntimeError, "Rust CADD branch is suspended"):
            installer.load_pins(path)
        with patch("sys.stderr", io.StringIO()), self.assertRaises(SystemExit):
            installer.arguments(["--config-json", "{}", "--standalone-program", "old-program"])

    def test_mutable_container_tag_is_rejected(self):
        pins = dict(self.pins, image="docker://example/cadd:latest")
        path = self.root / "mutable.json"
        path.write_text(json.dumps(pins))
        with self.assertRaisesRegex(RuntimeError, "immutable"):
            installer.load_pins(path)

    def test_existing_image_is_hashed_and_never_redownloaded(self):
        folder, image = self.seeded_runtime()
        with patch.object(installer, "run", side_effect=AssertionError("no transfer")):
            record, _ = installer.image_and_routes(folder, self.pins, self.root / "state",
                                                   self.root / "tmp", verify_only=True)
        self.assertEqual(record["sha256"], installer.sha256(image))

    def test_container_wrapper_binds_selected_resources_and_scratch_without_site_roots(self):
        folder, _ = self.seeded_runtime()
        tools = self.root / "tools"
        tools.mkdir()
        fake = tools / "apptainer"
        fake.write_text('#!/bin/sh\nprintf "%s\\n" "$@"\n')
        fake.chmod(0o755)
        source = self.root / "selected-cache/CADD/source"
        source.mkdir(parents=True)
        work = self.root / "work"
        work.mkdir()
        environment = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ["PATH"],
                           CADD=str(source), CADD_SOURCE_ROOT=str(source), TMPDIR=str(work))
        completed = subprocess.run([str(folder / "bin/cadd-python2"), "-c", "print(1)"],
                                   check=True, capture_output=True, text=True, env=environment, cwd=work)
        arguments = completed.stdout.splitlines()
        self.assertIn(str(source) + ":" + str(source) + ":ro", arguments)
        self.assertIn(str(work) + ":" + str(work) + ":rw", arguments)
        self.assertNotIn("/stornext/snfs170/next-gen/scratch/yangyxt:/stornext/snfs170/next-gen/scratch/yangyxt:ro", arguments)
        self.assertNotIn("/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade:/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade:rw", arguments)

    def test_changed_image_never_publishes_success(self):
        folder, image = self.seeded_runtime()
        image.write_bytes(b"different image\n")
        with patch.object(installer, "run", side_effect=AssertionError("no replacement")), \
                self.assertRaisesRegex(RuntimeError, "differs"):
            installer.image_and_routes(folder, self.pins, self.root / "state", self.root / "tmp")
        self.assertEqual(image.read_bytes(), b"different image\n")

    def test_interrupted_image_is_preserved_without_force_overwrite(self):
        folder = self.args.runtime_dir
        folder.mkdir()
        partial = folder / "cadd-upstream-0.1.1.sif.partial"
        partial.write_bytes(b"interrupted image")
        with patch.object(installer, "run", side_effect=AssertionError("no replacement")), \
                self.assertRaisesRegex(RuntimeError, "Retained incomplete"):
            installer.image_and_routes(folder, self.pins, self.root / "state", self.root / "tmp")
        self.assertEqual(partial.read_bytes(), b"interrupted image")

    def test_verify_only_refuses_missing_image_before_acquisition(self):
        with patch.object(installer, "run", side_effect=AssertionError("no acquisition")), \
                self.assertRaisesRegex(RuntimeError, "not installed"):
            installer.image_and_routes(self.args.runtime_dir, self.pins, self.root / "state",
                                       self.root / "tmp", verify_only=True)

    def invoke_stage(self, stage, annotation_failure=False, candidate_failure=False,
                     backend='container', profile_failure=False):
        self.args.stage = stage
        self.args.backend = backend
        image = self.root / "fixture.sif"
        image.write_bytes(b"fixture image")
        image_record = {"path": str(image), "identity": installer.cache.identity(image)}
        annotation = patch.object(installer, "annotation_resources", return_value={"verified": True})
        if annotation_failure:
            annotation = patch.object(installer, "annotation_resources", side_effect=RuntimeError("unfinished archive"))
        candidate = patch.object(installer.transcripts, "install", return_value={"vep_release": "116", "cadd_enabled": False})
        if candidate_failure:
            candidate = patch.object(installer.transcripts, "install", side_effect=ValueError("unfinished transcript candidate"))
        deployment = patch.object(installer.native, "install_profile", return_value={"status": "installed_not_scientifically_qualified"})
        if profile_failure:
            deployment = patch.object(installer.native, "install_profile", side_effect=RuntimeError("unfinished native profile"))
        assets = {"annotation_assets": {"GRCh38": {"status": "readable_not_scientifically_qualified"}}}
        with patch.object(installer.shutil, "which", return_value="available"), \
                patch.object(installer, "checkout"), \
                patch.object(installer, "image_and_routes", return_value=(image_record, {})), \
                patch.object(installer, "native_environments", return_value=(image_record, {})), \
                deployment as deployment_call, \
                annotation as annotation_call, \
                candidate as candidate_call, \
                patch.object(installer.runtime, "check_cadd_environments", return_value={}) as checks, \
                patch.object(installer.runtime, "cadd_asset_readiness", return_value=assets), \
                patch.object(installer, "calculation_resources", return_value=[{"original": True}]) as models, \
                patch.object(installer, "python_calculation", return_value={"implementation": "all_four_calculation_stages_upstream_python"}) as bridge:
            receipt = installer.install(self.config, self.args, self.pins)
        self.candidate_call = candidate_call
        self.deployment_call = deployment_call
        return receipt, annotation_call, checks, models, bridge

    def test_native_both_requires_unqualified_profile_publication(self):
        receipt, *_ = self.invoke_stage("both", backend="native")
        self.deployment_call.assert_called_once_with(self.config, receipt, False)
        self.assertEqual(json.loads(receipt.read_text())["status"], "installed_not_scientifically_qualified")

    def test_native_profile_failure_propagates_after_completed_component_receipt(self):
        with self.assertRaisesRegex(RuntimeError, "unfinished native profile"):
            self.invoke_stage("both", backend="native", profile_failure=True)
        completed_components = list((self.root / "state").glob("*.ready.json"))
        self.assertEqual(len(completed_components), 1)
        self.assertFalse(json.loads(completed_components[0].read_text())["qualification"]["scientific_profile_activated"])

    def test_partial_native_components_do_not_publish_a_complete_driver_profile(self):
        for stage in ("annotate", "calculate"):
            self.invoke_stage(stage, backend="native")
            self.deployment_call.assert_not_called()

    def test_calculation_only_never_fetches_annotation_or_checks_feature_models(self):
        receipt, annotation, checks, models, bridge = self.invoke_stage("calculate")
        self.assertIsNone(json.loads(receipt.read_text())["selected_transcript_candidate"])
        self.candidate_call.assert_not_called()
        annotation.assert_not_called()
        self.assertEqual(checks.call_args.kwargs, {"names": ("cadd-python2",), "include_vep": False})
        models.assert_called_once()
        bridge.assert_called_once()
        record = json.loads(receipt.read_text())
        self.assertEqual(record["status"], "installed_not_scientifically_qualified")
        self.assertFalse(record["qualification"]["shared_vep_qualified"])
        self.assertFalse(record["qualification"]["scientific_profile_activated"])
        self.assertEqual(record["qualification"]["legacy_transcript_release"], 110)
        self.assertEqual(record["qualification"]["shared_priva_vep_release"], "116")
        self.assertEqual(record["calculation_runtime"]["implementation"], "all_four_calculation_stages_upstream_python")
        self.assertNotIn("standalone", record)

    def test_annotation_only_avoids_calculation_program(self):
        receipt, annotation, checks, models, bridge = self.invoke_stage("annotate")
        annotation.assert_called_once()
        models.assert_not_called()
        bridge.assert_not_called()
        self.assertIn("cadd-mmsplice-python", checks.call_args.kwargs["names"])
        self.assertFalse(checks.call_args.kwargs["include_vep"])
        self.candidate_call.assert_called_once_with(self.config, verify_only=False)
        self.assertEqual(json.loads(receipt.read_text())["selected_transcript_candidate"],
                         {"vep_release": "116", "cadd_enabled": False})

    def test_receipts_separate_configured110_labels_from_coordinate_resource_releases(self):
        for assembly, release in (("hg19", 75), ("hg38", 110)):
            with self.subTest(assembly=assembly):
                self.config["assembly"] = assembly
                receipt, _, _, _, _ = self.invoke_stage("calculate")
                record = json.loads(receipt.read_text())
                evidence = record["qualification"]
                self.assertEqual(evidence["legacy_transcript_release"], 110)
                self.assertEqual(evidence["legacy_coordinate_annotation_release"], release)
                self.assertEqual(evidence["legacy_transcript_release_evidence"]["config_key"], "EnsemblDB")
                self.assertTrue(evidence["legacy_transcript_release_evidence"]["esm_peptide_path"].endswith("/esm/pep.110.fa"))
                expected_assembly = "GRCh37" if assembly == "hg19" else "GRCh38"
                self.assertTrue(evidence["legacy_coordinate_annotation_evidence"]["config_path"].endswith(f"references_{expected_assembly}_v1.7.cfg"))
                self.assertFalse(record["qualification"]["scientific_profile_activated"])

    def test_both_selects_union(self):
        _, annotation, _, models, bridge = self.invoke_stage("both")
        annotation.assert_called_once()
        models.assert_called_once()
        bridge.assert_called_once()

    def test_auto_backend_is_native_and_explicit_image_keeps_container(self):
        args = installer.arguments(["--config-json", "{}"])
        folder = self.root / "fresh-runtime"
        self.assertEqual(installer.select_backend(args, folder), "native")
        folder.mkdir()
        (folder / "cadd-upstream-0.1.1.sif").write_bytes(b"existing reference")
        self.assertEqual(installer.select_backend(args, folder), "native")
        args.runtime_dir = folder
        self.assertEqual(installer.select_backend(args, folder), "container")
        args.backend = "native"
        self.assertEqual(installer.select_backend(args, folder), "native")

    def test_explicit_container_failure_does_not_attempt_native_substitution(self):
        with patch.object(installer.shutil, "which", return_value="/fixture/tool"), \
                patch.object(installer, "checkout"), \
                patch.object(installer, "image_and_routes", side_effect=RuntimeError("reference image unavailable")), \
                patch.object(installer, "native_environments", side_effect=AssertionError("no native fallback")), \
                self.assertRaisesRegex(RuntimeError, "reference image unavailable"):
            installer.install(self.config, self.args, self.pins)

    def test_native_plan_creates_once_then_reuses_completed_prefix(self):
        source = self.root / "source"
        env_file = source / "envs/environment_minimal.yml"
        env_file.parent.mkdir(parents=True)
        env_file.write_text("dependencies: [python=2.7.13]\n")
        folder = self.root / "native-runtime"
        prefix = self.mamba / "envs/cadd-python2"
        def install_prefix(command, **kwargs):
            (prefix / "bin").mkdir(parents=True)
            (prefix / "bin/python").write_text("fixture interpreter")
            (prefix / "conda-meta").mkdir()
            (prefix / "conda-meta/history").write_text("managed history")
            (prefix / "conda-meta/python.json").write_text('{"name":"python","version":"2.7.13"}')
        with patch.object(installer, "run", side_effect=install_prefix) as run:
            record, _ = installer.native_environments(folder, source, ("cadd-python2",),
                self.root / "tmp", self.pins, manager="/selected/micromamba")
            self.assertEqual(record["backend"], "native")
            self.assertEqual(run.call_count, 1)
            self.assertIn("python=2.7.13=0", run.call_args.args[0])
            self.assertIn(str(self.mamba), run.call_args.args[0])
            self.assertEqual(run.call_args.kwargs["env"]["PIP_CACHE_DIR"], str(folder / "pip-cache"))
        with patch.object(installer, "run", side_effect=AssertionError("no package solve")):
            installer.native_environments(folder, source, ("cadd-python2",), self.root / "tmp",
                                          self.pins, verify_only=True)
        (prefix / "conda-meta/python.json").write_text('{"name":"python","version":"different"}')
        with self.assertRaisesRegex(RuntimeError, "metadata changed"):
            installer.native_environments(folder, source, ("cadd-python2",), self.root / "tmp",
                                          self.pins, verify_only=True)

    def test_native_installer_preserves_unmanaged_prefix(self):
        source = self.root / "source"
        env_file = source / "envs/environment_minimal.yml"
        env_file.parent.mkdir(parents=True)
        env_file.write_text("dependencies: []\n")
        prefix = self.root / "native/envs/cadd-python2"
        prefix.mkdir(parents=True)
        marker = prefix / "user-file"
        marker.write_bytes(b"user content")
        with patch.object(installer, "run", side_effect=AssertionError("no mutation")), \
                self.assertRaisesRegex(RuntimeError, "Unmanaged existing"):
            installer.native_environments(self.root / "native", source, ("cadd-python2",),
                                          self.root / "tmp", self.pins)
        self.assertEqual(marker.read_bytes(), b"user content")

    def test_native_prefix_reuses_byte_identical_source_at_new_cache_path(self):
        source = self.root / "old-source"
        original = source / "envs/environment_minimal.yml"
        original.parent.mkdir(parents=True)
        original.write_text("dependencies: [python=2.7.13]\n")
        folder = self.root / "runtime"
        prefix = self.mamba / "envs/cadd-python2"
        def install_prefix(command, **kwargs):
            (prefix / "bin").mkdir(parents=True)
            (prefix / "bin/python").write_text("fixture interpreter")
            (prefix / "conda-meta").mkdir()
            (prefix / "conda-meta/history").write_text("managed history")
            (prefix / "conda-meta/python.json").write_text('{"name":"python","version":"2.7.13"}')
        with patch.object(installer, "run", side_effect=install_prefix):
            before, _ = installer.native_environments(folder, source, ("cadd-python2",),
                self.root / "tmp", self.pins, manager="/selected/micromamba")
        selected = self.root / "cache/CADD/CADD-scripts-1.7.2"
        adopted = selected / "envs/environment_minimal.yml"
        adopted.parent.mkdir(parents=True)
        adopted.write_bytes(original.read_bytes())
        ready = folder / "native-state/cadd-python2.ready.json"
        identity = installer.cache.identity(ready)
        with patch.object(installer, "run", side_effect=AssertionError("no package solve or changes")):
            after, _ = installer.native_environments(folder, selected, ("cadd-python2",),
                                                    self.root / "tmp", self.pins, verify_only=True)
        self.assertEqual(after, before)
        self.assertEqual(installer.cache.identity(ready), identity)
        self.assertEqual(json.loads(ready.read_text())["plan"]["environment_file"], str(original))

    def esm_source(self, source, original=None):
        folder = source / "data/annotations/GRCh37_v1.7/esm"
        folder.mkdir(parents=True)
        for name in installer.runtime.ESM_MODELS:
            target = folder / (name + ".pt")
            if original:
                os.link(original / target.name, target)
            else:
                target.write_bytes(b"original model fixture\n")
        return folder

    def test_esm_checkpoint_handoff_retargets_same_model_inode_and_reuses_links(self):
        source = self.root / "old-source"
        original = self.esm_source(source)
        folder = self.root / "runtime"
        folder.mkdir()
        installer.runtime.cadd_asset_readiness(source, folder, assemblies=("GRCh37",), stages=("annotate",))
        selected = self.root / "cache/CADD/CADD-scripts-1.7.2"
        adopted = self.esm_source(selected, original)
        with patch("sys.stdout", io.StringIO()):
            expected = installer.runtime.cadd_asset_readiness(selected, folder, assemblies=("GRCh37",), stages=("annotate",))
        links = [folder / "esm-hub/checkpoints" / (name + ".pt") for name in installer.runtime.ESM_MODELS]
        identities = [installer.cache.identity(link) for link in links]
        for link in links:
            self.assertEqual(os.readlink(link), str(adopted / link.name))
            self.assertTrue(link.samefile(original / link.name))
        actual = installer.runtime.cadd_asset_readiness(selected, folder, assemblies=("GRCh37",), stages=("annotate",))
        self.assertEqual(actual, expected)
        self.assertEqual([installer.cache.identity(link) for link in links], identities)

    def test_esm_checkpoint_handoff_rejects_different_model_objects_before_any_retarget(self):
        source = self.root / "old-source"
        original = self.esm_source(source)
        folder = self.root / "runtime"
        folder.mkdir()
        installer.runtime.cadd_asset_readiness(source, folder, assemblies=("GRCh37",), stages=("annotate",))
        selected = self.root / "cache/CADD/CADD-scripts-1.7.2"
        adopted = self.esm_source(selected, original)
        last = adopted / (installer.runtime.ESM_MODELS[-1] + ".pt")
        payload = last.read_bytes()
        last.unlink()
        last.write_bytes(payload)  # Equal bytes alone do not authorize link replacement.
        links = [folder / "esm-hub/checkpoints" / (name + ".pt") for name in installer.runtime.ESM_MODELS]
        before = [(os.readlink(link), installer.cache.identity(link)) for link in links]
        with self.assertRaisesRegex(RuntimeError, "different target"):
            installer.runtime.cadd_asset_readiness(selected, folder, assemblies=("GRCh37",), stages=("annotate",))
        self.assertEqual([(os.readlink(link), installer.cache.identity(link)) for link in links], before)

    def test_native_manager_keeps_initialized_executable_outside_path(self):
        manager = self.root / "private-tools/micromamba"
        manager.parent.mkdir()
        manager.write_text("#!/bin/sh\nexit 0\n")
        manager.chmod(0o755)
        with patch.dict(os.environ, {"MAMBA_EXE": str(manager)}), \
                patch.object(installer.shutil, "which", return_value=None):
            self.assertEqual(installer.native_package_manager(), str(manager))

    def test_conda_native_install_preserves_original_pip_ownership_and_private_caches(self):
        source = self.root / "source"
        env_file = source / "envs/esm.yml"
        env_file.parent.mkdir(parents=True)
        env_file.write_text("channels: [pytorch, nvidia, bioconda, conda-forge, nodefaults]\n"
                            "dependencies:\n  - python\n  - numpy\n  - pip\n"
                            "  - pip: [fair-esm==2.0.0]\n")
        folder = self.root / "runtime"
        prefix = self.mamba / "envs/cadd-esm-python"
        calls = []
        def install_prefix(command, **kwargs):
            calls.append(command)
            self.assertNotIn("PYTHONHOME", kwargs["env"])
            self.assertNotIn("PYTHONPATH", kwargs["env"])
            self.assertEqual(kwargs["env"]["CONDA_PKGS_DIRS"], str(self.mamba / "pkgs"))
            self.assertEqual(kwargs["env"]["CONDA_CHANNEL_PRIORITY"], "flexible")
            if len(calls) == 1:
                self.assertNotIn("--no-rc", command)
                self.assertNotIn("--root-prefix", command)
                self.assertNotIn("--channel-priority", command)
                self.assertIn("--override-channels", command)
                self.assertIn("numpy=2.1.3", command)
                self.assertNotIn("fair-esm==2.0.0", command)
                self.assertNotIn("nodefaults", command)
                (prefix / "bin").mkdir(parents=True)
                (prefix / "bin/python").write_text("fixture interpreter")
                (prefix / "conda-meta").mkdir()
                (prefix / "conda-meta/history").write_text("managed history")
            else:
                self.assertEqual(command[:4], [str(prefix / "bin/python"), "-m", "pip", "install"])
                self.assertEqual(kwargs["env"]["CONDA_PREFIX"], str(prefix))
                self.assertEqual(kwargs["env"]["PATH"].split(os.pathsep)[0], str(prefix / "bin"))
                requirements = Path(command[-1])
                self.assertTrue(requirements.is_relative_to(folder))
                self.assertEqual(requirements.read_text(), "fair-esm==2.0.0\n")
        with patch.dict(os.environ, PYTHONHOME="/foreign/python", PYTHONPATH="/foreign/modules"), \
                patch.object(installer, "run", side_effect=install_prefix):
            installer.native_environments(folder, source, ("cadd-esm-python",),
                                          self.root / "tmp", self.pins, manager="/selected/conda")
        self.assertEqual(len(calls), 2)
        ready = json.loads((folder / "native-state/cadd-esm-python.ready.json").read_text())
        self.assertEqual(ready["pip_command"], calls[1])

    def test_native_feature_yaml_stays_under_runtime_and_retains_reference_distribution_pin(self):
        source = self.root / "source"
        env_file = source / "envs/regulatorySequence.yml"
        env_file.parent.mkdir(parents=True)
        env_file.write_text("channels: [bioconda, conda-forge, nodefaults]\ndependencies: [seqiolib]\n")
        folder = self.root / "runtime"
        prefix = self.mamba / "envs/cadd-regseq-python"
        def install_prefix(command, **kwargs):
            staged = Path(command[command.index("--file") + 1])
            self.assertTrue(staged.is_relative_to(folder))
            self.assertEqual(staged.read_bytes(), env_file.read_bytes())
            # seqiolib 0.2.4 exposes a stale module __version__ of 0.2.0.
            self.assertIn("seqiolib=0.2.4", command)
            self.assertNotIn("seqiolib=0.2.0", command)
            self.assertIn("tensorflow-gpu=2.6.2", command)
            (prefix / "bin").mkdir(parents=True)
            (prefix / "bin/python").write_text("fixture interpreter")
            (prefix / "conda-meta").mkdir()
            (prefix / "conda-meta/history").write_text("managed history")
        before = env_file.stat()
        with patch.object(installer, "run", side_effect=install_prefix):
            installer.native_environments(folder, source, ("cadd-regseq-python",),
                                          self.root / "tmp", self.pins, manager="/selected/micromamba")
        self.assertEqual(env_file.stat().st_mtime_ns, before.st_mtime_ns)
        self.assertEqual(list(env_file.parent.iterdir()), [env_file])

    def external_annotation_fixture(self, typed=False):
        source = self.root / "old-cache/CADD/CADD-scripts-1.7.2"
        target = source / "data/annotations/GRCh37_v1.7"
        (target / "vep").mkdir(parents=True)
        asset = target / "vep/fixture.txt"
        asset.write_bytes(b"unchanged original annotation")
        archive = self.root / "old-cache/.downloads/GRCh37_v1.7.tar.gz"
        archive.parent.mkdir()
        archive.write_bytes(b"publisher archive fixture")
        receipt_root = self.root / "external-receipts"
        receipt_root.mkdir()
        md5 = installer.cache.checksum(archive)
        config = {"cadd_version": "v1.7", "cadd_GRCh37_anno_url": "https://example.invalid/GRCh37.tar.gz",
                  "cadd_GRCh37_anno_md5": md5}
        record = {"schema": 1, "assembly": "GRCh37", "version": "v1.7",
                  "archive_url": config["cadd_GRCh37_anno_url"], "archive_md5": md5,
                  "archive_identity": installer.cache.identity(archive), "members": {
                      "vep": {"kind": "directory"}, "vep/fixture.txt": {"kind": "file",
                          "identity": installer.cache.identity(asset), "sha256": installer.sha256(asset)}}}
        if typed:
            record["archive_identity"] = installer.cache.receipt_identity(archive)
            record["members"]["vep/fixture.txt"]["identity"] = installer.cache.receipt_identity(asset)
        receipt = receipt_root / "GRCh37.members.ready.json"
        receipt.write_text(json.dumps(record))
        return source, target, asset, archive, receipt_root, config, receipt

    def test_external_annotation_receipt_allows_runtime_install_without_source_writes(self):
        source, target, asset, archive, receipt_root, config, receipt = self.external_annotation_fixture()
        before = installer.cache.identity(asset)
        ready_before = receipt.read_bytes(), installer.cache.identity(receipt)
        with patch.object(installer.cache, "annotations", side_effect=AssertionError("no source writes/acquisition")), \
                patch.object(installer.cache, "checksum", side_effect=self.metadata_hash_only(receipt)):
            result = installer.annotation_resources(config, "GRCh37", source, self.root / "new-cache",
                                                    receipt_root=receipt_root)
        self.assertEqual(result["archive"], str(archive))
        self.assertEqual(result["archive_identity"], installer.cache.receipt_identity(archive))
        self.assertEqual(installer.cache.identity(asset), before)
        self.assertEqual((receipt.read_bytes(), installer.cache.identity(receipt)), ready_before)
        self.assertFalse((target / installer.cache.RECEIPT).exists())

    def metadata_hash_only(self, receipt):
        original = installer.cache.checksum
        def checksum(path, algorithm="md5"):
            self.assertEqual((path, algorithm), (receipt, "sha256"), "annotation payload must not be rehashed")
            return original(path, algorithm)
        return checksum

    def test_verify_only_reuses_typed_archive_and_member_identities_without_writes(self):
        source, target, asset, archive, receipt_root, config, receipt = self.external_annotation_fixture(typed=True)
        before = {path:installer.cache.identity(path) for path in (archive, asset, receipt)}
        ready_bytes = receipt.read_bytes()
        with patch.object(installer.cache, "annotations", side_effect=AssertionError("verify-only cannot publish")), \
                patch.object(installer.cache, "checksum", side_effect=self.metadata_hash_only(receipt)):
            result = installer.annotation_resources(config, "GRCh37", source, self.root / "new-cache",
                                                    verify_only=True, receipt_root=receipt_root)
        self.assertEqual(result["archive_identity"], json.loads(ready_bytes)["archive_identity"])
        self.assertEqual({path:installer.cache.identity(path) for path in before}, before)
        self.assertEqual(receipt.read_bytes(), ready_bytes)
        self.assertFalse((target/installer.cache.RECEIPT).exists())

    def test_unknown_archive_identity_requires_md5_in_read_only_runtime_verification(self):
        source, target, asset, archive, receipt_root, config, receipt = self.external_annotation_fixture(typed=True)
        value = json.loads(receipt.read_text())
        value["archive_identity"] = {"kind":"unknown-future-format"}
        receipt.write_text(json.dumps(value))
        before = receipt.read_bytes(), installer.cache.identity(receipt)
        with patch.object(installer.cache, "checksum", wraps=installer.cache.checksum) as hashed:
            installer.annotation_resources(config, "GRCh37", source, self.root/"new-cache",
                                           verify_only=True, receipt_root=receipt_root)
        self.assertEqual([value.args for value in hashed.call_args_list], [(archive,), (receipt, "sha256")])
        self.assertEqual((receipt.read_bytes(), installer.cache.identity(receipt)), before)

    def test_unfinished_annotation_never_publishes_stage_receipt(self):
        with self.assertRaisesRegex(RuntimeError, "unfinished archive"):
            self.invoke_stage("both", annotation_failure=True)
        self.assertEqual(list((self.root / "state").glob("*.ready.json")), [])

    def test_unfinished_selected_transcripts_never_publish_stage_receipt(self):
        with self.assertRaisesRegex(ValueError, "unfinished transcript candidate"):
            self.invoke_stage("both", candidate_failure=True)
        self.assertEqual(list((self.root / "state").glob("*.ready.json")), [])

    def test_existing_receipt_is_immutable_and_identical_reuse_keeps_inode(self):
        path = self.root / "ready.json"
        installer.publish_json(path, {"status": "installed", "checked_utc": "old"})
        before = installer.cache.identity(path)
        installer.publish_json(path, {"status": "installed", "checked_utc": "new"})
        self.assertEqual(installer.cache.identity(path), before)
        self.assertEqual(json.loads(path.read_text())["checked_utc"], "old")
        with self.assertRaisesRegex(RuntimeError, "differs"):
            installer.publish_json(path, {"status": "different"})
        self.assertEqual(installer.cache.identity(path), before)

    def test_shared_writer_lock_stops_duplicate_installation(self):
        lock = self.root / "shared.lock"
        with installer.locked(lock), self.assertRaisesRegex(RuntimeError, "writer is active"):
            with installer.locked(lock):
                self.fail("entered active writer lock")

    def test_python_route_is_regular_frozen_file_without_rust_clone_build_or_binary(self):
        source = self.root / "upstream-source"
        program = self.root / "installed/bin/cadd-python-calculate"
        with patch.object(installer, "checkout", side_effect=AssertionError("no second source checkout")), \
                patch.object(installer.subprocess, "run", side_effect=AssertionError("no build or binary probe")):
            record = installer.python_calculation(program, source, self.args.runtime_dir, "GRCh38")
        self.assertTrue(program.is_file())
        self.assertFalse(program.is_symlink())
        self.assertTrue(os.access(program, os.X_OK))
        self.assertEqual(record["sha256"], self.pins["calculation"]["wrapper_sha256"])
        self.assertEqual(record["upstream_source_commit"], installer.SOURCE_COMMIT)
        self.assertEqual(record["environment"]["CADD_PYTHON2"], str(self.args.runtime_dir / "bin/cadd-python2"))
        self.assertTrue(record["environment"]["CADD_MODEL"].endswith("/GRCh38/CADDv1.7-GRCh38.mod"))
        self.assertEqual(list(program.parent.iterdir()), [program])
        before = installer.cache.identity(program)
        installer.python_calculation(program, source, self.args.runtime_dir, "GRCh38", verify_only=True)
        self.assertEqual(installer.cache.identity(program), before)

    def test_python_route_does_not_overwrite_existing_different_file_or_symlink(self):
        program = self.root / "selected-program"
        program.write_bytes(b"manual executable")
        program.chmod(0o755)
        with self.assertRaisesRegex(RuntimeError, "differs"):
            installer.runtime.install_python_calculation(program)
        self.assertEqual(program.read_bytes(), b"manual executable")
        link = self.root / "selected-link"
        link.symlink_to(program)
        with self.assertRaisesRegex(RuntimeError, "regular frozen file"):
            installer.runtime.install_python_calculation(link)
        self.assertTrue(link.is_symlink())

    def test_python_route_verify_only_requires_existing_exact_executable(self):
        program = self.root / "uninstalled/program"
        with self.assertRaisesRegex(RuntimeError, "not installed"):
            installer.runtime.install_python_calculation(program, verify_only=True)
        self.assertFalse(program.parent.exists())

    def test_python_route_failed_publication_leaves_no_partial_program(self):
        program = self.root / "installed/program"
        with patch.object(installer.runtime.os, "link", side_effect=OSError("publication failed")), \
                self.assertRaisesRegex(OSError, "publication failed"):
            installer.runtime.install_python_calculation(program)
        self.assertFalse(program.exists())
        self.assertEqual(list(program.parent.iterdir()), [])

    def test_missing_calculation_model_is_not_replaced_by_annotation_readiness(self):
        with self.assertRaisesRegex(RuntimeError, "Missing original calculation asset"):
            installer.calculation_resources(self.root / "source", "GRCh38")

    def test_verify_only_requires_completed_member_receipt(self):
        config = {"cadd_version": "v1.7", "cadd_GRCh38_anno_url": "https://example.invalid/GRCh38.tar.gz",
                  "cadd_GRCh38_anno_md5": "1" * 32}
        with patch.object(installer.cache, "annotations", side_effect=AssertionError("no acquisition")), \
                self.assertRaisesRegex(RuntimeError, "member verification is incomplete"):
            installer.annotation_resources(config, "GRCh38", self.root / "source",
                                           self.root / "cache", verify_only=True)

    def test_fresh_reviewed_archive_uses_pinned_git_without_data_or_image_install(self):
        source = self.root / "fresh-source"
        config = json.dumps({"cadd_zip_download_url": self.pins["source_archive_url"]})
        with patch.object(installer, "checkout", return_value=source) as checkout, \
                patch.object(installer, "install", side_effect=AssertionError("source preparation only")), \
                patch("sys.stdout", io.StringIO()):
            status = installer.main(["--config-json", config, "--source-only", "--source", str(source)])
        self.assertEqual(status, 0)
        checkout.assert_called_once_with(source, self.pins["source_url"], self.pins["source_commit"],
                                         installer.SOURCE_PATHS, False)

    def test_custom_archive_is_not_silently_replaced_by_reviewed_git(self):
        source = self.root / "custom-source"
        config = json.dumps({"cadd_zip_download_url": "https://example.invalid/custom-cadd.zip"})
        with patch.object(installer, "checkout") as checkout, patch("sys.stderr", io.StringIO()) as errors:
            status = installer.main(["--config-json", config, "--source-only", "--source", str(source)])
        self.assertEqual(status, 1)
        checkout.assert_not_called()
        self.assertFalse(source.exists())
        self.assertIn("custom archive/source bytes are retained", errors.getvalue())

    def test_existing_non_git_archive_tree_is_retained_without_reinitialization(self):
        source = self.root / "unpacked-original"
        source.mkdir()
        script = source / "CADD.sh"
        script.write_bytes(b"retained non-Git original source\n")
        before = installer.cache.identity(script)
        error = subprocess.CalledProcessError(128, ["git", "rev-parse", "HEAD"])
        with patch.object(installer.subprocess, "check_output", side_effect=error), \
                patch.object(installer, "run") as mutation, self.assertRaises(subprocess.CalledProcessError):
            installer.checkout(source, self.pins["source_url"], self.pins["source_commit"], installer.SOURCE_PATHS)
        mutation.assert_not_called()
        self.assertFalse((source / ".git").exists())
        self.assertEqual(installer.cache.identity(script), before)


if __name__ == "__main__":
    unittest.main()
