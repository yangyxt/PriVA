"""CADD runtime publication gates, offline checks and unchanged asset views."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("cadd_runtime_resources", ROOT / "scripts/deploy_cluster_resources.py")
resources = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(resources)


class CaddRuntimeDeployment(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="cadd_runtime_",
                                                     dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.runtime = self.root / "runtime"
        self.runtime.mkdir()
        self.source = self.root / "source"
        for assembly in ("GRCh37", "GRCh38"):
            folder = self.source / f"data/models/{assembly}"
            folder.mkdir(parents=True)
            for name in (f"CADDv1.7-{assembly}.mod", f"conversionTable_CADDv1.7-{assembly}.txt"):
                (folder / name).write_bytes(b"original small calculation asset\n")

    def seed_grch37_annotation_assets(self):
        folder = self.source / "data/annotations/GRCh37_v1.7"
        names = [*(f"esm/{name}.pt" for name in resources.ESM_MODELS), "esm/pep.110.fa",
                 "regseq/Hyperopt400InclNegatives.json", "regseq/Hyperopt400InclNegatives.h5",
                 "regseq/reference.fa", "regseq/reference.fa.fai", "regseq/reference.fa.genome"]
        for name in names:
            path = folder / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(("original " + name + "\n").encode())

    def test_entry_points_install_exact_source_and_reuse_without_rewriting(self):
        record = resources.install_cadd_entry_points(self.runtime)
        wrapper = self.runtime / "bin/cadd-container-exec"
        self.assertEqual(wrapper.read_bytes(), (ROOT / "scripts/cadd_container_exec.sh").read_bytes())
        self.assertTrue(os.access(wrapper, os.X_OK))
        for name in resources.CADD_ENTRY_POINTS:
            self.assertEqual(os.readlink(self.runtime / "bin" / name), wrapper.name)
        before = {path.name: path.lstat().st_mtime_ns for path in (self.runtime / "bin").iterdir()}
        self.assertEqual(resources.install_cadd_entry_points(self.runtime), record)
        self.assertEqual(before, {path.name: path.lstat().st_mtime_ns for path in (self.runtime / "bin").iterdir()})

    def test_different_existing_wrapper_is_retained_and_rejected(self):
        (self.runtime / "bin").mkdir()
        wrapper = self.runtime / "bin/cadd-container-exec"
        wrapper.write_text("unreviewed wrapper\n")
        with self.assertRaisesRegex(RuntimeError, "differs"):
            resources.install_cadd_entry_points(self.runtime)
        self.assertEqual(wrapper.read_text(), "unreviewed wrapper\n")

    def test_conflicting_route_is_rejected_before_any_entry_point_is_published(self):
        binary = self.runtime / "bin"
        binary.mkdir()
        (binary / "cadd-esm-python").symlink_to("different-wrapper")
        with self.assertRaisesRegex(RuntimeError, "different route"):
            resources.install_cadd_entry_points(self.runtime)
        self.assertFalse((binary / "cadd-container-exec").exists())
        self.assertEqual(os.readlink(binary / "cadd-esm-python"), "different-wrapper")

    def test_regular_file_cannot_impersonate_a_named_entry_point(self):
        binary = self.runtime / "bin"
        binary.mkdir()
        (binary / "cadd-python2").write_text("unrelated executable\n")
        with self.assertRaisesRegex(RuntimeError, "expected symlink"):
            resources.install_cadd_entry_points(self.runtime)
        self.assertFalse((binary / "cadd-container-exec").exists())

    def test_asset_inventory_keeps_missing_grch38_separate_from_available_models(self):
        self.seed_grch37_annotation_assets()
        record = resources.cadd_asset_readiness(self.source, self.runtime)
        self.assertEqual(len(record["calculation_models"]), 4)
        self.assertEqual(record["annotation_assets"]["GRCh37"]["status"], "readable_not_scientifically_qualified")
        self.assertEqual(record["annotation_assets"]["GRCh38"]["status"], "incomplete")
        self.assertTrue(any(path.endswith("mmsplice/homo_sapiens.110.gtf")
                            for path in record["annotation_assets"]["GRCh38"]["missing"]))
        self.assertEqual(record["esm_hub"]["status"], "ready")
        self.assertFalse(record["esm_hub"]["downloaded"])
        for item in record["esm_hub"]["checkpoints"]:
            self.assertEqual(os.readlink(item["link"]), item["target"])
            self.assertEqual(Path(item["link"]).stat().st_ino, Path(item["target"]).stat().st_ino)
        self.assertEqual(resources.cadd_asset_readiness(self.source, self.runtime), record)

    def test_absent_annotation_assets_do_not_manufacture_hub_or_readiness(self):
        record = resources.cadd_asset_readiness(self.source, self.runtime)
        self.assertEqual(record["esm_hub"]["status"], "unavailable")
        self.assertFalse((self.runtime / "esm-hub").exists())
        self.assertTrue(all(group["status"] == "incomplete" for group in record["annotation_assets"].values()))

    def test_existing_hub_checkpoint_cannot_be_repointed(self):
        self.seed_grch37_annotation_assets()
        checkpoints = self.runtime / "esm-hub/checkpoints"
        checkpoints.mkdir(parents=True)
        link = checkpoints / f"{resources.ESM_MODELS[0]}.pt"
        link.symlink_to(self.root / "different-model.pt")
        with self.assertRaisesRegex(RuntimeError, "different target"):
            resources.cadd_asset_readiness(self.source, self.runtime)
        self.assertEqual(os.readlink(link), str(self.root / "different-model.pt"))

    def test_empty_required_model_is_never_readable(self):
        target = self.source / "data/models/GRCh37/CADDv1.7-GRCh37.mod"
        target.write_bytes(b"")
        with self.assertRaisesRegex(RuntimeError, "empty or unreadable"):
            resources.cadd_asset_readiness(self.source, self.runtime)

    def fake_environment(self, command, **kwargs):
        name = Path(command[0]).name
        self.assertLessEqual(kwargs["timeout"], 180)
        if name == "cadd-vep":
            self.assertEqual(command[1:], ["--help"])
            return SimpleNamespace(stdout="ensembl-vep: 110.1\n", stderr="")
        prefix, major, modules = json.loads(command[-1])
        self.assertEqual((prefix, major, tuple(modules)), resources.CADD_ENVIRONMENTS[name])
        self.assertIn("socket.socket.connect = no_network", command[2])
        self.assertNotIn("load_model_and_alphabet(", command[2])
        return SimpleNamespace(stdout=json.dumps({"executable": prefix + "/bin/python", "conda_prefix": prefix}), stderr="")

    def test_every_named_environment_is_bounded_offline_and_separate_from_inference(self):
        with patch.object(resources.subprocess, "run", side_effect=self.fake_environment) as run:
            result = resources.check_cadd_environments(self.runtime, {})
        self.assertEqual(set(result), set(resources.CADD_ENTRY_POINTS))
        self.assertEqual(run.call_count, 5)
        self.assertEqual(result["cadd-vep"]["role"], "reference_only_not_shared_VEP113")

    def test_wrong_named_interpreter_route_is_rejected(self):
        bad = SimpleNamespace(stdout=json.dumps({"executable": "/other/python", "conda_prefix": "/other"}), stderr="")
        with patch.object(resources.subprocess, "run", return_value=bad), self.assertRaisesRegex(RuntimeError, "different route"):
            resources.check_cadd_environments(self.runtime, {})

    def test_vep_requires_its_own_110_banner_not_an_incidental_110_in_help(self):
        def fake(command, **kwargs):
            if Path(command[0]).name == "cadd-vep":
                return SimpleNamespace(stdout="ensembl-vep: 113.2\nDocumentation also mentions110\n", stderr="")
            return self.fake_environment(command, **kwargs)
        with patch.object(resources.subprocess, "run", side_effect=fake), self.assertRaisesRegex(RuntimeError, "VEP110"):
            resources.check_cadd_environments(self.runtime, {})

    def test_nonzero_vep_help_is_not_accepted(self):
        def fake(command, **kwargs):
            if Path(command[0]).name == "cadd-vep":
                raise subprocess.CalledProcessError(255, command, stderr="unsupported option\n")
            return self.fake_environment(command, **kwargs)
        with patch.object(resources.subprocess, "run", side_effect=fake), self.assertRaisesRegex(RuntimeError, "help/version check"):
            resources.check_cadd_environments(self.runtime, {})

    def test_failed_import_stops_without_trying_online_fallback(self):
        with patch.object(resources.subprocess, "run", side_effect=subprocess.CalledProcessError(1, ["python"])) as run:
            with self.assertRaisesRegex(RuntimeError, "offline import/route checks"):
                resources.check_cadd_environments(self.runtime, {})
        self.assertEqual(run.call_count, 1)

    def run_main(self, stage_result):
        state = self.root / "state"
        state.mkdir()
        old = state / "cadd-runtime.ready.json"
        old.write_text('{"slurm_job_id":"1514578","historical":true}\n')
        manifest = self.root / "manifest.json"
        manifest.write_text('{"assembly":"GRCh37"}\n')
        new = state / "cadd-runtime.verified.ready.json"
        arguments = ["deploy_cluster_resources.py", "cadd-runtime", "--manifest", str(manifest),
                     "--cache-root", str(self.root / "cache"), "--software-root", str(self.root / "software"),
                     "--state-dir", str(state), "--tmp-root", str(self.root / "tmp"),
                     "--runtime-receipt", str(new)]
        with patch("sys.argv", arguments), patch.object(resources, "source"), \
                patch.object(resources, "cadd_runtime", **stage_result):
            resources.main()
        return old, new

    def test_successful_recheck_preserves_historical_receipt(self):
        old, new = self.run_main({"return_value": {"environments": {"verified": True}}})
        self.assertEqual(old.read_text(), '{"slurm_job_id":"1514578","historical":true}\n')
        self.assertEqual(json.loads(new.read_text())["status"], "deployed")

    def test_failed_recheck_never_publishes_new_receipt(self):
        with self.assertRaisesRegex(RuntimeError, "missing model"):
            self.run_main({"side_effect": RuntimeError("missing model")})
        self.assertFalse((self.root / "state/cadd-runtime.verified.ready.json").exists())
        self.assertEqual((self.root / "state/cadd-runtime.ready.json").read_text(),
                         '{"slurm_job_id":"1514578","historical":true}\n')


if __name__ == "__main__":
    unittest.main()
