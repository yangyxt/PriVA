"""Named CADD software layout and subprocess boundaries; no solver or models.

The tiny interpreter fixtures exercise routing only. The retained old wrapper
is pinned source, not a benchmark result or a scientific acceptance artifact.
"""
import importlib.util
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("cadd_named_installer", ROOT / "scripts/cadd_runtime_install.py")
I = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(I)
R = I.runtime


class NamedCaddEnvironments(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="cadd_named_",
            dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.manager_root = self.root / "Mamba with spaces"
        self.main_prefix = self.manager_root / "envs/priva_acmg"
        self.folder = self.manager_root / "share/priva/cadd-1.7.2"
        self.source = self.root / "cache/CADD/CADD-scripts-1.7.2"
        self.config = {"cache_root": str(self.root / "cache"), "cadd_base_dir": str(self.source)}
        self.pins = I.load_pins(ROOT / "scripts/cadd_runtime_pins.json")
        environment = patch.dict(os.environ, CONDA_PREFIX=str(self.main_prefix),
                                 MAMBA_ROOT_PREFIX=str(self.manager_root))
        environment.start()
        self.addCleanup(environment.stop)
        for item in self.pins["native"]["environments"].values():
            path = self.source / item["environment_file"]
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("channels: [bioconda, conda-forge]\ndependencies: [python]\n")

    def interpreter(self, prefix):
        (prefix / "bin").mkdir(parents=True)
        (prefix / "bin/python").symlink_to(sys.executable)
        metadata = prefix / "conda-meta"
        metadata.mkdir()
        (metadata / "history").write_text("tiny routing fixture; no packages installed\n")
        (metadata / "fixture.json").write_text('{"test_fixture_only":true}\n')

    def complete(self, folder, name="cadd-python2", legacy=False):
        if not legacy:
            I.prepare_native_prefixes(folder)
        prefix = R.native_environment_prefix(folder, name)
        self.interpreter(prefix)
        state = folder / "native-state"
        state.mkdir(parents=True, exist_ok=True)
        plan = I.native_environment_plan(folder, self.source, name, self.pins)
        I.publish_json(state / (name + ".plan.json"), plan)
        I.publish_json(state / (name + ".ready.json"), {
            "plan": plan, "conda_metadata_sha256": {"fixture.json": I.sha256(prefix / "conda-meta/fixture.json")},
            "status": "installed_not_scientifically_qualified"})
        if legacy:
            wrapper = folder / "bin/cadd-native-exec"
            wrapper.parent.mkdir()
            wrapper.write_bytes((ROOT / "tests/fixtures/cadd_native_exec_legacy.sh").read_bytes())
            wrapper.chmod(0o755)
        R.install_cadd_entry_points(folder, native=True, names=(name,))
        return prefix

    def identity(self, folder):
        return {str(path.relative_to(folder)): (path.lstat().st_ino, path.lstat().st_mtime_ns,
                os.readlink(path) if path.is_symlink() else path.read_bytes())
                for path in folder.rglob("*") if path.is_file() or path.is_symlink()}

    def test_fresh_prefixes_are_named_beside_priva_and_keep_all_four_tools(self):
        I.prepare_native_prefixes(self.folder)
        for name in R.CADD_ENVIRONMENTS:
            self.assertEqual(R.native_environment_prefix(self.folder, name), self.manager_root / "envs" / name)
        self.assertFalse((self.root / "cache/runtime").exists())
        self.assertFalse((self.folder / "mamba").exists())
        self.assertEqual(I.native_root_prefix(), self.manager_root)

    def test_selected_manager_info_is_read_only_fallback_without_active_root(self):
        with patch.dict(os.environ, {"PATH": os.environ["PATH"]}, clear=True), \
                patch.object(I.subprocess, "run", return_value=subprocess.CompletedProcess(
                    [], 0, json.dumps({"root_prefix": str(self.manager_root)}), "")) as query:
            self.assertEqual(I.native_root_prefix("/selected/micromamba"), self.manager_root)
        self.assertEqual(query.call_args.args[0], ["/selected/micromamba", "info", "--json"])
        self.assertEqual(query.call_args.kwargs["timeout"], 30)

    def test_new_runtime_under_cache_manager_root_is_refused(self):
        nested = self.root / "cache/mamba/envs/priva_acmg"
        with patch.dict(os.environ, CONDA_PREFIX=str(nested)), \
                self.assertRaisesRegex(RuntimeError, "outside cache_root"):
            I.installation_paths(self.config, I.arguments(["--config-json", "{}"]))
        self.assertFalse(self.folder.exists())

    def test_completed_legacy_runtime_is_selected_and_reused_without_mutation(self):
        legacy = self.root / "cache/runtime/cadd-1.7.2"
        prefix = self.complete(legacy, legacy=True)
        before = self.identity(legacy)
        with patch.object(I, "native_root_prefix", side_effect=AssertionError("no new layout")), \
                patch.object(I, "run", side_effect=AssertionError("no solve or install")):
            args = I.arguments(["--config-json", "{}"])
            self.assertEqual(I.installation_paths(self.config, args)[2], legacy)
            record, _ = I.native_environments(legacy, self.source, ("cadd-python2",), self.root / "tmp",
                                               self.pins, verify_only=True)
            exported = I.installed_prefixes(self.config, args, self.pins)
        self.assertEqual(record["environments"]["cadd-python2"]["prefix"], str(prefix))
        self.assertEqual(exported["environments"]["cadd-python2"]["prefix"], str(prefix))
        self.assertEqual(before, self.identity(legacy))
        self.assertFalse((legacy / "native-prefixes").exists())
        self.assertEqual(record["wrappers"]["sha256"], R.LEGACY_NATIVE_WRAPPER_SHA256)

    def test_unowned_and_corrupt_old_runtimes_never_hide_behind_fresh_layout(self):
        legacy = self.root / "cache/runtime/cadd-1.7.2"
        legacy.mkdir(parents=True)
        marker = legacy / "user-file"
        marker.write_bytes(b"user data")
        args = I.arguments(["--config-json", "{}"])
        with self.assertRaisesRegex(RuntimeError, "not an owned"):
            I.installation_paths(self.config, args)
        state = legacy / "native-state"
        state.mkdir()
        (state / "cadd-python2.ready.json").write_text('{"plan":{}}\n')
        self.assertEqual(I.installation_paths(self.config, args)[2], legacy)
        with patch.object(I, "run", side_effect=AssertionError("no solve")), \
                self.assertRaisesRegex(RuntimeError, "plan differs"):
            I.native_environments(legacy, self.source, ("cadd-python2",), self.root / "tmp", self.pins)
        self.assertEqual(marker.read_bytes(), b"user data")
        self.assertFalse(self.folder.exists())

    def test_unmanaged_named_environment_is_retained_without_solver(self):
        prefix = self.manager_root / "envs/cadd-python2"
        self.interpreter(prefix)
        before = self.identity(prefix)
        with patch.object(I, "run", side_effect=AssertionError("no solve")), \
                self.assertRaisesRegex(RuntimeError, "Unmanaged existing"):
            I.native_environments(self.folder, self.source, ("cadd-python2",), self.root / "tmp", self.pins)
        self.assertEqual(before, self.identity(prefix))

    def test_prefixes_cli_reports_only_complete_real_routes_without_writes(self):
        prefix = self.complete(self.folder)
        before = self.identity(self.folder)
        output = io.StringIO()
        with patch("sys.stdout", output), patch.object(I, "run", side_effect=AssertionError("no install")):
            code = I.main(["prefixes", "--config-json", json.dumps(self.config)])
        self.assertEqual(code, 0)
        self.assertEqual(json.loads(output.getvalue()), {"runtime_dir": str(self.folder), "backend": "native", "environments": {
            "cadd-python2": {"prefix": str(prefix),
                "ready_receipt": str(self.folder / "native-state/cadd-python2.ready.json")}}})
        self.assertEqual(self.identity(self.folder), before)

    def test_prefixes_refuses_changed_package_plan_and_unfinished_tools(self):
        prefix = self.complete(self.folder)
        args = I.arguments(["--config-json", "{}"])
        metadata = prefix / "conda-meta/fixture.json"
        original = metadata.read_bytes()
        metadata.write_bytes(b"changed")
        with self.assertRaisesRegex(RuntimeError, "metadata changed"):
            I.installed_prefixes(self.config, args, self.pins)
        metadata.write_bytes(original)
        planned = self.folder / "native-state/cadd-esm-python.plan.json"
        planned.write_text("{}\n")
        with self.assertRaisesRegex(RuntimeError, "unfinished plan"):
            I.installed_prefixes(self.config, args, self.pins)

    def test_explicit_existing_container_reports_no_host_envs_without_solver(self):
        folder = self.root / "legacy-container"
        folder.mkdir()
        image = folder / "cadd-upstream-0.1.1.sif"
        image.write_bytes(b"tiny image fixture; no scientific content\n")
        (folder / "image.origin.json").write_text(json.dumps({
            "image_source": self.pins["image"], "image_sha256": I.sha256(image)}))
        R.install_cadd_entry_points(folder)
        before = self.identity(folder)
        args = I.arguments(["--config-json", "{}", "--runtime-dir", str(folder)])
        with patch.object(I, "native_root_prefix", side_effect=AssertionError("no native layout")), \
                patch.object(I, "run", side_effect=AssertionError("no solve or runtime check")):
            self.assertEqual(I.installed_prefixes(self.config, args, self.pins),
                {"runtime_dir": str(folder), "backend": "container", "environments": {}})
        self.assertEqual(before, self.identity(folder))
        (folder / "image.origin.json").unlink()
        with self.assertRaisesRegex(RuntimeError, "origin receipt"):
            I.installed_prefixes(self.config, args, self.pins)

    def test_prefix_metadata_tampering_and_missing_route_never_fall_back(self):
        self.complete(self.folder)
        metadata = self.folder / "native-prefixes/cadd-python2.prefix"
        original = metadata.read_bytes()
        metadata.chmod(0o644)  # Intentional corruption of this test fixture only.
        for payload in (b"relative/envs/cadd-python2\n", original + b"extra\n",
                        b"/tmp/../envs/cadd-python2\n", original.rstrip(b"\n")):
            metadata.write_bytes(payload)
            with self.subTest(payload=payload), self.assertRaisesRegex(RuntimeError, "Invalid"):
                R.native_environment_prefix(self.folder, "cadd-python2")
            result = subprocess.run([str(self.folder / "bin/cadd-python2"), "-V"], capture_output=True)
            self.assertEqual(result.returncode, 69)
        metadata.unlink()
        with self.assertRaisesRegex(RuntimeError, "incomplete"):
            R.native_environment_prefix(self.folder, "cadd-python2")

    def test_legacy_wrapper_cannot_be_used_for_named_prefixes(self):
        self.complete(self.folder)
        wrapper = self.folder / "bin/cadd-native-exec"
        wrapper.write_bytes((ROOT / "tests/fixtures/cadd_native_exec_legacy.sh").read_bytes())
        with self.assertRaisesRegex(RuntimeError, "wrapper differs"):
            I.native_environments(self.folder, self.source, ("cadd-python2",), self.root / "tmp", self.pins,
                                  verify_only=True)

    def parent(self, child_code, background=False):
        prefix = self.complete(self.folder)
        script = r'''
            before_prefix=$CONDA_PREFIX before_path=$PATH before_home=$PYTHONHOME before_pythonpath=$PYTHONPATH
            "$1" -u -c "$2" ''' + (r'''&
            child=$!
            wait "$child"
        ''' if background else "\n") + r'''
            status=$?
            [[ $CONDA_PREFIX == "$before_prefix" && $PATH == "$before_path" &&
               $PYTHONHOME == "$before_home" && $PYTHONPATH == "$before_pythonpath" ]] || exit 90
            printf 'PARENT_UNCHANGED:%s\n' "$status"
            exit "$status"
        '''
        environment = dict(os.environ, CONDA_PREFIX=str(self.main_prefix),
            PATH=str(self.main_prefix / "bin") + os.pathsep + os.environ["PATH"],
            PYTHONHOME="/parent-only/python", PYTHONPATH="/parent-only/modules")
        process = subprocess.Popen(["bash", "--noprofile", "--norc", "-c", script, "fixture",
            str(self.folder / "bin/cadd-python2"), child_code], env=environment,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return process, prefix

    def test_child_success_and_failure_restore_parent_without_activation(self):
        for code in (0, 71):
            with self.subTest(code=code):
                if self.folder.exists():
                    # Reuse the completed routing fixture in a second parent.
                    with patch.object(self, "complete", return_value=self.manager_root / "envs/cadd-python2"):
                        process, prefix = self.parent(self.child_code(code))
                else:
                    process, prefix = self.parent(self.child_code(code))
                output, error = process.communicate(timeout=15)
                self.assertEqual(process.returncode, code, error)
                record, parent = output.splitlines()
                record = json.loads(record)
                self.assertEqual(record["conda_prefix"], str(prefix))
                self.assertEqual(record["path_first"], str(prefix / "bin"))
                self.assertIsNone(record["pythonhome"])
                self.assertIsNone(record["pythonpath"])
                self.assertTrue(record["no_user_site"])
                self.assertTrue(record["dont_write_bytecode"])
                self.assertEqual(parent, "PARENT_UNCHANGED:" + str(code))

    @staticmethod
    def child_code(code):
        return ("import json,os,sys; print(json.dumps({'conda_prefix':os.environ.get('CONDA_PREFIX'),"
                "'path_first':os.environ['PATH'].split(os.pathsep)[0],"
                "'pythonhome':os.environ.get('PYTHONHOME'),'pythonpath':os.environ.get('PYTHONPATH'),"
                "'no_user_site':sys.flags.no_user_site,'dont_write_bytecode':sys.dont_write_bytecode}));"
                "sys.exit(" + str(code) + ")")

    def test_child_cancellation_returns_signal_status_and_preserves_parent(self):
        process, _ = self.parent("import os,time; print(os.getpid(),flush=True); time.sleep(60)", background=True)
        child = None
        try:
            child = int(process.stdout.readline().strip())
            os.kill(child, signal.SIGTERM)  # Only this test's explicitly owned child.
            output, error = process.communicate(timeout=15)
            self.assertEqual(process.returncode, 143, error)
            self.assertEqual(output.strip(), "PARENT_UNCHANGED:143")
        finally:
            if process.poll() is None:
                if child is not None:
                    try:
                        os.kill(child, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                process.terminate()
                process.communicate(timeout=10)


if __name__ == "__main__":
    unittest.main()
