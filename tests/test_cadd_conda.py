"""Exercise frontend routing with fixture executables, never package transactions."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("cadd_conda", ROOT / "scripts/cadd_conda.py")
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)


def write_conda(prefix, version="24.9.2", activate=True):
    """A test-only protocol peer for --version/info/config; it cannot install."""
    (prefix / "bin").mkdir(parents=True, exist_ok=True)
    source = f'''#!{sys.executable}
import json, sys
if sys.argv[1:] == ['--version']:
    print('conda {version}')
elif sys.argv[1:] == ['info', '--json']:
    print(json.dumps({{'conda_prefix': {str(prefix)!r}, 'platform': 'linux-64'}}))
elif sys.argv[1:] == ['config', '--get', 'channel_priority', '--json']:
    print('{{"success": true, "get": {{"channel_priority": "flexible"}}}}')
else:
    raise SystemExit(91)
'''
    executable = prefix / "bin/conda"
    executable.write_text(source)
    executable.chmod(0o755)
    if activate:
        (prefix / "bin/activate").write_text("# fixture activation entry\n")
    return executable


class CondaFrontendTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="cadd-conda-fixture-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.tools = self.root / "tools"
        self.tools.mkdir()
        (self.tools / "python").symlink_to(sys.executable)
        self.frontend = self.root / "cache/CADD/conda_frontend"
        self.env = dict(os.environ, PATH=str(self.tools), CONDA_EXE="", MAMBA_EXE="")
        self.manager = self.tools / "micromamba"
        # The fixture manager accepts only create; it records the exact argv.
        self.manager.write_text(f'''#!{sys.executable}
import importlib.util, json, pathlib, sys
pathlib.Path({str(self.root / 'manager.json')!r}).write_text(json.dumps(sys.argv[1:]))
assert sys.argv[1] == 'create', sys.argv
spec=importlib.util.spec_from_file_location('fixture', {str(Path(__file__).resolve())!r})
module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
module.write_conda(pathlib.Path(sys.argv[sys.argv.index('--prefix')+1]))
''')
        self.manager.chmod(0o755)

    def run_helper(self, procedure="check", **kwargs):
        command = [sys.executable, str(ROOT / "scripts/cadd_conda.py"), procedure,
                   "--frontend-prefix", str(self.frontend)]
        if procedure == "install":
            command += ["--manager", str(self.manager)]
        result = subprocess.run(command, env=self.env, capture_output=True, text=True, timeout=30)
        if kwargs.get("success", True):
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
        return result

    def selected(self, result):
        return json.loads(result.stdout)

    def test_compatible_user_conda_on_path_is_reused_without_special_provenance(self):
        executable = write_conda(self.root / "user", version="25.11.0")
        self.env["PATH"] += os.pathsep + str(executable.parent)
        result = self.selected(self.run_helper("install"))
        self.assertEqual(result["executable"], str(executable))
        self.assertEqual(result["path"], self.env["PATH"])
        self.assertFalse((self.root / "manager.json").exists())
        self.assertFalse(self.frontend.parent.exists())

    def test_conda_exe_outside_path_is_discovered_without_python_shadowing(self):
        executable = write_conda(self.root / "user")
        self.env["CONDA_EXE"] = str(executable)
        (executable.parent / "python").write_text("not the scientific Python")
        result = self.selected(self.run_helper())
        self.assertEqual(shutil.which("conda", path=result["path"]), str(executable))
        self.assertEqual(shutil.which("python", path=result["path"]), str(self.tools / "python"))

    def test_managed_frontend_is_found_read_only_without_historical_path(self):
        executable = write_conda(self.frontend)
        result = self.selected(self.run_helper())
        self.assertEqual(result["executable"], str(executable))
        self.assertFalse((self.frontend.parent / ".conda_frontend.lock").exists())

    def test_missing_manager_check_never_creates_or_downloads(self):
        result = self.run_helper(success=False)
        self.assertIn("needs a conda executable", result.stderr)
        self.assertFalse(self.frontend.parent.exists())
        self.assertFalse((self.root / "manager.json").exists())

    def test_old_or_broken_frontend_is_preserved_and_rejected_before_create(self):
        for version, activate in (("23.9.0", True), ("24.9.2", False)):
            with self.subTest(version=version, activate=activate):
                write_conda(self.frontend, version=version, activate=activate)
                activation = self.frontend / "bin/activate"
                if not activate and activation.exists():
                    activation.unlink()
                before = (self.frontend / "bin/conda").read_bytes()
                result = self.run_helper("install", success=False)
                self.assertIn("retained without modification", result.stderr)
                self.assertEqual((self.frontend / "bin/conda").read_bytes(), before)
                self.assertFalse((self.root / "manager.json").exists())

    def test_missing_frontend_create_preserves_all_existing_environments_and_exports(self):
        scientific = self.root / "scientific"
        for name in ("priva_acmg", "minimal", "esm", "mmsplice", "regulatory", "vep"):
            prefix = scientific / name
            (prefix / "conda-meta").mkdir(parents=True)
            (prefix / "conda-meta/package.json").write_text('{"name":"science","version":"1","build":"exact"}')
            (prefix / "conda-meta/history").write_text("original history\n")
            (prefix / "installed.yml").write_text("dependencies:\n  - science=1=exact\n")
            (prefix / "model").write_bytes(b"unchanged model")
        def inventory():
            return {str(p.relative_to(scientific)): (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns)
                    for p in scientific.rglob("*") if p.is_file()}
        before = inventory()
        self.env["CONDA_PREFIX"] = str(scientific / "priva_acmg")
        selected = self.selected(self.run_helper("install"))
        command = json.loads((self.root / "manager.json").read_text())
        self.assertEqual(command, ["create", "--yes", "--prefix", str(self.frontend),
                                  "--override-channels", "--channel", "conda-forge", "conda=24.9.2"])
        self.assertEqual(selected["prefix"], str(self.frontend))
        self.assertEqual(inventory(), before)
        (self.root / "manager.json").unlink()
        self.run_helper("install")
        self.assertFalse((self.root / "manager.json").exists())
        self.assertEqual(inventory(), before)

    def test_frontend_creation_failure_does_not_fall_back_to_scientific_install(self):
        self.manager.write_text(f"#!{sys.executable}\nraise SystemExit(17)\n")
        self.run_helper("install", success=False)
        self.assertFalse(self.frontend.exists())

    def test_creation_uses_the_selected_manager_or_wrapper_for_all_supported_kinds(self):
        for name in ("micromamba-custom", "mamba", "conda"):
            with self.subTest(manager=name):
                manager = self.root / "managers" / name
                manager.parent.mkdir(exist_ok=True)
                manager.write_bytes(self.manager.read_bytes())
                manager.chmod(0o755)
                frontend = self.root / name / "conda_frontend"
                result = subprocess.run([sys.executable, str(ROOT / "scripts/cadd_conda.py"), "install",
                                         "--frontend-prefix", str(frontend), "--manager", str(manager)],
                                        env=self.env, capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertTrue((frontend / "bin/conda").is_file())

    def test_main_stops_before_annotation_when_frontend_is_incomplete(self):
        write_conda(self.frontend, activate=False)
        config = self.root / "config.yaml"
        config.write_text("fixture: true\n")
        shell = '''source "$1" >/dev/null 2>&1 || exit
priva_install_basic_paths() { PRIVA_INSTALL_TMPDIR="$TEST_ROOT"; }
conda_install_vep() { PRIVA_INSTALL_ENV_NAME=priva_acmg; }
priva_cache_install_begin() { PRIVA_CACHE_INSTALL_BUILD_CONFIG=$2; }
priva_manual_cache_inputs() { :; }
priva_prepare_install_config() { :; }
bundled_snapshot_install() { :; }
priva_find_env_manager() { PRIVA_ENV_MANAGER="$TEST_MANAGER"; }
read_yaml() { printf '%s' "$TEST_PARENT"; }
reference_genome_install() { printf 'UNEXPECTED_ANNOTATION'; }
main_install "$2"
'''
        env = dict(self.env, PATH=self.env["PATH"] + ":/usr/bin:/bin", TEST_ROOT=str(self.root),
                   TEST_PARENT=str(self.frontend.parent), TEST_MANAGER=str(self.manager))
        result = subprocess.run(["/bin/bash", "-c", shell, "bash",
                                 str(ROOT / "scripts/install_utils.sh"), str(config)],
                                env=env, capture_output=True, text=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("retained without modification", result.stderr)
        self.assertNotIn("UNEXPECTED_ANNOTATION", result.stdout)
        self.assertFalse((self.root / "manager.json").exists())

    def test_condabin_can_replace_old_path_conda_without_shadowing_python(self):
        old = write_conda(self.root / "old", "23.9.0")
        executable = write_conda(self.frontend)
        (self.frontend / "condabin").mkdir()
        (self.frontend / "condabin/conda").symlink_to(executable)
        (self.frontend / "bin/python").write_text("not the scientific Python")
        self.env["PATH"] += os.pathsep + str(old.parent)
        result = self.selected(self.run_helper())
        self.assertEqual(Path(shutil.which("conda", path=result["path"])).resolve(), executable)
        self.assertEqual(shutil.which("conda", path=result["path"]), result["executable"])
        self.assertEqual(shutil.which("python", path=result["path"]), str(self.tools / "python"))

    def test_shell_routing_keeps_python_and_original_conda_only_cadd_contract(self):
        executable = write_conda(self.frontend)
        shell = '''source "$1" >/dev/null 2>&1 || exit
priva_cadd_conda check --frontend-prefix "$2" || exit
command -v conda
command -v python
'''
        # Only the standard tools needed to source the shell function pool.
        env = dict(self.env, PATH=self.env["PATH"] + ":/usr/bin:/bin")
        result = subprocess.run(["/bin/bash", "-c", shell, "bash",
                                 str(ROOT / "scripts/install_utils.sh"), str(self.frontend)],
                                env=env, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines(), [str(executable), str(self.tools / "python")])


if __name__ == "__main__":
    unittest.main()
