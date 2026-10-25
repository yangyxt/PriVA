"""Installer dispatch fixtures, without Java, downloads or scientific work.

Exercise the real main_install and Nextflow_runtime_install functions. Only the
driver's subprocess and unrelated acquisition stages are stubbed. The separate
bootstrap suite checks driver installation/recovery; these tests establish its
place in main_install, config-derived destination and failure propagation.
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"
PRELUDE = r'''
source "$INSTALLER" >/dev/null || exit
record_call() {
    "$TEST_PYTHON" -c 'import json,os,sys;open(os.environ["CALLS"],"a").write(json.dumps(sys.argv[1:])+"\n")' "$@"
}
# JSON is a YAML subset. Full YAML parsing/bootstrap has separate tests.
read_yaml() {
    "$TEST_PYTHON" -c 'import json,sys; v=json.load(open(sys.argv[1])).get(sys.argv[2]); print("null" if v is None else v)' "$1" "$2"
}
# Intercept only the external driver helper. No real network/JVM in this suite.
bash() {
    [[ $1 == "$SCRIPT_DIR/bootstrap_nextflow_runtime.sh" ]] || return 98
    record_call driver "$@"
    return "${DRIVER_STATUS:-0}"
}
'''
MAIN_FIXTURES = r'''
priva_install_basic_paths() { PRIVA_INSTALL_TMPDIR="$TEST_FOLDER"; }
conda_install_vep() { PRIVA_INSTALL_ENV_NAME=priva_acmg; }
priva_manual_cache_inputs() { record_call manual "$@"; return "${MANUAL_STATUS:-0}"; }
priva_prepare_install_config() { record_call config "$@"; return "${CONFIG_STATUS:-0}"; }
bundled_snapshot_install() { record_call bundles "$@"; return "${BUNDLES_STATUS:-0}"; }
# Stop at the next acquisition stage; no VEP/reference/model assets are used.
reference_genome_install() { record_call reference "$@"; return 73; }
main_install "$CONFIG"
'''


class NextflowDriverInstallTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.cache = self.folder / "cache"
        self.cache.mkdir()
        self.config = self.folder / "selected profile.json"
        self.settings = dict(cache_root=str(self.cache), assembly="hg38",
                             vep_cache_dir=str(self.folder), vep_plugins_dir=str(self.folder),
                             vep_plugins_cachedir=str(self.folder), ref_genome="unused.fa")
        self.calls = self.folder / "calls.jsonl"

    def run_fixture(self, body, **overrides):
        self.config.write_text(json.dumps(self.settings))
        before = self.config.read_bytes()
        self.calls.write_text("")
        env = {key: value for key, value in os.environ.items()
               if not key.startswith(("PRIVA_SMOKE_", "BASH_FUNC_"))
               and key not in ("BASH_ENV", "ENV")}
        env.update(INSTALLER=str(INSTALLER), TEST_PYTHON=sys.executable,
                   TEST_FOLDER=str(self.folder), CONFIG=str(self.config), CALLS=str(self.calls))
        env.update(overrides)
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c", PRELUDE + body],
                                env=env, capture_output=True, text=True, timeout=20)
        self.assertEqual(self.config.read_bytes(), before, "driver must not rewrite the profile")
        calls = [json.loads(line) for line in self.calls.read_text().splitlines()]
        return result, calls

    def expected_driver_call(self):
        return ["driver", str(ROOT / "scripts/bootstrap_nextflow_runtime.sh"),
                "install", "--prefix", str(self.cache / "runtime/priva-nextflow-25.10.4")]

    def test_direct_install_uses_config_root_and_named_helper_interface(self):
        for suffix in ("", "/"):
            with self.subTest(trailing_slash=bool(suffix)):
                self.settings["cache_root"] = str(self.cache) + suffix
                result, calls = self.run_fixture('Nextflow_runtime_install "$CONFIG"\n')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(calls, [self.expected_driver_call()])
                self.assertIn("Nextflow driver ready", result.stderr + result.stdout)

    def test_missing_non_directory_and_relative_root_stop_before_helper(self):
        regular = self.folder / "ordinary-file"
        regular.write_text("retained")
        missing = self.folder / "missing-root"
        for value in (str(missing), str(regular), "relative-cache", None):
            with self.subTest(cache_root=value):
                self.settings["cache_root"] = value
                result, calls = self.run_fixture('Nextflow_runtime_install "$CONFIG"\n')
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(calls, [])
                self.assertFalse(missing.exists())
                self.assertEqual(regular.read_text(), "retained")

    def test_existing_root_alias_uses_its_selected_destination_without_config_changes(self):
        alias = self.folder / "selected-cache-root"
        alias.symlink_to(self.cache, target_is_directory=True)
        self.settings["cache_root"] = str(alias)
        result, calls = self.run_fixture('Nextflow_runtime_install "$CONFIG"\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(calls, [self.expected_driver_call()])
        self.assertTrue(alias.is_symlink())
        self.assertEqual(alias.readlink(), self.cache)

    def test_filesystem_root_alias_stops_before_helper(self):
        alias = self.folder / "filesystem-root"
        alias.symlink_to("/", target_is_directory=True)
        self.settings["cache_root"] = str(alias)
        result, calls = self.run_fixture('Nextflow_runtime_install "$CONFIG"\n')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(calls, [])
        self.assertIn("dedicated storage directory", result.stderr + result.stdout)
        self.assertEqual(alias.readlink(), Path("/"))

    def test_parent_components_follow_existing_cache_resolver_lexical_order(self):
        external = self.folder / "different-physical-parent"
        target = external / "target"
        target.mkdir(parents=True)
        (external / "cache").mkdir()
        alias = self.folder / "linked-child"
        alias.symlink_to(target, target_is_directory=True)
        # Both paths exist. Resolve normpath first, like installer_config_paths,
        # rather than silently selecting different-physical-parent/cache.
        self.settings["cache_root"] = str(alias / ".." / "cache")
        result, calls = self.run_fixture('Nextflow_runtime_install "$CONFIG"\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(calls, [self.expected_driver_call()])

    def test_missing_config_or_wrong_argument_count_does_not_start_helper(self):
        for invocation in ('Nextflow_runtime_install',
                           'Nextflow_runtime_install "$CONFIG.missing"',
                           'Nextflow_runtime_install "$CONFIG" extra'):
            with self.subTest(invocation=invocation):
                result, calls = self.run_fixture(invocation)
                self.assertEqual(result.returncode, 2)
                self.assertEqual(calls, [])

    def test_main_orders_driver_after_preflight_and_before_reference(self):
        result, calls = self.run_fixture(MAIN_FIXTURES)
        self.assertEqual(result.returncode, 1)  # Reference fixture intentionally stops.
        self.assertEqual([row[0] for row in calls], ["manual", "config", "bundles", "driver", "reference"])
        self.assertEqual(calls[3], self.expected_driver_call())
        self.assertEqual(calls[-1][1:], [str(self.config)])

    def test_driver_failure_stops_main_before_reference_or_success(self):
        result, calls = self.run_fixture(MAIN_FIXTURES, DRIVER_STATUS="73")
        self.assertEqual(result.returncode, 1)
        self.assertEqual([row[0] for row in calls], ["manual", "config", "bundles", "driver"])
        self.assertIn("Failed to install the pinned Nextflow driver", result.stderr + result.stdout)
        self.assertNotIn("Nextflow driver ready", result.stderr + result.stdout)

    def test_preflight_failure_prevents_driver_and_resource_work(self):
        for variable, expected in (("MANUAL_STATUS", ["manual"]),
                                   ("CONFIG_STATUS", ["manual", "config"]),
                                   ("BUNDLES_STATUS", ["manual", "config", "bundles"])):
            with self.subTest(failed_stage=variable):
                result, calls = self.run_fixture(MAIN_FIXTURES, **{variable: "73"})
                self.assertEqual(result.returncode, 1)
                self.assertEqual([row[0] for row in calls], expected)


if __name__ == "__main__":
    unittest.main()
