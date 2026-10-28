"""Shell wiring fixtures; no image acquisition, model inference or cache access."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class CaddRuntimeShellTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="cadd_shell_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temp.cleanup)
        self.config = Path(self.temp.name) / "config with spaces.json"
        self.config.write_text('{"assembly":"hg38","cache_root":"/fixture/cache with spaces"}')

    def shell(self, body, **variables):
        script = r'''
            source "$1" >/dev/null || exit
            yq() { cat "$2"; }
            python() {
                printf '%s\n' "$@" > "$TEST_ARGUMENTS"
                return "${TEST_HELPER_EXIT:-0}"
            }
        ''' + body
        args = self.config.parent / "arguments.txt"
        result = subprocess.run(
            ["bash", "--noprofile", "--norc", "-c", script, "test",
             str(ROOT / "scripts/install_utils.sh"), str(self.config)],
            env=dict(os.environ, TEST_ARGUMENTS=str(args), **variables),
            text=True, capture_output=True, timeout=20,
        )
        return result, args.read_text().splitlines() if args.exists() else []

    def test_default_both_preserves_json_and_paths(self):
        result, args = self.shell('CADD_runtime_install "$2"')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(args[0], str(ROOT / "scripts/cadd_runtime_install.py"))
        self.assertEqual(args[1], "--config-json")
        self.assertEqual(json.loads(args[2]), json.loads(self.config.read_text()))
        self.assertEqual(args[3:], ["--stage", "both"])

    def test_independent_stage_routes(self):
        for stage in ("annotate", "calculate"):
            with self.subTest(stage=stage):
                result, args = self.shell('CADD_runtime_install "$2" "$TEST_STAGE"', TEST_STAGE=stage)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(args[-2:], ["--stage", stage])

    def test_helper_failure_propagates_without_errexit(self):
        result, args = self.shell('CADD_runtime_install "$2" calculate', TEST_HELPER_EXIT="71")
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue(args)

    def test_bad_stage_or_config_does_not_call_helper(self):
        for body in ('CADD_runtime_install "$2" unknown', 'CADD_runtime_install /absent/config'):
            with self.subTest(body=body):
                result, args = self.shell(body)
                self.assertEqual(result.returncode, 2)
                self.assertEqual(args, [])

    def test_yaml_failure_does_not_call_helper(self):
        result, args = self.shell('yq() { return 72; }; CADD_runtime_install "$2"')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(args, [])

    def main_route(self, code, profile="/fixture/qualified profile.json"):
        return self.shell(r'''
            # Substitute expensive installer stages; retain real main ordering
            # and runtime hook. This is a shell fixture, not an installation.
            while read -r _ _ name; do
                case "$name" in
                    main_install|CADD_runtime_install|python|yq) ;;
                    *) eval "$name() { return 0; }" ;;
                esac
            done < <(declare -F)
            PRIVA_INSTALL_TMPDIR="$(dirname "$2")"
            PRIVA_INSTALL_ENV_NAME=fixture
            PRIVA_SMOKE_SKIP_VEP=1
            read_yaml() {
                case "$2" in
                    ref_genome) printf '%s\n' "$1" ;;
                    assembly) printf 'hg38\n' ;;
                    cadd_profile) printf '%s\n' "$TEST_CADD_PROFILE" ;;
                esac
            }
            CADD_install() { printf 'CACHE_STAGE\n'; }
            ClinGen_deploy() { printf 'AFTER_RUNTIME\n'; }
            main_install "$2"
        ''', TEST_HELPER_EXIT=str(code), TEST_CADD_PROFILE=profile)

    def test_main_calls_both_after_cache_install(self):
        result, args = self.main_route(0)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("CACHE_STAGE", result.stdout)
        self.assertIn("AFTER_RUNTIME", result.stdout)
        self.assertEqual(args[-2:], ["--stage", "both"])

    def test_main_skips_runtime_without_selected_profile(self):
        # Upstream CADD.sh needs none of the shared-VEP runtime; yq reports an
        # absent key as null. A failing helper would stop main_install.
        for profile in ("", "null"):
            with self.subTest(profile=profile):
                result, args = self.main_route(73, profile)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("AFTER_RUNTIME", result.stdout)
                self.assertEqual(args, [])

    def test_main_stops_on_runtime_install_failure(self):
        result, args = self.main_route(73)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("CACHE_STAGE", result.stdout)
        self.assertNotIn("AFTER_RUNTIME", result.stdout)
        self.assertTrue(args)


if __name__ == "__main__":
    unittest.main()
