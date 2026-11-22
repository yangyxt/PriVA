"""Isolated fresh-start boundaries and real local-HTTP acquisition recovery.

Manual files and HTTP payloads are synthetic fixtures, not publisher resources
or a full installation. The manual gate, PTY, YAML writer, wget, BGZF/TBI
tools, and MaveDB installer are real. Only costly main-install boundaries are
stubbed, with the precise stop boundary asserted in each test.
"""
import importlib.util
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import os
from pathlib import Path
import pty
import select
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "scripts/installer_prerequisites.py"
INSTALLER = ROOT / "scripts/install_utils.sh"
SPEC = importlib.util.spec_from_file_location("fresh_prerequisites", HELPER)
PREREQS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PREREQS)
KEYS = PREREQS.MANUAL_INPUT_KEYS


class ManualInputIntegration(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_fresh_fixture_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.cache = self.folder / "cache"
        self.cache.mkdir()
        self.config = self.folder / "config.yaml"
        self.selected = {}
        for key in KEYS:
            path = self.folder / (key + ".fixture")
            # Deliberately opaque: this gate must not perform scientific parsing.
            path.write_bytes(b"nonempty synthetic manual-file fixture\n")
            self.selected[key] = str(path)
        self.env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", TEST_CONFIG=str(self.config),
                        TEST_PYTHON=sys.executable, TEST_FOLDER=str(self.folder))
        self.env.pop("PRIVA_SMOKE_SKIP_VEP", None)
        self.env["PRIVA_INSTALL_INSPECT_ROOTS"] = os.pathsep.join((str(ROOT), str(self.folder)))
        self.write_config()

    def write_config(self, selected=None):
        # The full distributed template, with only ordinary fresh-install
        # selections changed; sparse YAML is not an installation interface.
        updates = {"cache_root": str(self.cache), "tmp_dir": str(self.folder),
                   "assembly": "hg38", "vep_installed_version": "116"}
        updates.update({key: "" for key in KEYS})
        updates.update(selected or {})
        text = PREREQS.CONFIG_PATHS.render_config((ROOT / "config.yaml").read_text(), updates)
        self.config.write_text(text + "# unrelated trailing comment stays byte-identical\n")

    def cli(self, *arguments):
        return subprocess.run([sys.executable, str(HELPER), str(self.config),
                               "--stage", "manual-inputs", *arguments], env=self.env,
                              capture_output=True, text=True, timeout=20)

    def assert_unchanged(self, before):
        self.assertEqual(self.config.read_bytes(), before)
        self.assertEqual(list(self.cache.iterdir()), [])
        self.assertEqual(list(self.folder.glob(".config.yaml.*")), [])

    def prompt(self, answers, on_prompt=None):
        """Exercise a real terminal stdin while capturing the helper output."""
        master, slave = pty.openpty()
        process = subprocess.Popen([sys.executable, str(HELPER), str(self.config),
                                    "--stage", "manual-inputs", "--prompt"],
                                   stdin=slave, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   env=self.env)
        os.close(slave)
        output = b""
        deadline = time.monotonic() + 20
        try:
            for key, answer in answers:
                marker = (key + ": ").encode()
                while marker not in output:
                    if process.poll() is not None:
                        output += process.stdout.read()
                        self.fail(f"Helper exited before prompt {key}: {output.decode()}")
                    if time.monotonic() > deadline:
                        self.fail(f"Timed out awaiting prompt {key}: {output.decode()}")
                    if select.select([process.stdout], [], [], 0.1)[0]:
                        output += os.read(process.stdout.fileno(), 65536)
                if on_prompt is not None:
                    on_prompt(key)
                os.write(master, b"\x04" if answer is None else (answer + "\n").encode())
            remaining, _ = process.communicate(timeout=max(1, deadline - time.monotonic()))
            return process.returncode, (output + remaining).decode()
        finally:
            os.close(master)
            if process.poll() is None:
                process.kill()
                process.communicate()
            process.stdout.close()

    def test_missing_batch_inputs_fail_before_writes_and_list_all_keys(self):
        before = self.config.read_bytes()
        result = self.cli()
        self.assertNotEqual(result.returncode, 0)
        for key in KEYS:
            self.assertIn(key, result.stderr)
        self.assert_unchanged(before)

    def test_prepared_inputs_are_read_only_even_when_prompt_requested_without_tty(self):
        self.write_config(self.selected)
        before = self.config.read_bytes(), self.config.stat().st_ino, self.config.stat().st_mtime_ns
        result = self.cli("--prompt")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.config.read_bytes(), self.config.stat().st_ino,
                          self.config.stat().st_mtime_ns), before)
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_prompt_never_reads_piped_input(self):
        before = self.config.read_bytes()
        result = subprocess.run([sys.executable, str(HELPER), str(self.config),
                                 "--stage", "manual-inputs", "--prompt"], input="ignored\n" * 3,
                                text=True, capture_output=True, env=self.env, timeout=20)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("interactive terminal", result.stderr)
        self.assert_unchanged(before)

    def test_explicit_missing_path_is_not_replaced_in_prompt_mode(self):
        selected = dict(self.selected, primateai_prescore=str(self.folder / "absent.fixture"))
        self.write_config(selected)
        before = self.config.read_bytes()
        result = self.cli("--prompt")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("absent.fixture", result.stderr)
        self.assertNotIn("Select completed", result.stderr)
        self.assert_unchanged(before)

    def test_empty_directory_unreadable_and_unfinished_selected_files_fail(self):
        empty = self.folder / "empty.fixture"
        empty.touch()
        unreadable = self.folder / "unreadable.fixture"
        unreadable.write_bytes(b"fixture")
        unreadable.chmod(0)
        self.addCleanup(unreadable.chmod, 0o600)
        partial = self.folder / "input.partial.tsv.gz"
        partial.write_bytes(b"incomplete fixture")
        disguised = self.folder / "completed.fixture"
        disguised.symlink_to(partial)
        for path in (empty, self.folder, unreadable, partial, disguised):
            with self.subTest(path=path.name):
                self.write_config(dict(self.selected, spliceai_snv_prescore=str(path)))
                before = self.config.read_bytes()
                result = self.cli()
                self.assertNotEqual(result.returncode, 0, str(path))
                self.assertIn("spliceai_snv_prescore", result.stderr)
                self.assert_unchanged(before)

    def test_duplicate_key_is_rejected_without_changing_configuration(self):
        self.write_config(self.selected)
        with self.config.open("a") as handle:
            handle.write("primateai_prescore: ignored.fixture\n")
        before = self.config.read_bytes()
        result = self.cli()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Duplicate configuration key", result.stderr)
        self.assert_unchanged(before)

    def test_prompt_flag_cannot_change_other_stages(self):
        before = self.config.read_bytes()
        result = subprocess.run([sys.executable, str(HELPER), str(self.config),
                                 "--stage", "main", "--prompt"], env=self.env,
                                capture_output=True, text=True, timeout=20)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("restricted", result.stderr)
        self.assert_unchanged(before)

    def test_interactive_collection_saves_only_unset_keys_together_and_keeps_comments(self):
        self.write_config({KEYS[0]: self.selected[KEYS[0]]})
        self.config.chmod(0o640)
        before = self.config.read_text()
        snapshots = []
        code, output = self.prompt([(key, self.selected[key]) for key in KEYS[1:]],
                                   lambda key: snapshots.append(self.config.read_text()))
        self.assertEqual(code, 0, output)
        self.assertEqual(snapshots, [before, before])
        self.assertNotIn(KEYS[0] + ": ", output)
        expected = PREREQS.CONFIG_PATHS.render_config(before, {key: self.selected[key] for key in KEYS[1:]})
        self.assertEqual(self.config.read_text(), expected)
        self.assertEqual(self.config.stat().st_mode & 0o777, 0o640)
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_interactive_invalid_later_answer_saves_none(self):
        before = self.config.read_bytes()
        answers = [(key, self.selected[key]) for key in KEYS]
        answers[-1] = (KEYS[-1], str(self.folder / "invalid-last.fixture"))
        code, output = self.prompt(answers)
        self.assertNotEqual(code, 0)
        self.assertIn("invalid-last.fixture", output)
        self.assert_unchanged(before)

    def test_interactive_eof_saves_none(self):
        before = self.config.read_bytes()
        code, output = self.prompt([(KEYS[0], self.selected[KEYS[0]]), (KEYS[1], None)])
        self.assertNotEqual(code, 0)
        self.assertIn("collection ended", output)
        self.assert_unchanged(before)

    def test_interactive_explicit_invalid_selection_fails_before_any_prompt(self):
        self.write_config({KEYS[0]: str(self.folder / "explicit-missing.fixture")})
        before = self.config.read_bytes()
        code, output = self.prompt([])
        self.assertNotEqual(code, 0)
        self.assertIn("explicit-missing.fixture", output)
        self.assertNotIn("Select completed", output)
        self.assert_unchanged(before)

    def test_interactive_collection_preserves_a_concurrent_config_edit(self):
        changed = self.config.read_text() + "# concurrent writer fixture\n"

        def concurrent_edit(key):
            if key == KEYS[0]:
                self.config.write_text(changed)

        code, output = self.prompt([(key, self.selected[key]) for key in KEYS], concurrent_edit)
        self.assertNotEqual(code, 0)
        self.assertIn("Configuration changed", output)
        self.assertEqual(self.config.read_text(), changed)
        self.assertEqual(list(self.cache.iterdir()), [])
        self.assertEqual(list(self.folder.glob(".config.yaml.*")), [])

    def test_interactive_collection_rechecks_previously_selected_file_before_save(self):
        self.write_config({KEYS[0]: self.selected[KEYS[0]]})
        before = self.config.read_bytes()

        def remove_prepared_file(key):
            if key == KEYS[1]:
                Path(self.selected[KEYS[0]]).unlink()

        code, output = self.prompt([(key, self.selected[key]) for key in KEYS[1:]], remove_prepared_file)
        self.assertNotEqual(code, 0)
        self.assertIn(KEYS[0], output)
        self.assert_unchanged(before)

    def main_boundary(self):
        # Explicit boundary stubs: no environment creation, reference/cache
        # preparation, VEP, or main installation occurs in this test.
        script = r'''
source "$1"
priva_install_basic_paths() { printf 'BOUNDARY:basic\n'; PRIVA_INSTALL_TMPDIR="$TEST_FOLDER"; }
conda_install_vep() { printf 'BOUNDARY:environment\n'; PRIVA_INSTALL_ENV_NAME=priva_acmg; }
python() { "$TEST_PYTHON" "$@"; }
priva_prepare_install_config() { printf 'BOUNDARY:prepare\n'; return 73; }
reference_genome_install() { printf 'UNEXPECTED_REFERENCE\n'; return 99; }
VEP_plugins_install() { printf 'UNEXPECTED_WHOLE_TABLE_VEP\n'; return 99; }
main_install "$TEST_CONFIG"
'''
        return subprocess.run(["bash", "--noprofile", "--norc", "-c", script, "fixture", str(INSTALLER)],
                              env=self.env, capture_output=True, text=True, timeout=20)

    def test_main_missing_manual_stops_after_environment_before_prepare(self):
        before = self.config.read_bytes()
        result = self.main_boundary()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual([line for line in result.stdout.splitlines() if line.startswith("BOUNDARY:")],
                         ["BOUNDARY:basic", "BOUNDARY:environment"])
        self.assertIn("Manual publisher inputs", result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assert_unchanged(before)

    def test_main_invalid_manual_stops_before_prepare_without_replacement(self):
        self.write_config(dict(self.selected, primateai_prescore=str(self.folder / "invalid.fixture")))
        before = self.config.read_bytes()
        result = self.main_boundary()
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("BOUNDARY:prepare", result.stdout)
        self.assertIn("invalid.fixture", result.stderr)
        self.assert_unchanged(before)

    def test_main_prepared_manual_reaches_prepare_after_environment(self):
        self.write_config(self.selected)
        before = self.config.read_bytes()
        result = self.main_boundary()
        self.assertNotEqual(result.returncode, 0)  # Intentional stop at prepare.
        self.assertEqual([line for line in result.stdout.splitlines() if line.startswith("BOUNDARY:")],
                         ["BOUNDARY:basic", "BOUNDARY:environment", "BOUNDARY:prepare"])
        self.assertNotIn("UNEXPECTED", result.stdout)
        self.assert_unchanged(before)


@unittest.skipUnless(all(shutil.which(tool) for tool in ("bgzip", "tabix", "wget", "yq")),
                     "Real wget, yq, bgzip and tabix are required")
class LocalHttpMaveDBRecovery(unittest.TestCase):
    def test_real_tools_retry_failed_index_acquisition_without_publishing_partial_data(self):
        # Synthetic BGZF and loopback HTTP only: real acquisition mechanics,
        # explicitly not a publisher download or scientific installation.
        with tempfile.TemporaryDirectory(prefix="priva_http_fixture_", dir=os.environ.get(
                "PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR"))) as temporary:
            folder = Path(temporary)
            served = folder / "synthetic.tsv.gz"
            rows = b"#chrom\tstart\tend\nchr1\t1\t1\nchr1\t2\t2\n"
            with served.open("wb") as handle:
                subprocess.run(["bgzip", "-c"], input=rows, stdout=handle, check=True)
            subprocess.run(["tabix", "-s", "1", "-b", "2", "-e", "3", str(served)], check=True)
            state = {"fail_index": True, "requests": []}

            class FixtureServer(BaseHTTPRequestHandler):
                def do_GET(self):
                    state["requests"].append(self.path)
                    if self.path == "/synthetic.tsv.gz.tbi" and state["fail_index"]:
                        self.send_error(503, "Synthetic first-attempt index failure")
                        return
                    if self.path not in ("/synthetic.tsv.gz", "/synthetic.tsv.gz.tbi"):
                        self.send_error(404)
                        return
                    body = Path(str(served) + (".tbi" if self.path.endswith(".tbi") else "")).read_bytes()
                    self.send_response(200)
                    self.send_header("Content-Length", str(len(body)))
                    self.end_headers()
                    self.wfile.write(body)

                def log_message(self, *args):
                    pass

            server = ThreadingHTTPServer(("127.0.0.1", 0), FixtureServer)
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            try:
                cache = folder / "managed"
                cache.mkdir()
                config = folder / "config.yaml"
                url = f"http://127.0.0.1:{server.server_port}/synthetic.tsv.gz"
                text = PREREQS.CONFIG_PATHS.render_config((ROOT / "config.yaml").read_text(), {
                    "cache_root": str(cache), "tmp_dir": str(folder), "assembly": "hg38",
                    "vep_installed_version": "116", "mavedb_file": "", "mavedb_url": url,
                    "vep_plugins_cachedir": str(cache)})
                config.write_text(text)
                env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", TMPDIR=str(folder),
                           PRIVA_INSTALL_INSPECT_ROOTS=os.pathsep.join((str(ROOT), str(folder))))
                command = ["bash", "--noprofile", "--norc", "-c",
                           'source "$1"; MaveDB_install "$2" "$3"',
                           "fixture", str(INSTALLER), str(config), str(cache)]
                before = config.read_bytes()
                first = subprocess.run(command, env=env, text=True, capture_output=True, timeout=30)
                self.assertNotEqual(first.returncode, 0)
                self.assertIn("Failed to stage", first.stderr + first.stdout)
                target = cache / "MaveDB/MaveDB_variants.hg38.tsv.gz"
                self.assertFalse(target.exists())
                self.assertFalse(Path(str(target) + ".tbi").exists())
                self.assertEqual(config.read_bytes(), before)
                failed_stages = list(target.parent.glob(".download.*"))
                self.assertEqual(len(failed_stages), 1)
                retained = failed_stages[0] / target.name
                self.assertEqual(retained.read_bytes(), served.read_bytes())
                state["fail_index"] = False
                second = subprocess.run(command, env=env, text=True, capture_output=True, timeout=30)
                self.assertEqual(second.returncode, 0, second.stderr + second.stdout)
                self.assertEqual(target.read_bytes(), served.read_bytes())
                self.assertEqual(subprocess.check_output(["tabix", str(target), "chr1:1-2"]),
                                 b"chr1\t1\t1\nchr1\t2\t2\n")
                updated = yaml.safe_load(config.read_text())
                self.assertEqual(updated["mavedb_file"], str(target))
                self.assertEqual(updated["mavedb_url"], url)
                old = yaml.safe_load(before)
                self.assertEqual({key: value for key, value in updated.items() if key != "mavedb_file"},
                                 {key: value for key, value in old.items() if key != "mavedb_file"})
                self.assertEqual(state["requests"], ["/synthetic.tsv.gz", "/synthetic.tsv.gz.tbi"] * 2)
                self.assertEqual(list(target.parent.glob(".download.*")), failed_stages)
                self.assertEqual(retained.read_bytes(), served.read_bytes())
                # A warm rerun reuses exact selected bytes/index without HTTP.
                requests = list(state["requests"])
                identity = target.stat().st_ino, target.stat().st_mtime_ns, target.read_bytes()
                third = subprocess.run(command, env=env, text=True, capture_output=True, timeout=30)
                self.assertEqual(third.returncode, 0, third.stderr + third.stdout)
                self.assertEqual(state["requests"], requests)
                self.assertEqual((target.stat().st_ino, target.stat().st_mtime_ns, target.read_bytes()), identity)
            finally:
                server.shutdown()
                thread.join(timeout=5)
                server.server_close()


if __name__ == "__main__":
    unittest.main()
