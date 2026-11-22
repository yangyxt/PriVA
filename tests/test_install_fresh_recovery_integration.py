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

import fake_priva_binary


ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "scripts/installer_prerequisites.py"
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


if __name__ == "__main__":
    unittest.main()
