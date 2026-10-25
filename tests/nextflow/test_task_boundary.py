"""Synthetic task-boundary preflight tests; no genomic data or models."""
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    'priva_task_boundary', Path(__file__).resolve().parents[2] / 'bin/priva_task_boundary.py')
BOUNDARY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BOUNDARY)


class BoundaryTests(unittest.TestCase):
    def setUp(self):
        self.original_cwd = Path.cwd()
        self.fixture = tempfile.TemporaryDirectory()
        self.root = Path(self.fixture.name)
        self.task = self.root / 'task'
        self.task.mkdir()
        for name in ('.command.sh', '.command.run'):
            (self.task / name).write_text('true\n')
        os.chdir(self.task)

    def tearDown(self):
        os.chdir(self.original_cwd)
        self.fixture.cleanup()

    def test_controls_required(self):
        (self.task / '.command.run').unlink()
        with self.assertRaises(BOUNDARY.BoundaryError):
            BOUNDARY.check_task(self.task)

    def test_wrong_working_directory_refused(self):
        os.chdir(self.root)
        with self.assertRaises(BOUNDARY.BoundaryError):
            BOUNDARY.check_task(self.task)

    def test_root_symlink_refused(self):
        link = self.root / 'link'
        link.symlink_to(self.task, target_is_directory=True)
        with self.assertRaises(BOUNDARY.BoundaryError):
            BOUNDARY.check_task(link)

    def test_staged_symlinks_are_not_walked(self):
        outside = self.root / 'outside'
        outside.mkdir()
        os.mkfifo(outside / 'unopened-fifo')
        (self.task / 'input').symlink_to(outside, target_is_directory=True)
        self.assertEqual(BOUNDARY.check_task(self.task), self.task)

    def test_hardlink_alias_refused(self):
        source = self.root / 'cache'
        source.write_text('unchanged\n')
        os.link(source, self.task / 'input')
        with self.assertRaisesRegex(BOUNDARY.BoundaryError, 'hard-linked'):
            BOUNDARY.check_task(self.task)
        self.assertEqual(source.read_text(), 'unchanged\n')

    def test_control_symlink_refused(self):
        control = self.task / '.command.sh'
        control.unlink()
        control.symlink_to(self.task / '.command.run')
        with self.assertRaises(BOUNDARY.BoundaryError):
            BOUNDARY.check_task(self.task)

    def test_fifo_refused(self):
        os.mkfifo(self.task / 'fifo')
        with self.assertRaises(BOUNDARY.BoundaryError):
            BOUNDARY.check_task(self.task)

    def test_argv_retains_exact_command(self):
        command = ['/bin/bash', '-e', '-o', 'pipefail', '.command.sh']
        result = BOUNDARY.command_argv(self.task, command, '/usr/bin/bwrap')
        self.assertEqual(result[-len(command):], command)
        self.assertIn('--unshare-net', result)
        self.assertIn('--unshare-pid', result)
        self.assertIn('--ro-bind', result)

    def test_scratch_symlink_refused(self):
        (self.task / '.priva_task_tmp').symlink_to(self.root, target_is_directory=True)
        with self.assertRaises(BOUNDARY.BoundaryError):
            BOUNDARY.command_argv(self.task, ['true'], '/usr/bin/bwrap')

    def test_empty_command_refused(self):
        with self.assertRaises(BOUNDARY.BoundaryError):
            BOUNDARY.command_argv(self.task, [], '/usr/bin/bwrap')


if __name__ == '__main__':
    unittest.main()
