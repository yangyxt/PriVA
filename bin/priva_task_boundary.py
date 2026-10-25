#!/usr/bin/env python3
"""Run a Nextflow task with only its task directory writable on the host.

  python3 bin/priva_task_boundary.py check --task-dir /absolute/task
  python3 bin/priva_task_boundary.py exec --task-dir . -- /bin/bash -e .command.sh

24 October 2026 CT: integration of the synthetic read-only boundary pilot.
This opt-in Linux/Bubblewrap interface does not qualify annotation, scoring or
family results. It never retries a denied write outside the task directory.
Nextflow's driver, scheduler wrapper and publication remain outside this view.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import stat
import sys


class BoundaryError(ValueError):
    """The requested writable view cannot safely represent a task directory."""


def check_task(task_dir: Path) -> Path:
    """Check real Nextflow controls and reject pre-existing writable aliases.

    A read-only cache bind does not protect a hard link reached through a
    different writable bind. Stage symlinks or independent copies instead.
    Walk only this task's physical directories; never follow staged symlinks.
    Other writers must not change the task directory during launch.
    """
    if task_dir.is_symlink():
        raise BoundaryError('Task directory must not be a symlink')
    task = task_dir.resolve(strict=True)
    if not task.is_dir() or task != Path.cwd().resolve():
        raise BoundaryError('Task directory must be the current physical directory')
    for name in ('.command.sh', '.command.run'):
        path = task / name
        try:
            value = path.lstat()
        except FileNotFoundError as error:
            raise BoundaryError('Missing Nextflow task control: ' + name) from error
        if not stat.S_ISREG(value.st_mode) or value.st_nlink != 1:
            raise BoundaryError('Nextflow task control must be an independent regular file: ' + name)
    for current, directories, files in os.walk(task, followlinks=False):
        for name in (*directories, *files):
            path = Path(current) / name
            value = path.lstat()
            if stat.S_ISLNK(value.st_mode) or stat.S_ISDIR(value.st_mode):
                continue
            if not stat.S_ISREG(value.st_mode):
                raise BoundaryError('Unsupported writable task entry: ' + str(path))
            if value.st_nlink != 1:
                raise BoundaryError('Writable task contains a hard-linked file: ' + str(path))
    return task


def private_directory(task: Path, name: str) -> Path:
    """Create private scratch without following an input's directory symlink."""
    path = task / name
    if path.is_symlink() or (path.exists() and not path.is_dir()):
        raise BoundaryError('Private task scratch is not a physical directory: ' + name)
    path.mkdir(exist_ok=True)
    return path


def command_argv(task: Path, command: list[str], executable: str) -> list[str]:
    if not command:
        raise BoundaryError('Pass a command after --')
    temp = private_directory(task, '.priva_task_tmp')
    cache = private_directory(task, '.priva_task_cache')
    # Keep scientific inputs, HOME/configuration and package/model paths intact.
    # Only disposable scratch/cache defaults are redirected. No model is fetched
    # at runtime: network access is disabled inside this task view.
    return [executable, '--die-with-parent', '--new-session', '--unshare-pid',
            '--unshare-net', '--cap-drop', 'ALL', '--ro-bind', '/', '/',
            '--bind', str(task), str(task), '--tmpfs', '/tmp', '--proc', '/proc',
            '--dev', '/dev', '--chdir', str(task), '--setenv', 'TMPDIR', str(temp),
            '--setenv', 'XDG_CACHE_HOME', str(cache), '--setenv',
            'PYTHONDONTWRITEBYTECODE', '1', '--setenv', 'PRIVA_TASK_BOUNDARY',
            'bwrap-readonly-v1', '--', *command]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    procedures = parser.add_subparsers(dest='procedure', required=True)
    for name in ('check', 'exec'):
        command = procedures.add_parser(name)
        command.add_argument('--task-dir', type=Path, required=True)
        if name == 'exec':
            command.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    try:
        task = check_task(args.task_dir)
        if args.procedure == 'check':
            return 0
        command = args.command
        if command[:1] == ['--']:
            command = command[1:]
        executable = shutil.which('bwrap')
        if executable is None:
            raise BoundaryError('Bubblewrap (bwrap) is required for the task_readonly profile')
        os.execv(executable, command_argv(task, command, executable))
    except (BoundaryError, OSError) as error:
        print('PriVA task boundary refused: ' + str(error), file=sys.stderr)
        return 125
    return 125


if __name__ == '__main__':
    raise SystemExit(main())
