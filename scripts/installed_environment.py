#!/usr/bin/env python3
"""Export the packages actually installed in PriVA and CADD environments.

export --config CONFIG --manager MAMBA_OR_CONDA [--name priva_acmg] [--prefix PREFIX]
The portable YAML is written to cache_root/<name>.installed.yml. It describes
the selected environment only; it makes no claim about other CADD environments.
export-installation exports PriVA and every completed environment that
upstream CADD.sh's Snakemake workflow built under cadd_base_dir/envs/conda.
"""
import argparse
import hashlib
import os
from pathlib import Path
import re
import subprocess
import stat
import sys
import tempfile
import uuid
from urllib.parse import urlsplit

import yaml


def package_document(document, name, selected_prefix=None):
    if not isinstance(document, dict):
        raise ValueError('The package manager exported a different environment')
    reported_name = document.get('name', name)
    if selected_prefix is not None:
        # Micromamba exports a valid unregistered explicit prefix with name:null.
        # Bind its reported prefix before assigning a portable recreation name.
        reported_prefix = document.get('prefix')
        if not isinstance(reported_prefix, str) or Path(reported_prefix).resolve() != Path(selected_prefix).resolve():
            raise ValueError('The package manager exported a different environment prefix')
    if reported_name != name and not (reported_name is None and selected_prefix is not None):
        raise ValueError('The package manager exported a different environment')
    dependencies = document.get('dependencies')
    if not isinstance(dependencies, list) or not dependencies:
        raise ValueError('The installed environment export has no dependencies')
    for entry in dependencies:
        if isinstance(entry, str) and entry.strip():
            continue
        if (isinstance(entry, dict) and set(entry) == {'pip'}
                and isinstance(entry['pip'], list)
                and all(isinstance(item, str) and item.strip() for item in entry['pip'])):
            continue
        raise ValueError('The installed environment export has an unsupported dependency')
    channels = document.get('channels', [])
    if not isinstance(channels, list) or not all(isinstance(item, str) for item in channels):
        raise ValueError('The installed environment export has invalid channels')
    for channel in channels:
        parsed = urlsplit(channel)
        if (parsed.username or parsed.password or parsed.query or parsed.fragment
                or parsed.scheme == 'file' or channel.startswith('/')):
            raise ValueError('A channel cannot be exported as a portable public package source')
    # 2026-10-25 CT, user: "output the env yaml file only for installed packages".
    # Omit host-specific prefix, environment variables and solver input history.
    return dict(name=name, channels=channels, dependencies=dependencies)


def _stat_identity(value):
    return (value.st_dev, value.st_ino, value.st_mode, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns, value.st_nlink)


def _export_git(directory, *args):
    env = dict(os.environ, GIT_OPTIONAL_LOCKS='0', GIT_TERMINAL_PROMPT='0')
    for key in ('GIT_DIR', 'GIT_COMMON_DIR', 'GIT_WORK_TREE', 'GIT_INDEX_FILE'):
        env.pop(key, None)
    result = subprocess.run(['git', '-C', str(directory), '--literal-pathspecs', *args],
                            capture_output=True, env=env, timeout=30)
    if result.returncode:
        raise ValueError('Installed-package link has no committed git-annex provenance')
    return result.stdout


def _managed_export_snapshot(output):
    """Read a committed local annex link; never unlock or publish its target."""
    try:
        source = output.parent.resolve() / output.name
        before = _stat_identity(source.lstat())
        if not stat.S_ISLNK(before[2]):
            raise ValueError('Export link changed')
        literal = os.readlink(source)
        repo = Path(os.fsdecode(_export_git(source.parent, 'rev-parse', '--show-toplevel')).strip()).resolve()
        common = Path(os.fsdecode(_export_git(repo, 'rev-parse', '--git-common-dir')).strip())
        common = (common if common.is_absolute() else repo / common).resolve()
        annex_uuid = os.fsdecode(_export_git(repo, 'config', '--local', '--get', 'annex.uuid')).strip()
        uuid.UUID(annex_uuid)
        objects = common / 'annex/objects'
        target = Path(os.path.abspath(source.parent / literal))
        if (objects.is_symlink() or not objects.is_dir() or not target.is_relative_to(objects)
                or target.resolve(strict=True) != target):
            raise ValueError('Foreign export target')
        key = target.name
        match = re.fullmatch(r'(MD5E|SHA256E)-s([0-9]+)--([0-9a-f]+)(?:[.][^/]*)?', key)
        if not match or target.parent.name != key:
            raise ValueError('Unsupported annex content key')
        relative = str(source.relative_to(repo))
        entry = _export_git(repo, 'ls-tree', '-z', 'HEAD', '--', relative)
        meta, tracked = entry.rstrip(b'\0').split(b'\t')
        mode, kind, oid = meta.split()
        if (mode != b'120000' or kind != b'blob' or os.fsdecode(tracked) != relative
                or _export_git(repo, 'cat-file', 'blob', oid.decode('ascii')) != os.fsencode(literal)):
            raise ValueError('Export pointer is not the committed literal link')
        target_before = _stat_identity(target.lstat())
        if not stat.S_ISREG(target_before[2]) or target_before[3] > 1024**2:
            raise ValueError('Invalid installed-package object')
        with os.fdopen(os.open(target, os.O_RDONLY | os.O_NOFOLLOW), 'rb') as stream:
            if _stat_identity(os.fstat(stream.fileno())) != target_before:
                raise ValueError('Export object changed before read')
            payload = stream.read(1024**2 + 1)
            if _stat_identity(os.fstat(stream.fileno())) != target_before:
                raise ValueError('Export object changed during read')
        digest = (hashlib.md5(payload, usedforsecurity=False).hexdigest() if match[1] == 'MD5E'
                  else hashlib.sha256(payload).hexdigest())
        current_entry = _export_git(repo, 'ls-tree', '-z', 'HEAD', '--', relative)
        current_owner = os.fsdecode(_export_git(repo, 'config', '--local', '--get', 'annex.uuid')).strip()
        if (len(payload) != int(match[2]) or digest != match[3]
                or current_entry != entry or current_owner != annex_uuid
                or _stat_identity(source.lstat()) != before or os.readlink(source) != literal
                or _stat_identity(target.lstat()) != target_before):
            raise ValueError('Export link, key, object or committed pointer changed')
        # Other dataset metadata may be committed during this read-only export.
        # Bind this path's committed pointer and owner, rather than all of HEAD.
        return (str(source), before, literal, str(target), target_before,
                str(repo), str(common), entry, annex_uuid, payload)
    except (OSError, ValueError, UnicodeError, subprocess.SubprocessError) as error:
        raise ValueError('Refusing a linked installed-package export destination without stable committed git-annex content') from error


def export_environment(*, config_file, manager, name='priva_acmg', prefix=None):
    if not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*', name):
        raise ValueError('The environment name must be a simple name')
    config = yaml.safe_load(Path(config_file).read_text())
    selected = config.get('cache_root') if isinstance(config, dict) else None
    if not isinstance(selected, str) or not Path(selected).is_absolute():
        raise ValueError('Set cache_root to the existing absolute cache directory')
    root = Path(selected)
    if not root.is_dir() or root.resolve() == Path('/'):
        raise ValueError('cache_root must be an existing dedicated directory')
    output = root / (name + '.installed.yml')
    linked_before = _managed_export_snapshot(output) if output.is_symlink() else None
    selection = ['--name', name]
    if prefix is not None:
        prefix = Path(prefix)
        if not prefix.is_absolute() or not (prefix / 'conda-meta/history').is_file():
            raise ValueError('Select a completed absolute Mamba environment prefix')
        selection = ['--prefix', str(prefix)]
    completed = subprocess.run([str(manager), 'env', 'export', *selection],
                               capture_output=True, text=True)
    if completed.returncode:
        raise ValueError('The package manager could not export the installed environment')
    document = package_document(yaml.safe_load(completed.stdout), name, selected_prefix=prefix)
    payload = yaml.safe_dump(document, sort_keys=False)
    # 2026-10-26 CT: DataLad can record these small exports as annex links.
    # Reuse only the same committed/key-verified, byte-identical snapshot;
    # preserve every linked target, including changed or foreign destinations.
    if linked_before is not None:
        linked_after = _managed_export_snapshot(output)
        previous = yaml.safe_load(linked_after[-1])
        if (linked_after != linked_before or not isinstance(previous, dict)
                or set(previous) != {'name', 'channels', 'dependencies'} or previous.get('name') != name
                or linked_after[-1] != payload.encode('utf-8')):
            raise ValueError('Linked installed-package export changed or differs from the installed packages; retained')
        return output
    if output.is_symlink():
        raise ValueError('Refusing a linked installed-package export destination')
    if output.exists():
        previous = yaml.safe_load(output.read_text())
        if (not isinstance(previous, dict) or set(previous) != {'name', 'channels', 'dependencies'}
                or previous.get('name') != name):
            raise ValueError('Existing file is not a matching installed-package export; retained')
    if output.exists() and output.read_text() == payload:
        return output
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', prefix='.' + output.name + '.',
                                         dir=root, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return output


def cadd_environments(config):
    """Yield (name, prefix) for each completed environment of upstream CADD.sh.

    CADD.sh lets Snakemake build one environment per recipe in envs/ under
    envs/conda. Snakemake names each folder after the MD5 of that folder's real
    path and the recipe's bytes, and writes FOLDER.env_setup_done when the
    environment is complete.
    """
    base = config.get('cadd_base_dir') if isinstance(config, dict) else None
    if not isinstance(base, str) or not base or base == 'null':
        return
    recipes = Path(base) / 'envs'
    root = os.path.realpath(recipes / 'conda')
    for recipe in sorted(recipes.glob('*.yml')):
        prefix = recipes / 'conda' / (hashlib.md5(root.encode() + recipe.read_bytes()).hexdigest() + '_')
        if (prefix / 'conda-meta/history').is_file() and Path(str(prefix) + '.env_setup_done').exists():
            yield 'cadd-' + recipe.stem, prefix


def export_installation(*, config_file, manager):
    """Export PriVA and each completed CADD.sh environment without activating them."""
    config = yaml.safe_load(Path(config_file).read_text())
    exports = [export_environment(config_file=config_file, manager=manager)]
    for name, prefix in cadd_environments(config):
        exports.append(export_environment(config_file=config_file, manager=manager, name=name, prefix=prefix))
    return exports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='procedure', required=True)
    command = commands.add_parser('export')
    command.add_argument('--config', type=Path, required=True)
    command.add_argument('--manager', required=True)
    command.add_argument('--name', default='priva_acmg')
    command.add_argument('--prefix', type=Path, help='Export an installed prefix without activating it')
    installation = commands.add_parser('export-installation', help='Export PriVA and the completed CADD.sh environments')
    installation.add_argument('--config', type=Path, required=True)
    installation.add_argument('--manager', required=True)
    args = parser.parse_args()
    try:
        if args.procedure == 'export-installation':
            for path in export_installation(config_file=args.config, manager=args.manager):
                print(path)
        else:
            print(export_environment(config_file=args.config, manager=args.manager, name=args.name, prefix=args.prefix))
    except (OSError, ValueError, yaml.YAMLError) as error:
        parser.exit(1, 'ERROR: ' + str(error) + '\n')


if __name__ == '__main__':
    main()
