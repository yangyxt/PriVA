#!/usr/bin/env python3
"""Export the packages actually installed in PriVA and CADD environments.

export --config CONFIG --manager MAMBA_OR_CONDA [--name priva_acmg] [--prefix PREFIX]
The portable YAML is written to cache_root/<name>.installed.yml. It describes
the selected environment only; it makes no claim about other CADD environments.
export-installation exports PriVA and every completed native CADD environment.
"""
import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
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
    output = root / (name + '.installed.yml')
    if output.is_symlink():
        raise ValueError('Refusing a linked installed-package export destination')
    if output.exists():
        previous = yaml.safe_load(output.read_text())
        if (not isinstance(previous, dict) or set(previous) != {'name', 'channels', 'dependencies'}
                or previous.get('name') != name):
            raise ValueError('Existing file is not a matching installed-package export; retained')
    payload = yaml.safe_dump(document, sort_keys=False)
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


def export_installation(*, config_file, manager):
    """Export PriVA and each completed CADD environment without activating them."""
    helper = Path(__file__).with_name('cadd_runtime_install.py')
    completed = subprocess.run([sys.executable, str(helper), 'prefixes', '--config', str(config_file),
                                '--package-manager', str(manager)], capture_output=True, text=True)
    if completed.returncode:
        raise ValueError('Could not resolve the installed CADD environment prefixes')
    import json
    layout = json.loads(completed.stdout)
    backend = layout.get('backend', 'native')
    if backend not in ('native', 'container'):
        raise ValueError('Unknown installed CADD environment backend')
    environments = layout.get('environments', {})
    core = environments.get('cadd-python2', {})
    if backend == 'native' and (not core.get('prefix') or not Path(core['prefix'], 'conda-meta/history').is_file()):
        raise ValueError('The separate CADD scorer environment is not installed')
    exports = [export_environment(config_file=config_file, manager=manager)]
    for name, record in environments.items():
        prefix = Path(record['prefix'])
        if (prefix / 'conda-meta/history').is_file():
            exports.append(export_environment(config_file=config_file, manager=manager,
                                               name=name, prefix=prefix))
    return exports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='procedure', required=True)
    command = commands.add_parser('export')
    command.add_argument('--config', type=Path, required=True)
    command.add_argument('--manager', required=True)
    command.add_argument('--name', default='priva_acmg')
    command.add_argument('--prefix', type=Path, help='Export an installed prefix without activating it')
    installation = commands.add_parser('export-installation', help='Export PriVA and completed separate CADD environments')
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
