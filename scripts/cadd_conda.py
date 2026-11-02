#!/usr/bin/env python3
"""Find the real Conda frontend required by upstream CADD's Snakemake.

check is read-only. install may create only the dedicated frontend prefix;
it never installs into priva_acmg or CADD's scientific rule environments.
Use --print-path to obtain the child process PATH without shell evaluation.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from packaging.version import Version

MINIMUM_VERSION = Version("24.7.1")  # Snakemake 8.25.2 deployment/conda.py.
FRONTEND_SPEC = "conda=24.9.2"


def command_json(command):
    result = subprocess.run(command, check=True, capture_output=True, text=True, timeout=60)
    value = json.loads(result.stdout)
    if not isinstance(value, dict):
        raise ValueError("expected an object from " + " ".join(command))
    return value


def inspect_conda(executable):
    """Check the same interfaces Snakemake uses, not a particular installation origin."""
    executable = Path(executable).absolute()
    version = subprocess.run([str(executable), "--version"], check=True,
                             capture_output=True, text=True, timeout=60).stdout.strip()
    match = re.fullmatch(r"conda (\S+)", version)
    if not match or Version(match.group(1)) < MINIMUM_VERSION:
        raise ValueError(f"{executable}: CADD requires Conda >= 24.7.1; found {version!r}")
    info = command_json([str(executable), "info", "--json"])
    if not isinstance(info.get("platform"), str) or not info["platform"]:
        raise ValueError(f"{executable}: Conda info lacks its platform")
    prefix = Path(info["conda_prefix"])
    if not prefix.is_absolute() or not (prefix / "bin/activate").is_file():
        raise ValueError(f"{executable}: reported Conda prefix lacks bin/activate")
    config = command_json([str(executable), "config", "--get", "channel_priority", "--json"])
    if config.get("success") is False or not isinstance(config.get("get"), dict):
        raise ValueError(f"{executable}: cannot read Conda channel priority")
    return {"executable": str(executable), "prefix": str(prefix), "version": version}


def same_file(left, right):
    return bool(left and right and Path(left).samefile(right))


def child_path(selected):
    """Expose Conda while retaining the caller's Python and other tool precedence."""
    original = os.environ.get("PATH", "")
    executable = selected["executable"]
    if same_file(shutil.which("conda", path=original), executable):
        return original
    # Appending bin preserves scientific tools, including Python. A standard
    # condabin also permits replacing an incompatible earlier Conda without
    # putting the frontend's separate Python ahead of priva_acmg's Python.
    candidate = str(Path(executable).parent)
    path = original + os.pathsep + candidate
    if not same_file(shutil.which("conda", path=path), executable):
        condabin = Path(selected["prefix"]) / "condabin"
        alternate = inspect_conda(condabin / "conda")
        if alternate["prefix"] != selected["prefix"] or (condabin / "python").exists():
            raise ValueError("cannot expose the selected Conda without changing scientific Python")
        selected.update(alternate)
        path = str(condabin) + os.pathsep + original
    if shutil.which("python", path=path) != shutil.which("python", path=original):
        raise ValueError("the Conda frontend must not replace the active scientific Python")
    return path


def discover(frontend):
    candidates = [shutil.which("conda"), os.environ.get("CONDA_EXE"), str(frontend / "bin/conda")]
    failures, seen = [], set()
    for candidate in candidates:
        if not candidate or candidate in seen:
            continue
        seen.add(candidate)
        if not Path(candidate).exists():
            continue
        try:
            selected = inspect_conda(candidate)
            selected["path"] = child_path(selected)
            return selected
        except (OSError, ValueError, KeyError, subprocess.SubprocessError) as exc:
            failures.append(f"{candidate}: {exc}")
    if failures:
        print("Unusable Conda candidates: " + "; ".join(failures), file=sys.stderr)
    return None


def ensure(frontend, manager=None):
    selected = discover(frontend)
    if selected:
        return selected
    if manager is None:
        raise ValueError("CADD.sh needs a conda executable (>= 24.7.1); run main_install or CADD_conda_install")
    # A reused scientific environment may have exact package/export receipts.
    # Installing a manager into it would invalidate those receipts even if no
    # scientific package changed. Reuse the dedicated CADD frontend location.
    frontend.parent.mkdir(parents=True, exist_ok=True)
    with (frontend.parent / ".conda_frontend.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        selected = discover(frontend)
        if selected:
            return selected
        if frontend.exists() or frontend.is_symlink():
            raise ValueError(f"incomplete or incompatible frontend at {frontend}; retained without modification")
        command = [str(manager), "create", "--yes", "--prefix", str(frontend),
                   "--override-channels", "--channel", "conda-forge", FRONTEND_SPEC]
        print("Creating CADD's tool-only Conda frontend: " + " ".join(command), file=sys.stderr)
        # Preserve the selected manager/wrapper and normal package verification.
        # Never rename a Conda prefix: installed entry points contain its path.
        subprocess.run(command, check=True, stdout=sys.stderr)
        selected = inspect_conda(frontend / "bin/conda")
        selected["path"] = child_path(selected)
        return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="procedure", required=True)
    for name in ("check", "install"):
        command = commands.add_parser(name)
        command.add_argument("--frontend-prefix", type=Path, required=True)
        command.add_argument("--print-path", action="store_true")
        if name == "install":
            command.add_argument("--manager", type=Path, required=True)
    args = parser.parse_args()
    try:
        if not args.frontend_prefix.is_absolute() or args.frontend_prefix.name != "conda_frontend":
            raise ValueError("use the absolute cadd_parent_dir/conda_frontend prefix")
        selected = ensure(args.frontend_prefix, getattr(args, "manager", None))
        print(f"Using {selected['version']}: {selected['executable']}", file=sys.stderr)
        print(selected["path"] if args.print_path else json.dumps(selected, sort_keys=True))
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
