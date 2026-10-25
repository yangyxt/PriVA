"""Relay only a labelled fixture's scientific command across the new subprocess.

The real publication helper and receipts run unchanged. The exported Bash
fixture_python function retains each existing test's own producer assertions.
This driver is never used by an installer.
"""
import argparse
import importlib.util
from pathlib import Path
import subprocess
import sys
from unittest import mock


def run(args):
    spec = importlib.util.spec_from_file_location("fixture_mapped_publication", args.publisher)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    native_run = subprocess.run

    def scientific(argv, **kwargs):
        return native_run(["bash", "--noprofile", "--norc", "-c", 'fixture_python "$@"',
                           "fixture", *argv[1:]], **kwargs)

    arguments = args.arguments[1:] if args.arguments[:1] == ["--"] else args.arguments
    # This bridge covers mapping/publication, with an explicitly labelled base
    # prerequisite stub. The new base-binding suite exercises the real option.
    if "--annotation-base" in arguments:
        position = arguments.index("--annotation-base")
        arguments = arguments[:position] + arguments[position + 2:]
    sys.argv = [str(args.publisher), *arguments]
    with mock.patch.object(module.subprocess, "run", side_effect=scientific):
        return module.main()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    child = commands.add_parser("run")
    child.add_argument("--publisher", type=Path, required=True)
    child.add_argument("arguments", nargs=argparse.REMAINDER)
    child.set_defaults(function=run)
    args = parser.parse_args()
    return args.function(args)


if __name__ == "__main__":
    raise SystemExit(main())
