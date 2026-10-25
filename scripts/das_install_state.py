#!/usr/bin/env python3
"""Completion metadata for the installer's multi-output DAS build.

This records successful execution/publication, not scientific qualification or
content validation. A nonempty leftover is not evidence of a completed build.
2026-10-24 CT: source CONTENT under a logical role is independent of checkout
path/inode. Scientific inputs/outputs retain their strict metadata identities.
Known bundled resources are checked only for existence/nonzero size; their
contents are never read. Named parameters retain their supplied order.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import sys
import tempfile

VERSION = 3
IDENTITY_POLICY = ("data=path,target,inode,size,mtime_ns; node_local_st_dev_excluded; "
                   "source=role,sha256; bundled=role,exists_nonempty; parameters=ordered")


def identities(paths):
    result = []
    for name in paths:
        path = Path(name).absolute()
        info = path.stat()
        if not stat.S_ISREG(info.st_mode) or info.st_size == 0:
            raise ValueError(f"Missing/nonregular/empty DAS file: {path}")
        result.append({"path": str(path), "target": str(path.resolve()),
                       "inode": info.st_ino,
                       "size": info.st_size, "mtime_ns": info.st_mtime_ns})
    return result


def named(values):
    result = []
    for value in values:
        name, separator, selected = value.partition('=')
        if not separator or not name or not selected:
            raise ValueError('Use a nonempty logical NAME=VALUE: ' + value)
        result.append([name, selected])
    return result


def sources(values):
    result = []
    for role, selected in named(values):
        path = Path(selected)
        before = identities([path])
        before_ctime = path.stat().st_ctime_ns
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(block)
        if identities([path]) != before or path.stat().st_ctime_ns != before_ctime:
            raise ValueError('DAS source changed while hashing: ' + str(path))
        result.append({'role': role, 'sha256': digest.hexdigest()})
    return result


def bundled(values):
    result = []
    for role, selected in named(values):
        path = Path(selected)
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError('Missing/nonregular/empty bundled DAS file: ' + str(path))
        result.append({'role': role, 'exists_nonempty': True})
    return result


def input_state(args):
    return {"profile": args.profile, "identity_policy": IDENTITY_POLICY,
            "inputs": identities(args.input),
            "sources": sources(getattr(args, 'source', [])),
            "bundled": bundled(getattr(args, 'bundled', [])),
            "parameters": named(getattr(args, 'parameter', []))}


def record(args):
    current = input_state(args)
    if current != json.loads(args.expected_inputs):
        raise ValueError("DAS inputs changed during generation; no completion recorded")
    value = {"version": VERSION, **current, "outputs": identities(args.output)}
    # Publish only the small receipt atomically. The scientific writer's order,
    # filenames, formats and output bytes remain unchanged.
    path = Path(args.receipt)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(value, handle, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def check(args):
    expected = {"version": VERSION, **input_state(args), "outputs": identities(args.output)}
    with open(args.receipt) as handle:
        if json.load(handle) != expected:
            raise ValueError("DAS completion metadata no longer matches selected files")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("capture", "check", "record"):
        command = commands.add_parser(name)
        command.add_argument("--profile", required=True)
        command.add_argument("--input", action="append", required=True)
        command.add_argument("--source", action="append", default=[], metavar="NAME=PATH")
        command.add_argument("--bundled", action="append", default=[], metavar="NAME=PATH")
        command.add_argument("--parameter", action="append", default=[], metavar="NAME=VALUE")
        if name != "capture":
            command.add_argument("--receipt", required=True)
            command.add_argument("--output", action="append", required=True)
        if name == "record":
            command.add_argument("--expected-inputs", required=True)
    args = parser.parse_args()
    try:
        if args.command == "capture":
            print(json.dumps(input_state(args)))
        elif args.command == "check":
            check(args)
        else:
            record(args)
    except (OSError, ValueError) as error:
        print(f"DAS install completion: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
