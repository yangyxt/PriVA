#!/usr/bin/env python3
"""A stand-in for curl and wget that serves files from a folder.

Installed as FOLDER/bin/curl and FOLDER/bin/wget (``fake_net.py install
FOLDER``), first on PATH, it answers URL with the file FAKE_NET_ROOT/HOST/PATH
(the query string left out), so an install step's old shell function and its
process can be run on the same bytes without the network. A missing file is
the server's 404: curl -f exits 22, wget exits 8. Every request is appended to
FAKE_NET_ROOT/requests.log.
"""
import os
import shutil
import sys
from pathlib import Path
from urllib.parse import urlsplit

CURL_VALUES = {"-o", "--output", "--retry", "--retry-delay", "--connect-timeout",
               "--max-time", "-m", "--max-filesize", "-C", "--continue-at", "-H",
               "--header", "-A", "--user-agent", "-w", "--write-out", "-D", "--dump-header"}
WGET_VALUES = {"-O", "--output-document", "-t", "--tries", "-T", "--timeout", "-P"}


def served(url: str) -> Path:
    parts = urlsplit(url)
    return Path(os.environ["FAKE_NET_ROOT"]) / parts.netloc / parts.path.lstrip("/")


def answer(url: str, output: str | None) -> bool:
    root = Path(os.environ["FAKE_NET_ROOT"])
    with open(root / "requests.log", "a") as log:
        log.write(url + "\n")
    source = served(url)
    if not source.is_file():
        print(f"fake_net: 404 {url}", file=sys.stderr)
        return False
    if output in (None, "-"):
        sys.stdout.buffer.write(source.read_bytes())
    else:
        shutil.copyfile(source, output)
    return True


def curl(args: list[str]) -> int:
    output, urls, at = None, [], 0
    while at < len(args):
        arg = args[at]
        if arg in CURL_VALUES:
            if arg in ("-o", "--output"):
                output = args[at + 1]
            at += 2
            continue
        if arg.startswith("-"):
            at += 1
            continue
        urls.append(arg)
        at += 1
    ok = all(answer(url, output) for url in urls)
    return 0 if ok else 22


def wget(args: list[str]) -> int:
    output, urls, at = None, [], 0
    while at < len(args):
        arg = args[at]
        if arg in WGET_VALUES:
            if arg in ("-O", "--output-document"):
                output = args[at + 1]
            at += 2
            continue
        if arg.startswith("-O") and len(arg) > 2:
            output = arg[2:]
        elif not arg.startswith("-"):
            urls.append(arg)
        at += 1
    for url in urls:
        if not answer(url, output or os.path.basename(urlsplit(url).path)):
            if output and os.path.exists(output) and os.path.getsize(output) == 0:
                os.remove(output)
            return 8
    return 0


def install(folder: Path) -> None:
    bin_dir = folder / "bin"
    bin_dir.mkdir(parents=True, exist_ok=True)
    for name in ("curl", "wget"):
        target = bin_dir / name
        target.write_text(f"#!/bin/sh\nexec {sys.executable} {Path(__file__).resolve()} {name} \"$@\"\n")
        target.chmod(0o755)


if __name__ == "__main__":
    tool, rest = sys.argv[1], sys.argv[2:]
    if tool == "install":
        install(Path(rest[0]))
    elif tool == "curl":
        sys.exit(curl(rest))
    elif tool == "wget":
        sys.exit(wget(rest))
    else:
        raise SystemExit(f"fake_net: unknown tool {tool}")
