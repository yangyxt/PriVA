#!/usr/bin/env python3
"""A stand-in for curl and wget that serves files from a folder.

Installed as FOLDER/bin/curl and FOLDER/bin/wget (``fake_net.py install
FOLDER``), first on PATH, it answers URL with the file FAKE_NET_ROOT/HOST/PATH
(the query string left out; a PATH that ends in / is PATH/_index), so an
install step's old shell function and its process can be run on the same
bytes without the network. A missing file is
the server's 404: curl -f exits 22, wget exits 8. curl --include (-i) writes
the reply's status line and headers before the body (Content-Length,
Last-Modified from the file's time, an ETag from its bytes; a Range header is
ignored, as some servers do), and a 404 is then a reply too: curl exits 0
unless -f is given. Every request is appended to FAKE_NET_ROOT/requests.log.
"""
import email.utils
import hashlib
import os
import shutil
import sys
from pathlib import Path
from urllib.parse import urlsplit

CURL_VALUES = {"-o", "--output", "--retry", "--retry-delay", "--connect-timeout",
               "--max-time", "-m", "--max-filesize", "-C", "--continue-at", "-H",
               "--header", "-A", "--user-agent", "-w", "--write-out", "-D", "--dump-header",
               "-x", "--proxy", "--max-redirs", "--speed-limit", "--speed-time"}
WGET_VALUES = {"-O", "--output-document", "-t", "--tries", "-T", "--timeout", "-P"}


def served(url: str) -> Path:
    parts = urlsplit(url)
    path = parts.path.lstrip("/")
    # A folder's own URL (a path that ends in /) is its _index file.
    if path.endswith("/") or not path:
        path += "_index"
    return Path(os.environ["FAKE_NET_ROOT"]) / parts.netloc / path


def head(source: Path) -> bytes:
    """The status line and headers of a reply with SOURCE (None: a 404)."""
    if not source.is_file():
        return b"HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\n\r\n"
    body = source.read_bytes()
    modified = email.utils.formatdate(source.stat().st_mtime, usegmt=True)
    tag = hashlib.sha256(body).hexdigest()[:16]
    return (f"HTTP/1.1 200 OK\r\nContent-Length: {len(body)}\r\nLast-Modified: {modified}\r\n"
            f"ETag: \"{tag}\"\r\n\r\n").encode()


def answer(url: str, output: str | None, include: bool = False, fail: bool = True) -> bool:
    root = Path(os.environ["FAKE_NET_ROOT"])
    with open(root / "requests.log", "a") as log:
        log.write(url + "\n")
    source = served(url)
    found = source.is_file()
    if not found:
        print(f"fake_net: 404 {url}", file=sys.stderr)
        if fail or not include:
            return False
    reply = (head(source) if include else b"") + (source.read_bytes() if found else b"")
    if output in (None, "-"):
        sys.stdout.buffer.write(reply)
        sys.stdout.buffer.flush()
    elif include:
        Path(output).write_bytes(reply)
    else:
        shutil.copyfile(source, output)
    return True


def curl(args: list[str]) -> int:
    output, urls, at, flags = None, [], 0, set()
    while at < len(args):
        arg = args[at]
        if arg in CURL_VALUES:
            if arg in ("-o", "--output"):
                output = args[at + 1]
            at += 2
            continue
        if arg.startswith("--"):
            flags.add(arg)
        elif arg.startswith("-"):
            flags.update(f"-{letter}" for letter in arg[1:])
        else:
            urls.append(arg)
        at += 1
    include = bool(flags & {"--include", "-i"})
    fail = bool(flags & {"--fail", "-f"})
    ok = all(answer(url, output, include, fail or not include) for url in urls)
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
