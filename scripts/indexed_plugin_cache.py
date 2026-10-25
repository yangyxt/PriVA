#!/usr/bin/env python3
"""Prepare external indexed plugin tables without rewriting selected source data.

``ensure-index`` repairs only the adjacent TBI, preserving readable publisher
coordinate/header settings. ``format-primate`` streams the original PrimateAI
formatting steps into a managed, source- and assembly-bound cache directory.
Neither entry point downloads or substitutes a scientific source.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile


DEFAULT_INDEX = {
    "mavedb": (0, 1, 2, 3, ord("#"), 0),
    "splicevault": (0, 1, 2, 2, ord("#"), 0),
    "primateai": (0, 1, 2, 2, ord("#"), 0),
}


def sha256(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def stat_identity(path: Path) -> tuple:
    value = path.stat()
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)


def run(argv: list[str], **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(argv, check=True, **kwargs)


def index_settings(index: Path) -> tuple | None:
    """Read the actual generic TBI coordinate flags, columns, comment and skip."""
    try:
        with gzip.open(index, "rb") as stream:
            header = stream.read(36)
        if len(header) != 36 or header[:4] != b"TBI\x01":
            return None
        n_ref, fmt, seq, begin, end, meta, skip, names = struct.unpack("<8i", header[4:])
        if (n_ref < 1 or fmt < 0 or min(seq, begin, end) < 1
                or not 0 <= meta <= 127 or skip < 0 or names < 1):
            return None
        return fmt, seq, begin, end, meta, skip
    except (OSError, EOFError, ValueError):
        return None


def first_record(path: Path, settings: tuple) -> tuple[str, int, int, str]:
    fmt, seq, begin, end, meta, skip = settings
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for number in range(4096):
            line = stream.readline(1024 * 1024)
            if not line:
                break
            if not line.endswith("\n") and len(line) >= 1024 * 1024:
                raise ValueError("Plugin table header/record exceeds the bounded probe")
            if number < skip or not line.strip() or (meta and line.startswith(chr(meta))):
                continue
            fields = line.rstrip("\r\n").split("\t")
            try:
                chromosome = fields[seq - 1]
                start, finish = int(fields[begin - 1]), int(fields[end - 1])
            except (IndexError, ValueError) as error:
                raise ValueError("Plugin table has no valid first coordinate record") from error
            if not chromosome or start < (0 if fmt & 0x10000 else 1) or finish < start:
                raise ValueError("Plugin table coordinates disagree with the TBI coordinate settings")
            return chromosome, start + bool(fmt & 0x10000), max(finish, start + 1), line.rstrip("\r\n")
    raise ValueError("Plugin table contains no bounded first data record")


def boundary_queries(path: Path, record: tuple) -> dict[str, bytes]:
    chromosome, start, end, _ = record
    positions = sorted({max(1, start - 1), max(1, start), start + 1, max(1, end), end + 1})
    return {f"{chromosome}:{pos}-{pos}":
            run(["tabix", str(path), f"{chromosome}:{pos}-{pos}"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout for pos in positions}


def readable_queries(path: Path, record: tuple) -> dict[str, bytes] | None:
    try:
        run(["tabix", "-l", str(path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        result = boundary_queries(path, record)
        if any(record[3].encode() in value.splitlines() for value in result.values()):
            return result
    except subprocess.CalledProcessError:
        pass
    return None


def ensure_index(*, input_file: Path, kind: str) -> dict:
    """Validate BGZF and repair only the index, retaining publisher TBI semantics."""
    path = input_file.absolute()
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty selected {kind} table: {path}")
    before = stat_identity(path)
    # gzip -t accepts ordinary gzip, whereas tabix requires a real BGZF member.
    with path.open("rb") as stream:
        prefix = stream.read(18)
    if (len(prefix) < 18 or prefix[:4] != b"\x1f\x8b\x08\x04"
            or prefix[12:16] != b"BC\x02\x00"):
        raise ValueError(f"Selected {kind} table is not BGZF: {path}")
    run(["bgzip", "-t", str(path)])
    index = Path(str(path) + ".tbi")
    settings = index_settings(index) or DEFAULT_INDEX[kind]
    if settings[0] & ~0x10000:
        raise ValueError(f"Selected {kind} TBI has an unsupported non-generic format")
    if settings[1:4] != DEFAULT_INDEX[kind][1:4]:
        raise ValueError(f"Selected {kind} TBI columns disagree with the plugin's columns")
    record = first_record(path, settings)
    old_queries = readable_queries(path, record) if index.is_file() else None
    if old_queries is not None and index.stat().st_mtime_ns >= path.stat().st_mtime_ns:
        if stat_identity(path) != before:
            raise ValueError("Selected table changed during index validation")
        return {"path": str(path), "index": str(index), "repaired": False}
    # tabix has no output-index option. Index a temporary alias, then publish
    # only its independent TBI; the selected table is never opened for writing.
    with tempfile.TemporaryDirectory(prefix=f".{path.name}.index-", dir=path.parent) as temporary:
        staged = Path(temporary) / path.name
        staged.symlink_to(path)
        fmt, seq, begin, end, meta, skip = settings
        argv = ["tabix", "-s", str(seq), "-b", str(begin), "-e", str(end),
                "-c", chr(meta) if meta else "", "-S", str(skip)]
        if fmt & 0x10000:
            argv.append("-0")
        run(argv + [str(staged)])
        new_queries = readable_queries(staged, record)
        if new_queries is None:
            raise ValueError("Repaired index cannot retrieve the first data record")
        if old_queries is not None and new_queries != old_queries:
            raise ValueError("Repaired index changed bounded coordinate-boundary query results")
        if index_settings(Path(str(staged) + ".tbi")) != settings:
            raise ValueError("Repaired TBI did not preserve publisher coordinate/header settings")
        if stat_identity(path) != before:
            raise ValueError("Selected table changed during index repair")
        os.replace(str(staged) + ".tbi", index)
    return {"path": str(path), "index": str(index), "repaired": True}


def primate_format(path: Path) -> str:
    with gzip.open(path, "rb") as stream:
        lines = [stream.readline(1024 * 1024) for _ in range(12)]
    header = lines[11].split(b"\t", 1)[0].lower().strip()
    if header in (b"chr", b"chrom", b"chromosome"):
        return "raw"
    first_record(path, DEFAULT_INDEX["primateai"])
    return "prepared"


@contextmanager
def locked(path: Path):
    with path.open("a") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        yield


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def stream_primate(*, source: Path, output: Path, threads: int, tmp_dir: Path) -> None:
    """Exact gunzip/sed/sed/awk+sort rules from the existing upstream command.

    Physical line 12 is prefixed before blank removal. The first eleven
    remaining records bypass sort; all others use the original GNU sort keys
    and full-row tie ordering. sort/bgzip failures propagate independently.
    """
    body = output.parent / "sorted.body"
    header = bytearray()
    with body.open("wb") as sorted_stream:
        sorter = subprocess.Popen(["sort", "-k1,1", "-k", "2,2n", "-V",
                                   f"--parallel={threads}", "--buffer-size=256M", "-T", str(tmp_dir)],
                                  stdin=subprocess.PIPE, stdout=sorted_stream)
        try:
            count = 0
            with gzip.open(source, "rb") as stream:
                for number, line in enumerate(stream, 1):
                    line = line.removesuffix(b"\n")
                    if number == 12:
                        line = b"#" + line
                    if not line:
                        continue
                    count += 1
                    if count < 12:
                        header.extend(line + b"\n")
                    else:
                        sorter.stdin.write(line + b"\n")
            sorter.stdin.close()
            if sorter.wait() != 0:
                raise ValueError("PrimateAI sort failed")
        finally:
            if sorter.poll() is None:
                sorter.kill()
                sorter.wait()
    with output.open("wb") as compressed:
        compressor = subprocess.Popen(["bgzip", "-@", str(threads), "-c"],
                                      stdin=subprocess.PIPE, stdout=compressed)
        try:
            compressor.stdin.write(header)
            with body.open("rb") as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    compressor.stdin.write(block)
            compressor.stdin.close()
            if compressor.wait() != 0:
                raise ValueError("PrimateAI BGZF compression failed")
        finally:
            if compressor.poll() is None:
                compressor.kill()
                compressor.wait()
    body.unlink()


def format_primate(*, input_file: Path, cache_dir: Path, assembly: str,
                   threads: int = 1, tmp_dir: Path | None = None) -> dict:
    if assembly not in ("GRCh37", "hg19", "GRCh38", "hg38"):
        raise ValueError("An explicit supported PrimateAI source assembly is required")
    source = input_file.absolute()
    if not source.is_file() or not source.stat().st_size:
        raise ValueError(f"Missing or empty selected PrimateAI source: {source}")
    if primate_format(source) == "prepared":
        return ensure_index(input_file=source, kind="primateai")
    before = stat_identity(source)
    identity = {"path": str(source), "sha256": sha256(source), "size": before[2], "assembly": assembly}
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:20]
    parent = cache_dir.absolute() / assembly
    parent.mkdir(parents=True, exist_ok=True)
    destination = parent / key
    name = source.name.removesuffix(".gz").removesuffix(".tsv") + ".sorted.tsv.bgz"
    with locked(parent / f".{key}.lock"):
        if destination.exists():
            receipt = json.loads((destination / "source.ready.json").read_text())
            output = destination / name
            if (receipt.get("source") != identity or receipt.get("algorithm") != "upstream-primate-format-v1"
                    or receipt.get("output_sha256") != sha256(output)):
                raise ValueError(f"Conflicting managed PrimateAI cache: {destination}")
            result = ensure_index(input_file=output, kind="primateai")
            if receipt.get("index_sha256") != sha256(Path(str(output) + ".tbi")):
                raise ValueError(f"Conflicting managed PrimateAI index: {output}.tbi")
            if stat_identity(source) != before:
                raise ValueError("Selected PrimateAI source changed during cache reuse")
            return result
        with tempfile.TemporaryDirectory(prefix=f".{key}.build-", dir=parent) as temporary:
            stage = Path(temporary)
            output = stage / name
            temporary_sort = tmp_dir or stage
            temporary_sort.mkdir(parents=True, exist_ok=True)
            stream_primate(source=source, output=output, threads=threads, tmp_dir=temporary_sort)
            ensure_index(input_file=output, kind="primateai")
            if stat_identity(source) != before:
                raise ValueError("Selected PrimateAI source changed during formatting")
            write_json(stage / "source.ready.json", {
                "source": identity, "algorithm": "upstream-primate-format-v1",
                "output_sha256": sha256(output), "index_sha256": sha256(Path(str(output) + ".tbi")),
                "scientific_qualified": False,
            })
            os.replace(stage, destination)
    return {"path": str(destination / name), "index": str(destination / (name + ".tbi")), "formatted": True}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    index = commands.add_parser("ensure-index", help="Validate BGZF and atomically repair its TBI")
    index.add_argument("--input", type=Path, required=True)
    index.add_argument("--kind", choices=tuple(DEFAULT_INDEX), required=True)
    primate = commands.add_parser("format-primate", help="Prepare a selected PrimateAI source under managed cache")
    primate.add_argument("--input", type=Path, required=True)
    primate.add_argument("--cache-dir", type=Path, required=True)
    primate.add_argument("--assembly", required=True, choices=("hg19", "GRCh37", "hg38", "GRCh38"))
    primate.add_argument("--threads", type=int, default=1)
    primate.add_argument("--tmp-dir", type=Path)
    primate.add_argument("--print-path", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.command == "ensure-index":
            result = ensure_index(input_file=args.input, kind=args.kind)
        else:
            if args.threads < 1:
                raise ValueError("threads must be positive")
            result = format_primate(input_file=args.input, cache_dir=args.cache_dir, assembly=args.assembly,
                                    threads=args.threads, tmp_dir=args.tmp_dir)
        print(result["path"] if getattr(args, "print_path", False) else json.dumps(result, sort_keys=True))
        return 0
    except (OSError, ValueError, EOFError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
