#!/usr/bin/env python3
"""Compare the standalone MD5 probe with Python hashlib, retaining every result.

Compile tests/support/md5_probe.rs with rustc --edition=2024 -O. Run:
  python md5_compare.py compare --binary PROBE --output RESULTS [--long-streams]
Fixtures are deterministic synthetic bytes generated in RESULTS, never patient
data. --long-streams hashes virtual messages across 512 MiB and 4 GiB boundaries
without writing or retaining those messages. No Python MD5 implementation is used.
"""
from __future__ import annotations

import argparse
import errno
import hashlib
import io
import json
import math
from pathlib import Path
import random
import re
import subprocess
import sys
import time


RFC_VECTORS = [
    (b"", "d41d8cd98f00b204e9800998ecf8427e"),
    (b"a", "0cc175b9c0f1b6a831c399e269772661"),
    (b"abc", "900150983cd24fb0d6963f7d28e17f72"),
    (b"message digest", "f96b697d7cb7938d525a2f31aaf161d0"),
    (b"abcdefghijklmnopqrstuvwxyz", "c3fcd3d76192e4007dfb496cca67e13b"),
    (b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789", "d174ab98d277d9f5a5611c2c9f419d9f"),
    (b"1234567890" * 8, "57edf4a22be3c955ac49da2e2107b67a"),
]


def md5(data=b""):
    return hashlib.md5(data, usedforsecurity=False)


def identity(path: Path) -> dict:
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest}


def messages(seed: int):
    for index, (data, expected) in enumerate(RFC_VECTORS):
        assert md5(data).hexdigest() == expected
        yield f"rfc_{index}", data
    yield "million_a", b"a" * 1_000_000
    pattern = bytes(range(256))
    lengths = set(range(261)) | {447, 448, 449, 511, 512, 513, 1023, 1024, 1025,
                                      4095, 4096, 4097, (1 << 20)-1, 1 << 20, (1 << 20)+1,
                                      (2 << 20)-1, 2 << 20, (2 << 20)+1}
    for length in sorted(lengths):
        yield f"binary_{length}", (pattern * (length // 256 + 1))[:length]
    generator = random.Random(seed)
    for index in range(128):
        yield f"random_{index}", generator.randbytes(generator.randrange(0, 65537))
    yield "space name", b"\x00\xff\x80binary\x00\n\r\t" * 31


def repeated_digest(pattern: bytes, length: int) -> str:
    # The block length is an exact multiple of the pattern, so each starts at 0.
    block = pattern * max(1, (1 << 20) // len(pattern))
    result = md5()
    while length:
        count = min(length, len(block))
        result.update(block[:count])
        length -= count
    return result.hexdigest()


class RaisingReader(io.RawIOBase):
    """Independent Python stream error, consumed by hashlib.file_digest."""
    def __init__(self, data: bytes, fail_at: int):
        self.data, self.fail_at, self.position = data, fail_at, 0

    def readable(self):
        return True

    def readinto(self, buffer):
        if self.position >= self.fail_at:
            raise OSError(errno.EIO, "injected read failure")
        count = min(len(buffer), self.fail_at - self.position, len(self.data) - self.position)
        buffer[:count] = self.data[self.position:self.position + count]
        self.position += count
        return count


def compare(args: argparse.Namespace) -> None:
    started = time.monotonic()
    args.output.mkdir(exist_ok=False)
    fixtures = args.output / "fixtures"
    fixtures.mkdir()
    pins = {"binary": identity(args.binary), "module": identity(args.module), "comparator": identity(Path(__file__))}
    source = args.module.read_text()
    constants = [int(value, 16) for value in re.findall(r"0x[0-9a-fA-F]+", source.split("const CONSTANTS:", 1)[1].split("];", 1)[0])]
    derived = [int((1 << 32) * abs(math.sin(index))) for index in range(1, 65)]
    counts = {"cases": 0, "digest_comparisons": 0, "io_failure_cases": 0, "interrupted_successes": 0, "mismatches": 0}
    records = []
    process = None
    try:
        with (args.output / "native.stderr.log").open("x") as stderr, (args.output / "results.jsonl").open("x") as output, (args.output / "mismatches.jsonl").open("x") as ledger:
            process = subprocess.Popen([str(args.binary), "serve"], text=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=stderr)

            def check(name, arguments, expected, *, fail=False, error_kind=None, consumed=None):
                counts["cases"] += 1
                if any("\t" in token or "\n" in token for token in arguments):
                    raise ValueError("probe server requests cannot contain tabs/newlines")
                process.stdin.write("\t".join(arguments) + "\n"); process.stdin.flush()
                raw = process.stdout.readline()
                if not raw:
                    raise RuntimeError("probe ended without a result; see native.stderr.log")
                actual = json.loads(raw)
                problems = []
                if fail:
                    counts["io_failure_cases"] += 1
                    if actual.get("ok") is not False or actual.get("error_kind") != error_kind:
                        problems.append("missing/wrong I/O error")
                    if any(key in actual for key in ("bytes", "stream", "reader", "file")) or re.search(r"[0-9a-f]{32}", raw):
                        problems.append("partial digest emitted on I/O failure")
                else:
                    if actual.get("ok") is not True:
                        problems.append("unexpected failure")
                    for key in expected:
                        counts["digest_comparisons"] += 1
                        if actual.get(key) != expected[key]:
                            problems.append("digest differs: " + key)
                if consumed is not None and actual.get("consumed") != consumed:
                    problems.append("stream consumed an unexpected number of bytes")
                row = {"name": name, "argv": arguments, "expected": expected, "actual": actual, "problems": problems}
                json.dump(row, output); output.write("\n")
                if problems:
                    counts["mismatches"] += 1
                    json.dump(row, ledger); ledger.write("\n")

            if constants != derived:
                row = {"name": "RFC1321_sine_constants", "expected": derived, "actual": constants}
                json.dump(row, ledger); ledger.write("\n"); counts["mismatches"] += 1
            for number, (name, data) in enumerate(messages(args.seed)):
                path = fixtures / f"{number}_{name}.bin"
                with path.open("xb") as handle:
                    handle.write(data)
                expected = md5(data).hexdigest()
                with path.open("rb") as handle:
                    assert hashlib.file_digest(handle, md5).hexdigest() == expected
                records.append(identity(path))
                check(name + "_file", ["file", "--input", str(path)], {"file": expected})
                chunks = [1, 2, 3, 7, 55, 56, 57, 63, 64, 65, 127, 128, 129, 255, 256, 257] if len(data) <= 260 else [1, 63, 64, 65, 1048576]
                for chunk in chunks:
                    check(name + f"_chunks_{chunk}", ["bytes", "--input", str(path), "--chunks", str(chunk), "--read-chunk", str(chunk)],
                          {"bytes": expected, "stream": expected, "reader": expected})
                check(name + "_mixed", ["bytes", "--input", str(path), "--chunks", "1,7,55,56,63,64,65,1025", "--read-chunk", "127"],
                      {"bytes": expected, "stream": expected, "reader": expected})
            for path, kind in ((fixtures / "missing.bin", "NotFound"), (fixtures, "IsADirectory")):
                try:
                    with path.open("rb") as handle:
                        hashlib.file_digest(handle, md5)
                except (FileNotFoundError, IsADirectoryError):
                    pass
                else:
                    raise RuntimeError("the independent filesystem failure did not fail")
                check("file_failure_" + kind, ["file", "--input", str(path)], {}, fail=True, error_kind=kind)
            pattern = bytes(range(256))
            short = pattern * 4
            for at in (0, 1, 55, 56, 63, 64, 65, 127, 128, 1023, 1024):
                python_stream = RaisingReader(short, at)
                try:
                    hashlib.file_digest(python_stream, md5)
                except OSError:
                    assert python_stream.position == at
                else:
                    raise RuntimeError("Python stream unexpectedly returned a partial digest")
                for chunk in (1, 7, 64):
                    check(f"fault_{at}_{chunk}", ["fault", "--pattern", pattern.hex(), "--length", "1024", "--fail-at", str(at), "--read-chunk", str(chunk)],
                          {}, fail=True, error_kind="InvalidData", consumed=at)
            for at in (0, 1, 55, 56, 63, 64, 65, 127, 1024):
                check(f"interrupt_{at}", ["repeat", "--pattern", pattern.hex(), "--length", "1024", "--interrupt-at", str(at), "--read-chunk", "65"],
                      {"reader": md5(short).hexdigest()}, consumed=1024)
                counts["interrupted_successes"] += 1
            check("interrupt_then_failure", ["fault", "--pattern", pattern.hex(), "--length", "1024", "--interrupt-at", "63", "--fail-at", "129", "--read-chunk", "64"],
                  {}, fail=True, error_kind="InvalidData", consumed=129)
            virtual = []
            if args.long_streams:
                for length in ((1 << 29)-1, (1 << 29)+65, (1 << 32)+65):
                    expected = repeated_digest(pattern, length)
                    check(f"virtual_{length}", ["repeat", "--pattern", pattern.hex(), "--length", str(length), "--read-chunk", "1048576"], {"reader": expected}, consumed=length)
                    virtual.append({"bytes": length, "md5": expected})
            process.stdin.close()
            returncode = process.wait(timeout=10)
            process = None
            if returncode:
                raise RuntimeError("probe returned a nonzero status")
    finally:
        if process is not None:
            process.kill(); process.wait()
    unchanged = all(identity(Path(row["path"])) == row for row in pins.values()) and all(identity(Path(row["path"])) == row for row in records)
    result = {**counts, "identities": pins, "fixtures": records, "virtual_streams": virtual,
              "constants_match_rfc1321_formula": constants == derived, "rfc_vectors": len(RFC_VECTORS),
              "inputs_unchanged": unchanged, "python": sys.version, "seed": args.seed,
              "wall_seconds": time.monotonic() - started,
              "results": identity(args.output / "results.jsonl"), "mismatch_ledger": identity(args.output / "mismatches.jsonl")}
    result["passed"] = counts["mismatches"] == 0 and unchanged and constants == derived
    with (args.output / "comparison.json").open("x") as handle:
        json.dump(result, handle, sort_keys=True, indent=2); handle.write("\n")
    print(json.dumps({key: value for key, value in result.items() if key not in ("identities", "fixtures", "python")}), flush=True)
    if not result["passed"]:
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("compare")
    command.add_argument("--binary", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    command.add_argument("--module", type=Path, default=Path(__file__).resolve().parents[1] / "src/md5.rs")
    command.add_argument("--long-streams", action="store_true")
    command.add_argument("--seed", type=lambda value: int(value, 0), default=0x4D44355F5354524D)
    command.set_defaults(function=compare)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
