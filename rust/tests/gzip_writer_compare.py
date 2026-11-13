#!/usr/bin/env python3
"""Check the standalone gzip writer against Python gzip/zlib and pandas.

Run `compare --binary PROBE --output DIRECTORY`. Output is newly created and
retains every synthetic input, output, result and mismatch. Compressed payload
bytes may differ between backends; decompressed bytes and selected gzip metadata
must match exactly. No serialization, pandas formatting or BGZF is implemented.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import random
import struct
import subprocess
import sys
import time
import zlib

import pandas as pd


def identity(path):
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest}


def header(data):
    if len(data) < 10 or data[:3] != b"\x1f\x8b\x08" or data[3] not in (0, 8):
        raise ValueError("unexpected gzip header")
    end = data.index(b"\0", 10) if data[3] & 8 else 9
    return {"flags": data[3], "mtime": int.from_bytes(data[4:8], "little"),
            "xfl": data[8], "os": data[9], "filename_hex": data[10:end].hex() if data[3] & 8 else ""}


def python_gzip(data, *, level=9, mtime=0, filename=b""):
    output = io.BytesIO()
    with gzip.GzipFile(fileobj=output, mode="wb", compresslevel=level, mtime=mtime, filename=filename) as writer:
        writer.write(data)
    return output.getvalue()


def messages(seed):
    yield "empty", b""
    yield "unicode", "gene\tβ\t中\né\tU0001f9ec\r\n".encode()
    yield "json", json.dumps({"gene": "β", "null": None, "rows": ["x\ny", "\0", "\\"]}, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    pattern = bytes(range(256))
    for length in (1, 2, 7, 15, 16, 31, 32, 33, 55, 56, 57, 63, 64, 65, 127, 128, 129,
                   255, 256, 257, 1023, 1024, 1025, 65535, 65536, 65537, 1048575, 1048576, 1048577):
        yield f"binary_{length}", (pattern * (length // 256 + 1))[:length]
    randomizer = random.Random(seed)
    for index in range(16):
        yield f"random_{index}", randomizer.randbytes(randomizer.randrange(1, 65537))


def compare(args):
    started = time.monotonic()
    args.output.mkdir(exist_ok=False)
    inputs, outputs = args.output / "inputs", args.output / "outputs"
    inputs.mkdir(); outputs.mkdir()
    pins = {"binary": identity(args.binary), "module": identity(args.module), "comparator": identity(Path(__file__))}
    counts = {"cases": 0, "payload_checks": 0, "header_checks": 0, "trailer_checks": 0,
              "deterministic_repeats": 0, "compressed_backend_differences": 0, "corruption_rejections": 0,
              "failure_cases": 0, "interrupted_successes": 0, "pandas_cases": 0,
              "plain_pandas_cases": 0, "python_io_failures": 0, "mismatches": 0}
    fixture_pins = []
    process = None
    with (args.output / "native.stderr.log").open("x") as stderr, (args.output / "results.jsonl").open("x") as results, (args.output / "mismatches.jsonl").open("x") as failures:
        process = subprocess.Popen([str(args.binary), "serve"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=stderr, text=True)

        def record(row):
            json.dump(row, results); results.write("\n")
            if row.get("problems"):
                counts["mismatches"] += 1
                json.dump(row, failures); failures.write("\n")

        def check(name, path, data, *, level=9, mtime=0, filename=None, chunks="65536", flush_every=0,
                  options=(), error_kind=None, expected_python=None, repeated=None, corrupt=False):
            counts["cases"] += 1
            target = outputs / f"{counts['cases']:05d}.gz"
            arguments = ["encode", "--input", str(path), "--output", str(target), "--level", str(level),
                         "--mtime", str(mtime), "--chunks", chunks, "--flush-every", str(flush_every), *options]
            if filename is not None:
                arguments += ["--name-hex", filename.hex()]
            if any("\n" in token or "\t" in token for token in arguments):
                raise ValueError("probe server paths cannot contain tabs/newlines")
            process.stdin.write("\t".join(arguments) + "\n"); process.stdin.flush()
            reply = process.stdout.readline()
            if not reply:
                raise RuntimeError("probe ended; inspect native.stderr.log")
            actual = json.loads(reply)
            compressed = target.read_bytes()
            problems = []
            if "--fault-at" in options and not actual["fault_fired"]:
                problems.append("requested byte-position fault was not exercised")
            if "interrupt" in options and actual["interruptions"] != 1:
                problems.append("requested interrupted I/O was not exercised exactly once")
            if error_kind:
                counts["failure_cases"] += 1
                if actual["ok"] or actual["finish_ok"] or actual["error_kind"] != error_kind:
                    problems.append("I/O or option failure reported success or the wrong error kind")
                if error_kind == "InvalidInput" and compressed:
                    problems.append("invalid metadata wrote bytes")
            else:
                if not actual["ok"] or not actual["finish_ok"]:
                    problems.append("unexpected native failure")
                expected_python = expected_python if expected_python is not None else python_gzip(data, level=level, mtime=mtime, filename=filename or b"")
                counts["compressed_backend_differences"] += compressed != expected_python
                try:
                    counts["header_checks"] += 1
                    if header(compressed) != header(expected_python):
                        problems.append("selected header metadata differs")
                    counts["payload_checks"] += 2
                    if gzip.decompress(compressed) != data or zlib.decompress(compressed, wbits=31) != data:
                        problems.append("decompressed bytes differ")
                    decoder = zlib.decompressobj(wbits=31)
                    decoder.decompress(compressed)
                    if not decoder.eof or decoder.unused_data or decoder.unconsumed_tail:
                        problems.append("not exactly one complete gzip member")
                    counts["trailer_checks"] += 1
                    if compressed[-8:] != struct.pack("<II", zlib.crc32(data), len(data) & 0xffffffff):
                        problems.append("CRC32 or ISIZE differs")
                except (ValueError, OSError, EOFError, zlib.error) as error:
                    problems.append("decode/header failure: " + type(error).__name__)
                if repeated is not None:
                    counts["deterministic_repeats"] += 1
                    if compressed != repeated:
                        problems.append("same inputs/options/write schedule changed compressed bytes")
                if corrupt:
                    broken = {"crc": compressed[:-8] + bytes([compressed[-8] ^ 1]) + compressed[-7:],
                              "size": compressed[:-4] + bytes([compressed[-4] ^ 1]) + compressed[-3:]}
                    for length in (1, 9, 10, len(compressed) - 9, len(compressed) - 8, len(compressed) - 1):
                        broken[f"truncated_{length}"] = compressed[:length]
                    for reason, blob in broken.items():
                        (outputs / f"{counts['cases']:05d}_{reason}.gz").write_bytes(blob)
                        for decoder_name, decoder in (("gzip", gzip.decompress), ("zlib", lambda value: zlib.decompress(value, wbits=31))):
                            try:
                                decoder(blob)
                            except (OSError, EOFError, zlib.error):
                                counts["corruption_rejections"] += 1
                            else:
                                problems.append(decoder_name + " accepted " + reason)
            row = {"name": name, "argv": arguments, "actual": actual, "expected_error": error_kind,
                   "input": identity(path), "output": identity(target), "problems": problems}
            record(row)
            return compressed

        try:
            for index, (name, data) in enumerate(messages(args.seed)):
                path = inputs / f"{index}_{name}.bin"
                path.write_bytes(data); fixture_pins.append(identity(path))
                for level in (0, 1, 6, 9):
                    for chunks in ("1", "7,63,64,65,1025", "65536"):
                        baseline = check(f"{name}_L{level}_{chunks}", path, data, level=level, chunks=chunks,
                                         corrupt=level == 9 and chunks == "65536" and index < 5)
                    if level == 9:
                        check(name + "_repeat", path, data, repeated=baseline)
                check(name + "_flush", path, data, chunks="65536", flush_every=1)
            path = inputs / "metadata.bin"
            data = b"a\0b\n" * 257
            path.write_bytes(data); fixture_pins.append(identity(path))
            for mtime in (0, 1, 0xffffffff):
                for filename in (None, b"", b"records.tsv", b"caf\xe9.json"):
                    check("explicit_metadata", path, data, mtime=mtime, filename=filename)

            frame = pd.DataFrame({"gene": ["β", "é", "x\ty", None], "value": [1.5, 0.0, float("nan"), -2.0]})
            for filename in ("records.tsv.gz", "café.tsv.gz", "测试.tsv.gz"):
                python_path = outputs / ("python_" + filename)
                frame.to_csv(python_path, sep="\t", index=False)
                compressed = python_path.read_bytes()
                data, metadata = gzip.decompress(compressed), header(compressed)
                path = inputs / (filename + ".plain")
                path.write_bytes(data); fixture_pins.append(identity(path))
                check("pandas_inferred_gzip", path, data, filename=bytes.fromhex(metadata["filename_hex"]),
                      mtime=metadata["mtime"], expected_python=compressed)
                counts["pandas_cases"] += 1
            for filename in ("table.tsv.tmp", "table.tsv.gz.tmp"):
                path = outputs / filename
                frame.to_csv(path, sep="\t", index=False)
                expected = frame.to_csv(None, sep="\t", index=False).encode("utf-8")
                problems = [] if path.read_bytes() == expected and not expected.startswith(b"\x1f\x8b") else ["pandas temporary suffix compression changed"]
                record({"name": "pandas_plain_tmp", "output": identity(path), "problems": problems})
                counts["plain_pandas_cases"] += 1

            data = random.Random(args.seed + 1).randbytes(4096)
            path = inputs / "faults.bin"
            path.write_bytes(data); fixture_pins.append(identity(path))
            baseline = check("fault_baseline", path, data)
            for at in (0, 3, 9, 10, 11, len(baseline) - 8, len(baseline) - 1):
                for kind, error in (("error", "BrokenPipe"), ("zero", "WriteZero"), ("interrupt", None)):
                    check(f"fault_{at}_{kind}", path, data, options=("--fault-at", str(at), "--fault-kind", kind, "--sink-chunk", "3"), error_kind=error)
                    if kind == "interrupt": counts["interrupted_successes"] += 1
            check("final_flush_error", path, data, options=("--flush-fault", "error"), error_kind="PermissionDenied")
            check("earlier_flush_error", path, data, flush_every=1, options=("--flush-fault", "error"), error_kind="PermissionDenied")
            check("flush_interrupted", path, data, options=("--flush-fault", "interrupt"))
            counts["interrupted_successes"] += 1
            check("invalid_level", path, data, level=10, error_kind="InvalidInput")
            check("invalid_filename", path, data, filename=b"x\0y", error_kind="InvalidInput")

            class BrokenWriter(io.BytesIO):
                def write(self, value):
                    raise OSError("injected Python writer failure")
            class BrokenFlush(io.BytesIO):
                def flush(self):
                    raise OSError("injected Python flush failure")
            for name, target in (("write", BrokenWriter()), ("flush", BrokenFlush())):
                try:
                    with gzip.GzipFile(fileobj=target, mode="wb", mtime=0) as writer:
                        writer.write(data)
                        writer.flush()
                except OSError:
                    counts["python_io_failures"] += 1
                    record({"name": "python_" + name + "_failure", "expected": "OSError", "problems": []})
                else:
                    record({"name": "python_" + name + "_failure", "problems": ["Python did not propagate the injected I/O failure"]})
            process.stdin.close()
            code = process.wait(timeout=10)
            process = None
            if code: raise RuntimeError("probe exited unsuccessfully")
        finally:
            if process is not None:
                process.kill(); process.wait()
    stable = all(identity(Path(pin["path"])) == pin for pin in [*pins.values(), *fixture_pins])
    result = {**counts, "passed": counts["mismatches"] == 0 and stable, "inputs_unchanged": stable,
              "identities": pins, "fixtures": fixture_pins, "python": sys.version,
              "gzip_module": identity(Path(gzip.__file__)), "pandas_version": pd.__version__,
              "zlib_version": zlib.ZLIB_VERSION, "zlib_runtime_version": zlib.ZLIB_RUNTIME_VERSION,
              "seed": args.seed, "wall_seconds": time.monotonic() - started,
              "results": identity(args.output / "results.jsonl"), "mismatch_ledger": identity(args.output / "mismatches.jsonl"),
              "contract": "Exact decompressed bytes, selected header metadata, CRC32/ISIZE and explicit failure reporting; cross-backend compressed bytes are not required to match"}
    with (args.output / "comparison.json").open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True); handle.write("\n")
    print(json.dumps({key: value for key, value in result.items() if key not in ("fixtures", "identities", "python")}), flush=True)
    if not result["passed"]: raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("compare")
    command.add_argument("--binary", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    command.add_argument("--module", type=Path, default=Path(__file__).resolve().parents[1] / "src/gzip_writer.rs")
    command.add_argument("--seed", type=lambda value: int(value, 0), default=0x475A49505F535452)
    command.set_defaults(function=compare)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
