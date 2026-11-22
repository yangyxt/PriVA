#!/usr/bin/env python3
"""Compare the standalone exon-AF reductions with the unchanged Python producer.

python exon_af_stats_compare.py compare --probe PATH --source PATH --output DIR
All operands and four result bit patterns are retained, including every mismatch.
"""
from __future__ import annotations
import argparse
import importlib.util
import itertools
import json
import math
from pathlib import Path
import random
import struct
import subprocess
import sys
import warnings
sys.dont_write_bytecode = True


def bits(value):
    return struct.pack(">d", float(value)).hex()


def floating(value):
    return struct.unpack(">d", bytes.fromhex(value))[0]


def cases():
    rng = random.Random(30631)
    lengths = list(range(1, 35)) + [63, 64, 65, 127, 128, 129, 255, 256, 257,
                                  511, 512, 513, 1023, 1024, 1025, 1803, 2049,
                                  8191, 8192, 8193, 16383, 16384, 16385, 65537]
    for n in lengths:
        for name, values in [
            ("zero", [0.0] * n),
            ("one", [1.0] * n),
            ("ordered", [float(i) for i in range(n)]),
            ("reversed", [float(i) for i in range(n - 1, -1, -1)]),
            ("organ_pipe", [float(min(i, n - i)) for i in range(n)]),
            ("alternating", [float(i % 2) for i in range(n)]),
            ("f32", [float(struct.unpack("f", struct.pack("f", rng.random()))[0]) for _ in range(n)]),
            ("f64", [rng.random() for _ in range(n)]),
            ("cancellation", [rng.choice([1e100, -1e100, 1e-100, -1e-100, 1.0]) for _ in range(n)]),
        ]:
            yield f"{name}_{n}", values
        for name, special in [("positive_nan", float("nan")), ("negative_nan", -float("nan")),
                              ("positive_inf", math.inf), ("negative_inf", -math.inf)]:
            for at in sorted({0, n // 2, n - 1}):
                values = [float(i % 17) / 17 for i in range(n)]
                values[at] = special
                yield f"{name}_{n}_{at}", values
        values = [rng.choice([0.0, -0.0]) for _ in range(n)]
        yield f"signed_zero_{n}", values
        values = [rng.choice([1.0, float("nan"), -float("nan")]) for _ in range(n)]
        yield f"mixed_nan_{n}", values
    # Short traversal/tie reductions, including distinct quiet-NaN payloads.
    choices = [0.0, -0.0, 1.0, -1.0, math.inf, -math.inf,
               floating("7ff8000000000001"), floating("fff8000000000002")]
    for n in range(1, 5):
        for index, values in enumerate(itertools.product(choices, repeat=n)):
            yield f"short_{n}_{index}", values
    for index in range(250):
        n = rng.randrange(5, 400)
        yield f"ties_{index}_{n}", [rng.choice(choices) for _ in range(n)]


def compare(args):
    import numpy as np
    args.output.mkdir(parents=True, exist_ok=False)
    spec = importlib.util.spec_from_file_location("exon_source", args.source)
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    records = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name, values in cases():
            expected = source.calculate_stats({"T": {"1/1": list(values)}})["T"]["1/1"]
            records.append({"case": name, "input": [bits(x) for x in values],
                            "python": [bits(x) for x in expected[:4]]})
    request = "".join(" ".join(row["input"]) + "\n" for row in records)
    (args.output / "request.txt").write_text(request)
    result = subprocess.run([str(args.probe), "calculate"], input=request,
                            text=True, capture_output=True, timeout=120)
    (args.output / "native.stdout").write_text(result.stdout)
    (args.output / "native.stderr").write_text(result.stderr)
    lines = result.stdout.splitlines()
    mismatches = []
    for at, record in enumerate(records):
        record["native"] = lines[at].split() if at < len(lines) else None
        record["passed"] = record["native"] == record["python"]
        if not record["passed"]:
            mismatches.append(record)
    for filename, rows in [("results.jsonl", records), ("mismatches.jsonl", mismatches)]:
        with (args.output / filename).open("x") as stream:
            for row in rows:
                stream.write(json.dumps(row, allow_nan=False) + "\n")
    fields = ["median", "mean", "maximum", "minimum"]
    counts = {name: sum(row["native"] is None or row["native"][i] != row["python"][i]
                        for row in mismatches) for i, name in enumerate(fields)}
    report = {"cases": len(records), "mismatches": len(mismatches), "by_field": counts,
              "exit": result.returncode, "numpy": np.__version__,
              "cpu_features": {k: v for k, v in np.core._multiarray_umath.__cpu_features__.items() if k in ["AVX2", "AVX512F", "AVX512_SKX"]},
              "passed": not mismatches and result.returncode == 0 and len(lines) == len(records)}
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)
    return 0 if report["passed"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("compare")
    for name in ("probe", "source", "output"):
        command.add_argument("--" + name, type=Path, required=True)
    command.set_defaults(function=compare)
    args = parser.parse_args()
    return args.function(args)


if __name__ == "__main__":
    raise SystemExit(main())
