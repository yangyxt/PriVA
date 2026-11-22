#!/usr/bin/env python3
"""Compare native string hashes and set iteration with CPython for explicit seeds.

  python python_set_order_compare.py check --binary PROBE --work-root NEW_DIR
The standalone probe compiles the production module directly, without a Python
binding. Cases cover collisions, resizes, the large-set growth rule and Unicode
character widths; all expected values come from independent Python processes.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


def cases():
    for size in [0, 1, 4, 5, 18, 19, 76, 77, 306, 307, 1228, 1229, 50000, 50001, 78643]:
        values = [f"ENSG:DOMAIN_{index}" for index in range(size)]
        yield f"forward_{size}", values
        yield f"reverse_{size}", values[::-1]
    yield "unicode", ["", "a", "é", "α", "𐀀", "😀", "Aé", "Aα", "A😀", "\0", "é", "a"]
    yield "duplicate_resize", [f"ENST_{index%97}" for index in range(1000)]


def reference(args):
    if not (sys.version_info[:2] == (3, 12) and sys.byteorder == "little" and sys.hash_info.algorithm == "siphash13"
            and sys.hash_info.width == 64 and sys.hash_info.cutoff == 0):
        raise RuntimeError("reference requires the declared CPython 3.12, 64-bit little-endian SipHash13 profile")
    args.work_root.mkdir()
    for name, values in cases():
        (args.work_root/(name+".input")).write_text("" if not values else "\n".join(values)+"\n")
        ordered = list(set(values))
        (args.work_root/(name+".order")).write_text("" if not ordered else "\n".join(ordered)+"\n")
        (args.work_root/(name+".hashes")).write_text("".join(str(hash(value))+"\n" for value in values))


def check(args):
    args.work_root.mkdir(parents=True, exist_ok=False)
    rows = []
    for seed in [0, 1, 42, 4294967295]:
        folder = args.work_root/str(seed)
        subprocess.run([sys.executable, __file__, "reference", "--work-root", str(folder)],
                       env=dict(os.environ, PYTHONHASHSEED=str(seed)), check=True)
        for name, values in cases():
            output, hashes = folder/(name+".native.order"), folder/(name+".native.hashes")
            result = subprocess.run([str(args.binary), "--input", str(folder/(name+".input")), "--output", str(output),
                                     "--hashes", str(hashes), "--seed", str(seed)], capture_output=True)
            (folder/(name+".stderr.log")).write_bytes(result.stderr)
            equal_order = result.returncode == 0 and output.read_bytes() == (folder/(name+".order")).read_bytes()
            equal_hashes = result.returncode == 0 and hashes.read_bytes() == (folder/(name+".hashes")).read_bytes()
            rows.append({"seed": seed, "case": name, "values": len(values), "returncode": result.returncode,
                         "equal_order": equal_order, "equal_hashes": equal_hashes})
    report = {"cases": len(rows), "hash_values": sum(row["values"] for row in rows),
              "mismatches": sum(not (row["equal_order"] and row["equal_hashes"]) for row in rows), "results": rows}
    (args.work_root/"comparison.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({key: value for key, value in report.items() if key != "results"}))
    if report["mismatches"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name, function in [("reference", reference), ("check", check)]:
        sub = commands.add_parser(name)
        sub.add_argument("--work-root", type=Path, required=True)
        if name == "check":
            sub.add_argument("--binary", type=Path, required=True)
        sub.set_defaults(function=function)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
