#!/usr/bin/env python3
"""Independent CPython differential for the standalone JSON writer probe.

Compile support/json_text_probe.rs against current and pre-change json_text.rs, using
--cfg legacy for the latter. Run:
  python json_text_compare.py compare --binary PROBE --legacy-binary OLD_PROBE \
      --output RESULTS --random-floats 5000
All requests, responses and differences are retained beneath RESULTS. No PriVA
or scientific Python module is imported. Explicit pre-existing parser refusals
and legacy decimal-format differences are recorded, with exact unchanged legacy
behavior checked separately from the new writer's exact CPython output.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import struct
import subprocess
import sys
import time


def identity(path: Path) -> dict:
    with path.open("rb") as handle:
        sha = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha}


def fixed_inputs():
    strings = ["", "plain", 'quote" slash/ backslash\\', "".join(map(chr, range(32))) + chr(127),
               "é e\u0301 β 中文 😀 𝄞", "\u2028\u2029\ufeff\uffff\U0010ffff"]
    for index, value in enumerate([None, True, False, 0, -1, 2**63, -(2**127), 2**127-1,
                                    [], {}, [[], {}, None], *strings]):
        yield f"scalar_{index}", json.dumps(value), False
    for index, source in enumerate([
        '{"2":2,"10":10,"01":1,"1":1,"0":0,"-1":-1}',
        '{"β":1,"a":2,"é":3,"😀":4,"中":5,"𝄞":6,"e\\u0301":7}',
        '{"z":[{"b":2,"a":1},[],{}],"a":{"2":0,"10":1}}',
        '{"b":0,"a":1,"b":2}',
        '[NaN,Infinity,-Infinity]', '{"z":{"value":NaN},"a":Infinity}',
        'NaN', 'Infinity', '-Infinity', '1e400', '-1e400', '-0.0', '-0',
        '[1e-5,1e-4,1e15,1e16,1e23,0.30000000000000004,5e-324,1.7976931348623157e308]',
    ]):
        yield f"structured_{index}", source, False
    for depth in (1, 2, 20, 100, 399, 400):
        yield f"depth_{depth}", "[" * depth + "0" + "]" * depth, False
    for index, source in enumerate(["", " ", "[1,]", '{"a":1,}', "01", "[1 2]", '{"a" 1}',
                                    '"unterminated', '"\\q"', '"\\u12G4"', '"a\nb"', "{} x",
                                    '\ufeff{}', '{\n"é": β\n}', "[true,NaN,]"]):
        yield f"parse_error_{index}", source, False
    # These are unchanged, explicit limits of json_text's existing parser.
    for name, source in [("integer_beyond_i128", str(2**127)), ("lone_surrogate", '"\\ud800"'),
                         *[(f"depth_{depth}", "[" * depth + "0" + "]" * depth) for depth in (401, 600, 1100)]]:
        yield "preexisting_refusal_" + name, source, True


def modes():
    for indent, sort_keys, ensure_ascii, compact, allow_nan in itertools.product((None, 2), (False, True), (False, True), (False, True), (False, True)):
        yield {"indent": indent, "sort_keys": sort_keys, "ensure_ascii": ensure_ascii,
               "compact": compact, "allow_nan": allow_nan}


def requests(random_floats: int, seed: int):
    all_modes = list(modes())
    for name, source, refusal in fixed_inputs():
        for number, options in enumerate(all_modes):
            yield {"name": name + f"_mode_{number}", "input": source, "mode": "options", "options": options,
                   "preexisting_refusal": refusal}
        for mode in ("legacy_default", "legacy_sorted", "legacy_indented"):
            yield {"name": name + "_" + mode, "input": source, "mode": mode, "preexisting_refusal": refusal}
    # Exact f64 bit patterns, including negative zero and all exponent ranges.
    generator = random.Random(seed)
    samples = [0, 1 << 63, 1, (1 << 63) | 1, 0x0010000000000000, 0x7FEFFFFFFFFFFFFF,
               0x43192587AD7A8149, 0x430710BD04F20BE2, 0xC3192587AD7A8149, 0xC30710BD04F20BE2]
    samples += [generator.getrandbits(64) for _ in range(random_floats)]
    for index, bits in enumerate(samples):
        value = struct.unpack(">d", struct.pack(">Q", bits))[0]
        if not math.isfinite(value):
            continue
        source = json.dumps(value)
        for number in (0, 7, 15, 16, 23, 31):
            yield {"name": f"float_{index}_{bits:016x}_mode_{number}", "input": source,
                   "mode": "options", "options": all_modes[number], "preexisting_refusal": False}
        for mode in ("legacy_default", "legacy_sorted", "legacy_indented"):
            yield {"name": f"float_{index}_{bits:016x}_{mode}", "input": source, "mode": mode, "preexisting_refusal": False}
    # Duplicate rejection is a pre-existing caller-selectable parser contract.
    yield {"name": "reject_nested_duplicate", "input": '{"a":{"z":1,"z":2}}',
           "mode": "options", "options": {}, "reject_duplicates": True, "preexisting_refusal": False}


def python_result(case: dict) -> dict:
    def unique(pairs):
        output = {}
        for key, value in pairs:
            if key in output:
                raise ValueError("Duplicate manifest key: " + key)
            output[key] = value
        return output
    try:
        value = json.loads(case["input"], object_pairs_hook=unique if case.get("reject_duplicates") else None)
    except Exception as error:
        return {"ok": False, "phase": "parse", "error": str(error), "refusal": False}
    mode = case["mode"]
    if mode == "legacy_default":
        options = {}
    elif mode == "legacy_sorted":
        options = {"sort_keys": True}
    elif mode == "legacy_indented":
        options = {"indent": 2}
    else:
        options = {key: value for key, value in case["options"].items() if key != "compact"}
        if case["options"].get("compact"):
            options["separators"] = (",", ":")
    try:
        return {"ok": True, "output": json.dumps(value, **options)}
    except Exception as error:
        return {"ok": False, "phase": "dump", "error": str(error), "refusal": False}


def send(process: subprocess.Popen, case: dict) -> dict:
    process.stdin.write(json.dumps(case, ensure_ascii=True) + "\n")
    process.stdin.flush()
    line = process.stdout.readline()
    if not line:
        raise RuntimeError("native probe ended before returning a response; see its stderr log")
    return json.loads(line)


def compare(args: argparse.Namespace) -> None:
    started = time.monotonic()
    args.output.mkdir(exist_ok=False)
    pins = {key: identity(path) for key, path in (("binary", args.binary), ("legacy_binary", args.legacy_binary), ("comparator", Path(__file__)))}
    counts = {"cases": 0, "supported_cases": 0, "supported_mismatches": 0,
              "preexisting_refusal_differences": 0, "unexpected_refusal_behavior": 0,
              "legacy_comparisons": 0, "legacy_mismatches": 0, "legacy_cpython_differences": 0}
    processes = []
    try:
        with (args.output / "native.stderr.log").open("x") as native_log, (args.output / "legacy.stderr.log").open("x") as legacy_log, \
                (args.output / "requests.jsonl").open("x") as inputs, (args.output / "responses.jsonl").open("x") as outputs, \
                (args.output / "mismatches.jsonl").open("x") as ledger:
            native = subprocess.Popen([str(args.binary), "serve"], text=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=native_log)
            legacy = subprocess.Popen([str(args.legacy_binary), "serve"], text=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=legacy_log)
            processes = [native, legacy]
            for case in requests(args.random_floats, args.seed):
                json.dump(case, inputs, ensure_ascii=True); inputs.write("\n")
                expected = python_result(case)
                actual = send(native, case)
                counts["cases"] += 1
                old = None
                if case["mode"].startswith("legacy_") or case["preexisting_refusal"]:
                    old_case = dict(case, mode=case["mode"] if case["mode"].startswith("legacy_") else "legacy_default")
                    old = send(legacy, old_case)
                    counts["legacy_comparisons"] += 1
                    if old != actual:
                        counts["legacy_mismatches"] += 1
                        json.dump({"kind": "legacy_changed", "case": case, "expected": old, "actual": actual}, ledger); ledger.write("\n")
                equal = expected == actual
                if case["preexisting_refusal"]:
                    if not actual.get("refusal") or actual.get("phase") != "parse" or old != actual:
                        counts["unexpected_refusal_behavior"] += 1
                    if not equal:
                        counts["preexisting_refusal_differences"] += 1
                elif case["mode"].startswith("legacy_"):
                    if not equal:
                        counts["legacy_cpython_differences"] += 1
                else:
                    counts["supported_cases"] += 1
                    if not equal:
                        counts["supported_mismatches"] += 1
                if not equal:
                    kind = ("preexisting_parser_refusal" if case["preexisting_refusal"] else
                            "preexisting_legacy_format" if case["mode"].startswith("legacy_") and old == actual else
                            "cpython_difference")
                    json.dump({"kind": kind,
                               "case": case, "expected": expected, "actual": actual}, ledger, ensure_ascii=True); ledger.write("\n")
                json.dump({"name": case["name"], "actual": actual, "expected": expected, "equal": equal}, outputs, ensure_ascii=True); outputs.write("\n")
            for process in processes:
                process.stdin.close()
            for process in processes:
                if process.wait(timeout=10) != 0:
                    raise RuntimeError("a native probe returned a nonzero status")
            processes = []
    finally:
        for process in processes:
            process.kill(); process.wait()
    stable = all(identity(Path(row["path"])) == row for row in pins.values())
    result = {**counts, "identities": pins, "inputs_unchanged": stable, "python": sys.version,
              "wall_seconds": time.monotonic() - started, "random_floats_requested": args.random_floats, "seed": args.seed,
              "mismatches": str(args.output / "mismatches.jsonl"), "requests": identity(args.output / "requests.jsonl"),
              "responses": identity(args.output / "responses.jsonl")}
    result["passed"] = stable and all(counts[key] == 0 for key in ("supported_mismatches", "unexpected_refusal_behavior", "legacy_mismatches"))
    result["all_cpython_cases_exact"] = result["passed"] and counts["preexisting_refusal_differences"] == 0 and counts["legacy_cpython_differences"] == 0
    with (args.output / "comparison.json").open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True); handle.write("\n")
    print(json.dumps({key: value for key, value in result.items() if key not in ("identities", "python", "requests", "responses")}), flush=True)
    if not result["passed"]:
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("compare")
    command.add_argument("--binary", type=Path, required=True)
    command.add_argument("--legacy-binary", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    command.add_argument("--random-floats", type=int, default=5000)
    command.add_argument("--seed", type=lambda value: int(value, 0), default=0x50726956414A534F)
    command.set_defaults(function=compare)
    args = parser.parse_args()
    if args.random_floats < 0:
        parser.error("--random-floats must be nonnegative")
    args.function(args)


if __name__ == "__main__":
    main()
