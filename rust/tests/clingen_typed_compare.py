#!/usr/bin/env python3
"""Compare typed ClinGen values, old archive bytes, and actual criteria consumers.

  python clingen_typed_compare.py check --binary NEW --baseline-binary OLD \
      --reference-pickle PYTHON.pkl.gz --native-v1 NATIVE.rkyv \
      --converted-v1 CONVERTED.rkyv --work-root NEW_DIRECTORY

No VCF producer is rerun. Existing v1 archives are read only. New synthetic
pickles exercise actual converter/export and ACMG consumer routes; the Python
summary functions are independently source-pinned by differential_scoring.py.
"""
from __future__ import annotations
import argparse
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import pickle
import subprocess
from types import SimpleNamespace

import pandas as pd
import acmg_step_oracle_common as frames
import cache_values_compare as values
import differential_scoring as scoring


def identity(path):
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": values.sha256(path)}


def invoke(binary, args, root, name):
    with (root / f"{name}.stdout.log").open("x") as out, (root / f"{name}.stderr.log").open("x") as err:
        result = subprocess.run([str(binary), *map(str, args)], stdout=out, stderr=err, check=False)
    if result.returncode:
        raise AssertionError(f"{name} returned {result.returncode}; see its saved diagnostics")


def compare(binary, source, archive, root, name):
    with redirect_stdout(io.StringIO()):
        result = values.compare(SimpleNamespace(binary=binary, cache="clingen_map", pickle=source,
            archive=archive, report=root / f"{name}.json", mismatches=root / f"{name}.mismatches.jsonl"))
    if not result["passed"]:
        raise AssertionError(f"{name}: {result['mismatches']} unequal values")
    return result


def consumer(binary, archive, reference_pickle, keys, root, scope, strength, enabled):
    root.mkdir()
    original = pd.DataFrame({"variant_id": keys})
    criteria = {"PVS1": [strength] * len(keys), "PM2": [2 if i % 2 else 0 for i in range(len(keys))]}
    expected_table, expected = scope["summarize_acmg_criteria"](
        original.copy(), criteria, str(reference_pickle), apply_clingen_override=enabled)
    expected["__ACMG_criteria"] = expected_table["ACMG_criteria"]
    frames.write_frame(expected, root / "expected.frame")
    input_frame = original.copy()
    for name, column in criteria.items():
        input_frame["criterion:" + name] = column
    frames.write_frame(input_frame, root / "input.frame")
    status, stderr = frames.run_rust(binary, "acmg_scoring", "summarize_acmg_criteria", root,
        {"clingen_map_pkl": str(archive), "apply_clingen_override": enabled})
    (root / "native.stderr.log").write_text(stderr)
    differences = [] if status else frames.compare_frames(root / "expected.frame", root / "observed.frame", show=10000)
    (root / "differences.json").write_text(json.dumps(differences, indent=2) + "\n")
    if status or differences:
        raise AssertionError(f"consumer check failed: status={status}, differences={len(differences)}")
    return {"rows": len(keys), "columns": len(expected.columns), "mismatches": 0,
            "expected": identity(root / "expected.frame"), "actual": identity(root / "observed.frame")}


def check(args):
    root = args.work_root
    root.mkdir(parents=True, exist_ok=False)
    inputs = {name: identity(getattr(args, name)) for name in
              ("binary", "baseline_binary", "reference_pickle", "native_v1", "converted_v1")}
    result = {"inputs": inputs, "python_reference": frames.manifest([
        Path(__file__).resolve(), Path(frames.__file__).resolve(),
        Path(values.__file__).resolve(), Path(scoring.__file__).resolve(),
        scoring.ROOT / "scripts/acmg_scoring.py",
    ])}
    new_v1 = root / "reconverted.v1.rkyv"
    invoke(args.binary, ["cache-convert", "clingen_map", "--pickle", args.reference_pickle, "--output", new_v1], root, "reconvert")
    if new_v1.read_bytes() != args.converted_v1.read_bytes():
        raise AssertionError("the same full Python result no longer produces identical version-1 archive bytes")
    result["full_v1_conversion_bytes_identical"] = True
    result["v1_value_comparisons"] = {}
    for name, archive in (("native", args.native_v1), ("converted", args.converted_v1), ("reconverted", new_v1)):
        result["v1_value_comparisons"][name] = compare(args.binary, args.reference_pickle, archive, root, name)
        streams = []
        for label, binary in (("old", args.baseline_binary), ("new", args.binary)):
            stream = root / f"{name}.{label}.jsonl"
            invoke(binary, ["cache-export", "clingen_map", "--input", archive, "--output", stream], root, name + "." + label)
            streams.append(stream)
        if streams[0].read_bytes() != streams[1].read_bytes():
            raise AssertionError("version-1 exported reader behavior changed")

    codes = {"text": "PS3, PM2_Supporting", "empty": "", "space": " \t", "integer": 17,
             "negative_integer": -(2**63), "unsigned_integer": 2**64 - 1,
             "float": 0.25, "positive_zero": 0.0, "negative_zero": -0.0,
             "positive_infinity": float("inf"), "negative_infinity": -float("inf"),
             "true": True, "false": False, "missing": float("nan")}
    typed_pickle = root / "typed.pkl"
    typed_pickle.write_bytes(pickle.dumps(codes, protocol=5))
    typed_archive = root / "typed.v2.rkyv"
    invoke(args.binary, ["cache-convert", "clingen_map", "--pickle", typed_pickle, "--output", typed_archive], root, "typed")
    if int.from_bytes(typed_archive.read_bytes()[16:20], "little") != 2:
        raise AssertionError("typed scalars did not select explicit schema version 2")
    result["typed_value_comparison"] = compare(args.binary, typed_pickle, typed_archive, root, "typed")
    # Control v1 keeps all keys but contains no replacement text for non-strings,
    # exactly matching Python's existing isinstance(..., str) override predicate.
    control_pickle = root / "consumer-control.pkl"
    control_pickle.write_bytes(pickle.dumps({key: code if isinstance(code, str) else float("nan") for key, code in codes.items()}, protocol=5))
    control = root / "consumer-control.v1.rkyv"
    invoke(args.baseline_binary, ["cache-convert", "clingen_map", "--pickle", control_pickle, "--output", control], root, "control")
    with redirect_stdout(io.StringIO()):
        scope = scoring.reference_functions()
    keys = [*codes, "absent", "text", "integer"]
    result["consumer_comparisons"] = []
    for strength in (0, 1, 4):
        for enabled in (False, True):
            label = f"strength-{strength}.override-{enabled}"
            new = consumer(args.binary, typed_archive, typed_pickle, keys, root / (label + ".new"), scope, strength, enabled)
            old = consumer(args.baseline_binary, control, typed_pickle, keys, root / (label + ".old"), scope, strength, enabled)
            if Path(new["actual"]["path"]).read_bytes() != Path(old["actual"]["path"]).read_bytes():
                raise AssertionError("typed codes changed the old string-only consumer behavior")
            result["consumer_comparisons"].extend([new, old])
    result["inputs_unchanged"] = all(values.sha256(Path(row["path"])) == row["sha256"] for row in inputs.values())
    result["passed"] = result["inputs_unchanged"]
    (root / "comparison.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"passed": result["passed"], "v1_keys": result["v1_value_comparisons"]["native"]["reference_entries"],
                      "typed_keys": len(codes), "consumer_checks": len(result["consumer_comparisons"]), "report": str(root / "comparison.json")}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("check")
    for name in ("binary", "baseline-binary", "reference-pickle", "native-v1", "converted-v1", "work-root"):
        command.add_argument("--" + name, type=Path, required=True)
    command.set_defaults(function=check)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
