#!/usr/bin/env python3
"""Compare pinned Rust condition masks and metadata with a closed actual-Python oracle.

python condition_masks_compare.py compare --binary PATH --binary-sha256 SHA
  --oracle DATA_ROOT/oracle --oracle-manifest-sha256 SHA --output DATA_ROOT/comparison
All individual failures remain in comparison.json and per-case stderr/stdout logs.
"""
import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.reader(stream))


def metadata_rows(case):
    encode = lambda value: ["null", ""] if value is None else ["str", value.encode().hex()]
    expected = [["property", "position", "kind", "value"]]
    for i, name in enumerate(case["index_names"]):
        expected.append(["index_name", str(i), *encode(name)])
    expected.append(["index_dtype", "0", *encode(case["index_dtype"])])
    for i, item in enumerate(case["series_metadata"]):
        assert item["index_names"] == case["index_names"] and item["index_dtype"] == case["index_dtype"]
        expected.extend([["series_name", str(i), *encode(item["name"])],
                         ["series_dtype", str(i), *encode(item["dtype"])]])
    return expected


def command(args, source, output, metadata, case):
    cmd = [str(args.binary), "condition-masks", "--input", str(source), "--output", str(output),
           "--metadata-output", str(metadata), "--index-dtype", case["index_dtype"]]
    if case["index_names"][0] is not None:
        cmd.extend(["--index-name", case["index_names"][0]])
    for position, dtype in enumerate(case["dtypes"]):
        cmd.extend(["--column-dtype", f"{position}={dtype}"])
    return cmd


def compare(args):
    args.output.mkdir(parents=True, exist_ok=False)
    manifest_path = args.oracle / "manifest.json"
    assert digest(manifest_path) == args.oracle_manifest_sha256, "oracle manifest pin changed"
    manifest = json.loads(manifest_path.read_text())

    def verify():
        assert digest(args.binary) == args.binary_sha256, "binary pin changed"
        for name, sha in manifest.items():
            assert digest(args.oracle / name) == sha, name

    verify()
    oracle = json.loads((args.oracle / "oracle.json").read_text())
    results = []
    for case in oracle["cases"]:
        name = case["name"]
        source = args.oracle / (name + ".input.csv")
        output = args.output / (name + ".actual.csv")
        metadata = args.output / (name + ".metadata.csv")
        cmd = command(args, source, output, metadata, case)
        entry = {"name": name, "rows": case["rows"], "argv": cmd, "passed": False}
        try:
            process = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            (args.output / (name + ".stderr.log")).write_text(process.stderr)
            (args.output / (name + ".stdout.log")).write_text(process.stdout)
            entry["returncode"] = process.returncode
            if "exception" in case:
                error = case["exception"]
                expected = f"priva: {error['type']}: {error['message']}"
                entry["passed"] = (process.returncode == 2 and process.stderr.strip() == expected
                                   and not output.exists() and not metadata.exists())
                entry["expected_exception"] = error
                if not entry["passed"]:
                    entry["actual_stderr"] = process.stderr
            else:
                expected = rows(args.oracle / (name + ".expected.csv"))
                actual = rows(output) if output.exists() else None
                actual_metadata = rows(metadata) if metadata.exists() else None
                entry["values_equal"] = actual == expected
                entry["metadata_equal"] = actual_metadata == metadata_rows(case)
                entry["passed"] = process.returncode == 0 and entry["values_equal"] and entry["metadata_equal"]
                if not entry["passed"]:
                    entry["actual_stderr"] = process.stderr
                    entry["expected_metadata"] = metadata_rows(case)
                    entry["actual_metadata"] = actual_metadata
                    entry["first_different_row"] = next((i for i, (a, e) in enumerate(zip(actual or [], expected)) if a != e), None)
                    entry["actual_row_count"] = len(actual) if actual is not None else None
                    entry["expected_row_count"] = len(expected)
        except Exception as error:
            entry["comparison_failure"] = {"type": type(error).__name__, "message": str(error)}
        results.append(entry)
    # Boundary safeguards are distinct from original-Python semantics.
    valid = next(c for c in oracle["cases"] if "exception" not in c)
    source = args.oracle / (valid["name"] + ".input.csv")
    sentinel = args.output / "existing-output.csv"
    sentinel.write_bytes(b"preserved existing output\n")
    for label, extra, output in [
        ("existing_output", [], sentinel),
        ("input_alias", [], source),
        ("unsupported_dtype", ["--column-dtype", "0=category"], args.output / "unsupported.csv"),
        ("unpaired_option", ["--index-name"], args.output / "unpaired.csv"),
    ]:
        meta = args.output / (label + ".metadata.csv")
        cmd = [str(args.binary), "condition-masks", "--input", str(source), "--output", str(output),
               "--metadata-output", str(meta), *extra]
        process = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        (args.output / (label + ".stderr.log")).write_text(process.stderr)
        results.append({"name": "transport_" + label, "passed": process.returncode == 2 and not meta.exists()
                        and (output in [sentinel, source] or not output.exists()), "returncode": process.returncode})
    assert sentinel.read_bytes() == b"preserved existing output\n"
    verify()
    result = {"completed": all(c["passed"] for c in results), "cases": results,
              "binary_sha256": args.binary_sha256, "oracle_manifest_sha256": args.oracle_manifest_sha256,
              "oracle_sha256": digest(args.oracle / "oracle.json"),
              "comparator_sha256": digest(Path(__file__)), "scientific_acceptance": False,
              "production_backend_selected": False, "scope": oracle["scope"]}
    (args.output / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"completed": result["completed"], "passed": sum(c["passed"] for c in results),
                      "cases": len(results), "failed": [c["name"] for c in results if not c["passed"]]}))
    if not result["completed"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    command = sub.add_parser("compare")
    for name in ["binary", "oracle", "output"]:
        command.add_argument("--" + name, type=Path, required=True)
    for name in ["binary-sha256", "oracle-manifest-sha256"]:
        command.add_argument("--" + name, required=True)
    compare(parser.parse_args())


if __name__ == "__main__":
    main()
