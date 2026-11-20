#!/usr/bin/env python3
"""Compare the original AlphaMissense projection and native entry point.

  python alphamissense_vcf_compare.py check --producer SCRIPT --binary PRIVA \
      --work-root NEW_DIRECTORY
  python alphamissense_vcf_compare.py probe --producer SCRIPT --work-root NEW_DIRECTORY

The original shell invocation is used only on controlled, ordinary test paths.
Native literal-path checks never pass shell syntax to the original script.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pysam

PRODUCER_SHA256 = "da52271f22c7187be3dcaba76a20520a2cb70b4029c0032965579c99ae286cbf"
HEADER = "# publisher metadata\n#CHROM\tPOS\tREF\tALT\tSTRAND\tUNIPROT\tTRANSCRIPT\tPVAR\tSCORE\tCLASS\n"
ROW = "chr1\t10\tA\tT\t+\tP1\tENST1\tp.A1T\t1e-06\tlikely_benign\n"
ROW2 = "chr2\t20\tG\tC\t-\tP2\tENST2\tp.G2C\t0.99000\tlikely_pathogenic\n"


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def identity(path):
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest(path)}


def cases():
    return {
        "ordinary": {}, "duplicates": {"text": HEADER + ROW + ROW + ROW2},
        "literal_scores": {"text": HEADER + ROW.replace("1e-06", "-0.000") + ROW2.replace("0.99000", "NaN")},
        "extra_fields": {"text": HEADER + ROW.rstrip("\n") + "\textra\tunused\n"},
        "skip_first_data_line": {"text": ROW + ROW2},
        "skip_first_blank_line": {"text": "\n" + ROW + ROW2},
        "comment_between_records": {"text": HEADER + ROW + "# ignored\n" + ROW2},
        "whitespace": {"text": HEADER + " \t" + ROW.replace("\t", "\u001c\u00a0") + ROW2.replace("\t", "  ")},
        "unicode_whitespace": {"text": HEADER + ROW.replace("\t", "\u00a0\u2003") + ROW2.replace("\t", "  ")},
        "late_control_whitespace": {"text": HEADER + "#" + "padding" * 160 + "\n" + ROW.replace("\t", "\u001c")},
        "crlf": {"text": (HEADER + ROW + ROW2).replace("\n", "\r\n")},
        "cr_only": {"text": (HEADER + ROW + ROW2).replace("\n", "\r")},
        "no_final_newline": {"text": (HEADER + ROW + ROW2).rstrip("\n")},
        "blank_record": {"text": HEADER + ROW + "\n" + ROW2},
        "space_before_comment": {"text": HEADER + " # not a comment\n" + ROW},
        "short_record": {"text": HEADER + ROW + "chr2\t20\tG\tC\n"},
        "empty_stream": {"text": ""}, "only_first_line": {"text": "discarded\n"},
        "no_index": {"index": "missing"}, "only_csi": {"index": "csi"},
        "csi_as_tbi": {"index": "csi_as_tbi"}, "corrupt_index": {"index": "corrupt"},
        "ordinary_gzip": {"gzip": True}, "not_gzip": {"plain": True},
        "bad_utf8": {"bytes": (HEADER + ROW).encode() + b"\xff\n"},
        "utf8_cross_block": {"text": "#" + "a" * 8190 + "β\n" + ROW},
        "bad_output_suffix": {"suffix": ".vcf"},
        "bgzip_failure": {"bgzip": 17}, "tabix_failure": {"tabix": 19},
        "bgzip_missing": {"bgzip": "missing"}, "tabix_missing": {"tabix": "missing"},
        "output_is_directory": {"output_directory": True},
        "unsorted_records": {"text": HEADER + ROW.replace("10", "20") + ROW.replace("10", "10")},
    }


def prepare(root, case, tools):
    root.mkdir(parents=True)
    plain, source = root / "source.tsv", root / "source.tsv.gz"
    plain.write_text(HEADER + ROW + ROW2)
    pysam.tabix_compress(str(plain), str(source), force=True)
    pysam.tabix_index(str(source), seq_col=0, start_col=1, end_col=1, force=True)
    if case.get("index") in ("csi", "csi_as_tbi"):
        pysam.tabix_index(str(source), seq_col=0, start_col=1, end_col=1, force=True, csi=True)
        Path(str(source) + ".tbi").unlink()
        if case["index"] == "csi_as_tbi":
            Path(str(source) + ".csi").rename(Path(str(source) + ".tbi"))
    elif case.get("index") == "missing":
        Path(str(source) + ".tbi").unlink()
    elif case.get("index") == "corrupt":
        Path(str(source) + ".tbi").write_bytes(b"not an index")
    body = case.get("bytes", case.get("text", HEADER + ROW + ROW2).encode())
    plain.write_bytes(body)
    if case.get("plain"):
        source.write_bytes(body)
    elif case.get("gzip"):
        with gzip.open(source, "wb") as stream:
            stream.write(body)
    else:
        pysam.tabix_compress(str(plain), str(source), force=True)
    child_path = root / "tools"
    child_path.mkdir()
    for name in ("bgzip", "tabix"):
        choice = case.get(name)
        if choice == "missing":
            continue
        path = child_path / name
        if choice is None:
            path.symlink_to(tools[name])
        else:
            path.write_text(f"#!/bin/sh\nexit {choice}\n")
            path.chmod(0o755)
    return source, {**os.environ, "PATH": str(child_path)}


def invoke(command, root, label, environment):
    with (root / (label + ".stdout.log")).open("x") as out, (root / (label + ".stderr.log")).open("x") as err:
        return subprocess.run(list(map(str, command)), env=environment, stdout=out, stderr=err, check=False).returncode


def artifact_state(output):
    result = {}
    for name, path in (("output", output), ("temporary", Path(str(output) + ".tmp.vcf")),
                       ("index", Path(str(output) + ".tbi"))):
        result[name] = {"exists": path.exists(), "directory": path.is_dir()}
        if path.is_file():
            result[name].update(bytes=path.stat().st_size, sha256=digest(path))
            if name == "output":
                try:
                    with gzip.open(path, "rb") as stream:
                        body = stream.read()
                    result[name].update(decompressed_sha256=hashlib.sha256(body).hexdigest(), lines=body.count(b"\n"))
                except Exception as error:
                    result[name]["gzip_error"] = type(error).__name__
    return result


def check(args):
    if digest(args.producer) != PRODUCER_SHA256:
        raise RuntimeError("original Python producer identity changed")
    args.work_root.mkdir(parents=True, exist_ok=False)
    tools = {name: shutil.which(name) for name in ("bgzip", "tabix")}
    if not all(tools.values()):
        raise RuntimeError("the existing bgzip and tabix executables are required")
    results = {}
    for name, case in cases().items():
        root = args.work_root / name
        source, environment = prepare(root, case, tools)
        row = {}
        for label in (["python", "native"] if args.binary else ["python"]):
            output = root / (label + case.get("suffix", ".vcf.gz"))
            if case.get("output_directory"):
                output.mkdir()
            else:
                output.write_bytes(b"previous output")
            command = ([sys.executable, args.producer, source, output] if label == "python" else
                       [args.binary, "alphamissense-vcf", "build", "--input", source, "--output", output])
            row[label] = {"status": invoke(command, root, label, environment), "artifacts": artifact_state(output)}
        if args.binary:
            row["passed"] = row["python"] == row["native"]
        results[name] = row
    report = {"producer": identity(args.producer), "binary": identity(args.binary) if args.binary else None,
              "cases": results, "count": len(results), "passed": all(row.get("passed", True) for row in results.values())}
    if args.binary:
        root = args.work_root / "literal_paths"
        source, environment = prepare(root, {}, tools)
        output = root / "literal $name;{}.vcf.gz"
        status = invoke([args.binary, "alphamissense-vcf", "build", "--input", source, "--output", output], root, "native", environment)
        observed = artifact_state(output)
        report["literal_paths"] = {"status": status, "passed": status == 0 and observed == results["ordinary"]["python"]["artifacts"],
            "scope": "literal paths are supported without shell interpretation; original shell expansion is deliberately not reproduced"}
        report["passed"] = report["passed"] and report["literal_paths"]["passed"]
    (args.work_root / "comparison.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"passed": report["passed"], "count": report["count"],
                      "outcomes": {name: {label: value["status"] for label, value in row.items() if label in ("python", "native")} for name, row in results.items()},
                      "failed": [name for name, row in results.items() if not row.get("passed", True)]}, indent=2))
    if not report["passed"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("probe", "check"):
        command = commands.add_parser(name)
        command.add_argument("--producer", type=Path, required=True)
        command.add_argument("--binary", type=Path, required=name == "check")
        command.add_argument("--work-root", type=Path, required=True)
        command.set_defaults(function=check)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
