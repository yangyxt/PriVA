#!/usr/bin/env python3
"""Compare `priva alphamissense-cache` with its unchanged Python source.

check --binary PRIVA --repo CHECKOUT --work-root NEW_DIRECTORY --hg19 PUBLISHER_FILE [--only CASE ...]
oracle --producer SCRIPT --exception JSON -- ARGUMENTS (internal reference launcher)

The publisher size and MD5 are pinned, so success paths need the real hg19
prescore (PUBLISHER_FILE, read only). Each case builds two identical fixture
folders; the real file enters them as a hard link (a copy when linking is not
possible). No case writes into it: the script and the native command only
copy from it or rename the link. Setup steps that need an earlier
installation run the unchanged script in each folder. No case reaches the
network: download cases stop at their tool checks.

The Python arm runs CHECKOUT/scripts/alphamissense_cache.py through the oracle
launcher, which records the class and text of the exception its main()
reports through parser.exit; the native arm runs `PRIVA alphamissense-cache`
with the same arguments. Compared exactly: exit status, stdout, the progress
lines on stderr, the exception category and message (case folder written as
<ARM>, retained-artifact timestamps as <STAMP>), and the fixture tree after
the run: every entry's kind, mode and content (files over 64 KiB by size and
the SHA-256 of their first and last MiB; receipts with identities equal to
the current file's written as <CURRENT>). Argument errors exit 2 in both.
Writes WORK_ROOT/complete.json, removes each case's large fixtures, and exits
1 on any mismatch.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
PRODUCER_SHA256 = "0daae6928374ada5402e31f2bdea4fb01757fe4ac60932b7552bcb940206d854"
PROGRESS = ("Reusing ", "Recovering ", "Retained ", "Downloading ")
STAMP = re.compile(r"\.invalid-\d{8}T\d{6}\.\d{6}Z")
HG19_BYTES = 622293310


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def oracle(args):
    """Run the unchanged main() and record the exception it reports."""
    module = load(args.producer, "original_alphamissense_cache")
    record = {}
    original_exit = module.argparse.ArgumentParser.exit

    def observed_exit(parser, status=0, message=None):
        # main() reports a caught failure with parser.exit(1, ...) inside its
        # except clause, where sys.exc_info() is that exception.
        kind, value, _ = sys.exc_info()
        if status == 1 and kind is not None:
            record.update(category=kind.__name__, message=str(value))
            write_json(args.exception, record)
        return original_exit(parser, status, message)

    module.argparse.ArgumentParser.exit = observed_exit
    sys.argv = [str(args.producer), *args.arguments]
    try:
        module.main()
    except SystemExit:
        raise
    except BaseException as exc:
        record.update(category=type(exc).__name__, message=str(exc), uncaught=True)
        write_json(args.exception, record)
        raise
    raise SystemExit(0)


# ---------------------------------------------------------------- fixtures

def vep_vcf(header='##VEP="v116" API="v116" assembly="GRCh38.p14" cache="x"', extra=""):
    return f"##fileformat=VCFv4.2\n{header}\n{extra}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t5\t.\tA\tC\t.\t.\t.\n"


def place_real(real, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(real, target)
    except OSError:
        shutil.copyfile(real, target)


def materialize(arm, case, real):
    arm.mkdir(parents=True)
    for relative, spec in case.get("files", {}).items():
        path = arm / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if spec == "REAL":
            place_real(real, path)
        elif isinstance(spec, tuple) and spec[0] == "sparse":
            with path.open("wb") as handle:
                handle.truncate(spec[1])
        elif isinstance(spec, tuple) and spec[0] == "gzip":
            path.write_bytes(gzip.compress(spec[1].encode(), mtime=0))
        elif isinstance(spec, tuple) and spec[0] == "directory":
            path.mkdir()
        else:
            path.write_text(spec)


def tools(work_root, names):
    folder = work_root / ("tools-" + "-".join(names) if names else "tools-none")
    if not folder.exists():
        folder.mkdir()
        for name in names:
            found = shutil.which(name)
            if found:
                (folder / name).symlink_to(found)
    return folder


def catalogue():
    cache = "{arm}/cache"
    install = ["install", "--cache-dir", cache, "--assembly", "hg19"]
    real_source = {"src/AlphaMissense_hg19.tsv.gz": "REAL"}
    managed = "cache/AlphaMissense/AlphaMissense_hg19.tsv.gz"
    annotation = ["check-annotation", "--input", "{arm}/a.vcf", "--assembly", "hg38", "--vep-major", "116"]
    setup_install = [["install", "--cache-dir", cache, "--assembly", "hg19", "--source",
                      "{arm}/src/AlphaMissense_hg19.tsv.gz"]]
    return {
        "argument_unknown_command": {"args": ["fetch"]},
        "argument_install_without_assembly": {"args": ["install", "--cache-dir", cache]},
        "validate_unknown_assembly": {"args": ["validate", "--input", "{arm}/x.tsv.gz", "--assembly", "T2T"]},
        "validate_missing": {"args": ["validate", "--input", "{arm}/x.tsv.gz", "--assembly", "hg19"]},
        "validate_empty": {"files": {"x.tsv.gz": ""}, "args": ["validate", "--input", "{arm}/x.tsv.gz", "--assembly", "GRCh37"]},
        "validate_partial_name": {"files": {"x.Partial/a.tsv.gz": "x"},
                                  "args": ["validate", "--input", "{arm}/x.Partial/a.tsv.gz", "--assembly", "hg19"]},
        "validate_size_mismatch": {"files": {"x.tsv.gz": "x"}, "args": ["validate", "--input", "{arm}/x.tsv.gz", "--assembly", "hg38"]},
        "validate_md5_mismatch": {"files": {"x.tsv.gz": ("sparse", HG19_BYTES)},
                                  "args": ["validate", "--input", "x.tsv.gz", "--assembly", "hg19"]},
        "validate_real": {"files": real_source, "real": True,
                          "args": ["validate", "--input", "{arm}/src/AlphaMissense_hg19.tsv.gz", "--assembly", "hg19"]},
        "validate_real_other_assembly": {"files": real_source, "real": True,
                                         "args": ["validate", "--input", "{arm}/src/AlphaMissense_hg19.tsv.gz",
                                                  "--assembly", "hg38"]},
        "annotation_plain": {"files": {"a.vcf": vep_vcf()}, "args": annotation},
        "annotation_gzip": {"files": {"a.vcf": ("gzip", vep_vcf())}, "args": annotation},
        "annotation_other_major": {"files": {"a.vcf": vep_vcf()}, "args": [*annotation[:-1], "115"]},
        "annotation_major_text": {"files": {"a.vcf": vep_vcf()}, "args": [*annotation[:-1], "v116"]},
        "annotation_other_assembly": {"files": {"a.vcf": vep_vcf()}, "args": [*annotation[:4], "GRCh37", *annotation[5:]]},
        "annotation_without_vep": {"files": {"a.vcf": vep_vcf("##source=x")}, "args": annotation},
        "annotation_two_vep": {"files": {"a.vcf": vep_vcf(extra='##VEP="v116" assembly="GRCh38"\n')}, "args": annotation},
        "annotation_missing": {"args": annotation},
        "annotation_no_column_header": {"files": {"a.vcf": "##fileformat=VCFv4.2\n1\t5\n"}, "args": annotation},
        "install_partial_cache": {"args": ["install", "--cache-dir", "{arm}/c.partial", "--assembly", "hg19"]},
        "install_parent_reference": {"args": ["install", "--cache-dir", "{arm}/a/../c", "--assembly", "hg19"]},
        "install_unknown_assembly": {"args": [*install[:4], "T2T"]},
        "install_without_tools": {"path": [], "args": install},
        "install_without_curl": {"path": ["bgzip", "tabix"], "args": install},
        "install_source_too_small": {"files": {"src/x.tsv.gz": "x"},
                                     "args": [*install, "--source", "{arm}/src/x.tsv.gz"]},
        "install_external_real": {"files": real_source, "real": True,
                                  "args": [*install, "--source", "{arm}/src/AlphaMissense_hg19.tsv.gz"]},
        "install_reuse": {"files": real_source, "real": True, "setup": setup_install,
                          "args": [*install, "--source", "{arm}/src/AlphaMissense_hg19.tsv.gz"]},
        "install_reuse_without_source": {"files": real_source, "real": True, "setup": setup_install,
                                         "args": install},
        "install_existing_valid_reindex": {"files": {managed: "REAL", managed + ".tbi": "stale index"},
                                           "real": True, "args": install},
        "install_existing_invalid_without_curl": {"files": {managed: "not the publisher file"},
                                                  "path": ["bgzip", "tabix"], "args": install},
    }


# ---------------------------------------------------------------- runs

def normalize(text, arm):
    if not isinstance(text, str):
        return text
    for spelling in sorted({str(arm), os.path.realpath(arm)}, key=len, reverse=True):
        text = text.replace(spelling, "<ARM>")
    return STAMP.sub(".invalid-<STAMP>", text)


def identity(path):
    value = path.stat()
    return {"inode": value.st_ino, "bytes": value.st_size,
            "mtime_ns": value.st_mtime_ns, "ctime_ns": value.st_ctime_ns}


def receipt_state(data, root):
    try:
        body = json.loads(data)
    except ValueError:
        return normalize(data.decode("utf-8", "replace"), root)
    if isinstance(body, dict) and isinstance(body.get("path"), str):
        for key, path in (("file_identity", Path(body["path"])), ("index_identity", Path(body["path"] + ".tbi"))):
            if path.exists() and body.get(key) == identity(path):
                body[key] = "<CURRENT>"
    return normalize(json.dumps(body, indent=2), root)


def tree(root):
    state = {}
    for folder, directories, names in os.walk(root, followlinks=False):
        for name in sorted(directories + names):
            path = Path(folder) / name
            relative = normalize(str(path.relative_to(root)), root)
            if path.is_symlink():
                state[relative] = ["link", normalize(os.readlink(path), root)]
                continue
            mode = oct(path.stat().st_mode & 0o7777)
            if path.is_dir():
                state[relative] = ["directory", mode]
                continue
            size = path.stat().st_size
            if name.endswith(".ready.json"):
                state[relative] = ["receipt", receipt_state(path.read_bytes(), root), mode]
            elif size <= 65536:
                data = path.read_bytes()
                try:
                    state[relative] = ["text", normalize(data.decode("utf-8"), root), mode]
                except UnicodeDecodeError:
                    state[relative] = ["bytes", hashlib.sha256(data).hexdigest(), mode]
            else:
                with path.open("rb") as handle:
                    head = handle.read(1 << 20)
                    handle.seek(max(size - (1 << 20), 0))
                    tail = handle.read()
                state[relative] = ["large", size, hashlib.sha256(head + tail).hexdigest(), mode]
    return state


def run(command, folder, label, environment, cwd):
    with (folder / (label + ".stdout")).open("wb") as out, (folder / (label + ".stderr")).open("wb") as err:
        process = subprocess.run(list(map(str, command)), cwd=cwd, stdout=out, stderr=err,
                                 stdin=subprocess.DEVNULL, env=environment, timeout=3600)
    return process.returncode


def outcome(folder, label, status, arm):
    stdout = (folder / (label + ".stdout")).read_bytes().decode("utf-8", "replace")
    stderr = (folder / (label + ".stderr")).read_bytes().decode("utf-8", "replace")
    # validate prints the file's identity; links share it, copies do not.
    stdout = re.sub(r'"inode": \d+, "bytes": \d+, "mtime_ns": \d+, "ctime_ns": \d+', "<IDENTITY>", stdout)
    progress = [normalize(line, arm) for line in stderr.splitlines() if line.startswith(PROGRESS)]
    result = {"exit": status, "stdout": normalize(stdout, arm), "progress": progress, "tree": tree(arm)}
    if label == "python":
        recorded = folder / "exception.json"
        value = json.loads(recorded.read_text()) if recorded.exists() else {}
        if value:
            result.update(category=value.get("category"), message=value.get("message"))
    else:
        matched = re.search(r"(?:\A|\n)priva: ([A-Za-z_][A-Za-z_0-9]*): ((?:(?!\npriva: ).)*)\n\Z", stderr, re.S)
        if matched:
            result.update(category=matched[1], message=matched[2])
        elif status not in (0, 2):
            result.update(category=None, message=stderr)
    if "message" in result:
        result["message"] = normalize(result["message"], arm)
    return result


def compare(python, native, before):
    if python["exit"] == 2 or native["exit"] == 2:
        return python["exit"] == native["exit"] == 2 and python["tree"] == native["tree"] == before
    left, right = dict(python), dict(native)
    left.pop("uncaught", None)
    return left == right


def check(args):
    producer = args.repo / "scripts/alphamissense_cache.py"
    if hashlib.sha256(producer.read_bytes()).hexdigest() != PRODUCER_SHA256:
        raise SystemExit("alphamissense_cache.py differs from the qualified source; update the pin")
    cases = catalogue()
    unknown = set(args.only or ()) - cases.keys()
    if unknown:
        raise SystemExit("Unknown cases: " + ", ".join(sorted(unknown)))
    real = args.hg19.resolve()
    if real.stat().st_size != HG19_BYTES:
        raise SystemExit("--hg19 is not the pinned hg19 publisher file")
    args.work_root.mkdir(parents=True, exist_ok=False)
    rows = {}
    for name, case in cases.items():
        if args.only and name not in args.only:
            continue
        folder = args.work_root / name
        folder.mkdir()
        results, befores = {}, {}
        for label in ("python", "native"):
            arm = (folder / label).absolute()
            materialize(arm, case, real)
            environment = {key: value for key, value in os.environ.items() if not key.startswith("PRIVA_")}
            environment["PYTHONDONTWRITEBYTECODE"] = "1"
            if "path" in case:
                environment["PATH"] = str(tools(args.work_root, case["path"]))
            for step in case.get("setup", []):
                values = [value.replace("{arm}", str(arm)) for value in step]
                done = subprocess.run([sys.executable, producer, *values], cwd=arm, env=environment,
                                      capture_output=True, text=True)
                if done.returncode != 0:
                    raise SystemExit(f"setup failed in {arm}: {done.stderr}")
            befores[label] = tree(arm)
            arguments = [value.replace("{arm}", str(arm)) for value in case["args"]]
            if label == "python":
                command = [sys.executable, Path(__file__).resolve(), "oracle", "--producer", producer,
                           "--exception", folder / "exception.json", "--", *arguments]
            else:
                command = [args.binary, "alphamissense-cache", *arguments]
            status = run(command, folder, label, environment, arm)
            results[label] = outcome(folder, label, status, arm)
        passed = befores["python"] == befores["native"] and \
            compare(results["python"], results["native"], befores["native"])
        rows[name] = {"passed": passed, "python": results["python"], "native": results["native"]}
        write_json(folder / "comparison.json", rows[name])
        if case.get("real") or "sparse" in json.dumps(case.get("files", {})):
            for label in ("python", "native"):
                shutil.rmtree(folder / label)
    failed = sorted(name for name, row in rows.items() if not row["passed"])
    report = {
        "producer": str(producer), "producer_sha256": PRODUCER_SHA256,
        "binary": str(args.binary),
        "binary_sha256": hashlib.sha256(Path(args.binary).read_bytes()).hexdigest(),
        "cases": len(rows), "python_successes": sum(row["python"]["exit"] == 0 for row in rows.values()),
        "failed": failed, "rows": rows,
    }
    write_json(args.work_root / "complete.json", report)
    print(json.dumps({key: report[key] for key in ("cases", "python_successes", "failed")}))
    if failed:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="procedure", required=True)
    p = commands.add_parser("check")
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--work-root", type=Path, required=True)
    p.add_argument("--hg19", type=Path, required=True)
    p.add_argument("--only", action="append")
    p.set_defaults(function=check)
    p = commands.add_parser("oracle")
    p.add_argument("--producer", type=Path, required=True)
    p.add_argument("--exception", type=Path, required=True)
    p.add_argument("arguments", nargs=argparse.REMAINDER)
    p.set_defaults(function=oracle)
    args = parser.parse_args()
    if args.procedure == "oracle" and args.arguments[:1] == ["--"]:
        args.arguments = args.arguments[1:]
    args.function(args)


if __name__ == "__main__":
    main()
