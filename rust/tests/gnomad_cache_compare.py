#!/usr/bin/env python3
"""Compare `priva gnomad-cache` with its unchanged Python source.

check --binary PRIVA --repo CHECKOUT --work-root NEW_DIRECTORY [--only CASE ...]
oracle --producer SCRIPT --exception JSON -- ARGUMENTS (internal reference launcher)

Each case builds two identical fixture folders holding a generated
24-chromosome family (headers with every caller-required INFO field, one
record per member, bgzip plus a tabix or CSI index written by pysam), and
optionally a reference .fai. No gnomAD data is read. The Python arm runs
CHECKOUT/scripts/gnomad_cache.py through the oracle launcher, which records
the class and text of the exception main() reports; the native arm runs
`PRIVA gnomad-cache` with the same arguments, both from the case folder.
Compared exactly: exit status, stdout and the exception category and message
(case folder written as <ARM>), and the fixture tree before and after (kind,
permission bits, content digest and modification time), which must not
change. Argument errors exit 2 in both; argparse's usage text is not
compared. Writes WORK_ROOT/complete.json and exits 1 on any mismatch.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys

sys.dont_write_bytecode = True
PRODUCER_SHA256 = "2a1991059f944a85c78a4b86c37fe317124085041fa168848e3c6b21db2619a2"
FIXED_NS = 1_700_000_000_000_000_000


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def oracle(args):
    """Run the unchanged main() and record the exception it reports."""
    module = load(args.producer, "original_gnomad_cache")
    record = {}
    builtin_print = print

    def observed_print(*values, **options):
        if options.get("file") is sys.stderr and values and str(values[0]).startswith("ERROR: "):
            kind, value, _ = sys.exc_info()
            if kind is not None:
                record.update(category=kind.__name__, message=str(value))
        return builtin_print(*values, **options)

    module.print = observed_print
    sys.argv = [str(args.producer), *args.arguments]
    try:
        status = module.main()
    except SystemExit:
        raise
    except BaseException as exc:
        record.update(category=type(exc).__name__, message=str(exc), uncaught=True)
        write_json(args.exception, record)
        raise
    write_json(args.exception, record)
    raise SystemExit(status)


# ---------------------------------------------------------------- fixtures

def family(arm, script, assembly="hg38", variant=None, name="selected.chrX.vcf.bgz"):
    import pysam
    members = script.family_paths(arm / "gnomad" / name)
    lengths = {chromosome: 100000 + index for index, chromosome in enumerate(script.CHROMOSOMES)}
    lengths["chr1"] = script.CHR1_LENGTH[assembly]
    for chromosome, path in members.items():
        text = "##fileformat=VCFv4.2\n" + "".join(
            f"##contig=<ID={name},length={length + (variant == 'wrong_length' and name == chromosome == 'chr7')}>\n"
            for name, length in lengths.items())
        for field in script.REQUIRED_INFO:
            if variant == "missing_info" and chromosome == "chr7" and field in ("AF_joint", "AN_joint"):
                continue
            kind = "Float" if field.startswith("AF_") else "Integer"
            if variant == "wrong_type" and chromosome == "chr7" and field == "AF_joint":
                kind = "String"
            text += f'##INFO=<ID={field},Number=1,Type={kind},Description="fixture">\n'
        text += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        if variant != "empty_contigs":
            text += f"{chromosome}\t{1 if variant != 'beyond_length' else lengths[chromosome] + 5}\t.\tA\tC\t.\tPASS\t.\n"
        plain = arm / "gnomad/plain.vcf"
        plain.parent.mkdir(parents=True, exist_ok=True)
        plain.write_text(text)
        pysam.tabix_compress(str(plain), str(path), force=True)
        pysam.tabix_index(str(path), preset="vcf", force=True, csi=variant == "csi")
        plain.unlink()
    if variant in ("missing_member", "missing_index", "corrupt_index", "empty_member"):
        target = members["chr7"] if variant in ("missing_member", "empty_member") else Path(str(members["chr7"]) + ".tbi")
        if variant == "corrupt_index":
            target.write_bytes(b"broken")
        elif variant == "empty_member":
            target.write_bytes(b"")
        else:
            target.unlink()
    return members


def reference(arm, rows):
    path = arm / "reference.fa"
    path.write_text(">chr1\nA\n")
    if rows is not None:
        Path(str(path) + ".fai").write_text("".join(f"{row}\n" for row in rows))
    return path


def catalogue():
    chr_x = "{arm}/gnomad/selected.chrX.vcf.bgz"
    validate = ["validate", "--chr-x", chr_x, "--assembly", "GRCh38"]
    with_reference = [*validate, "--reference", "{arm}/reference.fa"]
    cases = {
        "argument_unknown_command": {"args": ["inspect"]},
        "argument_members_without_path": {"args": ["members"]},
        "argument_validate_without_assembly": {"args": ["validate", "--chr-x", chr_x]},
        "members_bgz": {"args": ["members", "--chr-x", chr_x]},
        "members_gz_relative": {"args": ["members", "--chr-x", "cache/./a.b.chrX.vcf.gz"]},
        "members_other_template": {"args": ["members", "--chr-x", "{arm}/gnomad/selected.chr1.vcf.bgz"]},
        "members_partial_folder": {"args": ["members", "--chr-x", "{arm}/x.PARTIAL/selected.chrX.vcf.bgz"]},
        "members_building_folder": {"args": ["members", "--chr-x", "{arm}/.building-1/selected.chrX.vcf.bgz"]},
        "default_hg19": {"args": ["default-path", "--cache-dir", "{arm}/g//x/", "--assembly", "hg19"]},
        "default_grch38_relative": {"args": ["default-path", "--cache-dir", "rel/../g", "--assembly", "GRCh38"]},
        "default_unknown_assembly": {"args": ["default-path", "--cache-dir", "/g", "--assembly", "hg39"]},
        "validate_hg38": {"family": {}, "args": validate},
        "validate_hg19": {"family": {"assembly": "hg19"}, "args": [*validate[:4], "hg19"]},
        "validate_csi": {"family": {"variant": "csi"}, "args": validate},
        "validate_empty_contigs": {"family": {"variant": "empty_contigs"}, "args": validate},
        "validate_wrong_assembly": {"family": {}, "args": [*validate[:4], "GRCh37"]},
        "validate_unknown_assembly": {"family": {}, "args": [*validate[:4], "T2T"]},
        "validate_missing_member": {"family": {"variant": "missing_member"}, "args": validate},
        "validate_empty_member": {"family": {"variant": "empty_member"}, "args": validate},
        "validate_missing_index": {"family": {"variant": "missing_index"}, "args": validate},
        "validate_corrupt_index": {"family": {"variant": "corrupt_index"}, "args": validate},
        "validate_missing_info": {"family": {"variant": "missing_info"}, "args": validate},
        "validate_wrong_type": {"family": {"variant": "wrong_type"}, "args": validate},
        "validate_wrong_length": {"family": {"variant": "wrong_length"}, "args": validate},
        "validate_record_beyond_length": {"family": {"variant": "beyond_length"}, "args": validate},
        "validate_other_template": {"family": {}, "args": ["validate", "--chr-x", "{arm}/gnomad/selected.chr2.vcf.bgz",
                                                           "--assembly", "hg38"]},
        "validate_reference_matching": {"family": {}, "reference": ["chr1\t248956422\t6\t1\t2", "2\t100001\t1\t1\t2",
                                                                   "chrM\t16569\t1\t1\t2"],
                                        "args": with_reference},
        "validate_reference_overrides_header": {"family": {}, "reference": ["chr1\t248956422\t6\t1\t2",
                                                                           "chr3\t5\t1\t1\t2"],
                                                "args": with_reference},
        "validate_reference_other_assembly": {"family": {}, "reference": ["1\t249250621\t6\t1\t2"],
                                              "args": with_reference},
        "validate_reference_malformed": {"family": {}, "reference": ["chr1\t248956422", "chr2"],
                                         "args": with_reference},
        "validate_reference_conflicting": {"family": {}, "reference": ["chr1\t248956422\t6\t1\t2",
                                                                      "1\t248956423\t6\t1\t2"],
                                           "args": with_reference},
        "validate_reference_crlf": {"family": {}, "reference": ["chr1\t248956422\t6\t1\t2\r"],
                                    "args": with_reference},
        "validate_reference_without_index": {"family": {}, "reference": None, "args": with_reference},
    }
    return cases


# ---------------------------------------------------------------- runs

def normalize(text, arm):
    if not isinstance(text, str):
        return text
    for spelling in sorted({str(arm), os.path.realpath(arm)}, key=len, reverse=True):
        text = text.replace(spelling, "<ARM>")
    return text


def tree(root):
    state = {}
    for folder, directories, names in os.walk(root, followlinks=False):
        for name in sorted(directories + names):
            path = Path(folder) / name
            relative = str(path.relative_to(root))
            if path.is_symlink():
                state[relative] = ["link", normalize(os.readlink(path), root)]
                continue
            value = path.stat()
            if path.is_dir():
                state[relative] = ["directory", oct(value.st_mode & 0o7777)]
                continue
            state[relative] = ["file", hashlib.sha256(path.read_bytes()).hexdigest(),
                               oct(value.st_mode & 0o7777), value.st_mtime_ns == FIXED_NS]
    return state


def run(command, folder, label, environment, cwd):
    with (folder / (label + ".stdout")).open("wb") as out, (folder / (label + ".stderr")).open("wb") as err:
        process = subprocess.run(list(map(str, command)), cwd=cwd, stdout=out, stderr=err,
                                 stdin=subprocess.DEVNULL, env=environment, timeout=600)
    return process.returncode


def outcome(folder, label, status, arm):
    stdout = (folder / (label + ".stdout")).read_bytes().decode("utf-8", "replace")
    stderr = (folder / (label + ".stderr")).read_bytes().decode("utf-8", "replace")
    result = {"exit": status, "stdout": normalize(stdout, arm), "tree": tree(arm)}
    if label == "python":
        recorded = folder / "exception.json"
        value = json.loads(recorded.read_text()) if recorded.exists() else {}
        if value:
            result.update(category=value.get("category"), message=value.get("message"))
    else:
        matched = re.search(r"(?:\A|\n)priva: ([A-Za-z_][A-Za-z_0-9]*): ((?:(?!\npriva: ).)*)\n\Z", stderr, re.S)
        if matched:
            result.update(category=matched[1], message=matched[2])
        elif stderr or status not in (0, 2):
            result.update(category=None, message=stderr)
    if "message" in result:
        result["message"] = normalize(result["message"], arm)
    return result


def compare(python, native, before):
    if python["exit"] == 2 or native["exit"] == 2:
        return python["exit"] == native["exit"] == 2 and python["tree"] == native["tree"] == before
    left, right = dict(python), dict(native)
    left.pop("uncaught", None)
    return left == right and right["tree"] == before


def check(args):
    producer = args.repo / "scripts/gnomad_cache.py"
    if hashlib.sha256(producer.read_bytes()).hexdigest() != PRODUCER_SHA256:
        raise SystemExit("gnomad_cache.py differs from the qualified source; update the pin")
    script = load(producer, "qualified_gnomad_cache")
    cases = catalogue()
    unknown = set(args.only or ()) - cases.keys()
    if unknown:
        raise SystemExit("Unknown cases: " + ", ".join(sorted(unknown)))
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
            arm.mkdir()
            if "family" in case:
                family(arm, script, **case["family"])
            if "reference" in case:
                reference(arm, case["reference"])
            for path in arm.rglob("*"):
                if path.is_file():
                    os.utime(path, ns=(FIXED_NS, FIXED_NS))
            befores[label] = tree(arm)
            environment = {key: value for key, value in os.environ.items() if not key.startswith("PRIVA_")}
            environment["PYTHONDONTWRITEBYTECODE"] = "1"
            arguments = [value.replace("{arm}", str(arm)) for value in case["args"]]
            if label == "python":
                command = [sys.executable, Path(__file__).resolve(), "oracle", "--producer", producer,
                           "--exception", folder / "exception.json", "--", *arguments]
            else:
                command = [args.binary, "gnomad-cache", *arguments]
            status = run(command, folder, label, environment, arm)
            results[label] = outcome(folder, label, status, arm)
        passed = befores["python"] == befores["native"] and \
            compare(results["python"], results["native"], befores["native"])
        rows[name] = {"passed": passed, "python": results["python"], "native": results["native"]}
        write_json(folder / "comparison.json", rows[name])
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
