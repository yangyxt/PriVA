#!/usr/bin/env python3
"""Compare the native PM1 builder with an unchanged prepare_pm1_regions.py.

python pm1_builder_compare.py compare --binary PRIVA --repo FROZEN_SOURCE \
    --work-root EMPTY_DIRECTORY [--hg38-das-pickle TRUSTED_PICKLE]

The Python producer is copied byte-for-byte into each disposable reference tree
so its normal bundled-path resolution is exercised. All successful outputs use
the independent complete-value cache comparator: scalar types, float bits,
metadata order and every set member are exact. Only the cache format's declared
lookup-key/set order is canonicalized. No pickle writer enters production.

Failure checks compare exit status, exception category and output preservation;
traceback/log text is retained but is not the native CLI's text contract. Timing
and RSS include each executable's startup and are single measurements, not a
production speedup claim. The hg19 complete bundled-table case uses an explicitly
synthetic DAS set; only --hg38-das-pickle enables a real full PM1 input profile.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import pickle
import random
import re
import shutil
import subprocess
import sys
import time

MAPPING_HEADER = "Gene stable ID\tTranscript stable ID\tHGNC symbol\tUniProtKB/Swiss-Prot ID\tRefSeq mRNA ID\n"
MAPPING = MAPPING_HEADER + (
    "G1\tENST2\tGENE\tP1\tNM_1.2\n"
    "G1\tENST1\tGENE\tP1\tNM_1.3\n"
    "G1\tENST2\tGENE\tP1\tNM_1\n"
    "G2\tENST3\tOTHER\tP2\tNM_2\n"
)
HC_HEADER = "type\tchrom\tstart\tend\tisoform\taa_start\taa_end\n"
HOTSPOTS = HC_HEADER + (
    "hotspot\t1\t1\t9\tNM_1.4\t1\t4\n"
    "hotspot\t1\t3\t9\tNM_1\t3\t5\n"
    "hotspot\t1\t1\t9\tNM_2\t8\t6\n"
    "coldspot\t1\t1\t9\tNM_1\tbad\tbad\n"
    "hotspot\t1\t1\t9\tUNKNOWN\t1\t4\n"
)
RMC_HEADER = "transcript\tstart_aa_raw\tend_aa_raw\toe\tp\n"
RMC = RMC_HEADER + (
    "ENST2\tArg1\tGly4\t0.3\t0.001\n"
    "ENST2\tGly3\tMet7\t0.2\t0.002\n"
    "ENST3\tMet9\tGly7\t0.1\t0.001\n"
    "ENST1\tArg8\tMet9\t0.4\t0.001\n"
)


def sha(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path: Path, data: object) -> None:
    with path.open("x") as output:
        json.dump(data, output, indent=2, sort_keys=True, allow_nan=False)
        output.write("\n")


def measured(command: list[str], directory: Path, name: str) -> dict:
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    start = time.monotonic()
    with (directory / f"{name}.stdout").open("x") as output, (directory / f"{name}.stderr").open("x") as errors:
        # wait4 reports this child's peak RSS, unlike the cumulative
        # RUSAGE_CHILDREN maximum. Do not depend on a separately installed time.
        process = subprocess.Popen(command, cwd=directory, env=env, stdout=output, stderr=errors)
        try:
            while True:
                child, status, usage = os.wait4(process.pid, os.WNOHANG)
                if child:
                    process.returncode = os.waitstatus_to_exitcode(status)
                    break
                if time.monotonic() - start > 300:
                    raise TimeoutError("comparison child exceeded 300 seconds")
                time.sleep(0.01)
        except BaseException:
            process.kill()
            _, status, _ = os.wait4(process.pid, 0)
            process.returncode = os.waitstatus_to_exitcode(status)
            raise
    timing = dict(wall_s=time.monotonic() - start, peak_rss_kib=usage.ru_maxrss)
    save(directory / f"{name}.time.json", timing)
    return dict(command=command, returncode=process.returncode, **timing)


def fixtures() -> list[dict]:
    cases = []

    def add(name: str, **changes) -> None:
        cases.append(dict(name=name, **changes))

    add("mapping-multiplicity-overlap-and-empty-range")
    add("all-empty", mapping="", hotspots="", rmc="", das=set())
    add("header-only", mapping=MAPPING_HEADER, hotspots=HC_HEADER, rmc=RMC_HEADER)
    add("missing-rmc-columns-with-no-records", rmc="unrecognized\n")
    add("empty-enst-and-nm-version", mapping=MAPPING + "G\t \tS\tP\tNM_1.19\n")
    add("non-nm-and-short-mapping-ignored", mapping=MAPPING + "G\tENST99\tS\tP\tXM_1\nshort\n")
    add("mapping-unconditional-header-skip", mapping=MAPPING.removeprefix(MAPPING_HEADER))
    add("hc-unconditional-header-skip", hotspots=HOTSPOTS.removeprefix(HC_HEADER))
    add("hotspot-whitespace-negative-and-underscore-int",
        hotspots=HC_HEADER + " hotspot \t1\t1\t9\t NM_1.2 \t -2 \t+0_2\n")
    add("unicode-mapping-and-integer",
        mapping=MAPPING + "G\t\u2003ENSTユ\u2003\tS\tP\t\u2003NM_3.2\u2003\n",
        hotspots=HC_HEADER + "hotspot\t1\t1\t9\tNM_3\t１２\t١٤\n")
    add("short-hc-and-nonhotspot-ignored", hotspots=HOTSPOTS + "short\nother\t1\t1\t9\tNM_1\tx\ty\n")
    add("unmapped-hotspot-still-parses", hotspots=HOTSPOTS + "hotspot\t1\t1\t9\tUNKNOWN\tx\t2\n", failure="ValueError")
    add("hc-failure-before-missing-rmc", hotspots=HC_HEADER + "hotspot\t1\t1\t9\tNM_1\tx\t2\n", missing=["rmc"], failure="ValueError")
    add("mapping-decode-failure-before-missing-das", mapping=b"h\n\xff\n", missing=["das"], failure="UnicodeDecodeError")
    add("hc-negative-empty-range", hotspots=HC_HEADER + "hotspot\t1\t1\t9\tNM_1\t-1\t-3\n")
    add("hc-invalid-float-int", hotspots=HC_HEADER + "hotspot\t1\t1\t9\tNM_1\t1.0\t2\n", failure="ValueError")
    add("hc-invalid-superscript-int", hotspots=HC_HEADER + "hotspot\t1\t1\t9\tNM_1\t²\t2\n", failure="ValueError")
    for label, ending in (("crlf", "\r\n"), ("bare-cr", "\r"), ("no-final-newline", "")):
        convert = (lambda value: value.rstrip("\n")) if not ending else (lambda value: value.replace("\n", ending))
        add(label, mapping=convert(MAPPING), hotspots=convert(HOTSPOTS), rmc=convert(RMC))
    add("rmc-first-unicode-decimal-run", rmc=RMC_HEADER + " ENST1 \tArg９７then100\t-Arg١٠٠.5\t٠.٣\t٠.٠٠١\n")
    add("rmc-empty-transcript", rmc=RMC_HEADER + " \tArg1\tArg2\t0.1\t0.001\n")
    add("rmc-duplicate-column-last-wins", rmc=RMC_HEADER.rstrip("\n") + "\tp\nENST1\tA1\tA2\t0.2\tBAD\t0.001\n")
    add("rmc-columns-reordered", rmc="p\tend_aa_raw\tstart_aa_raw\ttranscript\toe\n0.001\tA2\tA1\tENST1\t0.2\n")
    add("rmc-empty-range", rmc=RMC_HEADER + "ENST1\tArg9\tArg2\t0.1\t0.001\n")
    add("rmc-strict-oe-boundary", rmc=RMC_HEADER + "ENST1\tA1\tA2\t0.4\t0.001\n")
    add("rmc-strict-bh-boundary", rmc=RMC_HEADER + "ENST1\tA1\tA2\t0.1\t0.05\n")
    add("rmc-oe-nan-inf", rmc=RMC_HEADER + "ENST1\tA1\tA2\tnan\t0.001\nENST2\tA1\tA2\t-inf\t0.002\nENST3\tA1\tA2\tinf\t0.001\n")
    for name, value in (("nan", "nan"), ("positive-inf", "inf"), ("negative-inf", "-inf"), ("signed-zero", "-0.0"), ("large", "2.0")):
        add("threshold-" + name, oe=value, p=value)
    add("singleton-negative-zero", rmc=RMC_HEADER + "ENST1\tA1\tA2\t0.1\t-0.0\n", p="0.0")
    for name, value in (("nan", "nan"), ("infinity", "inf"), ("negative", "-0.01"), ("above-one", "1.01")):
        add("invalid-p-" + name, rmc=RMC_HEADER + f"ENST1\tA1\tA2\t0.1\t{value}\n", failure="ValueError")
    for name, row, category in (
        ("missing-column", "name\nENST1\n", "KeyError"),
        ("short-row", RMC_HEADER + "ENST1\tA1\n", "IndexError"),
        ("blank-record", RMC_HEADER + "\n", "IndexError"),
        ("bad-aa", RMC_HEADER + "ENST1\tA\tA2\t0.1\t0.001\n", "ValueError"),
        ("superscript-aa", RMC_HEADER + "ENST1\tA²\tA2\t0.1\t0.001\n", "ValueError"),
        ("bad-oe-before-aa", RMC_HEADER + "ENST1\tA\tA2\tbad\t0.001\n", "ValueError"),
        ("bad-p", RMC_HEADER + "ENST1\tA1\tA2\t0.1\tbad\n", "ValueError"),
    ):
        add(name, rmc=row, failure=category)
    for name in ("mapping", "hotspots", "rmc"):
        add("missing-" + name, missing=[name], failure="FileNotFoundError")
    add("replace-existing-complete-output", prior=True)
    rng = random.Random(73219)
    for number in range(24):
        rows = []
        for _ in range(rng.randrange(1, 101)):
            start, end = rng.randrange(0, 50), rng.randrange(0, 50)
            rows.append(f"ENST{rng.randrange(8)}\tA{start}\tB{end}\t{rng.choice([0.1,0.4,0.9])!r}\t{rng.choice([0.0,0.001,0.025,0.05,0.5,1.0])!r}\n")
        add(f"seeded-rmc-{number:02d}", rmc=RMC_HEADER + "".join(rows))
    return cases


def error_category(path: Path) -> str | None:
    value = path.read_text(errors="replace")
    categories = re.findall(r"(?:^|: )(ValueError|IndexError|KeyError|UnicodeDecodeError|FileNotFoundError|IsADirectoryError):", value, re.M)
    if categories:
        return categories[-1]
    if "No such file or directory" in value:
        return "FileNotFoundError"
    if "valid UTF-8" in value:
        return "UnicodeDecodeError"
    return None


def one_case(args: argparse.Namespace, spec: dict) -> dict:
    root = args.work_root / spec["name"]
    root.mkdir()
    assembly = spec.get("assembly", "hg38")
    reference = root / "reference"
    producer = reference / "scripts/prepare_pm1_regions.py"
    producer.parent.mkdir(parents=True)
    shutil.copy2(args.repo / "scripts/prepare_pm1_regions.py", producer)
    paths = dict(mapping=root / "mapping.tsv", hotspots=reference / f"data/HCSeeker/HC_spots.{assembly}.tsv",
                 rmc=reference / f"data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.{assembly}.tsv",
                 das=root / "das.pkl", native_das=root / "das.rkyv")
    for name, default in (("mapping", MAPPING), ("hotspots", HOTSPOTS), ("rmc", RMC)):
        target = paths[name]
        target.parent.mkdir(parents=True, exist_ok=True)
        value = spec.get(name, default)
        if isinstance(value, Path):
            shutil.copy2(value, target)
        else:
            target.write_bytes(value if isinstance(value, bytes) else value.encode())
    das = spec.get("das", {"ENSG1:Pfam:PF1", "ENSG2:SMART:SM1", "ENSG2:Pfam:PF2"})
    if isinstance(das, Path):
        shutil.copy2(das, paths["das"])
    else:
        with paths["das"].open("xb") as output:
            pickle.dump(das, output, protocol=pickle.HIGHEST_PROTOCOL)
    conversion = measured([str(args.binary), "cache-convert", "all_intolerant_domains", "--pickle", str(paths["das"]), "--output", str(paths["native_das"])], root, "das-convert")
    if conversion["returncode"]:
        raise RuntimeError("fixture DAS conversion failed: " + spec["name"])
    for name in spec.get("missing", []):
        paths[name].unlink()
        if name == "das":
            paths["native_das"].unlink()
    config = root / "config.yaml"
    config.write_text("all_intolerant_domains: " + json.dumps(str(paths["das"])) + "\n")
    output_py, output_rs = root / "reference.pkl", root / "native.rkyv"
    prior = b"preserve this prior output until all sources finish\n"
    if spec.get("failure") or spec.get("prior"):
        output_py.write_bytes(prior)
        output_rs.write_bytes(prior)
    tracked = [producer, config, *[p for p in paths.values() if p.exists()]]
    before = {str(path): sha(path) for path in tracked}
    reference_command = [sys.executable, str(producer), "--config", str(config), "--output", str(output_py),
                         "--assembly", assembly, "--skip_das", "--nm_to_enst_tsv", str(paths["mapping"])]
    native_command = [str(args.binary), "pm1-regions", "build", "--assembly", assembly, "--output", str(output_rs),
                      "--nm-to-enst", str(paths["mapping"]), "--intolerant-domains", str(paths["native_das"]),
                      "--hcseeker", str(paths["hotspots"]), "--rmc", str(paths["rmc"])]
    for key, python_name, native_name in (("oe", "rmc_oe_threshold", "rmc-oe-threshold"), ("p", "rmc_p_bh_threshold", "rmc-p-bh-threshold")):
        if key in spec:
            reference_command.append(f"--{python_name}={spec[key]}")
            native_command.append(f"--{native_name}={spec[key]}")
    python = measured(reference_command, root, "python")
    native = measured(native_command, root, "native")
    mismatch = []
    if python["returncode"] != native["returncode"]:
        mismatch.append(dict(kind="exit_status", expected=python["returncode"], actual=native["returncode"]))
    result = dict(name=spec["name"], scope=spec.get("scope", "synthetic-edge-case"), python=python, native=native, das_conversion=conversion)
    if spec.get("failure"):
        classes = dict(python=error_category(root / "python.stderr"), native=error_category(root / "native.stderr"))
        result["failure_categories"] = classes
        if classes != {"python": spec["failure"], "native": spec["failure"]}:
            mismatch.append(dict(kind="failure_category", expected=spec["failure"], actual=classes))
        if not python["returncode"] or not native["returncode"]:
            mismatch.append(dict(kind="expected_failure", expected="nonzero", actual=[python["returncode"], native["returncode"]]))
        for name, output in (("python", output_py), ("native", output_rs)):
            if output.read_bytes() != prior:
                mismatch.append(dict(kind="prior_output_changed", arm=name))
    elif python["returncode"] or native["returncode"]:
        mismatch.append(dict(kind="unexpected_failure", actual=[python["returncode"], native["returncode"]]))
    else:
        comparison = measured([sys.executable, str(args.repo / "rust/tests/cache_values_compare.py"), "compare", "--binary", str(args.binary),
                               "--cache", "pm1_regions_pkl", "--pickle", str(output_py), "--archive", str(output_rs),
                               "--report", str(root / "values.json"), "--mismatches", str(root / "values.mismatches.jsonl")], root, "compare")
        result["comparison"] = json.loads((root / "values.json").read_text()) if (root / "values.json").exists() else comparison
        if comparison["returncode"]:
            mismatch.append(dict(kind="all_value_comparison", report=str(root / "values.json")))
    after = {str(path): sha(path) for path in tracked}
    if after != before:
        mismatch.append(dict(kind="input_mutation", before=before, after=after))
    result.update(inputs=before, inputs_unchanged=before == after, mismatches=mismatch, passed=not mismatch)
    save(root / "case.json", result)
    return result


def compare(args: argparse.Namespace) -> None:
    args.work_root.mkdir(parents=True)
    cases = fixtures()
    for assembly, mapping_name in (("hg19", "ensembl_GRCh37_mart116_nm_to_enst.tsv"), ("hg38", "ensembl_GRCh38_mart116_ftp_candidate.tsv")):
        spec = dict(name="complete-bundled-tables-" + assembly, assembly=assembly,
                    mapping=args.repo / "data/hgnc" / mapping_name,
                    hotspots=args.repo / f"data/HCSeeker/HC_spots.{assembly}.tsv",
                    rmc=args.repo / f"data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.{assembly}.tsv",
                    scope="complete-bundled-tables-with-synthetic-DAS")
        if assembly == "hg38" and args.hg38_das_pickle:
            spec.update(das=args.hg38_das_pickle, scope="complete-hg38-prepared-PM1-inputs")
        cases.append(spec)
    pins = {str(path): sha(path) for path in (args.binary, Path(__file__), args.repo / "scripts/prepare_pm1_regions.py", args.repo / "rust/tests/cache_values_compare.py")}
    results = []
    with (args.work_root / "mismatches.jsonl").open("x") as ledger:
        for case in cases:
            result = one_case(args, case)
            results.append(result)
            for mismatch in result["mismatches"]:
                ledger.write(json.dumps(dict(case=case["name"], **mismatch)) + "\n")
                ledger.flush()
            print(json.dumps(dict(case=case["name"], passed=result["passed"], mismatches=len(result["mismatches"]))), flush=True)
    stable = pins == {path: sha(Path(path)) for path in pins}
    report = dict(cases=len(results), passed=stable and all(case["passed"] for case in results),
                  matching_successes=sum(not case["python"]["returncode"] and case["passed"] for case in results),
                  matching_failures=sum(bool(case["python"]["returncode"]) and case["passed"] for case in results),
                  mismatch_count=sum(len(case["mismatches"]) for case in results),
                  identities=pins, identities_unchanged=stable,
                  results=results, limits=["Not a DAS regeneration test", "hg19 uses synthetic DAS", "Exception/log text and formats differ by interface", "Residues must fit the frozen signed-64-bit cache schema", "No new installer caller switch"])
    save(args.work_root / "report.json", report)
    print(json.dumps({key: value for key, value in report.items() if key not in ("identities", "results")}), flush=True)
    if not report["passed"]:
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    command = sub.add_parser("compare")
    for name in ("binary", "repo", "work-root"):
        command.add_argument("--" + name, type=Path, required=True)
    command.add_argument("--hg38-das-pickle", type=Path)
    command.set_defaults(function=compare)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
