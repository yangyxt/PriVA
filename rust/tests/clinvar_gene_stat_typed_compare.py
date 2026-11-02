#!/usr/bin/env python3
"""Compare schema-5 v2 with existing source-produced pickles and unchanged consumers.

  python clinvar_gene_stat_typed_compare.py check --binary NEW --baseline-binary OLD \
      --source-probes EXISTING_FIXTURES --scripts SCRIPTS --work-root NEW_DIRECTORY

No VCF producer is rerun. Existing source-produced pickles stay read-only. Additional
synthetic caches exercise scalar bits, consumer branches and rejected pickle forms.
"""
from __future__ import annotations
import argparse
import copy
import gzip
import importlib.util
import json
from pathlib import Path
import pickle
import struct
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import acmg_step_oracle_common as oracle
import cache_values_compare as values

CACHE = "clinvar_gene_stat"
CONSUMER_SHA256 = "6852e46ac99f10e39c1ffb428028d446da651fa1a8e39fc46724f1610f5d198e"
SET_NAMES = ("clinvar_variants", "clinvar_pathogenic", "clinvar_not_pathogenic",
             "putative_nmd_variants", "small_aachange", "large_aachange", "no_aachange")
SENTINEL = b"existing archive must survive rejected input\n"


def run(command, directory, name):
    with (directory/(name+".stdout.log")).open("x") as out, (directory/(name+".stderr.log")).open("x") as err:
        return subprocess.run(list(map(str, command)), stdout=out, stderr=err).returncode


def converted(args, source, directory, name, version, baseline=False):
    archive = directory/(name+".rkyv")
    code = run([args.binary, "cache-convert", CACHE, "--pickle", source, "--output", archive], directory, name)
    result = {"exit": code, "passed": code == 0}
    if code:
        return result, archive
    header = struct.unpack("<8sIII", archive.read_bytes()[:20])
    result["header"] = [header[0].decode(), *header[1:]]
    result["passed"] &= header == (b"PRIVARKV", 1, 5, version)
    result["comparison"] = values.compare(SimpleNamespace(binary=args.binary, cache=CACHE,
        pickle=source, archive=archive, report=directory/(name+".values.json"), mismatches=directory/(name+".mismatches.jsonl")))
    result["passed"] &= result["comparison"]["passed"]
    if baseline:
        old = directory/(name+".old.rkyv")
        status = run([args.baseline_binary, "cache-convert", CACHE, "--pickle", source, "--output", old], directory, name+".old")
        result["v1_bytes_equal"] = status == 0 and old.read_bytes() == archive.read_bytes()
        result["passed"] &= result["v1_bytes_equal"]
    return result, archive


def consumers(args, source, archive, data, directory, module):
    directory.mkdir()
    keys = list(data)+["ABSENT"]
    frame = pd.DataFrame({"Gene": pd.Series(keys, dtype=object)})
    oracle.write_frame(frame, directory/"input.frame")
    expected = [module.analyze_bp1_pp2(data.get(key, {"ensembl_id": key})) for key in keys]
    expected = pd.DataFrame({"ensembl_id": [item[0] for item in expected],
                             "pp2": [bool(item[1][0]) for item in expected],
                             "bp1": [bool(item[1][1]) for item in expected]})
    oracle.write_frame(expected, directory/"expected.frame")
    status, stderr = oracle.run_rust(args.binary, "acmg_pm1_pp2_bp1_missense", "analyze_bp1_pp2", directory,
                                   {"clinvar_gene_stat_pkl": str(archive)})
    (directory/"native.stderr.log").write_text(stderr)
    diff = [] if status else oracle.compare_frames(directory/"expected.frame", directory/"observed.frame", show=1000000)
    result = {"exit": status, "differences": diff, "passed": status == 0 and not diff}
    combined = directory/"combined"
    combined.mkdir()
    table = pd.DataFrame({"Gene": keys*2+[None],
        "Consequence": ["missense_variant"]*len(keys)+["synonymous_variant"]*len(keys)+[""],
        "Gene_avg_AM_score": [0.0]*(len(keys)*2)+[float("nan")]})
    domains = combined/"domains.tsv"
    domains.write_text("domain\tis_more_tolerant\nOTHER:Pfam:X\tFalse\n")
    pp2, bp1 = module.PP2_BP1_criteria(table.copy(), data, str(domains), threads=1)
    oracle.write_frame(pd.DataFrame({"PP2": pp2, "BP1": bp1}), combined/"expected.frame")
    oracle.write_frame(table, combined/"input.frame")
    code, stderr = oracle.run_rust(args.binary, "acmg_pm1_pp2_bp1_missense", "PP2_BP1_criteria", combined,
        {"clinvar_gene_stat_pkl": str(archive), "am_intol_domains_tsv": str(domains)})
    (combined/"native.stderr.log").write_text(stderr)
    diff = [] if code else oracle.compare_frames(combined/"expected.frame", combined/"observed.frame", show=1000000)
    result["combined"] = {"exit": code, "differences": diff, "passed": code == 0 and not diff}
    result["passed"] &= result["combined"]["passed"]
    return result


def legacy_none_boundary(args, module, data, archive):
    """Retain the existing Python/native input-column limitation separately."""
    folder = args.work_root/"legacy_object_none"
    folder.mkdir()
    keys = list(data)
    table = pd.DataFrame({"Gene": [keys[0], None], "Consequence": ["missense_variant", None],
                          "Gene_avg_AM_score": [0.0, float("nan")]})
    domains = folder/"domains.tsv"
    domains.write_text("domain\tis_more_tolerant\nOTHER:Pfam:X\tFalse\n")
    module.PP2_BP1_criteria(table.copy(), data, str(domains), threads=1)
    statuses = {}
    for name, binary in [("baseline", args.baseline_binary), ("updated", args.binary)]:
        run_root = folder/name
        run_root.mkdir()
        oracle.write_frame(table, run_root/"input.frame")
        code, stderr = oracle.run_rust(binary, "acmg_pm1_pp2_bp1_missense", "PP2_BP1_criteria", run_root,
            {"clinvar_gene_stat_pkl": str(archive), "am_intol_domains_tsv": str(domains)})
        (run_root/"stderr.log").write_text(stderr)
        statuses[name] = {"exit": code, "error": stderr}
    return {"python_success": True, **statuses,
            "unchanged_existing_limitation": all(row["exit"] != 0 for row in statuses.values())
                and statuses["baseline"]["error"] == statuses["updated"]["error"]}


def record(key, symbol, index):
    sets = {name: set() for name in SET_NAMES}
    variants = [f"v{number}" for number in range(30)]
    sets["clinvar_variants"] = set(variants)
    sets["clinvar_pathogenic"] = set(variants[:[0, 2, 8, 20][index % 4]])
    sets["clinvar_not_pathogenic"] = set(variants)-sets["clinvar_pathogenic"]
    sets["small_aachange"] = set(variants[::2])
    sets["large_aachange"] = set(variants[1::2])
    sets["putative_nmd_variants"] = set(variants[:index % 30])
    return {"gene_symbol": symbol, "ensembl_id": key, **sets,
            "transcripts": {"ENST_Z": copy.deepcopy(sets), "ENST_A": copy.deepcopy(sets)},
            "domains": {key+":Pfam:Z": copy.deepcopy(sets), key+":Pfam:A": copy.deepcopy(sets)}}


def check(args):
    args.work_root.mkdir(parents=True)
    source = args.scripts/"acmg_pm1_pp2_bp1_missense.py"
    if values.sha256(source) != CONSUMER_SHA256:
        raise RuntimeError("consumer source differs from pinned identity")
    sys.path.insert(0, str(args.scripts))
    import acmg_pm1_pp2_bp1_missense as module
    report = {"inputs": {"binary": values.sha256(args.binary), "baseline_binary": values.sha256(args.baseline_binary),
                         "consumer": CONSUMER_SHA256}, "cases": {}, "failures": {}}
    for name in ["text", "all_missing", "mixed_missing_text", "integers", "floats", "mixed_numeric_text",
                 "nullable_integer", "nullable_float", "boolean_text"]:
        folder = args.work_root/name
        folder.mkdir()
        original = args.source_probes/name/"reference.pkl.gz"
        with gzip.open(original, "rb") as stream:
            data = pickle.load(stream)
        v1 = all(type(entry["gene_symbol"]) is str for entry in data.values())
        row, archive = converted(args, original, folder, "convert", 1 if v1 else 2, baseline=v1)
        if row["passed"]:
            row["consumer"] = consumers(args, original, archive, data, folder/"consumer", module)
            row["passed"] &= row["consumer"]["passed"]
        report["cases"][name] = row
        if name == "text":
            report["existing_consumer_boundary"] = legacy_none_boundary(args, module, data, archive)

    symbols = ["TEXT", pd.NA, np.int32(-2147483648), np.int32(2147483647), -0.0, float("inf"), float("-inf")]
    symbols += [struct.unpack(">d", bytes.fromhex(bits))[0] for bits in ["7ff8000000000042", "7ff0000000000042"]]
    symbols += [np.frombuffer(struct.pack("<I", bits), dtype="<f4")[0]
                for bits in [0, 1, 0x80000000, 0x7f800000, 0xff800000, 0x7fc00042, 0x7f800042]]
    data = {f"GENE_{index}": record(f"GENE_{index}", symbol, index) for index, symbol in enumerate(symbols)}
    for name, typed in [("scalar_bits", True), ("text_control", False)]:
        folder = args.work_root/name
        folder.mkdir()
        current = copy.deepcopy(data)
        if not typed:
            for gene in current.values():
                gene["gene_symbol"] = "TEXT"
        source = folder/"source.pkl"
        source.write_bytes(pickle.dumps(current, protocol=5))
        row, archive = converted(args, source, folder, "convert", 2 if typed else 1, baseline=not typed)
        if row["passed"]:
            row["consumer"] = consumers(args, source, archive, current, folder/"consumer", module)
            row["passed"] &= row["consumer"]["passed"]
        report["cases"][name] = row
    # Unapproved globals are never unpickled by this harness. Only native decoding
    # sees these payloads; the destination must remain intact on rejection.
    invalid_symbols = [None, True, 1, np.int64(1), np.float64(1), np.uint32(1), np.array([1]), [], {}, set()]
    for index, value in enumerate(invalid_symbols):
        folder = args.work_root/f"invalid_{index}"
        folder.mkdir()
        source, archive = folder/"bad.pkl", folder/"destination.rkyv"
        source.write_bytes(pickle.dumps({"G": record("G", value, 1)}, protocol=5))
        archive.write_bytes(SENTINEL)
        code = run([args.binary, "cache-convert", CACHE, "--pickle", source, "--output", archive], folder, "convert")
        report["failures"][str(index)] = {"exit": code, "passed": code != 0 and archive.read_bytes() == SENTINEL}
    for index, value in enumerate([pd.NA, np.int32(1), np.float32(-0.0)]):
        folder = args.work_root/f"default_reader_{index}"
        folder.mkdir()
        source = folder/"scalar.pkl"
        source.write_bytes(pickle.dumps(value, protocol=5))
        code = run([args.binary, "acmg-step", "inspect", "pickle", source], folder, "inspect")
        report["failures"]["default_"+str(index)] = {"exit": code, "passed": code != 0}
    report["passed"] = all(row["passed"] for row in [*report["cases"].values(), *report["failures"].values()])
    (args.work_root/"comparison.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"passed": report["passed"], "cases": len(report["cases"]), "rejection_controls": len(report["failures"]),
        "failed": [name for name, row in {**report["cases"], **report["failures"]}.items() if not row["passed"]]}))
    if not report["passed"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    sub = commands.add_parser("check")
    for name in ("binary", "baseline-binary", "source-probes", "scripts", "work-root"):
        sub.add_argument("--"+name, type=Path, required=True)
    sub.set_defaults(function=check)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
