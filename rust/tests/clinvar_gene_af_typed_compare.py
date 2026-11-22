#!/usr/bin/env python3
"""Check typed gene caches and complete unchanged Python PM2 -> BS1 behavior.

  python clinvar_gene_af_typed_compare.py check --scripts SCRIPTS_DIR --binary PRIVA --work-root NEW_DIR

The independent Python decoder compares all archive values and ordering. The
consumer comparison includes every frame cell/type and all three BS1 route masks.
Only the gene cache representation changes; the Python consumer source is pinned.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import pickle
from pathlib import Path
import struct
import subprocess
import sys
import traceback
from types import SimpleNamespace

import numpy as np
import pandas as pd

import acmg_step_oracle_common as oracle
import cache_values_compare as cache_compare

PINNED = {
    "acmg_pm2_bs1_ba1_frequency.py": "ef2e70115ddfe7087b4d24002ce4af95c4f4b91fbf38a06e8a52275eecca4fb3",
    "acmg_variant_mechanism.py": "c7e8da80ca0ec47fb38a5fd1d965641c1791189cb0d25a651253f44508d0bf26",
}
SENTINEL = b"existing destination must survive failure\n"


def execute(command, folder, name):
    with (folder / (name + ".stdout.log")).open("x") as out, (folder / (name + ".stderr.log")).open("x") as err:
        return subprocess.run(list(map(str, command)), stdout=out, stderr=err).returncode


def convert(binary, cache, source, output, folder, name):
    return execute([binary, "cache-convert", cache, "--pickle", source, "--output", output], folder, name)


def values(binary, folder):
    folder.mkdir()
    source, archive = folder / "typed.pkl", folder / "typed.rkyv"
    data = {"Z": {"acc_patho_af": False, "af": True, "pos": 17, "chrom": "chr1",
                  "impact": -0.0, "hgvsp": struct.unpack(">d", bytes.fromhex("7ff8000000000042"))[0]},
            "NONE": None, "EMPTY": {}}
    source.write_bytes(pickle.dumps(data, protocol=4))
    status = convert(binary, "clinvar_patho_af_stat", source, archive, folder, "convert")
    result = {"convert_exit": status, "passed": status == 0}
    if status:
        return result
    result["comparison"] = cache_compare.compare(SimpleNamespace(binary=binary, cache="clinvar_patho_af_stat", pickle=source,
        archive=archive, report=folder/"values.json", mismatches=folder/"mismatches.jsonl"))
    header = struct.unpack("<8sIII", archive.read_bytes()[:20])
    result["container_version"], result["schema_id"], result["schema_version"] = header[1:]
    result["passed"] &= result["comparison"]["passed"] and header == (b"PRIVARKV", 1, 6, 2)
    before = archive.read_bytes()
    result["failures"] = []
    for number, bad in enumerate([None, [], {"x": 1}, (1,), set()]):
        source.write_bytes(pickle.dumps({"Z": {"af": bad}}, protocol=4))
        status = convert(binary, "clinvar_patho_af_stat", source, archive, folder, f"invalid_{number}")
        passed = status != 0 and archive.read_bytes() == before
        result["failures"].append({"exit": status, "preserved": archive.read_bytes() == before, "passed": passed})
        result["passed"] &= passed
    return result


def frame(genes):
    n = len(genes)
    def repeat(items):
        return [items[i % len(items)] for i in range(n)]
    return pd.DataFrame({
        "Gene": pd.Series(genes, dtype=object),
        "Feature": pd.Series(["ENST1"] * n, dtype=object),
        "EXON": pd.Series(["1/2"] * n, dtype=object),
        "chrom": pd.Series(repeat(["chr1", "chrX", "chr2", "chrY", "chrM"]), dtype=object),
        "variant_inheritance": pd.Series(repeat(["dominant", "x_linked_recessive", "recessive", "y_linked", "mitochondrial"]), dtype=object),
        "variant_penetrance": pd.Series(["complete"] * n, dtype=object),
        "CLNSIG": pd.Series(["Pathogenic"] * n, dtype=object),
        "CLNREVSTAT": pd.Series(["reviewed_by_expert_panel"] * n, dtype=object),
        "gnomAD_joint_AF": repeat([0.0, 0.00001, 0.2, np.nan, 0.01]),
        "gnomAD_joint_AF_max": repeat([0.0, 0.0001, 0.3, np.nan, 0.02]),
        "gnomAD_joint_AN": repeat([100000.0, 100.0, 0.0, 10000.0, 80000.0]),
        "gnomAD_joint_AN_max": repeat([50000.0, 1000.0, 0.0, 10000.0, 40000.0]),
        "gnomAD_nhomalt_max": repeat([0.0, 1.0, 5.0, 0.0, 100.0]),
        "gnomAD_nhomalt_XX": repeat([0.0, 0.0, 3.0, 0.0, 10.0]),
        "gnomAD_nhomalt_XY": repeat([0.0, 1.0, 4.0, 0.0, 20.0]),
        "gnomAD_joint_AF_XY": repeat([0.0, 0.001, 0.2, 0.0, 0.01]),
        "gnomAD_joint_AF_XX": repeat([0.0, 0.001, 0.2, 0.0, 0.01]),
        "gnomAD_joint_AN_XY": [10000.0] * n,
        "gnomAD_joint_AN_XX": [10000.0] * n,
    })


def consumers(args, module):
    specs = {
        "all_boolean": ([True, False, True, False, True], []),
        "boolean_integer": ([True, 1, False, 0, 2], []),
        "boolean_float": ([True, 0.02, False, -0.0, float("nan")], []),
        "boolean_defaults": ([True], ["EMPTY", "NONE", "ABSENT", "G0"]),
        "unmatched_boolean": ([0.1, 0, -0.0, 1, 0.2], []),
        "v1_numeric": ([0.1, 0, -0.0, 1, float("nan")], []),
        "text_af_failure": ([True, "bad"], []),
    }
    results = {}
    for name, (numbers, extras) in specs.items():
        folder = args.work_root / name
        folder.mkdir()
        genes = [f"G{i}" for i in range(len(numbers))] + extras
        cache = {f"G{i}": {"af": value} for i, value in enumerate(numbers)}
        cache.update(EMPTY={}, NONE=None)
        if name == "unmatched_boolean":
            cache["UNMATCHED"] = {"af": True}
        source, archive = folder/"gene.pkl", folder/"gene.rkyv"
        source.write_bytes(pickle.dumps(cache, protocol=4))
        exon_source, exon_archive = folder/"exon.pkl", folder/"exon.rkyv"
        exon_source.write_bytes(pickle.dumps({"ENST1":{"1/2": (0.01,0.02,0.03,0.0,(0.01,0.03))}}, protocol=4))
        codes = [convert(args.binary, "clinvar_patho_af_stat", source, archive, folder, "gene_convert"),
                 convert(args.binary, "clinvar_patho_exon_af_stat", exon_source, exon_archive, folder, "exon_convert")]
        if any(codes):
            results[name] = {"passed":False, "conversion_exits":codes}
            continue
        original = frame(genes)
        pm2_dir = folder / "pm2"
        pm2_dir.mkdir()
        oracle.write_frame(original, pm2_dir/"input.frame")
        python_frame = original.copy(deep=True)
        def pm2_frame(value):
            python_frame["__PM2"] = value
            return python_frame
        pm2 = oracle.capture(pm2_dir, lambda: module.PM2_criteria(python_frame, str(source), str(exon_source)), pm2_frame)
        (pm2_dir/"observed.frame").write_bytes(SENTINEL)
        status, error = oracle.run_rust(args.binary, "acmg_pm2_bs1_ba1_frequency", "PM2_criteria", pm2_dir,
            {"clinvar_patho_af_stat":str(archive), "clinvar_patho_exon_af_stat":str(exon_archive)})
        (pm2_dir/"native.stderr.log").write_text(error)
        expected_ok = (pm2_dir/"expected.frame").exists()
        result = {"python_success":expected_ok, "native_exit":status, "passed":expected_ok == (status == 0)}
        if not expected_ok:
            result["python_error"] = json.loads((pm2_dir/"expected_error.json").read_text())
        if expected_ok and status == 0:
            result["pm2_differences"] = oracle.compare_frames(pm2_dir/"expected.frame", pm2_dir/"observed.frame", show=1000000)
            result["mapped_dtype"] = str(python_frame["clinvar_patho_gene_max_af"].dtype)
            result["passed"] &= not result["pm2_differences"]
            for mode in ["pm2", "no_pm2"]:
                bs1_dir = folder / ("bs1_"+mode)
                bs1_dir.mkdir()
                # Feed each arm its own PM2 result, preserving mutation order.
                (bs1_dir/"input.frame").write_bytes((pm2_dir/"observed.frame").read_bytes())
                py_bs1 = python_frame.copy(deep=True)
                pm2_values = pm2 if mode == "pm2" else None
                def bs1_frame(value):
                    criteria, components = value
                    py_bs1["__BS1"] = criteria
                    for key, values in components.items():
                        py_bs1["__"+key] = values
                    return py_bs1
                oracle.capture(bs1_dir, lambda: module.BS1_criteria(py_bs1, pm2_criteria=pm2_values, return_frequency_components=True), bs1_frame)
                (bs1_dir/"observed.frame").write_bytes(SENTINEL)
                code, error = oracle.run_rust(args.binary, "acmg_pm2_bs1_ba1_frequency", "BS1_criteria", bs1_dir,
                    {"pm2_criteria":None if pm2_values is None else [int(v) for v in pm2_values]})
                (bs1_dir/"native.stderr.log").write_text(error)
                py_ok = (bs1_dir/"expected.frame").exists()
                phase = {"python_success":py_ok,"native_exit":code,"passed":py_ok == (code == 0)}
                if not py_ok:
                    phase["python_error"] = json.loads((bs1_dir/"expected_error.json").read_text())
                if py_ok and code == 0:
                    phase["differences"] = oracle.compare_frames(bs1_dir/"expected.frame", bs1_dir/"observed.frame", show=1000000)
                    phase["passed"] &= not phase["differences"]
                elif not py_ok:
                    phase["output_preserved"] = (bs1_dir/"observed.frame").read_bytes() == SENTINEL
                    phase["passed"] &= phase["output_preserved"]
                result[mode] = phase
                result["passed"] &= phase["passed"]
        elif not expected_ok:
            result["output_preserved"] = (pm2_dir/"observed.frame").read_bytes() == SENTINEL
            result["passed"] &= result["output_preserved"]
        results[name] = result
    return results


def check(args):
    args.work_root.mkdir(parents=True, exist_ok=False)
    for name, wanted in PINNED.items():
        if oracle.sha256_file(args.scripts/name) != wanted:
            raise RuntimeError("consumer source does not match its pinned identity")
    provenance = oracle.manifest([args.scripts/name for name in PINNED])
    sys.path.insert(0, str(args.scripts))
    import acmg_pm2_bs1_ba1_frequency as module
    report = {"provenance":provenance, "values":values(args.binary,args.work_root/"values")}
    report["consumers"] = consumers(args,module)
    report["passed"] = report["values"]["passed"] and all(v["passed"] for v in report["consumers"].values())
    (args.work_root/"comparison.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"passed":report["passed"],"consumer_cases":len(report["consumers"]),
                      "failed":[k for k,v in report["consumers"].items() if not v["passed"]]}))
    if not report["passed"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command",required=True)
    command = commands.add_parser("check")
    command.add_argument("--scripts",type=Path,required=True)
    command.add_argument("--binary",type=Path,required=True)
    command.add_argument("--work-root",type=Path,required=True)
    args = parser.parse_args()
    check(args)


if __name__ == "__main__":
    main()
