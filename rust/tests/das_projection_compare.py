#!/usr/bin/env python3
"""Compare DAS score arrays and common-missense summaries with the original procedures.

check --binary PRIVA --source-dir SOURCE --work-root NEW
Each generated VCF is collected twice: by the unchanged Python statistics script
(final pickle) and by the native builder (typed Parquet). The native projection
commands are compared with `collect_domain_data` plus the dictionary assignment
in `visualize_domain_distribution`, and with `filter_common_missense`. Order,
float64 bits, scalar types and removal flags are exact. The current producer
always writes residue dictionaries, so the flat-array fallback is out of reach
here; das_common_compare.py covers it.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib
import importlib.util
import json
import logging
import math
from pathlib import Path
import pickle
import subprocess
import sys

sys.dont_write_bytecode = True


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def cases(factory):
    def rows(specs):
        # specs: (gene, domains, residue, alt, score, af); positions stay sorted.
        return [factory.row(pos=10 + i, score=score, af=af, variant=f"A{residue}{alt}",
                            csq=factory.annotation(Gene=gene, DOMAINS=domains, Protein_position=f"{residue}/400"))
                for i, (gene, domains, residue, alt, score, af) in enumerate(specs)]

    mixed = []
    for residue in range(1, 9):  # eight residues: scored
        for alt, af in [("T", "0.2"), ("G", "0.06" if residue % 3 else "0.005")]:
            mixed.append(("ENSG_A", "Pfam:PF1&PANTHER:P1:SF2", residue, alt, f"0.{residue}5", af))
    for residue in range(20, 23):  # three residues: below the six-value minimum
        mixed.append(("ENSG_A", "Pfam:PF2", residue, "T", "0.9", "NaN"))
    for residue in range(30, 37):
        mixed.append(("ENSG_B", "Gene3D:G1&Pfam:PF1", residue, "C", "0.375", "0.0"))
    edges = []
    for residue in range(1, 21):  # 20 residues; exactly one or two fully common
        af = ("0.051" if residue == 1 else "0.05") if residue <= 2 else "0.0"
        edges.append(("ENSG_E", "Pfam:ONE", residue, "T", "0.5", af))
        edges.append(("ENSG_E", "Pfam:TWO", residue, "T", "0.5", "0.06" if residue <= 2 else "0.0"))
    for residue in range(1, 6):  # five residues only
        edges.append(("ENSG_E", "Pfam:FIVE", residue, "T", "0.5", "0.011"))
    for residue in range(1, 7):  # rare minimum exactly at the boundary
        edges.append(("ENSG_E", "Pfam:RARE", residue, "T", "0.5", "0.06"))
        edges.append(("ENSG_E", "Pfam:RARE", residue, "G", "0.5", "0.01"))
    nan = []
    for residue in range(1, 8):
        nan.append(("ENSG_N", "Pfam:FIRST_NAN", residue, "T", "0.5", "NaN"))
        nan.append(("ENSG_N", "Pfam:FIRST_NAN", residue, "G", "0.5", "0.2"))
        nan.append(("ENSG_N", "Pfam:LATER_NAN", residue, "T", "0.5", "0.2"))
        nan.append(("ENSG_N", "Pfam:LATER_NAN", residue, "G", "0.5", "NaN"))
    collision = []
    for residue in range(1, 8):
        # "ENSG_C:Pfam" + "X:PF9" and "ENSG_C" + "Pfam:X:PF9" join to one path.
        collision.append(("ENSG_C", "Pfam:X:PF9", residue, "T", "0.25", "0.2"))
    for residue in range(1, 10):
        collision.append(("ENSG_C:Pfam", "X:PF9", residue, "T", "0.75", "0.0"))
    return {
        "mixed": rows(mixed),
        "threshold_edges": rows(edges),
        "nan_order": rows(nan),
        "joined_path_collision": rows(collision),
        "empty": [],
    }


def tree_paths(value, path, found):
    for key, child in value.items():
        if isinstance(child, dict) and "distribution" not in key:
            joined = ":".join([*path, key])
            found.append(joined)
            tree_paths(child, [*path, key], found)


def expected_scores(module, model):
    data = []
    for gene, gene_data in model.items():
        for db, db_data in gene_data.items():
            data.extend(module.collect_domain_data(db_data, "unused", 6, path=[gene, db]))
    scores = {}
    for path, values, _ in data:
        scores[path] = values
    import numpy as np
    return [[path, len(values), hashlib.sha256(np.asarray(values, dtype=np.float64).astype(">f8").tobytes()).hexdigest()]
            for path, values in scores.items()]


def typed(value):
    if isinstance(value, bool):
        return ["bool", value]
    if isinstance(value, int):
        return ["int", value]
    if isinstance(value, float):
        return ["float", value.hex() if not math.isnan(value) else "nan"]
    return [type(value).__name__, repr(value)]


def run(command, folder, name):
    with (folder / (name + ".stdout")).open("wb") as out, (folder / (name + ".stderr")).open("wb") as err:
        return subprocess.run(list(map(str, command)), stdout=out, stderr=err, timeout=600).returncode


def check(args):
    logging.disable(logging.CRITICAL)
    args.work_root.mkdir(parents=True, exist_ok=False)
    factory = load(args.source_dir / "rust/tests/alpha_domain_compare.py", "projection_cases")
    sys.path.insert(0, str(args.source_dir / "scripts"))
    das = importlib.import_module("am_pick_intolerant_domains")
    results = []
    for name, case_rows in cases(factory).items():
        folder = args.work_root / name
        folder.mkdir()
        vcf = factory.make_case(folder, {"rows": case_rows})
        reference, native = folder / "reference.pkl", folder / "statistics.parquet"
        python = run([sys.executable, args.source_dir / "rust/tests/alpha_domain_compare.py", "oracle", "--procedure", "statistics",
                      "--source", args.source_dir / "scripts/stat_protein_domain_amscores.py", "--vcf", vcf,
                      "--output", reference, "--threads", "1"], folder, "python_statistics")
        built = run([args.binary, "alphamissense-domain", "build", "--vcf", vcf, "--output", native, "--threads", "1"],
                    folder, "native_statistics")
        record = dict(case=name, python_statistics_exit=python, native_statistics_exit=built)
        if python or built:
            record.update(passed=False, problem="statistics producer failed")
            results.append(record)
            print(json.dumps(record), flush=True)
            continue
        import pyarrow.parquet as pq
        builder = pq.ParquetFile(native).metadata.metadata[b"priva.builder_source_sha256"].decode()
        with reference.open("rb") as stream:
            model = pickle.load(stream)
        paths = []
        tree_paths(model, [], paths)
        requested = list(dict.fromkeys(paths)) + ["ABSENT:Pfam:NONE"]
        (folder / "domains.jsonl").write_text("".join(json.dumps(p) + "\n" for p in requested))
        scores_exit = run([args.binary, "das-domains", "scores", "--statistics", native, "--source-builder-sha256", builder],
                          folder, "native_scores")
        common_exit = run([args.binary, "das-domains", "common", "--statistics", native, "--source-builder-sha256", builder,
                           "--domains", folder / "domains.jsonl"], folder, "native_common")
        record.update(scores_exit=scores_exit, common_exit=common_exit)
        problems = []
        if scores_exit or common_exit:
            problems.append("native projection failed")
        else:
            native_scores = [json.loads(line) for line in (folder / "native_scores.stdout").read_text().splitlines()]
            reference_scores = expected_scores(das, model)
            record.update(score_domains=len(reference_scores), native_score_domains=len(native_scores))
            if native_scores != reference_scores:
                problems.append("ordered score arrays differ")
            summary = das.filter_common_missense(set(requested), model, 0.05)[2]
            native_common = [json.loads(line) for line in (folder / "native_common.stdout").read_text().splitlines()]
            fields = ["n_fully_common_residues", "n_total_residues", "frac_common", "removed"]
            expected_rows = [[domain, [typed(summary[domain][f]) for f in fields]] for domain in requested]
            actual_rows = [[row["domain"], [typed(row[f]) for f in fields]] for row in native_common]
            record.update(common_domains=len(expected_rows), removed=sum(summary[d]["removed"] for d in requested))
            if actual_rows != expected_rows:
                problems.append("common-missense summaries differ")
                (folder / "common_difference.json").write_text(json.dumps(dict(expected=expected_rows, actual=actual_rows), indent=1))
            if native_common and list(native_common[0]) != ["domain", *fields]:
                problems.append("common-missense field order differs")
        record.update(passed=not problems, problems=problems)
        results.append(record)
        print(json.dumps(record), flush=True)
    summary = dict(passed=all(r["passed"] for r in results), cases=len(results), results=results,
                   flat_fallback_covered_here=False)
    (args.work_root / "complete.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "results"}))
    if not summary["passed"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    test = commands.add_parser("check")
    for name in ("binary", "source-dir", "work-root"):
        test.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    globals()[args.command](args)


if __name__ == "__main__":
    main()
