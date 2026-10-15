#!/usr/bin/env python3
"""Differential-test source-mapped scoring against the actual frozen Python AST.

This qualifies already-assigned integer criteria, default posterior parameters,
and pre-resolved ClinGen strings. It does not qualify criterion assignment,
pickle conversion, ranking penalties or any real-family scientific output.
"""
import argparse
import ast
import csv
import gzip
import hashlib
import io
import itertools
import logging
import os
from pathlib import Path
import pickle
import random
import subprocess
import tempfile
from typing import Dict, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SOURCE_SHA256 = "9111a76e3fc94a064e14d97b535bb1dd2132f7c61836f8e8892477d4c5e69646"
CRITERIA = ["PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
            "PP1", "PP2", "PP3", "PP4", "PP5", "BA1", "BS1", "BS2", "BS3", "BS4", "BP1",
            "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"]


def reference_functions():
    source = (ROOT / "scripts/acmg_scoring.py").read_bytes()
    actual_sha = hashlib.sha256(source).hexdigest()
    if actual_sha != SOURCE_SHA256:
        raise AssertionError(f"reference source changed: expected {SOURCE_SHA256}, found {actual_sha}")
    names = {"fix_clingen_term", "translate_strength_level", "create_criteria_summary",
             "summarize_acmg_criteria", "calculate_posterior_probability"}
    nodes = [node for node in ast.parse(source).body
             if isinstance(node, ast.FunctionDef) and node.name in names]
    assert {node.name for node in nodes} == names
    namespace = {"pd": pd, "np": np, "Dict": Dict, "Tuple": Tuple, "os": os,
                 "gzip": gzip, "pickle": pickle, "logger": logging.getLogger("frozen-scoring")}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), "frozen-acmg-scoring", "exec"), namespace)
    return namespace


def rows_for_fixture():
    rng = random.Random(271828)
    rows = []
    # Exhaust every strength across two strong, moderate, supporting and benign
    # slots: includes 2-PS adjustment, BA override and the original ignored
    # pathogenic strength 5 / benign strength 4 behavior.
    for strengths in itertools.product(range(6), repeat=5):
        row = dict.fromkeys(CRITERIA, 0)
        row.update(zip(["PS1", "PS2", "PM1", "PP3", "BS1"], strengths))
        rows.append(row)
    for _ in range(2048):
        rows.append({name: rng.randrange(6) if rng.random() < 0.45 else 0 for name in CRITERIA})
    # Source's special PVS + PM branch counts only one PVS and requires NO PP.
    for pvs_count, pm_count, pp_count in itertools.product(range(4), range(4), range(3)):
        row = dict.fromkeys(CRITERIA, 0)
        for name in CRITERIA[:pvs_count]: row[name] = 4
        for name in ["PM1", "PM2", "PM3"][:pm_count]: row[name] = 2
        for name in ["PP1", "PP2"][:pp_count]: row[name] = 1
        rows.append(row)
    evidence_values = [None, float("nan"), 17, "", " \t\n", "\x1c\x1f",
                       "PVS1_Very, PM2_Supporting", "BA1_Stand", "PS1,PS1_Supporting",
                       "PS1_Supporting,PS1_Strong", "PP3_Strong, BP4", ", ,",
                       "  PM2_Moderate, \x1cPP3\x1f ", "PVS1_Very_extra, BA1_Standalone",
                       "PM2_Supporting_extra", "UNKNOWN", "PS1;PM2_Supporting"]
    evidence_map = {}
    for index, row in enumerate(rows):
        row["variant_id"] = f"v{index}"
        row["Feature"] = f"ENST{index % 5}"
        row["unchanged_note"] = 'tab\tquote";unicode β\nsecond line'
        evidence = evidence_values[index % len(evidence_values)]
        evidence_map[row["variant_id"]] = evidence
        # This boundary consumes pre-resolved strings. Original non-string map
        # entries have no replacement effect and export as an empty field.
        row["resolved_clingen"] = evidence if isinstance(evidence, str) else ""
        row["ACMG_class"] = "replace old class"
        row["ACMG_quant_score"] = "replace old score"
        row["ACMG_criteria"] = "replace old summary"
    # Duplicate IDs and row order must survive; resolved evidence is identical.
    rows.extend([dict(rows[3]), dict(rows[20]), dict(rows[3])])
    return pd.DataFrame(rows), evidence_map


def oracle(frame, evidence_map, enabled, scope, folder):
    path = folder / "clingen.pkl"
    with path.open("wb") as stream:
        pickle.dump(evidence_map, stream)
    result, matrix = scope["summarize_acmg_criteria"](
        frame.copy(), {name: frame[name].to_numpy() for name in CRITERIA},
        str(path), apply_clingen_override=enabled)
    result[CRITERIA] = matrix[CRITERIA]
    if len(matrix):
        probabilities, classes = zip(*matrix.apply(scope["calculate_posterior_probability"], axis=1))
        result["ACMG_quant_score"] = probabilities
        result["ACMG_class"] = classes
    else:
        result["ACMG_quant_score"] = pd.Series(dtype=float)
        result["ACMG_class"] = pd.Series(dtype=str)
    return result


def parsed(text):
    return list(csv.reader(io.StringIO(text), delimiter="\t"))


def invoke(binary, folder, frame, enabled=True):
    input_path, output_path = folder / "input.tsv", folder / "output.tsv"
    if output_path.exists(): output_path.unlink()
    frame.to_csv(input_path, sep="\t", index=False)
    original = input_path.read_bytes()
    flag = ["--clingen-evidence-column", "resolved_clingen"] if enabled else ["--no-clingen-override"]
    process = subprocess.run([str(binary), "score-criteria", str(input_path), str(output_path), *flag],
                             text=True, capture_output=True, env={**os.environ, "POLARS_MAX_THREADS": "1"})
    assert input_path.read_bytes() == original, "input was mutated"
    return process, output_path


def run(binary):
    scope = reference_functions()
    frame, evidence_map = rows_for_fixture()
    comparisons = 0
    failures = 0
    with tempfile.TemporaryDirectory(prefix="priva-scoring-parity-", dir=os.environ.get("TMPDIR")) as temporary:
        folder = Path(temporary)
        for enabled in (False, True):
            expected = oracle(frame, evidence_map, enabled, scope, folder)
            process, output_path = invoke(binary, folder, frame, enabled)
            if process.returncode: raise AssertionError(process.stderr)
            actual_rows, expected_rows = parsed(output_path.read_text()), parsed(expected.to_csv(sep="\t", index=False))
            if actual_rows != expected_rows:
                for index, (actual, reference) in enumerate(zip(actual_rows, expected_rows)):
                    if actual != reference:
                        differences = [(expected_rows[0][i], a, r) for i, (a, r) in enumerate(zip(actual, reference)) if a != r]
                        raise AssertionError(f"row {index}, override={enabled}: {differences}")
                raise AssertionError("output row count differs")
            comparisons += len(frame)

        # All-BA tables infer an integer score column in pandas; mixed tables
        # infer float. Both serialize exactly, including 0 versus 0.0.
        for count in (1, 3):
            benign = frame.iloc[:count].copy()
            benign["BA1"] = 5
            expected = oracle(benign, {}, False, scope, folder)
            process, output_path = invoke(binary, folder, benign, False)
            assert process.returncode == 0, process.stderr
            assert parsed(output_path.read_text()) == parsed(expected.to_csv(sep="\t", index=False))
            comparisons += count

        # Preserve the source's empty DataFrame.apply failure rather than
        # silently invent a successful empty scientific result.
        empty = frame.iloc[:0]
        try:
            oracle(empty, evidence_map, False, scope, folder)
        except ValueError:
            pass
        else:
            raise AssertionError("empty source behavior changed; requalify explicitly")
        process, output_path = invoke(binary, folder, empty, False)
        assert process.returncode != 0 and not output_path.exists()
        assert "ValueError:" in process.stderr
        failures += 1

        # Original prefix-match quirks/malformed strengths must fail with the
        # same Python exception category, never be normalized or ignored.
        for evidence in ("PS1_Unknown", "XPS1", "PS10", "PS1_", "PS1_ Supporting", "PVS1_Unknown"):
            bad = frame.iloc[:1].copy()
            bad.loc[0, "resolved_clingen"] = evidence
            try:
                oracle(bad, {"v0": evidence}, True, scope, folder)
            except Exception as error:
                category = type(error).__name__
            else:
                raise AssertionError(f"fixture unexpectedly valid: {evidence}")
            process, output_path = invoke(binary, folder, bad)
            assert process.returncode != 0 and not output_path.exists()
            assert f"{category}:" in process.stderr, process.stderr
            failures += 1

        # The reference rebuilds each column across all rows before the next
        # criterion. Its PVS1 KeyError in row 2 precedes row 1's PS1 IndexError.
        bad = frame.iloc[:2].copy()
        bad.loc[0, "resolved_clingen"] = "XPS1"
        bad.loc[1, "resolved_clingen"] = "PVS1_Unknown"
        try:
            oracle(bad, {"v0": "XPS1", "v1": "PVS1_Unknown"}, True, scope, folder)
        except KeyError:
            pass
        else:
            raise AssertionError("source override-column order changed")
        process, output_path = invoke(binary, folder, bad)
        assert process.returncode != 0 and not output_path.exists()
        assert "KeyError:" in process.stderr, process.stderr
        failures += 1

        # Reject non-typed/missing strengths explicitly; do not silently apply
        # pandas inference or treat unknown evidence as zero.
        for value in ("", "NaN", "1.5", "1.0", "6", "-1"):
            bad = frame.iloc[:1].astype(object).copy()
            bad.loc[0, "PS1"] = value
            process, output_path = invoke(binary, folder, bad)
            assert process.returncode != 0 and not output_path.exists(), value
            failures += 1

        bad = frame.iloc[:1].drop(columns=["PM2"])
        process, output_path = invoke(binary, folder, bad)
        assert process.returncode != 0 and not output_path.exists()
        failures += 1

        # create_new refuses overwrite, including accidentally choosing input.
        good = frame.iloc[:1]
        process, output_path = invoke(binary, folder, good, False)
        assert process.returncode == 0, process.stderr
        original = output_path.read_bytes()
        process = subprocess.run([str(binary), "score-criteria", str(folder / "input.tsv"), str(output_path),
                                  "--no-clingen-override"], capture_output=True)
        assert process.returncode != 0 and output_path.read_bytes() == original
        failures += 1

    print(f"PASS: {comparisons} exact parsed-cell differential rows (ClinGen OFF/ON), {failures} failure/overwrite cases including empty source failure")
    print(f"source_commit=09223ef5a8b0ceb4c6d61055c6eefe62b51ba903 source_sha256={SOURCE_SHA256}")
    print(f"pandas={pd.__version__} numpy={np.__version__}; default model only; NOT criterion-assignment or family parity")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    run(parser.parse_args().binary.resolve())
