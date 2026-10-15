#!/usr/bin/env python3
"""Compare Rust with the ACTUAL terminal block of legacy sort_and_rank_variants.

Only the terminal block (group max, sort, factorize) is under test. The earlier
SciPy/control, expression and inheritance penalties and ACMG judgement are not
implemented by this Rust module and are not certified by this test.
"""
import argparse
import ast
import csv
import io
import os
from pathlib import Path
import random
import subprocess
import tempfile

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def legacy_terminal_block():
    tree = ast.parse((ROOT / "scripts/acmg_scoring.py").read_text())
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                    and node.name == "sort_and_rank_variants")
    start = next(index for index, node in enumerate(function.body)
                 if isinstance(node, ast.Assign)
                 and any(isinstance(target, ast.Name) and target.id == "variant_groups" for target in node.targets))
    function.body = function.body[start:]
    function.args = ast.arguments(posonlyargs=[], args=[ast.arg(arg="df")], kwonlyargs=[],
                                  kw_defaults=[], defaults=[])
    function.returns = None
    scope = {"pd": pd}
    exec(compile(ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[])),
                 "legacy-terminal-order", "exec"), scope)
    return scope["sort_and_rank_variants"]


def parsed(text):
    return list(csv.reader(io.StringIO(text), delimiter="\t"))


def run(binary):
    oracle = legacy_terminal_block()
    rng = random.Random(161803)
    cases = 0
    with tempfile.TemporaryDirectory(prefix="priva-rust-parity-", dir=os.environ.get("TMPDIR")) as folder:
        folder = Path(folder)
        for chrom_type in ("string", "integer"):
            for count in (0, 1, 2, 8, 40, 200):
                for replicate in range(5):
                    rows = []
                    for index in range(count):
                        variant = rng.randrange(max(1, count // 3))
                        chromosome = ("chr2", "chr10", "chrX")[variant % 3] if chrom_type == "string" else (2, 10, 22)[variant % 3]
                        rows.append({"chrom": chromosome, "pos": 1 + variant // 3,
                                     "ref": "A", "alt": "G", "variant_id": f"v{variant}",
                                     "Feature": f"ENST{index}",
                                     "sort_index": rng.choice([0.0, 0.95, 0.8, 0.1, 0.9500000000000001]),
                                     "ACMG_quant_score": rng.choice([0.1, 0.8, 0.95]),
                                     "ACMG_criteria": "PM2_Supporting;PP3" if index % 2 else "",
                                     "quoted_note": 'literal\ttab and "quote"',
                                     "variant_rank": 999})
                    frame = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "variant_id", "Feature",
                                                       "sort_index", "ACMG_quant_score", "ACMG_criteria", "quoted_note", "variant_rank"])
                    # Empty object columns still permit the legacy terminal block.
                    expected = oracle(frame.copy()).to_csv(sep="\t", index=False)
                    input_path, output_path = folder / "input.tsv", folder / "output.tsv"
                    frame.to_csv(input_path, sep="\t", index=False)
                    result = subprocess.run([str(binary), "final-order", str(input_path), str(output_path),
                                             "--chrom-type", chrom_type], text=True, capture_output=True)
                    if result.returncode:
                        raise AssertionError(result.stderr)
                    actual = output_path.read_text()
                    if parsed(actual) != parsed(expected):
                        raise AssertionError(f"Rust mismatch: chrom={chrom_type}, rows={count}, replicate={replicate}")
                    output_path.unlink()
                    cases += 1
        # Invalid/nonfinite keys must fail; they are outside this first module's contract.
        for field, value in (("sort_index", "NaN"), ("ACMG_quant_score", "inf"),
                             ("variant_id", ""), ("pos", "1.5"), ("chrom", "")):
            row = {"chrom": "chr1", "pos": "1", "ref": "A", "alt": "G", "variant_id": "v1",
                   "sort_index": "0.9", "ACMG_quant_score": "0.9"}
            row[field] = value
            pd.DataFrame([row]).to_csv(input_path, sep="\t", index=False)
            result = subprocess.run([str(binary), "final-order", str(input_path), str(output_path),
                                     "--chrom-type", "string"], text=True, capture_output=True)
            if result.returncode == 0 or output_path.exists():
                raise AssertionError(f"invalid input accepted: {field}={value!r}")
    print(f"PASS: {cases} seeded differential cases against actual legacy terminal block; 5 rejection cases")
    print(f"pandas={pd.__version__}; this is NOT whole-function or single-family scientific parity")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    run(parser.parse_args().binary.resolve())
