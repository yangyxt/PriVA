#!/usr/bin/env python3
"""Execute the unchanged pinned Python function with actual pandas, not a reimplementation.

Use: python condition_masks_oracle.py capture --source-file scripts/acmg_variant_mechanism.py
     --output DATA_ROOT/oracle
Only primitive typed cells and a single string/None index name cross the CSV boundary.
The original function AST is compiled verbatim; unrelated hub imports/model inputs are not run.
"""
import argparse
import ast
import csv
import hashlib
import json
import platform
import random
import re
import struct
from pathlib import Path

import numpy as np
import pandas as pd

SOURCE_SHA = "c7e8da80ca0ec47fb38a5fd1d965641c1791189cb0d25a651253f44508d0bf26"
REFERENCE_COMMIT = "09223ef5a8b0ceb4c6d61055c6eefe62b51ba903"
INHERITANCE, PENETRANCE = "variant_inheritance", "variant_penetrance"


def encode(value):
    if value is None or value is pd.NA:
        return ["null", ""]
    if isinstance(value, (bool, np.bool_)):
        return ["bool", str(int(value))]
    if isinstance(value, (int, np.integer)):
        return ["uint" if value > 2**63 - 1 else "int", str(value)]
    if isinstance(value, (float, np.floating)):
        return ["float", struct.pack(">d", value).hex()]
    if isinstance(value, str):
        return ["str", value.encode().hex()]
    raise TypeError(f"outside typed primitive transport: {type(value)}")


def frame(columns, *, index=None, index_name="gapped_duplicate"):
    """Keep per-position dtypes and duplicate column labels, without dict collapse."""
    size = len(columns[0][1]) if columns else (len(index) if index is not None else 1)
    if columns:
        df = pd.concat([pd.Series(values, dtype=dtype) for _, values, dtype in columns], axis=1)
        df.columns = [name for name, _, _ in columns]
    else:
        df = pd.DataFrame(index=range(size))
    df.index = (pd.Index([7 if i % 3 else -2 for i in range(size)], name=index_name)
                if index is None else index)
    df.columns.name = "retained_input_column_axis"
    df.attrs["fixture"] = "condition masks; unchanged source"
    return df


def pair(inheritance, penetrance=None, *, dtype="object", index=None, index_name="gapped_duplicate"):
    return frame([(INHERITANCE, inheritance, dtype),
                  (PENETRANCE, [""] * len(inheritance) if penetrance is None else penetrance, "object")],
                 index=index, index_name=index_name)


def cases():
    inheritance_tokens = ["recessive", "dominant", "x_linked_recessive", "x_linked_dominant",
                          "x_linked_unspecified", "y_linked", "mitochondrial", "polygenic",
                          "digenic", "oligogenic", "non_mendelian"]
    templates = ["{}", ";{};", ";;{};;", "x;{};y", "{}\n", "{}\n\n", "{}\r\n",
                 "\n{}", ";{}\n", "{}\n;", "{};\n", "x;{}\n", "{}\u2028",
                 " {}", "{} ", "\t{}", "{}\t", "x{}", "{}x", "{}_LOF", "{}+dominant"]
    values = [template.format(token) for token in inheritance_tokens for template in templates]
    values += ["", ";", "dominantly_inherited", "non-monogenic", "non_mendelian_like",
               "recessive\0", "recessive;\0;dominant", "dominant;recessive;polygenic",
               ";".join(inheritance_tokens), "\n;recessive;dominant\n", "recessive\n;dominant"]
    yield "tokens_and_regex_boundaries", pair(values, ["complete;incomplete"] * len(values))
    values = [template.format(token) for token in ["complete", "incomplete", "high", "unknown"]
              for template in templates]
    yield "penetrance_boundaries", pair(["non_mendelian"] * len(values), values)
    values = ["DOMINANT", "X_LINKED_RECESSIVE", "X_LIN\u212aED_DOMINANT", "M\u0130TOCHONDRIAL",
              "POLYGENIC;DIGENIC;OLIGOGENIC", "\u00a0recessive", "recessive\u001c", "\u200brecessive",
              "\uff38_LINKED_DOMINANT", "non_MENDELIAN", "\u03a3;RECESSIVE", "mitochondrial"]
    yield "lowercase_without_strip", pair(values, ["INCOMPLETE;COMPLETE"] * len(values))
    scalars = [None, pd.NA, float("nan"), float("inf"), -float("inf"), True, False, -1, 0,
               2**63 - 1, 2**63, 2**64 - 1, 0.0, -0.0, 2.5, "None", "NaN", "<NA>", ""]
    yield "object_missing_and_scalars", pair(scalars, list(reversed(scalars)))
    rng = random.Random(267340)
    choices = inheritance_tokens + [t.upper() for t in inheritance_tokens] + ["", "x", "recessive\n", "complete"]
    values = [";".join(rng.choices(choices, k=rng.randrange(1, 9))) for _ in range(600)]
    yield "ordered_mixed_histories", pair(values, [rng.choice(["complete", "incomplete", "complete;incomplete", "", None]) for _ in values])
    dtype_values = {
        "object": [[None, pd.NA, float("nan")], ["dominant", "RECESSIVE"], [1, True], []],
        "string": [[None, pd.NA], ["dominant", "INCOMPLETE"], ["", "none"], []],
        "Int64": [[1, 2], [1, pd.NA], [pd.NA, pd.NA], []],
        "Float64": [[1.0, float("inf")], [1.0, pd.NA], [pd.NA, pd.NA], []],
        "boolean": [[True, False], [True, pd.NA], [pd.NA, pd.NA], []],
        "float64": [[1.0, float("nan")], [float("inf"), -float("inf")], []],
        "int64": [[-1, 2], []], "uint64": [[2**63, 2**64 - 1], []], "bool": [[True, False], []],
    }
    for dtype, groups in dtype_values.items():
        for number, values in enumerate(groups):
            for phase, name in enumerate([INHERITANCE, PENETRANCE]):
                other = PENETRANCE if phase == 0 else INHERITANCE
                yield f"dtype_{dtype}_{number}_{phase}", frame([(name, values, dtype), (other, [""] * len(values), "object")])
    for name, columns in [("both", []), ("inheritance", [(PENETRANCE, [pd.NA], "boolean")]),
                          ("penetrance", [(INHERITANCE, [pd.NA], "Int64")])]:
        yield "missing_" + name, frame(columns)
    yield "inheritance_before_penetrance", frame([(PENETRANCE, [pd.NA], "boolean"),
                                                  (INHERITANCE, [pd.NA], "Int64")])
    yield "blank_inheritance_does_not_skip_penetrance", frame([(INHERITANCE, [""], "object"),
                                                               (PENETRANCE, [pd.NA], "Float64")])
    duplicates = [
        ("inheritance", [(INHERITANCE, [""], "object"), (INHERITANCE, [""], "object"), (PENETRANCE, [pd.NA], "boolean")]),
        ("penetrance", [(INHERITANCE, [""], "object"), (PENETRANCE, [""], "object"), (PENETRANCE, [""], "object")]),
        ("inheritance_nullable", [(INHERITANCE, [""], "object"), (INHERITANCE, [pd.NA], "Int64"), (PENETRANCE, [pd.NA], "boolean")]),
        ("penetrance_nullable", [(INHERITANCE, [""], "object"), (PENETRANCE, [""], "object"), (PENETRANCE, [pd.NA], "boolean")]),
        ("two_nullable", [(PENETRANCE, [pd.NA], "boolean"), (INHERITANCE, [pd.NA], "Float64"), (INHERITANCE, [pd.NA], "Int64")]),
        ("required_before_duplicate", [(INHERITANCE, [pd.NA], "Int64"), (INHERITANCE, [pd.NA], "boolean")]),
        ("unrelated", [("ignored", [pd.NA], "boolean"), ("ignored", [pd.NA], "Int64"), (PENETRANCE, ["complete"], "object"), (INHERITANCE, ["dominant"], "object")]),
        ("empty_inheritance", [(INHERITANCE, [], "Int64"), (INHERITANCE, [], "boolean"), (PENETRANCE, [], "object")]),
    ]
    for name, columns in duplicates:
        yield "duplicate_" + name, frame(columns)
    for name, index in [
        ("object", pd.Index(["same", "", "same", None], dtype=object, name="")),
        ("float", pd.Index([1.0, float("nan"), -0.0, 1.0], name="floating_labels")),
        ("uint", pd.Index([2**64 - 1, 2**63, 0, 2**63], dtype="uint64", name="unsigned_labels")),
        ("boolean", pd.Index([True, False, True, False], name=None)),
        ("string", pd.Index(["a", None, "a", ""], dtype="string", name="nullable_labels")),
    ]:
        yield "index_" + name, pair(["recessive", "dominant", "", None], ["incomplete", "complete", "", None], index=index)


def capture(args):
    source = args.source_file.read_bytes()
    assert hashlib.sha256(source).hexdigest() == SOURCE_SHA, "original source content changed"
    tree = ast.parse(source, filename=str(args.source_file))
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "_variant_condition_masks")
    namespace = {"pd": pd, "re": re}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(args.source_file), "exec"), namespace)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "original_acmg_variant_mechanism.py").write_bytes(source)
    result = {
        "harness_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reference_source_commit": REFERENCE_COMMIT, "source_path": str(args.source_file),
        "source_sha256": SOURCE_SHA, "function_ast_sha256": hashlib.sha256(ast.dump(function).encode()).hexdigest(),
        "runtime": {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__},
        "execution": "unchanged function AST with actual pandas and Python re", "cases": [],
        "scientific_acceptance": False, "production_backend_selected": False,
        "scope": "typed primitive cells, declared pandas dtype, opaque primitive index labels and one string/None index name",
    }
    for name, df in cases():
        before = df.copy(deep=True)
        entry = {"name": name, "rows": len(df), "columns": list(df.columns),
                 "dtypes": [str(v) for v in df.dtypes], "index_names": list(df.index.names), "index_dtype": str(df.index.dtype)}
        with (args.output / (name + ".input.csv")).open("x", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["index_kind", "index_value"] + [v for c in df.columns for v in (c + "_kind", c)])
            for i in range(len(df)):
                row = [df.iloc[i, j] for j in range(len(df.columns))]
                writer.writerow(encode(df.index[i]) + [x for v in row for x in encode(v)])
        try:
            masks = namespace["_variant_condition_masks"](df)
        except Exception as error:
            entry["exception"] = {"type": type(error).__name__, "message": str(error)}
        else:
            entry["keys"] = list(masks)
            entry["series_metadata"] = [{"name": v.name, "dtype": str(v.dtype),
                "index_names": list(v.index.names), "index_dtype": str(v.index.dtype)} for v in masks.values()]
            assert len(masks) == 13 and all(v.index.equals(df.index) and str(v.dtype) == "bool" for v in masks.values())
            with (args.output / (name + ".expected.csv")).open("x", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(["index_kind", "index_value", *masks])
                for i in range(len(df)):
                    writer.writerow(encode(df.index[i]) + [int(v.iloc[i]) for v in masks.values()])
        pd.testing.assert_frame_equal(before, df, check_exact=True)
        assert before.attrs == df.attrs
        entry["input_unchanged"] = True
        result["cases"].append(entry)
    assert hashlib.sha256(args.source_file.read_bytes()).hexdigest() == SOURCE_SHA
    (args.output / "oracle.json").write_text(json.dumps(result, indent=2) + "\n")
    manifest = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(args.output.iterdir()) if p.is_file()}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"cases": len(result["cases"]), "rows": sum(c["rows"] for c in result["cases"]),
                      "expected_failures": sum("exception" in c for c in result["cases"]),
                      "manifest_sha256": hashlib.sha256((args.output / "manifest.json").read_bytes()).hexdigest()}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    command = sub.add_parser("capture")
    command.add_argument("--source-file", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    capture(parser.parse_args())


if __name__ == "__main__":
    main()
