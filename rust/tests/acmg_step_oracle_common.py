"""Shared Python side of the group 3 comparisons (the ACMG step, Python to Rust).

A thread captures cases by running the real Python function in the priva_acmg
environment and writing its input and output frames; `priva acmg-step` then runs
the Rust function on the same input, and `compare_frames` checks every cell
(floats by their bits, no tolerance).

    import acmg_step_oracle_common as oracle
    oracle.write_frame(df, case / "input.frame")
    oracle.capture(case, lambda: PM2_criteria(df, ...), to_frame=...)
    oracle.run_rust(priva, "acmg_pm2_bs1_ba1_frequency", "PM2_criteria", case)
    oracle.compare_frames(case / "expected.frame", case / "observed.frame")

The frame file format is described in rust/src/acmg_step/frame_file.rs.
"""
import hashlib
import json
import math
import platform
import struct
import subprocess
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

MAGIC = "priva-frame"
DTYPES = {"int64", "float64", "bool", "object", "Int64", "boolean"}
PINNED = {"python": "3.12.3", "pandas": "2.2.3", "numpy": "1.26.4"}


def _hex(text):
    return text.encode("utf-8").hex()


def encode_cell(value, dtype):
    """One cell of the frame file. Raises TypeError for values the format does not carry."""
    if value is pd.NA:
        if dtype not in ("Int64", "boolean"):
            raise TypeError(f"pd.NA in a {dtype} column")
        return "A"
    if value is None:
        return "N"
    if isinstance(value, (bool, np.bool_)):
        return "b1" if value else "b0"
    if isinstance(value, (int, np.integer)):
        if not -(2 ** 63) <= int(value) < 2 ** 63:
            raise TypeError(f"integer {value} beyond 64 bits")
        return f"i{int(value)}"
    if isinstance(value, (float, np.floating)):
        if math.isnan(value):
            return "n"
        return "f" + struct.pack(">d", float(value)).hex()
    if isinstance(value, str):
        return "s" + _hex(value)
    raise TypeError(f"value of type {type(value).__name__} is not carried by frame files")


def frame_text(df):
    if not isinstance(df.index, pd.RangeIndex) or df.index.start != 0 or df.index.step != 1:
        raise ValueError("frames carry the default index 0..n; reset_index(drop=True) first")
    dtypes = [str(dtype) for dtype in df.dtypes]
    unknown = sorted(set(dtypes) - DTYPES)
    if unknown:
        raise TypeError(f"dtypes {unknown} are not carried by frame files")
    lines = [f"{MAGIC}\t1\t{len(df)}", "\t".join(_hex(str(name)) for name in df.columns), "\t".join(dtypes)]
    columns = [df.iloc[:, position].tolist() if dtypes[position] in ("Int64", "boolean")
               else list(df.iloc[:, position].to_numpy(dtype=object))
               for position in range(df.shape[1])]
    for row in range(len(df)):
        lines.append("\t".join(encode_cell(column[row], dtype) for column, dtype in zip(columns, dtypes)))
    return "\n".join(lines) + "\n"


def write_frame(df, path):
    Path(path).write_text(frame_text(df), encoding="utf-8")


def _decode_cell(cell):
    kind, payload = cell[:1], cell[1:]
    if cell == "n":
        return np.nan
    if cell == "N":
        return None
    if cell == "A":
        return pd.NA
    if kind == "b":
        return payload == "1"
    if kind == "i":
        return int(payload)
    if kind == "f":
        return struct.unpack(">d", bytes.fromhex(payload))[0]
    if kind == "s":
        return bytes.fromhex(payload).decode("utf-8")
    raise ValueError(f"bad frame cell {cell!r}")


def read_frame(path):
    lines = Path(path).read_text(encoding="utf-8").split("\n")
    magic, version, rows = lines[0].split("\t")
    if magic != MAGIC or version != "1":
        raise ValueError(f"{path}: not a version 1 frame file")
    names = [bytes.fromhex(name).decode("utf-8") for name in lines[1].split("\t")] if lines[1] else []
    dtypes = lines[2].split("\t") if lines[2] else []
    cells = [line.split("\t") if names else [] for line in lines[3:3 + int(rows)]]
    data = {}
    for position, (name, dtype) in enumerate(zip(names, dtypes)):
        values = [_decode_cell(row[position]) for row in cells]
        data[position] = pd.Series(values, dtype=dtype if dtype != "object" else object)
    df = pd.DataFrame(data, index=pd.RangeIndex(int(rows)))
    df.columns = names
    return df


def compare_frames(expected_path, observed_path, show=12):
    """Cell-by-cell comparison of the frame texts. Returns a list of differences (empty when equal)."""
    expected = Path(expected_path).read_text(encoding="utf-8").split("\n")
    observed = Path(observed_path).read_text(encoding="utf-8").split("\n")
    differences = []
    if expected[:3] != observed[:3]:
        differences.append(f"header lines differ: {expected[:3]!r} vs {observed[:3]!r}")
        return differences
    names = [bytes.fromhex(name).decode() for name in expected[1].split("\t")] if expected[1] else []
    for number, (left, right) in enumerate(zip(expected[3:], observed[3:])):
        if left == right:
            continue
        for name, a, b in zip(names, left.split("\t"), right.split("\t")):
            if a != b and len(differences) < show:
                differences.append(f"row {number} column {name}: expected {_decode_cell(a)!r} observed {_decode_cell(b)!r}")
    if len(expected) != len(observed):
        differences.append(f"line counts differ: {len(expected)} vs {len(observed)}")
    return differences


def capture(case_dir, function, to_frame):
    """Runs `function()`; writes expected.frame (through `to_frame(result)`) or expected_error.json."""
    case_dir = Path(case_dir)
    try:
        result = function()
    except Exception as error:
        (case_dir / "expected_error.json").write_text(json.dumps({
            "type": type(error).__name__, "message": str(error),
            "where": traceback.extract_tb(error.__traceback__)[-1].name,
        }, indent=1, sort_keys=True))
        return None
    write_frame(to_frame(result), case_dir / "expected.frame")
    return result


def run_rust(priva, module, function, case_dir, options=None):
    """`priva acmg-step MODULE FUNCTION` on case_dir/input.frame; returns (exit code, stderr)."""
    case_dir = Path(case_dir)
    command = [str(priva), "acmg-step", module, function, "--input", str(case_dir / "input.frame"),
               "--output", str(case_dir / "observed.frame")]
    if options is not None:
        command += ["--options", json.dumps(options, sort_keys=True)]
    done = subprocess.run(command, capture_output=True, text=True)
    return done.returncode, done.stderr


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def manifest(scripts, resources=(), seed=None):
    """Versions, hashes of the Python scripts and resources a capture used. Refuses other versions."""
    versions = {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__}
    try:
        import scipy
        versions["scipy"] = scipy.__version__
    except ImportError:
        versions["scipy"] = None
    for name, wanted in PINNED.items():
        if versions[name] != wanted:
            sys.exit(f"{name} {versions[name]} is not the pinned {wanted}; refusing to capture")
    return {
        "versions": versions,
        "scripts": {str(path): sha256_file(path) for path in scripts},
        "resources": {str(path): sha256_file(path) for path in resources},
        "seed": seed,
    }
