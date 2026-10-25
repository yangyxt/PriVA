#!/usr/bin/env python3
"""Read original gzip bundles and regenerated plain InterPro pickle bytes.

Consumers share this byte-magic decision instead of inferring compression from
the historical .pkl.gz filename. No mapping values, order or multiplicity are
changed; ordinary pickle decoding errors propagate without another file choice
or additional stream integrity checks.
"""
from __future__ import annotations

import argparse
from collections.abc import Sized
import gzip
import json
from pathlib import Path
import pickle
import sys


def load_interpro_mapping(path: str | Path):
    """Deserialize the selected file without rewriting it or guessing a fallback."""
    with open(path, "rb") as original:
        magic = original.read(2)
        original.seek(0)
        if magic == b"\x1f\x8b":
            with gzip.GzipFile(fileobj=original, mode="rb") as compressed:
                # Preserve the ordinary scientific read for bundled gzip data.
                return pickle.load(compressed)
        return pickle.load(original)


def inspect_mapping(path: str | Path) -> dict:
    """Return bounded metadata for the CLI; never print mapping contents."""
    value = load_interpro_mapping(path)
    return {"input": str(path), "value_type": type(value).__name__,
            "entries": len(value) if isinstance(value, Sized) else None}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    inspect = commands.add_parser("inspect", help="Read the selected mapping and print metadata only")
    inspect.add_argument("--input", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        print(json.dumps(inspect_mapping(args.input), sort_keys=True))
    except (OSError, EOFError, ValueError, pickle.PickleError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
