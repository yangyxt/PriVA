#!/usr/bin/env python3
"""Check newly acquired source formats without changing scientific readers.

Usage: python mechanism_source_input.py validate --kind clingen_erepo --input FILE
Checks apply only to downloads/explicit URL overrides, never genuine bundles.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import tempfile
from urllib.parse import urlsplit
import xml.etree.ElementTree as ET


def source_url(value: str) -> str:
    parsed = urlsplit(value)
    if parsed.scheme not in ("http", "https") or not parsed.netloc or any(c.isspace() for c in value):
        raise ValueError("Source URL must be an absolute HTTP(S) URL")
    return value


def validate_source_payload(path: Path, kind: str) -> None:
    """Require the existing reader's source columns, or complete Product 6 XML."""
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"Downloaded {kind} resource is missing or empty: {path}")
    if kind == "orphadata":
        root = None
        has_disorders = False
        for event, element in ET.iterparse(path, events=("start", "end")):
            if "}" in element.tag:
                raise ValueError("The original Orphadata reader requires unnamespaced Product 6 XML")
            tag = element.tag.rsplit("}", 1)[-1]
            if root is None:
                root = tag
            has_disorders |= tag == "DisorderList"
            if event == "end":
                element.clear()
        if root != "JDBOR" or not has_disorders:
            raise ValueError("Expected Orphadata Product 6 JDBOR/DisorderList XML")
        return
    groups = {
        "clingen_erepo": (("clinvar_variation_id",), ("applied_evidence_codes_met",)),
        "clingen_dosage": (("gene_symbol", "symbol"), ("haploinsufficiency_score", "hi_score"),
                            ("triplosensitivity_score", "ts_score")),
        "clingen_dosage_runtime": (("gene_symbol",), ("haploinsufficiency_score",),
                                    ("triplosensitivity_score",)),
        "g2p_ddg2p": (("gene_symbol",), ("disease_name", "disease_label", "disease"),
                      ("molecular_mechanism", "mutation_consequence")),
    }
    with path.open(encoding="utf-8-sig", newline="") as handle:
        lines = (line for line in handle if not (kind.startswith("clingen_dosage")
                 and line.startswith("#") and "\t" not in line))
        rows = csv.reader(lines, delimiter="," if kind == "g2p_ddg2p" else "\t", strict=True)
        header = next(rows, [])
        exact = {"clingen_erepo": {"ClinVar Variation Id", "Applied Evidence Codes (Met)"},
                 "clingen_dosage_runtime": {"#Gene Symbol", "Genomic Location", "Haploinsufficiency Score",
                     "Haploinsufficiency Description", "Haploinsufficiency Disease ID", "Triplosensitivity Score"}}
        if kind in exact and not exact[kind].issubset(header):
            raise ValueError(f"Downloaded {kind} lacks the original runtime reader's exact columns: {path}")
        normalized = {re.sub(r"[^a-z0-9]+", "_", field.lower()).strip("_") for field in header}
        if any(not normalized.intersection(group) for group in groups[kind]):
            raise ValueError(f"Downloaded {kind} has incompatible source columns: {path}")
        count = 0
        for row in rows:
            if not any(row):
                continue
            if len(row) != len(header):
                raise ValueError(f"Downloaded {kind} contains an incomplete source row: {path}")
            # These two existing readers cast their nonmissing source values
            # to integers. Accept finite numerics, including their float-to-int
            # behavior; do not add biological score or variant-ID rules.
            numeric = {"clingen_erepo": "ClinVar Variation Id", "clingen_dosage_runtime": "Haploinsufficiency Score"}
            if kind in numeric:
                value = row[header.index(numeric[kind])].strip()
                missing = {"", "NaN", "nan", "NA", "N/A", "NULL", "null", "None", "<NA>", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "n/a"}
                if value not in missing:
                    try:
                        valid = math.isfinite(float(value))
                    except ValueError:
                        valid = False
                    if not valid:
                        raise ValueError(f"Downloaded {kind} has a nonnumeric {numeric[kind]}: {path}")
            count += 1
        if not count:
            raise ValueError(f"Downloaded {kind} contains no source rows: {path}")


def source_receipt(path: Path, url: str, *, record: bool = False, receipt: Path | None = None,
                   metadata_only: bool = False) -> bool:
    """Bind completed acquired bytes to their selected URL, without relabeling old data."""
    sidecar = receipt or Path(str(path) + ".source.json")
    if metadata_only:
        try:
            value = json.loads(sidecar.read_text())
        except (OSError, ValueError):
            return False
        return (isinstance(value, dict) and value.get("kind") == "priva_clingen_source_v1"
                and value.get("url") == source_url(url)
                and isinstance(value.get("sha256"), str) and re.fullmatch(r"[0-9a-f]{64}", value["sha256"]) is not None)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    value = {"kind": "priva_clingen_source_v1", "url": source_url(url), "sha256": digest.hexdigest()}
    try:
        previous = json.loads(sidecar.read_text())
    except (OSError, ValueError):
        previous = {}
    if not isinstance(previous, dict):
        previous = {}
    if record:
        value["download_sha256"] = previous.get("download_sha256", value["sha256"]) if previous.get("url") == url else value["sha256"]
        with tempfile.NamedTemporaryFile(mode="w", dir=sidecar.parent, prefix="." + sidecar.name + ".", delete=False) as stream:
            json.dump(value, stream, sort_keys=True)
            stream.write("\n")
            temporary = Path(stream.name)
        temporary.replace(sidecar)
        return True
    return all(previous.get(key) == item for key, item in value.items())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("validate", help="Validate a newly acquired source before publication")
    validate.add_argument("--input", required=True, type=Path)
    validate.add_argument("--kind", required=True, choices=("clingen_erepo", "clingen_dosage", "clingen_dosage_runtime", "g2p_ddg2p", "orphadata"))
    for command in ("record-source", "check-source"):
        receipt = commands.add_parser(command, help="Record/check acquired bytes against the selected source URL")
        receipt.add_argument("--input", required=True, type=Path)
        receipt.add_argument("--url", required=True, type=source_url)
        if command == "check-source":
            receipt.add_argument("--receipt", type=Path)
            receipt.add_argument("--metadata-only", action="store_true",
                                 help="Check owned pending metadata before attempting interrupted publication recovery")
    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            validate_source_payload(args.input, args.kind)
        elif not source_receipt(args.input, args.url, record=args.command == "record-source",
                                receipt=getattr(args, "receipt", None), metadata_only=getattr(args, "metadata_only", False)):
            return 1
    except (ValueError, OSError, csv.Error, ET.ParseError) as error:
        parser.exit(1, f"{error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
