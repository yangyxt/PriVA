#!/usr/bin/env python3
"""Associate derived DAS metadata with the selected VCF and implementation.

This manages an acquired/derived intermediate, not a repository bundle. VCF
identity uses metadata only; no whole-table hash or scientific validation is
added. Pickle values, insertion order and the original default protocol stay
unchanged. Concurrent builds of the same cache must remain serialized by the
caller, as in main_install.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
from pathlib import Path
import pickle
import stat
import sys
import tempfile


def canonical_assembly(assembly: str) -> str:
    """Treat the installer's two accepted names for each assembly equally."""
    aliases = {"hg19": "hg19", "GRCh37": "hg19", "hg38": "hg38", "GRCh38": "hg38"}
    try:
        return aliases[assembly]
    except KeyError as error:
        raise ValueError(f"Unsupported metadata-cache assembly: {assembly}") from error


def file_identity(path: str | Path, *, source_code: bool = False) -> dict:
    """Describe a selected file without mount-local device numbers."""
    selected = Path(path).absolute()
    info = selected.stat()
    if not stat.S_ISREG(info.st_mode) or info.st_size == 0:
        raise ValueError(f"Missing/nonregular/empty metadata-cache file: {selected}")
    value = {"path": str(selected), "target": str(selected.resolve()),
             "inode": info.st_ino, "size": info.st_size, "mtime_ns": info.st_mtime_ns}
    if source_code:
        value["sha256"] = hashlib.sha256(selected.read_bytes()).hexdigest()
    return value


def input_state(*, vcf_path: str | Path, assembly: str, source_paths) -> dict:
    """Capture the exact selected input and small extraction source files."""
    sources = list(dict.fromkeys(str(Path(path).absolute()) for path in
                                [*source_paths, Path(__file__).resolve()]))
    if len(sources) < 2:
        raise ValueError("The metadata extractor source must be supplied")
    return {"assembly": canonical_assembly(assembly), "vcf": file_identity(vcf_path),
            "implementation": [file_identity(path, source_code=True) for path in sources]}


def receipt_path(cache_path: str | Path) -> Path:
    cache = Path(cache_path)
    return cache.with_name(cache.name + ".install.json")


def association_matches(*, cache_path: str | Path, expected_inputs: dict) -> bool:
    """Check completed publication by metadata; never deserialize the cache here."""
    try:
        expected = {"version": 1, "inputs": expected_inputs, "output": file_identity(cache_path)}
        with receipt_path(cache_path).open() as handle:
            return json.load(handle) == expected
    except (OSError, ValueError):
        return False


def _existing_mode(path: Path):
    try:
        return stat.S_IMODE(path.stat().st_mode)
    except FileNotFoundError:
        return None


def _write_pickle_and_receipt(*, value, cache_path: Path, expected_inputs: dict,
                              vcf_path, assembly, source_paths) -> None:
    """Publish complete data, then its association; leave no false ready pair."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    sidecar = receipt_path(cache_path)
    old_mode = _existing_mode(cache_path)
    receipt_mode = _existing_mode(sidecar)
    with tempfile.TemporaryDirectory(prefix=f".{cache_path.name}.write-", dir=cache_path.parent) as temporary:
        staged = Path(temporary) / cache_path.name
        # Normal open honors umask for a new file; existing final mode is kept.
        with staged.open("wb") as handle:
            pickle.dump(value, handle)
            handle.flush()
            os.fsync(handle.fileno())
        file_identity(staged)
        if old_mode is not None:
            os.chmod(staged, old_mode)
        if input_state(vcf_path=vcf_path, assembly=assembly, source_paths=source_paths) != expected_inputs:
            raise ValueError("Selected VCF or extraction source changed during metadata generation")
        os.replace(staged, cache_path)
        completed = {"version": 1, "inputs": expected_inputs, "output": file_identity(cache_path)}
        staged_receipt = Path(temporary) / sidecar.name
        with staged_receipt.open("w") as handle:
            json.dump(completed, handle, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        if receipt_mode is not None:
            os.chmod(staged_receipt, receipt_mode)
        if input_state(vcf_path=vcf_path, assembly=assembly, source_paths=source_paths) != expected_inputs:
            raise ValueError("Selected VCF or extraction source changed before metadata association")
        os.replace(staged_receipt, sidecar)


def load_or_extract_domain_metadata(*, cache_path: str | Path, vcf_path: str | Path,
                                    assembly: str, extractor, source_paths,
                                    logger=None):
    """Reuse associated metadata or regenerate from this selected VCF only.

    Unassociated or mismatched derived files stay intact until a complete new
    pickle is ready. A failed sidecar publication propagates; its old output
    identity cannot validate the new pickle, so retry regenerates. Ordinary
    pickle errors on an otherwise associated cache also propagate unchanged.
    """
    log = logger if logger is not None else logging.getLogger(__name__)
    cache = Path(cache_path)
    source_paths = tuple(source_paths)
    expected = input_state(vcf_path=vcf_path, assembly=assembly, source_paths=source_paths)
    if association_matches(cache_path=cache, expected_inputs=expected):
        output_before = file_identity(cache)
        log.info(f"[DOMAIN_META] Loading associated metadata from {cache}")
        with cache.open("rb") as handle:
            value = pickle.load(handle)
        if (input_state(vcf_path=vcf_path, assembly=assembly, source_paths=source_paths) != expected
                or file_identity(cache) != output_before
                or not association_matches(cache_path=cache, expected_inputs=expected)):
            raise ValueError("Selected VCF, extraction source or metadata output changed during cache loading")
        # Preserve the existing prot_length regeneration rule exactly.
        sample_meta = next(iter(value.values()), {}) if value else {}
        if sample_meta.get("prot_length") is not None:
            return value
        log.warning("[DOMAIN_META] Stale cache (no prot_length), regenerating...")
    log.info(f"[DOMAIN_META] Extracting from selected AM VCF: {vcf_path}")
    value = extractor()
    _write_pickle_and_receipt(value=value, cache_path=cache, expected_inputs=expected,
                              vcf_path=vcf_path, assembly=assembly, source_paths=source_paths)
    log.info(f"[DOMAIN_META] Cached {len(value)} domains to {cache}")
    return value


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("check", help="Check selected-source association using metadata only")
    check.add_argument("--cache", type=Path, required=True)
    check.add_argument("--vcf", type=Path, required=True)
    check.add_argument("--assembly", required=True)
    check.add_argument("--source-file", type=Path, action="append", required=True)
    args = parser.parse_args(argv)
    try:
        expected = input_state(vcf_path=args.vcf, assembly=args.assembly, source_paths=args.source_file)
        if not association_matches(cache_path=args.cache, expected_inputs=expected):
            raise ValueError("Metadata cache has no matching completed source association")
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
