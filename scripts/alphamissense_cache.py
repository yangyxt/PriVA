#!/usr/bin/env python3
"""Install the selected original AlphaMissense table into the managed cache.

Examples:
  python alphamissense_cache.py install --assembly hg38 --cache-dir /cache/VEP_plugin_cache
  python alphamissense_cache.py install --assembly hg19 --cache-dir /cache/VEP_plugin_cache --source /downloads/AlphaMissense_hg19.tsv.gz
  python alphamissense_cache.py validate --assembly hg38 --input /cache/VEP_plugin_cache/AlphaMissense/AlphaMissense_hg38.tsv.gz
  python alphamissense_cache.py check-annotation --assembly hg38 --vep-major 116 --input /cache/VEP_plugin_cache/AlphaMissense/AlphaMissense_hg38.vep.vcf.gz

Only external full prescores are handled here. Repository-bundled per-gene mean
scores are neither read nor validated. Downloads use immutable publisher object
generations and stay partial until their original bytes, BGZF stream and rebuilt
tabix index pass. Explicit sources are never modified. This changes no scores,
rows, order, assembly coordinates, scientific filters or VEP annotation policy.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile


# Original publisher objects, confirmed using the public Google Storage metadata
# API on 2026-10-18; also match the existing HGSC acquisition receipts.
SOURCES = {
    "hg19": {"bytes": 622293310, "md5": "3f14fba08c60b09a90a2ae14b7802e4c",
             "generation": "1691073495414540"},
    "hg38": {"bytes": 642961469, "md5": "9fd167735f16a1b87da6eb3e4c25fcb5",
             "generation": "1691073413649109"},
}
HEADER = "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class"
ANNOTATION_HEADER_BYTES = 4 * 1024 * 1024
ANNOTATION_HEADER_LINES = 50000


def require(condition, message):
    if not condition:
        raise ValueError(message)


def assembly_name(value):
    try:
        return {"hg19": "hg19", "GRCh37": "hg19", "hg38": "hg38", "GRCh38": "hg38"}[value]
    except KeyError as error:
        raise ValueError("AlphaMissense requires hg19/GRCh37 or hg38/GRCh38") from error


def file_identity(path):
    """Local change detection, without node-dependent filesystem device IDs."""
    value = path.stat()
    return {"inode": value.st_ino, "bytes": value.st_size,
            "mtime_ns": value.st_mtime_ns, "ctime_ns": value.st_ctime_ns}


class AnnotationHeaderChangedError(ValueError):
    """A concurrent change must never be classified as an obsolete cache."""


def _annotation_header_identity(path):
    def fields(info):
        return tuple(getattr(info, 'st_' + name) for name in
                     ('dev', 'ino', 'mode', 'nlink', 'size', 'mtime_ns', 'ctime_ns'))
    return fields(path.lstat()), fields(path.stat())


def check_vep_annotation(*, input_file, assembly, vep_major, resource):
    """Check only bounded derived-VCF provenance, never raw/bundled payloads.

    Official VEP records its software release and assembly together in one
    ##VEP line. Patch assembly labels retain the same selected coordinate
    family (for example GRCh38.p14 and hg38). Missing or duplicate metadata
    cannot establish reuse. Stop at #CHROM without reading variant records.
    """
    selected_assembly = assembly_name(assembly)
    requested = str(vep_major)
    require(re.fullmatch(r"[1-9][0-9]*", requested), "Requested VEP major must be a positive integer")
    path = Path(input_file).absolute()
    require(path.is_file() and path.stat().st_size > 0, "Missing/empty " + resource + " annotation: " + str(path))
    before = _annotation_header_identity(path)
    # 2026-10-23 CT: check identity on every exit, including malformed metadata
    # and mismatches. A changed read cannot authorize ClinVar cache retirement.
    try:
        return _check_vep_annotation_header(path, assembly, selected_assembly, requested, resource)
    finally:
        try:
            unchanged = _annotation_header_identity(path) == before
        except OSError as error:
            raise AnnotationHeaderChangedError(resource + " annotation changed during its header check") from error
        if not unchanged:
            raise AnnotationHeaderChangedError(resource + " annotation changed during its header check")


def _check_vep_annotation_header(path, assembly, selected_assembly, requested, resource):
    headers = []
    remaining = ANNOTATION_HEADER_BYTES
    with path.open("rb") as raw:
        compressed = raw.read(2) == b"\x1f\x8b"
    with (gzip.open(path, "rb") if compressed else path.open("rb")) as stream:
        for _ in range(ANNOTATION_HEADER_LINES):
            line = stream.readline(remaining + 1)
            remaining -= len(line)
            require(line and remaining >= 0 and line.endswith(b"\n"),
                    resource + " annotation has no complete bounded VCF header")
            if line.startswith(b"##VEP="):
                headers.append(line.rstrip(b"\r\n").decode("utf-8"))
            if line.startswith(b"#CHROM\t"):
                break
            require(line.startswith(b"##"), resource + " annotation has no complete VCF column header")
        else:
            raise ValueError(resource + " annotation VCF header exceeds the bounded line limit")
    require(len(headers) == 1, resource + " annotation requires exactly one unambiguous ##VEP header")
    fields = [token.partition("=") for token in shlex.split(headers[0])]
    versions = [value for key, separator, value in fields if key == "##VEP" and separator]
    assemblies = [value for key, separator, value in fields if key == "assembly" and separator]
    require(len(versions) == 1 and len(assemblies) == 1,
            resource + " annotation requires unambiguous VEP release and assembly metadata")
    version = re.fullmatch(r"v([1-9][0-9]*)(?:\.[0-9]+)*", versions[0])
    require(version is not None and version.group(1) == requested,
            f"{resource} annotation VEP release {versions[0]} differs from requested major {requested}")
    observed_assembly = re.sub(r"^(GRCh(?:37|38))\.p[1-9][0-9]*$", r"\1", assemblies[0])
    require(assembly_name(observed_assembly) == selected_assembly,
            f"{resource} annotation assembly {assemblies[0]} differs from requested {assembly}")
    return {"path": str(path), "assembly": selected_assembly,
            "vep_major": requested, "vep_version": versions[0][1:]}


def check_annotation(*, input_file, assembly, vep_major):
    """Preserve the public AlphaMissense interface and diagnostics."""
    return check_vep_annotation(input_file=input_file, assembly=assembly,
                                vep_major=vep_major, resource="AlphaMissense")


def checksum(path, algorithm="md5"):
    digest = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(command, **kwargs):
    return subprocess.run([str(value) for value in command], check=True, **kwargs)


def first_record(path, assembly):
    """Bounded schema/assembly probe; pinned full-file MD5 establishes all bytes."""
    header_seen = False
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for number, line in enumerate(stream, 1):
            require(number <= 100, "AlphaMissense header/first record was not found within 100 lines")
            line = line.rstrip("\r\n")
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    require(line == HEADER, "Unexpected AlphaMissense raw TSV header")
                    header_seen = True
                continue
            fields = line.split("\t")
            require(header_seen and len(fields) == 10, "AlphaMissense raw TSV has no expected ten-column record")
            require(fields[4] == assembly, "AlphaMissense raw TSV assembly differs from requested assembly")
            require(fields[1].isdigit() and int(fields[1]) > 0, "Invalid AlphaMissense position")
            return line, fields
    raise ValueError("AlphaMissense raw TSV contains no records")


def validate(path, assembly):
    """Verify complete original publisher bytes and mature-tool BGZF integrity."""
    assembly = assembly_name(assembly)
    path = Path(path).absolute()
    require(not any(".partial" in part.lower() for part in path.resolve().parts),
            "An explicit AlphaMissense source must not resolve to a partial file")
    require(path.is_file() and path.stat().st_size > 0, "Missing/empty AlphaMissense source: " + str(path))
    before = file_identity(path)
    expected = SOURCES[assembly]
    require(before["bytes"] == expected["bytes"], "AlphaMissense publisher size mismatch: " + str(path))
    require(checksum(path) == expected["md5"], "AlphaMissense publisher MD5 mismatch: " + str(path))
    run(["bgzip", "-t", path], stdout=subprocess.DEVNULL)
    first_record(path, assembly)
    require(file_identity(path) == before, "AlphaMissense source changed during validation")
    return before


def check_index(path, assembly):
    index = Path(str(path) + ".tbi")
    require(index.is_file() and index.stat().st_size > 0, "Missing/empty AlphaMissense tabix index")
    record, fields = first_record(path, assembly)
    contigs = subprocess.check_output(["tabix", "-l", str(path)], text=True).splitlines()
    require(contigs and fields[0] in contigs, "AlphaMissense index has no source chromosome")
    query = subprocess.check_output(["tabix", str(path), f"{fields[0]}:{fields[1]}-{fields[1]}"], text=True)
    require(record in query.splitlines(), "AlphaMissense index does not recover its first source record")
    return contigs


def preserve(path):
    """Retain replaced managed artifacts; never delete a supplied source."""
    if path.exists() or path.is_symlink():
        suffix = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        retained = path.with_name(path.name + ".invalid-" + suffix)
        path.rename(retained)
        print("Retained previous managed artifact: " + str(retained), file=sys.stderr, flush=True)


def download(path, assembly):
    """Resume an immutable object; retain failed or corrupt attempts for recovery."""
    expected = SOURCES[assembly]
    url = (f"https://storage.googleapis.com/dm_alphamissense/AlphaMissense_{assembly}.tsv.gz"
           f"?generation={expected['generation']}")
    for attempt in range(2):
        if path.exists() and path.stat().st_size >= expected["bytes"]:
            if path.stat().st_size == expected["bytes"] and checksum(path) == expected["md5"]:
                return url
            preserve(path)
        print("Downloading selected AlphaMissense " + assembly + " to " + str(path), file=sys.stderr, flush=True)
        run(["curl", "--fail", "--location", "--retry", "3", "--connect-timeout", "60",
             "--continue-at", "-", "--output", path, url], stdout=sys.stderr)
        if path.is_file() and path.stat().st_size == expected["bytes"] and checksum(path) == expected["md5"]:
            return url
        preserve(path)
    raise ValueError("AlphaMissense download failed publisher byte verification after one clean retry")


def receipt_valid(receipt, path, assembly):
    try:
        data = json.loads(receipt.read_text())
        return (data.get("schema_version") == 1 and data.get("assembly") == assembly
                and data.get("publisher") == SOURCES[assembly]
                and data.get("file_identity") == file_identity(path)
                and data.get("index_identity") == file_identity(Path(str(path) + ".tbi")))
    except (ValueError, OSError, AttributeError):
        return False


def install(*, cache_dir, assembly, source=None):
    """Adopt/download one selected prescore and publish a verified managed index."""
    assembly = assembly_name(assembly)
    cache_dir = Path(cache_dir).absolute()
    require(not any(".partial" in part.lower() or part == ".." for part in cache_dir.parts),
            "AlphaMissense cache directory must be complete and canonical")
    folder = cache_dir / "AlphaMissense"
    folder.mkdir(parents=True, exist_ok=True)
    destination = folder / f"AlphaMissense_{assembly}.tsv.gz"
    index = Path(str(destination) + ".tbi")
    receipt = folder / f"AlphaMissense_{assembly}.ready.json"
    supplied = Path(source).absolute() if source else None
    external = supplied is not None and supplied != destination
    for name in ("bgzip", "tabix"):
        require(shutil.which(name), "Missing required AlphaMissense mature tool: " + name)
    with (folder / f".AlphaMissense_{assembly}.install.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        source_identity = validate(supplied, assembly) if external else None
        # A completed, unchanged receipt avoids rereading the large stream. The
        # current data/index pair still has to pass a real indexed record probe.
        if not destination.is_symlink() and receipt_valid(receipt, destination, assembly):
            try:
                check_index(destination, assembly)
                print("Reusing verified managed AlphaMissense " + assembly, file=sys.stderr, flush=True)
                return destination
            except (ValueError, OSError, subprocess.SubprocessError):
                pass
        existing_valid = False
        if destination.is_file() and not destination.is_symlink():
            try:
                validate(destination, assembly)
                existing_valid = True
            except (ValueError, OSError, subprocess.SubprocessError) as error:
                print("Recovering managed AlphaMissense cache: " + str(error), file=sys.stderr, flush=True)
        original_source = str(supplied) if external else None
        source_url = f"https://storage.googleapis.com/dm_alphamissense/AlphaMissense_{assembly}.tsv.gz"
        # All intermediate bytes and indexes stay under the configured root.
        with tempfile.TemporaryDirectory(prefix=".indexing-", dir=folder) as temporary:
            staged = Path(temporary) / destination.name
            if existing_valid:
                staged.symlink_to(destination.resolve())
            else:
                partial = folder / (destination.name + ".partial")
                if external:
                    preserve(partial)
                    shutil.copyfile(supplied, partial)
                    require(file_identity(supplied) == source_identity, "Supplied AlphaMissense source changed while copying")
                elif destination.is_symlink() and destination.is_file():
                    # A prior compatibility view is an external source, not a
                    # bundled resource. Adopt its verified bytes as a real file.
                    try:
                        validate(destination, assembly)
                    except (ValueError, OSError, subprocess.SubprocessError):
                        # Recover the managed view without changing its old
                        # external target, even when that target is invalid.
                        download(partial, assembly)
                    else:
                        original_source = str(destination.resolve())
                        preserve(partial)
                        shutil.copyfile(destination, partial)
                else:
                    require(shutil.which("curl"), "Missing curl for automatic AlphaMissense acquisition")
                    download(partial, assembly)
                require(partial.stat().st_size == SOURCES[assembly]["bytes"]
                        and checksum(partial) == SOURCES[assembly]["md5"],
                        "Staged AlphaMissense source failed publisher byte verification")
                # Validation rejects partial input names; the staging alias is
                # created only after exact complete publisher bytes are present.
                shutil.copyfile(partial, staged)
                validate(staged, assembly)
            run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", "-S", "1", staged])
            contigs = check_index(staged, assembly)
            if not existing_valid:
                preserve(destination)
                staged.replace(destination)
                partial.unlink()
            preserve(index)
            Path(str(staged) + ".tbi").replace(index)
        data = {"schema_version": 1, "resource": "AlphaMissense", "assembly": assembly,
                "publisher": SOURCES[assembly], "source_url": source_url,
                "source_generation": SOURCES[assembly]["generation"], "supplied_source": original_source,
                "path": str(destination), "file_identity": file_identity(destination),
                "index_identity": file_identity(index), "contigs": contigs,
                "scientific_acceptance": False,
                "scope": "Exact original publisher prescore bytes and working tabix index; no VEP/statistics."}
        with tempfile.NamedTemporaryFile(mode="w", dir=folder, prefix=".ready-", delete=False) as stream:
            staged_receipt = Path(stream.name)
            json.dump(data, stream, indent=2)
            stream.write("\n")
        staged_receipt.replace(receipt)
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    install_parser = commands.add_parser("install", help="Download/adopt one selected assembly into its managed cache")
    install_parser.add_argument("--cache-dir", type=Path, required=True)
    install_parser.add_argument("--assembly", required=True)
    install_parser.add_argument("--source", type=Path)
    validate_parser = commands.add_parser("validate", help="Verify an external original publisher prescore without modifying it")
    validate_parser.add_argument("--input", type=Path, required=True)
    validate_parser.add_argument("--assembly", required=True)
    annotation_parser = commands.add_parser("check-annotation", help="Check the bounded derived VCF header for its selected VEP major and assembly")
    annotation_parser.add_argument("--input", type=Path, required=True)
    annotation_parser.add_argument("--assembly", required=True)
    annotation_parser.add_argument("--vep-major", required=True)
    args = parser.parse_args()
    try:
        if args.command == "install":
            print(install(cache_dir=args.cache_dir, assembly=args.assembly, source=args.source))
        elif args.command == "check-annotation":
            print(json.dumps(check_annotation(input_file=args.input, assembly=args.assembly, vep_major=args.vep_major)))
        else:
            print(json.dumps(validate(args.input, args.assembly)))
    except (OSError, ValueError, EOFError, subprocess.SubprocessError) as error:
        parser.exit(1, "AlphaMissense installation failed: " + str(error) + "\n")


if __name__ == "__main__":
    main()
