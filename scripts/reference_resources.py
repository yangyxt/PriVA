#!/usr/bin/env python3
"""Prepare pinned UCSC references without replacing supplied scientific inputs.

Only missing standard destinations are downloaded. Supplied FASTAs keep their
sequence bytes, contig names, masking and provenance; their .fai may be repaired.
The original UCSC bigZips files are deliberate pins, not the moving latest/ or
analysisSet releases, and are not substitutes for a selected DRAGEN reference.
"""
from __future__ import annotations

import argparse
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import yaml


BLOCK = 4 * 1024 * 1024
CHR1_LENGTH = {"hg19": 249250621, "hg38": 248956422}
# Publisher manifests and HTTP lengths verified 2026-10-16. These identify the
# exact archives used by the previous installer, including original hg19 chrM.
SOURCES = {
    "hg19": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
        "manifest": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/md5sum.txt",
        "md5": "806c02398f5ac5da8ffd6da2d1d5d1a9", "bytes": 948731419,
        "identity": "UCSC hg19 original bigZips FASTA; original chrM; not latest/analysisSet",
    },
    "hg38": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        "manifest": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt",
        "md5": "1c9dcaddfa41027f17cd8f7a82c7293b", "bytes": 983659424,
        "identity": "UCSC hg38 original bigZips FASTA; not latest/analysisSet",
    },
    "hg19_to_hg38": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        "manifest": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/md5sum.txt",
        "md5": "35887f73fe5e2231656504d1f6430900", "bytes": 227698,
        "identity": "UCSC forward hg19-to-hg38 over.chain; target=hg19, query=hg38",
    },
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def selected_assembly(value):
    assembly = {"hg19": "hg19", "GRCh37": "hg19", "hg38": "hg38", "GRCh38": "hg38"}.get(value)
    require(assembly is not None, "Set assembly to hg19/GRCh37 or hg38/GRCh38; reference filenames do not select an assembly")
    return assembly


def checked_path(value):
    require(isinstance(value, str) and value.strip() == value and value and not any(c in value for c in "\t\r\n"),
            "Reference resource paths must be nonempty path strings without tabs/newlines")
    require(not value.startswith("/paedyl01/"), f"Unresolved old-host reference resource path: {value}; supply its actual local path")
    path = Path(os.path.abspath(value))
    roots = [Path(os.path.abspath(root)) for root in os.environ.get("PRIVA_INSTALL_INSPECT_ROOTS", "").split(os.pathsep) if root]
    if roots:
        require(any(path == root or root in path.parents for root in roots),
                f"Reference resource path outside permitted inspection roots: {value}")
        # Check each symlink before resolving it; do not inspect an outside target.
        cursor = Path(path.anchor)
        parts = list(path.parts[1:])
        links = 0
        while parts:
            cursor /= parts.pop(0)
            if cursor.is_symlink():
                links += 1
                require(links <= 40, f"Too many symlinks in reference resource path: {value}")
                target = Path(os.path.abspath(cursor.parent / os.readlink(cursor)))
                require(any(target == root or root in target.parents for root in roots),
                        f"Reference resource symlink leaves permitted inspection roots: {value}")
                parts = list(target.parts[1:]) + parts
                cursor = Path(target.anchor)
    return path


def present(value):
    return value not in (None, "", "null")


def resource_plan(config, data_dir, scope="gofcards", *, use_cache_root=False):
    require(isinstance(config, dict), "The configuration must be a YAML mapping")
    require(scope in ("selected", "both", "gofcards"), f"Unknown reference scope: {scope}")
    selected = selected_assembly(config.get("assembly"))
    data_dir = checked_path(str(data_dir))
    managed_dir = data_dir
    if use_cache_root and present(config.get("cache_root")):
        root = config["cache_root"]
        require(isinstance(root, str) and re.fullmatch(r"/[A-Za-z0-9_./-]+", root),
                "cache_root must be an absolute path without spaces or shell characters")
        managed_dir = checked_path(root)
        require(managed_dir != Path(managed_dir.anchor), "cache_root must name a dedicated cache directory")
    legacy_defaults = {a: data_dir / "reference" / f"ucsc.{a}.fasta" for a in CHR1_LENGTH}
    defaults = {a: managed_dir / "reference" / f"ucsc.{a}.fasta" for a in CHR1_LENGTH}
    resources = {}
    updates = {}

    def destination(value, old_default, new_default):
        path = checked_path(value if present(value) else str(new_default))
        if path == old_default and path != new_default and not os.path.lexists(path):
            # Redirect only an absent exact managed destination. Do not abandon
            # an interrupted acquisition or a concurrent writer at its old path.
            for suffix in (".fai", ".fai.partial", ".gzi", ".gz", ".gz.partial", ".partial",
                           ".partial.fai", ".partial.fai.partial", ".source.json",
                           ".source.json.partial", ".install.lock"):
                artifact = checked_path(str(path) + suffix)
                require(not os.path.lexists(artifact),
                        f"Legacy managed reference has retained index/source/staging state: {artifact}; "
                        "resolve the selected original destination explicitly before redirecting it")
            path = new_default
        # Existing old defaults keep their publisher/receipt verification;
        # changing cache_root must never relabel them as custom inputs.
        return path, path in (old_default, new_default)

    def add(key, value, assembly):
        path, managed = destination(value, legacy_defaults[assembly], defaults[assembly])
        updates[key] = str(path)
        if path in resources:
            require(resources[path]["assembly"] == assembly,
                    f"GoFCards hg19 and hg38 FASTAs must be distinct assembly-correct inputs: {path}")
        else:
            resources[path] = {"path": path, "assembly": assembly, "key": key,
                               "managed": managed, "kind": "fasta"}
        return path

    chosen = add("ref_genome", config.get("ref_genome"), selected)
    if scope != "selected":
        for assembly in CHR1_LENGTH:
            key = f"gofcards_{assembly}_fasta"
            value = config.get(key)
            if not present(value) and assembly == selected:
                value = str(chosen)
            add(key, value, assembly)
        # Catch distinct symlinks/hardlinks selecting the same wrong assembly.
        paths = list(resources)
        for number, first in enumerate(paths):
            for second in paths[number + 1:]:
                if first.exists() and second.exists() and os.path.samefile(first, second):
                    require(resources[first]["assembly"] == resources[second]["assembly"],
                            "GoFCards hg19 and hg38 FASTAs must be distinct assembly-correct inputs")
    if scope == "gofcards":
        legacy_chain = data_dir / "liftOver" / "hg19ToHg38.over.chain"
        default_chain = managed_dir / "liftOver" / "hg19ToHg38.over.chain"
        key = "gofcards_hg19_to_hg38_chain"
        path, managed = destination(config.get(key), legacy_chain, default_chain)
        require(path not in resources, "The forward-chain destination must be distinct from every FASTA")
        resources[path] = {"path": path, "assembly": "hg19_to_hg38", "key": key,
                           "managed": managed, "kind": "chain"}
        updates[key] = str(path)
    return list(resources.values()), updates


def identity(path):
    stat = path.stat()
    return [stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]


def digest(path, algorithm="sha256"):
    before = identity(path)
    result = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(BLOCK), b""):
            result.update(block)
    require(identity(path) == before, f"Reference resource changed during verification: {path}; wait for its writer")
    return result.hexdigest()


def source_digest(archive, output=None):
    before = identity(archive)
    result = hashlib.sha256()
    with gzip.open(archive, "rb") as handle:
        for block in iter(lambda: handle.read(BLOCK), b""):
            result.update(block)
            if output is not None:
                output.write(block)
    require(identity(archive) == before, f"UCSC source archive changed during decompression: {archive}")
    return result.hexdigest()


def fasta_chr1(path):
    """Count actual chr1 bases when a supplied index needs repair; never infer from its name."""
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    found = False
    length = 0
    with opener(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if found:
                    return length
                found = line[1:].split()[0] in (b"1", b"chr1")
            elif found:
                length += len(line.strip())
    return length if found else None


def valid_index(path, assembly):
    index = checked_path(str(path) + ".fai")
    if not index.is_file() or index.stat().st_mtime_ns < path.stat().st_mtime_ns:
        return False
    try:
        rows = [line.rstrip().split("\t") for line in index.read_text().splitlines()]
        require(rows and all(len(row) == 5 and all(part.isdigit() for part in row[1:]) for row in rows),
                f"Invalid FASTA index: {index}")
        require(len({row[0] for row in rows}) == len(rows), f"Duplicate FASTA index contigs: {index}")
        lengths = {row[0]: int(row[1]) for row in rows}
        return lengths.get("chr1", lengths.get("1")) == CHR1_LENGTH[assembly]
    except (OSError, ValueError):
        return False


def validate_fasta(path, assembly):
    with path.open("rb") as handle:
        first = handle.read(2)
    if first == b"\x1f\x8b":
        with gzip.open(path, "rb") as handle:
            first = handle.read(1)
    require(first.startswith(b">"), f"Not a FASTA reference: {path}; existing bytes retained")
    if not valid_index(path, assembly):
        require(fasta_chr1(path) == CHR1_LENGTH[assembly],
                f"Reference does not identify the required {assembly} chromosome 1 length: {path}; existing bytes retained")


def validate_chain(path):
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8") as handle:
        header = next((line.split() for line in handle if line.startswith("chain ")), [])
    require(len(header) == 13 and header[2].removeprefix("chr") == "1" and header[7].removeprefix("chr") == "1"
            and (header[3], header[8]) == (str(CHR1_LENGTH["hg19"]), str(CHR1_LENGTH["hg38"])),
            f"Expected the forward UCSC hg19-to-hg38 chain: {path}; existing bytes retained")


def receipt_path(path):
    return checked_path(str(path) + ".source.json")


def verify_managed(resource):
    path = resource["path"]
    source = SOURCES[resource["assembly"]]
    receipt = receipt_path(path)
    if receipt.exists():
        metadata = json.loads(receipt.read_text())
        require(metadata.get("source") == source and metadata.get("output") == str(path),
                f"Pinned UCSC source identity mismatch: {receipt}; existing bytes retained")
        require(re.fullmatch(r"[0-9a-f]{64}", str(metadata.get("sha256", ""))), f"Invalid UCSC source receipt: {receipt}")
        require(metadata.get("identity") == identity(path) or metadata["sha256"] == digest(path),
                f"Pinned UCSC reference checksum mismatch: {path}; existing bytes retained")
        return metadata["sha256"]
    # Recover a publication interrupted before its receipt, or an earlier
    # installation whose exact pinned compressed source is still available.
    archive = checked_path(str(path) + ".gz")
    require(archive.is_file() and archive.stat().st_size == source["bytes"] and digest(archive, "md5") == source["md5"],
            f"Existing default reference lacks verified pinned UCSC provenance: {path}. "
            "Select an original reference at its own supplied path or move this file aside explicitly; it will not be replaced.")
    expected = source_digest(archive)
    require(digest(path) == expected, f"Existing default differs from its pinned UCSC source: {path}; existing bytes retained")
    return expected


def preflight(resources):
    """Read-only: check all supplied inputs before any lock, mkdir or download."""
    for resource in resources:
        path = resource["path"]
        if os.path.lexists(path):
            require(path.is_file() and path.stat().st_size > 0, f"Missing, empty, or non-file reference resource: {path}; existing bytes retained")
            if resource["managed"]:
                verify_managed(resource)
            if resource["kind"] == "fasta":
                validate_fasta(path, resource["assembly"])
            else:
                validate_chain(path)
        else:
            require(resource["managed"],
                    f"Explicit reference resource is missing ({resource.get('key', resource['assembly'])}): {path}; "
                    "restore the selected original input. A UCSC download will not replace it.")
            if resource["kind"] == "fasta":
                index = checked_path(str(path) + ".fai")
                require(not os.path.lexists(index), f"Orphan FASTA index exists without its reference: {index}; resolve it explicitly")


def staging(path):
    path = checked_path(str(path))
    require(not path.is_symlink() and (not path.exists() or path.is_file()), f"Reference staging path is not a regular file: {path}")
    return path


def publish_new(partial, destination):
    # link fails if a concurrent process supplied the final file; never replace it.
    os.link(partial, destination)
    partial.unlink()


def download(source, destination):
    destination = staging(destination)
    if destination.exists():
        require(destination.stat().st_size == source["bytes"] and digest(destination, "md5") == source["md5"],
                f"Existing UCSC source archive checksum mismatch: {destination}; existing bytes retained")
        return destination
    partial = staging(str(destination) + ".partial")
    offset = partial.stat().st_size if partial.exists() else 0
    require(offset <= source["bytes"], f"UCSC partial is larger than the pinned source: {partial}; move it aside explicitly")
    if offset < source["bytes"]:
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        try:
            response = urlopen(Request(source["url"], headers=headers), timeout=60)
        except HTTPError as error:
            raise ValueError(f"UCSC download failed: HTTP {error.code}; resumable partial retained: {partial}") from error
        with response:
            status = response.getcode()
            if status == 206:
                content_range = response.headers.get("Content-Range", "")
                require(content_range == f"bytes {offset}-{source['bytes'] - 1}/{source['bytes']}",
                        f"Incorrect UCSC download Content-Range: {content_range}; partial retained")
                mode = "ab"
            else:
                require(status == 200, f"Unexpected UCSC HTTP status: {status}")
                mode = "wb"  # Server ignored Range: restart this staged download.
            print(f"Downloading pinned UCSC resource: {source['url']} -> {partial}", file=sys.stderr)
            with partial.open(mode) as output:
                shutil.copyfileobj(response, output, BLOCK)
    require(partial.stat().st_size == source["bytes"] and digest(partial, "md5") == source["md5"],
            f"UCSC source archive checksum/size mismatch: {partial}; partial retained; move it aside explicitly before restarting")
    publish_new(partial, destination)
    return destination


def index_fasta(path, assembly):
    if valid_index(path, assembly):
        return
    index = checked_path(str(path) + ".fai")
    partial = staging(str(index) + ".partial")
    partial.unlink(missing_ok=True)
    before = identity(path)
    subprocess.run(["samtools", "faidx", "--fai-idx", str(partial), str(path)], check=True)
    require(identity(path) == before, f"Reference changed while indexing: {path}; staged index retained")
    require(valid_index_at(partial, assembly), f"samtools index does not identify required {assembly} chromosome 1: {path}")
    os.replace(partial, index)  # Repair an index only; the FASTA bytes are unchanged.


def valid_index_at(index, assembly):
    with index.open() as handle:
        lengths = {row[0]: int(row[1]) for line in handle if len(row := line.rstrip().split("\t")) == 5}
    return lengths.get("chr1", lengths.get("1")) == CHR1_LENGTH[assembly]


def install_resource(resource):
    path = resource["path"]
    source = SOURCES[resource["assembly"]]
    path.parent.mkdir(parents=True, exist_ok=True)
    lock = staging(str(path) + ".install.lock")
    with lock.open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        preflight([resource])  # Recheck after waiting for a concurrent writer.
        if path.exists():
            output_sha = verify_managed(resource) if resource["managed"] else None
        else:
            archive = download(source, str(path) + ".gz")
            partial = staging(str(path) + ".partial")
            with partial.open("wb") as output:
                output_sha = source_digest(archive, output)
            if resource["kind"] == "fasta":
                staging(str(partial) + ".fai").unlink(missing_ok=True)
                index_fasta(partial, resource["assembly"])
            else:
                validate_chain(partial)
            publish_new(partial, path)
            if resource["kind"] == "fasta":
                publish_new(Path(str(partial) + ".fai"), Path(str(path) + ".fai"))
        if resource["kind"] == "fasta":
            index_fasta(path, resource["assembly"])
        if resource["managed"]:
            receipt = receipt_path(path)
            partial_receipt = staging(str(receipt) + ".partial")
            metadata = {"source": source, "output": str(path), "sha256": output_sha, "identity": identity(path)}
            serialized = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
            if not receipt.exists() or receipt.read_text() != serialized:
                partial_receipt.write_text(serialized)
                os.replace(partial_receipt, receipt)
        print(f"Validated {resource['assembly']} {resource['kind']}: {path}", file=sys.stderr)


def render_config_values(original, updates):
    """Quote path scalars; preserve unrelated YAML text, ordering and comments."""
    if original.lstrip().startswith("{"):
        config = json.loads(original)
        config.update(updates)
        return json.dumps(config, indent=2) + "\n"
    updated = original
    for key, value in updates.items():
        pattern = re.compile(r"^" + re.escape(key) + r"[ \t]*:[^\n]*$", re.MULTILINE)
        matches = list(pattern.finditer(updated))
        require(len(matches) <= 1, f"Duplicate reference configuration key: {key}")
        comment = ""
        if matches:
            scalar = matches[0].group(0).split(":", 1)[1].strip()
            parsed = re.fullmatch(r'''(?:"(?:\\.|[^"\\])*"|'(?:''|[^'])*'|[^#]*?)(\s+#.*)?''', scalar)
            if parsed and parsed.group(1):
                comment = parsed.group(1)
        replacement = key + ": " + json.dumps(value) + comment
        if matches:
            updated = pattern.sub(lambda _: replacement, updated)
        else:
            updated = updated.rstrip("\n") + "\n" + replacement + "\n"
    return updated


def write_config_values(path, original, updated):
    require(path.read_text() == original, f"Configuration changed during reference preparation: {path}; existing configuration retained")
    if updated != original:
        partial = staging(str(path) + ".reference.partial")
        partial.write_text(updated)
        shutil.copymode(path, partial)
        os.replace(partial, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("sources", "preflight", "install"))
    parser.add_argument("config", type=Path, nargs="?")
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parents[1] / "data")
    parser.add_argument("--scope", choices=("selected", "both", "gofcards"), default="gofcards")
    parser.add_argument("--use-cache-root", action="store_true",
                        help="Installer mode: route new managed defaults through configured cache_root; keep existing selections")
    parser.add_argument("--write-config", action="store_true", help="Install only: save prepared paths after every requested resource succeeds")
    args = parser.parse_args()
    if args.command == "sources":
        print(json.dumps(SOURCES, indent=2, sort_keys=True))
        return 0
    try:
        require(args.config is not None, "Supply a configuration file")
        require(not args.write_config or args.command == "install", "--write-config requires install")
        config_path = checked_path(str(args.config))
        original = config_path.read_text()
        resources, updates = resource_plan(yaml.safe_load(original), args.data_dir, args.scope,
                                          use_cache_root=args.use_cache_root)
        updated = render_config_values(original, updates) if args.write_config else original
        preflight(resources)
        if args.command == "install":
            require(shutil.which("samtools") is not None, "samtools is required to install/index references")
            for resource in resources:
                install_resource(resource)
            if args.write_config:
                write_config_values(config_path, original, updated)
            for key, value in updates.items():
                print(f"{key}\t{value}")
    except (OSError, ValueError, EOFError, yaml.YAMLError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
