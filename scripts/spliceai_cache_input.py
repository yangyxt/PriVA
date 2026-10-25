#!/usr/bin/env python3
"""Inspect existing SpliceAI inputs; convert only confirmed native SpliceAI2 tables.

Header classification reads bounded headers, not variant rows or whole-file CRCs.
Explicit prepared-cache adoption streams original/copy SHA256 once and preserves
the original receipt; it does not rerun the converter or qualify scientific use.
Fresh installation guidance is assembly-specific; existing compatible v1 VCFs
also remain usable on GRCh38. Neither classification nor registration qualifies
SpliceAI2 for the unchanged ACMG/CADD scientific contracts.
"""
import argparse
import contextlib
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import zlib

MAX_LINE = 64 * 1024
MAX_HEADER = 1024 * 1024
MAX_LINES = 4096
MANIFEST_NAME = "SpliceAI2_GRCh38.downloads-complete.json"
REVISION = "cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
PUBLISHED = {(kind, chrom) for kind in ("snv", "indel") for chrom in CHROMS
             if (kind, chrom) != ("indel", "chrY")}
VCF_COLUMNS = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
CACHE_KEYS = ("spliceai_snv_prescore", "spliceai_indel_prescore")
RELOCATION_KIND = "priva_spliceai2_prepared_relocation_v1"
COPY_BLOCK = 8 * 1024 * 1024


def native_columns():
    columns = ["chrom", "pos", "ref", "alt", "strand"]
    for effect in ("donor_gain", "donor_loss", "acceptor_gain", "acceptor_loss", "jxn_gain", "jxn_loss"):
        fields = ["delta_score", "ref_score", "alt_score"]
        fields += ["donor_dist", "acceptor_dist"] if effect.startswith("jxn") else ["dist"]
        columns.extend(f"{effect}_{field}_{index}" for field in fields for index in range(3))
    return columns + ["spliceai2_summary_score"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical_assembly(assembly):
    aliases = {"hg19": "GRCh37", "grch37": "GRCh37", "hg38": "GRCh38", "grch38": "GRCh38"}
    result = aliases.get(str(assembly).lower())
    require(result is not None, f"Unsupported assembly: {assembly!r}")
    return result


def _identity(path):
    value = path.stat()
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)


def _complete_path(path):
    path = Path(path).absolute()
    def unfinished(part):
        part = part.lower()
        return ".partial" in part or part == ".building" or part.startswith(".building-")

    require(not any(unfinished(part) for part in path.parts), f"Incomplete/staging input path: {path}")
    resolved = path.resolve(strict=True)
    require(not any(unfinished(part) for part in resolved.parts), f"Incomplete/staging input target: {path}")
    require(resolved.is_file() and resolved.stat().st_size > 0, f"Input is not a nonempty regular file: {path}")
    return resolved


def _assembly_in(value):
    return {canonical_assembly(match) for match in
            re.findall(r"(?<![A-Za-z0-9])(?:GRCh37|GRCh38|hg19|hg38)(?![A-Za-z0-9])", value, re.I)}


def _metadata(line):
    match = re.fullmatch(r"##([A-Za-z0-9_.-]+)=(.+)", line)
    require(match is not None, "Invalid VCF metadata line")
    key, value = match.groups()
    # Structured metadata must close its brackets and quoted strings. Values may
    # contain commas inside quotes (INFO descriptions routinely do).
    if value.startswith("<"):
        require(value.endswith(">"), "Truncated structured VCF metadata")
        quoted, escaped = False, False
        for char in value[1:-1]:
            if escaped:
                escaped = False
            elif char == "\\" and quoted:
                escaped = True
            elif char == '"':
                quoted = not quoted
            elif char in "<>" and not quoted:
                raise ValueError("Invalid structured VCF metadata brackets")
        require(not quoted and not escaped, "Unclosed VCF metadata quotation")
    return key, value


def inspect_header(path):
    """Return {'format': ..., 'assembly': ...} or raise ValueError/OSError.

    Gzip/BGZF is detected by magic, independent of filename. Only the bounded
    decompressed header and one lookahead byte are read: corruption encountered
    there fails, but later body corruption needs the downstream full-file check.
    Unmarked legacy-compatible VCF means ``vcf_v1``; scores cannot prove version.
    """
    resolved = _complete_path(path)
    before = _identity(resolved)
    assemblies, versions = set(), set()
    total = 0
    result = None
    try:
        with resolved.open("rb") as raw:
            compressed = raw.read(2) == b"\x1f\x8b"
            raw.seek(0)
            stream = gzip.GzipFile(fileobj=raw) if compressed else raw
            try:
                for number in range(1, MAX_LINES + 1):
                    line = stream.readline(MAX_LINE + 1)
                    require(line, "Missing complete VCF/SpliceAI2 header")
                    total += len(line)
                    require(len(line) <= MAX_LINE and total <= MAX_HEADER, "Cache header exceeds bounded size limit")
                    require(line.endswith(b"\n"), "Truncated cache header line (missing newline)")
                    line = line[:-1]
                    if line.endswith(b"\r"):
                        line = line[:-1]
                    line = line.decode("utf-8")
                    require(not any((ord(char) < 32 and char != "\t") or ord(char) == 127 for char in line), "Control byte in cache header")
                    if number == 1:
                        if line.split("\t") == native_columns():
                            result = {"format": "spliceai2_tsv", "assembly": "GRCh38"}
                            break
                        require(re.fullmatch(r"##fileformat=VCFv4\.[0-5]", line) is not None,
                                "Unsupported cache header: expected VCF or exact 84-column SpliceAI2 TSV")
                        continue
                    if line.startswith("##"):
                        key, value = _metadata(line)
                        lower = key.lower()
                        require(lower != "fileformat", "Duplicate VCF fileformat line")
                        if lower in ("reference", "assembly"):
                            assemblies.update(_assembly_in(value))
                        elif lower == "contig":
                            for assembly in re.findall(r'(?:^|,)assembly=("[^"\n]*"|[^,>]+)', value[1:-1], re.I):
                                assemblies.update(_assembly_in(assembly))
                        if lower.startswith("spliceai2") or (lower == "source" and
                                re.search(r"spliceai(?:2|[ _-]+v?2)(?:\b|_)", value, re.I)):
                            versions.add(2)
                        if lower == "spliceai2version":
                            require(re.fullmatch(r"2(?:\.\d+)*", value) is not None,
                                    "Unsupported explicit SpliceAI2 version")
                        if lower in ("spliceaiversion", "spliceai_version"):
                            require(re.fullmatch(r"[12](?:\.\d+)*", value) is not None,
                                    "Unsupported explicit SpliceAI version")
                            versions.add(int(value[0]))
                        continue
                    columns = line.split("\t")
                    require(columns[:8] == VCF_COLUMNS, "Missing valid complete #CHROM VCF header")
                    require(len(columns) == 8 or (len(columns) >= 10 and columns[8] == "FORMAT" and
                            all(name and not re.search(r"\s", name) for name in columns[9:]) and
                            len(set(columns[9:])) == len(columns[9:])), "Invalid VCF sample columns")
                    require(len(assemblies) <= 1, "Conflicting assembly metadata in VCF header")
                    require(len(versions) <= 1, "Conflicting SpliceAI version metadata")
                    v2 = 2 in versions
                    require(not v2 or assemblies != {"GRCh37"}, "SpliceAI2 VCF conflicts with GRCh37 assembly")
                    result = {"format": "vcf_spliceai2" if v2 else "vcf_v1",
                              "assembly": next(iter(assemblies), "GRCh38" if v2 else None)}
                    break
                require(result is not None, "Cache header exceeds bounded line limit")
                # Also discovers a truncated gzip footer for header-only fixtures;
                # this is not a scan or validation of any variant record.
                stream.read(1)
            finally:
                if compressed:
                    stream.close()
    except (UnicodeError, EOFError, zlib.error) as exc:
        raise ValueError(f"Corrupt/truncated cache header: {resolved}: {exc}") from exc
    require(_identity(resolved) == before and Path(path).resolve(strict=True) == resolved,
            f"Input changed during header inspection: {path}")
    return result


def validate_pair(assembly, snv, indel):
    """Inspect both inputs and enforce matching format/version and assembly.

    Return {'format', 'assembly', 'snv': header_dict, 'indel': header_dict}.
    This does not discover/verify native source manifests or scientific policy.
    """
    assembly = canonical_assembly(assembly)
    headers = {"snv": inspect_header(snv), "indel": inspect_header(indel)}
    require(headers["snv"]["format"] == headers["indel"]["format"], "Mixed SpliceAI cache formats/versions are not supported")
    for kind, header in headers.items():
        require(header["assembly"] in (None, assembly), f"{kind} cache assembly conflicts with {assembly}")
    return {"format": headers["snv"]["format"], "assembly": assembly, **headers}


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"Duplicate manifest key: {key}")
        result[key] = value
    return result


def discover_manifest(snv, indel):
    """Cheap full-47 preflight; full content hashes remain the driver's job."""
    selected = {"snv": _complete_path(snv), "indel": _complete_path(indel)}
    candidates = {path.parent / MANIFEST_NAME for path in selected.values()}
    require(len(candidates) == 1, "Selected native tables must share an adjacent downloads-complete manifest")
    path = _complete_path(candidates.pop())
    before = _identity(path)
    require(before[2] <= MAX_HEADER, "Source manifest exceeds bounded size limit")
    manifest = json.loads(path.read_text(), object_pairs_hook=_unique_object)
    require(isinstance(manifest, dict), "Source manifest must be an object")
    require(type(manifest.get("schema_version")) is int and manifest["schema_version"] == 1 and
            manifest.get("status") == "downloads_complete" and isinstance(manifest.get("confirmed_by"), str) and
            manifest["confirmed_by"].strip(), "A confirmed downloads_complete manifest is required")
    require(manifest.get("resource") == "SpliceAI2" and manifest.get("version") == "2.0" and
            manifest.get("assembly") == "GRCh38" and manifest.get("source_revision") == REVISION,
            "Wrong SpliceAI2 source version/assembly/revision")
    files = manifest.get("files")
    require(isinstance(files, list) and len(files) == 47, "Full manifest requires all 47 published shards")
    roles, paths = set(), {}
    for item in files:
        require(isinstance(item, dict), "Invalid source shard entry")
        require(isinstance(item.get("kind"), str) and isinstance(item.get("chrom"), str), "Invalid source shard role")
        role = (item.get("kind"), item.get("chrom"))
        require(role in PUBLISHED and role not in roles, "Duplicate/unsupported source shard role")
        roles.add(role)
        name = item.get("path")
        require(isinstance(name, str) and Path(name).is_absolute(), "Manifest source paths must be absolute")
        require(Path(name).name == f"SpliceAI2_GRCh38_{role[0]}_{role[1]}.tsv.gz", "Source filename conflicts with manifest role")
        target = _complete_path(name)
        require(target not in paths, "One source path cannot represent multiple shard roles")
        require(type(item.get("bytes")) is int and item["bytes"] > 0 and target.stat().st_size == item["bytes"],
                "Source shard size does not match complete manifest")
        require(isinstance(item.get("sha256"), str) and re.fullmatch(r"[0-9a-f]{64}", item["sha256"]),
                "Source shard needs a pinned SHA256")
        paths[target] = role
    require(roles == PUBLISHED, "Full manifest requires all 47 published shards")
    for kind, target in selected.items():
        require(target in paths and paths[target][0] == kind, f"Selected {kind} input is absent from its manifest role")
    require(_identity(path) == before, "Source manifest changed during inspection")
    return path


@contextlib.contextmanager
def _completed_projection(snv, indel, *, require_existing_lock=False):
    """Cheap publication check, not a repeat of the driver's full SHA checks."""
    selected = {"snv": _complete_path(snv), "indel": _complete_path(indel)}
    folders = {path.parent for path in selected.values()}
    require(len(folders) == 1, "SpliceAI2 VCF pair must share one completed adapter receipt")
    folder = folders.pop()
    receipt = _complete_path(folder / "adapter.ready.json")
    require(receipt == folder / "adapter.ready.json", "Adapter receipt must be a regular adjacent file, not a symlink")
    # Hold the driver's own lock across the two-key registration. A completed
    # receipt is published last; a final VCF alone is not a completion signal.
    lock_path = folder / ".install.lock"
    require(not lock_path.is_symlink(), "Adapter installation lock must not be a symlink")
    if require_existing_lock:
        require(lock_path.is_file(), "Original completed adapter must already contain its installation lock")
    with lock_path.open("r" if require_existing_lock else "a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("SpliceAI2 preparation is active; refusing VCF registration") from exc
        require(not any(path.name == ".building" or path.name.startswith(".building-") for path in folder.iterdir()),
                "SpliceAI2 output has a pending .building directory")
        before = _identity(receipt)
        require(before[2] <= MAX_HEADER, "Adapter receipt exceeds bounded size limit")
        data = json.loads(receipt.read_text(), object_pairs_hook=_unique_object)
        require(isinstance(data, dict), "Adapter receipt must be an object")
        expected = {"schema_version": 1, "resource": "SpliceAI2", "version": "2.0", "assembly": "GRCh38",
                    "status": "prepared_annotation_only", "coverage": "complete_published",
                    "scientific_qualified": False, "source_revision": REVISION}
        require(all(type(data.get(key)) is type(value) and data[key] == value for key, value in expected.items()),
                "SpliceAI2 VCF requires a completed full GRCh38 adapter receipt")
        require(isinstance(data.get("source_manifest_sha256"), str) and
                re.fullmatch(r"[0-9a-f]{64}", data["source_manifest_sha256"]), "Adapter receipt lacks source manifest identity")
        outputs = data.get("outputs")
        require(isinstance(outputs, dict) and set(outputs) == {"snv", "indel"}, "Adapter receipt needs both output roles")
        tracked, seen = [(str(receipt), receipt, before)], set()
        for role, target in selected.items():
            item = outputs[role]
            require(isinstance(item, dict), "Invalid adapter output entry")
            for key, checksum in (("path", "sha256"), ("index_path", "index_sha256")):
                name = item.get(key)
                require(isinstance(name, str) and Path(name).is_absolute(), "Adapter output paths must be absolute")
                path = _complete_path(name)
                require(path.parent == folder and path not in seen, "Adapter outputs must be distinct files beside their receipt")
                require(isinstance(item.get(checksum), str) and re.fullmatch(r"[0-9a-f]{64}", item[checksum]),
                        "Adapter output lacks recorded SHA256")
                if key == "path":
                    require(path == target, f"Selected {role} VCF is not its completed receipt output")
                    require(type(item.get("bytes")) is int and item["bytes"] == path.stat().st_size,
                            "Adapter VCF size does not match receipt")
                else:
                    require(name == item["path"] + ".tbi", "Adapter index does not match its VCF")
                seen.add(path)
                tracked.append((name, path, _identity(path)))
        if "relocation" in data:
            tracked.extend(_relocation_provenance(folder, receipt, data))
        _unchanged_inputs(tracked)
        yield tracked


def _sha256_file(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(COPY_BLOCK), b""):
            value.update(block)
    return value.hexdigest()


def _relocated_receipt(original, destination, relocation):
    value = json.loads(json.dumps(original))
    for item in value["outputs"].values():
        item["path"] = str(destination / Path(item["path"]).name)
        item["index_path"] = str(destination / Path(item["index_path"]).name)
    value["gene_context"]["path"] = str(destination / Path(value["gene_context"]["path"]).name)
    value["relocation"] = relocation
    return value


def _relocation_provenance(folder, receipt, data):
    """Check small relocation records while retaining every original v2 gate."""
    relocation = data["relocation"]
    require(isinstance(relocation, dict) and relocation.get("kind") == RELOCATION_KIND,
            "Unknown prepared SpliceAI2 relocation provenance")
    tracked = []
    related = {}
    for key, basename in (("original_receipt", "adapter.original.ready.json"),
                          ("record", "relocation.ready.json")):
        require(relocation.get(key) == str(folder / basename), "Relocation metadata must be adjacent to the copied cache")
        path = _complete_path(relocation[key])
        require(path == folder / basename and path.stat().st_size <= MAX_HEADER,
                "Relocation metadata must be a bounded regular adjacent file")
        tracked.append((str(path), path, _identity(path)))
        related[key] = path
    original = related["original_receipt"].read_bytes()
    require(hashlib.sha256(original).hexdigest() == relocation.get("original_receipt_sha256"),
            "Copied original adapter receipt identity changed")
    original_data = json.loads(original, object_pairs_hook=_unique_object)
    require(_relocated_receipt(original_data, folder, relocation) == data,
            "Relocated receipt changes more than cache locations and relocation provenance")
    record = json.loads(related["record"].read_text(), object_pairs_hook=_unique_object)
    require(record.get("kind") == RELOCATION_KIND and record.get("status") == "complete"
            and record.get("scientific_qualified") is False
            and record.get("original_receipt_sha256") == relocation["original_receipt_sha256"]
            and record.get("relocated_receipt_sha256") == _sha256_file(receipt),
            "Relocation completion record does not bind both adapter receipts")
    return tracked


def _copy_verified_member(source, target, expected_sha):
    """Resume an owned partial prefix; hash source and read back copied bytes."""
    require(isinstance(expected_sha, str) and re.fullmatch(r"[0-9a-f]{64}", expected_sha),
            "Prepared member lacks its original SHA256")
    before = _identity(source)
    require(not target.is_symlink(), "Refusing symlinked adoption staging member")
    require(not target.exists() or target.stat().st_nlink == 1,
            "Refusing hardlinked adoption staging member")
    partial = target.with_name(target.name + ".copying")
    require(not partial.is_symlink(), "Refusing symlinked adoption partial member")
    require(not partial.exists() or partial.stat().st_nlink == 1,
            "Refusing hardlinked adoption partial member")
    if target.exists():
        require(target.is_file() and target.stat().st_size == before[2]
                and _sha256_file(source) == expected_sha and _sha256_file(target) == expected_sha,
                "Completed adoption staging bytes differ from original receipt")
    else:
        require(not partial.exists() or partial.is_file(), "Invalid adoption partial member")
        value = hashlib.sha256()
        with source.open("rb") as original, partial.open("r+b" if partial.exists() else "x+b") as copied:
            remaining = os.fstat(copied.fileno()).st_size
            require(remaining <= before[2], "Adoption partial is longer than its original input")
            while remaining:
                block = original.read(min(COPY_BLOCK, remaining))
                require(block and copied.read(len(block)) == block,
                        "Retained adoption partial prefix differs from the original input")
                value.update(block)
                remaining -= len(block)
            for block in iter(lambda: original.read(COPY_BLOCK), b""):
                copied.write(block)
                value.update(block)
            copied.flush()
            os.fsync(copied.fileno())
        require(value.hexdigest() == expected_sha, "Original prepared member SHA256 differs from its receipt")
        require(partial.stat().st_size == before[2] and _sha256_file(partial) == expected_sha,
                "Copied prepared member SHA256 differs from its receipt")
        shutil.copystat(source, partial)
        partial.replace(target)
    require(_identity(source) == before, "Original prepared member changed during adoption")
    return {"source": str(source), "destination": str(target), "bytes": before[2],
            "sha256": expected_sha, "source_identity": list(before), "copied_identity": list(_identity(target))}


def _write_adoption_metadata(path, value):
    require(not path.is_symlink() and (not path.exists() or path.stat().st_nlink == 1),
            "Refusing linked adoption metadata")
    descriptor, temporary = tempfile.mkstemp(prefix=".adoption-metadata-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def adopt_completed_projection(*, source_dir, cache_dir, cache_root):
    """Physically adopt a completed v2 publication without touching its originals.

    The original adapter receipt remains byte-exact in the destination. Only
    output/context locations change in the new ready receipt, which explicitly
    records relocation provenance and retains scientific_qualified=False.
    """
    source, selected_root, selected_destination = map(Path, (source_dir, cache_root, cache_dir))
    require(all(path.is_absolute() for path in (source, selected_root, selected_destination)),
            "Prepared cache adoption paths must be absolute")
    source = source.resolve(strict=True)
    require(all(path == Path(os.path.abspath(path)) for path in (selected_root, selected_destination))
            and not any(".partial" in part.lower() or part == ".building" or part.startswith(".building-")
                        for part in selected_destination.parts), "Prepared cache adoption needs complete canonical destination paths")
    require(selected_destination.is_relative_to(selected_root) and selected_destination != selected_root,
            "Prepared cache adoption requires a dedicated destination inside cache_root")
    require(not selected_root.is_symlink() or selected_root.is_dir(), "Selected cache_root alias must name an existing directory")
    # A deliberately selected root alias is supported by the installer. Resolve
    # that boundary once; never follow links created beneath the managed root.
    root = selected_root.resolve()
    require(root != Path(root.anchor), "Prepared cache adoption requires a dedicated cache_root directory")
    destination = root / selected_destination.relative_to(selected_root)
    require(not any(part.is_symlink() for part in (destination, *destination.parents) if part.is_relative_to(root)),
            "Prepared cache adoption must not traverse destination symlinks")
    require(not destination.is_relative_to(source) and not source.is_relative_to(destination),
            "Prepared cache source and destination must be separate directories")
    _complete_path(source / "adapter.ready.json")
    receipt_path = source / "adapter.ready.json"
    require(receipt_path.stat().st_size <= MAX_HEADER, "Original adapter receipt exceeds bounded size limit")
    original_before = _identity(receipt_path)
    original = receipt_path.read_bytes()
    data = json.loads(original, object_pairs_hook=_unique_object)
    require(isinstance(data, dict) and isinstance(data.get("outputs"), dict)
            and set(data["outputs"]) == {"snv", "indel"}
            and all(isinstance(data["outputs"][role], dict)
                    and isinstance(data["outputs"][role].get("path"), str) for role in ("snv", "indel")),
            "Prepared source receipt must declare both VCF output roles")
    snv, indel = (data["outputs"][role]["path"] for role in ("snv", "indel"))
    require(validate_pair("GRCh38", snv, indel)["format"] == "vcf_spliceai2",
            "Prepared adoption requires completed SpliceAI2 VCFs")
    with _completed_projection(snv, indel, require_existing_lock=True) as tracked:
        require(_identity(receipt_path) == original_before, "Original adapter receipt changed before adoption")
        require(Path(snv).resolve().parent == source and Path(indel).resolve().parent == source,
                "Original prepared receipt belongs to a different source directory")
        members = []
        for role in ("snv", "indel"):
            item = data["outputs"][role]
            for key, checksum in (("path", "sha256"), ("index_path", "index_sha256")):
                members.append((_complete_path(item[key]), item[checksum]))
        context = data.get("gene_context", {})
        require(isinstance(context, dict) and isinstance(context.get("path"), str)
                and isinstance(context.get("sha256"), str), "Prepared source receipt lacks gene-context provenance")
        context_path = _complete_path(context["path"])
        require(context_path.parent == source, "Prepared gene context must be beside its original receipt")
        members.append((context_path, context["sha256"]))
        require(len({path.name for path, _ in members}) == 5, "Prepared cache members must have distinct filenames")
        tracked.append((str(context_path), context_path, _identity(context_path)))
        original_sha = hashlib.sha256(original).hexdigest()
        destination.parent.mkdir(parents=True, exist_ok=True)
        with (destination.parent / ("." + destination.name + ".adoption.lock")).open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            if destination.exists():
                require(destination.is_dir(), "Prepared adoption destination is occupied")
                copied = json.loads((destination / "adapter.ready.json").read_text())
                require(copied.get("relocation", {}).get("original_receipt_sha256") == original_sha,
                        "Prepared adoption destination belongs to a different original receipt")
                with _completed_projection(*(copied["outputs"][role]["path"] for role in ("snv", "indel"))):
                    for path, checksum in members:
                        _copy_verified_member(path, destination / path.name, checksum)
                _unchanged_inputs(tracked)
                return {"kind": RELOCATION_KIND, "status": "reused", "destination": str(destination),
                        "original_receipt_sha256": original_sha}
            stage = destination.parent / ("." + destination.name + ".adoption-" + original_sha[:16])
            require(not stage.is_symlink(), "Refusing symlinked adoption workspace")
            stage.mkdir(exist_ok=True)
            owner = stage / "adoption.plan.json"
            planned = {"kind": RELOCATION_KIND, "source": str(source), "destination": str(destination),
                       "original_receipt_sha256": original_sha}
            if owner.exists():
                require(json.loads(owner.read_text()) == planned, "Occupied adoption workspace belongs to another source")
            else:
                require(not any(stage.iterdir()), "Occupied adoption workspace has no owned plan")
                _write_adoption_metadata(owner, (json.dumps(planned, indent=2) + "\n").encode())
            records = [_copy_verified_member(path, stage / path.name, checksum) for path, checksum in members]
            for item in records:
                item["destination"] = str(destination / Path(item["destination"]).name)
            # Small metadata are independent files; never hardlink a receipt.
            _write_adoption_metadata(stage / "adapter.original.ready.json", original)
            relocation = {"kind": RELOCATION_KIND, "original_receipt_sha256": original_sha,
                          "original_receipt": str(destination / "adapter.original.ready.json"),
                          "record": str(destination / "relocation.ready.json")}
            relocated = _relocated_receipt(data, destination, relocation)
            ready_bytes = (json.dumps(relocated, indent=2) + "\n").encode()
            record = dict(planned, status="complete", scientific_qualified=False, members=records,
                          relocated_receipt_sha256=hashlib.sha256(ready_bytes).hexdigest(),
                          original_receipt_preserved=True, converter_executed=False)
            _write_adoption_metadata(stage / "relocation.ready.json", (json.dumps(record, indent=2) + "\n").encode())
            _write_adoption_metadata(stage / "adapter.ready.json", ready_bytes)
            _unchanged_inputs(tracked)
            require(not destination.exists() and not destination.is_symlink(), "Prepared adoption destination became occupied")
            stage.rename(destination)
            with _completed_projection(*(relocated["outputs"][role]["path"] for role in ("snv", "indel"))):
                pass
            return record


def _config(path):
    import yaml  # Existing installer dependency; header inspection is stdlib-only.
    path = _complete_path(path)
    before = _identity(path)
    with path.open(newline="") as stream:
        text = stream.read()
    try:
        node = yaml.compose(text, Loader=yaml.SafeLoader)
        data = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid YAML config: {exc}") from exc
    require(isinstance(node, yaml.MappingNode) and not node.flow_style, "Expected a block-mapping YAML config")
    nodes = {}
    for key, value in node.value:
        require(isinstance(key, yaml.ScalarNode) and key.value not in nodes, "Duplicate/unsupported YAML config key")
        nodes[key.value] = (key, value)
    require(isinstance(data, dict), "Expected a YAML config mapping")
    require(_identity(path) == before, "Config changed during reading")
    return path, before, text, nodes, data


def _register_vcfs(config, snv, indel, expected_identity, inputs):
    """Replace just two scalar values, preserving unrelated YAML and comments."""
    import yaml
    path, before, text, nodes, data = _config(config)
    require(before == expected_identity, "Config changed after cache validation")
    selected = dict(zip(CACHE_KEYS, (str(Path(snv).absolute()), str(Path(indel).absolute()))))
    _unchanged_inputs(inputs)
    if all(data.get(key) == value for key, value in selected.items()):
        return  # Existing VCF reuse must leave even the quoting/mtime untouched.
    changes, additions = [], []
    for key, value in selected.items():
        value = json.dumps(value, ensure_ascii=False)
        if key not in nodes:
            additions.append(f"{key}: {value}\n")
            continue
        key_node, old = nodes[key]
        require(isinstance(old, yaml.ScalarNode) and key_node.start_mark.column == 0 and
                old.start_mark.line == key_node.start_mark.line and old.end_mark.line == key_node.start_mark.line,
                f"Expected a single-line scalar for {key}")
        # Null fields have an empty span immediately after ':'. Ensure spacing
        # without consuming the existing comment or any following line.
        prefix = " " if old.start_mark.index and not text[old.start_mark.index - 1].isspace() else ""
        changes.append((old.start_mark.index, old.end_mark.index, prefix + value))
    for start, end, value in sorted(changes, reverse=True):
        text = text[:start] + value + text[end:]
    if additions:
        text += ("" if text.endswith("\n") else "\n") + "".join(additions)
    # Validate the generated YAML before publishing either key.
    updated = yaml.safe_load(text)
    expected = {**data, **selected}
    require(updated == expected, "Config update would modify unrelated YAML values")
    fd, temporary = tempfile.mkstemp(prefix=".spliceai-config-", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as stream:
            os.fchmod(stream.fileno(), stat.S_IMODE(path.stat().st_mode))
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        require(_identity(path) == before, "Config changed before atomic registration")
        _unchanged_inputs(inputs)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _unchanged_inputs(inputs):
    for requested, resolved, before in inputs:
        require(Path(requested).resolve(strict=True) == resolved and _identity(resolved) == before,
                f"Cache input changed after validation: {requested}")


def _validate_v1_indexes(snv, indel):
    """Bounded index/query checks; neither full-file CRC nor whole-cache identity.

    Keep original licensed VCF/index bytes. Return identities so registration
    also rejects a file/index replaced while validation is in progress.
    """
    try:
        import pysam
    except ImportError as error:
        raise ValueError("SpliceAI index validation requires the bundled pysam package") from error
    inputs = []
    for name in (snv, indel):
        path = _complete_path(name)
        index_name = str(name) + ".tbi"
        try:
            index = _complete_path(index_name)
        except (OSError, ValueError) as error:
            raise ValueError(f"Missing/incomplete SpliceAI tabix index: {index_name}; supply the matching .tbi or run tabix -p vcf on the completed original VCF") from error
        inputs.extend([(name, path, _identity(path)), (index_name, index, _identity(index))])
        try:
            # Read only one original record, then demand its exact indexed
            # retrieval. A readable index alone does not establish pairing.
            with gzip.open(path, "rt") as stream:
                for _ in range(MAX_LINES + 1):
                    line = stream.readline(MAX_LINE + 1)
                    require(len(line) <= MAX_LINE, f"SpliceAI record exceeds validation bound: {path}")
                    if line and not line.startswith("#"):
                        break
                    require(line, f"SpliceAI VCF contains no record for index validation: {path}")
                else:
                    raise ValueError(f"SpliceAI header exceeds validation bound: {path}")
            fields = line.rstrip("\r\n").split("\t")
            require(len(fields) >= 8, f"Malformed first SpliceAI VCF record: {path}")
            position = int(fields[1])
            require(position > 0, f"Invalid first SpliceAI VCF position: {path}")
            with pysam.TabixFile(str(path), index=str(index)) as indexed:
                require(fields[0] in indexed.contigs, f"SpliceAI index lacks first-record contig: {index}")
                found = False
                for number, row in enumerate(indexed.fetch(fields[0], position - 1, position)):
                    require(number < 10000, f"SpliceAI index query exceeds validation bound: {path}")
                    if row == line.rstrip("\r\n"):
                        found = True
                        break
                require(found, f"SpliceAI index does not retrieve the original first record: {index}")
        except (OSError, ValueError, EOFError) as error:
            raise ValueError(f"SpliceAI VCF/index validation failed for {path}: {error}") from error
    _unchanged_inputs(inputs)
    return inputs


def install(config, cache_dir, snv, indel):
    config_path, config_before, _, _, data = _config(config)
    inputs = []
    for name in (snv, indel):
        resolved = _complete_path(name)
        inputs.append((name, resolved, _identity(resolved)))
    report = validate_pair(data.get("assembly"), snv, indel)
    if report["format"] == "spliceai2_tsv":
        manifest = discover_manifest(snv, indel)
        command = ["bash", str(Path(__file__).resolve().with_name("install_utils.sh")), "SpliceAI2_install",
                   str(Path(config).absolute()), str(manifest),
                   str(Path(cache_dir).absolute() / "SpliceAI2" / "GRCh38")]
        require(Path(config).resolve(strict=True) == config_path and _identity(config_path) == config_before,
                "Config changed after cache validation")
        _unchanged_inputs(inputs)
        subprocess.run(command, check=True)
    else:
        require(Path(config).resolve(strict=True) == config_path, "Config target changed after cache validation")
        with contextlib.ExitStack() as contexts:
            if report["format"] == "vcf_spliceai2":
                inputs.extend(contexts.enter_context(_completed_projection(snv, indel)))
            else:
                inputs.extend(_validate_v1_indexes(snv, indel))
            _register_vcfs(config, snv, indel, config_before, inputs)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("install")
    for name in ("config", "cache-dir", "snv", "indel"):
        command.add_argument("--" + name, required=True)
    command = commands.add_parser("adopt-prepared", help="Copy a completed SpliceAI2 cache into managed storage without conversion")
    for name in ("source-dir", "output-dir", "cache-root"):
        command.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "adopt-prepared":
            report = adopt_completed_projection(source_dir=args.source_dir, cache_dir=args.output_dir,
                                                cache_root=args.cache_root)
        else:
            report = install(args.config, args.cache_dir, args.snv, args.indel)
    except subprocess.CalledProcessError as exc:
        print(f"SpliceAI cache preparation failed (exit {exc.returncode})", file=sys.stderr)
        return exc.returncode if exc.returncode > 0 else 1
    except (ValueError, OSError) as exc:
        print(f"SpliceAI cache input rejected: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
