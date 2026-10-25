#!/usr/bin/env python3
"""Install configured CADD annotation/prescore bytes without publishing partials.

This helper prepares caches only. URLs and separate data/index MD5s come from
the caller's configuration; it does not select a model, mirror, or replacement.
"""
import argparse
import fcntl
import hashlib
import http.client
import json
import os
from pathlib import Path, PurePosixPath
import posixpath
import re
import shutil
import sys
import tarfile
import tempfile
from urllib.parse import unquote, urlsplit
from urllib.error import HTTPError
from urllib.request import Request, urlopen


RECEIPT = ".priva-cadd-ready.json"
BLOCK = 8 * 1024 * 1024


def require(condition, message):
    if not condition:
        raise ValueError(message)


def identity(path):
    value = path.lstat()
    return [value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns,
            value.st_ctime_ns, value.st_mode]


def checksum(path, algorithm="md5"):
    before = path.stat()
    digest = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(BLOCK), b""):
            digest.update(block)
    after = path.stat()
    require((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns)
            == (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns),
            f"CADD cache changed during verification: {path}; wait for its writer to finish")
    return digest.hexdigest()


def configured_resources(config, assembly):
    require(isinstance(config, dict), "CADD configuration must be a mapping")
    assembly = {"hg19": "GRCh37", "hg38": "GRCh38"}.get(assembly, assembly)
    require(assembly in ("GRCh37", "GRCh38"), f"Unsupported CADD assembly: {assembly}")
    version = config.get("cadd_version")
    require(isinstance(version, str) and re.fullmatch(r"v[0-9]+(?:\.[0-9]+)*", version),
            "Set cadd_version to the matching CADD resource version (for example v1.7)")
    prefix = f"cadd_{assembly}_"

    def url(key):
        value = config.get(key)
        require(isinstance(value, str) and value.strip() == value and value,
                f"Set {key} to the approved URL matching this assembly/version and checksum")
        parsed = urlsplit(value)
        require(parsed.scheme in ("https", "http", "file") and parsed.path and not parsed.fragment,
                f"Invalid CADD resource URL in {key}")
        return value

    def md5(key):
        value = config.get(key)
        require(isinstance(value, str) and re.fullmatch(r"[0-9a-fA-F]{32}", value),
                f"Set {key} to the approved 32-digit MD5 of that exact resource; "
                "data and .tbi indexes need separate checksums")
        return value.lower()

    annotation = (url(prefix + "anno_url"), md5(prefix + "anno_md5"))
    prescores = []
    for kind in ("snv", "indel"):
        key = prefix + kind + "_anno_"
        source = url(key + "url")
        filename = unquote(PurePosixPath(urlsplit(source).path).name)
        require(re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", filename)
                and ".partial" not in filename.lower(), f"Invalid CADD filename in {key}url")
        data_md5 = md5(key + "md5")
        index_md5 = md5(key + "tbi_md5")
        if key + "tbi_url" in config:
            index_url = url(key + "tbi_url")
        else:
            require(not urlsplit(source).query,
                    f"Set {key}tbi_url explicitly for a data URL containing a query")
            index_url = source + ".tbi"  # Preserve the installer's companion-index convention.
        prescores.extend([(source, filename, data_md5),
                          (index_url, filename + ".tbi", index_md5)])
    require(len({item[1] for item in prescores}) == 4,
            "CADD SNV/indel data and indexes must have distinct destination filenames")
    return {"assembly": assembly, "version": version,
            "annotation": annotation, "prescores": prescores}


def download(source, destination, md5, verified_identity=None):
    """Verify final bytes, or resume a sibling partial and publish after MD5.

    The optional identity is from a prior verified annotation receipt. Existing
    final mismatches are never replaced; failed staged bytes remain inspectable.
    """
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".partial")
    if os.path.lexists(destination):
        require(destination.is_file(), f"Existing CADD cache is not a regular file: {destination}")
        require((not destination.is_symlink() and identity(destination) == verified_identity)
                or checksum(destination) == md5,
                f"CADD checksum mismatch: {destination}; existing bytes retained. "
                "Check the configured source/checksum and move the file aside explicitly before retrying.")
        return destination
    require(not partial.is_symlink() and (not partial.exists() or partial.is_file()),
            f"CADD staging path is not a regular file: {partial}")
    if not partial.exists() or checksum(partial) != md5:
        offset = partial.stat().st_size if partial.exists() else 0
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        try:
            response = urlopen(Request(source, headers=headers), timeout=60)
        except HTTPError as error:
            error.close()
            if error.code == 416 and offset:
                raise ValueError(f"CADD partial does not match the configured MD5 and HTTP resume was refused: {partial}; "
                                 "staged bytes retained. Check source/checksum and move the partial aside explicitly before retrying.") from error
            raise
        with response:
            resumed = response.status == 206
            if resumed:
                content_range = response.headers.get("Content-Range", "")
                require(re.fullmatch(rf"bytes {offset}-[0-9]+/(?:[0-9]+|\*)", content_range),
                        f"Invalid HTTP resume response for {partial}; staged bytes retained")
            # HTTP servers that ignore Range and file:// sources restart safely.
            with partial.open("ab" if resumed else "wb") as output:
                received = 0
                while True:
                    try:
                        block = response.read(BLOCK)
                    except http.client.IncompleteRead as error:
                        output.write(error.partial)
                        raise
                    if not block:
                        break
                    output.write(block)
                    received += len(block)
                length = response.headers.get("Content-Length")
                require(length is None or received == int(length),
                        f"Incomplete CADD transfer: {partial}; staged bytes retained for retry")
        require(checksum(partial) == md5,
                f"CADD checksum mismatch: {partial}; staged bytes retained, no final file published. "
                "Check the configured source/checksum before retrying.")
    # Linking is atomic and refuses to overwrite a concurrently created final.
    os.link(partial, destination)
    partial.unlink()
    print(f"Verified CADD resource: {destination}", flush=True)
    return destination


def _member_path(member, root_name):
    path = PurePosixPath(member.name)
    require(not path.is_absolute() and ".." not in path.parts
            and path.parts and path.parts[0] == root_name,
            f"CADD archive has an unexpected member path: {member.name}")
    return PurePosixPath(*path.parts[1:])


def verify_or_extract(archive_path, target, extract):
    """Compare every legacy member, or extract into an unpublished empty tree."""
    members = {}
    with tarfile.open(archive_path, "r|gz") as archive:
        for member in archive:
            relative = _member_path(member, target.name)
            name = str(relative)
            if name == ".":
                require(member.isdir(), "CADD archive root must be a directory")
                continue
            require(name != RECEIPT, "CADD archive conflicts with its verification receipt")
            path = target.joinpath(*relative.parts)
            require(name not in members, f"Duplicate CADD archive member: {member.name}")
            # Never traverse a previously extracted symlink while writing.
            require(all(not target.joinpath(*parent.parts).is_symlink() for parent in relative.parents),
                    f"CADD annotation member traverses a symlink: {member.name}")
            if member.isdir():
                if extract:
                    path.mkdir(parents=True, exist_ok=True)
                require(path.is_dir() and not path.is_symlink(),
                        f"Incomplete existing CADD annotations: {path}")
                members[name] = {"kind": "directory"}
            elif member.issym():
                linked = posixpath.normpath(posixpath.join(str(relative.parent), member.linkname))
                require(member.linkname and not PurePosixPath(member.linkname).is_absolute()
                        and linked != ".." and not linked.startswith("../"),
                        f"CADD archive symlink escapes its annotation directory: {member.name}")
                if extract:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.symlink_to(member.linkname)
                require(path.is_symlink() and os.readlink(path) == member.linkname,
                        f"Existing CADD annotation link differs from the archive: {path}")
                members[name] = {"kind": "symlink", "link": member.linkname, "identity": identity(path)}
            elif member.isfile() or member.islnk():
                if not extract:
                    require(path.is_file() and not path.is_symlink()
                            and (member.islnk() or path.stat().st_size == member.size),
                            f"Incomplete/different existing CADD annotation: {path}; bytes retained")
                if member.islnk():
                    linked = _member_path(tarfile.TarInfo(member.linkname), target.name)
                    linked_path = target.joinpath(*linked.parts)
                    require(str(linked) in members and members[str(linked)]["kind"] == "file",
                            f"Unsupported CADD archive hardlink: {member.name}")
                    expected = members[str(linked)]["sha256"]
                    if extract:
                        path.parent.mkdir(parents=True, exist_ok=True)
                        os.link(linked_path, path)
                else:
                    digest = hashlib.sha256()
                    with archive.extractfile(member) as source:
                        if extract:
                            path.parent.mkdir(parents=True, exist_ok=True)
                            with path.open("xb") as output:
                                for block in iter(lambda: source.read(BLOCK), b""):
                                    output.write(block)
                                    digest.update(block)
                            path.chmod(member.mode & 0o777)
                        else:
                            for block in iter(lambda: source.read(BLOCK), b""):
                                digest.update(block)
                    expected = digest.hexdigest()
                require(path.is_file() and not path.is_symlink(),
                        f"Incomplete existing CADD annotations: {path}")
                before = identity(path)
                require(extract or checksum(path, "sha256") == expected,
                        f"Existing CADD annotation differs from the verified archive: {path}; bytes retained")
                require(identity(path) == before, f"CADD annotation changed during verification: {path}")
                members[name] = {"kind": "file", "sha256": expected, "identity": before}
            else:
                raise ValueError(f"Unsupported CADD archive member: {member.name}")
    require((target / "vep").is_dir() and any(item["kind"] == "file" for item in members.values()),
            f"Incomplete CADD annotation archive/tree: {target}; expected populated VEP resources")
    # A later archive hardlink can change an earlier inode's ctime.
    for name, entry in members.items():
        if entry["kind"] == "file":
            entry["identity"] = identity(target / name)
    return members


def verify_receipt(target, receipt):
    members = receipt.get("members")
    require(isinstance(members, dict) and members, f"Invalid CADD annotation receipt: {target / RECEIPT}")
    for name, entry in members.items():
        relative = PurePosixPath(name)
        require(not relative.is_absolute() and ".." not in relative.parts,
                f"Invalid member in CADD annotation receipt: {name}")
        require(isinstance(entry, dict), f"Invalid member in CADD annotation receipt: {name}")
        path = target.joinpath(*relative.parts)
        kind = entry.get("kind")
        if kind == "directory":
            require(path.is_dir() and not path.is_symlink(), f"Incomplete existing CADD annotations: {path}")
        elif kind == "symlink":
            require(path.is_symlink() and os.readlink(path) == entry.get("link"),
                    f"Existing CADD annotation link changed: {path}; bytes retained")
            entry["identity"] = identity(path)
        else:
            require(kind == "file" and path.is_file() and not path.is_symlink(),
                    f"Incomplete existing CADD annotations: {path}")
            require(isinstance(entry.get("sha256"), str)
                    and re.fullmatch(r"[0-9a-f]{64}", entry["sha256"]),
                    f"Invalid member checksum in CADD annotation receipt: {name}")
            before = identity(path)
            if before != entry.get("identity"):
                require(checksum(path, "sha256") == entry.get("sha256"),
                        f"Existing CADD annotation changed: {path}; bytes retained")
                require(identity(path) == before, f"CADD annotation changed during verification: {path}")
                entry["identity"] = before
    require((target / "vep").is_dir() and any(item["kind"] == "file" for item in members.values()),
            f"Incomplete existing CADD annotations: {target}")
    return members


def annotation_archive(resources, folder, cache_root=None):
    """Reuse the old cluster archive in place; never bypass a present canonical file."""
    canonical = folder / (resources["assembly"] + "_" + resources["version"] + ".tar.gz")
    if not os.path.lexists(canonical) and cache_root:
        legacy = Path(cache_root) / ".downloads" / canonical.name
        if legacy.is_file():
            return legacy
    return canonical


def annotations(resources, folder, archive_path=None):
    root_name = resources["assembly"] + "_" + resources["version"]
    target = folder / root_name
    archive_path = Path(archive_path) if archive_path is not None else folder / (root_name + ".tar.gz")
    source, md5 = resources["annotation"]
    receipt = None
    if os.path.lexists(target):
        require(target.is_dir() and not target.is_symlink() and (target / "vep").is_dir(),
                f"Incomplete existing CADD annotations: {target}; bytes retained. "
                "Verify or move this directory aside explicitly before retrying.")
    if (target / RECEIPT).exists():
        receipt = json.loads((target / RECEIPT).read_text())
        require(isinstance(receipt, dict), f"Invalid CADD annotation receipt: {target / RECEIPT}")
        if any(receipt.get(key) != value for key, value in
               (("schema", 1), ("assembly", resources["assembly"]), ("version", resources["version"]),
                ("archive_md5", md5), ("archive_url", source))):
            receipt = None
    known_archive = receipt.get("archive_identity") if receipt else None
    download(source, archive_path, md5, known_archive)
    folder.mkdir(parents=True, exist_ok=True)
    if os.path.lexists(target):
        require(target.is_dir() and not target.is_symlink(),
                f"Existing CADD annotation path is not a directory: {target}; bytes retained")
        members = verify_receipt(target, receipt) if receipt else verify_or_extract(archive_path, target, False)
        published = target
        staging = None
    else:
        staging = Path(tempfile.mkdtemp(prefix=".cadd-extract-", dir=folder))
        published = staging / root_name
        published.mkdir()
        try:
            members = verify_or_extract(archive_path, published, True)
        except Exception:
            shutil.rmtree(staging)
            raise
    ready = {"schema": 1, "assembly": resources["assembly"], "version": resources["version"],
             "archive_url": source, "archive_md5": md5,
             "archive_identity": identity(archive_path), "members": members}
    try:
        with tempfile.NamedTemporaryFile(mode="w", prefix=RECEIPT + ".", dir=published, delete=False) as output:
            json.dump(ready, output, indent=2, sort_keys=True)
            output.write("\n")
            receipt_path = Path(output.name)
        receipt_path.replace(published / RECEIPT)
        if staging:
            require(not os.path.lexists(target), f"CADD annotation target appeared during extraction: {target}")
            published.rename(target)
    finally:
        if staging:
            shutil.rmtree(staging)
    print(f"Verified CADD annotations: {target}", flush=True)


def install(config, assembly, base_dir):
    resources = configured_resources(config, assembly)  # Fail before any cache transfer.
    data = Path(base_dir) / "data"
    data.mkdir(parents=True, exist_ok=True)
    with (data / ".priva-cadd.install.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        folder = data / "annotations"
        archive = annotation_archive(resources, folder, config.get("cache_root"))
        annotations(resources, folder, archive_path=archive)
        folder = data / "prescored" / (resources["assembly"] + "_" + resources["version"]) / "incl_anno"
        for source, filename, md5 in resources["prescores"]:
            download(source, folder / filename, md5)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config-json", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--base-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        install(json.loads(args.config_json), args.assembly, args.base_dir)
    except (OSError, ValueError, EOFError, tarfile.TarError, http.client.HTTPException) as error:
        print(f"ERROR: CADD cache installation failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
