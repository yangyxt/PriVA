#!/usr/bin/env python3
"""Inspect existing SpliceAI inputs; convert only confirmed native SpliceAI2 tables.

The importable API reads bounded headers, not variant rows or whole-file CRCs.
Fresh installation guidance is assembly-specific; existing compatible v1 VCFs
also remain usable on GRCh38. Neither classification nor registration qualifies
SpliceAI2 for the unchanged ACMG/CADD scientific contracts.
"""
import argparse
import contextlib
import fcntl
import gzip
import json
import os
from pathlib import Path
import re
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
def _completed_projection(snv, indel):
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
    with lock_path.open("a") as lock:
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
        _unchanged_inputs(tracked)
        yield tracked


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
                   str(Path(cache_dir).absolute() / "SpliceAI2" / "GRCh38"), "priva"]
        require(Path(config).resolve(strict=True) == config_path and _identity(config_path) == config_before,
                "Config changed after cache validation")
        _unchanged_inputs(inputs)
        subprocess.run(command, check=True)
    else:
        require(Path(config).resolve(strict=True) == config_path, "Config target changed after cache validation")
        with contextlib.ExitStack() as contexts:
            if report["format"] == "vcf_spliceai2":
                inputs.extend(contexts.enter_context(_completed_projection(snv, indel)))
            _register_vcfs(config, snv, indel, config_before, inputs)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("install")
    for name in ("config", "cache-dir", "snv", "indel"):
        command.add_argument("--" + name, required=True)
    args = parser.parse_args()
    try:
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
