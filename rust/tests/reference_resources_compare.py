#!/usr/bin/env python3
"""Compare `priva reference-resources` with its unchanged Python source.

check --binary PRIVA --repo CHECKOUT --work-root NEW_DIRECTORY [--only CASE ...]
oracle --producer SCRIPT --exception JSON -- ARGUMENTS (internal reference launcher)

Each case builds two identical fixture folders. The Python arm runs
CHECKOUT/scripts/reference_resources.py through the oracle launcher, which
records the class and text of the exception main() reports (or does not
catch) without changing main() or any check. The native arm runs
`PRIVA reference-resources` with the same arguments. Both run with a fixture
`samtools` first on PATH, which writes a .fai from the FASTA it is given (a
`samtools=fail|signal|short|garbage` word in a FASTA header selects a failure).
Compared exactly: exit status, stdout and the `Validated`/`Downloading`
progress lines of stderr (with the case folder written as <ARM>), the
exception category and message, and the complete fixture tree after the run:
every entry's kind, permission bits and content (text with the case folder
written as <ARM>; source receipts with an identity equal to the current
`os.stat` identity of their output written as <CURRENT>).

No case reaches the network: every install case has its resources present,
or fails on a pinned-archive check before a download would start. The
FASTAs are tiny with hand-written .fai indexes carrying the real chromosome
1 lengths; the cases that count bases without an index use one generated
gzip FASTA whose chr1 has the full hg19 length. Declared differences, each
named in its case:
- YAML syntax errors: PyYAML's subclass and text against the native
  YAMLError; only the category family is compared.
- refused profiles: the native command must fail with the named category
  and leave the fixture unchanged (anchors and a flow/JSON configuration).
- argument errors exit 2 in both; argparse's usage text is not compared.
Writes WORK_ROOT/complete.json and exits 1 on any mismatch.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys

sys.dont_write_bytecode = True
PRODUCER_SHA256 = "4a5276d3e812225e6254a1d9cd02d4af6c6c755069cb5f1c15856c50a464c4ec"
YAML_ERRORS = {"YAMLError", "MarkedYAMLError", "ReaderError", "ScannerError", "ParserError",
               "ComposerError", "ConstructorError"}
CHR1 = {"hg19": 249250621, "hg38": 248956422}
PROGRESS = ("Validated ", "Downloading ")
CHAIN_BYTES = 227698

FAKE_SAMTOOLS = r'''#!PYTHON
"""Fixture samtools: only `faidx --fai-idx OUTPUT FASTA`."""
import gzip, os, sys
args = sys.argv[1:]
if len(args) != 4 or args[:2] != ["faidx", "--fai-idx"]:
    sys.exit("fixture samtools: unsupported arguments " + repr(args))
output, fasta = args[2], args[3]
with open(fasta, "rb") as handle:
    opener = gzip.open if handle.read(2) == b"\x1f\x8b" else open
rows, name, length, mode = [], None, 0, ""
with opener(fasta, "rb") as handle:
    for line in handle:
        if line.startswith(b">"):
            if name is not None:
                rows.append((name, length))
            words = line[1:].split()
            name, length = words[0].decode(), 0
            mode = next((word[9:].decode() for word in words[1:] if word.startswith(b"samtools=")), mode)
        else:
            length += len(line.strip())
if name is not None:
    rows.append((name, length))
if mode == "fail":
    sys.exit(3)
if mode == "signal":
    os.kill(os.getpid(), 9)
with open(output, "w") as out:
    for name, length in rows:
        if mode == "garbage":
            out.write(f"{name}\tabc\t1\t1\t1\n")
        else:
            out.write(f"{name}\t{length - (mode == 'short')}\t0\t60\t61\n")
'''


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_json(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")


def oracle(args):
    """Run the unchanged main() and record the exception it reports."""
    module = load(args.producer, "original_reference_resources")
    record = {}
    builtin_print = print

    def observed_print(*values, **options):
        # main() reports a caught failure with print(f"ERROR: {error}") inside
        # its except clause, where sys.exc_info() is that exception.
        if options.get("file") is sys.stderr and values and str(values[0]).startswith("ERROR: "):
            kind, value, _ = sys.exc_info()
            if kind is not None:
                record.update(category=kind.__name__, message=str(value))
        return builtin_print(*values, **options)

    module.print = observed_print
    sys.argv = [str(args.producer), *args.arguments]
    try:
        status = module.main()
    except SystemExit:
        raise
    except BaseException as exc:
        record.update(category=type(exc).__name__, message=str(exc), uncaught=True)
        write_json(args.exception, record)
        raise
    write_json(args.exception, record)
    raise SystemExit(status)


# ---------------------------------------------------------------- fixtures

def fasta(extra=""):
    return f">chr1{extra}\nACGT\n>chr2\nAC\n"


def fai(assembly, length=None):
    return f"chr1\t{length or CHR1[assembly]}\t6\t60\t61\nchr2\t2\t17\t60\t61\n"


def chain(source=CHR1["hg19"], target=CHR1["hg38"], prefix="chr"):
    return f"chain 1 {prefix}1 {source} + 0 1 {prefix}1 {target} + 0 1 1\n1\n"


def big_body():
    """A gzip member holding an hg19-length chr1 in 1 MB lines (no header)."""
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", mtime=0) as handle:
        line = b"A" * 1_000_000 + b"\n"
        whole, rest = divmod(CHR1["hg19"], 1_000_000)
        for _ in range(whole):
            handle.write(line)
        handle.write(b"A" * rest + b"\n>chr2\nAC\n")
    return buffer.getvalue()


def big(header=">chr1"):
    """A two-member gzip FASTA: the header, then the shared chr1 body."""
    return ("big", gzip.compress((header + "\n").encode(), mtime=0))


def base_files():
    return {
        "orig/hg38.fa": fasta(),
        "orig/hg38.fa.fai": fai("hg38"),
        "orig/hg19.fa": fasta(),
        "orig/hg19.fa.fai": fai("hg19"),
        "orig/forward.chain": chain(),
    }


def base_config():
    return {
        "assembly": "hg38",
        "ref_genome": "{arm}/orig/hg38.fa",
        "gofcards_hg19_fasta": "{arm}/orig/hg19.fa",
        "gofcards_hg19_to_hg38_chain": "{arm}/orig/forward.chain",
        "cache_root": "{arm}/cache",
    }


def managed(folder="data", receipts=None):
    """Managed defaults under FOLDER with pinned-source receipts."""
    receipts = receipts or {}
    files = {}
    for assembly in ("hg19", "hg38"):
        name = f"{folder}/reference/ucsc.{assembly}.fasta"
        files[name] = fasta()
        files[name + ".fai"] = fai(assembly)
    files[f"{folder}/liftOver/hg19ToHg38.over.chain"] = chain()
    for name, assembly in ((f"{folder}/reference/ucsc.hg19.fasta", "hg19"),
                           (f"{folder}/reference/ucsc.hg38.fasta", "hg38"),
                           (f"{folder}/liftOver/hg19ToHg38.over.chain", "hg19_to_hg38")):
        files[name + ".source.json"] = ("receipt", name, assembly, receipts.get(assembly, {}))
    return files


def yaml_value(value):
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, list):
        return "[" + ", ".join(value) + "]"
    return json.dumps(value, ensure_ascii=False)


def config_text(case, arm):
    if "raw" in case:
        return case["raw"].replace("{arm}", str(arm))
    config = base_config()
    for key in case.get("unset", ()):
        config.pop(key, None)
    config.update(case.get("set", {}))
    lines = [f"{key}: {yaml_value(value)}" for key, value in config.items()]
    text = "# fixture configuration\n" + "\n".join(lines + case.get("append", [])) + "\n"
    return text.replace("{arm}", str(arm))


def stat_identity(path):
    stat = os.stat(path)
    return [stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]


def receipt(arm, relative, assembly, overrides, sources):
    output = arm / relative
    body = {"source": json.loads(json.dumps(sources[assembly])), "output": str(output),
            "sha256": hashlib.sha256(output.read_bytes()).hexdigest(), "identity": [0, 0, 0, 0, 0]}
    for key, value in overrides.items():
        if value == "<CURRENT>":
            value = stat_identity(output)
        elif isinstance(value, str):
            value = value.replace("{arm}", str(arm))
        target = body
        *parents, last = key.split(".")
        for name in parents:
            target = target[name]
        target[last] = value
    return json.dumps(body, indent=2, sort_keys=True) + "\n"


def materialize(arm, case, sources, big_member):
    arm.mkdir(parents=True)
    files = base_files() if case.get("base", True) else {}
    files.update(case.get("files", {}))
    receipts = []
    for relative, content in files.items():
        path = arm / relative
        if content is None:
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        kind = content[0] if isinstance(content, tuple) else None
        if kind == "link":
            path.symlink_to(content[1].replace("{arm}", str(arm)))
        elif kind == "hardlink":
            os.link(arm / content[1], path)
        elif kind == "directory":
            path.mkdir()
        elif kind == "sparse":
            with path.open("wb") as handle:
                handle.truncate(content[1])
        elif kind == "big":
            path.write_bytes(content[1] + big_member)
        elif kind == "receipt":
            receipts.append((path, content))
        elif isinstance(content, bytes):
            path.write_bytes(content)
        else:
            path.write_text(content.replace("{arm}", str(arm)))
    for path, (_, relative, assembly, overrides) in receipts:
        if isinstance(overrides, str):
            path.write_text(overrides)
        else:
            path.write_text(receipt(arm, relative, assembly, overrides, sources))
    for relative, seconds in case.get("touch", {}).items():
        stat = os.stat(arm / relative)
        os.utime(arm / relative, ns=(stat.st_atime_ns, stat.st_mtime_ns + seconds * 10**9))
    config = arm / case.get("config_rel", "config.yaml")
    config.parent.mkdir(parents=True, exist_ok=True)
    text = config_text(case, arm)
    if case.get("crlf"):
        text = text.replace("\n", "\r\n")
    config.write_bytes(case["config_bytes"] if "config_bytes" in case else text.encode())
    if "config_mode" in case:
        config.chmod(case["config_mode"])
    if case.get("config_directory"):
        config.unlink()
        config.mkdir()


def normalize(text, arm):
    if not isinstance(text, str):
        return text
    for spelling in sorted({str(arm), os.path.realpath(arm)}, key=len, reverse=True):
        text = text.replace(spelling, "<ARM>")
    return text


def receipt_state(data, arm):
    """A receipt's content, its current identity written as <CURRENT>, and
    whether its text is exactly json.dumps(..., indent=2, sort_keys=True)."""
    try:
        body = json.loads(data)
    except ValueError:
        return normalize(data.decode("utf-8", "replace"), arm)
    canonical = data == (json.dumps(body, indent=2, sort_keys=True) + "\n").encode()
    if isinstance(body, dict) and isinstance(body.get("output"), str):
        output = Path(body["output"])
        if output.exists() and body.get("identity") == stat_identity(output):
            body["identity"] = "<CURRENT>"
    return [canonical, normalize(json.dumps(body, sort_keys=True, indent=2), arm)]


def tree(root):
    state = {}
    for folder, directories, names in os.walk(root, followlinks=False):
        for name in sorted(directories + names):
            path = Path(folder) / name
            relative = str(path.relative_to(root))
            if path.is_symlink():
                state[relative] = ["link", normalize(os.readlink(path), root)]
                continue
            mode = oct(path.stat().st_mode & 0o7777)
            if path.is_dir():
                state[relative] = ["directory", mode]
                continue
            data = path.read_bytes()
            if name.endswith(".source.json"):
                state[relative] = ["receipt", receipt_state(data, root), mode]
            elif len(data) <= 65536:
                try:
                    state[relative] = ["text", normalize(data.decode("utf-8"), root), mode]
                except UnicodeDecodeError:
                    state[relative] = ["bytes", hashlib.sha256(data).hexdigest(), mode]
            else:
                state[relative] = ["bytes", hashlib.sha256(data).hexdigest(), mode]
    return state


# ---------------------------------------------------------------- cases

def catalogue():
    gz_fasta = gzip.compress(fasta().encode(), mtime=0)
    cache_flag = ["--use-cache-root"]
    hg19_only = {"set": {"assembly": "hg19", "ref_genome": "{arm}/orig/hg19.fa",
                         "gofcards_hg19_fasta": "{arm}/orig/big.fa.gz",
                         "gofcards_hg38_fasta": "{arm}/orig/hg38.fa"}}
    defaults_only = "# defaults\nassembly: hg38\n"
    cases = {
        # sources
        "sources": {"command": "sources", "config_arg": None, "data_dir": None},
        "sources_ignores_config": {"command": "sources", "config_arg": "absent.yaml", "data_dir": None,
                                   "args": ["--scope", "both"]},
        # preflight of supplied originals
        "originals_gofcards": {},
        "originals_selected": {"args": ["--scope", "selected"]},
        "originals_both": {"args": ["--scope", "both"]},
        "originals_with_cache_root_flag": {"args": cache_flag},
        "alias_grch38": {"set": {"assembly": "GRCh38"}},
        "alias_grch37": {"set": {"assembly": "GRCh37", "ref_genome": "{arm}/orig/hg19.fa",
                                 "gofcards_hg19_fasta": None, "gofcards_hg38_fasta": "{arm}/orig/hg38.fa"}},
        "relative_paths": {"set": {"ref_genome": "orig/../orig/hg38.fa"}, "data_dir": "data"},
        "explicit_hg38_key": {"set": {"gofcards_hg38_fasta": "{arm}/orig/hg38.fa"}},
        "gzip_fasta_indexed": {"set": {"ref_genome": "{arm}/orig/hg38.fa.gz"},
                               "files": {"orig/hg38.fa.gz": gz_fasta, "orig/hg38.fa.gz.fai": fai("hg38")}},
        "index_alias_one": {"files": {"orig/hg38.fa.fai": f"1\t{CHR1['hg38']}\t6\t60\t61\n"}},
        "index_unicode_digits": {"files": {"orig/hg38.fa.fai":
                                           "chr1\t２48956422\t6\t²\t61\r\nchr2\t2\t17\t60\t61\n"}},
        "chain_gzip": {"set": {"gofcards_hg19_to_hg38_chain": "{arm}/orig/forward.chain.gz"},
                       "files": {"orig/forward.chain.gz": gzip.compress(chain().encode(), mtime=0)}},
        "chain_unprefixed_late_header": {"files": {"orig/forward.chain":
                                                   "# note\r\nchainx\n" + chain(prefix="")}},
        # preflight of absent managed defaults
        "defaults_absent_data": {"raw": defaults_only},
        "defaults_absent_selected": {"raw": defaults_only, "args": ["--scope", "selected"]},
        "defaults_absent_cache_root": {"raw": "assembly: hg19\ncache_root: {arm}/cache\n", "args": cache_flag},
        "defaults_empty_and_null_strings": {"raw": "assembly: hg38\nref_genome: ''\ngofcards_hg19_fasta: 'null'\n"},
        "cache_root_ignored_without_flag": {"raw": "assembly: hg38\ncache_root: relative\n"},
        "orphan_default_index": {"raw": defaults_only,
                                 "files": {"data/reference/ucsc.hg19.fasta.fai": fai("hg19")}},
        # managed defaults with receipts
        "managed_receipts": {"raw": defaults_only, "files": managed()},
        "managed_receipt_current_identity": {"raw": defaults_only, "files": managed(receipts={
            "hg38": {"identity": "<CURRENT>", "sha256": "0" * 64}})},
        "managed_receipt_float_bytes": {"raw": defaults_only, "files": managed(receipts={
            "hg19": {"source.bytes": float(948731419)}})},
        "managed_receipt_wrong_source": {"raw": defaults_only, "files": managed(receipts={
            "hg19": {"source.url": "https://example.invalid/hg19.fa.gz"}})},
        "managed_receipt_wrong_output": {"raw": defaults_only, "files": managed(receipts={
            "hg38": {"output": "{arm}/elsewhere.fa"}})},
        "managed_receipt_invalid_sha": {"raw": defaults_only, "files": managed(receipts={
            "hg38": {"sha256": "ABC"}})},
        "managed_receipt_integer_sha": {"raw": defaults_only, "files": managed(receipts={
            "hg38": {"sha256": 12}})},
        "managed_receipt_sha_mismatch": {"raw": defaults_only, "files": managed(receipts={
            "hg19_to_hg38": {"sha256": "f" * 64}})},
        "managed_receipt_not_json": {"raw": defaults_only, "files": managed(receipts={"hg38": "{\n"})},
        "managed_receipt_list": {"raw": defaults_only, "files": managed(receipts={"hg38": "[1]\n"})},
        "managed_receipt_invalid_utf8": {"raw": defaults_only, "files": {
            **managed(), "data/reference/ucsc.hg38.fasta.source.json": b"{\"a\": \"\xff\"}"}},
        "managed_without_receipt": {"raw": defaults_only, "files": {
            **managed(), "data/reference/ucsc.hg19.fasta.source.json": None}},
        "managed_without_receipt_wrong_archive": {"raw": defaults_only, "files": {
            **managed(), "data/reference/ucsc.hg19.fasta.source.json": None,
            "data/reference/ucsc.hg19.fasta.gz": gzip.compress(fasta().encode(), mtime=0)}},
        "managed_receipt_symlinked_default": {"raw": defaults_only, "files": {
            "store/hg38.fa": fasta(), **managed(),
            "data/reference/ucsc.hg38.fasta": ("link", "{arm}/store/hg38.fa")}},
        # legacy defaults and cache_root routing
        "legacy_redirect_absent": {"raw": "assembly: hg38\ncache_root: {arm}/cache\n"
                                          "ref_genome: {arm}/data/reference/ucsc.hg38.fasta\n",
                                   "args": cache_flag},
        "legacy_retained_lock": {"raw": "assembly: hg38\ncache_root: {arm}/cache\n"
                                        "ref_genome: {arm}/data/reference/ucsc.hg38.fasta\n",
                                 "args": cache_flag,
                                 "files": {"data/reference/ucsc.hg38.fasta.install.lock": ""}},
        "legacy_retained_partial_link": {"raw": "assembly: hg38\ncache_root: {arm}/cache\n"
                                                "ref_genome: {arm}/data/reference/ucsc.hg38.fasta\n",
                                         "args": cache_flag,
                                         "files": {"data/reference/ucsc.hg38.fasta.gz.partial":
                                                   ("link", "{arm}/absent")}},
        "legacy_existing_stays_managed": {"raw": "assembly: hg38\ncache_root: {arm}/cache\n"
                                                 "ref_genome: {arm}/data/reference/ucsc.hg38.fasta\n",
                                          "args": cache_flag, "files": managed()},
        "cache_root_relative": {"set": {"cache_root": "relative"}, "args": cache_flag},
        "cache_root_filesystem_root": {"set": {"cache_root": "/"}, "args": cache_flag},
        "cache_root_space": {"set": {"cache_root": "{arm}/with space"}, "args": cache_flag},
        "cache_root_bool": {"set": {"cache_root": True}, "args": cache_flag},
        "cache_root_old_host": {"set": {"cache_root": "/paedyl01/cache"}, "args": cache_flag},
        # configuration values
        "assembly_missing": {"unset": ["assembly"]},
        "assembly_unknown": {"set": {"assembly": "hg20"}},
        "assembly_integer": {"set": {"assembly": 38}},
        "assembly_list": {"set": {"assembly": ["hg38"]}},
        "filename_does_not_select": {"raw": "ref_genome: ucsc.hg19.fasta\n"},
        "path_integer": {"set": {"ref_genome": 5}},
        "path_bool": {"set": {"gofcards_hg19_fasta": True}},
        "path_padded": {"set": {"ref_genome": " {arm}/orig/hg38.fa"}},
        "path_tab": {"set": {"gofcards_hg19_to_hg38_chain": "{arm}/orig/a\tb"}},
        "path_old_host": {"set": {"ref_genome": "/paedyl01/ref/hg38.fa"}},
        "explicit_missing_reference": {"set": {"ref_genome": "{arm}/orig/absent.fa"}},
        "explicit_missing_hg19": {"set": {"gofcards_hg19_fasta": "{arm}/orig/absent.fa"}},
        "explicit_missing_chain": {"set": {"gofcards_hg19_to_hg38_chain": "{arm}/orig/absent.chain"}},
        "same_path_two_assemblies": {"set": {"gofcards_hg19_fasta": "{arm}/orig/hg38.fa"}},
        "hardlink_two_assemblies": {"files": {"orig/hg19.fa": ("hardlink", "orig/hg38.fa")}},
        "symlink_two_assemblies": {"files": {"orig/hg19.fa": ("link", "hg38.fa")}},
        "chain_equals_fasta": {"set": {"gofcards_hg19_to_hg38_chain": "{arm}/orig/hg38.fa"}},
        # FASTA and index checks
        "empty_original": {"files": {"orig/hg38.fa": ""}},
        "directory_original": {"files": {"orig/hg38.fa": ("directory",), "orig/hg38.fa.fai": None}},
        "dangling_original": {"files": {"orig/hg38.fa": ("link", "{arm}/absent.fa")}},
        "not_fasta": {"files": {"orig/hg38.fa": "hello\n", "orig/hg38.fa.fai": fai("hg38")}},
        "bad_gzip_method": {"files": {"orig/hg38.fa": b"\x1f\x8bjunk-after-magic"}},
        "truncated_gzip": {"files": {"orig/hg38.fa": gz_fasta[:15]}},
        "gzip_not_fasta": {"files": {"orig/hg38.fa": gzip.compress(b"hello\n", mtime=0)}},
        "wrong_index_length": {"files": {"orig/hg19.fa.fai": fai("hg19", 1000)}},
        "stale_index": {"touch": {"orig/hg38.fa": 10}},
        "index_four_columns": {"files": {"orig/hg38.fa.fai": f"chr1\t{CHR1['hg38']}\t6\t60\n"}},
        "index_duplicate_names": {"files": {"orig/hg38.fa.fai": fai("hg38") + "chr2\t2\t1\t1\t1\n"}},
        "index_invalid_utf8": {"files": {"orig/hg38.fa.fai": fai("hg38").encode() + b"\xff\n"}},
        "index_superscript_length": {"files": {"orig/hg38.fa.fai": "chr1\t²48956422\t6\t60\t61\n"}},
        "index_missing_counts_bases": {"files": {"orig/hg38.fa.fai": None}},
        "header_without_name": {"files": {"orig/hg38.fa": ">\nACGT\n", "orig/hg38.fa.fai": None}},
        "alias_one_counted": {"files": {"orig/hg38.fa": ">1 desc\nAC GT \r\n", "orig/hg38.fa.fai": None}},
        "chain_reverse": {"files": {"orig/forward.chain": chain(CHR1["hg38"], CHR1["hg19"])}},
        "chain_no_header": {"files": {"orig/forward.chain": "1\n"}},
        "chain_twelve_fields": {"files": {"orig/forward.chain": "chain 1 chr1 249250621 + 0 1 chr1 248956422 + 0 1\n"}},
        "chain_zero_padded_length": {"files": {"orig/forward.chain": chain("0249250621")}},
        "chain_invalid_utf8": {"files": {"orig/forward.chain": b"# \xff\n" + chain().encode()}},
        "chain_bad_gzip": {"set": {"gofcards_hg19_to_hg38_chain": "{arm}/orig/forward.chain.gz"},
                           "files": {"orig/forward.chain.gz": b"\x1f\x8b\x08\x00broken"}},
        # bases counted in a full-length chr1 without an index
        "big_unindexed_preflight": {**hg19_only, "files": {"orig/big.fa.gz": big()}},
        # inspection-root policy
        "policy_inside": {"env": {"PRIVA_INSTALL_INSPECT_ROOTS": "{arm}"}},
        "policy_outside": {"env": {"PRIVA_INSTALL_INSPECT_ROOTS": "{arm}/orig:{arm}/data"},
                           "config_rel": "orig/config.yaml",
                           "set": {"gofcards_hg19_to_hg38_chain": "{arm}/chains/forward.chain"},
                           "files": {"chains/forward.chain": chain()}},
        "policy_config_outside": {"env": {"PRIVA_INSTALL_INSPECT_ROOTS": "{arm}/orig:{arm}/data"}},
        "policy_link_escape": {"env": {"PRIVA_INSTALL_INSPECT_ROOTS": "{arm}/orig:{arm}/data"},
                               "config_rel": "orig/config.yaml",
                               "set": {"gofcards_hg19_fasta": "{arm}/orig/linked.fa"},
                               "files": {"outside/hg19.fa": fasta(), "outside/hg19.fa.fai": fai("hg19"),
                                         "orig/linked.fa": ("link", "{arm}/outside/hg19.fa")}},
        "policy_link_between_roots": {"env": {"PRIVA_INSTALL_INSPECT_ROOTS": "{arm}/orig:{arm}/data:{arm}/outside"},
                                      "config_rel": "orig/config.yaml",
                                      "set": {"gofcards_hg19_fasta": "{arm}/orig/linked.fa"},
                                      "files": {"outside/hg19.fa": fasta(),
                                                "orig/linked.fa": ("link", "../outside/hg19.fa"),
                                                "orig/linked.fa.fai": fai("hg19")}},
        "policy_link_loop": {"env": {"PRIVA_INSTALL_INSPECT_ROOTS": "{arm}/orig:{arm}/data"},
                             "config_rel": "orig/config.yaml",
                             "set": {"ref_genome": "{arm}/orig/a.fa"},
                             "files": {"orig/a.fa": ("link", "b.fa"), "orig/b.fa": ("link", "a.fa")}},
        "policy_link_outside_unchecked_target": {
            "env": {"PRIVA_INSTALL_INSPECT_ROOTS": "{arm}/orig:{arm}/data"}, "config_rel": "orig/config.yaml",
            "set": {"ref_genome": "{arm}/orig/private.fa"},
            "files": {"orig/private.fa": ("link", "/not-permitted/private-reference")}},
        # the configuration file and arguments
        "config_missing": {"config_arg": "{arm}/absent.yaml"},
        "config_directory": {"config_directory": True},
        "config_invalid_utf8": {"config_bytes": b"assembly: hg38\n# \xff\n"},
        "config_list": {"raw": "- a\n- b\n"},
        "config_empty": {"raw": "# nothing\n"},
        "config_braces": {"raw": "{}\n"},
        "config_tab": {"raw": "assembly:\thg38\n\tx: 1\n", "yaml_error": True},
        "config_anchor": {"raw": "assembly: &a hg38\n", "refused": "NotImplementedError"},
        "config_json": {"raw": "{\"assembly\": \"hg38\"}\n", "refused": "NotImplementedError"},
        "config_relative_argument": {"config_arg": "config.yaml"},
        "config_duplicate_keys_last_wins": {"append": ["assembly: hg19", "ref_genome: {arm}/orig/hg19.fa",
                                                       "gofcards_hg38_fasta: {arm}/orig/hg38.fa"]},
        "no_config": {"config_arg": None},
        "write_config_needs_install": {"args": ["--write-config"]},
        "admission_command": {"command": "fetch"},
        "admission_scope": {"args": ["--scope", "all"]},
        "admission_flag_value": {"args": ["--write-config=yes"]},
        "admission_extra_positional": {"args": ["extra"]},
        # install
        "install_originals": {"command": "install"},
        "install_selected": {"command": "install", "args": ["--scope", "selected"]},
        "install_originals_cache_root": {"command": "install", "args": cache_flag},
        "install_write_config": {
            "command": "install", "args": ["--write-config"],
            "raw": "# fixture configuration\nassembly: hg38  # selected\n"
                   "ref_genome: {arm}/orig/hg38.fa   # original\n"
                   "gofcards_hg19_fasta: '{arm}/orig/hg19.fa'\n"
                   "gofcards_hg19_to_hg38_chain: \"{arm}/orig/forward.chain\" # chain\nthreads: 7  # keep\n"},
        "install_write_config_unchanged": {
            "command": "install", "args": ["--write-config"],
            "raw": "assembly: hg38\nref_genome: \"{arm}/orig/hg38.fa\"\n"
                   "gofcards_hg19_fasta: \"{arm}/orig/hg19.fa\"\ngofcards_hg38_fasta: \"{arm}/orig/hg38.fa\"\n"
                   "gofcards_hg19_to_hg38_chain: \"{arm}/orig/forward.chain\"\n"},
        "install_write_config_crlf": {"command": "install", "args": ["--write-config"], "crlf": True},
        "install_write_config_mode": {"command": "install", "args": ["--write-config"], "config_mode": 0o640},
        "install_write_config_unicode": {"command": "install", "args": ["--write-config"],
                                         "set": {"ref_genome": "{arm}/orig/réf.fa"},
                                         "files": {"orig/réf.fa": fasta(), "orig/réf.fa.fai": fai("hg38")}},
        "install_write_config_duplicate_key": {"command": "install", "args": ["--write-config"],
                                               "append": ["ref_genome: {arm}/orig/hg38.fa"]},
        "install_managed_receipts": {"command": "install", "raw": defaults_only, "files": managed()},
        "install_managed_current_receipts": {"command": "install", "raw": defaults_only, "files": managed(receipts={
            name: {"identity": "<CURRENT>"} for name in ("hg19", "hg38", "hg19_to_hg38")})},
        "install_cache_root_redirect": {
            "command": "install", "args": [*cache_flag, "--write-config"], "files": managed("cache"),
            "raw": "assembly: hg38\ncache_root: {arm}/cache\nref_genome: {arm}/data/reference/ucsc.hg38.fasta\n"},
        "install_no_samtools": {"command": "install", "path": "empty"},
        "install_lock_directory": {"command": "install", "files": {"orig/hg19.fa.install.lock": ("directory",)}},
        "install_lock_symlink": {"command": "install", "files": {"orig/hg38.fa.install.lock": ("link", "{arm}/x")}},
        "install_receipt_partial_symlink": {"command": "install", "raw": defaults_only, "files": {
            **managed(), "data/reference/ucsc.hg38.fasta.source.json.partial": ("link", "{arm}/x")}},
        "install_index_partial_symlink": {
            **hg19_only, "command": "install",
            "files": {"orig/big.fa.gz": big(), "orig/big.fa.gz.fai.partial": ("link", "{arm}/x")}},
        "install_existing_wrong_chain_archive": {
            "command": "install", "unset": ["gofcards_hg19_to_hg38_chain", "cache_root"],
            "files": {"data/liftOver/hg19ToHg38.over.chain.gz": "wrong\n"}},
        "install_chain_partial_oversized": {
            "command": "install", "unset": ["gofcards_hg19_to_hg38_chain", "cache_root"],
            "files": {"data/liftOver/hg19ToHg38.over.chain.gz.partial": ("sparse", CHAIN_BYTES + 1)}},
        "install_chain_partial_exact_mismatch": {
            "command": "install", "unset": ["gofcards_hg19_to_hg38_chain", "cache_root"],
            "files": {"data/liftOver/hg19ToHg38.over.chain.gz.partial": ("sparse", CHAIN_BYTES)}},
        "install_chain_archive_symlink": {
            "command": "install", "unset": ["gofcards_hg19_to_hg38_chain", "cache_root"],
            "files": {"data/liftOver/hg19ToHg38.over.chain.gz": ("link", "{arm}/absent.gz")}},
        "install_big_unindexed_repairs": {**hg19_only, "command": "install", "args": ["--write-config"],
                                          "files": {"orig/big.fa.gz": big()}},
        "install_big_stale_index_replaced": {**hg19_only, "command": "install",
                                             "files": {"orig/big.fa.gz.fai": fai("hg19", 5),
                                                       "orig/big.fa.gz": big()}},
        "install_big_samtools_fails": {**hg19_only, "command": "install",
                                       "files": {"orig/big.fa.gz": big(">chr1 samtools=fail")}},
        "install_big_samtools_signal": {**hg19_only, "command": "install",
                                        "files": {"orig/big.fa.gz": big(">chr1 samtools=signal")}},
        "install_big_samtools_short": {**hg19_only, "command": "install",
                                       "files": {"orig/big.fa.gz": big(">chr1 samtools=short")}},
        "install_big_samtools_garbage": {**hg19_only, "command": "install",
                                         "files": {"orig/big.fa.gz": big(">chr1 samtools=garbage")}},
    }
    for case in cases.values():
        case.setdefault("command", "preflight")
    return cases


# ---------------------------------------------------------------- runs

def arguments(case, arm):
    values = [case["command"]]
    config = case.get("config_arg", "{arm}/" + case.get("config_rel", "config.yaml"))
    if config is not None:
        values.append(config)
    data = case.get("data_dir", "{arm}/data")
    if data is not None:
        values += ["--data-dir", data]
    values += case.get("args", [])
    return [value.replace("{arm}", str(arm)) for value in values]


def run(command, folder, label, environment):
    with (folder / (label + ".stdout")).open("wb") as out, (folder / (label + ".stderr")).open("wb") as err:
        process = subprocess.run(list(map(str, command)), cwd=folder / label, stdout=out, stderr=err,
                                 stdin=subprocess.DEVNULL, env=environment, timeout=600)
    return process.returncode


def outcome(folder, label, status, arm):
    stdout = (folder / (label + ".stdout")).read_bytes().decode("utf-8", "replace")
    stderr = (folder / (label + ".stderr")).read_bytes().decode("utf-8", "replace")
    progress = [normalize(line, arm) for line in stderr.splitlines() if line.startswith(PROGRESS)]
    result = {"exit": status, "stdout": normalize(stdout, arm), "progress": progress, "tree": tree(arm)}
    if label == "python":
        recorded = folder / "exception.json"
        value = json.loads(recorded.read_text()) if recorded.exists() else {}
        if value:
            result.update(category=value.get("category"), message=value.get("message"))
    else:
        matched = re.search(r"(?:\A|\n)priva: ([A-Za-z_][A-Za-z_0-9]*): ((?:(?!\npriva: ).)*)\n\Z", stderr, re.S)
        if matched:
            result.update(category=matched[1], message=matched[2])
        elif status not in (0, 2):
            result.update(category=None, message=stderr)
    if "message" in result:
        result["message"] = normalize(result["message"], arm)
    return result


def compare(case, python, native, before):
    if case.get("refused"):
        return native["exit"] == 1 and native.get("category") == case["refused"] and native["tree"] == before
    if python["exit"] == 2 or native["exit"] == 2:
        # Argument admission: argparse usage text differs by design.
        return python["exit"] == native["exit"] == 2 and python["tree"] == native["tree"] == before
    left, right = dict(python), dict(native)
    if case.get("yaml_error") and left.get("category") in YAML_ERRORS and right.get("category") == "YAMLError":
        left["category"] = "YAMLError"
        left.pop("message"), right.pop("message")
    left.pop("uncaught", None)
    return left == right


def check(args):
    producer = args.repo / "scripts/reference_resources.py"
    if hashlib.sha256(producer.read_bytes()).hexdigest() != PRODUCER_SHA256:
        raise SystemExit("reference_resources.py differs from the qualified source; update the pin")
    sources = load(producer, "qualified_reference_resources").SOURCES
    cases = catalogue()
    unknown = set(args.only or ()) - cases.keys()
    if unknown:
        raise SystemExit("Unknown cases: " + ", ".join(sorted(unknown)))
    if not re.fullmatch(r"/[A-Za-z0-9_./-]+", str(args.work_root.absolute())):
        # Cases select cache_root under the case folder; the script admits only these characters.
        raise SystemExit("Use a work root whose absolute path has only letters, digits and _./-")
    args.work_root.mkdir(parents=True, exist_ok=False)
    tools = args.work_root / "tools"
    tools.mkdir()
    (tools / "samtools").write_text(FAKE_SAMTOOLS.replace("#!PYTHON", "#!" + sys.executable, 1))
    (tools / "samtools").chmod(0o755)
    empty = args.work_root / "empty_path"
    empty.mkdir()
    big_member = big_body()
    rows = {}
    for name, case in cases.items():
        if args.only and name not in args.only:
            continue
        folder = args.work_root / name
        folder.mkdir()
        results = {}
        befores = {}
        for label in ("python", "native"):
            arm = (folder / label).absolute()
            materialize(arm, case, sources, big_member)
            befores[label] = tree(arm)
            environment = {key: value for key, value in os.environ.items()
                           if key != "PRIVA_INSTALL_INSPECT_ROOTS"}
            environment.update({key: value.replace("{arm}", str(arm))
                                for key, value in case.get("env", {}).items()})
            environment["PYTHONDONTWRITEBYTECODE"] = "1"
            environment["PATH"] = str(empty) if case.get("path") == "empty" else \
                str(tools) + os.pathsep + environment.get("PATH", "/usr/bin:/bin")
            if label == "python":
                command = [sys.executable, Path(__file__).resolve(), "oracle", "--producer", producer,
                           "--exception", folder / "exception.json", "--", *arguments(case, arm)]
            else:
                command = [args.binary, "reference-resources", *arguments(case, arm)]
            status = run(command, folder, label, environment)
            results[label] = outcome(folder, label, status, arm)
        passed = befores["python"] == befores["native"] and \
            compare(case, results["python"], results["native"], befores["native"])
        rows[name] = {"passed": passed, "python": results["python"], "native": results["native"],
                      "refused": case.get("refused")}
        write_json(folder / "comparison.json", rows[name])
    failed = sorted(name for name, row in rows.items() if not row["passed"])
    report = {
        "producer": str(producer), "producer_sha256": PRODUCER_SHA256,
        "binary": str(args.binary),
        "binary_sha256": hashlib.sha256(Path(args.binary).read_bytes()).hexdigest(),
        "cases": len(rows), "python_successes": sum(row["python"]["exit"] == 0 for row in rows.values()),
        "failed": failed, "rows": rows,
    }
    write_json(args.work_root / "complete.json", report)
    print(json.dumps({key: report[key] for key in ("cases", "python_successes", "failed")}))
    if failed:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="procedure", required=True)
    p = commands.add_parser("check")
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--work-root", type=Path, required=True)
    p.add_argument("--only", action="append")
    p.set_defaults(function=check)
    p = commands.add_parser("oracle")
    p.add_argument("--producer", type=Path, required=True)
    p.add_argument("--exception", type=Path, required=True)
    p.add_argument("arguments", nargs=argparse.REMAINDER)
    p.set_defaults(function=oracle)
    args = parser.parse_args()
    if args.procedure == "oracle" and args.arguments[:1] == ["--"]:
        args.arguments = args.arguments[1:]
    args.function(args)


if __name__ == "__main__":
    main()
