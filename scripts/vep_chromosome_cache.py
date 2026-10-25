#!/usr/bin/env python3
"""Inventory and preserve chromosome-wise VEP cache work without publication.

User, 2026-10-23 America/Chicago: "Always remember to parallel these
installations by chromosomes." A later request targets ten CPUs per chromosome
and twelve concurrent chromosomes, remembering those already completed.

An output coordinate is only an inventory hint. `preserve-prefix` checks every
ordered input occurrence before certifying a completed chromosome. It copies
only the completed prefix, then independently rereads that exact byte extent;
the original VEP writer may continue appending beyond it. No command here stops
a job, replaces a canonical cache, or asserts whole-installer qualification.
All new files stay under the supplied cache root. Run substantial work in Slurm.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import stat
import subprocess


# 2026-10-23 CT: preserve the user's original scientific execution while
# recovering chromosome work. These are the installed VEP 116.3 parser and
# effective-command-header sources, not a new policy for skipping variants.
VEP116_OMISSION_SOURCES = {
    "modules/Bio/EnsEMBL/VEP/Parser/VCF.pm": "b80cd16203133b757d63406e488f3ef7adfc453340822f50128cb7a68222048a",
    "Bio/EnsEMBL/IO/Parser/BaseVCF4.pm": "66cd60c39c64d1f55f7797ea8198eddd25216e2de8796f2d02fdd382e92b88d7",
    "modules/Bio/EnsEMBL/VEP/Config.pm": "523ef2348fc0d4e16258f4a053b1cb35ce2dd038bc349430e84c9fafa8bfa6b9",
    "modules/Bio/EnsEMBL/VEP/OutputFactory/VCF.pm": "a3f1371cd2c6fc1758708a12a5c8eab058cc012a7845d9b0f2f9d104ba11958f",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def now():
    return datetime.now(timezone.utc).isoformat()


def identity(path):
    path = Path(path).absolute()
    s = path.lstat()
    require(stat.S_ISREG(s.st_mode), f"Regular file required: {path}")
    return dict(path=str(path), inode=s.st_ino, size=s.st_size, mode=s.st_mode,
                nlink=s.st_nlink, mtime_ns=str(s.st_mtime_ns),
                ctime_ns=str(s.st_ctime_ns), observed_device=s.st_dev)


def same_identity(left, right):
    # Mount device numbers are node-local. Keep every other field exact; this
    # does not alter any historical native-controller preservation policy.
    keys = set(left) - {"observed_device"}
    return keys == set(right) - {"observed_device"} and all(left[k] == right[k] for k in keys)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def owned(path, cache_root):
    path, root = Path(path).absolute(), Path(cache_root).absolute()
    require(path.is_relative_to(root) and ".." not in path.parts, "Destination escapes cache_root")
    for ancestor in [path, *path.parents]:
        require(not ancestor.is_symlink(), f"Managed symlink refused: {ancestor}")
        if ancestor == root:
            break
    return path


def write_new(path, value):
    with Path(path).open("x") as f:
        json.dump(value, f, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")


def open_vcf(path):
    with Path(path).open("rb") as f:
        compressed = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rb") if compressed else Path(path).open("rb")


def read_header(stream, limit=8 * 1024 * 1024):
    result = []
    size = 0
    while True:
        line = stream.readline(limit - size + 1)
        require(line.startswith(b"#") and line.endswith(b"\n"), "Incomplete VCF header")
        size += len(line)
        require(size <= limit, "VCF header exceeds bounded limit")
        result.append(line)
        if line.startswith(b"#CHROM\t"):
            require(len(line.rstrip(b"\n").split(b"\t")) >= 8, "Malformed column header")
            return b"".join(result)


def index_inventory(path, bcftools):
    before = identity(path)
    indexes = [Path(str(path) + suffix) for suffix in (".tbi", ".csi")]
    present = [p for p in indexes if p.is_file() and not p.is_symlink()]
    require(present, "An existing input index is required; inventory never builds one")
    pins = [dict(identity=identity(p), sha256=digest(p)) for p in present]
    p = subprocess.run([str(bcftools), "index", "--stats", str(path)], check=True,
                       text=True, capture_output=True, timeout=30)
    rows = [line.split("\t") for line in p.stdout.splitlines() if line]
    require(rows and all(len(r) == 3 and r[2].isdigit() for r in rows), "Malformed index counts")
    names = [r[0] for r in rows]
    require(len(set(names)) == len(names), "Repeated index contig")
    with open_vcf(path) as f:
        header = read_header(f)
    header_names = [v.decode() for v in re.findall(rb"^##contig=<ID=([^,>]+)", header, re.M)]
    require(len(set(header_names)) == len(header_names), "Repeated header contig")
    require([c for c in header_names if c in names] == names, "Header/index order differs")
    require(same_identity(before, identity(path)), "Input changed during inventory")
    require(all(same_identity(v["identity"], identity(v["identity"]["path"]))
                and digest(v["identity"]["path"]) == v["sha256"] for v in pins), "Input index changed")
    return dict(input=before, indexes=pins, header_sha256=hashlib.sha256(header).hexdigest(),
                contig_order=names, records_by_contig={r[0]: int(r[2]) for r in rows},
                total_records=sum(int(r[2]) for r in rows))


def tail_boundary(path, order, limit=1024 * 1024):
    before = identity(path)
    with Path(path).open("rb") as f:
        raw_header = read_header(f)
        offset = max(len(raw_header), before["size"] - limit)
        f.seek(offset)
        raw = f.read(min(limit, before["size"] - offset))
    if offset > len(raw_header):
        cut = raw.find(b"\n")
        require(cut >= 0, "No complete line in bounded output tail")
        offset += cut + 1
        raw = raw[cut + 1:]
    witness = None
    for line in raw.splitlines(keepends=True):
        if line.endswith(b"\n") and not line.startswith(b"#"):
            fields = line.rstrip(b"\n").split(b"\t")
            require(len(fields) >= 8 and fields[0].decode() in order, "Malformed output tail record")
            witness = dict(contig=fields[0].decode(), position=int(fields[1]),
                           offset=offset, length=len(line), sha256=hashlib.sha256(line).hexdigest())
        offset += len(line)
    require(witness is not None, "No complete output row in bounded tail")
    after = identity(path)
    require(all(before[k] == after[k] for k in ("path", "inode", "mode", "nlink"))
            and after["size"] >= before["size"], "Active output replaced or shrank")
    return dict(before=before, after=after, witness=witness,
                header_sha256=hashlib.sha256(raw_header).hexdigest())


def inventory(args):
    out = owned(args.output_json, args.cache_root)
    require(out.parent.is_dir(), "Inventory parent must already exist")
    indexed = index_inventory(args.input_vcf, args.bcftools)
    boundary = tail_boundary(args.serial_output, indexed["contig_order"])
    rank = indexed["contig_order"].index(boundary["witness"]["contig"])
    candidates = indexed["contig_order"][:rank]
    value = dict(kind="priva-vep-chromosome-inventory-v1", created_utc=now(),
                 helper_source=dict(path=str(Path(__file__).absolute()), sha256=digest(__file__)),
                 producer_job=str(args.producer_job), input_inventory=indexed,
                 serial_output=boundary, candidate_completed_contigs=candidates,
                 validated_completed_contigs=[], active_contig=indexed["contig_order"][rank],
                 remaining_contigs=indexed["contig_order"][rank:],
                 states={c: "observed_serial_prefix_NOT_VALIDATED" if c in candidates
                         else "in_progress" if i == rank else "pending"
                         for i, c in enumerate(indexed["contig_order"])},
                 target_cpus_per_chromosome=args.cpus_per_chromosome,
                 target_max_parallel=args.max_parallel,
                 prefix_body_validated=False, canonical_publication=False,
                 installation_accepted=False, scientific_qualification=False)
    write_new(out, value)
    print(json.dumps(dict(inventory=str(out), sha256=digest(out),
                         candidate_completed_contigs=candidates,
                         active_contig=value["active_contig"],
                         validated_completed_contigs=[])))


def row_invariant(line, annotated=False):
    require(line.endswith(b"\n"), "Incomplete VCF occurrence")
    fields = line[:-1].split(b"\t")
    require(len(fields) >= 8 and fields[1].isdigit() and int(fields[1]) > 0, "Malformed VCF occurrence")
    tokens = [] if fields[7] == b"." else fields[7].split(b";")
    require(all(token and token != b"." for token in tokens),
            "Empty/missing INFO token mixed with a nonmissing INFO list")
    csq = [v for v in tokens if v.split(b"=", 1)[0] == b"CSQ"]
    if annotated:
        require(len(csq) == 1 and csq[0].startswith(b"CSQ=") and len(csq[0]) > 4,
                "Annotated occurrence needs exactly one nonempty CSQ")
        tokens = [v for v in tokens if v.split(b"=", 1)[0] != b"CSQ"]
    else:
        require(not csq, "Source unexpectedly contains existing CSQ")
    return fields[:7] + [b";".join(tokens) if tokens else b"."] + fields[8:]


def vep116_omission_policy(header, source_root=None):
    """Bind the narrow omission rule to pinned code and recorded effective flags.

    Config.pm:693-708 loads config files/environment before :731-749 builds
    full_command; OutputFactory/VCF.pm:281 records it in the VCF header.
    Config.pm:410-425 subsequently lets individual/individual_zyg enable
    allow_non_variant. Without this proof, retain the original strict matcher.
    """
    if source_root is None:
        return dict(kind="strict-no-parser-omissions", nonvariant_omissions_authorized=False)
    pins = []
    for relative, expected in VEP116_OMISSION_SOURCES.items():
        path = Path(source_root) / relative
        before = identity(path)
        require(digest(path) == expected, f"Pinned VEP116 omission source differs: {relative}")
        require(same_identity(before, identity(path)), "VEP116 source changed during pinning")
        pins.append(dict(relative_path=relative, identity=before, sha256=expected))
    versions = re.findall(rb'^##VEP="(v[^"\n]+)"', header, re.M)
    commands = re.findall(rb"^##VEP-command-line='([^\n]*)'\n", header, re.M)
    require(len(versions) == 1 and re.fullmatch(rb"v116(?:\.[0-9]+)*", versions[0]),
            "Exactly one recorded VEP116 version is required for parser omissions")
    require(len(commands) == 1, "Exactly one effective VEP command header is required for parser omissions")
    command = commands[0].decode("ascii")
    tokens = shlex.split(command)
    require(tokens and tokens[0] == "vep", "Malformed effective VEP command header")
    flags = {}
    for flag in ("allow_non_variant", "individual", "individual_zyg"):
        positions = [i for i, token in enumerate(tokens) if token == "--" + flag]
        require(len(positions) <= 1 and not any(token.startswith("--" + flag + "=") for token in tokens),
                "Ambiguous effective nonvariant option")
        if positions:
            i = positions[0] + 1
            flags[flag] = tokens[i] if i < len(tokens) and not tokens[i].startswith("--") else "1"
    # Perl's false string is exactly "0" (or empty/undef). The individual
    # options become arrayrefs before their implication is applied, so their
    # presence is conservatively treated as enabling nonvariant output.
    allow = flags.get("allow_non_variant", "0") not in ("", "0") or any(
        flag in flags for flag in ("individual", "individual_zyg"))
    return dict(kind="pinned-vep116-nonvariant-omissions-v1", source_pins=pins,
                effective_command=command, effective_command_sha256=hashlib.sha256(commands[0]).hexdigest(),
                nonvariant_flags=flags, allow_non_variant=allow,
                nonvariant_omissions_authorized=not allow,
                predicate_source="Parser/VCF.pm:233-266; BaseVCF4.pm:486-490,576-590")


def vep116_nonvariant_omitted(invariant, *, allow_non_variant):
    """Exact early-skip predicate; all other omissions remain mismatches.

    Parser/VCF.pm:233-242 dispatches SVs before :259-266 tests first ALT '.'.
    BaseVCF4.pm drops trailing empty ALT components, parses only the second
    '=' component of INFO values, and uses the last duplicate INFO key.
    dont_skip and reference-validation flags are downstream of this branch.
    """
    if allow_non_variant:
        return False
    ref, raw_alts, raw_info = invariant[3], invariant[4], invariant[7]
    alts = raw_alts.split(b",") if raw_alts not in (b"", b"0") else []
    while alts and not alts[-1]:
        alts.pop()
    if not alts or alts[0] != b".":
        return False
    svtype = None
    for token in raw_info.split(b";"):
        key, *values = token.split(b"=")
        if key == b"SVTYPE":
            svtype = values[0] if values else None
    structural = not re.fullmatch(rb"[ACGT]+", ref + b"".join(alts)) and (
        svtype not in (None, b"", b"0") or
        re.search(rb"[<\[\]][^\*]+[>\]\[]|^\.\w+|\w+\.$", b",".join(alts)) is not None)
    return not structural


def preserve_prefix(args):
    require(digest(args.inventory) == args.inventory_sha256, "Inventory checksum differs")
    prior = json.loads(Path(args.inventory).read_text())
    require(prior["kind"] == "priva-vep-chromosome-inventory-v1", "Wrong inventory kind")
    indexed = prior["input_inventory"]
    source, serial = Path(indexed["input"]["path"]), Path(prior["serial_output"]["after"]["path"])
    require(same_identity(indexed["input"], identity(source)), "Inventoried source changed")
    old = prior["serial_output"]["after"]
    current = identity(serial)
    require(all(old[k] == current[k] for k in ("path", "inode", "mode", "nlink"))
            and current["size"] >= old["size"], "Original output replaced or shrank")
    order = indexed["contig_order"]
    completed = prior["candidate_completed_contigs"]
    require(completed and completed == order[:len(completed)], "A nonempty whole-contig prefix is required")
    root = owned(args.output_dir, args.cache_root)
    require(not os.path.lexists(root), "Fresh preservation directory required")
    root.mkdir(mode=0o700)
    kept, omitted = Counter({c: 0 for c in completed}), Counter({c: 0 for c in completed})
    counts, records = Counter(), []
    prefix_hash = hashlib.sha256()
    omission_hash = hashlib.sha256()
    omission_path = root / "omitted-source-occurrences.jsonl"
    active = None
    writer = None
    try:
        with open_vcf(source) as src, serial.open("rb") as original, omission_path.open("xb") as omission_file:
            source_header, vep_header = read_header(src), read_header(original)
            require(hashlib.sha256(source_header).hexdigest() == indexed["header_sha256"]
                    and hashlib.sha256(vep_header).hexdigest() == prior["serial_output"]["header_sha256"],
                    "Inventoried header differs")
            require(source_header.splitlines()[-1] == vep_header.splitlines()[-1], "VCF column/sample header differs")
            for line in source_header.splitlines(keepends=True):
                if line.startswith((b"##fileformat=", b"##INFO=", b"##contig=")):
                    require(source_header.splitlines(keepends=True).count(line)
                            == vep_header.splitlines(keepends=True).count(line), "Original scientific header changed")
            require(len(re.findall(rb"^##INFO=<ID=CSQ,", vep_header, re.M)) == 1, "Missing/duplicate CSQ schema")
            policy = vep116_omission_policy(vep_header, getattr(args, "vep_source_root", None))
            prefix_hash.update(vep_header)
            previous_rank, previous_pos = -1, 0
            output_ordinal = 0
            for source_ordinal, line in enumerate(src, 1):
                invariant = row_invariant(line)
                chrom = invariant[0].decode()
                require(chrom in indexed["records_by_contig"], "Source contig absent from index")
                rank, position = order.index(chrom), int(invariant[1])
                require(rank >= previous_rank and (rank != previous_rank or position >= previous_pos), "Source order regressed")
                previous_rank, previous_pos = rank, position
                counts[chrom] += 1
                if chrom not in completed:
                    continue
                if chrom != active:
                    if writer is not None:
                        writer.close()
                        records[-1].update(identity=identity(records[-1]["path"]), sha256=output_hash.hexdigest())
                    path = root / f"contig-{rank:03d}.vep.vcf"
                    writer = path.open("xb")
                    output_hash = hashlib.sha256(vep_header)
                    writer.write(vep_header)
                    records.append(dict(contig=chrom, path=str(path), source_ordinal_start=sum(indexed["records_by_contig"][c] for c in order[:rank]) + 1))
                    active = chrom
                if policy["nonvariant_omissions_authorized"] and vep116_nonvariant_omitted(
                        invariant, allow_non_variant=policy["allow_non_variant"]):
                    omission = dict(source_ordinal=source_ordinal, contig=chrom,
                                    source_contig_ordinal=counts[chrom], next_output_ordinal=output_ordinal + 1,
                                    source_row_sha256=hashlib.sha256(line).hexdigest(), source_row_bytes=len(line),
                                    reason="vep116_first_alt_dot_allow_non_variant_false_nonstructural")
                    encoded = (json.dumps(omission, sort_keys=True, allow_nan=False) + "\n").encode()
                    omission_file.write(encoded)
                    omission_hash.update(encoded)
                    omitted[chrom] += 1
                    continue
                annotated = original.readline()
                require(row_invariant(annotated, True) == invariant, "Missing/extra/reordered/changed prefix occurrence")
                writer.write(annotated)
                output_hash.update(annotated)
                prefix_hash.update(annotated)
                kept[chrom] += 1
                output_ordinal += 1
            if writer is not None:
                writer.close()
                writer = None
                records[-1].update(identity=identity(records[-1]["path"]), sha256=output_hash.hexdigest())
            end = original.tell()
            boundary = original.readline()
            future = row_invariant(boundary, True)[0].decode()
            require(future in order[len(completed):], "Extra occurrence after certified prefix")
            witness = dict(offset=end, length=len(boundary), sha256=hashlib.sha256(boundary).hexdigest(), contig=future)
        require(dict(counts) == indexed["records_by_contig"], "Full source/index counts differ")
        require(sum(counts.values()) == indexed["total_records"], "Full source/index total differs")
        require(list(kept) == completed and all(kept[c] + omitted[c] == counts[c] for c in completed),
                "Completed-contig occurrence accounting differs")
        # Independent bounded-extent reread proves the copied prefix stayed
        # byte-identical while the writer was allowed to append elsewhere.
        h, remaining = hashlib.sha256(), end
        with serial.open("rb") as f:
            while remaining:
                block = f.read(min(1024 * 1024, remaining))
                require(block, "Serial output shrank during preservation")
                h.update(block)
                remaining -= len(block)
            require(f.read(witness["length"]) == boundary, "Boundary occurrence changed")
        require(h.hexdigest() == prefix_hash.hexdigest(), "Serial prefix changed during preservation")
        require(same_identity(indexed["input"], identity(source)), "Source changed during preservation")
        after = identity(serial)
        require(all(current[k] == after[k] for k in ("path", "inode", "mode", "nlink"))
                and after["size"] >= current["size"], "Serial output replaced or shrank")
        for record in records:
            record["records"] = kept[record["contig"]]
            record["source_records"] = counts[record["contig"]]
            record["omitted_records"] = omitted[record["contig"]]
            record["source_ordinal_end"] = record["source_ordinal_start"] + record["source_records"] - 1
            # Verify the closed destination bytes, independently of the bytes
            # supplied to write(); do not certify a corrupt/changed copy.
            copy_before = identity(record["path"])
            require(same_identity(copy_before, record["identity"]), "Preserved copy identity changed")
            require(digest(record["path"]) == record["sha256"], "Preserved chromosome copy checksum differs")
            require(same_identity(copy_before, identity(record["path"])), "Preserved copy changed during checksum")
        omission_before = identity(omission_path)
        require(digest(omission_path) == omission_hash.hexdigest(), "Omitted-occurrence ledger checksum differs")
        require(same_identity(omission_before, identity(omission_path)), "Omitted-occurrence ledger changed during checksum")
        value = dict(kind="priva-vep-chromosome-preserved-prefix-v1", created_utc=now(),
                     helper_source=dict(path=str(Path(__file__).absolute()), sha256=digest(__file__)),
                     inventory=str(Path(args.inventory).absolute()), inventory_sha256=args.inventory_sha256,
                     input_inventory=indexed, serial_output_before=current, serial_output_after=after,
                     serial_prefix_bytes=end, serial_prefix_sha256=h.hexdigest(), boundary=witness,
                     validated_completed_contigs=completed, records_by_contig=dict(kept), files=records,
                     source_records_by_contig={c: counts[c] for c in completed}, omitted_records_by_contig=dict(omitted),
                     omitted_occurrences=dict(identity=omission_before, sha256=omission_hash.hexdigest(),
                                              records=sum(omitted.values())), parser_omission_policy=policy,
                     source_ordinal_occurrences_preserved=not any(omitted.values()),
                     source_ordinal_occurrences_accounted=True, non_vep_fields_and_info_exact=True,
                     annotations_preserved_without_recalculation=True, serial_writer_may_still_be_running=True,
                     canonical_publication=False, installation_accepted=False, scientific_qualification=False)
        proof = root / "preserved-prefix.json"
        write_new(proof, value)
        print(json.dumps(dict(proof=str(proof), sha256=digest(proof), validated_completed_contigs=completed,
                             records=sum(kept.values()), canonical_publication=False)))
    finally:
        if writer is not None:
            writer.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    p = sub.add_parser("inventory", help="Persist provisional chromosome progress from index and bounded output")
    p.add_argument("--input-vcf", required=True, type=Path)
    p.add_argument("--serial-output", required=True, type=Path)
    p.add_argument("--output-json", required=True, type=Path)
    p.add_argument("--cache-root", required=True, type=Path)
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--producer-job", required=True)
    p.add_argument("--cpus-per-chromosome", type=int, default=10)
    p.add_argument("--max-parallel", type=int, default=12)
    p.set_defaults(func=inventory)
    p = sub.add_parser("preserve-prefix", help="Validate/copy every completed ordered occurrence without publication")
    p.add_argument("--inventory", required=True, type=Path)
    p.add_argument("--inventory-sha256", required=True)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--cache-root", required=True, type=Path)
    p.add_argument("--vep-source-root", type=Path,
                   help="Pinned VEP116.3 source root; required to authorize its explicit nonvariant omissions")
    p.set_defaults(func=preserve_prefix)
    args = parser.parse_args()
    if args.operation == "inventory":
        require(args.cpus_per_chromosome > 1 and args.max_parallel > 0, "Positive resource limits required")
    args.func(args)


if __name__ == "__main__":
    main()
