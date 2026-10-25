#!/usr/bin/env python3
"""Inspect the exact indexed gnomAD chromosome family consumed by PriVA.

No score, allele, coordinate, filename, row or source byte is rewritten. The
existing caller substitutes a literal .chrX.vcf token to select chr1–22/X/Y.
These read-only checks establish its required headers and indexed access, not
scientific qualification of a different release or changed source dataset.

  python gnomad_cache.py members --chr-x /cache/selected.chrX.vcf.bgz
  python gnomad_cache.py validate --chr-x /cache/selected.chrX.vcf.bgz --assembly GRCh38 --reference /cache/reference.fa
  python gnomad_cache.py default-path --cache-dir /cache/gnomAD --assembly hg19
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys


CHROMOSOMES = tuple("chr" + str(number) for number in range(1, 23)) + ("chrX", "chrY")
CHR1_LENGTH = {"hg19": 249250621, "hg38": 248956422}
POPULATIONS = ("nfe", "eas", "afr", "amr", "asj", "fin", "sas", "mid", "remaining")
REQUIRED_INFO = (
    "AC_joint", "AN_joint", "AF_joint", "nhomalt_joint",
    "AC_grpmax_joint", "AF_grpmax_joint", "AN_grpmax_joint", "nhomalt_grpmax_joint",
    "AC_joint_XX", "AN_joint_XX", "AC_joint_XY", "AN_joint_XY", "AF_joint_XX", "AF_joint_XY",
    "nhomalt_joint_XX", "nhomalt_joint_XY",
    *(field for population in POPULATIONS for field in
      (f"AF_joint_{population}_XX", f"AF_joint_{population}_XY", f"nhomalt_joint_{population}_XX")),
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def assembly_name(value):
    require(value in ("hg19", "GRCh37", "hg38", "GRCh38"), "Select hg19/GRCh37 or hg38/GRCh38")
    return "hg19" if value in ("hg19", "GRCh37") else "hg38"


def default_path(*, cache_dir, assembly):
    assembly = assembly_name(assembly)
    name = ("gnomad.joint.v4.1.sites.hg19.chrX.vcf.gz" if assembly == "hg19"
            else "gnomad.joint.v4.1.sites.chrX.vcf.bgz")
    return Path(cache_dir).absolute() / name


def family_paths(chr_x):
    """Enumerate only the 24 names already required by the caller's template."""
    path = Path(chr_x).absolute()
    matched = re.fullmatch(r"(.+)\.chrX\.vcf(\.(?:gz|bgz))", path.name)
    require(matched is not None,
            "Selected gnomAD filename must retain the caller's .chrX.vcf.gz or .chrX.vcf.bgz template; refusing to substitute another dataset")
    require(not any(".partial" in part.lower() or part.startswith(".building") for part in path.parts),
            "Selected gnomAD family must not be a partial/building cache")
    prefix, suffix = matched.groups()
    return {chromosome: path.with_name(prefix + "." + chromosome + ".vcf" + suffix) for chromosome in CHROMOSOMES}


def identity(path):
    value = path.stat()
    return (value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)


def reference_lengths(reference):
    """Read an available selected reference index; never create/repair it here."""
    if reference is None:
        return {}
    index = Path(str(reference) + ".fai")
    if not index.is_file():
        return {}
    lengths = {}
    with index.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            require(len(fields) >= 2 and fields[1].isdigit(), "Malformed selected reference index: " + str(index))
            name = "chr" + fields[0].removeprefix("chr")
            if name in CHROMOSOMES:
                require(name not in lengths or lengths[name] == int(fields[1]), "Conflicting selected reference contig lengths")
                lengths[name] = int(fields[1])
    return lengths


def validate_family(*, chr_x, assembly, reference=None):
    """Check all exact members without downloading, indexing or changing bytes.

    Use mature htslib parsing/indexed fetches. Header contig lengths must match
    the selected reference where available and one coherent family header;
    chromosome 1 identifies the requested assembly. These are bounded indexed
    probes, not a full stream integrity scan or a new release qualification.
    """
    import pysam

    assembly = assembly_name(assembly)
    family = family_paths(chr_x)
    tracked, records = [], []
    for chromosome, path in family.items():
        require(path.is_file() and path.stat().st_size > 0, f"Selected gnomAD family member is missing/empty: {path}; restore that selected family")
        tracked.append((path, identity(path)))
        indexes = [Path(str(path) + suffix) for suffix in (".tbi", ".csi") if Path(str(path) + suffix).is_file()]
        require(indexes and all(index.stat().st_size > 0 for index in indexes), f"Selected gnomAD family requires its nonempty matching index: {path}")
        tracked.extend((index, identity(index)) for index in indexes)

    selected_lengths = reference_lengths(reference)
    require(not selected_lengths or selected_lengths.get("chr1") == CHR1_LENGTH[assembly],
            "Selected reference index does not identify the requested gnomAD assembly")
    try:
        with pysam.VariantFile(str(family["chr1"])) as anchor:
            family_lengths = {name: anchor.header.contigs[name].length for name in CHROMOSOMES
                              if name in anchor.header.contigs and anchor.header.contigs[name].length}
        require(family_lengths.get("chr1") == CHR1_LENGTH[assembly],
                "Selected gnomAD chr1 header does not identify the requested assembly")
        expected_lengths = dict(family_lengths, **selected_lengths)
        for chromosome, path in family.items():
            with pysam.VariantFile(str(path)) as variants:
                missing = [name for name in REQUIRED_INFO if name not in variants.header.info]
                require(not missing, f"Selected gnomAD header lacks caller-required INFO fields in {path}: {','.join(missing)}")
                for name in REQUIRED_INFO:
                    expected_type = "Float" if name.startswith("AF_") else "Integer"
                    require(variants.header.info[name].type == expected_type,
                            f"Selected gnomAD INFO/{name} has incompatible type in {path}")
                require(chromosome in variants.header.contigs and chromosome in expected_lengths,
                        f"Selected gnomAD header/reference lacks {chromosome} assembly length: {path}")
                require(variants.header.contigs[chromosome].length == expected_lengths[chromosome],
                        f"Selected gnomAD {chromosome} header differs from the selected assembly/reference: {path}")
                # htslib selects the same adjacent index as bcftools. A first
                # indexed record also catches unusable/truncated index formats.
                first = next(variants.fetch(chromosome), None)
                require(first is None or (first.contig == chromosome and 0 < first.pos <= expected_lengths[chromosome]),
                        f"Selected gnomAD indexed record is incompatible with its chromosome: {path}")
                records.append({"chromosome": chromosome, "path": str(path), "bytes": path.stat().st_size})
    except (OSError, ValueError, KeyError) as error:
        raise ValueError("Cannot reuse selected gnomAD family: " + str(error)) from error
    require(all(identity(path) == before for path, before in tracked), "Selected gnomAD family changed during inspection")
    return {"assembly": assembly, "chr_x": str(family["chrX"]), "members": records,
            "scientific_acceptance": False, "scope": "Original caller headers, coherent assembly and indexed access; no data/release substitution."}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    members = commands.add_parser("members", help="List the exact caller-required chromosome paths")
    members.add_argument("--chr-x", type=Path, required=True)
    validate = commands.add_parser("validate", help="Read-only header/assembly/index validation of the selected complete family")
    validate.add_argument("--chr-x", type=Path, required=True)
    validate.add_argument("--assembly", required=True)
    validate.add_argument("--reference", type=Path)
    default = commands.add_parser("default-path", help="Print the original managed release-4.1 chrX destination")
    default.add_argument("--cache-dir", type=Path, required=True)
    default.add_argument("--assembly", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        if args.command == "members":
            print(json.dumps({key: str(value) for key, value in family_paths(args.chr_x).items()}))
        elif args.command == "default-path":
            print(default_path(cache_dir=args.cache_dir, assembly=args.assembly))
        else:
            print(json.dumps(validate_family(chr_x=args.chr_x, assembly=args.assembly, reference=args.reference)))
    except (OSError, ValueError, TypeError) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
