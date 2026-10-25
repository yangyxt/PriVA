#!/usr/bin/env python3
"""Validate a Conservation BigWig and optionally publish a verified download.

Uses the already-declared perl-bio-bigfile reader, not an optional Python
scientific package. Resource checks do not establish scientific qualification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


ASSEMBLIES = {"hg19": "GRCh37", "GRCh37": "GRCh37", "hg38": "GRCh38", "GRCh38": "GRCh38"}
CHR1_LENGTH = {"GRCh37": 249250621, "GRCh38": 248956422}
GRCH38_DEFAULT_URL = (
    "https://ftp.ensembl.org/pub/release-114/compara/conservation_scores/"
    "91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
)
# Exact publisher size and BSD CHECKSUMS recorded by acquisition 1515153's
# frozen plan; BSD is not a cryptographic publisher checksum.
GRCH38_DEFAULT_PIN = {"bytes": 9611134022, "bsd_sum": [30633, 9385874]}
PERL_READER = r'''
use strict;
use warnings;
use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;
use JSON::PP;
my $path = shift @ARGV;
my $wig = Bio::DB::BigFile->bigWigFileOpen($path);
die "Cannot open BigWig\n" unless $wig;
my $list = $wig->chromList();
die "BigWig has no chromosome dictionary\n" unless $list;
my %chroms;
for (my $entry = $list->head; $entry; $entry = $entry->next) {
    my ($name, $size) = ($entry->name, $entry->size);
    die "Invalid BigWig chromosome\n" unless defined($name) && length($name) && $size > 0;
    die "Duplicate BigWig chromosome\n" if exists $chroms{$name};
    $chroms{$name} = 0 + $size;
    # Exercise actual signal/index reads without materializing whole tracks.
    my $window = $size < 4096 ? $size : 4096;
    $wig->bigWigSingleSummary($name, 0, $window, bbiSumMean);
    $wig->bigWigSingleSummary($name, $size - $window, $size, bbiSumMean) if $size > $window;
}
die "Empty BigWig chromosome dictionary\n" unless keys %chroms;
print JSON::PP->new->canonical->encode({chroms => \%chroms,
    reader => 'Bio::DB::BigFile', reader_module => $INC{'Bio/DB/BigFile.pm'}});
'''


def require(condition, message):
    if not condition:
        raise ValueError(message)


def identity(path):
    value = Path(path).stat()
    return {"device": value.st_dev, "inode": value.st_ino, "bytes": value.st_size,
            "mtime_ns": value.st_mtime_ns, "ctime_ns": value.st_ctime_ns}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized_chromosome(name):
    value = name[3:] if name.startswith("chr") else name
    return "MT" if value == "M" else value


def normalized_lengths(chromosomes, label):
    result = {}
    for name, length in chromosomes.items():
        key = normalized_chromosome(name)
        require(key not in result or result[key] == length,
                f"Conflicting normalized chromosome lengths in {label}: {name}")
        result[key] = length
    return result


def reference_lengths(reference, assembly):
    reference = Path(reference)
    require(reference.is_file() and reference.stat().st_size > 0,
            f"Missing or empty configured reference: {reference}")
    index = Path(str(reference) + ".fai")
    require(index.is_file() and index.stat().st_size > 0,
            f"Missing configured reference FASTA index: {index}; index the selected reference first")
    before = identity(index)
    chromosomes = {}
    with index.open(encoding="utf-8") as stream:
        for line in stream:
            columns = line.rstrip("\r\n").split("\t")
            require(len(columns) >= 5 and columns[0] and all(value.isdigit() for value in columns[1:5])
                    and int(columns[1]) > 0 and int(columns[3]) > 0 and int(columns[4]) >= int(columns[3]),
                    f"Malformed configured reference FASTA index: {index}")
            require(columns[0] not in chromosomes, f"Duplicate reference chromosome: {columns[0]}")
            chromosomes[columns[0]] = int(columns[1])
    normalized = normalized_lengths(chromosomes, "reference")
    require(normalized.get("1") == CHR1_LENGTH[assembly],
            f"Configured reference chr1 length does not identify {assembly}: {reference}")
    record = {"path": str(reference.absolute()), "identity": identity(reference),
              "fai": str(index.absolute()), "fai_sha256": sha256(index)}
    require(identity(index) == before, "Configured reference index changed during validation")
    return normalized, record


def bigwig_chromosomes(path):
    # Prefer the Perl that belongs to the running Python's installed environment.
    perl = Path(sys.executable).parent / "perl"
    executable = str(perl) if perl.is_file() else shutil.which("perl")
    require(executable, "Install the declared perl-bio-bigfile runtime before Conservation validation")
    result = subprocess.run([executable, "-e", PERL_READER, str(path)],
                            text=True, capture_output=True, timeout=120)
    require(result.returncode == 0,
            f"Cannot validate Conservation BigWig with declared Bio::DB::BigFile: {path}: "
            + result.stderr.strip())
    value = json.loads(result.stdout)
    require(isinstance(value.get("chroms"), dict) and value["chroms"], "Invalid BigWig chromosome dictionary")
    return value


def bsd_sum(path):
    result = subprocess.run(["sum", "-r", "--", str(path)], check=True,
                            capture_output=True, text=True, env=dict(os.environ, LC_ALL="C"))
    columns = result.stdout.split()
    require(len(columns) >= 2 and all(value.isdigit() for value in columns[:2]), "Invalid BSD checksum output")
    return list(map(int, columns[:2]))


def validate(path, assembly, reference, source_url=None):
    assembly = ASSEMBLIES.get(assembly)
    require(assembly, "Conservation requires hg19/GRCh37 or hg38/GRCh38")
    path = Path(path)
    require(path.is_file() and path.stat().st_size > 0, f"Missing or empty Conservation BigWig: {path}")
    before = identity(path)
    known_default = source_url == GRCH38_DEFAULT_URL
    if known_default:
        require(assembly == "GRCh38", "Archived release114 Conservation requires GRCh38")
        require(before["bytes"] == GRCH38_DEFAULT_PIN["bytes"],
                "Archived release114 Conservation publisher byte count mismatch; resource retained")
    lengths, reference_record = reference_lengths(reference, assembly)
    bigwig = bigwig_chromosomes(path)
    normalized = normalized_lengths(bigwig["chroms"], "BigWig")
    require(normalized.get("1") == CHR1_LENGTH[assembly],
            f"Conservation BigWig chr1 length does not identify {assembly}: {path}")
    shared = sorted(normalized.keys() & lengths.keys())
    for name in shared:
        require(normalized[name] == lengths[name],
                f"Conservation/reference chromosome length mismatch for {name}: "
                f"{normalized[name]} != {lengths[name]}")
    measured_bsd = bsd_sum(path) if known_default else None
    if known_default:
        require(measured_bsd == GRCH38_DEFAULT_PIN["bsd_sum"],
                "Archived release114 Conservation publisher BSD checksum mismatch; resource retained")
    digest = sha256(path)
    require(identity(path) == before, "Conservation resource changed during validation; resource retained")
    return {"path": str(path.absolute()), "assembly": assembly, "sha256": digest, "identity": before,
            "reader": bigwig["reader"], "reader_module": bigwig["reader_module"],
            "chromosome_count": len(bigwig["chroms"]), "shared_reference_chromosomes": shared,
            "extra_non_reference_chromosomes": sorted(normalized.keys() - lengths.keys()),
            "reference": reference_record, "source_url": source_url,
            "source_kind": "download" if source_url else "configured_existing_identity",
            "release": 114 if known_default else None,
            "analysis": "91_mammals.gerp_conservation_score" if known_default else None,
            "publisher_byte_count": GRCH38_DEFAULT_PIN["bytes"] if known_default else None,
            "publisher_bsd_sum": measured_bsd, "publisher_cryptographic_checksum_claimed": False,
            "scientific_qualified": False}


def publish(path, destination, assembly, reference, source_url):
    path, destination = Path(path), Path(destination)
    require(source_url, "Download publication requires its explicit source URL")
    require(path.absolute() != destination.absolute(), "Download staging and final paths must differ")
    require(not path.is_symlink(), "Conservation staging file must not be a symlink")
    record = validate(path, assembly, reference, source_url)
    require(identity(path) == record["identity"], "Conservation staging file changed before publication")
    # Atomic link refuses to overwrite any existing target. Errors retain staging.
    os.link(path, destination)
    path.unlink()
    record.update(path=str(destination.absolute()), identity=identity(destination))
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--source-url")
    parser.add_argument("--publish-to", type=Path)
    args = parser.parse_args(argv)
    try:
        result = (publish(args.input, args.publish_to, args.assembly, args.reference, args.source_url)
                  if args.publish_to else validate(args.input, args.assembly, args.reference, args.source_url))
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        print(f"ERROR: Conservation validation failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
