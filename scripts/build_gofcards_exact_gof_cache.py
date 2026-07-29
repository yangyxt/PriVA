#!/usr/bin/env python3
"""Build the PriVA GoFCards exact gain-of-function evidence cache.

GoFCards publishes GRCh37 coordinates only, so GRCh37 is the sole ground truth
here.  Everything else is derived from it:

  1. fetch_sources      public workbook (one download)
     + fetch_summary_annotations   GRCh38 position and ClinVar accession,
       one cached request per allele
  2. build_source_vcf   pad sparse indels, check REF against GRCh37, gate failures
  3. normalize_vcf      bcftools norm                  (reused for both assemblies)
  4. run_vep            Ensembl VEP                    (reused for both assemblies)
  5. liftover_vcf       CrossMap GRCh37 -> GRCh38, then normalize_vcf again
  6. run_vep on the lifted VCF
  7. build_cache        join, resolve HGNC symbols, apply reviews, export

Two checks carry the weight, and they are deliberately separate.  The reference
check in step 2 decides which records are *trustworthy*: an allele whose REF
does not match GRCh37 describes a variant that does not exist, so its HGVS is
meaningless no matter how confidently VEP reports it.  The liftover in step 5
decides only which records can *offer a GRCh38 coordinate key*; failing it
costs an allele one match route, not its eligibility.

HGVS is a property of the transcript sequence, not of the assembly.  The two
VEP runs exist because the GRCh37 and GRCh38 caches ship different transcript
catalogues, and PriVA must be able to match against whichever one it annotated
its own query variants with.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import requests
from pyfaidx import Fasta

sys.path.insert(0, str(Path(__file__).resolve().parent))

from clinvar_vcv import (  # noqa: E402
    clean_text,
    normalize_allele,
    normalize_chrom,
    normalize_int,
    variant_id_of,
)

BACKEND_BASE = "https://java.genemed.tech/admin-api/backend/data/hg19"
PUBLIC_EXCEL_URL = "https://download.genemed.tech/upload/GainFunCards/gofcards_data_download.xlsx"
# One request per allele. It is the only place GoFCards publishes its own GRCh38
# position and its ClinVar cross-reference, and it returns neither a GRCh38
# reference nor alternate allele -- only a position, which is why every position
# it gives has to be checked against the genome before it can be used.
SUMMARY_ENDPOINT = f"{BACKEND_BASE}/variantLevel/summary"
# Only two fields are kept. The position recovers variants the liftover
# cannot place; the accession is GoFCards' ClinVar VariationID, which links
# a variant to its ClinVar record by identity. The endpoint also returns
# CLNSIG, CLNDN, CLNREVSTAT and rsID, and nothing downstream reads any of
# them -- ClinVar significance is taken live from ClinVar itself.
SUMMARY_FIELDS = ("hg38_start", "Accession")
API_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
        "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125 Safari/537.36"
    ),
    "Accept": "application/json, text/plain, */*",
    "Referer": "https://www.genemed.tech/gofcards/",
    "Origin": "https://www.genemed.tech",
    "Admin-Work-Version": "x",
    "Content-Type": "application/json; charset=UTF-8",
}

VEP_FIELDS = (
    "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,"
    "HGVSc,HGVSp,MANE_SELECT,MANE_PLUS_CLINICAL,CANONICAL,SYMBOL"
)

ELIGIBLE = "eligible"
QUARANTINE_GENE = "quarantined_gene_discordance"
QUARANTINE_SIBLING = "quarantined_allele_gene_discordance"

# The canonical schema encodes the reviewed mechanism in the quarantine state
# itself, so a reader can tell LOF from mixed without opening the review table.
REVIEWED_QUARANTINE = {
    "LOF": "quarantined_reviewed_lof",
    "MIXED": "quarantined_reviewed_mixed",
    "DOMINANT_NEGATIVE": "quarantined_reviewed_dominant_negative",
    "DN": "quarantined_reviewed_dominant_negative",
    "UNCERTAIN": "quarantined_reviewed_uncertain",
    "EXCLUDE": "quarantined_reviewed_exclude",
}
QUARANTINE_REVIEW_REQUIRED = "quarantined_mechanism_review_required"
# Reference-check failures never reach the cache at all; they are written to the
# rejects table in step 2, before any annotation happens.

REVIEW_COLUMNS = [
    "source_order", "source_gene", "allele_key", "source_mechanism",
    "reviewed_mechanism", "review_status", "gof_eligibility", "reason_code",
    "evidence_summary", "pmids", "reviewer", "reviewed_at", "review_version",
]


def log(message: str) -> None:
    print(f"[gofcards] {message}", file=sys.stderr, flush=True)


def run_command(cmd: list[str], log_path: Path) -> None:
    """Run a command, merging stdout and stderr into one log file."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log(f"$ {' '.join(cmd)}")
    with log_path.open("w", encoding="utf-8") as handle:
        result = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}); see {log_path}: {' '.join(cmd)}")


def require_executable(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise RuntimeError(
            f"Required executable {name!r} is not on PATH. Install it into the PriVA "
            f"environment rather than working around it."
        )
    return path


# ---------------------------------------------------------------------------
# Shared normalization of the raw GoFCards fields
# ---------------------------------------------------------------------------

# clean_text, normalize_chrom, normalize_allele, normalize_int and
# variant_id_of are imported from clinvar_vcv, which owns the cache contract.
# The identifier this module mints is the same one every consumer looks a
# variant up by, so it is defined once, in the module both sides share.


def locus_key_of(row: Any) -> str:
    """Type-free allele identity: chrom|start|end|ref|alt, all GRCh37.

    This is the real identity.  The `variant_type` label is derived, not
    published consistently: GoFCards' backend calls multi-base substitutions
    SNV while the public workbook calls them Indel, so the label must never be
    part of a key used to join anything.
    """
    return "|".join([
        normalize_chrom(row.get("Chr")),
        normalize_int(row.get("Start")),
        normalize_int(row.get("End")),
        normalize_allele(row.get("Ref")),
        normalize_allele(row.get("Alt")),
    ])


def allele_key_of(row: Any) -> str:
    """Published allele key, `type|chrom|start|end|ref|alt`.

    Retained in this exact shape because PriVA deduplicates GoFCards matches on
    it. Use `locus_key_of` for joining, never this.
    """
    return f"{clean_text(row.get('variant_type'))}|{locus_key_of(row)}"


# ---------------------------------------------------------------------------
# Step 1 - fetch
# ---------------------------------------------------------------------------

PUBLIC_COLUMN_MAP = {
    "genesymbol": "Gene_Symbol", "transcript": "Transcript", "chr": "Chr",
    "hg19start": "Start", "hg19end": "End", "ref": "Ref", "alt": "Alt",
    "function": "Function", "pathways proteins involved": "Pathways_proteins_involved",
    "disorder involved": "Phenotype", "pmid": "PMID",
    "animal model": "Animal_model", "cell model": "Cell_model",
    "pscore": "Pscore", "order numbe": "source_order",
}


def fetch_summary_annotations(records: pd.DataFrame, cache_jsonl: Path,
                              timeout_seconds: int = 300) -> dict[str, dict]:
    """Query the per-allele summary endpoint once for every unique allele.

    This is 2,033 sequential HTTP requests and the endpoint has been observed
    both slow and entirely unreachable, so every response is appended to a
    JSONL cache keyed by allele. A rerun re-fetches only what is missing, and an
    outage part-way through loses nothing already collected.

    It supplies two things available nowhere else in GoFCards: the database's own
    GRCh38 position, and its ClinVar cross-reference (``Accession`` plus the
    clinical significance, condition names and review status).
    """
    cache: dict[str, dict] = {}
    if cache_jsonl.exists():
        with cache_jsonl.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if line:
                    row = json.loads(line)
                    cache[row["allele_key"]] = row.get("data") or {}
        log(f"summary cache already holds {len(cache)} alleles")

    unique = records.drop_duplicates("locus_key")
    session = requests.Session()
    cache_jsonl.parent.mkdir(parents=True, exist_ok=True)
    pending = [r for _, r in unique.iterrows() if r["allele_key"] not in cache]
    if pending:
        log(f"querying the summary endpoint for {len(pending)} uncached alleles")
    with cache_jsonl.open("a", encoding="utf-8") as handle:
        for index, row in enumerate(pending, start=1):
            params = {
                "projectCode": "GoFCards",
                "variantLevelType": "Indel" if row["variant_type"] == "Indel" else "SNV",
                "chr": normalize_chrom(row.get("Chr")),
                "start": normalize_int(row.get("Start")),
                "end": normalize_int(row.get("End")),
                "ref": normalize_allele(row.get("Ref")),
                "alt": normalize_allele(row.get("Alt")),
            }
            data: dict = {}
            for attempt in range(1, 5):
                try:
                    reply = session.get(SUMMARY_ENDPOINT, headers=API_HEADERS,
                                        params=params, timeout=timeout_seconds)
                    reply.raise_for_status()
                    data = (reply.json() or {}).get("data") or {}
                    break
                except Exception as exc:
                    if attempt == 4:
                        raise RuntimeError(
                            f"summary endpoint failed for {row['allele_key']} "
                            f"after 4 attempts: {exc}") from exc
                    time.sleep(5 * attempt)
            kept = {field: data[field] for field in SUMMARY_FIELDS if data.get(field)}
            cache[row["allele_key"]] = kept
            handle.write(json.dumps({"allele_key": row["allele_key"], "data": kept},
                                    ensure_ascii=False, sort_keys=True) + "\n")
            handle.flush()
            if index == 1 or index % 200 == 0:
                log(f"  summary {index}/{len(pending)}")
    with_hg38 = sum(1 for v in cache.values() if v.get("hg38_start"))
    with_clinvar = sum(1 for v in cache.values() if v.get("Accession"))
    log(f"summary endpoint: {len(cache)} alleles; {with_hg38} carry a GRCh38 position; "
        f"{with_clinvar} carry a ClinVar accession")
    return cache


def fetch_sources(workdir: Path, public_url: str,
                  timeout_seconds: int = 300) -> tuple[pd.DataFrame, dict]:
    """Download the public GoFCards workbook.

    One request. The workbook is the source of record: it is the citable
    artifact, it is hash-stable, and it carries every coordinate and every
    evidence field the cache needs.

    The backend SNV/Indel table endpoints are deliberately not called. They add
    only ANNOVAR's protein change, which nothing downstream reads -- our own
    protein change comes from VEP, on transcripts whose versions match the ones
    PriVA annotates with.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    session = requests.Session()

    public_path = workdir / "gofcards_public.xlsx"
    log(f"Downloading public workbook -> {public_path}")
    resp = session.get(public_url, headers={"User-Agent": "Mozilla/5.0"},
                       timeout=timeout_seconds)
    resp.raise_for_status()
    public_path.write_bytes(resp.content)

    raw = pd.read_excel(public_path, sheet_name=0)
    public = raw.rename(columns={c: PUBLIC_COLUMN_MAP[str(c).strip().lower()]
                                 for c in raw.columns
                                 if str(c).strip().lower() in PUBLIC_COLUMN_MAP})
    public["variant_type"] = public.apply(
        lambda r: "SNV" if len(normalize_allele(r.get("Ref"))) == 1
        and len(normalize_allele(r.get("Alt"))) == 1 else "Indel", axis=1)
    public["allele_key"] = public.apply(allele_key_of, axis=1)
    public["locus_key"] = public.apply(locus_key_of, axis=1)

    provenance = {
        "public_excel_url": public_url,
        "public_excel_last_modified": resp.headers.get("Last-Modified", ""),
        "public_workbook_rows": int(len(public)),
        "public_unique_alleles": int(public["locus_key"].nunique()),
    }
    log(f"public workbook: {len(public)} rows, "
        f"{public['locus_key'].nunique()} unique alleles")
    return public, provenance


# ---------------------------------------------------------------------------
# Step 2 - build the GRCh37 VCF and gate on the reference allele
# ---------------------------------------------------------------------------

def fasta_contig_name(fasta: Fasta, chrom: str) -> str | None:
    """Return the contig name this FASTA actually uses for a bare chromosome.

    References differ in whether contigs are `1` or `chr1`, and every downstream
    tool -- bcftools, CrossMap, VEP -- needs the VCF to agree with the FASTA it
    is given.  Resolving the name once here keeps that agreement in one place.
    """
    for candidate in (chrom, f"chr{chrom}", chrom.removeprefix("chr")):
        if candidate in fasta:
            return candidate
    return None


def _fasta_slice(fasta: Fasta, chrom: str, start_1based: int, length: int) -> str:
    name = fasta_contig_name(fasta, chrom)
    if name is None:
        raise KeyError(f"Chromosome {chrom!r} not present in FASTA")
    if start_1based < 1:
        raise ValueError(f"Invalid coordinate {start_1based}")
    return str(fasta[name][start_1based - 1:start_1based - 1 + length]).upper()


def _vcf_representation(fasta: Fasta, chrom: str, start: int, ref: str, alt: str) -> dict:
    """Convert one GoFCards allele into a valid VCF record, or explain why not.

    GoFCards writes indels in the sparse ANNOVAR style with an empty REF or ALT,
    so those need an anchor base from the genome.  Substitutions are checked
    directly.  A REF that disagrees with GRCh37 is reported, never silently
    replaced: replacing it would invent a variant nobody observed.
    """
    try:
        if ref and alt:
            observed = _fasta_slice(fasta, chrom, start, max(len(ref), 1))
            if observed != ref:
                return {"status": "ref_mismatch", "detail": f"GRCh37 has {observed}, source claims {ref}"}
            return {"status": "ok", "kind": "raw_ref_alt", "pos": start, "ref": ref, "alt": alt}

        if ref and not alt:  # deletion, needs a preceding anchor base
            anchor_pos = start - 1
            anchor = _fasta_slice(fasta, chrom, anchor_pos, 1)
            observed = _fasta_slice(fasta, chrom, start, max(len(ref), 1))
            if observed != ref:
                return {"status": "ref_mismatch", "detail": f"GRCh37 has {observed}, source claims {ref}"}
            return {"status": "ok", "kind": "deletion_padded_ref_match",
                    "pos": anchor_pos, "ref": anchor + ref, "alt": anchor}

        if alt and not ref:  # insertion, anchored on the base at start
            anchor = _fasta_slice(fasta, chrom, start, 1)
            return {"status": "ok", "kind": "insertion_padded",
                    "pos": start, "ref": anchor, "alt": anchor + alt}

        return {"status": "missing_ref_and_alt", "detail": "source has neither REF nor ALT"}
    except Exception as exc:
        return {"status": "fasta_lookup_error", "detail": str(exc)}


def build_source_vcf(source_records: pd.DataFrame, hg19_fasta: Path, out_vcf: Path,
                     rejects_tsv: Path) -> pd.DataFrame:
    """Write one GRCh37 VCF record per unique allele; reject what fails the check.

    Every record carries a short synthetic ID that survives bcftools, CrossMap
    and VEP unchanged, so downstream tables join back by identity rather than
    by position.
    """
    fasta = Fasta(str(hg19_fasta), sequence_always_upper=True, rebuild=False)
    unique = source_records.drop_duplicates("allele_key").reset_index(drop=True)

    kept: list[dict] = []
    rejected: list[dict] = []
    for index, row in unique.iterrows():
        chrom = normalize_chrom(row.get("Chr"))
        start_text = normalize_int(row.get("Start"))
        if not chrom or not start_text.isdigit():
            rejected.append({"allele_key": row["allele_key"], "reject_reason": "missing_coordinate",
                             "detail": f"chrom={chrom!r} start={start_text!r}"})
            continue
        rep = _vcf_representation(fasta, chrom, int(start_text),
                                  normalize_allele(row.get("Ref")), normalize_allele(row.get("Alt")))
        if rep["status"] != "ok":
            rejected.append({"allele_key": row["allele_key"],
                             "reject_reason": rep["status"], "detail": rep.get("detail", "")})
            continue
        kept.append({
            "vcf_id": f"gof{index:06d}",
            "allele_key": row["allele_key"],
            "hg19_chrom": chrom,
            "vcf_contig": fasta_contig_name(fasta, chrom),
            "hg19_pos": rep["pos"],
            "hg19_ref": rep["ref"],
            "hg19_alt": rep["alt"],
            "hg19_vcf_status": rep["kind"],
            "source_hg19_pos": start_text,
            "source_hg19_ref": normalize_allele(row.get("Ref")),
            "source_hg19_alt": normalize_allele(row.get("Alt")),
        })

    core = pd.DataFrame(kept)
    # Write contigs under the names this FASTA uses, with a header entry each, so
    # bcftools and CrossMap can resolve every record against the reference.
    contigs = {rec["vcf_contig"] for rec in kept}
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    with out_vcf.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for name in sorted(contigs):
            handle.write(f"##contig=<ID={name},length={len(fasta[name])}>\n")
        handle.write('##INFO=<ID=GOFID,Number=1,Type=String,Description="GoFCards record id">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for rec in sorted(kept, key=lambda r: (r["vcf_contig"], int(r["hg19_pos"]))):
            handle.write(f"{rec['vcf_contig']}\t{rec['hg19_pos']}\t{rec['vcf_id']}\t"
                         f"{rec['hg19_ref']}\t{rec['hg19_alt']}\t.\t.\tGOFID={rec['vcf_id']}\n")

    pd.DataFrame(rejected).to_csv(rejects_tsv, sep="\t", index=False)
    log(f"reference check: {len(kept)} alleles kept, {len(rejected)} rejected -> {rejects_tsv}")
    return core


# ---------------------------------------------------------------------------
# Step 3 and 5b - one normalizer, used for both assemblies
# ---------------------------------------------------------------------------

def normalize_vcf(in_vcf: Path, fasta: Path, out_vcf: Path, log_dir: Path) -> Path:
    """Left-align and split multi-allelic records against the given genome.

    Needed once per assembly: an indel that is left-aligned in GRCh37 is not
    guaranteed to be left-aligned in GRCh38, because the surrounding repeat
    structure can differ between the two.
    """
    require_executable("bcftools")
    run_command(
        ["bcftools", "norm", "-f", str(fasta), "-c", "e", "-m-any", "-Ov", "-o", str(out_vcf), str(in_vcf)],
        log_dir / f"bcftools_norm.{out_vcf.stem}.log",
    )
    return out_vcf


# ---------------------------------------------------------------------------
# Step 4 and 6 - one VEP runner, used for both assemblies
# ---------------------------------------------------------------------------

def run_vep(in_vcf: Path, assembly: str, cache_dir: Path, cache_version: str,
            fasta: Path, out_tsv: Path, log_dir: Path, merged: bool = True) -> Path:
    """Annotate one VCF with Ensembl VEP, keeping every transcript.

    All transcripts are retained on purpose.  Transcripts disagree with each
    other about the protein change for most variants, so keeping only one would
    make a match depend on PriVA having picked the same transcript we did.
    """
    vep = require_executable("vep")
    cmd = [
        vep, "--input_file", str(in_vcf), "--output_file", str(out_tsv),
        "--format", "vcf", "--tab", "--force_overwrite",
        "--assembly", assembly, "--fasta", str(fasta),
        "--symbol", "--hgvs", "--mane", "--canonical", "--transcript_version",
        "--fields", VEP_FIELDS,
        "--cache", "--offline", "--dir_cache", str(cache_dir),
    ]
    if cache_version:
        cmd += ["--cache_version", str(cache_version)]
    if merged:
        cmd.append("--merged")
    run_command(cmd, log_dir / f"vep.{assembly}.log")
    return out_tsv


def parse_vep(tsv: Path, assembly: str) -> pd.DataFrame:
    """Read a VEP tab output into one row per (record, transcript)."""
    rows: list[dict] = []
    with tsv.open("r", encoding="utf-8") as handle:
        lines = [line for line in handle if not line.startswith("##")]
    for rec in csv.DictReader(lines, delimiter="\t"):
        feature = clean_text(rec.get("Feature"))
        if not feature:
            continue
        base, _, version = feature.partition(".")
        hgvsc = clean_text(rec.get("HGVSc"))
        hgvsp = clean_text(rec.get("HGVSp"))
        rows.append({
            "vcf_id": clean_text(rec.get("#Uploaded_variation")),
            "VEP_assembly": assembly,
            "VEP_transcript": feature,
            "transcript_base": base,
            "transcript_version": version,
            "feature_type": clean_text(rec.get("Feature_type")),
            "consequence": clean_text(rec.get("Consequence")),
            "canonical_transcript": clean_text(rec.get("CANONICAL")),
            "mane_select": clean_text(rec.get("MANE_SELECT")),
            "VEP_HGNC_Symbol": clean_text(rec.get("SYMBOL")),
            "HGVSc": hgvsc,
            "HGVSp": hgvsp,
            "hgvsc_key": hgvsc.split(":")[-1].upper() if hgvsc else "",
            "hgvsp_key": re.sub(r"^P\.", "", hgvsp.split(":")[-1].upper()) if hgvsp else "",
        })
    frame = pd.DataFrame(rows)
    log(f"VEP {assembly}: {len(frame)} transcript rows over {frame['vcf_id'].nunique() if not frame.empty else 0} records")
    return frame


# ---------------------------------------------------------------------------
# Step 5 - liftover
# ---------------------------------------------------------------------------

def liftover_vcf(in_vcf: Path, chain: Path, target_fasta: Path, out_vcf: Path,
                 log_dir: Path) -> tuple[Path, set[str]]:
    """Map GRCh37 records onto GRCh38 with CrossMap.

    Records CrossMap cannot place are written to a companion .unmap file and
    reported here.  They lose their GRCh38 coordinate key and nothing else:
    their HGVS came from the GRCh37 run and is unaffected.
    """
    require_executable("CrossMap")
    run_command(["CrossMap", "vcf", str(chain), str(in_vcf), str(target_fasta), str(out_vcf)],
                log_dir / "crossmap.log")
    unmapped: set[str] = set()
    unmap_path = Path(f"{out_vcf}.unmap")
    if unmap_path.exists():
        with unmap_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 2 and fields[2]:
                    unmapped.add(fields[2])
    log(f"liftover: {len(unmapped)} records could not be placed on GRCh38")
    return out_vcf, unmapped


def endpoint_hg38_coordinates(core: pd.DataFrame, summary: dict[str, dict],
                              hg38_fasta: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    """Turn GoFCards' own GRCh38 positions into checked coordinates.

    The endpoint returns a position and nothing else -- no GRCh38 reference or
    alternate allele. The alleles therefore have to be carried over from GRCh37,
    which is only valid if the GRCh38 sequence at that position really does
    start with the same reference. That check is the whole safeguard: it is what
    rejects a position pointing at the wrong part of the chromosome.

    Anything that fails here gets no endpoint coordinate and falls through to the
    liftover instead.
    """
    fasta = Fasta(str(hg38_fasta), sequence_always_upper=True, rebuild=False)
    rows: list[dict] = []
    counts: Counter[str] = Counter()
    for _, record in core.iterrows():
        data = summary.get(record["allele_key"]) or {}
        position = normalize_int(data.get("hg38_start"))
        if not position.isdigit():
            counts["no_endpoint_position"] += 1
            continue
        chrom, ref, alt = record["hg19_chrom"], record["hg19_ref"], record["hg19_alt"]
        try:
            observed = _fasta_slice(fasta, chrom, int(position), len(ref))
        except Exception:
            counts["endpoint_position_off_the_contig"] += 1
            continue
        if observed != ref:
            counts["endpoint_position_rejected_by_reference_check"] += 1
            continue
        counts["endpoint_position_accepted"] += 1
        rows.append({"vcf_id": record["vcf_id"], "hg38_chrom": chrom,
                     "hg38_pos": position, "hg38_ref": ref, "hg38_alt": alt})
    return pd.DataFrame(rows), dict(counts)


def read_vcf_coordinates(vcf: Path, prefix: str) -> pd.DataFrame:
    """Read back CHROM/POS/REF/ALT keyed by record id."""
    rows: list[dict] = []
    with vcf.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or not fields[2]:
                continue
            rows.append({
                "vcf_id": fields[2],
                f"{prefix}_chrom": normalize_chrom(fields[0]),
                f"{prefix}_pos": fields[1],
                f"{prefix}_ref": fields[3].upper(),
                f"{prefix}_alt": fields[4].upper(),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# HGNC symbol resolution
# ---------------------------------------------------------------------------

class HgncResolver:
    """Offline resolver mapping aliases and withdrawn symbols to approved ones."""

    def __init__(self, table: Path | None):
        self.approved: set[str] = set()
        self.aliases: dict[str, set[str]] = {}
        self.ids: dict[str, str] = {}
        if table and Path(table).exists():
            self._load(Path(table))
            log(f"HGNC resolver loaded {len(self.approved)} approved symbols from {table}")
        else:
            raise RuntimeError(
                f"HGNC complete-set table is required but was not found: {table}. "
                f"Install it rather than skipping symbol normalization."
            )

    def _load(self, table: Path) -> None:
        with table.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                symbol = clean_text(row.get("symbol")).upper()
                if not symbol:
                    continue
                self.approved.add(symbol)
                hgnc_id = clean_text(row.get("hgnc_id"))
                if hgnc_id:
                    self.ids[symbol] = hgnc_id
                for field in ("alias_symbol", "prev_symbol"):
                    for alias in re.split(r"[|,]", clean_text(row.get(field))):
                        alias = alias.strip().upper()
                        if alias:
                            self.aliases.setdefault(alias, set()).add(symbol)

    def resolve(self, value: Any) -> str:
        query = clean_text(value).upper()
        if not query or query in self.approved:
            return query
        matches = self.aliases.get(query, set())
        return next(iter(matches)) if len(matches) == 1 else query

    def hgnc_id(self, symbol: Any) -> str | None:
        """Return the stable HGNC identifier for an approved symbol, if known."""
        return self.ids.get(clean_text(symbol).upper())


# ---------------------------------------------------------------------------
# Mechanism polarity reviews
# ---------------------------------------------------------------------------

def load_mechanism_reviews(path: Path | None) -> dict[tuple[str, str], dict]:
    """Load validated per-allele polarity decisions, keyed by (gene, allele)."""
    if not path or not Path(path).exists():
        log("no mechanism review table supplied; every allele keeps its source GOF claim")
        return {}
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = [c for c in REVIEW_COLUMNS if c not in frame.columns]
    if missing:
        raise RuntimeError(f"{path} is missing review columns: {', '.join(missing)}")
    reviews = {
        (clean_text(r["source_gene"]).upper(), clean_text(r["allele_key"])): r.to_dict()
        for _, r in frame.iterrows()
    }
    log(f"loaded {len(reviews)} reviewed alleles from {path}")
    return reviews


# ---------------------------------------------------------------------------
# Step 7 - assemble and export
# ---------------------------------------------------------------------------

def _aggregate_evidence(source_records: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, list[dict]]]:
    """Group source rows by allele, returning join fields and per-publication evidence.

    An allele can carry many independent literature records -- KIT p.Asp816Val
    has 33 -- so each is kept as its own entry rather than being flattened into
    one string. The frame carries only the few fields the pipeline itself needs
    to join on; the bulky text stays in the per-publication records, where it is
    stored once instead of once per transcript.
    """
    def join_unique(values: Iterable[Any]) -> str:
        seen: list[str] = []
        for value in values:
            text = clean_text(value)
            if text and text not in seen:
                seen.append(text)
        return "; ".join(seen)

    join_fields: list[dict] = []
    records: dict[str, list[dict]] = {}
    for key, group in source_records.groupby("allele_key", sort=False):
        join_fields.append({
            "allele_key": key,
            "GoFCards_source_symbol": join_unique(group.get("Gene_Symbol", pd.Series(dtype=str))),
            "GoFCards_transcript": join_unique(group.get("Transcript", pd.Series(dtype=str))),
        })
        records[key] = [{
            "source_order": clean_text(row.get("source_order")),
            "pmid": clean_text(row.get("PMID")),
            "pscore": clean_text(row.get("Pscore")),
            "disease": clean_text(row.get("Phenotype")),
            "function": clean_text(row.get("Function")),
            "pathway": clean_text(row.get("Pathways_proteins_involved")),
            "animal_model": clean_text(row.get("Animal_model")),
            "cell_model": clean_text(row.get("Cell_model")),
            "source_transcript": clean_text(row.get("Transcript")),
        } for _, row in group.iterrows()]
    return pd.DataFrame(join_fields), records


def _nest(merged: pd.DataFrame, evidence_records: dict[str, list[dict]],
          summary: dict[str, dict], resolver: HgncResolver,
          provenance: dict, unmapped: set[str]) -> dict:
    """Fold the joined table into gene -> assembly -> {transcripts, genomic}.

    A GoFCards record is established by one of two routes, and the structure is
    keyed on exactly those:

      route 1   assembly + gene + transcript + HGVSc  (HGVSp maps to HGVSc)
      route 2   assembly + chromosome, position, reference and alternate allele

    Both routes are assembly-bound: a transcript identifier carries a different
    version per genome build, and a coordinate is meaningless without naming the
    build it belongs to.  Assembly is therefore the second key, above everything
    it qualifies.  Nothing derived from one build is allowed to sit above it.

    ``by_hgvsc`` is keyed because (gene, assembly, transcript, HGVSc) is unique
    across the catalogue.  ``by_hgvsp`` maps to a LIST of HGVSc keys because a
    protein change is not unique: two different coding changes can produce it.

    Facts that do not depend on the build -- the eligibility verdict, the
    literature, and what GoFCards originally published -- live once per source
    record under ``records``, which both assembly views reference by key.
    """
    genes: dict[str, dict] = {}
    for symbol, gene_rows in merged.groupby("HGNC_Symbol", sort=True):
        if not symbol:
            continue
        gene = genes.setdefault(symbol, {"hgnc_id": resolver.hgnc_id(symbol), "variants": {}})

        for source_key, rows in gene_rows.groupby("allele_key", sort=True):
            first = rows.iloc[0]
            identifier = variant_id_of(first["hg19_chrom"], first["source_hg19_pos"],
                                       first["source_hg19_ref"], first["source_hg19_alt"])
            variant = gene["variants"].setdefault(identifier, {
                # Everything that does not depend on the genome build, stored once.
                "record": {
                    "eligibility": {
                        "status": first["match_eligibility"],
                        "gene_match_status": first["gene_match_status"],
                        "reason": first["reject_reason"] or None,
                        "vep_symbol": first["VEP_HGNC_Symbol"] or None,
                    },
                    "source": {
                        "gofcards_allele_key": source_key,
                        "variant_type_label": source_key.split("|")[0],
                        "assembly": "hg19",
                        "chrom": first["hg19_chrom"],
                        "start": first["source_hg19_pos"],
                        "ref": first["source_hg19_ref"] or "-",
                        "alt": first["source_hg19_alt"] or "-",
                    },
                    "liftover_status": "unmapped" if first["vcf_id"] in unmapped else "mapped",
                    # Only what a downstream consumer reads. The ClinVar
                    # variation id links a variant to its ClinVar record by
                    # identity and reaches the final ACMG output. GoFCards' own
                    # snapshot of ClinVar significance, conditions and review
                    # status is deliberately not kept: the injection step
                    # attaches the live ClinVar record, and nothing reads the
                    # snapshot.
                    "annotations": {
                        "clinvar_variation_id": first["gofcards_accession_id"] or None,
                    },
                    "evidence": evidence_records.get(source_key, []),
                },
                "assemblies": {},
            })

            for assembly, arows in rows.groupby("VEP_assembly", sort=True):
                if not assembly:
                    continue
                head = arows.iloc[0]
                chrom, pos = head[f"{assembly}_chrom"], head[f"{assembly}_pos"]
                if not chrom:
                    continue
                # One variant has exactly one position per build, so this is a
                # single object rather than a map.
                block = variant["assemblies"].setdefault(assembly, {
                    "genomic": {
                        "chrom": chrom, "pos": int(pos),
                        "ref": head[f"{assembly}_ref"], "alt": head[f"{assembly}_alt"],
                        "status": head["hg19_vcf_status"] if assembly == "hg19"
                                  else head["hg38_refalt_status"],
                    },
                    "transcripts": {},
                })
                for _, r in arows.iterrows():
                    transcript = r["VEP_transcript"]
                    if not transcript:
                        continue
                    hgvsc = r["HGVSc"].split(":")[-1] if r["HGVSc"] else ""
                    hgvsp = r["HGVSp"].split(":")[-1] if r["HGVSp"] else ""
                    if not hgvsc:
                        # Nothing to key on for the transcript route; the variant
                        # stays reachable through its coordinates.
                        continue
                    view = block["transcripts"].setdefault(
                        transcript, {"by_hgvsc": {}, "by_hgvsp": {}})
                    view["by_hgvsc"][hgvsc] = {
                        "hgvsp": hgvsp or None,
                        "consequence": r["consequence"],
                        "canonical": r["canonical_transcript"] == "YES",
                        "mane_select": r["mane_select"] or None,
                    }
                    if hgvsp:
                        coding = view["by_hgvsp"].setdefault(hgvsp, [])
                        if hgvsc not in coding:
                            coding.append(hgvsc)

    return {
        "metadata": {
            "source": "GoFCards", "mechanism": "GOF",
            "derived_on": date.today().isoformat(),
            "structure": "genes -> alleles (normalized GRCh37 chrom|pos|ref|alt) "
                         "-> assemblies -> transcripts",
            **provenance,
        },
        "genes": genes,
    }


def build_cache(core: pd.DataFrame, evidence: pd.DataFrame, evidence_records: dict[str, list[dict]],
                vep_tables: list[pd.DataFrame],
                hg19_coords: pd.DataFrame, hg38_coords: pd.DataFrame,
                unmapped: set[str], reviews: dict, summary: dict[str, dict],
                resolver: HgncResolver, provenance: dict, out_tsv: Path) -> dict:
    """Join every layer, decide eligibility, and write the compact cache."""
    annotation = pd.concat([t for t in vep_tables if not t.empty], ignore_index=True)
    frame = core.merge(evidence, on="allele_key", how="left")

    # bcftools norm left-aligns and trims, so a tandem duplication written
    # longhand as TGAT->TGATTGAT becomes C->CTGAT. The authoritative coordinates
    # for BOTH builds are therefore the ones read back from the normalized VCFs.
    # Taking them from only one build would describe the same variant two
    # different ways depending on the assembly.
    for coords, prefix in ((hg19_coords, "hg19"), (hg38_coords, "hg38")):
        columns = [f"{prefix}_{field}" for field in ("chrom", "pos", "ref", "alt")]
        frame = frame.drop(columns=[c for c in columns if c in frame.columns])
        if not coords.empty:
            frame = frame.merge(coords, on="vcf_id", how="left")
        for column in columns:
            if column not in frame.columns:
                frame[column] = ""
    frame = frame.fillna("")
    merged = annotation.merge(frame, on="vcf_id", how="right").fillna("")

    merged["GoFCards_HGNC_Symbol"] = merged["GoFCards_source_symbol"].map(resolver.resolve)
    merged["VEP_HGNC_Symbol"] = merged["VEP_HGNC_Symbol"].map(resolver.resolve)
    merged["HGNC_Symbol"] = merged["GoFCards_HGNC_Symbol"]

    # A row with no curated gene cannot be attributed to anything, so it never
    # reaches the cache.
    merged = merged.loc[merged["GoFCards_HGNC_Symbol"] != ""].copy()

    # Gene provenance is a property of the ALLELE, not of an individual
    # transcript row. VEP reports every feature overlapping a locus, so a variant
    # in JAK2 also yields rows for neighbouring features like INSL6. Those rows
    # say nothing about JAK2 and must not be allowed to condemn it. The question
    # to ask per allele is: is the curated gene present at these coordinates?
    curated_present = set(
        merged.loc[merged["VEP_HGNC_Symbol"] == merged["GoFCards_HGNC_Symbol"], "vcf_id"])
    any_symbol = set(merged.loc[merged["VEP_HGNC_Symbol"] != "", "vcf_id"])

    # For a concordant allele, only the curated gene's own transcripts describe
    # it; other features' rows are dropped rather than retained as noise.
    keep = (~merged["vcf_id"].isin(curated_present)) | \
           (merged["VEP_HGNC_Symbol"] == merged["GoFCards_HGNC_Symbol"])
    dropped = int((~keep).sum())
    merged = merged.loc[keep].copy()
    log(f"dropped {dropped} transcript rows belonging to other genes at the same loci")

    # Status describes the row; eligibility describes the allele. A row carrying
    # no VEP symbol is not itself evidence of a disagreement, so it is only ever
    # a collateral casualty of one.
    def row_status(row: Any) -> str:
        if not row["VEP_HGNC_Symbol"]:
            return "source_gene_only"
        return "gene_concordant" if row["VEP_HGNC_Symbol"] == row["GoFCards_HGNC_Symbol"] \
            else "gene_discordant"

    merged["gene_match_status"] = merged.apply(row_status, axis=1)

    def eligibility(row: Any) -> tuple[str, str]:
        if row["vcf_id"] in curated_present:
            return ELIGIBLE, ""
        if row["vcf_id"] not in any_symbol:
            return ELIGIBLE, ""            # VEP placed no gene here; coordinates only
        if row["gene_match_status"] == "gene_discordant":
            return QUARANTINE_GENE, "curated_gene_absent_from_locus"
        return QUARANTINE_SIBLING, "sibling_of_gene_discordant_assertion"

    verdict = merged.apply(eligibility, axis=1, result_type="expand")
    merged["match_eligibility"], merged["reject_reason"] = verdict[0], verdict[1]

    # Polarity review applies only to rows that are otherwise eligible; a
    # gene-discordant row keeps its gene-provenance verdict.
    if reviews:
        def apply_review(row: Any) -> tuple[str, str]:
            if row["match_eligibility"] != ELIGIBLE:
                return row["match_eligibility"], row["reject_reason"]
            review = reviews.get((row["GoFCards_HGNC_Symbol"], row["allele_key"]))
            if not review or clean_text(review.get("gof_eligibility")).lower() != "quarantined":
                return row["match_eligibility"], row["reject_reason"]
            mechanism = clean_text(review.get("reviewed_mechanism")).upper()
            state = REVIEWED_QUARANTINE.get(mechanism, QUARANTINE_REVIEW_REQUIRED)
            reason = clean_text(review.get("reason_code")) or "reviewed_non_gof"
            return state, f"{reason}:{mechanism or 'unspecified'}"
        applied = merged.apply(apply_review, axis=1, result_type="expand")
        merged["match_eligibility"], merged["reject_reason"] = applied[0], applied[1]

    # Coordinate keys are no longer precomputed here. Each one is only meaningful
    # once its assembly is named, so the exporter builds them inside the assembly
    # block they belong to.
    merged["liftover_status"] = merged["vcf_id"].map(
        lambda v: "unmapped" if v in unmapped else "mapped")

    merged["hg38_refalt_status"] = merged.apply(
        lambda r: "liftover_unmapped" if r["liftover_status"] == "unmapped"
        else ("lifted_ref_match" if r["hg38_chrom"] else "no_hg38_coordinate"), axis=1)
    # This field reaches the final ACMG output, and its name means what GoFCards
    # means by it: the ClinVar accession, i.e. the ClinVar VariationID.
    merged["gofcards_accession_id"] = merged["allele_key"].map(
        lambda key: clean_text((summary.get(key) or {}).get("Accession")))

    # Every transcript is kept. The transcript identifier is a key in the exported
    # structure, so collapsing transcripts that happen to share an HGVS string
    # would delete a lookup path: a MANE pair gives the same HGVSc from its
    # Ensembl and its RefSeq member, and PriVA may annotate with either one.
    log(f"retaining {len(merged)} transcript annotations across "
        f"{merged['allele_key'].nunique()} source records")

    output = _nest(merged, evidence_records, summary, resolver, provenance, unmapped)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_tsv, "wt", encoding="utf-8") as handle:
        json.dump(output, handle, ensure_ascii=False, sort_keys=True, separators=(",", ":"))

    variants = [v for gene in output["genes"].values() for v in gene["variants"].values()]
    blocks = [(name, b) for v in variants for name, b in v["assemblies"].items()]
    views = [view for _, b in blocks for view in b["transcripts"].values()]
    stats = {
        "genes": len(output["genes"]),
        "variants": len(variants),
        "eligible_variants": sum(1 for v in variants
                                 if v["record"]["eligibility"]["status"] == ELIGIBLE),
        "quarantined_variants": sum(1 for v in variants
                                    if v["record"]["eligibility"]["status"] != ELIGIBLE),
        "evidence_entries": sum(len(v["record"]["evidence"]) for v in variants),
        "variants_on_hg19": sum(1 for v in variants if "hg19" in v["assemblies"]),
        "variants_on_hg38": sum(1 for v in variants if "hg38" in v["assemblies"]),
        # Transcript route: assembly + transcript + HGVSc (HGVSp maps to HGVSc)
        "transcript_views": len(views),
        "hgvsc_keys": sum(len(v["by_hgvsc"]) for v in views),
        "hgvsp_keys": sum(len(v["by_hgvsp"]) for v in views),
        "liftover_unmapped": len(unmapped),
    }
    log("cache written: " + "; ".join(f"{k}={v}" for k, v in stats.items()))
    return stats


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--out-json", "--out-tsv", dest="out_json", required=True, type=Path,
                        help="gzipped nested JSON cache")
    parser.add_argument("--hg19-fasta", required=True, type=Path)
    parser.add_argument("--hg38-fasta", required=True, type=Path)
    parser.add_argument("--chain", required=True, type=Path, help="UCSC hg19ToHg38 chain file")
    parser.add_argument("--vep-cache-dir", required=True, type=Path)
    parser.add_argument("--vep-cache-version", default="")
    parser.add_argument("--hgnc-table", required=True, type=Path)
    parser.add_argument("--mechanism-review-tsv", type=Path, default=None)
    parser.add_argument("--public-excel-url", default=PUBLIC_EXCEL_URL)
    parser.add_argument("--api-timeout", type=int, default=300,
                        help="seconds per backend request; the endpoint routinely takes ~50s")
    parser.add_argument("--stats-json", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    work = args.workdir
    logs = work / "logs"
    work.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)

    # Step 1
    source_records, provenance = fetch_sources(work, args.public_excel_url,
                                               timeout_seconds=args.api_timeout)
    summary = fetch_summary_annotations(source_records, work / "gofcards_summary_cache.jsonl",
                                        timeout_seconds=args.api_timeout)

    # Step 2 - GRCh37 is the only assembly GoFCards asserts, so it is the only
    # one that can adjudicate the reference allele.
    source_vcf = build_source_vcf(source_records, args.hg19_fasta, work / "gofcards.hg19.vcf",
                                  work / "gofcards_reference_rejects.tsv")

    # Step 3
    norm19 = normalize_vcf(work / "gofcards.hg19.vcf", args.hg19_fasta,
                           work / "gofcards.hg19.norm.vcf", logs)

    # Step 4
    vep19 = parse_vep(run_vep(norm19, "GRCh37", args.vep_cache_dir, args.vep_cache_version,
                              args.hg19_fasta, work / "gofcards.hg19.vep.tsv", logs), "hg19")

    # Step 5
    lifted, unmapped = liftover_vcf(norm19, args.chain, args.hg38_fasta,
                                    work / "gofcards.hg38.vcf", logs)
    norm38 = normalize_vcf(lifted, args.hg38_fasta, work / "gofcards.hg38.norm.vcf", logs)

    # Step 6
    vep38 = parse_vep(run_vep(norm38, "GRCh38", args.vep_cache_dir, args.vep_cache_version,
                              args.hg38_fasta, work / "gofcards.hg38.vep.tsv", logs), "hg38")

    # Step 7
    provenance.update({
        "vep_cache_dir": str(args.vep_cache_dir),
        "vep_cache_version": str(args.vep_cache_version),
        "liftover_chain": str(args.chain),
    })
    # CrossMap produces the representative GRCh38 record. It is the only source
    # that supplies a reference and alternate allele as well as a position: the
    # GoFCards endpoint returns a position alone, so its coordinate could only
    # ever be completed by assuming the GRCh37 alleles carry over. The endpoint
    # position is used solely to recover variants CrossMap could not place, and
    # only when the GRCh38 sequence confirms the carried-over reference.
    lifted_coords = read_vcf_coordinates(norm38, "hg38")
    placed = set(lifted_coords["vcf_id"]) if not lifted_coords.empty else set()
    unplaced = source_vcf.loc[~source_vcf["vcf_id"].isin(placed)]
    recovered, endpoint_counts = endpoint_hg38_coordinates(
        unplaced, summary, args.hg38_fasta)
    hg38_coords = pd.concat([lifted_coords, recovered], ignore_index=True) \
        if not recovered.empty else lifted_coords
    log(f"GRCh38 coordinates: {len(lifted_coords)} from the liftover, "
        f"{len(recovered)} recovered from the GoFCards endpoint")
    provenance.update({f"hg38_{k}": v for k, v in endpoint_counts.items()})

    evidence_frame, evidence_records = _aggregate_evidence(source_records)
    stats = build_cache(
        core=source_vcf,
        evidence=evidence_frame,
        evidence_records=evidence_records,
        vep_tables=[vep19, vep38],
        hg19_coords=read_vcf_coordinates(norm19, "hg19"),
        hg38_coords=hg38_coords,
        unmapped=unmapped,
        reviews=load_mechanism_reviews(args.mechanism_review_tsv),
        summary=summary,
        resolver=HgncResolver(args.hgnc_table),
        provenance=provenance,
        out_tsv=args.out_json,
    )
    if args.stats_json:
        args.stats_json.write_text(json.dumps({**provenance, **stats}, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
