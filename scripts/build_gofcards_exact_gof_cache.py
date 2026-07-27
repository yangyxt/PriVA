#!/usr/bin/env python3
"""Build the PriVA GoFCards exact gain-of-function evidence cache.

GoFCards publishes GRCh37 coordinates only, so GRCh37 is the sole ground truth
here.  Everything else is derived from it:

  1. fetch_sources      public workbook + backend SNV/Indel tables
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
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import requests
from pyfaidx import Fasta

BACKEND_BASE = "https://java.genemed.tech/admin-api/backend/data/hg19"
PUBLIC_EXCEL_URL = "https://download.genemed.tech/upload/GainFunCards/gofcards_data_download.xlsx"
TABLE_ENDPOINTS = {
    "SNV": f"{BACKEND_BASE}/GainFunCards_SNV/geneSymbol/page",
    "Indel": f"{BACKEND_BASE}/GainFunCards_Indel/geneSymbol/page",
}
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
    "HGVSc,HGVSp,MANE_SELECT,MANE_PLUS_CLINICAL,CANONICAL,SYMBOL,Existing_variation"
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

def clean_text(value: Any) -> str:
    """Normalize a field to text, mapping every placeholder for "empty" to "".

    VEP writes `-` for an absent field and GoFCards writes `_`, `.` or `-`.
    Treating those as text would produce match keys like `hgvsp_key="-"`, which
    would collide across every variant that has no protein change.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "_", ".", "-"} else text


def normalize_chrom(value: Any) -> str:
    return re.sub(r"^chr", "", clean_text(value), flags=re.IGNORECASE)


def normalize_allele(value: Any) -> str:
    text = clean_text(value)
    return "" if text in {"-", "."} else text.upper()


def normalize_int(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    text = str(value).strip()
    return text[:-2] if text.endswith(".0") else text


def locus_key_of(row: Any) -> str:
    """Type-free allele identity: chrom|start|end|ref|alt, all GRCh37.

    This is the real identity.  The two GoFCards sources disagree about the
    `variant_type` label for multi-base substitutions -- the backend calls them
    SNV, the public parser calls them Indel -- so the label must never be part
    of the key used to join the two together.
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


def fetch_sources(workdir: Path, public_url: str, page_size: int = 5000,
                  timeout_seconds: int = 300) -> tuple[pd.DataFrame, dict]:
    """Download the public workbook, then enrich it from the backend tables.

    The public workbook is the source of record: it is the citable artifact, it
    is hash-stable, and it carries every coordinate and every evidence field.
    The backend tables contribute exactly one thing the workbook omits, the
    ANNOVAR protein change (`AAChange_refGene`), which is what lets us flag
    alleles where GoFCards' own protein numbering disagrees with ours.

    The backend endpoint is slow -- a single 2.4 MB page takes around 50 seconds
    -- so the timeout is generous by design.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    session = requests.Session()

    public_path = workdir / "gofcards_public.xlsx"
    log(f"Downloading public workbook -> {public_path}")
    resp = session.get(public_url, headers={"User-Agent": "Mozilla/5.0"}, timeout=120)
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

    records: list[dict] = []
    for variant_type, url in TABLE_ENDPOINTS.items():
        page, pages = 1, 1
        while page <= pages:
            log(f"Pulling backend {variant_type} table, page {page}")
            payload = {
                "reference": "hg19", "queryby": "geneSymbol", "terms": "",
                "page": page, "pageNo": page, "currentPage": page, "pageSize": page_size,
            }
            for attempt in range(1, 5):
                try:
                    reply = session.post(url, headers=API_HEADERS, json=payload,
                                         timeout=timeout_seconds)
                    reply.raise_for_status()
                    body = reply.json()
                    break
                except Exception as exc:
                    if attempt == 4:
                        raise RuntimeError(
                            f"GoFCards backend {variant_type} table page {page} failed "
                            f"after 4 attempts: {exc}") from exc
                    log(f"  attempt {attempt} failed ({exc}); retrying")
                    time.sleep(5 * attempt)
            data = body.get("data") or {}
            pages = int(data.get("pages") or 1)
            for rec in data.get("records") or []:
                rec["variant_type"] = variant_type
                records.append(rec)
            page += 1

    backend = pd.DataFrame(records)
    if backend.empty:
        raise RuntimeError("GoFCards backend returned no records")
    backend["locus_key"] = backend.apply(locus_key_of, axis=1)

    # Join on the type-free locus key: the two sources label multi-base
    # substitutions differently, so joining on allele_key would drop them.
    protein = (backend.loc[:, ["locus_key", "AAChange_refGene"]]
               .assign(AAChange_refGene=lambda d: d["AAChange_refGene"].map(clean_text))
               .query("AAChange_refGene != ''")
               .drop_duplicates("locus_key"))
    merged = public.merge(protein, on="locus_key", how="left").fillna({"AAChange_refGene": ""})

    provenance = {
        "public_excel_url": public_url,
        "public_excel_last_modified": resp.headers.get("Last-Modified", ""),
        "public_workbook_rows": int(len(public)),
        "public_unique_alleles": int(public["locus_key"].nunique()),
        "backend_rows": int(len(backend)),
        "backend_unique_alleles": int(backend["locus_key"].nunique()),
        "rows_with_source_protein_change": int((merged["AAChange_refGene"] != "").sum()),
    }
    log(f"public rows={len(public)} unique alleles={public['locus_key'].nunique()}; "
        f"backend rows={len(backend)}; "
        f"protein change attached to {(merged['AAChange_refGene'] != '').sum()} rows")
    return merged, provenance


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
            "existing_variation": clean_text(rec.get("Existing_variation")),
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
            "source_protein_change": join_unique(group.get("AAChange_refGene", pd.Series(dtype=str))),
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
            "source_protein_change": clean_text(row.get("AAChange_refGene")),
        } for _, row in group.iterrows()]
    return pd.DataFrame(join_fields), records


def _nest(merged: pd.DataFrame, evidence_records: dict[str, list[dict]],
          resolver: HgncResolver, provenance: dict, unmapped: set[str]) -> dict:
    """Fold the joined table into gene -> allele -> assembly -> transcript.

    The nesting follows what actually owns each fact.  Evidence and the
    eligibility verdict belong to the allele and are assembly-independent, so
    they sit above the assembly level and are stored once.  Coordinates are
    shared by every transcript of an allele in one assembly, so they sit at the
    assembly level.  Only HGVS and consequence genuinely vary per transcript.

    An allele that could not be lifted over simply has no ``hg38`` entry, which
    states the situation more plainly than empty coordinate fields would.
    """
    genes: dict[str, dict] = {}
    for (symbol, allele_id), group in merged.groupby(["HGNC_Symbol", "hg19_vcf_key"], sort=True):
        if not symbol or not allele_id:
            continue
        first = group.iloc[0]
        gene = genes.setdefault(symbol, {"hgnc_id": resolver.hgnc_id(symbol), "alleles": {}})

        assemblies: dict[str, dict] = {}
        for assembly, rows in group.groupby("VEP_assembly", sort=True):
            if not assembly:
                continue
            head = rows.iloc[0]
            chrom, pos = head[f"{assembly}_chrom"], head[f"{assembly}_pos"]
            if not chrom:
                continue
            transcripts = {}
            for _, r in rows.iterrows():
                if not r["VEP_transcript"]:
                    continue
                transcripts[r["VEP_transcript"]] = {
                    "hgvsc": r["HGVSc"].split(":")[-1] if r["HGVSc"] else "",
                    "hgvsp": r["HGVSp"].split(":")[-1] if r["HGVSp"] else "",
                    "consequence": r["consequence"],
                    "mane_select": r["mane_select"] or None,
                    "canonical": r["canonical_transcript"] == "YES",
                }
            assemblies[assembly] = {
                "coordinates": {"chrom": chrom, "pos": int(pos),
                                "ref": head[f"{assembly}_ref"], "alt": head[f"{assembly}_alt"]},
                "status": head["hg19_vcf_status"] if assembly == "hg19" else head["hg38_refalt_status"],
                "transcripts": transcripts,
            }

        source_key = first["allele_key"]
        gene["alleles"][allele_id] = {
            "source_allele": {
                "variant_type_label": source_key.split("|")[0],
                "chrom": first["hg19_chrom"], "start": first["source_hg19_pos"],
                "ref": first["source_hg19_ref"], "alt": first["source_hg19_alt"],
                "allele_key": source_key,
            },
            "eligibility": {
                "status": first["match_eligibility"],
                "gene_match_status": first["gene_match_status"],
                "reason": first["reject_reason"] or None,
                "vep_symbol": first["VEP_HGNC_Symbol"] or None,
            },
            "liftover_status": "unmapped" if first["vcf_id"] in unmapped else "mapped",
            # How the GRCh37 VCF record was constructed at step 2. This describes
            # the allele itself, so it must survive even when VEP produced no
            # annotation for an assembly.
            "vcf_construction": first["hg19_vcf_status"],
            "rsid": first["gofcards_accession_id"] or None,
            "source_protein_change": first["source_protein_change"] or None,
            "protein_change_agreement": first["match_status"],
            "match_keys": {
                "hgvsp": sorted({r for r in group["hgvsp_key"] if r}),
                "hgvsc": sorted({r for r in group["hgvsc_key"] if r}),
                "hg19_vcf": allele_id,
                "hg38_vcf": first["hg38_vcf_key"] or None,
                "hg19_genomic": first["hg19_genomic_key"] or None,
            },
            "assemblies": assemblies,
            "evidence": evidence_records.get(source_key, []),
        }

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
                hg38_coords: pd.DataFrame, unmapped: set[str], reviews: dict,
                resolver: HgncResolver, provenance: dict, out_tsv: Path) -> dict:
    """Join every layer, decide eligibility, and write the compact cache."""
    annotation = pd.concat([t for t in vep_tables if not t.empty], ignore_index=True)
    frame = core.merge(evidence, on="allele_key", how="left")
    if not hg38_coords.empty:
        frame = frame.merge(hg38_coords, on="vcf_id", how="left")
    for column in ("hg38_chrom", "hg38_pos", "hg38_ref", "hg38_alt"):
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

    merged["hg19_genomic_key"] = merged.apply(
        lambda r: "|".join([r["hg19_chrom"], r["source_hg19_pos"], r["source_hg19_ref"],
                            r["source_hg19_alt"]]) if r["hg19_chrom"] else "", axis=1)
    merged["hg19_vcf_key"] = merged.apply(
        lambda r: "|".join([r["hg19_chrom"], str(r["hg19_pos"]), r["hg19_ref"], r["hg19_alt"]])
        if r["hg19_chrom"] else "", axis=1)
    merged["hg38_vcf_key"] = merged.apply(
        lambda r: "|".join([r["hg38_chrom"], str(r["hg38_pos"]), r["hg38_ref"], r["hg38_alt"]])
        if r["hg38_chrom"] else "", axis=1)
    merged["hg38_genomic_key"] = merged["hg38_vcf_key"]
    merged["liftover_status"] = merged["vcf_id"].map(
        lambda v: "unmapped" if v in unmapped else "mapped")

    # Columns the downstream non-LOF and HPO cache builders read. The normalized
    # VCF representation is the only one produced here, so the *_vcf_* fields
    # mirror hg19_*/hg38_* rather than duplicating a second convention.
    merged["hg38_refalt_status"] = merged.apply(
        lambda r: "liftover_unmapped" if r["liftover_status"] == "unmapped"
        else ("lifted_ref_match" if r["hg38_chrom"] else "no_hg38_coordinate"), axis=1)
    for assembly in ("hg19", "hg38"):
        for field in ("pos", "ref", "alt"):
            merged[f"{assembly}_vcf_{field}"] = merged[f"{assembly}_{field}"]
    merged["raw_GoFCards_HGVS"] = merged["source_protein_change"]
    merged["match_key_types"] = merged.apply(
        lambda r: ";".join(
            ([f"hgvsp"] if r["hgvsp_key"] else [])
            + (["hg19_genomic", "hg19_vcf"] if r["hg19_vcf_key"] else [])
            + (["hg38_genomic", "hg38_vcf"] if r["hg38_vcf_key"] else [])
        ), axis=1)

    merged["source"] = "GoFCards"
    merged["mechanism"] = "GOF"
    merged["build"] = "hg19_and_hg38"
    merged["gofcards_variant_id"] = merged["allele_key"]
    merged["gofcards_accession_id"] = merged["existing_variation"].map(
        lambda v: next((p for p in str(v).split(",") if p.startswith("rs")), ""))
    merged["match_status"] = merged.apply(
        lambda r: "source_protein_change_absent" if not r["source_protein_change"]
        else ("protein_change_agrees" if r["hgvsp_key"] and r["hgvsp_key"].replace("%3D", "=")
              in r["source_protein_change"].upper().replace("P.", "")
              else "protein_change_differs"), axis=1)
    merged["derived_on"] = date.today().isoformat()

    # Transcripts that yield the same HGVSc and HGVSp are indistinguishable for
    # matching, so only one representative of each distinct key pair is kept.
    # MANE Select first, then canonical, so the retained row is the one a reader
    # would expect to see. This is not a coverage reduction: every distinct
    # protein and coding change an allele can present is still present.
    before = len(merged)
    merged = (merged
              .assign(_rank=lambda d: list(zip(d["mane_select"] == "",
                                               d["canonical_transcript"] != "YES",
                                               d["VEP_transcript"])))
              .sort_values("_rank")
              .drop_duplicates(["allele_key", "VEP_assembly", "hgvsp_key", "hgvsc_key"],
                               keep="first")
              .drop(columns=["_rank"]))
    log(f"collapsed {before} transcript rows to {len(merged)} distinct HGVS keys")

    output = _nest(merged, evidence_records, resolver, provenance, unmapped)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_tsv, "wt", encoding="utf-8") as handle:
        json.dump(output, handle, ensure_ascii=False, sort_keys=True, separators=(",", ":"))

    alleles = [a for gene in output["genes"].values() for a in gene["alleles"].values()]
    stats = {
        "genes": len(output["genes"]),
        "alleles": len(alleles),
        "annotation_entries": int(len(merged)),
        "eligible_alleles": sum(1 for a in alleles if a["eligibility"]["status"] == ELIGIBLE),
        "quarantined_alleles": sum(1 for a in alleles if a["eligibility"]["status"] != ELIGIBLE),
        "alleles_with_hgvsp": sum(1 for a in alleles if a["match_keys"]["hgvsp"]),
        "alleles_with_hg38": sum(1 for a in alleles if a["match_keys"]["hg38_vcf"]),
        "evidence_records": sum(len(a["evidence"]) for a in alleles),
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
    evidence_frame, evidence_records = _aggregate_evidence(source_records)
    stats = build_cache(
        core=source_vcf,
        evidence=evidence_frame,
        evidence_records=evidence_records,
        vep_tables=[vep19, vep38],
        hg38_coords=read_vcf_coordinates(norm38, "hg38"),
        unmapped=unmapped,
        reviews=load_mechanism_reviews(args.mechanism_review_tsv),
        resolver=HgncResolver(args.hgnc_table),
        provenance=provenance,
        out_tsv=args.out_json,
    )
    if args.stats_json:
        args.stats_json.write_text(json.dumps({**provenance, **stats}, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
