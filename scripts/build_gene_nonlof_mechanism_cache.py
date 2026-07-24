#!/usr/bin/env python3
"""Build PriVA's portable non-LOF pathogenic-mechanism cache.

The builder is called by ``scripts/install_utils.sh`` after the exact GoFCards
cache is available. It checks each upstream source according to that source's
refresh interval, records checksums, injects eligible ClinVar VCV evidence,
validates the complete result, and atomically publishes one canonical JSON.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd

from clinvar_vcv import (
    flatten_match_audit,
    load_exact_gofcards_lookup,
    stream_parse_clinvar_vcv,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CACHE_DIR = PROJECT_ROOT / "data" / "gene_pathogenic_mechanism"
DEFAULT_SHARED_RAW_DIR = DEFAULT_CACHE_DIR / "raw"
DEFAULT_GOFCARDS_EXACT_VARIANTS = (
    PROJECT_ROOT / "data" / "gofcards" / "gofcards_exact_gof_hgvsp.tsv.gz"
)
DEFAULT_HGNC_TABLE = PROJECT_ROOT / "data" / "hgnc" / "non_alt_loci_set.tsv"
NONLOF_ASSERTIONS_FILENAME = "gene_nonlof_mechanism_curated_assertions.json"
NONLOF_ASSERTIONS_SCHEMA_FILENAME = (
    "gene_nonlof_mechanism_curated_assertions.schema.json"
)
DEFAULT_NONLOF_ASSERTIONS_JSON = (
    DEFAULT_CACHE_DIR / "prepared" / NONLOF_ASSERTIONS_FILENAME
)
DEFAULT_OUTPUT_SCHEMA = (
    DEFAULT_CACHE_DIR
    / "schema"
    / NONLOF_ASSERTIONS_SCHEMA_FILENAME
)
NONLOF_SOURCE_MANIFEST_FILENAME = "nonlof_source_manifest.json"
NONLOF_SOURCE_MANIFEST_TSV_FILENAME = "nonlof_source_manifest.tsv"
NONLOF_RUN_SUMMARY_FILENAME = "nonlof_run_summary.json"

GOFCARDS_EXACT_COLUMNS = [
    "source",
    "mechanism",
    "build",
    "HGNC_Symbol",
    "VEP_assembly",
    "VEP_transcript",
    "feature_type",
    "consequence",
    "HGVSc",
    "HGVSp",
    "hgvsp_key",
    "match_status",
    "raw_GoFCards_HGVS",
    "GoFCards_transcript",
    "canonical_transcript",
    "hg19_chrom",
    "hg19_pos",
    "hg19_ref",
    "hg19_alt",
    "hg19_vcf_status",
    "hg38_chrom",
    "hg38_pos",
    "hg38_ref",
    "hg38_alt",
    "hg38_refalt_status",
    "gofcards_accession_id",
    "gofcards_variant_id",
    "disease",
    "pmids",
    "pscore",
    "function",
    "pathway",
    "derived_on",
    "allele_key",
    "hg19_genomic_key",
    "hg19_vcf_pos",
    "hg19_vcf_ref",
    "hg19_vcf_alt",
    "hg38_vcf_pos",
    "hg38_vcf_ref",
    "hg38_vcf_alt",
    "hg19_vcf_key",
    "hg38_genomic_key",
    "hg38_vcf_key",
    "match_key_types",
]

GOFCARDS_EXACT_MISSING_VALUES = {"", ".", "_", "-", "na", "n/a", "nan", "none"}

USER_AGENT = (
    "PriVA/0.1 (portable non-LOF pathogenic-mechanism cache)"
)

NON_LOF_MECHANISM_LABELS = {
    "GOF",
    "ACTIVATING",
    "DOMINANT_NEGATIVE",
    "TRIPLOSENSITIVITY",
    "INCREASED_DOSAGE",
}
PANELAPP_NON_LOF_LABEL = "PANELAPP_GREEN_NON_LOF_PATHO_HISTORY"
PROMPT_EXCEPTION_LABELS = NON_LOF_MECHANISM_LABELS | {PANELAPP_NON_LOF_LABEL}

@dataclass(frozen=True)
class SourceSpec:
    name: str
    url: str
    raw_filename: str
    refresh_days: int
    source_role: str
    official_update_note: str
    docs_url: str
    parser: str
    local_fallback_path: str = ""


SOURCES: list[SourceSpec] = [
    SourceSpec(
        name="g2p_ddg2p",
        url="https://www.ebi.ac.uk/gene2phenotype/api/panel/all/download",
        raw_filename="AllG2P.csv",
        refresh_days=30,
        source_role="gene-disease molecular mechanism, allelic requirement, confidence, PMIDs",
        official_update_note=(
            "G2P was redesigned as a Vue SPA in 2025. The old DDG2P.csv.gz "
            "static download URL is gone. Use the new API endpoint: "
            "/api/panel/{panel}/download (CSV, not gzipped). "
            "AllG2P covers all panels including DD, Eye, Skin, Cancer, Cardiac."
        ),
        docs_url="https://www.ebi.ac.uk/gene2phenotype/api/",
        parser="g2p_ddg2p",
    ),
    SourceSpec(
        name="clingen_dosage_gene_grch38",
        url="https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv",
        raw_filename="ClinGen_gene_curation_list_GRCh38.tsv",
        refresh_days=1,
        source_role="haploinsufficiency and triplosensitivity gene curation",
        official_update_note=(
            "ClinGen says dosage curated gene TSV files are updated daily; "
            "the FTP index is refreshed nightly."
        ),
        docs_url="https://search.clinicalgenome.org/kb/downloads",
        parser="clingen_dosage",
    ),
    SourceSpec(
        name="gofcards_gof_variants",
        url="https://download.genemed.tech/upload/GainFunCards/gofcards_data_download.xlsx",
        raw_filename="gofcards/gofcards_data_download.xlsx",
        refresh_days=180,
        source_role="curated gain-of-function variant database (3161 variants, 579 genes)",
        official_update_note=(
            "GoFCards (NAR 2025, 53:D976-D988, PMID:39578693) is a curated "
            "GOF variant database. Updates follow publication cycles."
        ),
        docs_url="https://academic.oup.com/nar/article/53/D1/D976/7907365",
        parser="gofcards",
    ),
]

PANELAPP_SPEC = SourceSpec(
    name="panelapp_all_panels",
    url="https://panelapp.genomicsengland.co.uk/api/v1/panels/",
    raw_filename="panelapp_all_panels.json",
    refresh_days=7,
    source_role="PanelApp mode of inheritance and mode of pathogenicity",
    official_update_note=(
        "PanelApp is curator/API-backed and does not publish one fixed bulk "
        "release cycle; panel versions expose version_created timestamps."
    ),
    docs_url="https://panelapp.genomicsengland.co.uk/api/v1/panels/",
    parser="panelapp",
)

CLINVAR_VCV_SPEC = SourceSpec(
    name="clinvar_vcv_weekly",
    url=(
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/weekly_release/"
        "ClinVarVCVRelease_00-latest_weekly.xml.gz"
    ),
    raw_filename="clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz",
    refresh_days=7,
    source_role=(
        "condition-specific germline classifications and submitted observation "
        "evidence for exact normalized GoFCards GOF alleles"
    ),
    official_update_note=(
        "ClinVar generates VCV XML weekly. Check the companion MD5 first and "
        "download the multi-gigabyte XML only when that checksum changes."
    ),
    docs_url="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/_README",
    parser="clinvar_vcv",
)

CLINVAR_XML_README_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/_README"
CLINVAR_XSD_INDEX_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/"

def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def iso_now() -> str:
    return utc_now().isoformat(timespec="seconds")


def parse_iso(value: str | None) -> datetime | None:
    if not value:
        return None
    try:
        return datetime.fromisoformat(value)
    except ValueError:
        return None


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path) -> str:
    digest = hashlib.md5()  # nosec: provenance checksum, not security
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_mtime_iso(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(
        timespec="seconds"
    )


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")
    tmp.replace(path)


class BuildLock:
    def __init__(self, path: Path, stale_hours: float) -> None:
        self.path = path
        self.stale_hours = stale_hours
        self.fd: int | None = None

    def __enter__(self) -> "BuildLock":
        if self.path.exists():
            age_hours = (time.time() - self.path.stat().st_mtime) / 3600
            if age_hours > self.stale_hours:
                self.path.unlink()
            else:
                raise SystemExit(
                    f"Lock exists: {self.path}. Another cache build may be running."
                )
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        os.write(self.fd, f"pid={os.getpid()}\nstarted_at={iso_now()}\n".encode())
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        if self.fd is not None:
            os.close(self.fd)
        try:
            self.path.unlink()
        except FileNotFoundError:
            pass


def source_is_fresh(
    path: Path, meta: dict[str, Any], refresh_days: int, force: bool
) -> bool:
    if force or not path.exists() or path.stat().st_size == 0:
        return False
    checked_at = parse_iso(meta.get("last_checked_at_utc"))
    if checked_at is None:
        checked_at = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc)
    age_days = (utc_now() - checked_at).total_seconds() / 86400
    return age_days < refresh_days


def urlopen_with_optional_proxy(
    request: urllib.request.Request,
    timeout: int,
    proxy_url: str = "",
):
    if proxy_url:
        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({"http": proxy_url, "https": proxy_url})
        )
        return opener.open(request, timeout=timeout)
    return urllib.request.urlopen(request, timeout=timeout)


def proxy_needs_curl(proxy_url: str) -> bool:
    return proxy_url.lower().startswith(("socks4://", "socks4h://", "socks5://", "socks5h://"))


def download_url_with_curl(
    url: str,
    out_path: Path,
    timeout: int,
    retries: int,
    proxy_url: str = "",
    max_time: int | None = None,
    resume: bool = False,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "curl",
        "-L",
        "--fail",
        "--show-error",
        "--silent",
        "--connect-timeout",
        str(timeout),
        "--max-time",
        str(max_time if max_time is not None else max(timeout * 4, timeout + 60)),
        "--retry",
        str(max(0, retries - 1)),
        "--user-agent",
        USER_AGENT,
        "--output",
        str(out_path),
    ]
    if resume:
        cmd.extend(["--continue-at", "-"])
    if proxy_url:
        cmd.extend(["--proxy", proxy_url])
    cmd.append(url)
    result = subprocess.run(cmd, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(f"curl failed for {url}: {detail}")


def download_url(
    url: str,
    out_path: Path,
    timeout: int,
    retries: int,
    proxy_url: str = "",
    download_tool: str = "auto",
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if download_tool == "curl" or (
        download_tool == "auto" and proxy_needs_curl(proxy_url)
    ):
        download_url_with_curl(
            url,
            out_path,
            timeout=timeout,
            retries=retries,
            proxy_url=proxy_url,
        )
        return
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urlopen_with_optional_proxy(
                request, timeout=timeout, proxy_url=proxy_url
            ) as response:
                with out_path.open("wb") as handle:
                    shutil.copyfileobj(response, handle)
            return
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            last_error = error
            if attempt < retries:
                time.sleep(min(30, 2 * attempt))
    raise RuntimeError(f"failed to download {url}: {last_error}")


def update_metadata_for_file(
    source: SourceSpec,
    raw_path: Path,
    status: str,
    old_sha256: str | None,
    error: str | None = None,
) -> dict[str, Any]:
    metadata = {
        "name": source.name,
        "url": source.url,
        "docs_url": source.docs_url,
        "raw_path": str(raw_path),
        "source_role": source.source_role,
        "official_update_note": source.official_update_note,
        "refresh_days": source.refresh_days,
        "parser": source.parser,
        "last_checked_at_utc": iso_now(),
        "status": status,
        "error": error or "",
    }
    if raw_path.exists() and raw_path.stat().st_size > 0:
        sha256 = sha256_file(raw_path)
        metadata.update(
            {
                "size_bytes": raw_path.stat().st_size,
                "sha256": sha256,
                "md5": md5_file(raw_path),
                "file_mtime_utc": file_mtime_iso(raw_path),
                "content_changed": bool(old_sha256 and old_sha256 != sha256),
            }
        )
        if status in {"downloaded", "updated"}:
            metadata["last_downloaded_at_utc"] = iso_now()
    return metadata


def fetch_static_source(
    source: SourceSpec,
    raw_dir: Path,
    previous_meta: dict[str, Any],
    force: bool,
    timeout: int,
    retries: int,
    proxy_url: str,
    download_tool: str,
) -> dict[str, Any]:
    raw_path = raw_dir / source.raw_filename
    old_meta = previous_meta.get(source.name, {})
    if old_meta.get("raw_path") != str(raw_path):
        old_meta = {}
    if source_is_fresh(raw_path, old_meta, source.refresh_days, force):
        return update_metadata_for_file(
            source=source,
            raw_path=raw_path,
            status="skipped_fresh",
            old_sha256=old_meta.get("sha256"),
        )

    tmp = raw_path.with_name(raw_path.name + ".download.tmp")
    try:
        download_url(
            source.url,
            tmp,
            timeout=timeout,
            retries=retries,
            proxy_url=proxy_url,
            download_tool=download_tool,
        )
        new_sha = sha256_file(tmp)
        old_sha = old_meta.get("sha256")
        if raw_path.exists() and old_sha == new_sha:
            tmp.unlink()
            status = "checked_same_checksum"
        else:
            tmp.replace(raw_path)
            status = "updated" if old_sha else "downloaded"
        return update_metadata_for_file(
            source=source,
            raw_path=raw_path,
            status=status,
            old_sha256=old_sha,
        )
    except Exception as error:  # noqa: BLE001 - preserve cache and report source error
        if tmp.exists():
            tmp.unlink()
        if raw_path.exists():
            meta = update_metadata_for_file(
                source=source,
                raw_path=raw_path,
                status="using_existing_after_download_error",
                old_sha256=old_meta.get("sha256"),
                error=str(error),
            )
            return meta
        if source.local_fallback_path:
            fb = Path(source.local_fallback_path)
            if fb.exists() and fb.stat().st_size > 0:
                shutil.copy2(fb, raw_path)
                return update_metadata_for_file(
                    source=source,
                    raw_path=raw_path,
                    status="copied_local_fallback_after_download_error",
                    old_sha256=old_meta.get("sha256"),
                    error=f"{error}; copied local fallback {fb}",
                )
        return update_metadata_for_file(
            source=source,
            raw_path=raw_path,
            status="missing_after_download_error",
            old_sha256=old_meta.get("sha256"),
            error=str(error),
        )


def _download_text(
    url: str,
    timeout: int,
    retries: int,
    proxy_url: str,
    download_tool: str,
) -> str:
    with tempfile.NamedTemporaryFile(prefix="clinvar.meta.", delete=False) as handle:
        path = Path(handle.name)
    try:
        download_url(
            url,
            path,
            timeout=timeout,
            retries=retries,
            proxy_url=proxy_url,
            download_tool=download_tool,
        )
        return path.read_text(encoding="utf-8")
    finally:
        path.unlink(missing_ok=True)


def _write_downloaded_text(path: Path, text: str) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".download.tmp")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(path)
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "md5": md5_file(path),
    }


def _latest_clinvar_xsd_url(index_text: str) -> tuple[str, str]:
    versions = re.findall(r"ClinVar_VCV_(\d+(?:\.\d+)*)\.xsd", index_text)
    if not versions:
        raise ValueError("ClinVar XSD index contains no ClinVar_VCV_<version>.xsd")
    version = max(versions, key=lambda value: tuple(int(x) for x in value.split(".")))
    filename = f"ClinVar_VCV_{version}.xsd"
    return version, CLINVAR_XSD_INDEX_URL + filename


def _parse_remote_md5(text: str) -> str:
    match = re.search(r"\b([0-9a-fA-F]{32})\b", text)
    if not match:
        raise ValueError(f"Could not parse ClinVar companion MD5: {text[:200]!r}")
    return match.group(1).lower()


def fetch_clinvar_vcv(
    raw_dir: Path,
    previous_meta: dict[str, Any],
    force: bool,
    timeout: int,
    retries: int,
    proxy_url: str,
    download_tool: str,
    max_download_seconds: int,
) -> dict[str, Any]:
    """Check the weekly MD5 and download VCV XML only when content changed."""
    source = CLINVAR_VCV_SPEC
    raw_path = raw_dir / source.raw_filename
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    old = previous_meta.get(source.name, {})
    if old.get("raw_path") != str(raw_path):
        old = {}

    last_remote_check = parse_iso(old.get("last_remote_md5_checked_at_utc"))
    if (
        not force
        and raw_path.exists()
        and raw_path.stat().st_size > 0
        and last_remote_check is not None
        and (utc_now() - last_remote_check).total_seconds()
        < source.refresh_days * 86400
    ):
        metadata = dict(old)
        metadata["status"] = "skipped_fresh"
        metadata["last_build_seen_at_utc"] = iso_now()
        return metadata

    try:
        remote_md5_text = _download_text(
            source.url + ".md5",
            timeout,
            retries,
            proxy_url,
            download_tool,
        )
        remote_md5 = _parse_remote_md5(remote_md5_text)

        readme_text = _download_text(
            CLINVAR_XML_README_URL,
            timeout,
            retries,
            proxy_url,
            download_tool,
        )
        xsd_index = _download_text(
            CLINVAR_XSD_INDEX_URL,
            timeout,
            retries,
            proxy_url,
            download_tool,
        )
        xsd_version, xsd_url = _latest_clinvar_xsd_url(xsd_index)
        xsd_text = _download_text(
            xsd_url,
            timeout,
            retries,
            proxy_url,
            download_tool,
        )
        format_dir = raw_path.parent / "format"
        readme_meta = _write_downloaded_text(format_dir / "NCBI_CLINVAR_XML_README", readme_text)
        xsd_meta = _write_downloaded_text(
            format_dir / f"ClinVar_VCV_{xsd_version}.xsd", xsd_text
        )

        local_md5 = (
            md5_file(raw_path)
            if raw_path.exists() and raw_path.stat().st_size > 0
            else ""
        )
        local_verified = local_md5 == remote_md5
        downloaded = False
        if not local_verified:
            part = raw_path.with_name(f"{raw_path.name}.{remote_md5}.part")
            download_url_with_curl(
                source.url,
                part,
                timeout=timeout,
                retries=retries,
                proxy_url=proxy_url,
                max_time=max_download_seconds,
                resume=True,
            )
            actual_md5 = md5_file(part)
            if actual_md5 != remote_md5:
                part.unlink(missing_ok=True)
                raise ValueError(
                    f"ClinVar VCV MD5 mismatch: expected {remote_md5}, got {actual_md5}"
                )
            part.replace(raw_path)
            for stale_part in raw_path.parent.glob(f"{raw_path.name}.*.part"):
                stale_part.unlink(missing_ok=True)
            downloaded = True

        sha256 = sha256_file(raw_path) if downloaded or not old.get("sha256") else old["sha256"]
        metadata = {
            "name": source.name,
            "url": source.url,
            "md5_url": source.url + ".md5",
            "docs_url": source.docs_url,
            "raw_path": str(raw_path),
            "source_role": source.source_role,
            "official_update_note": source.official_update_note,
            "refresh_days": source.refresh_days,
            "parser": source.parser,
            "last_checked_at_utc": iso_now(),
            "last_remote_md5_checked_at_utc": iso_now(),
            "last_downloaded_at_utc": iso_now()
            if downloaded
            else old.get("last_downloaded_at_utc", ""),
            "status": "updated" if downloaded and old else "downloaded"
            if downloaded
            else "checked_same_remote_md5",
            "error": "",
            "size_bytes": raw_path.stat().st_size,
            "sha256": sha256,
            "md5": remote_md5,
            "remote_md5": remote_md5,
            "remote_md5_text": remote_md5_text.strip(),
            "file_mtime_utc": file_mtime_iso(raw_path),
            "content_changed": bool(
                downloaded and old.get("md5") and old.get("md5") != remote_md5
            ),
            "format_metadata": {
                "readme_url": CLINVAR_XML_README_URL,
                "readme": readme_meta,
                "xsd_index_url": CLINVAR_XSD_INDEX_URL,
                "xsd_version": xsd_version,
                "xsd_url": xsd_url,
                "xsd": xsd_meta,
            },
        }
        (raw_path.with_suffix(raw_path.suffix + ".md5")).write_text(
            remote_md5_text, encoding="utf-8"
        )
        return metadata
    except Exception as error:  # noqa: BLE001 - retain last verified release
        if raw_path.exists() and raw_path.stat().st_size > 0:
            metadata = dict(old)
            metadata.update(
                {
                    "name": source.name,
                    "url": source.url,
                    "docs_url": source.docs_url,
                    "raw_path": str(raw_path),
                    "source_role": source.source_role,
                    "official_update_note": source.official_update_note,
                    "refresh_days": source.refresh_days,
                    "parser": source.parser,
                    "last_checked_at_utc": iso_now(),
                    "status": "using_existing_after_download_error",
                    "error": str(error),
                    "size_bytes": raw_path.stat().st_size,
                    "file_mtime_utc": file_mtime_iso(raw_path),
                }
            )
            return metadata
        return {
            "name": source.name,
            "url": source.url,
            "docs_url": source.docs_url,
            "raw_path": str(raw_path),
            "source_role": source.source_role,
            "official_update_note": source.official_update_note,
            "refresh_days": source.refresh_days,
            "parser": source.parser,
            "last_checked_at_utc": iso_now(),
            "status": "missing_after_download_error",
            "error": str(error),
        }


def carry_forward_source_metadata(
    source: SourceSpec,
    raw_dir: Path,
    previous_meta: dict[str, Any],
) -> dict[str, Any]:
    """Carry a non-ClinVar source through a ClinVar-only refresh unchanged."""
    raw_path = raw_dir / source.raw_filename
    old = dict(previous_meta.get(source.name, {}))
    old.update(
        {
            "name": source.name,
            "url": source.url,
            "docs_url": source.docs_url,
            "raw_path": str(raw_path),
            "source_role": source.source_role,
            "official_update_note": source.official_update_note,
            "refresh_days": source.refresh_days,
            "parser": source.parser,
            "status": "not_checked_clinvar_only_refresh",
            "last_build_seen_at_utc": iso_now(),
        }
    )
    if not raw_path.exists() or raw_path.stat().st_size == 0:
        raise FileNotFoundError(
            f"ClinVar-only refresh requires existing source file: {raw_path}"
        )
    return old


def fetch_json_with_curl(
    url: str,
    timeout: int,
    retries: int,
    proxy_url: str = "",
) -> Any:
    with tempfile.NamedTemporaryFile(prefix="panelapp.", suffix=".json", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        download_url_with_curl(
            url,
            tmp_path,
            timeout=timeout,
            retries=retries,
            proxy_url=proxy_url,
        )
        return json.loads(tmp_path.read_text(encoding="utf-8"))
    finally:
        try:
            tmp_path.unlink()
        except FileNotFoundError:
            pass


def fetch_json(
    url: str,
    timeout: int,
    retries: int,
    proxy_url: str = "",
    download_tool: str = "auto",
) -> Any:
    if download_tool == "curl" or (
        download_tool == "auto" and proxy_needs_curl(proxy_url)
    ):
        return fetch_json_with_curl(
            url,
            timeout=timeout,
            retries=retries,
            proxy_url=proxy_url,
        )
    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
    )
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urlopen_with_optional_proxy(
                request, timeout=timeout, proxy_url=proxy_url
            ) as response:
                return json.loads(response.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, OSError, json.JSONDecodeError) as error:
            last_error = error
            if attempt < retries:
                time.sleep(min(30, 2 * attempt))
    raise RuntimeError(f"failed to fetch JSON {url}: {last_error}")


def fetch_panelapp(
    raw_dir: Path,
    previous_meta: dict[str, Any],
    force: bool,
    timeout: int,
    retries: int,
    proxy_url: str,
    download_tool: str,
    max_panels: int | None,
) -> dict[str, Any]:
    source = PANELAPP_SPEC
    raw_path = raw_dir / source.raw_filename
    old_meta = previous_meta.get(source.name, {})
    if old_meta.get("raw_path") != str(raw_path):
        old_meta = {}
    if source_is_fresh(raw_path, old_meta, source.refresh_days, force):
        return update_metadata_for_file(
            source=source,
            raw_path=raw_path,
            status="skipped_fresh",
            old_sha256=old_meta.get("sha256"),
        )

    try:
        panels: list[dict[str, Any]] = []
        next_url: str | None = source.url
        while next_url:
            payload = fetch_json(
                next_url,
                timeout=timeout,
                retries=retries,
                proxy_url=proxy_url,
                download_tool=download_tool,
            )
            panels.extend(payload.get("results", []))
            next_url = payload.get("next")
            if max_panels is not None and len(panels) >= max_panels:
                panels = panels[:max_panels]
                break

        detailed: list[dict[str, Any]] = []
        for index, panel in enumerate(panels, start=1):
            panel_id = panel.get("id")
            if not panel_id:
                continue
            detail_url = (
                f"https://panelapp.genomicsengland.co.uk/api/v1/panels/"
                f"{panel_id}/?format=json"
            )
            try:
                detail = fetch_json(
                    detail_url,
                    timeout=timeout,
                    retries=retries,
                    proxy_url=proxy_url,
                    download_tool=download_tool,
                )
            except Exception as error:  # noqa: BLE001
                detail = {
                    **panel,
                    "genes": [],
                    "detail_fetch_error": str(error),
                    "detail_url": detail_url,
                }
            detailed.append(detail)
            if index % 50 == 0:
                print(f"panelapp_details_fetched={index}", file=sys.stderr)

        tmp = raw_path.with_name(raw_path.name + ".download.tmp")
        payload = {
            "downloaded_at_utc": iso_now(),
            "list_url": source.url,
            "panel_count": len(panels),
            "panels": detailed,
        }
        with tmp.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True)
            handle.write("\n")
        new_sha = sha256_file(tmp)
        old_sha = old_meta.get("sha256")
        if raw_path.exists() and old_sha == new_sha:
            tmp.unlink()
            status = "checked_same_checksum"
        else:
            tmp.replace(raw_path)
            status = "updated" if old_sha else "downloaded"
        return update_metadata_for_file(source, raw_path, status, old_sha)
    except Exception as error:  # noqa: BLE001
        if raw_path.exists():
            return update_metadata_for_file(
                source, raw_path, "using_existing_after_download_error", old_meta.get("sha256"), str(error)
            )
        return update_metadata_for_file(
            source, raw_path, "missing_after_download_error", old_meta.get("sha256"), str(error)
        )


def norm_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def resolve_column(df: pd.DataFrame, *candidates: str) -> str:
    """Return the first column in *df* whose normalized name matches a candidate."""
    col_map = {norm_key(c): c for c in df.columns}
    for candidate in candidates:
        actual = col_map.get(norm_key(candidate))
        if actual is not None:
            return actual
    return ""


def _val(row: pd.Series, col: str) -> str:
    """Extract a stripped string value from *row* at resolved column *col*."""
    if not col:
        return ""
    v = row.get(col)
    return str(v).strip() if pd.notna(v) else ""


def split_multi(value: str) -> list[str]:
    if not value:
        return []
    parts = re.split(r"\s*(?:\||;|,\s(?=[A-Z0-9]))\s*", value)
    return [part.strip() for part in parts if part and part.strip()]


def normalize_mechanism(raw: str, source: str = "") -> list[str]:
    text = raw.lower().replace("_", " ").replace("-", " ")
    labels: set[str] = set()

    if source == "PanelApp":
        if "do not cause" in text and any(
            term in text for term in ["loss of function", "loss function", "lof"]
        ):
            return [PANELAPP_NON_LOF_LABEL]
        return []

    if re.search(r"\bdominant\s+negative\b", text):
        labels.add("DOMINANT_NEGATIVE")
    if any(term in text for term in ["gain of function", "gain function", "activating"]):
        labels.add("GOF")
    if "activated" in text or "activation" in text or "increased activity" in text:
        labels.add("ACTIVATING")
    if any(term in text for term in ["loss of function", "loss function", "lof"]):
        labels.add("LOF")
    if any(term in text for term in ["absent gene product", "absence of gene product"]):
        labels.add("LOF")
    if any(term in text for term in ["reduced function", "decreased function", "partial loss"]):
        labels.add("HYPMORPHIC_LOF")
    if "haploinsufficiency" in text or "haplo insufficient" in text:
        labels.add("HAPLOINSUFFICIENCY")
    if "triplosensitivity" in text or "triplo sensitivity" in text:
        labels.add("TRIPLOSENSITIVITY")
    if any(term in text for term in ["increased gene dosage", "increased dosage"]):
        labels.add("INCREASED_DOSAGE")
    if any(term in text for term in ["decreased gene dosage", "decreased dosage"]):
        labels.add("HAPLOINSUFFICIENCY")
    if any(term in text for term in ["altered gene product", "altered protein"]):
        labels.add("ALTERED_PRODUCT")
    if any(term in text for term in ["all missense", "missense", "in frame", "inframe"]):
        labels.add("MISSENSE_OR_INFRAME")
    if not labels and raw.strip():
        labels.add("UNMAPPED_MECHANISM")
    return sorted(labels)


def mechanism_confidence_from_text(value: str) -> str:
    text = value.lower()
    if value.strip() == "3":
        return "high"
    if value.strip() == "2":
        return "moderate"
    if value.strip() == "1":
        return "limited"
    if any(term in text for term in ["definitive", "sufficient", "green", "strong"]):
        return "high"
    if any(term in text for term in ["moderate", "emerging"]):
        return "moderate"
    if any(term in text for term in ["limited", "amber", "supportive"]):
        return "limited"
    if any(term in text for term in ["refuted", "disputed", "red", "unlikely"]):
        return "conflicting_or_refuted"
    return ""


UNIFIED_COLUMNS = [
    "gene_symbol",
    "mechanism",
    "assertion_level",
    "source",
    "source_record_id",
    "disease",
    "inheritance",
    "confidence",
    "disease_confidence",
    "pmids",
    "patho_mode_raw",
    "assembly",
    "chromosome",
    "position",
    "ref",
    "alt",
    "hgvsc",
    "hgvsp",
    "consequence",
    "allele_key",
    "notes",
]


def make_unified_row(
    *,
    gene_symbol: str,
    mechanism: Iterable[str] = (),
    assertion_level: str = "gene_level",
    source: str,
    source_record_id: str = "",
    disease: str = "",
    inheritance: str = "",
    confidence: str = "",
    disease_confidence: str = "",
    pmids: str = "",
    patho_mode_raw: str = "",
    assembly: str = "",
    chromosome: str = "",
    position: str = "",
    ref: str = "",
    alt: str = "",
    hgvsc: str = "",
    hgvsp: str = "",
    consequence: str = "",
    allele_key: str = "",
    notes: str = "",
) -> dict[str, str]:
    return {
        "gene_symbol": gene_symbol.strip(),
        "mechanism": "|".join(sorted(set(mechanism))),
        "assertion_level": assertion_level,
        "source": source,
        "source_record_id": source_record_id,
        "disease": disease,
        "inheritance": inheritance,
        "confidence": confidence,
        "disease_confidence": disease_confidence,
        "pmids": pmids,
        "patho_mode_raw": patho_mode_raw,
        "assembly": assembly,
        "chromosome": chromosome,
        "position": position,
        "ref": ref,
        "alt": alt,
        "hgvsc": hgvsc,
        "hgvsp": hgvsp,
        "consequence": consequence,
        "allele_key": allele_key,
        "notes": notes,
    }


def parse_g2p(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=UNIFIED_COLUMNS)
    df = pd.read_csv(path, dtype=str, encoding="utf-8-sig", on_bad_lines="skip").fillna("")
    c_gene = resolve_column(df, "gene symbol", "gene_symbol")
    c_disease = resolve_column(df, "disease name", "disease_label", "disease")
    c_consequence = resolve_column(df, "molecular mechanism", "mutation consequence", "mutation_consequence")
    c_confidence = resolve_column(df, "confidence", "confidence category", "g2p relation label")
    c_inheritance = resolve_column(df, "allelic requirement", "allelic_requirement")
    c_pmids = resolve_column(df, "pmids", "publications")
    c_hgnc = resolve_column(df, "hgnc id", "hgnc_id")
    rows: list[dict[str, str]] = []
    for _, r in df.iterrows():
        gene = _val(r, c_gene)
        if not gene:
            continue
        consequence = _val(r, c_consequence)
        confidence = _val(r, c_confidence)
        rows.append(
            make_unified_row(
                gene_symbol=gene,
                mechanism=normalize_mechanism(consequence, "G2P"),
                source="G2P_DDG2P",
                source_record_id=_val(r, c_hgnc),
                disease=_val(r, c_disease),
                inheritance=_val(r, c_inheritance),
                confidence=mechanism_confidence_from_text(confidence),
                disease_confidence=confidence,
                pmids=_val(r, c_pmids),
                patho_mode_raw=consequence,
            )
        )
    return pd.DataFrame(rows, columns=UNIFIED_COLUMNS)


def parse_panelapp(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=UNIFIED_COLUMNS)
    payload = read_json(path)
    panels = payload.get("panels", [])
    rows: list[dict[str, str]] = []
    for panel in panels:
        panel_name = str(panel.get("name", ""))
        panel_id = str(panel.get("id", ""))
        panel_version = str(panel.get("version", ""))
        genes = panel.get("genes", []) or []
        if not isinstance(genes, list):
            continue
        for gene_entry in genes:
            gene_data = gene_entry.get("gene_data") or {}
            gene = (
                gene_data.get("gene_symbol")
                or gene_data.get("hgnc_symbol")
                or gene_entry.get("entity_name")
                or gene_entry.get("gene_symbol")
                or ""
            )
            gene = str(gene).strip()
            if not gene:
                continue
            mop = str(
                gene_entry.get("mode_of_pathogenicity")
                or gene_entry.get("mode_of_pathogenicity_description")
                or ""
            ).strip()
            moi = str(gene_entry.get("mode_of_inheritance") or "").strip()
            publications = gene_entry.get("publications") or []
            pmids = "|".join(str(x) for x in publications if x)
            confidence = str(
                gene_entry.get("confidence_level")
                or gene_entry.get("penetrance")
                or ""
            )
            panelapp_labels = normalize_mechanism(mop, "PanelApp")
            if confidence.strip() != "3" or PANELAPP_NON_LOF_LABEL not in panelapp_labels:
                continue
            rows.append(
                make_unified_row(
                    gene_symbol=gene,
                    mechanism=panelapp_labels,
                    source="PanelApp",
                    source_record_id=f"panel:{panel_id};version:{panel_version}",
                    disease=panel_name,
                    inheritance=moi,
                    confidence=mechanism_confidence_from_text(confidence),
                    disease_confidence=confidence,
                    pmids=pmids,
                    patho_mode_raw=mop,
                )
            )
    return pd.DataFrame(rows, columns=UNIFIED_COLUMNS)


def score_positive(score: str) -> bool:
    return score.strip() in {"3", "2"}


def parse_clingen_dosage(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=UNIFIED_COLUMNS)
    with open(path, encoding="utf-8-sig") as fh:
        skip = 0
        for line in fh:
            if not line.startswith("#") or "\t" in line:
                break
            skip += 1
    df = pd.read_csv(path, sep="\t", dtype=str, encoding="utf-8-sig", skiprows=skip, on_bad_lines="skip").fillna("")
    if df.columns[0].startswith("#"):
        df.rename(columns={df.columns[0]: df.columns[0].lstrip("#").strip()}, inplace=True)
    c_gene = resolve_column(df, "Gene Symbol", "#Gene Symbol", "gene_symbol", "symbol")
    c_disease = resolve_column(df, "Disease Name", "Disease", "disease_label")
    c_hi_score = resolve_column(df, "Haploinsufficiency Score", "HI Score")
    c_hi_desc = resolve_column(df, "Haploinsufficiency Description", "HI Evidence Strength", "Haploinsufficiency")
    c_ts_score = resolve_column(df, "Triplosensitivity Score", "TS Score")
    c_ts_desc = resolve_column(df, "Triplosensitivity Description", "TS Evidence Strength", "Triplosensitivity")
    c_gene_id = resolve_column(df, "Gene ID", "HGNC ID")
    pmid_cols = [c for c in df.columns if "pmid" in norm_key(c)]
    rows: list[dict[str, str]] = []
    for _, r in df.iterrows():
        gene = _val(r, c_gene)
        if not gene:
            continue
        disease = _val(r, c_disease)
        hi_score = _val(r, c_hi_score)
        hi_desc = _val(r, c_hi_desc)
        ts_score = _val(r, c_ts_score)
        ts_desc = _val(r, c_ts_desc)
        pmid_values = [str(r[c]).strip() for c in pmid_cols if str(r.get(c, "")).strip()]
        pmids = "|".join(sorted(set("|".join(pmid_values).split("|")))) if pmid_values else ""
        gene_id = _val(r, c_gene_id)
        if hi_score or hi_desc:
            labels = ["HAPLOINSUFFICIENCY"] if score_positive(hi_score) else []
            rows.append(
                make_unified_row(
                    gene_symbol=gene,
                    mechanism=labels,
                    source="ClinGen_Dosage",
                    source_record_id=gene_id,
                    disease=disease,
                    confidence=mechanism_confidence_from_text(hi_desc),
                    disease_confidence=hi_desc,
                    pmids=pmids,
                    patho_mode_raw=f"HI_score={hi_score}; {hi_desc}".strip("; "),
                )
            )
        if ts_score or ts_desc:
            labels = ["TRIPLOSENSITIVITY"] if score_positive(ts_score) else []
            rows.append(
                make_unified_row(
                    gene_symbol=gene,
                    mechanism=labels,
                    source="ClinGen_Dosage",
                    source_record_id=gene_id,
                    disease=disease,
                    confidence=mechanism_confidence_from_text(ts_desc),
                    disease_confidence=ts_desc,
                    pmids=pmids,
                    patho_mode_raw=f"TS_score={ts_score}; {ts_desc}".strip("; "),
                )
            )
    return pd.DataFrame(rows, columns=UNIFIED_COLUMNS)


def _sanitize_gofcards_placeholder(val: str) -> str:
    return "" if val.strip() in ("_", "-", "nan", "NA", "N/A") else val


def _gofcards_public_allele_key(row: pd.Series) -> str:
    """Build the normalizer's stable source-allele key from the public row."""
    ref = _sanitize_gofcards_placeholder(str(row.get("ref", ""))).strip().upper()
    alt = _sanitize_gofcards_placeholder(str(row.get("alt", ""))).strip().upper()
    variant_type = "SNV" if len(ref) == 1 and len(alt) == 1 else "Indel"
    chrom = re.sub(r"^chr", "", str(row.get("chr", "")).strip(), flags=re.I)

    def integer_text(value: Any) -> str:
        text = str(value or "").strip()
        return text[:-2] if text.endswith(".0") else text

    return "|".join(
        [
            variant_type,
            chrom,
            integer_text(row.get("hg19start", "")),
            integer_text(row.get("hg19end", "")),
            ref,
            alt,
        ]
    )


def parse_gofcards(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=UNIFIED_COLUMNS)
    df = pd.read_excel(path, sheet_name="total3161", engine="openpyxl", dtype=str).fillna("")
    rows: list[dict[str, str]] = []
    for _, r in df.iterrows():
        gene = str(r.get("genesymbol", "")).strip()
        if not gene:
            continue
        pscore_raw = str(r.get("Pscore", "")).strip()
        try:
            pscore = float(pscore_raw)
        except (ValueError, TypeError):
            pscore = 0.0
        if pscore >= 5.0:
            conf = "high"
        elif pscore >= 3.0:
            conf = "moderate"
        else:
            conf = "limited"
        function_text = _sanitize_gofcards_placeholder(str(r.get("Function", "")))
        disorder = _sanitize_gofcards_placeholder(str(r.get("Disorder involved", "")))
        pmid = str(r.get("PMID", "")).strip()
        if pmid in ("", "nan"):
            pmid = ""
        chrom_raw = str(r.get("chr", "")).strip()
        chrom = f"chr{chrom_raw}" if chrom_raw and not chrom_raw.startswith("chr") else chrom_raw
        pos = str(r.get("hg19start", "")).strip()
        ref_allele = str(r.get("ref", "")).strip()
        alt_allele = str(r.get("alt", "")).strip()
        transcript = str(r.get("transcript", "")).strip()
        if transcript and chrom and pos:
            hgvsc = f"{transcript}:{chrom}:g.{pos}{ref_allele}>{alt_allele}"
        else:
            hgvsc = ""
        if len(ref_allele) == 1 and len(alt_allele) == 1 and ref_allele not in ("-", "") and alt_allele not in ("-", ""):
            consequence = "missense_variant"
        else:
            consequence = "indel"
        pathway = _sanitize_gofcards_placeholder(str(r.get("Pathways proteins involved", "")))
        animal = _sanitize_gofcards_placeholder(str(r.get("Animal model", "")))
        cell = _sanitize_gofcards_placeholder(str(r.get("Cell model", "")))
        notes_parts = [p for p in [
            f"pathway:{pathway}" if pathway else "",
            f"animal:{animal}" if animal else "",
            f"cell:{cell}" if cell else "",
        ] if p]
        rows.append(
            make_unified_row(
                gene_symbol=gene,
                mechanism=["GOF"],
                assertion_level="variant_level",
                source="GoFCards",
                source_record_id=str(r.get("Order numbe", "")),
                disease=disorder,
                confidence=conf,
                disease_confidence=f"Pscore={pscore_raw}" if pscore_raw else "",
                pmids=pmid,
                patho_mode_raw=function_text or "curated GOF variant",
                assembly="hg19",
                chromosome=chrom,
                position=pos,
                ref=ref_allele,
                alt=alt_allele,
                hgvsc=hgvsc,
                consequence=consequence,
                allele_key=_gofcards_public_allele_key(r),
                notes="; ".join(notes_parts),
            )
        )
    return pd.DataFrame(rows, columns=UNIFIED_COLUMNS)


EVIDENCE_PARSERS: dict[str, Callable[[Path], pd.DataFrame]] = {
    "g2p_ddg2p": parse_g2p,
    "panelapp": parse_panelapp,
    "clingen_dosage": parse_clingen_dosage,
    "gofcards": parse_gofcards,
}




def write_tsv(
    path: Path,
    columns: list[str],
    data: pd.DataFrame | Iterable[dict[str, str]],
) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    if not isinstance(data, pd.DataFrame):
        data = pd.DataFrame(data)
    for col in columns:
        if col not in data.columns:
            data[col] = ""
    data[columns].to_csv(tmp, sep="\t", index=False)
    tmp.replace(path)
    return len(data)


def parse_all_sources(raw_dir: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for spec in [*SOURCES, PANELAPP_SPEC]:
        parser_fn = EVIDENCE_PARSERS.get(spec.parser)
        if parser_fn is None:
            continue
        raw_path = raw_dir / spec.raw_filename
        if not raw_path.exists() and spec.local_fallback_path:
            raw_path = Path(spec.local_fallback_path)
        try:
            frames.append(parser_fn(raw_path))
        except Exception as exc:  # noqa: BLE001
            print(f"parse_error source={spec.name}: {exc}", file=sys.stderr)
            continue
    unified = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=UNIFIED_COLUMNS)
    if unified.empty:
        return unified
    keep = unified["mechanism"].apply(
        lambda v: bool(set(v.split("|")) & PROMPT_EXCEPTION_LABELS) if v.strip() else False
    )
    dropped = len(unified) - keep.sum()
    if dropped:
        print(f"filter: dropped {dropped} rows without curated non-LOF mechanism", file=sys.stderr)
    return unified[keep].reset_index(drop=True)


def write_source_manifest_tsv(path: Path, manifest: dict[str, Any]) -> int:
    rows = []
    for name, meta in sorted(manifest.get("sources", {}).items()):
        rows.append(
            {
                "source": name,
                "status": meta.get("status", ""),
                "url": meta.get("url", ""),
                "docs_url": meta.get("docs_url", ""),
                "raw_path": meta.get("raw_path", ""),
                "refresh_days": str(meta.get("refresh_days", "")),
                "last_checked_at_utc": meta.get("last_checked_at_utc", ""),
                "last_downloaded_at_utc": meta.get("last_downloaded_at_utc", ""),
                "file_mtime_utc": meta.get("file_mtime_utc", ""),
                "size_bytes": str(meta.get("size_bytes", "")),
                "sha256": meta.get("sha256", ""),
                "md5": meta.get("md5", ""),
                "remote_md5": meta.get("remote_md5", ""),
                "xsd_version": meta.get("format_metadata", {}).get(
                    "xsd_version", ""
                ),
                "xsd_path": meta.get("format_metadata", {})
                .get("xsd", {})
                .get("path", ""),
                "content_changed": str(meta.get("content_changed", "")),
                "error": meta.get("error", ""),
                "official_update_note": meta.get("official_update_note", ""),
            }
        )
    return write_tsv(
        path,
        [
            "source",
            "status",
            "url",
            "docs_url",
            "raw_path",
            "refresh_days",
            "last_checked_at_utc",
            "last_downloaded_at_utc",
            "file_mtime_utc",
            "size_bytes",
            "sha256",
            "md5",
            "remote_md5",
            "xsd_version",
            "xsd_path",
            "content_changed",
            "error",
            "official_update_note",
        ],
        rows,
    )


def load_hgnc_mapping(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    columns = {column.lower(): column for column in df.columns}

    def find_column(*names: str, required: bool = False) -> str:
        for name in names:
            if name.lower() in columns:
                return columns[name.lower()]
        if required:
            raise ValueError(
                f"{path} lacks required HGNC column; expected one of {names}"
            )
        return ""

    hgnc_col = find_column("hgnc_id", "HGNC ID", required=True)
    symbol_col = find_column("symbol", "Approved symbol", required=True)
    alias_col = find_column("alias_symbol", "Alias symbols")
    previous_col = find_column("prev_symbol", "Previous symbols")
    entrez_col = find_column(
        "entrez_id", "NCBI Gene ID(supplied by NCBI)"
    )
    ensembl_col = find_column(
        "ensembl_gene_id", "Ensembl ID(supplied by Ensembl)"
    )

    lookup: dict[str, dict[str, str]] = {}
    for _, r in df.iterrows():
        entry = {
            "hgnc_id": str(r.get(hgnc_col, "")).strip(),
            "symbol": str(r.get(symbol_col, "")).strip(),
            "entrez_id": str(r.get(entrez_col, "")).strip().removesuffix(".0"),
            "ensembl_id": str(r.get(ensembl_col, "")).strip(),
        }
        sym = entry["symbol"].upper()
        if sym:
            lookup[sym] = entry
        for alias in re.split(r"[|,]", str(r.get(alias_col, ""))):
            a = alias.strip().upper()
            if a and a not in lookup:
                lookup[a] = entry
        for prev in re.split(r"[|,]", str(r.get(previous_col, ""))):
            p = prev.strip().upper()
            if p and p not in lookup:
                lookup[p] = entry
    return lookup


def _split_pmids(raw: str) -> list[str]:
    if not raw or not raw.strip():
        return []
    parts = re.split(r"[|;,\s]+", raw.strip())
    return [p.strip() for p in parts if p.strip() and p.strip().replace("PMID:", "").isdigit()]


def _nonempty(d: dict) -> dict:
    return {k: v for k, v in d.items() if v is not None and v != "" and v != []}


def _clean_gofcards_exact_record(row: dict[str, Any]) -> dict[str, Any]:
    """Convert one exact-cache row to canonical JSON without placeholders."""
    record: dict[str, Any] = {}
    for key in GOFCARDS_EXACT_COLUMNS:
        value = str(row.get(key, "")).strip()
        if value.lower() in GOFCARDS_EXACT_MISSING_VALUES:
            continue
        if key == "pmids":
            pmids = _split_pmids(value)
            if pmids:
                record[key] = pmids
        elif key == "pscore":
            try:
                record[key] = float(value)
            except ValueError:
                record[key] = value
        elif key == "match_key_types":
            values = sorted({part.strip() for part in value.split(";") if part.strip()})
            if values:
                record[key] = values
        else:
            record[key] = value
    return record


def load_gofcards_exact_records(
    path: Path,
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, int]]:
    """Load every normalized GoFCards row, grouped by stable source allele."""
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [column for column in GOFCARDS_EXACT_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"GoFCards exact cache lacks columns: {missing}")

    by_allele: dict[str, list[dict[str, Any]]] = defaultdict(list)
    seen: set[str] = set()
    duplicate_rows = 0
    for row in frame.to_dict(orient="records"):
        record = _clean_gofcards_exact_record(row)
        allele_key = str(record.get("allele_key", ""))
        symbol = str(record.get("HGNC_Symbol", ""))
        if not allele_key or not symbol:
            raise ValueError("Every GoFCards exact-cache row requires allele_key and HGNC_Symbol")
        signature = json.dumps(record, sort_keys=True, separators=(",", ":"))
        if signature in seen:
            duplicate_rows += 1
            continue
        seen.add(signature)
        by_allele[allele_key].append(record)

    for records in by_allele.values():
        records.sort(
            key=lambda item: (
                str(item.get("HGNC_Symbol", "")),
                str(item.get("VEP_assembly", "")),
                str(item.get("VEP_transcript", "")),
                str(item.get("HGVSc", "")),
            )
        )
    return dict(by_allele), {
        "cache_rows": len(frame),
        "unique_cache_rows": len(seen),
        "unique_alleles": len(by_allele),
        "duplicate_cache_rows_removed": duplicate_rows,
    }


def _g2p_assertion(row: dict) -> dict:
    return _nonempty({
        "mechanism": row.get("mechanism", ""),
        "mechanism_raw": row.get("patho_mode_raw", ""),
        "disease": row.get("disease", ""),
        "inheritance": row.get("inheritance", ""),
        "confidence": row.get("confidence", ""),
        "pmids": _split_pmids(row.get("pmids", "")),
    })


def _panelapp_assertion(row: dict) -> dict:
    rid = row.get("source_record_id", "")
    panel_id = ""
    if "panel:" in rid:
        panel_id = rid.split("panel:")[1].split(";")[0]
    return _nonempty({
        "mechanism": row.get("mechanism", ""),
        "panel": row.get("disease", ""),
        "panel_id": panel_id,
        "inheritance": row.get("inheritance", ""),
        "confidence": row.get("confidence", ""),
    })


def _clingen_assertion(row: dict) -> dict:
    patho = row.get("patho_mode_raw", "")
    score = ""
    for part in patho.split(";"):
        p = part.strip()
        if p.startswith("TS_score=") or p.startswith("HI_score="):
            score = p.split("=")[1].strip()
            break
    return _nonempty({
        "mechanism": row.get("mechanism", ""),
        "score": score,
        "score_description": row.get("disease_confidence", ""),
        "pmids": _split_pmids(row.get("pmids", "")),
    })


def _gofcards_assertion(row: dict) -> dict:
    pscore_raw = row.get("disease_confidence", "")
    pscore = None
    if pscore_raw.startswith("Pscore="):
        try:
            pscore = float(pscore_raw.split("=")[1])
        except ValueError:
            pass
    transcript = ""
    hgvsc = row.get("hgvsc", "")
    if hgvsc and ":" in hgvsc:
        transcript = hgvsc.split(":")[0]
    return _nonempty({
        "mechanism": row.get("mechanism", ""),
        "source_record_id": row.get("source_record_id", ""),
        "allele_key": row.get("allele_key", ""),
        "disease": row.get("disease", ""),
        "chr": row.get("chromosome", ""),
        "pos": row.get("position", ""),
        "ref": row.get("ref", ""),
        "alt": row.get("alt", ""),
        "transcript": transcript,
        "consequence": row.get("consequence", ""),
        "pscore": pscore,
        "pmids": _split_pmids(row.get("pmids", "")),
        "function": row.get("patho_mode_raw", "") if row.get("patho_mode_raw", "") != "curated GOF variant" else "",
        "pathway": "",
        "animal_model": "",
        "cell_model": "",
    })


def _parse_gofcards_notes(assertion: dict, notes: str) -> None:
    for part in notes.split(";"):
        p = part.strip()
        if p.startswith("pathway:"):
            assertion["pathway"] = p[len("pathway:"):].strip()
        elif p.startswith("animal:"):
            val = p[len("animal:"):].strip()
            if val.upper() == "Y":
                assertion["animal_model"] = True
        elif p.startswith("cell:"):
            val = p[len("cell:"):].strip()
            if val.upper() == "Y":
                assertion["cell_model"] = True


_SOURCE_BUILDERS = {
    "G2P_DDG2P": ("gene_level", _g2p_assertion),
    "PanelApp": ("gene_level", _panelapp_assertion),
    "ClinGen_Dosage": ("gene_level", _clingen_assertion),
    "GoFCards": ("variant_level", _gofcards_assertion),
}


def build_nonlof_assertions_json(
    nonlof_df: pd.DataFrame,
    hgnc_map: dict[str, dict[str, str]],
    gofcards_exact_by_allele: dict[str, list[dict[str, Any]]] | None = None,
) -> dict[str, Any]:
    gofcards_exact_by_allele = gofcards_exact_by_allele or {}
    genes: dict[str, dict[str, Any]] = {}
    unmapped: list[str] = []

    for row in nonlof_df.to_dict(orient="records"):
        symbol = row.get("gene_symbol", "").strip()
        if not symbol:
            continue
        hgnc = hgnc_map.get(symbol.upper())
        if not hgnc:
            unmapped.append(symbol)
            hgnc_id = f"SYMBOL:{symbol}"
            entry_ids = {"symbol": symbol}
        else:
            hgnc_id = hgnc["hgnc_id"]
            entry_ids = {
                "symbol": hgnc["symbol"],
                "entrez_id": hgnc["entrez_id"],
                "ensembl_id": hgnc["ensembl_id"],
            }

        if hgnc_id not in genes:
            genes[hgnc_id] = {**_nonempty(entry_ids), "mechanisms": set(), "gene_level": [], "variant_level": []}

        gene = genes[hgnc_id]
        source = row.get("source", "")
        level_key, builder = _SOURCE_BUILDERS.get(source, ("gene_level", None))
        if builder is None:
            continue

        assertion = builder(row)
        if source == "GoFCards":
            _parse_gofcards_notes(assertion, row.get("notes", ""))
            allele_key = str(assertion.get("allele_key", ""))
            exact_rows = gofcards_exact_by_allele.get(allele_key, [])
            matching_rows = []
            for exact_row in exact_rows:
                exact_symbol = str(exact_row.get("HGNC_Symbol", ""))
                exact_hgnc = hgnc_map.get(exact_symbol.upper())
                exact_gene_id = (
                    exact_hgnc["hgnc_id"]
                    if exact_hgnc
                    else f"SYMBOL:{exact_symbol.upper()}"
                )
                if exact_gene_id == hgnc_id:
                    matching_rows.append(exact_row)
            if matching_rows:
                assertion["exact_normalization_status"] = "matched_gene_concordant"
                assertion["exact_normalized_variants"] = matching_rows
            elif exact_rows:
                assertion["exact_normalization_status"] = "gene_discordant_coordinate_collision"
                assertion["exact_cache_gene_symbols"] = sorted(
                    {
                        str(exact_row.get("HGNC_Symbol", ""))
                        for exact_row in exact_rows
                        if exact_row.get("HGNC_Symbol")
                    }
                )
            else:
                assertion["exact_normalization_status"] = "unmatched_public_source_allele"
            assertion = _nonempty(assertion)

        mech = row.get("mechanism", "")
        for m in mech.split("|"):
            if m.strip():
                gene["mechanisms"].add(m.strip())

        gene[level_key].append({source: assertion})

    result: dict[str, Any] = {}
    result["_meta"] = {
        "version": "2.0",
        "built_at": iso_now(),
        "total_genes": len(genes),
        "unmapped_symbols": sorted(set(unmapped)) if unmapped else [],
        "sources": {
            "G2P_DDG2P": {
                "level": "gene_level",
                "keys": ["mechanism", "mechanism_raw", "disease", "inheritance", "confidence", "pmids"],
                "mechanism_raw": "original G2P molecular mechanism text",
                "confidence_values": "definitive | strong | moderate | limited | conflicting_or_refuted",
            },
            "GoFCards": {
                "level": "variant_level",
                "keys": [
                    "mechanism", "source_record_id", "allele_key",
                    "exact_normalization_status", "exact_cache_gene_symbols",
                    "exact_normalized_variants", "disease", "chr", "pos", "ref",
                    "alt", "transcript", "consequence", "pscore", "pmids",
                    "function", "pathway", "animal_model", "cell_model",
                ],
                "exact_normalized_variant_keys": GOFCARDS_EXACT_COLUMNS,
                "raw_public_allele_fields": {
                    "assembly": "hg19",
                    "keys": ["chr", "pos", "ref", "alt"],
                    "note": (
                        "Applies only to the legacy top-level GoFCards source "
                        "fields; exact_normalized_variants identify hg19 and "
                        "hg38 fields explicitly"
                    ),
                },
                "exact_match_methods": [
                    "normalized gene symbol plus exact HGVSp",
                    "normalized gene symbol plus exact hg19/hg38 genomic allele",
                ],
                "pscore": "GoFCards pathogenicity score (float); >=5 high, >=3 moderate, else limited",
                "animal_model": "true if animal model evidence exists",
                "cell_model": "true if cell model evidence exists",
            },
            "PanelApp": {
                "level": "gene_level",
                "keys": ["mechanism", "panel", "panel_id", "inheritance", "confidence"],
                "note": "non-LOF history marker; LOF variants do NOT cause phenotype in this panel",
            },
            "ClinGen_Dosage": {
                "level": "gene_level",
                "keys": ["mechanism", "score", "score_description", "pmids"],
                "score": "ClinGen dosage score; 3=sufficient, 2=emerging evidence",
            },
            "ClinVar_VCV": {
                "level": "variant_level",
                "scope": (
                    "condition-specific germline RCV classifications and linked "
                    "SCV observations for exact normalized GoFCards alleles"
                ),
                "minimum_review_stars": 2,
                "eligible_review_statuses": [
                    "criteria provided, multiple submitters, no conflicts",
                    "reviewed by expert panel",
                    "practice guideline",
                ],
                "mechanism_caveat": (
                    "ClinVar does not assign GOF here; GOF context comes from the "
                    "exact matched GoFCards allele."
                ),
                "inheritance_caveat": (
                    "Observed zygosity and submitted mode of inheritance are retained "
                    "as evidence and do not set allelic_requirement."
                ),
            },
        },
        "mechanism_values": ["GOF", "DOMINANT_NEGATIVE", "TRIPLOSENSITIVITY", "PANELAPP_GREEN_NON_LOF_PATHO_HISTORY"],
    }

    for hgnc_id in sorted(genes):
        gene = genes[hgnc_id]
        gene["mechanisms"] = sorted(gene["mechanisms"])
        if not gene["gene_level"]:
            del gene["gene_level"]
        if not gene["variant_level"]:
            del gene["variant_level"]
        result[hgnc_id] = gene

    return result


def load_or_parse_clinvar_matches(
    xml_path: Path,
    exact_gofcards_path: Path,
    min_review_stars: int,
) -> tuple[list[dict[str, Any]], dict[str, int], str]:
    """Parse ClinVar matches without creating a second JSON evidence store."""
    if not xml_path.exists() or xml_path.stat().st_size == 0:
        raise FileNotFoundError(
            "ClinVar VCV is required for canonical schema 2.0; refusing to "
            f"replace the canonical JSON without it: {xml_path}"
        )
    if not exact_gofcards_path.exists():
        raise FileNotFoundError(
            f"Exact GoFCards lookup is required for ClinVar matching: {exact_gofcards_path}"
        )

    exact_lookup = load_exact_gofcards_lookup(exact_gofcards_path)
    matches, stats = stream_parse_clinvar_vcv(
        xml_path,
        exact_lookup,
        min_review_stars=min_review_stars,
    )
    return matches, stats, "parsed"


def inject_clinvar_matches(
    result: dict[str, Any],
    matches: list[dict[str, Any]],
    hgnc_map: dict[str, dict[str, str]],
) -> dict[str, int]:
    """Append ClinVar evidence to the same canonical HGNC-keyed JSON."""
    counts: Counter[str] = Counter()
    seen: set[tuple[str, str, str, str]] = set()
    for match in matches:
        payload = match.get("ClinVar_VCV", {})
        variation = payload.get("variation", {})
        for symbol in match.get("gene_symbols", []):
            hgnc = hgnc_map.get(str(symbol).upper())
            if not hgnc:
                counts["skipped_unmapped_gene"] += 1
                continue
            hgnc_id = hgnc["hgnc_id"]
            gene = result.get(hgnc_id)
            if gene is None:
                gene = {
                    **_nonempty(
                        {
                            "symbol": hgnc.get("symbol", symbol),
                            "entrez_id": hgnc.get("entrez_id", ""),
                            "ensembl_id": hgnc.get("ensembl_id", ""),
                        }
                    ),
                    "mechanisms": ["GOF"],
                    "variant_level": [],
                }
                result[hgnc_id] = gene
                counts["created_gene_entries"] += 1
            marker = (
                hgnc_id,
                str(variation.get("vcv_accession", "")),
                str(variation.get("component_variation_id", "")),
                str(variation.get("matched_component_context", "")),
            )
            if marker in seen:
                counts["deduplicated_gene_matches"] += 1
                continue
            seen.add(marker)
            gene.setdefault("variant_level", []).append(
                {"ClinVar_VCV": json.loads(json.dumps(payload))}
            )
            if "GOF" not in gene.setdefault("mechanisms", []):
                gene["mechanisms"].append("GOF")
                gene["mechanisms"].sort()
            counts["injected_gene_matches"] += 1
            counts["injected_rcv_conditions"] += len(
                payload.get("condition_assertions", [])
            )

    result["_meta"]["total_genes"] = len(result) - 1
    result["_meta"]["clinvar_vcv_integration"] = {
        "exact_match_method": "normalized GRCh37/GRCh38 VCF allele key",
        "minimum_review_stars": 2,
        "match_count": len(matches),
        **dict(sorted(counts.items())),
    }
    return dict(sorted(counts.items()))


def validate_canonical_json(result: dict[str, Any], schema_path: Path) -> None:
    """Validate the complete output before atomically replacing production."""
    try:
        from jsonschema import Draft202012Validator
    except ImportError as error:
        raise RuntimeError(
            "jsonschema is required to validate canonical schema 2.0 output"
        ) from error
    schema = read_json(schema_path)
    if not schema:
        raise FileNotFoundError(f"Canonical JSON Schema is missing: {schema_path}")
    validator = Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(result), key=lambda item: list(item.path))
    if errors:
        examples = []
        for error in errors[:10]:
            location = "/".join(str(part) for part in error.absolute_path) or "<root>"
            examples.append(f"{location}: {error.message}")
        raise ValueError(
            f"Canonical JSON failed schema validation ({len(errors)} errors): "
            + " | ".join(examples)
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument(
        "--shared-raw-dir",
        default=str(DEFAULT_SHARED_RAW_DIR),
        help=(
            "Shared raw-data directory used for downloads and parsing. Prepared "
            "project-specific outputs still go under --cache-dir/prepared. "
            "Expected filenames include AllG2P.csv, "
            "ClinGen_gene_curation_list_GRCh38.tsv, the weekly ClinVar VCV XML, "
            "GoFCards, and panelapp_all_panels.json. Use --shared-raw-dir '' to store raw "
            "files under --cache-dir/raw instead."
        ),
    )
    parser.add_argument(
        "--gofcards-exact-variants",
        type=Path,
        default=DEFAULT_GOFCARDS_EXACT_VARIANTS,
        help="Normalized hg19/hg38 exact GoFCards variant lookup.",
    )
    parser.add_argument(
        "--hgnc-table",
        type=Path,
        default=DEFAULT_HGNC_TABLE,
        help="PriVA HGNC complete-set table used for stable gene identifiers.",
    )
    parser.add_argument(
        "--output-schema",
        type=Path,
        default=DEFAULT_OUTPUT_SCHEMA,
        help="Draft 2020-12 schema required before canonical JSON replacement.",
    )
    parser.add_argument(
        "--clinvar-min-review-stars",
        type=int,
        choices=[2, 3, 4],
        default=2,
        help="Minimum condition-specific germline RCV review stars (default: 2).",
    )
    parser.add_argument(
        "--clinvar-max-download-seconds",
        type=int,
        default=86400,
        help="Maximum wall time for the resumable weekly VCV XML transfer.",
    )
    parser.add_argument(
        "--clinvar-only-refresh",
        action="store_true",
        help=(
            "Check/download only ClinVar VCV, carry other source metadata "
            "forward, then rebuild the canonical JSON from all existing raw files."
        ),
    )
    parser.add_argument("--force", action="store_true", help="Check/download all sources.")
    parser.add_argument(
        "--max-panelapp-panels",
        type=int,
        default=None,
        help="Debug cap for PanelApp panels fetched.",
    )
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument(
        "--proxy-url",
        default="",
        help=(
            "Optional proxy URL used only by this downloader, e.g. "
            "http://host:port, socks5://localhost:8888, or "
            "socks5h://localhost:8888. SOCKS proxies require curl, which is "
            "selected automatically when --download-tool auto is used."
        ),
    )
    parser.add_argument(
        "--download-tool",
        choices=["auto", "urllib", "curl"],
        default="auto",
        help=(
            "Downloader backend. auto uses Python urllib except when "
            "--proxy-url is a SOCKS proxy, where curl is required. The curl "
            "backend also supports installer runs on systems requiring a proxy."
        ),
    )
    parser.add_argument(
        "--validate-only",
        type=Path,
        metavar="JSON",
        help="Validate an existing canonical non-LOF JSON and exit.",
    )
    parser.add_argument("--stale-lock-hours", type=float, default=24.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.validate_only:
        canonical = read_json(args.validate_only.resolve())
        if not canonical:
            raise FileNotFoundError(args.validate_only.resolve())
        validate_canonical_json(canonical, args.output_schema.resolve())
        print(
            "schema_validation=passed "
            f"genes={len(canonical) - 1} json={args.validate_only.resolve()}"
        )
        return 0

    cache_dir = args.cache_dir.resolve()
    raw_dir = (
        Path(args.shared_raw_dir).expanduser().resolve()
        if args.shared_raw_dir.strip()
        else cache_dir / "raw"
    )
    raw_dir_is_shared = raw_dir != cache_dir / "raw"
    prepared_dir = cache_dir / "prepared"
    metadata_dir = cache_dir / "metadata"
    manifest_path = metadata_dir / NONLOF_SOURCE_MANIFEST_FILENAME
    lock_path = raw_dir / ".build.lock" if raw_dir_is_shared else cache_dir / ".build.lock"

    with BuildLock(lock_path, stale_hours=args.stale_lock_hours):
        previous_manifest = read_json(manifest_path)
        previous_sources = previous_manifest.get("sources", {})
        source_meta: dict[str, Any] = {}

        if args.clinvar_only_refresh:
            for source in [*SOURCES, PANELAPP_SPEC]:
                print(f"carrying_source={source.name}", file=sys.stderr)
                source_meta[source.name] = carry_forward_source_metadata(
                    source,
                    raw_dir,
                    previous_sources,
                )
        else:
            for source in SOURCES:
                print(f"checking_source={source.name}", file=sys.stderr)
                source_meta[source.name] = fetch_static_source(
                    source=source,
                    raw_dir=raw_dir,
                    previous_meta=previous_sources,
                    force=args.force,
                    timeout=args.timeout,
                    retries=args.retries,
                    proxy_url=args.proxy_url,
                    download_tool=args.download_tool,
                )

            print(f"checking_source={PANELAPP_SPEC.name}", file=sys.stderr)
            source_meta[PANELAPP_SPEC.name] = fetch_panelapp(
                raw_dir=raw_dir,
                previous_meta=previous_sources,
                force=args.force,
                timeout=args.timeout,
                retries=args.retries,
                proxy_url=args.proxy_url,
                download_tool=args.download_tool,
                max_panels=args.max_panelapp_panels,
            )

        print(f"checking_source={CLINVAR_VCV_SPEC.name}", file=sys.stderr)
        source_meta[CLINVAR_VCV_SPEC.name] = fetch_clinvar_vcv(
            raw_dir=raw_dir,
            previous_meta=previous_sources,
            force=args.force,
            timeout=args.timeout,
            retries=args.retries,
            proxy_url=args.proxy_url,
            download_tool=args.download_tool,
            max_download_seconds=args.clinvar_max_download_seconds,
        )

        manifest = {
            "schema_version": "2.0",
            "cache_dir": str(cache_dir),
            "raw_dir": str(raw_dir),
            "raw_dir_mode": "shared" if raw_dir_is_shared else "cache_local",
            "lock_path": str(lock_path),
            "proxy_mode": "explicit_proxy_url" if args.proxy_url else "environment_or_none",
            "download_tool": args.download_tool,
            "expected_raw_filenames": [
                source.raw_filename for source in SOURCES
            ]
            + [PANELAPP_SPEC.raw_filename, CLINVAR_VCV_SPEC.raw_filename],
            "built_at_utc": iso_now(),
            "sources": source_meta,
        }
        gofcards_exact_path = args.gofcards_exact_variants.resolve()
        hgnc_path = args.hgnc_table.resolve()
        schema_path = args.output_schema.resolve()
        build_inputs = {
            "builder_sha256": sha256_file(Path(__file__).resolve()),
            "clinvar_parser_sha256": sha256_file(
                Path(__file__).with_name("clinvar_vcv.py")
            ),
            "schema_sha256": sha256_file(schema_path),
            "hgnc_sha256": sha256_file(hgnc_path),
            "gofcards_exact_sha256": sha256_file(gofcards_exact_path),
            "source_sha256": {
                name: str(metadata.get("sha256", ""))
                for name, metadata in sorted(source_meta.items())
            },
        }
        manifest["build_inputs"] = build_inputs
        manifest_tsv_path = prepared_dir / NONLOF_SOURCE_MANIFEST_TSV_FILENAME
        json_path = prepared_dir / NONLOF_ASSERTIONS_FILENAME
        if (
            not args.force
            and json_path.exists()
            and previous_manifest.get("build_inputs") == build_inputs
        ):
            validate_canonical_json(read_json(json_path), schema_path)
            write_json(manifest_path, manifest)
            manifest_tsv_count = write_source_manifest_tsv(
                manifest_tsv_path,
                manifest,
            )
            summary_path = prepared_dir / NONLOF_RUN_SUMMARY_FILENAME
            run_summary = read_json(summary_path)
            run_summary.update(
                {
                    "last_checked_at_utc": manifest["built_at_utc"],
                    "status": "checked_inputs_unchanged",
                    "source_manifest_rows": manifest_tsv_count,
                }
            )
            write_json(summary_path, run_summary)
            print(
                json.dumps(
                    {
                        "status": "checked_inputs_unchanged",
                        "nonlof_assertions_json": str(json_path),
                        "schema_validation": "passed",
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
            return 0

        nonlof_df = parse_all_sources(raw_dir=raw_dir)

        hgnc_map = load_hgnc_mapping(hgnc_path)
        print(f"hgnc_mapping: {len(hgnc_map)} entries", file=sys.stderr)

        gofcards_exact_by_allele, gofcards_exact_stats = load_gofcards_exact_records(
            gofcards_exact_path
        )
        result = build_nonlof_assertions_json(
            nonlof_df,
            hgnc_map,
            gofcards_exact_by_allele=gofcards_exact_by_allele,
        )
        gofcards_assertions = [
            assertion["GoFCards"]
            for gene_id, gene in result.items()
            if gene_id != "_meta"
            for assertion in gene.get("variant_level", [])
            if "GoFCards" in assertion
        ]
        gofcards_status_counts = Counter(
            assertion["exact_normalization_status"]
            for assertion in gofcards_assertions
        )
        result["_meta"]["gofcards_exact_integration"] = {
            "exact_cache_path": str(gofcards_exact_path),
            "exact_cache_sha256": sha256_file(gofcards_exact_path),
            **gofcards_exact_stats,
            "public_assertions": len(gofcards_assertions),
            "emitted_exact_normalized_records": sum(
                len(assertion.get("exact_normalized_variants", []))
                for assertion in gofcards_assertions
            ),
            "assertion_status_counts": dict(sorted(gofcards_status_counts.items())),
        }
        clinvar_xml_path = raw_dir / CLINVAR_VCV_SPEC.raw_filename
        clinvar_matches, clinvar_stats, clinvar_parse_status = (
            load_or_parse_clinvar_matches(
                xml_path=clinvar_xml_path,
                exact_gofcards_path=args.gofcards_exact_variants.resolve(),
                min_review_stars=args.clinvar_min_review_stars,
            )
        )
        clinvar_injection = inject_clinvar_matches(result, clinvar_matches, hgnc_map)
        result["_meta"]["clinvar_vcv_integration"].update(
            {
                "source_vcv_md5": source_meta[CLINVAR_VCV_SPEC.name].get(
                    "md5", ""
                ),
                "gofcards_exact_sha256": sha256_file(
                    args.gofcards_exact_variants.resolve()
                ),
            }
        )
        clinvar_audit_count = write_tsv(
            prepared_dir / "clinvar_vcv_gofcards_matches.tsv",
            [
                "gene_symbol",
                "vcv_accession",
                "variation_id",
                "variation_name",
                "classification_scope",
                "matched_component_context",
                "rcv_accession",
                "condition",
                "clinical_significance",
                "review_status",
                "review_stars",
                "matched_scv_count",
            ],
            flatten_match_audit(clinvar_matches),
        )
        print(
            "clinvar_vcv: "
            f"parse_status={clinvar_parse_status}, matches={len(clinvar_matches)}, "
            f"audit_rows={clinvar_audit_count}, stats={json.dumps(clinvar_stats, sort_keys=True)}",
            file=sys.stderr,
        )
        validate_canonical_json(result, schema_path)
        print(
            f"json_schema_validation=passed schema={schema_path}",
            file=sys.stderr,
        )
        write_json(json_path, result)
        # Advance provenance only after the new canonical JSON has passed
        # schema validation and been atomically published. If parsing fails,
        # the prior manifest remains in place and the next run must retry.
        write_json(manifest_path, manifest)
        manifest_tsv_count = write_source_manifest_tsv(
            manifest_tsv_path,
            manifest,
        )
        n_genes = result["_meta"]["total_genes"]
        n_unmapped = len(result["_meta"].get("unmapped_symbols", []))
        print(f"json_output: {n_genes} genes, {n_unmapped} unmapped to {json_path}", file=sys.stderr)

        run_summary = {
            "built_at_utc": manifest["built_at_utc"],
            "last_checked_at_utc": manifest["built_at_utc"],
            "status": "rebuilt",
            "cache_dir": str(cache_dir),
            "total_genes": n_genes,
            "unmapped_symbols": n_unmapped,
            "source_manifest_rows": manifest_tsv_count,
            "clinvar_vcv": {
                "parse_status": clinvar_parse_status,
                "minimum_review_stars": args.clinvar_min_review_stars,
                "match_count": len(clinvar_matches),
                "audit_rows": clinvar_audit_count,
                "parser_stats": clinvar_stats,
                "injection_stats": clinvar_injection,
            },
            "outputs": {
                "source_manifest_json": str(manifest_path),
                "source_manifest_tsv": str(
                    prepared_dir / NONLOF_SOURCE_MANIFEST_TSV_FILENAME
                ),
                "nonlof_assertions_json": str(json_path),
                "clinvar_vcv_match_audit": str(
                    prepared_dir / "clinvar_vcv_gofcards_matches.tsv"
                ),
            },
        }
        write_json(prepared_dir / NONLOF_RUN_SUMMARY_FILENAME, run_summary)
        print(json.dumps(run_summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
