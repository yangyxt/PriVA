#!/usr/bin/env python3
"""Build PriVA's condition-resolved pathogenic-mechanism evidence table.

The builder keeps raw source files, records checksums and download/check
timestamps, and emits ``gene_pathogenic_mechanism_evidence.tsv``. That table is
a G2P/Orphadata build input to the integrated HPO condition cache; ACMG runtime
does not read it directly. ClinGen dosage and exact GoF assertions enter that
cache through their dedicated inputs. Exact GoF and dominant-negative allele
assertions are built by ``build_gene_nonlof_mechanism_cache.py`` instead.

The script requires pandas and is invoked by ``install_utils.sh`` as part of
the complete PriVA installer refresh. When a SOCKS proxy is requested,
downloads are delegated to curl because Python urllib does not natively
support SOCKS proxies.
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
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable

import pandas as pd

from hpo_penetrance import (
    normalize_penetrance_evidence,
    recognized_penetrance_hpo_ids,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CACHE_DIR = PROJECT_ROOT / "data" / "patho_mechanism"
DEFAULT_SHARED_RAW_DIR = DEFAULT_CACHE_DIR / "raw"
DEFAULT_DISEASE_SCOPE_REGISTRY = PROJECT_ROOT / "data" / "mondo" / "disease_scope.tsv.gz"

USER_AGENT = (
    "PriVA/0.1 "
    "(local pathogenic-mechanism cache; contact: yangyxt)"
)

PANELAPP_NON_LOF_LABEL = "PANELAPP_GREEN_NON_LOF_PATHO_HISTORY"

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
        source_role=(
            "gene-disease molecular mechanism, allelic requirement, "
            "condition-linked penetrance, confidence, PMIDs"
        ),
        official_update_note=(
            "G2P was redesigned as a Vue SPA in 2025. The old DDG2P.csv.gz "
            "static download URL is gone. Use the new API endpoint: "
            "/api/panel/{panel}/download (CSV, not gzipped). "
            "AllG2P covers all panels including DD, Eye, Skin, Cancer, Cardiac."
        ),
        docs_url="https://www.ebi.ac.uk/gene2phenotype/api/",
        parser="g2p_ddg2p",
        local_fallback_path="",
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
        local_fallback_path="",
    ),
    SourceSpec(
        name="orphadata_gene_disease",
        url="https://www.orphadata.com/data/xml/en_product6.xml",
        raw_filename="orphadata/en_product6.xml",
        refresh_days=30,
        source_role=(
            "condition-specific assessed germline gain-of-function and "
            "loss-of-function gene associations"
        ),
        official_update_note=(
            "Orphadata Product 6 contains disorder-gene relationship types. "
            "PriVA retains only assessed relationships whose type explicitly "
            "states a germline gain-of-function or loss-of-function mechanism."
        ),
        docs_url="https://www.orphadata.com/genes/",
        parser="orphadata",
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

CONDITION_EVIDENCE_SOURCE_NAMES = {"g2p_ddg2p", "orphadata_gene_disease"}
CONDITION_EVIDENCE_SPECS = [
    source for source in SOURCES if source.name in CONDITION_EVIDENCE_SOURCE_NAMES
]

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
        str(max(timeout * 4, timeout + 60)),
        "--retry",
        str(max(0, retries - 1)),
        "--user-agent",
        USER_AGENT,
        "--output",
        str(out_path),
    ]
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
    "source_condition_id",
    "mondo_id",
    "disease_scope",
    "priva_scope",
    "scope_review_status",
    "disease",
    "inheritance",
    "penetrance_raw",
    "penetrance_hpo_ids",
    "normalized_penetrance",
    "phenotype_terms",
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
    "notes",
]


def make_unified_row(
    *,
    gene_symbol: str,
    mechanism: Iterable[str] = (),
    assertion_level: str = "gene_level",
    source: str,
    source_record_id: str = "",
    source_condition_id: str = "",
    mondo_id: str = "",
    disease_scope: str = "",
    priva_scope: str = "",
    scope_review_status: str = "",
    disease: str = "",
    inheritance: str = "",
    penetrance_raw: str = "",
    penetrance_hpo_ids: str = "",
    normalized_penetrance: str = "",
    phenotype_terms: str = "",
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
    notes: str = "",
) -> dict[str, str]:
    return {
        "gene_symbol": gene_symbol.strip(),
        "mechanism": "|".join(sorted(set(mechanism))),
        "assertion_level": assertion_level,
        "source": source,
        "source_record_id": source_record_id,
        "source_condition_id": source_condition_id,
        "mondo_id": mondo_id,
        "disease_scope": disease_scope,
        "priva_scope": priva_scope,
        "scope_review_status": scope_review_status,
        "disease": disease,
        "inheritance": inheritance,
        "penetrance_raw": penetrance_raw,
        "penetrance_hpo_ids": penetrance_hpo_ids,
        "normalized_penetrance": normalized_penetrance,
        "phenotype_terms": phenotype_terms,
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
        "notes": notes,
    }


def parse_g2p(path: Path) -> pd.DataFrame:
    """Normalize G2P gene-condition mechanism assertions without losing IDs.

    G2P supplies three distinct identifiers that must not be conflated:

    * ``g2p id`` identifies the curated assertion itself;
    * ``disease mim`` identifies the source disease in OMIM; and
    * ``disease MONDO`` provides the cross-source disease identity used to join
      this assertion to HPO and Orphadata records.

    Earlier PriVA output stored the HGNC gene number as ``source_record_id`` and
    discarded both disease identifiers. That made condition-specific joins
    impossible and encouraged downstream gene-wide evidence transfer.
    """
    if not path.exists():
        return pd.DataFrame(columns=UNIFIED_COLUMNS)
    df = pd.read_csv(path, dtype=str, encoding="utf-8-sig", on_bad_lines="skip").fillna("")
    c_gene = resolve_column(df, "gene symbol", "gene_symbol")
    c_disease = resolve_column(df, "disease name", "disease_label", "disease")
    c_consequence = resolve_column(df, "molecular mechanism", "mutation consequence", "mutation_consequence")
    c_confidence = resolve_column(df, "confidence", "confidence category", "g2p relation label")
    c_inheritance = resolve_column(df, "allelic requirement", "allelic_requirement")
    c_penetrance_modifier = resolve_column(
        df,
        "cross cutting modifier",
        "cross_cutting_modifier",
    )
    c_phenotypes = resolve_column(df, "phenotypes", "phenotype_terms")
    c_pmids = resolve_column(df, "pmids", "publications")
    c_g2p_id = resolve_column(df, "g2p id", "g2p_id")
    c_disease_mim = resolve_column(df, "disease mim", "disease_mim")
    c_disease_mondo = resolve_column(df, "disease MONDO", "disease_mondo")
    rows: list[dict[str, str]] = []
    for _, r in df.iterrows():
        gene = _val(r, c_gene)
        if not gene:
            continue
        consequence = _val(r, c_consequence)
        confidence = _val(r, c_confidence)
        penetrance_raw = _val(r, c_penetrance_modifier)
        phenotype_terms = _val(r, c_phenotypes)
        penetrance_hpo_ids = recognized_penetrance_hpo_ids(phenotype_terms)
        rows.append(
            make_unified_row(
                gene_symbol=gene,
                mechanism=normalize_mechanism(consequence, "G2P"),
                source="G2P_DDG2P",
                source_record_id=_val(r, c_g2p_id),
                source_condition_id=(
                    f"OMIM:{_val(r, c_disease_mim)}"
                    if _val(r, c_disease_mim)
                    else ""
                ),
                mondo_id=_val(r, c_disease_mondo),
                disease=_val(r, c_disease),
                inheritance=_val(r, c_inheritance),
                penetrance_raw=penetrance_raw,
                penetrance_hpo_ids=";".join(penetrance_hpo_ids),
                normalized_penetrance=normalize_penetrance_evidence(
                    hpo_ids=penetrance_hpo_ids,
                    raw_values=penetrance_raw,
                ),
                phenotype_terms=phenotype_terms,
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
                notes="; ".join(notes_parts),
            )
        )
    return pd.DataFrame(rows, columns=UNIFIED_COLUMNS)


def parse_orphadata(path: Path) -> pd.DataFrame:
    """Stream condition-resolved germline mechanism assertions from Product 6.

    Orphadata includes many relationship types that are not interchangeable
    with an inherited disease mechanism. For example, a gene can be recorded
    as a biomarker, a susceptibility factor, part of a fusion, or the target of
    a somatic mutation. PriVA therefore accepts a relationship only when all of
    the following are true:

    1. Orphadata marks the relationship as ``Assessed``;
    2. the relationship explicitly says ``Disease-causing germline``; and
    3. the relationship text explicitly names gain or loss of function.

    The parser uses ``iterparse`` so the full Product 6 XML document is never
    materialized in memory. Each retained row preserves the ORPHA disease ID,
    disease name, relationship text, validation publications, and a stable
    source-record key. MONDO identity is added separately from PriVA's pinned
    disease-scope registry because Product 6 itself does not provide MONDO IDs.
    """
    if not path.exists():
        return pd.DataFrame(columns=UNIFIED_COLUMNS)

    rows: list[dict[str, str]] = []
    for _, disorder in ET.iterparse(path, events=("end",)):
        if disorder.tag != "Disorder":
            continue

        orpha_code = (disorder.findtext("./OrphaCode") or "").strip()
        disease = (disorder.findtext("./Name") or "").strip()
        if not orpha_code:
            disorder.clear()
            continue

        for association in disorder.findall(
            "./DisorderGeneAssociationList/DisorderGeneAssociation"
        ):
            status = (
                association.findtext("./DisorderGeneAssociationStatus/Name") or ""
            ).strip()
            relationship = (
                association.findtext("./DisorderGeneAssociationType/Name") or ""
            ).strip()
            relationship_lower = relationship.lower()
            if status != "Assessed" or "disease-causing germline" not in relationship_lower:
                continue

            mechanisms = normalize_mechanism(relationship, "Orphadata")
            explicit_mechanisms = [
                mechanism for mechanism in mechanisms if mechanism in {"GOF", "LOF"}
            ]
            if not explicit_mechanisms:
                continue

            gene = (association.findtext("./Gene/Symbol") or "").strip()
            if not gene:
                continue
            validation = (association.findtext("./SourceOfValidation") or "").strip()
            pmids = ";".join(
                dict.fromkeys(re.findall(r"(\d+)\s*\[PMID\]", validation))
            )
            source_condition_id = f"ORPHA:{orpha_code}"
            rows.append(
                make_unified_row(
                    gene_symbol=gene,
                    mechanism=explicit_mechanisms,
                    source="Orphadata",
                    source_record_id=(
                        f"{source_condition_id}|{gene}|{norm_key(relationship)}"
                    ),
                    source_condition_id=source_condition_id,
                    disease=disease,
                    confidence="high",
                    disease_confidence=status,
                    pmids=pmids,
                    patho_mode_raw=relationship,
                    notes=f"association_status:{status}",
                )
            )

        disorder.clear()

    return pd.DataFrame(rows, columns=UNIFIED_COLUMNS)


EVIDENCE_PARSERS: dict[str, Callable[[Path], pd.DataFrame]] = {
    "g2p_ddg2p": parse_g2p,
    "panelapp": parse_panelapp,
    "clingen_dosage": parse_clingen_dosage,
    "gofcards": parse_gofcards,
    "orphadata": parse_orphadata,
}


def attach_mondo_condition_identity(
    unified: pd.DataFrame,
    disease_scope_registry: Path | None,
) -> pd.DataFrame:
    """Attach stable MONDO identity and PriVA disease scope by source ID.

    The pinned disease-scope registry is keyed by a source identifier such as
    ``OMIM:615355`` or ``ORPHA:648``. It supplies the corresponding MONDO ID and
    the decision to include, review, or exclude that disease in a germline PriVA
    run. Source-provided MONDO IDs, currently supplied by G2P, remain primary.
    Registry values fill missing IDs, which is especially important for
    Orphadata, but do not silently replace a conflicting source assertion.

    Rows whose disease is absent from the registry remain usable as curated
    mechanism history. Their scope fields stay empty so downstream code can
    distinguish "not represented in HPO/MONDO scope data" from an explicit
    ``review`` or ``exclude`` decision.
    """
    output = unified.copy()
    for column in ("mondo_id", "disease_scope", "priva_scope", "scope_review_status"):
        if column not in output.columns:
            output[column] = ""
    if (
        output.empty
        or disease_scope_registry is None
        or not disease_scope_registry.exists()
    ):
        return output

    registry = pd.read_csv(
        disease_scope_registry,
        sep="\t",
        dtype=str,
        low_memory=False,
    ).fillna("")
    required = {
        "disease_id",
        "mondo_id",
        "disease_scope",
        "priva_scope",
        "scope_review_status",
    }
    missing = sorted(required - set(registry.columns))
    if missing:
        raise ValueError(
            f"disease scope registry missing columns: {', '.join(missing)}"
        )
    registry = registry.drop_duplicates("disease_id").set_index("disease_id")

    source_ids = output["source_condition_id"].fillna("").astype(str)
    mapped_mondo = source_ids.map(registry["mondo_id"]).fillna("")
    supplied_mondo = output["mondo_id"].fillna("").astype(str)
    conflicts = supplied_mondo.ne("") & mapped_mondo.ne("") & supplied_mondo.ne(mapped_mondo)
    if conflicts.any():
        print(
            "condition_identity: retained source MONDO ID for "
            f"{int(conflicts.sum())} registry conflicts",
            file=sys.stderr,
        )
    output["mondo_id"] = supplied_mondo.mask(supplied_mondo.eq(""), mapped_mondo)

    for column in ("disease_scope", "priva_scope", "scope_review_status"):
        mapped = source_ids.map(registry[column]).fillna("")
        current = output[column].fillna("").astype(str)
        output[column] = current.mask(current.eq(""), mapped)
    return output




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
    if (
        path.exists()
        and path.stat().st_size == tmp.stat().st_size
        and sha256_file(path) == sha256_file(tmp)
    ):
        tmp.unlink()
    else:
        tmp.replace(path)
    return len(data)


EVIDENCE_COLUMNS = [
    "gene_symbol",
    "source",
    "source_record_id",
    "source_condition_id",
    "mondo_id",
    "disease_scope",
    "priva_scope",
    "scope_review_status",
    "disease_label",
    "inheritance",
    "penetrance_raw",
    "penetrance_hpo_ids",
    "normalized_penetrance",
    "patho_mode_raw",
    "normalized_mechanisms",
    "mechanism_confidence",
    "disease_confidence",
    "pmids",
    "evidence_url",
    "panel_or_submitter",
    "phenotype_terms",
    "notes",
]


def source_url_for_evidence(source: str) -> str:
    for spec in [*SOURCES, PANELAPP_SPEC]:
        if source == "G2P_DDG2P" and spec.name == "g2p_ddg2p":
            return spec.url
        if source == "ClinGen_Dosage" and spec.name == "clingen_dosage_gene_grch38":
            return spec.url
        if source == "GoFCards" and spec.name == "gofcards_gof_variants":
            return spec.url
        if source == "PanelApp" and spec.name == "panelapp_all_panels":
            return spec.url
        if source == "Orphadata" and spec.name == "orphadata_gene_disease":
            return spec.url
    return ""


def to_gene_pathogenic_mechanism_evidence(unified: pd.DataFrame) -> pd.DataFrame:
    if unified.empty or "source" not in unified.columns:
        return pd.DataFrame(columns=EVIDENCE_COLUMNS)
    unified = unified.loc[
        unified["source"].isin({"G2P_DDG2P", "Orphadata"})
    ].copy()
    if unified.empty:
        return pd.DataFrame(columns=EVIDENCE_COLUMNS)
    evidence = pd.DataFrame(
        {
            "gene_symbol": unified.get("gene_symbol", ""),
            "source": unified.get("source", ""),
            "source_record_id": unified.get("source_record_id", ""),
            "source_condition_id": unified.get("source_condition_id", ""),
            "mondo_id": unified.get("mondo_id", ""),
            "disease_scope": unified.get("disease_scope", ""),
            "priva_scope": unified.get("priva_scope", ""),
            "scope_review_status": unified.get("scope_review_status", ""),
            "disease_label": unified.get("disease", ""),
            "inheritance": unified.get("inheritance", ""),
            "penetrance_raw": unified.get("penetrance_raw", ""),
            "penetrance_hpo_ids": unified.get("penetrance_hpo_ids", ""),
            "normalized_penetrance": unified.get("normalized_penetrance", ""),
            "patho_mode_raw": unified.get("patho_mode_raw", ""),
            "normalized_mechanisms": unified.get("mechanism", "").astype(str).str.replace("|", ";", regex=False),
            "mechanism_confidence": unified.get("confidence", ""),
            "disease_confidence": unified.get("disease_confidence", ""),
            "pmids": unified.get("pmids", ""),
            "evidence_url": unified.get("source", "").astype(str).map(source_url_for_evidence),
            "panel_or_submitter": "",
            "phenotype_terms": unified.get("phenotype_terms", ""),
            "notes": unified.get("notes", ""),
        }
    )
    return evidence[EVIDENCE_COLUMNS].fillna("")


def write_gene_pathogenic_mechanism_evidence(path: Path, unified: pd.DataFrame) -> int:
    return write_tsv(
        path,
        EVIDENCE_COLUMNS,
        to_gene_pathogenic_mechanism_evidence(unified),
    )


def parse_all_sources(
    raw_dir: Path,
    disease_scope_registry: Path | None = DEFAULT_DISEASE_SCOPE_REGISTRY,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for spec in CONDITION_EVIDENCE_SPECS:
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
    return attach_mondo_condition_identity(
        unified.reset_index(drop=True),
        disease_scope_registry,
    )


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
                "content_changed": str(meta.get("content_changed", "")),
                "official_update_note": meta.get("official_update_note", ""),
                "error": meta.get("error", ""),
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
            "content_changed",
            "official_update_note",
            "error",
        ],
        rows,
    )


DEFAULT_HGNC_MAP = DEFAULT_SHARED_RAW_DIR / "hgnc_gene_id_map.txt"


def load_hgnc_mapping(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    lookup: dict[str, dict[str, str]] = {}
    for _, r in df.iterrows():
        entry = {
            "hgnc_id": r.iloc[0].strip(),
            "symbol": r.iloc[1].strip(),
            "entrez_id": r.iloc[4].strip(),
            "ensembl_id": r.iloc[5].strip(),
        }
        sym = entry["symbol"].upper()
        if sym:
            lookup[sym] = entry
        for alias in r.iloc[2].split("|"):
            a = alias.strip().upper()
            if a and a not in lookup:
                lookup[a] = entry
        for prev in r.iloc[3].split("|"):
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


def _condition_mechanism_assertion(row: dict) -> dict:
    """Serialize one gene-condition mechanism without collapsing its identity.

    G2P and Orphadata differ in their native fields, but after source parsing
    they share one contract. Keeping the source record, source disease ID, and
    MONDO ID together lets the runtime hub join only the matching HPO disease
    assertion. ``allelic_requirement`` is the source's requirement; HPO-derived
    inheritance and penetrance are intentionally attached later and remain
    separately attributable.
    """
    return _nonempty(
        {
            "source_record_id": row.get("source_record_id", ""),
            "source_condition_id": row.get("source_condition_id", ""),
            "mondo_id": row.get("mondo_id", ""),
            "disease": row.get("disease", ""),
            "mechanism": row.get("mechanism", ""),
            "mechanism_raw": row.get("patho_mode_raw", ""),
            "allelic_requirement": row.get("inheritance", ""),
            "disease_scope": row.get("disease_scope", ""),
            "priva_scope": row.get("priva_scope", ""),
            "scope_review_status": row.get("scope_review_status", ""),
            "confidence": row.get("confidence", ""),
            "disease_confidence": row.get("disease_confidence", ""),
            "pmids": _split_pmids(row.get("pmids", "")),
        }
    )


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
    "G2P_DDG2P": ("gene_level", _condition_mechanism_assertion),
    "Orphadata": ("gene_level", _condition_mechanism_assertion),
    "PanelApp": ("gene_level", _panelapp_assertion),
    "ClinGen_Dosage": ("gene_level", _clingen_assertion),
    "GoFCards": ("variant_level", _gofcards_assertion),
}


def build_unified_json(
    unified_df: pd.DataFrame,
    hgnc_map: dict[str, dict[str, str]],
) -> dict[str, Any]:
    genes: dict[str, dict[str, Any]] = {}
    unmapped: list[str] = []

    for row in unified_df.to_dict(orient="records"):
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
            assertion = _nonempty(assertion)

        mech = row.get("mechanism", "")
        for m in mech.split("|"):
            if m.strip():
                gene["mechanisms"].add(m.strip())

        gene[level_key].append({source: assertion})

    result: dict[str, Any] = {}
    result["_meta"] = {
        "version": "1.0",
        "built_at": iso_now(),
        "total_genes": len(genes),
        "unmapped_symbols": sorted(set(unmapped)) if unmapped else [],
        "sources": {
            "G2P_DDG2P": {
                "level": "gene_level",
                "keys": ["source_record_id", "source_condition_id", "mondo_id", "disease", "mechanism", "mechanism_raw", "allelic_requirement", "disease_scope", "priva_scope", "scope_review_status", "confidence", "disease_confidence", "pmids"],
                "mechanism_raw": "original G2P molecular mechanism text",
                "confidence_values": "definitive | strong | moderate | limited | conflicting_or_refuted",
            },
            "Orphadata": {
                "level": "gene_level",
                "keys": ["source_record_id", "source_condition_id", "mondo_id", "disease", "mechanism", "mechanism_raw", "allelic_requirement", "disease_scope", "priva_scope", "scope_review_status", "confidence", "disease_confidence", "pmids"],
                "scope": "assessed disease-causing germline relationships with explicit GOF or LOF mechanism",
            },
            "GoFCards": {
                "level": "variant_level",
                "keys": ["mechanism", "disease", "chr", "pos", "ref", "alt", "transcript", "consequence", "pscore", "pmids", "function", "pathway", "animal_model", "cell_model"],
                "assembly": "hg19",
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
        },
        "mechanism_values": ["LOF", "GOF", "DOMINANT_NEGATIVE", "TRIPLOSENSITIVITY", "PANELAPP_GREEN_NON_LOF_PATHO_HISTORY"],
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument(
        "--output-tsv",
        type=Path,
        default=None,
        help=(
            "Canonical condition-mechanism evidence TSV. Defaults to "
            "<cache-dir>/gene_pathogenic_mechanism_evidence.tsv."
        ),
    )
    parser.add_argument(
        "--disease-scope-registry",
        type=Path,
        default=DEFAULT_DISEASE_SCOPE_REGISTRY,
        help=(
            "Pinned OMIM/ORPHA-to-MONDO disease-scope table used to enrich "
            "condition-specific mechanism assertions."
        ),
    )
    parser.add_argument(
        "--shared-raw-dir",
        default=str(DEFAULT_SHARED_RAW_DIR),
        help=(
            "Shared raw-data directory used for downloads and parsing. Audit "
            "outputs are written flat into --cache-dir. "
            "Expected filenames are exactly: AllG2P.csv, "
            "ClinGen_gene_curation_list_GRCh38.tsv, "
            "gofcards/gofcards_data_download.xlsx, and panelapp_all_panels.json. "
            "Use --shared-raw-dir '' to store raw "
            "files under --cache-dir/raw instead."
        ),
    )
    parser.add_argument("--force", action="store_true", help="Check/download all sources.")
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
            "--proxy-url is a SOCKS proxy, where curl is required."
        ),
    )
    parser.add_argument("--stale-lock-hours", type=float, default=24.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cache_dir = args.cache_dir.resolve()
    raw_dir = (
        Path(args.shared_raw_dir).expanduser().resolve()
        if args.shared_raw_dir.strip()
        else cache_dir / "raw"
    )
    raw_dir_is_shared = raw_dir != cache_dir / "raw"
    # Every output lands directly in the cache directory; the former
    # prepared/ and metadata/ split only added a folder level.
    manifest_path = cache_dir / "source_manifest.json"
    lock_path = raw_dir / ".build.lock" if raw_dir_is_shared else cache_dir / ".build.lock"

    with BuildLock(lock_path, stale_hours=args.stale_lock_hours):
        previous_manifest = read_json(manifest_path)
        previous_sources = previous_manifest.get("sources", {})
        source_meta: dict[str, Any] = {}

        for source in CONDITION_EVIDENCE_SPECS:
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

        manifest = {
            "schema_version": "1.0",
            "cache_dir": str(cache_dir),
            "raw_dir": str(raw_dir),
            "raw_dir_mode": "shared" if raw_dir_is_shared else "cache_local",
            "lock_path": str(lock_path),
            "proxy_mode": "explicit_proxy_url" if args.proxy_url else "environment_or_none",
            "download_tool": args.download_tool,
            "expected_raw_filenames": [
                source.raw_filename for source in CONDITION_EVIDENCE_SPECS
            ],
            "built_at_utc": iso_now(),
            "sources": source_meta,
        }
        all_evidence_df = parse_all_sources(
            raw_dir=raw_dir,
            disease_scope_registry=args.disease_scope_registry,
        )
        evidence_path = (
            args.output_tsv.expanduser().resolve()
            if args.output_tsv is not None
            else cache_dir / "gene_pathogenic_mechanism_evidence.tsv"
        )
        evidence_count = write_gene_pathogenic_mechanism_evidence(
            evidence_path,
            all_evidence_df,
        )
        print(f"evidence_output: {evidence_count} rows to {evidence_path}", file=sys.stderr)

        # Do not advance successful-build provenance until the required output
        # has been generated. A failed run must remain visibly retryable.
        write_json(manifest_path, manifest)
        manifest_tsv_count = write_source_manifest_tsv(
            cache_dir / "source_manifest.tsv",
            manifest,
        )

        run_summary = {
            "built_at_utc": manifest["built_at_utc"],
            "cache_dir": str(cache_dir),
            "evidence_rows": evidence_count,
            "source_manifest_rows": manifest_tsv_count,
            "outputs": {
                "source_manifest_json": str(manifest_path),
                "source_manifest_tsv": str(cache_dir / "source_manifest.tsv"),
                "gene_pathogenic_mechanism_evidence": str(evidence_path),
            },
        }
        write_json(cache_dir / "run_summary.json", run_summary)
        print(json.dumps(run_summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
