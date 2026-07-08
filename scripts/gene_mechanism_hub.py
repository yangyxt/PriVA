#!/usr/bin/env python3
"""Gene-level mechanism and inheritance hub for PriVA.

This module normalizes a gene query to one current HGNC symbol, then reports:

1. Curated mechanism history from the existing mechanism JSON cache and the
   companion DDG2P/G2P evidence table.
2. Known inheritance/HI status using PriVA's existing inheritance decision
   function from ``acmg_criteria_assign.py``.

The mechanism history is gene/history context only. It must not be treated as a
variant-level GoF/DN assertion unless a separate variant-specific matcher proves
that the query variant corresponds to a curated functional variant.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PRIVA_ROOT = SCRIPT_DIR.parent
DATA_DIR = PRIVA_ROOT / "data"

DEFAULT_HGNC_TABLE = DATA_DIR / "hgnc" / "non_alt_loci_set.tsv"
DEFAULT_HPO_COLLAPSED = DATA_DIR / "hpo" / "genes_to_phenotype.collapse.tsv.gz"
DEFAULT_CLINGEN_DOSAGE = DATA_DIR / "clingen" / "gene_dosage_sensitivity.hg19.tsv"
DEFAULT_LOEUF_TABLE = DATA_DIR / "loeuf" / "loeuf_dataset.tsv.gz"
DEFAULT_MECHANISM_JSON = (
    DATA_DIR / "gene_pathogenic_mechanism" / "prepared" / "gene_mechanism_curated_assertions.json"
)
DEFAULT_DDG2P_MECHANISM_EVIDENCE = (
    DATA_DIR / "gene_pathogenic_mechanism" / "prepared" / "gene_pathogenic_mechanism_evidence.tsv"
)
DEFAULT_GOFCARDS_EXACT_GOF_HGVSP = DATA_DIR / "gofcards" / "gofcards_exact_gof_hgvsp.tsv.gz"
DEFAULT_GOFCARDS_STEP1_TSV = DATA_DIR / "gofcards" / "legacy" / "gofcards_gof.hg19.numb.anno.step1.tsv"
DEFAULT_GOFCARDS_ACTIVE_TSV = DATA_DIR / "gofcards" / "legacy" / "active_json_gof_dn_variant_level.tsv"
DEFAULT_GOFCARDS_RAW_XLSX = (
    DATA_DIR / "gene_pathogenic_mechanism" / "raw" / "gofcards" / "gofcards_data_download.xlsx"
)
DEFAULT_GOFCARDS_CONVERSION_AUDIT_TSV = (
    DATA_DIR / "gofcards" / "legacy" / "gofcards_gof.hg19.numb.vcf_conversion_audit.tsv"
)

LOOKUP_FIELD_PRIORITY = (
    "symbol",
    "hgnc_id",
    "ensembl_gene_id",
    "entrez_id",
    "prev_symbol",
    "alias_symbol",
    "refseq_accession",
    "uniprot_ids",
    "mane_select",
)

CANONICAL_MECHANISMS = {"GOF", "DOMINANT_NEGATIVE", "TRIPLOSENSITIVITY"}
DDG2P_LOF_RAW_VALUES = {
    "loss of function",
    "absent gene product",
    "decreased gene product level",
    "reduced gene product level",
}
DDG2P_USABLE_MECHANISM_CONFIDENCE = {"high", "moderate"}
DDG2P_USABLE_DISEASE_CONFIDENCE = {"definitive", "strong", "moderate"}
DDG2P_RECESSIVE_LOF_INHERITANCE = {
    "biallelic_autosomal",
    "monoallelic_X",
    "monoallelic_X_hemizygous",
}
DDG2P_DOMINANT_LOF_INHERITANCE = {
    "monoallelic_autosomal",
    "monoallelic_X_heterozygous",
    "monoallelic_Y_hemizygous",
}
GENE_MECHANISM_CATEGORY_ORDER = (
    "LoF_recessive",
    "dominant_ambiguous",
    "dominant_HI",
    "dominant_GOF",
    "dominant_DN",
)
logger = logging.getLogger(__name__)


def _clean(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "<na>"}:
        return ""
    return text


def _norm(value: Any) -> str:
    return _clean(value).upper()


def _norm_chrom(value: Any) -> str:
    text = _clean(value)
    if not text:
        return ""
    return text if text.startswith("chr") else f"chr{text}"


def _norm_chrom_key(value: Any) -> str:
    text = _clean(value)
    if not text:
        return ""
    if text.lower().startswith("chr"):
        text = text[3:]
    if text.upper() in {"M", "MT"}:
        return "MT"
    return text.upper()


def _norm_pos_key(value: Any) -> str:
    text = _clean(value)
    if not text:
        return ""
    try:
        as_float = float(text)
    except ValueError:
        return text
    if as_float.is_integer():
        return str(int(as_float))
    return text


def _norm_allele_key(value: Any) -> str:
    return _clean(value).upper()


def _genomic_match_key(chrom: Any, pos: Any, ref: Any, alt: Any) -> str:
    chrom_key = _norm_chrom_key(chrom)
    pos_key = _norm_pos_key(pos)
    if not chrom_key or not pos_key:
        return ""
    return "|".join(
        [
            chrom_key,
            pos_key,
            _norm_allele_key(ref),
            _norm_allele_key(alt),
        ]
    )


def _genomic_match_key_from_text(value: Any) -> str:
    text = _clean(value)
    if not text:
        return ""
    parts = text.split("|")
    if len(parts) != 4:
        return ""
    return _genomic_match_key(parts[0], parts[1], parts[2], parts[3])


def _norm_hgvsp(value: Any) -> str:
    text = _clean(value)
    if not text:
        return ""
    protein = text.split(":")[-1].strip()
    if protein.startswith("p."):
        protein = protein[2:]
    return protein.upper()


def _split_multi(value: Any) -> list[str]:
    text = _clean(value)
    if not text:
        return []
    return [part.strip() for part in text.split("|") if part.strip()]


def _split_pmids(value: Any) -> list[str]:
    text = _clean(value)
    if not text:
        return []
    return [part.strip() for part in re.split(r"[;|,]", text) if part.strip()]


def _safe_int(value: Any) -> int | None:
    text = _clean(value)
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _safe_float(value: Any) -> float:
    text = _clean(value)
    if not text:
        return np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


def classify_gene_mechanism_categories(
    *,
    recessive: bool,
    dominant: bool,
    haplo_insufficient: bool,
    has_gof_history: bool = False,
    has_dn_history: bool = False,
    has_triplosensitivity_history: bool = False,
) -> list[str]:
    """Classify gene-level pathogenic mechanism categories.

    Inputs are the already-computed PriVA inheritance outputs plus curated
    gene/history mechanism flags. This function does not infer inheritance and
    does not assign a query variant as GoF/DN.
    """
    categories: list[str] = []
    if recessive:
        categories.append("LoF_recessive")
    if dominant:
        dominant_categories: list[str] = []
        if haplo_insufficient:
            dominant_categories.append("dominant_HI")
        if has_gof_history or has_triplosensitivity_history:
            dominant_categories.append("dominant_GOF")
        if has_dn_history:
            dominant_categories.append("dominant_DN")
        if dominant_categories:
            categories.extend(dominant_categories)
        else:
            categories.append("dominant_ambiguous")
    return categories


def _ordered_gene_mechanism_categories(categories: set[str]) -> list[str]:
    """Return compact mechanism categories in stable downstream order."""
    return [category for category in GENE_MECHANISM_CATEGORY_ORDER if category in categories]


def _hpo_inheritance_flags(inheritance_modes: Any) -> dict[str, bool]:
    """Parse the HPO inheritance mode text used by the collapsed gene table."""
    modes = {_clean(part) for part in _clean(inheritance_modes).split(";") if _clean(part)}
    return {
        "recessive": bool(
            {
                "Autosomal recessive inheritance",
                "X-linked recessive inheritance",
                "X-linked inheritance",
            }
            & modes
        ),
        "dominant": bool(
            {
                "Autosomal dominant inheritance",
                "X-linked dominant inheritance",
                "Y-linked inheritance",
            }
            & modes
        ),
    }


class _LocalHgncResolver:
    """Small fallback resolver using PriVA's local HGNC table."""

    def __init__(self, hgnc_table: Path = DEFAULT_HGNC_TABLE) -> None:
        self.hgnc_table = Path(hgnc_table)
        self.rows: list[dict[str, str]] = []
        self.lookup: dict[str, list[tuple[int, dict[str, str]]]] = defaultdict(list)
        self._load()

    def _load(self) -> None:
        with open(self.hgnc_table, encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.rows = [dict(row) for row in reader]
        for row in self.rows:
            for priority, field in enumerate(LOOKUP_FIELD_PRIORITY):
                for value in _split_multi(row.get(field, "")):
                    key = _norm(value)
                    if key:
                        self.lookup[key].append((priority, row))

    def resolve(self, query: Any, *, warn_ambiguous: bool = True) -> str:
        cleaned = _clean(query)
        if not cleaned:
            return ""
        matches = self.lookup.get(_norm(cleaned), [])
        if not matches:
            return cleaned
        matches = sorted(matches, key=lambda item: (item[0], item[1].get("symbol", "")))
        best_priority = matches[0][0]
        best = [row for priority, row in matches if priority == best_priority]
        symbols = sorted({row.get("symbol", "") for row in best if row.get("symbol", "")})
        if len(symbols) > 1:
            if warn_ambiguous:
                logger.warning(
                    "Ambiguous HGNC query %s resolves to %s; keeping input symbol",
                    cleaned,
                    ",".join(symbols),
                )
            return cleaned
        return symbols[0] if symbols else cleaned


class GeneMechanismHub:
    """Normalize symbols, then summarize gene mechanism and inheritance history."""

    def __init__(
        self,
        *,
        mechanism_json: str | Path = DEFAULT_MECHANISM_JSON,
        ddg2p_evidence: str | Path = DEFAULT_DDG2P_MECHANISM_EVIDENCE,
        hpo_collapsed: str | Path = DEFAULT_HPO_COLLAPSED,
        clingen_dosage: str | Path = DEFAULT_CLINGEN_DOSAGE,
        loeuf_table: str | Path = DEFAULT_LOEUF_TABLE,
        hgnc_table: str | Path = DEFAULT_HGNC_TABLE,
        use_hgnc_package: bool = True,
    ) -> None:
        self.mechanism_json = Path(mechanism_json)
        self.ddg2p_evidence = Path(ddg2p_evidence) if ddg2p_evidence else Path("")
        self.hpo_collapsed = Path(hpo_collapsed)
        self.clingen_dosage = Path(clingen_dosage)
        self.loeuf_table = Path(loeuf_table)
        self.hgnc_table = Path(hgnc_table)

        self._resolver = self._build_resolver(use_hgnc_package)
        self._mechanism_by_symbol: dict[str, dict[str, Any]] | None = None
        self._ddg2p_lof_by_symbol: dict[str, list[dict[str, Any]]] | None = None
        self._hpo_by_symbol: dict[str, dict[str, str]] | None = None
        self._clingen_by_symbol: dict[str, dict[str, Any]] | None = None
        self._loeuf_by_symbol: dict[str, float] | None = None
        self._gofcards_by_symbol_hgvsp: dict[tuple[str, str], list[dict[str, str]]] | None = None
        self._gofcards_by_symbol_genomic: dict[tuple[str, str, str, str], list[dict[str, str]]] | None = None

    def _build_resolver(self, use_hgnc_package: bool) -> Any:
        if use_hgnc_package:
            try:
                from hgnc_symbol_resolver import HgncResolver

                return HgncResolver()
            except Exception:
                pass
        return _LocalHgncResolver(self.hgnc_table)

    def resolve_symbol(self, gene_symbol: Any) -> str:
        """Return one current HGNC symbol for a query string."""
        if hasattr(self._resolver, "resolve"):
            return self._resolver.resolve(gene_symbol)
        raise TypeError("HGNC resolver object does not expose resolve()")

    def _try_resolve_symbol(self, gene_symbol: Any) -> str:
        """Resolve cache-side symbols without failing on ambiguous legacy aliases."""
        try:
            if isinstance(self._resolver, _LocalHgncResolver):
                return self._resolver.resolve(gene_symbol, warn_ambiguous=False)
            return self.resolve_symbol(gene_symbol)
        except ValueError:
            return _clean(gene_symbol)

    def _resolved_symbol_key(self, gene_symbol: Any) -> str:
        return self.resolve_symbol(gene_symbol)

    def _load_mechanisms(self) -> dict[str, dict[str, Any]]:
        if self._mechanism_by_symbol is not None:
            return self._mechanism_by_symbol
        with open(self.mechanism_json, encoding="utf-8") as handle:
            raw = json.load(handle)
        by_symbol: dict[str, dict[str, Any]] = {}
        for info in raw.values():
            symbol = _clean(info.get("symbol"))
            if not symbol:
                continue
            by_symbol[symbol] = info
            resolved = self._try_resolve_symbol(symbol)
            if resolved and resolved not in by_symbol:
                by_symbol[resolved] = info
        self._mechanism_by_symbol = by_symbol
        return by_symbol

    @staticmethod
    def _is_strict_ddg2p_lof(row: pd.Series) -> bool:
        raw = _clean(row.get("patho_mode_raw")).lower()
        if raw.startswith("undetermined") or "non-loss" in raw:
            return False
        norm_tokens = {
            token.strip().upper()
            for token in re.split(r"[;,]", _clean(row.get("normalized_mechanisms")))
            if token.strip()
        }
        return raw in DDG2P_LOF_RAW_VALUES or norm_tokens == {"LOF"}

    def _load_ddg2p_lof(self) -> dict[str, list[dict[str, Any]]]:
        if self._ddg2p_lof_by_symbol is not None:
            return self._ddg2p_lof_by_symbol
        by_symbol: dict[str, list[dict[str, Any]]] = defaultdict(list)
        if not self.ddg2p_evidence or not self.ddg2p_evidence.exists():
            self._ddg2p_lof_by_symbol = {}
            return self._ddg2p_lof_by_symbol

        try:
            df = pd.read_csv(
                self.ddg2p_evidence,
                sep="\t",
                dtype=str,
                low_memory=False,
            ).fillna("")
        except Exception as exc:
            logger.warning("Failed to load DDG2P/G2P evidence table %s: %s", self.ddg2p_evidence, exc)
            self._ddg2p_lof_by_symbol = {}
            return self._ddg2p_lof_by_symbol

        if "source" in df.columns:
            df = df.loc[df["source"].eq("G2P_DDG2P")].copy()
        if df.empty or "gene_symbol" not in df.columns:
            self._ddg2p_lof_by_symbol = {}
            return self._ddg2p_lof_by_symbol

        for _, row in df.iterrows():
            if not self._is_strict_ddg2p_lof(row):
                continue
            mechanism_confidence = _clean(row.get("mechanism_confidence")).lower()
            disease_confidence = _clean(row.get("disease_confidence")).lower()
            if mechanism_confidence not in DDG2P_USABLE_MECHANISM_CONFIDENCE:
                continue
            if disease_confidence not in DDG2P_USABLE_DISEASE_CONFIDENCE:
                continue

            symbol = _clean(row.get("gene_symbol"))
            if not symbol:
                continue
            inheritance = _clean(row.get("inheritance"))
            record = {
                "level": "gene_level",
                "source": "G2P_DDG2P",
                "mechanism": "LOF",
                "pmids": _split_pmids(row.get("pmids", "")),
                "disease": _clean(row.get("disease_label")),
                "inheritance": inheritance,
                "mechanism_confidence": mechanism_confidence,
                "disease_confidence": disease_confidence,
                "source_record_id": _clean(row.get("source_record_id")),
            }
            resolved = self._try_resolve_symbol(symbol)
            keys = {symbol}
            if resolved:
                keys.add(resolved)
            for key in keys:
                by_symbol[key].append(record)

        self._ddg2p_lof_by_symbol = {
            symbol: records for symbol, records in by_symbol.items() if records
        }
        return self._ddg2p_lof_by_symbol

    def _load_hpo(self) -> dict[str, dict[str, str]]:
        if self._hpo_by_symbol is not None:
            return self._hpo_by_symbol
        df = pd.read_csv(self.hpo_collapsed, sep="\t", dtype=str, low_memory=False).fillna("")
        by_symbol: dict[str, dict[str, str]] = {}
        for _, row in df.iterrows():
            symbol = _clean(row.get("gene_symbol"))
            if not symbol:
                continue
            record = {
                "ncbi_gene_id": _clean(row.get("ncbi_gene_id")),
                "hpo_id": _clean(row.get("hpo_id")),
                "inheritance_modes": _clean(row.get("inheritance_modes")),
                "inheritance_disease_ids": _clean(row.get("inheritance_disease_ids")),
            }
            by_symbol.setdefault(symbol, record)
            resolved = self._try_resolve_symbol(symbol)
            if resolved:
                by_symbol.setdefault(resolved, record)
        self._hpo_by_symbol = by_symbol
        return by_symbol

    def _load_clingen(self) -> dict[str, dict[str, Any]]:
        if self._clingen_by_symbol is not None:
            return self._clingen_by_symbol
        df = pd.read_csv(self.clingen_dosage, sep="\t", dtype=str, low_memory=False).fillna("")
        by_symbol: dict[str, dict[str, Any]] = {}
        for _, row in df.iterrows():
            symbol = _clean(row.get("#Gene Symbol"))
            if not symbol:
                continue
            record = {
                "haploinsufficiency_score": _safe_int(row.get("Haploinsufficiency Score")),
                "haploinsufficiency_description": _clean(
                    row.get("Haploinsufficiency Description")
                ),
                "triplosensitivity_score": _safe_int(row.get("Triplosensitivity Score")),
                "triplosensitivity_description": _clean(
                    row.get("Triplosensitivity Description")
                ),
            }
            by_symbol.setdefault(symbol, record)
            resolved = self._try_resolve_symbol(symbol)
            if resolved:
                by_symbol.setdefault(resolved, record)
        self._clingen_by_symbol = by_symbol
        return by_symbol

    def _load_loeuf(self) -> dict[str, float]:
        if self._loeuf_by_symbol is not None:
            return self._loeuf_by_symbol
        df = pd.read_csv(
            self.loeuf_table,
            sep="\t",
            dtype=str,
            low_memory=False,
            usecols=["#gene", "canonical", "oe_lof_upper"],
        ).fillna("")
        df["__loeuf"] = df["oe_lof_upper"].map(_safe_float)
        out: dict[str, float] = {}
        canonical = df[df["canonical"].str.lower().eq("true")]
        for symbol, grp in canonical.groupby("#gene"):
            values = [v for v in grp["__loeuf"] if not math.isnan(v)]
            if values:
                out.setdefault(symbol, min(values))
        for symbol, grp in df.groupby("#gene"):
            if symbol in out:
                continue
            values = [v for v in grp["__loeuf"] if not math.isnan(v)]
            if values:
                out[symbol] = min(values)
        resolved_out = dict(out)
        for symbol, value in out.items():
            resolved = self._try_resolve_symbol(symbol)
            if resolved:
                resolved_out.setdefault(resolved, value)
        self._loeuf_by_symbol = resolved_out
        return resolved_out

    def _load_gofcards_variant_hgvsp(
        self,
        *,
        gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
        gofcards_step1_tsv: str | Path = DEFAULT_GOFCARDS_STEP1_TSV,
        gofcards_active_tsv: str | Path = DEFAULT_GOFCARDS_ACTIVE_TSV,
        gofcards_raw_xlsx: str | Path = DEFAULT_GOFCARDS_RAW_XLSX,
        gofcards_conversion_audit_tsv: str | Path = DEFAULT_GOFCARDS_CONVERSION_AUDIT_TSV,
    ) -> dict[tuple[str, str], list[dict[str, str]]]:
        """Return exact-match GoFCards GOF protein-change evidence.

        The key is normalized ``(HGNC symbol, protein change)``. This lookup is
        only for variant-level evidence; it must not be used to convert all
        variants in a gene into GOF assertions.
        """
        if self._gofcards_by_symbol_hgvsp is not None:
            return self._gofcards_by_symbol_hgvsp

        compact_path = Path(gofcards_exact_hgvsp_tsv)
        if compact_path.exists():
            compact = pd.read_csv(
                compact_path,
                sep="\t",
                dtype=str,
                low_memory=False,
            ).fillna("")
            by_symbol_hgvsp: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
            for _, row in compact.iterrows():
                cache_symbol = (
                    row.get("HGNC_Symbol")
                    or row.get("symbol")
                    or row.get("match_symbol")
                    or row.get("gofcards_symbol_resolved")
                    or row.get("gofcards_symbol")
                )
                symbol = self._try_resolve_symbol(cache_symbol)
                if not symbol:
                    symbol = _clean(cache_symbol)
                hgvsp_key = _clean(row.get("hgvsp_key")) or _norm_hgvsp(row.get("HGVSp"))
                if not symbol or not hgvsp_key:
                    continue
                match_status = _clean(row.get("match_status")) or _clean(
                    row.get("gofcards_hgvs_match_status")
                )
                by_symbol_hgvsp[(symbol, hgvsp_key)].append(
                    {
                        "source": _clean(row.get("source")),
                        "mechanism": _clean(row.get("mechanism")),
                        "build": _clean(row.get("build")),
                        "gofcards_variant_id": _clean(row.get("gofcards_variant_id")),
                        "gofcards_accession_id": _clean(row.get("gofcards_accession_id")),
                        "disease": _clean(row.get("disease")),
                        "pmids": _clean(row.get("pmids")),
                        "pscore": _clean(row.get("pscore")),
                        "function": _clean(row.get("function")),
                        "pathway": _clean(row.get("pathway")),
                        "transcript": _clean(row.get("GoFCards_transcript"))
                        or _clean(row.get("transcript")),
                        "symbol": symbol,
                        "match_symbol": _clean(row.get("HGNC_Symbol"))
                        or _clean(row.get("match_symbol")),
                        "gofcards_symbol": _clean(row.get("gofcards_symbol")),
                        "vep_symbol": _clean(row.get("vep_symbol")),
                        "gofcards_symbol_resolved": _clean(row.get("gofcards_symbol_resolved")),
                        "vep_symbol_resolved": _clean(row.get("vep_symbol_resolved")),
                        "source_refseq_transcript": _clean(row.get("GoFCards_transcript"))
                        or _clean(row.get("source_refseq_transcript")),
                        "vep_transcript": _clean(row.get("VEP_transcript"))
                        or _clean(row.get("vep_transcript")),
                        "chrom": _norm_chrom(row.get("hg19_chrom") or row.get("chrom")),
                        "pos": _clean(row.get("hg19_pos")) or _clean(row.get("pos")),
                        "ref": _clean(row.get("hg19_ref")) or _clean(row.get("ref")),
                        "alt": _clean(row.get("hg19_alt")) or _clean(row.get("alt")),
                        "hg19_genomic_key": _clean(row.get("hg19_genomic_key")),
                        "hg19_vcf_key": _clean(row.get("hg19_vcf_key")),
                        "hg38_genomic_key": _clean(row.get("hg38_genomic_key")),
                        "hg38_vcf_key": _clean(row.get("hg38_vcf_key")),
                        "match_key_types": _clean(row.get("match_key_types")),
                        "hg19_chrom": _clean(row.get("hg19_chrom")),
                        "hg19_start": _clean(row.get("hg19_pos"))
                        or _clean(row.get("hg19_start")),
                        "hg19_end": _clean(row.get("hg19_pos")) or _clean(row.get("hg19_end")),
                        "hg19_ref": _clean(row.get("hg19_ref")),
                        "hg19_alt": _clean(row.get("hg19_alt")),
                        "hg19_vcf_pos": _clean(row.get("hg19_pos"))
                        or _clean(row.get("hg19_vcf_pos")),
                        "hg19_vcf_ref": _clean(row.get("hg19_ref"))
                        or _clean(row.get("hg19_vcf_ref")),
                        "hg19_vcf_alt": _clean(row.get("hg19_alt"))
                        or _clean(row.get("hg19_vcf_alt")),
                        "hg19_vcf_status": _clean(row.get("hg19_vcf_status")),
                        "hg38_chrom": _clean(row.get("hg38_chrom")),
                        "hg38_start": _clean(row.get("hg38_pos"))
                        or _clean(row.get("hg38_start")),
                        "hg38_end": _clean(row.get("hg38_pos")) or _clean(row.get("hg38_end")),
                        "hg38_ref": _clean(row.get("hg38_ref")),
                        "hg38_alt": _clean(row.get("hg38_alt")),
                        "hg38_vcf_pos": _clean(row.get("hg38_pos"))
                        or _clean(row.get("hg38_vcf_pos")),
                        "hg38_vcf_ref": _clean(row.get("hg38_ref"))
                        or _clean(row.get("hg38_vcf_ref")),
                        "hg38_vcf_alt": _clean(row.get("hg38_alt"))
                        or _clean(row.get("hg38_vcf_alt")),
                        "hg38_refalt_status": _clean(row.get("hg38_refalt_status")),
                        "HGVSc": _clean(row.get("HGVSc")),
                        "HGVSp": _clean(row.get("HGVSp")),
                        "normalized_hgvsp": _clean(row.get("normalized_hgvsp"))
                        or _norm_hgvsp(row.get("HGVSp")),
                        "gofcards_hgvsc_key": _clean(row.get("gofcards_hgvsc_key")),
                        "gofcards_hgvsp_key": _clean(row.get("gofcards_hgvsp_key")),
                        "vep_hgvsc_key": _clean(row.get("vep_hgvsc_key")),
                        "vep_hgvsp_key": _clean(row.get("vep_hgvsp_key")),
                        "gofcards_gene_match": _clean(row.get("gofcards_gene_match")),
                        "gofcards_hgvsc_match": _clean(row.get("gofcards_hgvsc_match")),
                        "gofcards_hgvsp_match": _clean(row.get("gofcards_hgvsp_match")),
                        "gofcards_hgvs_match_status": match_status,
                        "match_status": match_status,
                        "canonical_transcript": _clean(row.get("canonical_transcript")),
                        "gofcards_AAChange_refGene": _clean(row.get("raw_GoFCards_HGVS"))
                        or _clean(row.get("gofcards_AAChange_refGene")),
                        "derived_on": _clean(row.get("derived_on")),
                    }
                )
            self._gofcards_by_symbol_hgvsp = dict(by_symbol_hgvsp)
            return self._gofcards_by_symbol_hgvsp

        step1 = pd.read_csv(
            gofcards_step1_tsv,
            sep="\t",
            dtype=str,
            low_memory=False,
            usecols=["chrom", "pos", "ref", "alt", "SYMBOL", "HGVSc", "HGVSp", "CANONICAL"],
        ).fillna("")
        active = pd.read_csv(
            gofcards_active_tsv,
            sep="\t",
            dtype=str,
            low_memory=False,
        ).fillna("")
        active = active[
            (active["source"].eq("GoFCards"))
            & (active["mechanism"].eq("GOF"))
        ].copy()
        active["_row_index"] = [str(i) for i in range(1, len(active) + 1)]
        active["_gofcards_variant_id"] = active["_row_index"].map(
            lambda value: f"GOFCARDS_GOF_{int(value):04d}"
        )

        raw = pd.read_excel(
            gofcards_raw_xlsx,
            sheet_name="total3161",
            engine="openpyxl",
            dtype=str,
        ).fillna("")
        raw["_coord_key"] = list(
            zip(
                raw["chr"].map(_norm_chrom),
                raw["hg19start"].map(_clean),
                raw["ref"].map(_clean),
                raw["alt"].map(_clean),
            )
        )
        accessions_by_coord: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
        for _, row in raw.iterrows():
            accession = _clean(row.get("Order numbe"))
            if accession:
                accessions_by_coord[row["_coord_key"]].add(accession)

        audit = pd.read_csv(
            gofcards_conversion_audit_tsv,
            sep="\t",
            dtype=str,
            low_memory=False,
            usecols=[
                "row_index",
                "status",
                "chrom",
                "vcf_chrom",
                "vcf_pos",
                "vcf_ref",
                "vcf_alt",
                "source_pos",
                "source_ref",
                "source_alt",
            ],
        ).fillna("")
        audit = audit[audit["status"].eq("written")].copy()
        audit["_row_index"] = audit["row_index"].map(_clean)
        audit["_vcf_coord_key"] = list(
            zip(
                audit["vcf_chrom"].map(_norm_chrom),
                audit["vcf_pos"].map(_clean),
                audit["vcf_ref"].map(_clean),
                audit["vcf_alt"].map(_clean),
            )
        )
        audit["_source_coord_key"] = list(
            zip(
                audit["chrom"].map(_norm_chrom),
                audit["source_pos"].map(_clean),
                audit["source_ref"].map(_clean),
                audit["source_alt"].map(_clean),
            )
        )
        audit_by_row = {
            row["_row_index"]: {
                "vcf_coord_key": row["_vcf_coord_key"],
                "source_coord_key": row["_source_coord_key"],
            }
            for _, row in audit.iterrows()
        }

        active_by_coord: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
        for _, row in active.iterrows():
            row_index = _clean(row.get("_row_index"))
            source_coord_key = (
                _norm_chrom(row.get("chr")),
                _clean(row.get("pos")),
                _clean(row.get("ref")),
                _clean(row.get("alt")),
            )
            audit_row = audit_by_row.get(row_index, {})
            vcf_coord_key = audit_row.get("vcf_coord_key", source_coord_key)
            source_coord_key = audit_row.get("source_coord_key", source_coord_key)
            accessions = sorted(accessions_by_coord.get(source_coord_key, set()))
            active_by_coord[vcf_coord_key].append(
                {
                    "gofcards_variant_id": _clean(row.get("_gofcards_variant_id")),
                    "gofcards_accession_id": ";".join(accessions),
                    "disease": _clean(row.get("disease")),
                    "pmids": _clean(row.get("pmids")),
                    "pscore": _clean(row.get("pscore")),
                    "function": _clean(row.get("function")),
                    "pathway": _clean(row.get("pathway")),
                    "transcript": _clean(row.get("transcript")),
                }
            )

        by_symbol_hgvsp: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for _, row in step1.iterrows():
            hgvsp_key = _norm_hgvsp(row.get("HGVSp"))
            if not hgvsp_key:
                continue
            symbol = self._try_resolve_symbol(row.get("SYMBOL"))
            if not symbol:
                symbol = _clean(row.get("SYMBOL"))
            coord_key = (
                _norm_chrom(row.get("chrom")),
                _clean(row.get("pos")),
                _clean(row.get("ref")),
                _clean(row.get("alt")),
            )
            for active_row in active_by_coord.get(coord_key, []):
                match = {
                    **active_row,
                    "symbol": symbol,
                    "chrom": coord_key[0],
                    "pos": coord_key[1],
                    "ref": coord_key[2],
                    "alt": coord_key[3],
                    "HGVSc": _clean(row.get("HGVSc")),
                    "HGVSp": _clean(row.get("HGVSp")),
                    "canonical_transcript": _clean(row.get("CANONICAL")),
                }
                by_symbol_hgvsp[(symbol, hgvsp_key)].append(match)

        for key, rows in by_symbol_hgvsp.items():
            seen: set[tuple[str, str, str, str, str, str]] = set()
            unique_rows: list[dict[str, str]] = []
            for row in sorted(
                rows,
                key=lambda item: (
                    item.get("canonical_transcript") != "YES",
                    item.get("gofcards_accession_id", ""),
                    item.get("gofcards_variant_id", ""),
                    item.get("HGVSp", ""),
                ),
            ):
                dedup_key = (
                    row.get("gofcards_variant_id", ""),
                    row.get("gofcards_accession_id", ""),
                    row.get("HGVSc", ""),
                    row.get("HGVSp", ""),
                    row.get("chrom", ""),
                    row.get("pos", ""),
                )
                if dedup_key in seen:
                    continue
                seen.add(dedup_key)
                unique_rows.append(row)
            by_symbol_hgvsp[key] = unique_rows

        self._gofcards_by_symbol_hgvsp = dict(by_symbol_hgvsp)
        return self._gofcards_by_symbol_hgvsp

    @staticmethod
    def _deduplicate_gofcards_matches(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        seen: set[tuple[str, str, str, str, str, str, str]] = set()
        unique_rows: list[dict[str, str]] = []
        for row in sorted(
            rows,
            key=lambda item: (
                item.get("canonical_transcript") != "YES",
                item.get("gofcards_accession_id", ""),
                item.get("gofcards_variant_id", ""),
                item.get("HGVSp", ""),
                item.get("hg38_vcf_key", ""),
                item.get("hg19_vcf_key", ""),
            ),
        ):
            dedup_key = (
                row.get("gofcards_variant_id", ""),
                row.get("gofcards_accession_id", ""),
                row.get("HGVSc", ""),
                row.get("HGVSp", ""),
                row.get("hg38_vcf_key", ""),
                row.get("hg19_vcf_key", ""),
                row.get("allele_key", ""),
            )
            if dedup_key in seen:
                continue
            seen.add(dedup_key)
            unique_rows.append(row)
        return unique_rows

    def _load_gofcards_variant_genomic(
        self,
        *,
        gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    ) -> dict[tuple[str, str, str, str], list[dict[str, str]]]:
        """Return exact-match GoFCards GOF genomic allele evidence.

        The key is ``(HGNC symbol, assembly, key_type, chrom|pos|ref|alt)``.
        ``key_type=vcf`` uses the VCF-padded allele representation intended for
        caller/VCF matching; ``key_type=genomic`` uses the sparse source genomic
        fields retained from GoFCards.
        """
        if self._gofcards_by_symbol_genomic is not None:
            return self._gofcards_by_symbol_genomic

        compact_path = Path(gofcards_exact_hgvsp_tsv)
        if not compact_path.exists():
            self._gofcards_by_symbol_genomic = {}
            return self._gofcards_by_symbol_genomic

        compact = pd.read_csv(
            compact_path,
            sep="\t",
            dtype=str,
            low_memory=False,
        ).fillna("")
        by_symbol_genomic: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
        key_columns = {
            ("hg19", "genomic"): "hg19_genomic_key",
            ("hg19", "vcf"): "hg19_vcf_key",
            ("hg38", "genomic"): "hg38_genomic_key",
            ("hg38", "vcf"): "hg38_vcf_key",
        }
        for _, row in compact.iterrows():
            cache_symbol = (
                row.get("HGNC_Symbol")
                or row.get("symbol")
                or row.get("match_symbol")
                or row.get("gofcards_symbol_resolved")
                or row.get("gofcards_symbol")
            )
            symbol = self._try_resolve_symbol(cache_symbol)
            if not symbol:
                symbol = _clean(cache_symbol)
            if not symbol:
                continue
            match = {str(col): _clean(row.get(col)) for col in compact.columns}
            match["symbol"] = symbol
            match.setdefault("match_status", _clean(row.get("gofcards_hgvs_match_status")))
            match.setdefault("gofcards_hgvs_match_status", _clean(row.get("match_status")))
            match["input_cache"] = str(compact_path)
            for (assembly, key_type), col in key_columns.items():
                key = _genomic_match_key_from_text(row.get(col))
                if not key:
                    continue
                by_symbol_genomic[(symbol, assembly, key_type, key)].append(match)

        self._gofcards_by_symbol_genomic = {
            key: self._deduplicate_gofcards_matches(rows)
            for key, rows in by_symbol_genomic.items()
        }
        return self._gofcards_by_symbol_genomic

    def match_gofcards_variant_gof(
        self,
        gene_symbol: Any,
        hgvsp: Any,
        *,
        gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
        gofcards_step1_tsv: str | Path = DEFAULT_GOFCARDS_STEP1_TSV,
        gofcards_active_tsv: str | Path = DEFAULT_GOFCARDS_ACTIVE_TSV,
        gofcards_raw_xlsx: str | Path = DEFAULT_GOFCARDS_RAW_XLSX,
        gofcards_conversion_audit_tsv: str | Path = DEFAULT_GOFCARDS_CONVERSION_AUDIT_TSV,
    ) -> dict[str, Any]:
        """Return variant-level GOF evidence for an exact GoFCards HGVSp match."""
        symbol = self._resolved_symbol_key(gene_symbol)
        hgvsp_key = _norm_hgvsp(hgvsp)
        if not symbol or not hgvsp_key:
            return {
                "input_symbol": _clean(gene_symbol),
                "symbol": symbol,
                "input_hgvsp": _clean(hgvsp),
                "variant_gof_tag": "",
                "gofcards_accession_id": "",
                "matches": [],
            }

        matches = self._load_gofcards_variant_hgvsp(
            gofcards_exact_hgvsp_tsv=gofcards_exact_hgvsp_tsv,
            gofcards_step1_tsv=gofcards_step1_tsv,
            gofcards_active_tsv=gofcards_active_tsv,
            gofcards_raw_xlsx=gofcards_raw_xlsx,
            gofcards_conversion_audit_tsv=gofcards_conversion_audit_tsv,
        ).get((symbol, hgvsp_key), [])
        accession_ids = sorted(
            {row.get("gofcards_accession_id", "") for row in matches if row.get("gofcards_accession_id", "")}
        )
        variant_ids = sorted(
            {row.get("gofcards_variant_id", "") for row in matches if row.get("gofcards_variant_id", "")}
        )
        return {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "input_hgvsp": _clean(hgvsp),
            "matched_hgvsp_key": hgvsp_key,
            "variant_gof_tag": "GOF" if matches else "",
            "gofcards_accession_id": ";".join(accession_ids),
            "gofcards_variant_id": ";".join(variant_ids),
            "matches": matches,
        }

    def match_gofcards_variant_gof_by_genomic(
        self,
        gene_symbol: Any,
        chrom: Any,
        pos: Any,
        ref: Any,
        alt: Any,
        *,
        assembly: str = "hg38",
        key_type: str = "auto",
        gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    ) -> dict[str, Any]:
        """Return variant-level GOF evidence for an exact GoFCards genomic match."""
        symbol = self._resolved_symbol_key(gene_symbol)
        normalized_key = _genomic_match_key(chrom, pos, ref, alt)
        assembly_key = _clean(assembly).lower()
        assembly_aliases = {
            "grch37": "hg19",
            "hg19": "hg19",
            "b37": "hg19",
            "grch38": "hg38",
            "hg38": "hg38",
            "b38": "hg38",
        }
        assembly_key = assembly_aliases.get(assembly_key, assembly_key)
        key_type_value = _clean(key_type).lower() or "auto"
        if key_type_value == "auto":
            key_types = ("vcf", "genomic")
        elif key_type_value in {"vcf", "genomic"}:
            key_types = (key_type_value,)
        else:
            key_types = ()

        empty = {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "input_assembly": _clean(assembly),
            "assembly": assembly_key,
            "input_chrom": _clean(chrom),
            "input_pos": _clean(pos),
            "input_ref": _clean(ref),
            "input_alt": _clean(alt),
            "matched_genomic_key": normalized_key,
            "matched_key_type": "",
            "variant_gof_tag": "",
            "gofcards_accession_id": "",
            "gofcards_variant_id": "",
            "matches": [],
        }
        if not symbol or not normalized_key or assembly_key not in {"hg19", "hg38"} or not key_types:
            return empty

        lookup = self._load_gofcards_variant_genomic(
            gofcards_exact_hgvsp_tsv=gofcards_exact_hgvsp_tsv,
        )
        matches: list[dict[str, str]] = []
        matched_key_types: list[str] = []
        for current_key_type in key_types:
            current_matches = lookup.get((symbol, assembly_key, current_key_type, normalized_key), [])
            if current_matches:
                matched_key_types.append(current_key_type)
                matches.extend(current_matches)
        matches = self._deduplicate_gofcards_matches(matches)
        accession_ids = sorted(
            {row.get("gofcards_accession_id", "") for row in matches if row.get("gofcards_accession_id", "")}
        )
        variant_ids = sorted(
            {row.get("gofcards_variant_id", "") for row in matches if row.get("gofcards_variant_id", "")}
        )
        return {
            **empty,
            "matched_key_type": ";".join(matched_key_types),
            "variant_gof_tag": "GOF" if matches else "",
            "gofcards_accession_id": ";".join(accession_ids),
            "gofcards_variant_id": ";".join(variant_ids),
            "matches": matches,
        }

    @staticmethod
    def _iter_mechanism_entries(info: dict[str, Any]) -> list[dict[str, Any]]:
        entries: list[dict[str, Any]] = []
        for level in ("gene_level", "variant_level"):
            for entry in info.get(level, []) or []:
                for source, data in entry.items():
                    mechanism = _clean(data.get("mechanism"))
                    if not mechanism:
                        continue
                    pmids = [_clean(p) for p in data.get("pmids", []) if _clean(p)]
                    entries.append(
                        {
                            "level": level,
                            "source": source,
                            "mechanism": mechanism,
                            "pmids": pmids,
                            "disease": _clean(data.get("disease")),
                            "consequence": _clean(data.get("consequence")),
                            "chrom": _clean(data.get("chr")),
                            "pos": _clean(data.get("pos")),
                            "ref": _clean(data.get("ref")),
                            "alt": _clean(data.get("alt")),
                            "transcript": _clean(data.get("transcript")),
                        }
                    )
        return entries

    def mechanism_history(
        self,
        gene_symbol: Any,
        *,
        include_entries: bool = False,
    ) -> dict[str, Any]:
        """Return curated GoF/DN/triplosensitivity history for a normalized gene."""
        symbol = self._resolved_symbol_key(gene_symbol)
        info = self._load_mechanisms().get(symbol, {})
        entries = self._iter_mechanism_entries(info) if info else []
        entries.extend(self._load_ddg2p_lof().get(symbol, []))
        mechanism_counts = Counter(entry["mechanism"] for entry in entries)
        pmids_by_mechanism: dict[str, set[str]] = defaultdict(set)
        variant_counts = Counter()
        gene_counts = Counter()
        for entry in entries:
            mechanism = entry["mechanism"]
            pmids_by_mechanism[mechanism].update(entry["pmids"])
            if entry["level"] == "variant_level":
                variant_counts[mechanism] += 1
            elif entry["level"] == "gene_level":
                gene_counts[mechanism] += 1

        out: dict[str, Any] = {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "mechanisms": sorted(mechanism_counts),
            "mechanism_counts": dict(sorted(mechanism_counts.items())),
            "gene_level_counts": dict(sorted(gene_counts.items())),
            "variant_level_counts": dict(sorted(variant_counts.items())),
            "pmids_by_mechanism": {
                key: sorted(value) for key, value in sorted(pmids_by_mechanism.items())
            },
            "has_lof_history": mechanism_counts.get("LOF", 0) > 0,
            "has_gof_history": mechanism_counts.get("GOF", 0) > 0,
            "has_dn_history": mechanism_counts.get("DOMINANT_NEGATIVE", 0) > 0,
            "has_triplosensitivity_history": mechanism_counts.get("TRIPLOSENSITIVITY", 0) > 0,
        }
        if include_entries:
            out["entries"] = entries
        return out

    def ddg2p_lof_history(
        self,
        gene_symbol: Any,
        *,
        include_entries: bool = False,
    ) -> dict[str, Any]:
        """Return DDG2P/G2P LoF mechanism history and allelic requirement flags."""
        symbol = self._resolved_symbol_key(gene_symbol)
        entries = self._load_ddg2p_lof().get(symbol, [])
        inheritance_counts = Counter(_clean(entry.get("inheritance")) for entry in entries)
        recessive = any(
            inheritance in DDG2P_RECESSIVE_LOF_INHERITANCE
            for inheritance in inheritance_counts
        )
        dominant = any(
            inheritance in DDG2P_DOMINANT_LOF_INHERITANCE
            for inheritance in inheritance_counts
        )
        out: dict[str, Any] = {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "has_ddg2p_lof_history": bool(entries),
            "ddg2p_lof_recessive": recessive,
            "ddg2p_lof_dominant": dominant,
            "ddg2p_lof_inheritance_counts": dict(sorted(inheritance_counts.items())),
            "ddg2p_lof_disease_count": len({entry.get("disease", "") for entry in entries if entry.get("disease", "")}),
        }
        if include_entries:
            out["entries"] = entries
        return out

    def classify_categories(
        self,
        gene_symbol: Any,
        *,
        recessive: bool,
        dominant: bool,
        haplo_insufficient: bool,
    ) -> list[str]:
        """Return mechanism categories for one gene using precomputed inheritance."""
        history = self.mechanism_history(gene_symbol)
        return classify_gene_mechanism_categories(
            recessive=bool(recessive),
            dominant=bool(dominant),
            haplo_insufficient=bool(haplo_insufficient),
            has_gof_history=history["has_gof_history"],
            has_dn_history=history["has_dn_history"],
            has_triplosensitivity_history=history["has_triplosensitivity_history"],
        )

    def known_inheritance_mode(
        self,
        gene_symbol: Any,
        *,
        gene_mean_am_score: float = np.nan,
    ) -> dict[str, Any]:
        """Return inheritance/HI calls using PriVA's existing inheritance function."""
        from acmg_criteria_assign import identify_inheritance_mode_per_row

        symbol = self._resolved_symbol_key(gene_symbol)
        hpo_record = self._load_hpo().get(symbol, {})
        clingen_record = self._load_clingen().get(symbol, {})
        loeuf = self._load_loeuf().get(symbol, np.nan)
        clingen_hi_score = clingen_record.get("haploinsufficiency_score")
        ddg2p_lof = self.ddg2p_lof_history(symbol)

        row_dict = {
            "Gene": symbol,
            "SYMBOL": symbol,
            "LOEUF": loeuf,
            "HPO_IDs": hpo_record.get("hpo_id", ""),
            "HPO_gene_inheritance": hpo_record.get("inheritance_modes") or None,
        }
        result = identify_inheritance_mode_per_row(
            row_dict,
            gene_mean_am_score,
            clingen_hi_score,
        )
        recessive, dominant, non_monogenic, non_mendelian, haplo_insufficient, incomplete = result
        effective_recessive = bool(recessive) or bool(ddg2p_lof["ddg2p_lof_recessive"])
        effective_dominant = bool(dominant) or bool(ddg2p_lof["ddg2p_lof_dominant"])
        effective_haplo = bool(haplo_insufficient) or bool(ddg2p_lof["ddg2p_lof_dominant"])
        return {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "recessive": effective_recessive,
            "dominant": effective_dominant,
            "non_monogenic": bool(non_monogenic),
            "non_mendelian": bool(non_mendelian),
            "haplo_insufficient": effective_haplo,
            "incomplete_penetrance": bool(incomplete),
            "hpo_inheritance_modes": hpo_record.get("inheritance_modes", ""),
            "hpo_inheritance_disease_ids": hpo_record.get("inheritance_disease_ids", ""),
            "ddg2p_lof_history": ddg2p_lof,
            "clingen_haploinsufficiency_score": clingen_hi_score,
            "clingen_haploinsufficiency_description": clingen_record.get(
                "haploinsufficiency_description", ""
            ),
            "clingen_triplosensitivity_score": clingen_record.get("triplosensitivity_score"),
            "clingen_triplosensitivity_description": clingen_record.get(
                "triplosensitivity_description", ""
            ),
            "loeuf": None if isinstance(loeuf, float) and math.isnan(loeuf) else loeuf,
            "decision_function": "acmg_criteria_assign.identify_inheritance_mode_per_row",
        }

    def gene_summary(self, gene_symbol: Any, *, include_entries: bool = False) -> dict[str, Any]:
        """Return normalized symbol, mechanism history, and known inheritance summary."""
        symbol = self._resolved_symbol_key(gene_symbol)
        return {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "mechanism_history": self.mechanism_history(symbol, include_entries=include_entries),
            "ddg2p_lof_history": self.ddg2p_lof_history(symbol, include_entries=include_entries),
            "known_inheritance_mode": self.known_inheritance_mode(symbol),
        }


def resolve_gene_symbol(gene_symbol: Any) -> str:
    """Convenience function returning one current HGNC symbol."""
    return GeneMechanismHub().resolve_symbol(gene_symbol)


def annotate_gene_mechanism_categories(
    df: pd.DataFrame,
    *,
    recessive: np.ndarray,
    dominant: np.ndarray,
    haplo_insufficient: np.ndarray,
    mechanism_json: str | Path = DEFAULT_MECHANISM_JSON,
    ddg2p_evidence: str | Path = DEFAULT_DDG2P_MECHANISM_EVIDENCE,
    symbol_col: str = "SYMBOL",
    hpo_inheritance_col: str = "HPO_gene_inheritance",
    output_col: str = "gene_mech_inher_history",
    use_hgnc_package: bool = False,
) -> pd.DataFrame:
    """Annotate a PriVA dataframe with compact gene-level mechanism history.

    The arrays come from PriVA's existing row-level inheritance call path. This
    function also consults the central gene-level inheritance hub so rows with
    missing local HPO inheritance do not collapse to LOEUF-only HI calls. The
    returned dataframe is a shallow copy with one added column:
    ``gene_mech_inher_history``. Values are semicolon-separated category labels
    such as ``LoF_recessive;dominant_HI``. This is gene/history context only and
    does not assign the query variant as GoF/DN.
    """
    if not (len(df) == len(recessive) == len(dominant) == len(haplo_insufficient)):
        raise ValueError("df and inheritance arrays must have the same length")
    if symbol_col not in df.columns:
        raise KeyError(f"missing symbol column: {symbol_col}")

    hub = GeneMechanismHub(
        mechanism_json=mechanism_json,
        ddg2p_evidence=ddg2p_evidence,
        use_hgnc_package=use_hgnc_package,
    )
    out = df.copy()
    category_values: list[str] = []
    known_cache: dict[str, dict[str, Any]] = {}
    row_hpo_values = (
        out[hpo_inheritance_col]
        if hpo_inheritance_col in out.columns
        else pd.Series("", index=out.index)
    )

    for gene, row_hpo, rec, dom, hi in zip(
        out[symbol_col],
        row_hpo_values,
        recessive,
        dominant,
        haplo_insufficient,
        strict=True,
    ):
        symbol = hub.resolve_symbol(gene)
        history = hub.mechanism_history(gene)
        if symbol not in known_cache:
            known_cache[symbol] = hub.known_inheritance_mode(symbol)
        known = known_cache[symbol]
        hub_hpo_flags = _hpo_inheritance_flags(known.get("hpo_inheritance_modes", ""))
        row_hpo_flags = _hpo_inheritance_flags(row_hpo)
        ddg2p_lof = known.get("ddg2p_lof_history", {})
        clingen_hi_score = known.get("clingen_haploinsufficiency_score")
        clingen_hi = clingen_hi_score == 3
        clingen_haplosufficient = clingen_hi_score in {30, 40}

        effective_recessive = (
            bool(rec)
            or row_hpo_flags["recessive"]
            or hub_hpo_flags["recessive"]
            or clingen_haplosufficient
            or bool(ddg2p_lof.get("ddg2p_lof_recessive"))
        )
        effective_dominant = (
            bool(dom)
            or row_hpo_flags["dominant"]
            or hub_hpo_flags["dominant"]
            or clingen_hi
            or bool(ddg2p_lof.get("ddg2p_lof_dominant"))
        )

        categories = set(
            classify_gene_mechanism_categories(
                recessive=effective_recessive,
                dominant=effective_dominant,
                haplo_insufficient=bool(hi) or clingen_hi or bool(ddg2p_lof.get("ddg2p_lof_dominant")),
                has_gof_history=history["has_gof_history"],
                has_dn_history=history["has_dn_history"],
                has_triplosensitivity_history=history["has_triplosensitivity_history"],
            )
        )

        # If the only dominant evidence came from a row-level LOEUF/AM fallback
        # and the hub's HPO/ClinGen view does not support dominance or HI, drop
        # the spurious HI label. This prevents genes with missing row HPO fields
        # from being treated as dominant_HI solely because LOEUF is low.
        if (
            "dominant_HI" in categories
            and bool(hi)
            and not bool(known.get("dominant"))
            and not clingen_hi
            and not row_hpo_flags["dominant"]
            and not hub_hpo_flags["dominant"]
        ):
            categories.discard("dominant_HI")

        category_values.append(";".join(_ordered_gene_mechanism_categories(categories)))

    out[output_col] = category_values
    return out


def match_gofcards_variant_gof(
    gene_symbol: Any,
    hgvsp: Any,
    *,
    gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    gofcards_step1_tsv: str | Path = DEFAULT_GOFCARDS_STEP1_TSV,
    gofcards_active_tsv: str | Path = DEFAULT_GOFCARDS_ACTIVE_TSV,
    gofcards_raw_xlsx: str | Path = DEFAULT_GOFCARDS_RAW_XLSX,
    gofcards_conversion_audit_tsv: str | Path = DEFAULT_GOFCARDS_CONVERSION_AUDIT_TSV,
    use_hgnc_package: bool = False,
) -> dict[str, Any]:
    """Convenience wrapper for exact GoFCards HGVSp variant-level GOF matching."""
    return GeneMechanismHub(use_hgnc_package=use_hgnc_package).match_gofcards_variant_gof(
        gene_symbol,
        hgvsp,
        gofcards_exact_hgvsp_tsv=gofcards_exact_hgvsp_tsv,
        gofcards_step1_tsv=gofcards_step1_tsv,
        gofcards_active_tsv=gofcards_active_tsv,
        gofcards_raw_xlsx=gofcards_raw_xlsx,
        gofcards_conversion_audit_tsv=gofcards_conversion_audit_tsv,
    )


def match_gofcards_variant_gof_by_genomic(
    gene_symbol: Any,
    chrom: Any,
    pos: Any,
    ref: Any,
    alt: Any,
    *,
    assembly: str = "hg38",
    key_type: str = "auto",
    gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    use_hgnc_package: bool = False,
) -> dict[str, Any]:
    """Convenience wrapper for exact GoFCards genomic variant-level GOF matching."""
    return GeneMechanismHub(use_hgnc_package=use_hgnc_package).match_gofcards_variant_gof_by_genomic(
        gene_symbol,
        chrom,
        pos,
        ref,
        alt,
        assembly=assembly,
        key_type=key_type,
        gofcards_exact_hgvsp_tsv=gofcards_exact_hgvsp_tsv,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene", nargs="+", help="Gene symbol/alias/HGNC/Ensembl/Entrez query")
    parser.add_argument("--include-entries", action="store_true")
    parser.add_argument("--mechanism-json", default=str(DEFAULT_MECHANISM_JSON))
    parser.add_argument("--ddg2p-evidence", default=str(DEFAULT_DDG2P_MECHANISM_EVIDENCE))
    parser.add_argument("--hpo-collapsed", default=str(DEFAULT_HPO_COLLAPSED))
    parser.add_argument("--clingen-dosage", default=str(DEFAULT_CLINGEN_DOSAGE))
    parser.add_argument("--loeuf-table", default=str(DEFAULT_LOEUF_TABLE))
    parser.add_argument("--hgnc-table", default=str(DEFAULT_HGNC_TABLE))
    args = parser.parse_args()

    hub = GeneMechanismHub(
        mechanism_json=args.mechanism_json,
        ddg2p_evidence=args.ddg2p_evidence,
        hpo_collapsed=args.hpo_collapsed,
        clingen_dosage=args.clingen_dosage,
        loeuf_table=args.loeuf_table,
        hgnc_table=args.hgnc_table,
    )
    result = {
        gene: hub.gene_summary(gene, include_entries=args.include_entries)
        for gene in args.gene
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
