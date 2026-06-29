#!/usr/bin/env python3
"""Gene-level mechanism and inheritance hub for PriVA.

This module normalizes a gene query to one current HGNC symbol, then reports:

1. Curated non-LoF mechanism history from the existing mechanism JSON cache.
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
DEFAULT_MECHANISM_JSON = Path(
    "/paedyl01/disk1/yangyxt/llm_gene_reranker/data/gene_pathogenic_mechanism/"
    "prepared/gene_mechanism_curated_assertions.json"
)
DEFAULT_GOFCARDS_STEP1_TSV = Path(
    "/paedyl01/disk1/yangyxt/llm_gene_reranker/results/gof_dn_variant_extraction/"
    "priva_output/gofcards_gof.hg19.numb.anno.step1.tsv"
)
DEFAULT_GOFCARDS_ACTIVE_TSV = Path(
    "/paedyl01/disk1/yangyxt/llm_gene_reranker/results/gof_dn_variant_extraction/"
    "active_json_gof_dn_variant_level.tsv"
)
DEFAULT_GOFCARDS_RAW_XLSX = Path(
    "/paedyl01/disk1/yangyxt/public_data/gene_pathogenic_mechanism/raw/gofcards/"
    "gofcards_data_download.xlsx"
)
DEFAULT_GOFCARDS_CONVERSION_AUDIT_TSV = Path(
    "/paedyl01/disk1/yangyxt/llm_gene_reranker/results/gof_dn_variant_extraction/"
    "gofcards_gof.hg19.numb.vcf_conversion_audit.tsv"
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
        if haplo_insufficient:
            categories.append("dominant_HI")
        if has_gof_history or has_triplosensitivity_history:
            categories.append("dominant_GOF")
        if has_dn_history:
            categories.append("dominant_DN")
    return categories


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

    def resolve(self, query: Any) -> str:
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
        hpo_collapsed: str | Path = DEFAULT_HPO_COLLAPSED,
        clingen_dosage: str | Path = DEFAULT_CLINGEN_DOSAGE,
        loeuf_table: str | Path = DEFAULT_LOEUF_TABLE,
        hgnc_table: str | Path = DEFAULT_HGNC_TABLE,
        use_hgnc_package: bool = True,
    ) -> None:
        self.mechanism_json = Path(mechanism_json)
        self.hpo_collapsed = Path(hpo_collapsed)
        self.clingen_dosage = Path(clingen_dosage)
        self.loeuf_table = Path(loeuf_table)
        self.hgnc_table = Path(hgnc_table)

        self._resolver = self._build_resolver(use_hgnc_package)
        self._mechanism_by_symbol: dict[str, dict[str, Any]] | None = None
        self._hpo_by_symbol: dict[str, dict[str, str]] | None = None
        self._clingen_by_symbol: dict[str, dict[str, Any]] | None = None
        self._loeuf_by_symbol: dict[str, float] | None = None
        self._gofcards_by_symbol_hgvsp: dict[tuple[str, str], list[dict[str, str]]] | None = None

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

    def match_gofcards_variant_gof(
        self,
        gene_symbol: Any,
        hgvsp: Any,
        *,
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
            "has_gof_history": mechanism_counts.get("GOF", 0) > 0,
            "has_dn_history": mechanism_counts.get("DOMINANT_NEGATIVE", 0) > 0,
            "has_triplosensitivity_history": mechanism_counts.get("TRIPLOSENSITIVITY", 0) > 0,
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
        return {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "recessive": bool(recessive),
            "dominant": bool(dominant),
            "non_monogenic": bool(non_monogenic),
            "non_mendelian": bool(non_mendelian),
            "haplo_insufficient": bool(haplo_insufficient),
            "incomplete_penetrance": bool(incomplete),
            "hpo_inheritance_modes": hpo_record.get("inheritance_modes", ""),
            "hpo_inheritance_disease_ids": hpo_record.get("inheritance_disease_ids", ""),
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
    symbol_col: str = "SYMBOL",
    output_col: str = "gene_mech_inher_history",
    use_hgnc_package: bool = False,
) -> pd.DataFrame:
    """Annotate a PriVA dataframe with compact gene-level mechanism history.

    The arrays must come from PriVA's existing inheritance call path. The
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
        use_hgnc_package=use_hgnc_package,
    )
    out = df.copy()
    category_values: list[str] = []

    for gene, rec, dom, hi in zip(
        out[symbol_col],
        recessive,
        dominant,
        haplo_insufficient,
        strict=True,
    ):
        history = hub.mechanism_history(gene)
        categories = classify_gene_mechanism_categories(
            recessive=bool(rec),
            dominant=bool(dom),
            haplo_insufficient=bool(hi),
            has_gof_history=history["has_gof_history"],
            has_dn_history=history["has_dn_history"],
            has_triplosensitivity_history=history["has_triplosensitivity_history"],
        )
        category_values.append(";".join(categories))

    out[output_col] = category_values
    return out


def match_gofcards_variant_gof(
    gene_symbol: Any,
    hgvsp: Any,
    *,
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
        gofcards_step1_tsv=gofcards_step1_tsv,
        gofcards_active_tsv=gofcards_active_tsv,
        gofcards_raw_xlsx=gofcards_raw_xlsx,
        gofcards_conversion_audit_tsv=gofcards_conversion_audit_tsv,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene", nargs="+", help="Gene symbol/alias/HGNC/Ensembl/Entrez query")
    parser.add_argument("--include-entries", action="store_true")
    parser.add_argument("--mechanism-json", default=str(DEFAULT_MECHANISM_JSON))
    parser.add_argument("--hpo-collapsed", default=str(DEFAULT_HPO_COLLAPSED))
    parser.add_argument("--clingen-dosage", default=str(DEFAULT_CLINGEN_DOSAGE))
    parser.add_argument("--loeuf-table", default=str(DEFAULT_LOEUF_TABLE))
    parser.add_argument("--hgnc-table", default=str(DEFAULT_HGNC_TABLE))
    args = parser.parse_args()

    hub = GeneMechanismHub(
        mechanism_json=args.mechanism_json,
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
