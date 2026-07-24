#!/usr/bin/env python3
"""Gene-condition and query-variant mechanism hub for PriVA.

This module normalizes a gene query to one current HGNC symbol, then reports:

1. Curated mechanism history from the existing mechanism JSON cache and the
   companion DDG2P/G2P evidence table.
2. Known inheritance/HI status using PriVA's existing inheritance decision
   function from ``acmg_criteria_assign.py``.

Gene history and query-variant effect are deliberately represented separately.
The dataframe annotator combines them only in explicit applicability fields, so
an unrelated mechanism elsewhere in the gene is not assigned to every variant.
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
PACKAGED_MECHANISM_JSON = (
    DATA_DIR / "gene_pathogenic_mechanism" / "prepared" / "gene_mechanism_curated_assertions.json"
)
SHARED_CANONICAL_MECHANISM_JSON = (
    PRIVA_ROOT.parent
    / "llm_gene_reranker"
    / "data"
    / "gene_pathogenic_mechanism"
    / "prepared"
    / "gene_mechanism_curated_assertions.json"
)
DEFAULT_MECHANISM_JSON = Path(
    os.environ.get(
        "PRIVA_GENE_MECHANISM_JSON",
        str(
            SHARED_CANONICAL_MECHANISM_JSON
            if SHARED_CANONICAL_MECHANISM_JSON.exists()
            else PACKAGED_MECHANISM_JSON
        ),
    )
)
DEFAULT_DDG2P_MECHANISM_EVIDENCE = (
    DATA_DIR / "gene_pathogenic_mechanism" / "prepared" / "gene_pathogenic_mechanism_evidence.tsv"
)
DEFAULT_GOFCARDS_EXACT_GOF_HGVSP = DATA_DIR / "gofcards" / "gofcards_exact_gof_hgvsp.tsv.gz"
DEFAULT_GOFCARDS_RAW_XLSX = (
    DATA_DIR / "gene_pathogenic_mechanism" / "raw" / "gofcards" / "gofcards_data_download.xlsx"
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

CANONICAL_MECHANISMS = {"LOF", "GOF", "DOMINANT_NEGATIVE", "TRIPLOSENSITIVITY"}
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
VARIANT_MECHANISM_OUTPUT_COLUMNS = (
    "gene_lof_evidence",
    "variant_effect",
    "variant_effect_evidence",
    "variant_effect_conflict",
    "variant_mechanism_applicable",
    "variant_mechanism_uncertain",
    "variant_mechanism_incompatible",
    "variant_mechanism_applicability_detail",
    "clinvar_vcv_accessions",
    "clinvar_rcv_conditions",
    "clinvar_vcv_max_review_stars",
    "clinvar_vcv_hgvs",
)
logger = logging.getLogger(__name__)


def _clean(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "<na>"}:
        return ""
    return text


def _bool_value(value: Any) -> bool:
    """Coerce common dataframe boolean encodings without treating NaN as true."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return value != 0
    return _clean(value).lower() in {"1", "true", "t", "yes", "y"}


def _split_annotation_tokens(value: Any) -> set[str]:
    return {
        token.strip().upper()
        for token in re.split(r"[&,;|]", _clean(value))
        if token.strip()
    }


def infer_query_variant_effect(row: dict[str, Any] | pd.Series) -> dict[str, str]:
    """Infer the query allele's primary effect and retain every evidence route.

    An exact curated GoFCards match takes precedence over generic consequence
    prediction. LOFTEE ``HC`` and ``OS`` are both high-confidence predicted
    loss-of-function evidence; ``OS`` is LOFTEE's other-splice category.
    """
    evidence: list[str] = []
    loftee_tokens = _split_annotation_tokens(row.get("LoF"))
    if "HC" in loftee_tokens:
        evidence.append("LOFTEE_HC")
    if "OS" in loftee_tokens:
        evidence.append("LOFTEE_OS")

    consequence = _clean(row.get("Consequence")).lower()
    nmd = _clean(row.get("NMD")).lower()
    lof_filter = _clean(row.get("LoF_filter")).upper()
    truncating = "stop_gained" in consequence or "frameshift" in consequence
    nmd_lof = truncating and "escaping" not in nmd and "END_TRUNC" not in lof_filter
    if nmd_lof:
        evidence.append("NMD_PREDICTED_LOF")
    if _bool_value(row.get("vep_consq_lof")):
        evidence.append("VEP_LOF")
    if _bool_value(row.get("splicing_lof")):
        evidence.append("PRIVA_SPLICE_LOF")
    if _bool_value(row.get("5UTR_lof")):
        evidence.append("PRIVA_5UTR_LOF")

    gof_tokens = _split_annotation_tokens(row.get("variant_gof_tag"))
    exact_gof = "GOF" in gof_tokens
    if exact_gof:
        evidence.append("GOFCARDS_EXACT")

    evidence = list(dict.fromkeys(evidence))
    lof_evidence = [item for item in evidence if item != "GOFCARDS_EXACT"]
    if exact_gof:
        effect = "exact_known_GOF"
        conflict = "predicted_LOF_vs_exact_GOF" if lof_evidence else ""
    elif lof_evidence:
        effect = "predicted_LOF_high_confidence"
        conflict = ""
    else:
        effect = "uncertain"
        conflict = ""
    return {
        "variant_effect": effect,
        "variant_effect_evidence": ";".join(evidence),
        "variant_effect_conflict": conflict,
    }


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


# The former gene-only category classifier is deliberately disabled. Mechanism
# and inheritance are now resolved together in variant-level assertions.


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
                "allelic_requirement": inheritance,
                "confidence": mechanism_confidence,
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
        gofcards_step1_tsv: str | Path | None = None,
        gofcards_active_tsv: str | Path | None = None,
        gofcards_raw_xlsx: str | Path | None = None,
        gofcards_conversion_audit_tsv: str | Path | None = None,
    ) -> dict[tuple[str, str], list[dict[str, str]]]:
        """Return exact-match GoFCards GOF protein-change evidence.

        The key is normalized ``(HGNC symbol, protein change)``. This lookup is
        only for variant-level evidence; it must not be used to convert all
        variants in a gene into GOF assertions.

        The compact PriVA TSV is the only supported runtime source. The older
        ``gofcards_step1_tsv``/``gofcards_active_tsv``/audit arguments are
        retained only so external callers do not break, but they are ignored.
        """
        if self._gofcards_by_symbol_hgvsp is not None:
            return self._gofcards_by_symbol_hgvsp

        compact_path = Path(gofcards_exact_hgvsp_tsv)
        if not compact_path.exists():
            logger.warning(
                "GoFCards exact GOF compact cache not found: %s. "
                "Run `bash scripts/install_utils.sh gofcards_exact_gof_cache_install config.yaml`.",
                compact_path,
            )
            self._gofcards_by_symbol_hgvsp = {}
            return self._gofcards_by_symbol_hgvsp

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
        gofcards_step1_tsv: str | Path | None = None,
        gofcards_active_tsv: str | Path | None = None,
        gofcards_raw_xlsx: str | Path | None = None,
        gofcards_conversion_audit_tsv: str | Path | None = None,
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
                            "allelic_requirement": _clean(
                                data.get("inheritance")
                                or data.get("allelic_requirement")
                            ),
                            "confidence": _clean(
                                data.get("confidence")
                                or data.get("mechanism_confidence")
                            ),
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
            "condition_mechanism_assertions": self.condition_mechanism_assertions(symbol),
            "ddg2p_lof_history": self.ddg2p_lof_history(symbol, include_entries=include_entries),
            "known_inheritance_mode": self.known_inheritance_mode(symbol),
        }

    def condition_mechanism_assertions(self, gene_symbol: Any) -> list[dict[str, Any]]:
        """Return condition-specific curated mechanism/allelic assertions."""
        history = self.mechanism_history(gene_symbol, include_entries=True)
        assertions: list[dict[str, Any]] = []
        seen: set[tuple[str, str, str, str, str]] = set()
        for entry in history.get("entries", []):
            mechanism = _clean(entry.get("mechanism")).upper()
            if mechanism not in CANONICAL_MECHANISMS:
                continue
            # Other curated variants in the same gene establish history, but
            # they are not query-variant condition assertions. GoFCards exact
            # matches are added from the query row below.
            if entry.get("level") == "variant_level":
                continue
            record = {
                "source": _clean(entry.get("source")),
                "disease": _clean(entry.get("disease")),
                "mechanism": mechanism,
                "allelic_requirement": _clean(entry.get("allelic_requirement")),
                "confidence": _clean(entry.get("confidence")),
            }
            key = tuple(record[field] for field in (
                "source",
                "disease",
                "mechanism",
                "allelic_requirement",
                "confidence",
            ))
            if key not in seen:
                seen.add(key)
                assertions.append(record)
        return assertions

    def matched_clinvar_vcv_for_gofcards(
        self,
        gene_symbol: Any,
        *,
        gofcards_variant_ids: Any = "",
        gofcards_accession_ids: Any = "",
    ) -> list[dict[str, Any]]:
        """Return schema-v2 ClinVar VCV records linked to exact GoFCards IDs."""
        variant_ids = {
            _clean(token)
            for token in re.split(r"[;,]", _clean(gofcards_variant_ids))
            if _clean(token)
        }
        accession_ids = {
            _clean(token)
            for token in re.split(r"[;,]", _clean(gofcards_accession_ids))
            if _clean(token)
        }
        if not variant_ids and not accession_ids:
            return []

        symbol = self._resolved_symbol_key(gene_symbol)
        info = self._load_mechanisms().get(symbol, {})
        matches: list[dict[str, Any]] = []
        for entry in info.get("variant_level", []) or []:
            data = entry.get("ClinVar_VCV") if isinstance(entry, dict) else None
            if not isinstance(data, dict):
                continue
            linked_records = data.get("match", {}).get("matched_gofcards_records", []) or []
            linked_variant_ids = {
                _clean(record.get("gofcards_variant_id"))
                for record in linked_records
                if isinstance(record, dict) and _clean(record.get("gofcards_variant_id"))
            }
            linked_accessions = {
                _clean(record.get("gofcards_accession_id"))
                for record in linked_records
                if isinstance(record, dict) and _clean(record.get("gofcards_accession_id"))
            }
            if not (
                variant_ids & linked_variant_ids
                or accession_ids & linked_accessions
            ):
                continue
            matches.append(data)
        return matches


def resolve_gene_symbol(gene_symbol: Any) -> str:
    """Convenience function returning one current HGNC symbol."""
    return GeneMechanismHub().resolve_symbol(gene_symbol)


def _hpo_allelic_requirements(inheritance_modes: Any) -> set[str]:
    """Collapse HPO inheritance to the supported compact assertion vocabulary.

    HPO establishes only a broad dominant or recessive history here. It does
    not establish autosomal/X-linked dosage or a molecular mechanism.
    """
    modes = {
        _clean(part)
        for part in _clean(inheritance_modes).split(";")
        if _clean(part)
    }
    requirements: set[str] = set()
    if {
        "Autosomal recessive inheritance",
        "X-linked inheritance",
        "X-linked recessive inheritance",
        "Pseudoautosomal recessive inheritance",
    } & modes:
        requirements.add("recessive")
    if {
        "Autosomal dominant inheritance",
        "Autosomal dominant inheritance with maternal imprinting",
        "X-linked dominant inheritance",
        "Y-linked inheritance",
    } & modes:
        requirements.add("dominant")
    if "Mitochondrial inheritance" in modes:
        requirements.add("mitochondrial")
    return requirements


def _compact_inheritance(allelic_requirement: Any) -> str:
    requirement = _clean(allelic_requirement)
    if requirement in {"recessive", "dominant", "mitochondrial"}:
        return requirement
    if requirement.startswith("biallelic_") or requirement in {
        "monoallelic_X",
        "monoallelic_X_hemizygous",
    }:
        return "recessive"
    if requirement.startswith("monoallelic_"):
        return "dominant"
    return ""


def _mechanism_profile_tag(assertion: dict[str, Any]) -> str:
    inheritance = _compact_inheritance(assertion.get("allelic_requirement"))
    mechanism = _clean(assertion.get("mechanism")).upper() or "UNRESOLVED"
    if mechanism == "UNRESOLVED":
        return inheritance or "uncertain"
    mechanism = "DN" if mechanism == "DOMINANT_NEGATIVE" else mechanism
    return f"{inheritance}_{mechanism}" if inheritance else mechanism


def _compact_profile_tags(assertions: list[dict[str, Any]]) -> list[str]:
    tags = list(dict.fromkeys(_mechanism_profile_tag(assertion) for assertion in assertions))
    inheritance_with_mechanism = {
        tag.split("_", 1)[0]
        for tag in tags
        if "_" in tag and tag.split("_", 1)[0] in {"recessive", "dominant"}
    }
    qualified_mechanisms = {
        tag.split("_", 1)[1]
        for tag in tags
        if "_" in tag and tag.split("_", 1)[0] in {"recessive", "dominant"}
    }
    return [
        tag
        for tag in tags
        if tag not in inheritance_with_mechanism
        and tag not in qualified_mechanisms
        and tag != "uncertain"
    ]


def _deduplicate_assertions(assertions: list[dict[str, Any]]) -> list[dict[str, Any]]:
    fields = ("source", "disease", "mechanism", "allelic_requirement", "confidence")
    seen: set[tuple[str, ...]] = set()
    output: list[dict[str, Any]] = []
    for assertion in assertions:
        normalized = {
            "source": _clean(assertion.get("source")),
            "disease": _clean(assertion.get("disease")),
            "mechanism": _clean(assertion.get("mechanism")).upper() or "UNRESOLVED",
            "allelic_requirement": _clean(assertion.get("allelic_requirement")),
            "confidence": _clean(assertion.get("confidence")),
        }
        key = tuple(normalized[field] for field in fields)
        if key not in seen:
            seen.add(key)
            output.append(normalized)
    return output


def _classify_variant_applicability(
    assertions: list[dict[str, Any]],
    effect: str,
) -> dict[str, Any]:
    groups: dict[str, list[str]] = {
        "applicable": [],
        "uncertain": [],
        "incompatible": [],
    }
    details: list[dict[str, str]] = []
    for assertion in assertions:
        mechanism = assertion["mechanism"]
        if mechanism == "LOF":
            if effect == "predicted_LOF_high_confidence":
                status, reason = "applicable", "query_effect_matches_LOF"
            elif effect == "uncertain":
                status, reason = "uncertain", "query_LOF_effect_not_established"
            else:
                status, reason = "incompatible", "exact_GOF_does_not_match_LOF"
        elif mechanism == "GOF":
            if effect == "exact_known_GOF":
                status, reason = "applicable", "exact_query_GOF_match"
            else:
                status, reason = "incompatible", "GOF_requires_exact_variant_match"
        elif mechanism == "DOMINANT_NEGATIVE":
            if effect == "uncertain":
                status, reason = "uncertain", "no_variant_level_DN_assertion"
            else:
                status, reason = "incompatible", "query_effect_does_not_support_DN"
        elif mechanism == "TRIPLOSENSITIVITY":
            status, reason = "uncertain", "sequence_variant_not_equivalent_to_copy_gain"
        else:
            status, reason = "uncertain", "inheritance_known_mechanism_unresolved"

        tag = _mechanism_profile_tag(assertion)
        if tag not in groups[status]:
            groups[status].append(tag)
        details.append(
            {
                **assertion,
                "tag": tag,
                "applicability": status,
                "reason": reason,
            }
        )
    return {
        "plausible": ";".join(
            _compact_profile_tags(
                [
                    detail
                    for detail in details
                    if detail["applicability"] != "incompatible"
                ]
            )
        ),
        "applicable": ";".join(groups["applicable"]),
        "uncertain": ";".join(groups["uncertain"]),
        "incompatible": ";".join(groups["incompatible"]),
        "detail": json.dumps(details, sort_keys=True, separators=(",", ":")),
    }


def annotate_gene_mechanism_categories(
    df: pd.DataFrame,
    *,
    clinvar_pathogenic_genes: set[str] | None = None,
    gene_to_am_score_map: dict[str, float] | None = None,
    mechanism_json: str | Path = DEFAULT_MECHANISM_JSON,
    ddg2p_evidence: str | Path = DEFAULT_DDG2P_MECHANISM_EVIDENCE,
    symbol_col: str = "SYMBOL",
    gene_col: str = "Gene",
    hpo_inheritance_col: str = "HPO_gene_inheritance",
    output_col: str = "var_plausible_patho_mechs",
    use_hgnc_package: bool = False,
    hpo_collapsed: str | Path = DEFAULT_HPO_COLLAPSED,
    clingen_dosage: str | Path = DEFAULT_CLINGEN_DOSAGE,
    loeuf_table: str | Path = DEFAULT_LOEUF_TABLE,
    hgnc_table: str | Path = DEFAULT_HGNC_TABLE,
) -> pd.DataFrame:
    """Annotate compact gene history and query-variant applicability.

    HPO contributes only ``recessive`` or ``dominant``. A generic inheritance
    assertion gains a LOF mechanism only when the gene has a pathogenic
    ClinVar variant reviewed at two stars or above, LOEUF < 0.35, or mean
    AlphaMissense score > 0.564.
    """
    if symbol_col not in df.columns:
        raise KeyError(f"missing symbol column: {symbol_col}")
    if gene_col not in df.columns:
        raise KeyError(f"missing gene column: {gene_col}")

    hub = GeneMechanismHub(
        mechanism_json=mechanism_json,
        ddg2p_evidence=ddg2p_evidence,
        hpo_collapsed=hpo_collapsed,
        clingen_dosage=clingen_dosage,
        loeuf_table=loeuf_table,
        hgnc_table=hgnc_table,
        use_hgnc_package=use_hgnc_package,
    )
    out = df.copy()
    plausible_mechanism_values: list[str] = []
    variant_outputs: dict[str, list[str]] = {
        column: [] for column in VARIANT_MECHANISM_OUTPUT_COLUMNS
    }
    known_cache: dict[str, dict[str, Any]] = {}
    assertion_cache: dict[str, list[dict[str, Any]]] = {}
    clinvar_pathogenic_genes = set(clinvar_pathogenic_genes or set())
    gene_to_am_score_map = gene_to_am_score_map or {}
    row_hpo_values = (
        out[hpo_inheritance_col]
        if hpo_inheritance_col in out.columns
        else pd.Series("", index=out.index)
    )

    for (_, row), row_hpo in zip(
        out.iterrows(),
        row_hpo_values,
        strict=True,
    ):
        gene = row[symbol_col]
        symbol = hub.resolve_symbol(gene)
        if symbol not in assertion_cache:
            assertion_cache[symbol] = hub.condition_mechanism_assertions(gene)
        if symbol not in known_cache:
            known_cache[symbol] = hub.known_inheritance_mode(symbol)
        known = known_cache[symbol]
        clingen_hi_score = known.get("clingen_haploinsufficiency_score")
        clingen_hi = clingen_hi_score == 3
        clingen_haplosufficient = clingen_hi_score in {30, 40}

        assertions = list(assertion_cache[symbol])
        hpo_requirements = _hpo_allelic_requirements(row_hpo)
        hpo_requirements.update(
            _hpo_allelic_requirements(known.get("hpo_inheritance_modes", ""))
        )
        for requirement in sorted(hpo_requirements):
            assertions.append(
                {
                    "source": "HPO_inheritance",
                    "disease": known.get("hpo_inheritance_disease_ids", ""),
                    "mechanism": "UNRESOLVED",
                    "allelic_requirement": requirement,
                    "confidence": "inheritance_only",
                }
            )

        gene_id = _clean(row.get(gene_col))
        lof_evidence: list[str] = []
        if gene_id in clinvar_pathogenic_genes:
            lof_evidence.append("ClinVar_pathogenic_2plus")
        loeuf = _safe_float(row.get("LOEUF"))
        if not math.isnan(loeuf) and loeuf < 0.35:
            lof_evidence.append("LOEUF_lt_0.35")
        gene_avg_am = _safe_float(
            row.get("Gene_avg_AM_score", gene_to_am_score_map.get(gene_id))
        )
        if not math.isnan(gene_avg_am) and gene_avg_am > 0.564:
            lof_evidence.append("GeneAvgAM_gt_0.564")
        if lof_evidence:
            lof_requirements = hpo_requirements or {""}
            for requirement in sorted(lof_requirements):
                assertions.append(
                    {
                        "source": "PriVA_gene_LOF_evidence",
                        "disease": known.get("hpo_inheritance_disease_ids", ""),
                        "mechanism": "LOF",
                        "allelic_requirement": requirement,
                        "confidence": ";".join(lof_evidence),
                    }
                )
        if clingen_hi:
            assertions.append(
                {
                    "source": "ClinGen_Dosage",
                    "disease": "",
                    "mechanism": "LOF",
                    "allelic_requirement": "monoallelic_autosomal",
                    "confidence": "sufficient_HI_evidence",
                }
            )
        elif clingen_haplosufficient:
            assertions.append(
                {
                    "source": "ClinGen_Dosage",
                    "disease": "",
                    "mechanism": "LOF",
                    "allelic_requirement": "biallelic_autosomal",
                    "confidence": "AR_gene_association",
                }
            )
        effect_call = infer_query_variant_effect(row)
        vcv_accessions: list[str] = []
        vcv_conditions: list[str] = []
        vcv_hgvs: list[str] = []
        vcv_review_stars: list[int] = []
        if effect_call["variant_effect"] == "exact_known_GOF":
            assertions.append(
                {
                    "source": "GoFCards",
                    "disease": "",
                    "mechanism": "GOF",
                    "allelic_requirement": "",
                    "confidence": "exact_variant_match",
                }
            )
            vcv_matches = hub.matched_clinvar_vcv_for_gofcards(
                symbol,
                gofcards_variant_ids=row.get("gofcards_variant_id", ""),
                gofcards_accession_ids=row.get("gofcards_accession_id", ""),
            )
            for vcv in vcv_matches:
                variation = vcv.get("variation", {}) or {}
                accession = _clean(variation.get("vcv_accession"))
                if accession:
                    vcv_accessions.append(accession)
                for hgvs in variation.get("hgvs", []) or []:
                    if isinstance(hgvs, dict) and _clean(hgvs.get("expression")):
                        vcv_hgvs.append(_clean(hgvs.get("expression")))
                for condition_assertion in vcv.get("condition_assertions", []) or []:
                    if not isinstance(condition_assertion, dict):
                        continue
                    classification = condition_assertion.get("germline_classification", {}) or {}
                    stars = classification.get("review_stars")
                    try:
                        vcv_review_stars.append(int(stars))
                    except (TypeError, ValueError):
                        pass
                    condition_names = [
                        _clean(condition.get("name"))
                        for condition in condition_assertion.get("conditions", []) or []
                        if isinstance(condition, dict) and _clean(condition.get("name"))
                    ]
                    vcv_conditions.extend(condition_names)
                    assertions.append(
                        {
                            "source": "GoFCards_exact+ClinVar_VCV",
                            "disease": ";".join(condition_names),
                            "mechanism": "GOF",
                            "allelic_requirement": "",
                            "confidence": (
                                f"ClinVar_{stars}_star" if stars is not None else "ClinVar_unrated"
                            ),
                        }
                    )
        assertions = _deduplicate_assertions(assertions)
        applicability = _classify_variant_applicability(
            assertions,
            effect_call["variant_effect"],
        )

        plausible_mechanism_values.append(applicability["plausible"])
        variant_outputs["gene_lof_evidence"].append(";".join(lof_evidence))
        variant_outputs["variant_effect"].append(effect_call["variant_effect"])
        variant_outputs["variant_effect_evidence"].append(effect_call["variant_effect_evidence"])
        variant_outputs["variant_effect_conflict"].append(effect_call["variant_effect_conflict"])
        variant_outputs["variant_mechanism_applicable"].append(applicability["applicable"])
        variant_outputs["variant_mechanism_uncertain"].append(applicability["uncertain"])
        variant_outputs["variant_mechanism_incompatible"].append(applicability["incompatible"])
        variant_outputs["variant_mechanism_applicability_detail"].append(applicability["detail"])
        variant_outputs["clinvar_vcv_accessions"].append(
            ";".join(dict.fromkeys(vcv_accessions))
        )
        variant_outputs["clinvar_rcv_conditions"].append(
            ";".join(dict.fromkeys(vcv_conditions))
        )
        variant_outputs["clinvar_vcv_max_review_stars"].append(
            str(max(vcv_review_stars)) if vcv_review_stars else ""
        )
        variant_outputs["clinvar_vcv_hgvs"].append(
            ";".join(dict.fromkeys(vcv_hgvs))
        )

    out[output_col] = plausible_mechanism_values
    for column, values in variant_outputs.items():
        out[column] = values
    return out


def match_gofcards_variant_gof(
    gene_symbol: Any,
    hgvsp: Any,
    *,
    gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    gofcards_step1_tsv: str | Path | None = None,
    gofcards_active_tsv: str | Path | None = None,
    gofcards_raw_xlsx: str | Path | None = None,
    gofcards_conversion_audit_tsv: str | Path | None = None,
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
