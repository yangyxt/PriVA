"""Lazy resource loading and the public GeneMechanismHub service.

The class composes condition and variant definitions without importing the
top-level facade. This keeps dependency direction one-way.
"""

from __future__ import annotations

import csv
import json
import logging
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from acmg_inheritance import identify_inheritance_mode_per_row, parse_hpo_inheritance
from clinvar_vcv import open_text
from gene_mechanism_common import (
    CANONICAL_MECHANISMS,
    CONDITION_MECHANISM_EVIDENCE_COLUMNS,
    CONDITION_MECHANISM_SOURCES,
    DDG2P_DOMINANT_LOF_INHERITANCE,
    DDG2P_LOF_RAW_VALUES,
    DDG2P_RECESSIVE_LOF_INHERITANCE,
    DDG2P_X_LINKED_DOMINANT_LOF_INHERITANCE,
    DDG2P_X_LINKED_UNSPECIFIED_LOF_INHERITANCE,
    DDG2P_Y_LINKED_LOF_INHERITANCE,
    DDG2P_USABLE_DISEASE_CONFIDENCE,
    DDG2P_USABLE_MECHANISM_CONFIDENCE,
    DEFAULT_CLINGEN_DOSAGE,
    DEFAULT_DDG2P_MECHANISM_EVIDENCE,
    DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    DEFAULT_HGNC_TABLE,
    DEFAULT_HPO_COLLAPSED,
    DEFAULT_HPO_CONDITION_MECHANISM_CACHE,
    DEFAULT_LOEUF_TABLE,
    DEFAULT_MECHANISM_JSON,
    HPO_CONDITION_MECHANISM_SCHEMA_VERSION,
    HPO_INHERITANCE_TERMS,
    HPO_ONSET_TERMS,
    HPO_PENETRANCE_TERMS,
    LOOKUP_FIELD_PRIORITY,
    _clean,
    _norm,
    _safe_float,
    _safe_int,
    _split_multi,
    _split_pmids,
)
from gene_mechanism_conditions import (
    _hpo_inheritance_flags,
    build_hpo_condition_index,
    condition_cache_mechanism_assertions,
    condition_cache_mechanism_entries,
    enrich_condition_mechanism_assertion,
    normalize_inheritance_histories,
)
from gene_mechanism_variants import ExactVariantMatcher


logger = logging.getLogger(__name__)


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
        condition_cache: str | Path = DEFAULT_HPO_CONDITION_MECHANISM_CACHE,
        mechanism_json: str | Path = DEFAULT_MECHANISM_JSON,
        ddg2p_evidence: str | Path = DEFAULT_DDG2P_MECHANISM_EVIDENCE,
        hpo_collapsed: str | Path = DEFAULT_HPO_COLLAPSED,
        clingen_dosage: str | Path = DEFAULT_CLINGEN_DOSAGE,
        loeuf_table: str | Path = DEFAULT_LOEUF_TABLE,
        hgnc_table: str | Path = DEFAULT_HGNC_TABLE,
        use_hgnc_package: bool = True,
    ) -> None:
        self.condition_cache = Path(condition_cache)
        self.mechanism_json = Path(mechanism_json)
        self.ddg2p_evidence = Path(ddg2p_evidence) if ddg2p_evidence else Path("")
        self.hpo_collapsed = Path(hpo_collapsed)
        self.clingen_dosage = Path(clingen_dosage)
        self.loeuf_table = Path(loeuf_table)
        self.hgnc_table = Path(hgnc_table)

        self._resolver = self._build_resolver(use_hgnc_package)
        self._condition_cache_by_symbol: dict[str, dict[str, Any]] | None = None
        self._condition_cache_meta: dict[str, Any] | None = None
        self._mechanism_by_symbol: dict[str, dict[str, Any]] | None = None
        self._condition_mechanism_by_symbol: dict[str, list[dict[str, Any]]] | None = None
        self._ddg2p_lof_by_symbol: dict[str, list[dict[str, Any]]] | None = None
        self._hpo_by_symbol: dict[str, dict[str, str]] | None = None
        self._hpo_condition_index: dict[tuple[str, str], dict[str, Any]] | None = None
        self._clingen_by_symbol: dict[str, dict[str, Any]] | None = None
        self._loeuf_by_symbol: dict[str, float] | None = None
        self._exact_variant_matcher: ExactVariantMatcher | None = None

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

    def _load_condition_cache(self) -> dict[str, dict[str, Any]]:
        """Load and normalize PriVA's integrated condition-mechanism cache.

        The builder writes one atomic JSON document with a versioned top-level
        contract. Runtime loading is deliberately strict: using a stale schema
        could silently move inheritance or penetrance between conditions, which
        is not acceptable for ACMG criteria. Gene keys are indexed under both
        their cache spelling and current HGNC symbol so historical aliases do
        not force callers to understand resource-specific naming.
        """
        if self._condition_cache_by_symbol is not None:
            return self._condition_cache_by_symbol
        if not self.condition_cache.is_file():
            raise FileNotFoundError(
                "Integrated HPO condition-mechanism cache not found: "
                f"{self.condition_cache}. Run "
                "`bash scripts/install_utils.sh "
                "hpo_condition_mechanism_cache_install config.yaml`."
            )

        with open_text(self.condition_cache) as handle:
            payload = json.load(handle)
        if not isinstance(payload, dict):
            raise ValueError(f"{self.condition_cache} must contain a JSON object")
        meta = payload.get("_meta")
        if not isinstance(meta, dict):
            raise ValueError(f"{self.condition_cache} is missing object _meta")
        schema_version = _clean(meta.get("schema_version"))
        if schema_version != HPO_CONDITION_MECHANISM_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported condition cache schema {schema_version!r}; expected "
                f"{HPO_CONDITION_MECHANISM_SCHEMA_VERSION!r}"
            )
        genes = payload.get("genes")
        if not isinstance(genes, dict):
            raise ValueError(f"{self.condition_cache} is missing object genes")

        by_symbol: dict[str, dict[str, Any]] = {}
        for cache_symbol, gene in genes.items():
            if not isinstance(gene, dict):
                raise ValueError(
                    f"Condition cache gene {cache_symbol!r} must be an object"
                )
            symbol = _clean(cache_symbol)
            if not symbol:
                continue
            by_symbol[symbol] = gene
            resolved = self._try_resolve_symbol(symbol)
            if resolved:
                by_symbol.setdefault(resolved, gene)

        self._condition_cache_meta = meta
        self._condition_cache_by_symbol = by_symbol
        return by_symbol

    def _load_mechanisms(self) -> dict[str, dict[str, Any]]:
        if self._mechanism_by_symbol is not None:
            return self._mechanism_by_symbol
        with open_text(self.mechanism_json) as handle:
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

    def _load_condition_mechanism_evidence(
        self,
    ) -> dict[str, list[dict[str, Any]]]:
        """Load condition-resolved G2P and Orphadata mechanism assertions.

        The configured mechanism JSON and this TSV have deliberately different
        responsibilities. The shared JSON preserves exact ClinVar/GoFCards
        variant relationships, whereas the TSV is the canonical runtime source
        for G2P and Orphadata's gene-condition-mechanism assertions. Reading the
        TSV here lets PriVA keep both resources without replacing the richer
        shared JSON or copying condition records into it.

        A complete identity schema is mandatory. Older evidence tables that do
        not contain source disease IDs, MONDO mappings, or PriVA disease-scope
        fields are rejected rather than being interpreted gene-wide. Such old
        tables remain usable by ``_load_ddg2p_lof`` for the narrowly defined
        legacy PVS1/inheritance audit, but they cannot transfer condition-level
        inheritance, onset, or penetrance.

        One source row can encode more than one normalized mechanism. Each
        canonical mechanism is emitted as its own record so downstream variant
        matching can select GOF or LOF without carrying an unrelated mechanism.
        No disease-name, PMID, or phenotype-similarity matching occurs here.
        """
        if self._condition_mechanism_by_symbol is not None:
            return self._condition_mechanism_by_symbol

        by_symbol: dict[str, list[dict[str, Any]]] = defaultdict(list)
        if not self.ddg2p_evidence or not self.ddg2p_evidence.exists():
            self._condition_mechanism_by_symbol = {}
            return self._condition_mechanism_by_symbol

        try:
            with open(self.ddg2p_evidence, encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                missing = sorted(
                    CONDITION_MECHANISM_EVIDENCE_COLUMNS
                    - set(reader.fieldnames or [])
                )
                if missing:
                    logger.warning(
                        "Condition mechanism evidence table %s has a stale schema; "
                        "missing columns: %s",
                        self.ddg2p_evidence,
                        ", ".join(missing),
                    )
                    self._condition_mechanism_by_symbol = {}
                    return self._condition_mechanism_by_symbol

                seen_by_symbol: dict[str, set[tuple[str, ...]]] = defaultdict(set)
                for row in reader:
                    source = _clean(row.get("source"))
                    symbol = _clean(row.get("gene_symbol"))
                    if source not in CONDITION_MECHANISM_SOURCES or not symbol:
                        continue

                    mechanisms = {
                        token.strip().upper()
                        for token in re.split(
                            r"[|;,]", _clean(row.get("normalized_mechanisms"))
                        )
                        if token.strip()
                    } & CANONICAL_MECHANISMS
                    if not mechanisms:
                        continue

                    keys = {symbol}
                    resolved = self._try_resolve_symbol(symbol)
                    if resolved:
                        keys.add(resolved)
                    for mechanism in sorted(mechanisms):
                        record = {
                            "level": "gene_level",
                            "source": source,
                            "mechanism": mechanism,
                            "mechanism_raw": _clean(row.get("patho_mode_raw")),
                            "pmids": _split_pmids(row.get("pmids", "")),
                            "source_record_id": _clean(
                                row.get("source_record_id")
                            ),
                            "source_condition_id": _clean(
                                row.get("source_condition_id")
                            ),
                            "mondo_id": _clean(row.get("mondo_id")),
                            "disease": _clean(row.get("disease_label")),
                            "disease_scope": _clean(row.get("disease_scope")),
                            "priva_scope": _clean(row.get("priva_scope")),
                            "scope_review_status": _clean(
                                row.get("scope_review_status")
                            ),
                            "allelic_requirement": _clean(
                                row.get("inheritance")
                            ),
                            "confidence": _clean(
                                row.get("mechanism_confidence")
                            ),
                            "mechanism_confidence": _clean(
                                row.get("mechanism_confidence")
                            ),
                            "disease_confidence": _clean(
                                row.get("disease_confidence")
                            ),
                        }
                        identity = (
                            source,
                            record["source_record_id"],
                            record["source_condition_id"],
                            record["mondo_id"],
                            mechanism,
                            record["allelic_requirement"],
                        )
                        for key in keys:
                            if identity in seen_by_symbol[key]:
                                continue
                            seen_by_symbol[key].add(identity)
                            by_symbol[key].append(record)
        except Exception as exc:
            logger.warning(
                "Failed to load condition mechanism evidence table %s: %s",
                self.ddg2p_evidence,
                exc,
            )
            self._condition_mechanism_by_symbol = {}
            return self._condition_mechanism_by_symbol

        self._condition_mechanism_by_symbol = {
            symbol: records for symbol, records in by_symbol.items() if records
        }
        return self._condition_mechanism_by_symbol

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
                "source_condition_id": _clean(row.get("source_condition_id")),
                "mondo_id": _clean(row.get("mondo_id")),
                "disease_scope": _clean(row.get("disease_scope")),
                "priva_scope": _clean(row.get("priva_scope")),
                "scope_review_status": _clean(row.get("scope_review_status")),
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
        assertion_columns = {
            "gene_symbol",
            "disease_id",
            "hpo_id",
            "frequency",
            "evidence",
            "reference",
        }
        if assertion_columns.issubset(df.columns):
            raw_condition_index = build_hpo_condition_index(df)
            condition_index: dict[tuple[str, str], dict[str, Any]] = {}
            for (symbol, condition_id), record in raw_condition_index.items():
                condition_index[(symbol, condition_id)] = record
                resolved = self._try_resolve_symbol(symbol)
                if resolved:
                    condition_index.setdefault(
                        (resolved.upper(), condition_id),
                        record,
                    )
            self._hpo_condition_index = condition_index
            hpo_records = []
            for symbol, group in df.groupby("gene_symbol", sort=False):
                if "priva_scope" in group.columns:
                    included_group = group.loc[group["priva_scope"].eq("include")]
                    review_diseases = list(
                        dict.fromkeys(
                            filter(
                                None,
                                map(
                                    _clean,
                                    group.loc[
                                        group["priva_scope"].eq("review"), "disease_id"
                                    ],
                                ),
                            )
                        )
                    )
                    excluded_diseases = list(
                        dict.fromkeys(
                            filter(
                                None,
                                map(
                                    _clean,
                                    group.loc[
                                        group["priva_scope"].eq("exclude"), "disease_id"
                                    ],
                                ),
                            )
                        )
                    )
                else:
                    included_group = group
                    review_diseases = []
                    excluded_diseases = []
                hpo_ids = list(
                    dict.fromkeys(filter(None, map(_clean, included_group["hpo_id"])))
                )
                inheritance_rows = included_group[
                    included_group["hpo_id"].isin(HPO_INHERITANCE_TERMS)
                ]
                inheritance_modes = list(
                    dict.fromkeys(
                        HPO_INHERITANCE_TERMS[hpo_id]
                        for hpo_id in inheritance_rows["hpo_id"]
                    )
                )
                inheritance_diseases = list(
                    dict.fromkeys(filter(None, map(_clean, inheritance_rows["disease_id"])))
                )
                hpo_records.append(
                    (
                        symbol,
                        {
                            "ncbi_gene_id": "",
                            "hpo_id": ";".join(hpo_ids),
                            "inheritance_modes": ";".join(inheritance_modes),
                            "inheritance_disease_ids": ";".join(inheritance_diseases),
                            "scope_review_required": bool(review_diseases),
                            "scope_review_disease_ids": ";".join(review_diseases),
                            "scope_excluded_disease_ids": ";".join(excluded_diseases),
                        },
                    )
                )
        else:
            self._hpo_condition_index = {}
            hpo_records = []
            for _, row in df.iterrows():
                symbol = _clean(row.get("gene_symbol"))
                if not symbol:
                    continue
                hpo_records.append(
                    (
                        symbol,
                        {
                            "ncbi_gene_id": _clean(row.get("ncbi_gene_id")),
                            "hpo_id": _clean(row.get("hpo_id")),
                            "inheritance_modes": _clean(row.get("inheritance_modes")),
                            "inheritance_disease_ids": _clean(
                                row.get("inheritance_disease_ids")
                            ),
                        },
                    )
                )

        for symbol, record in hpo_records:
            symbol = _clean(symbol)
            if not symbol:
                continue
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
                "genomic_location": _clean(row.get("Genomic Location")),
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

    # ------------------------------------------------------------------
    # Exact sequence-variant mechanism evidence
    # ------------------------------------------------------------------

    def _exact_matcher(self) -> ExactVariantMatcher:
        if self._exact_variant_matcher is None:
            self._exact_variant_matcher = ExactVariantMatcher(
                mechanism_json=self.mechanism_json,
                load_mechanisms=self._load_mechanisms,
                try_resolve_symbol=self._try_resolve_symbol,
                resolved_symbol_key=self._resolved_symbol_key,
            )
        return self._exact_variant_matcher

    def _load_canonical_exact_nonlof_rows(self) -> list[dict[str, Any]]:
        return self._exact_matcher()._load_canonical_exact_nonlof_rows()

    def _load_gofcards_variant_hgvsp(
        self,
        *,
        gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    ) -> dict[tuple[str, str], list[dict[str, str]]]:
        return self._exact_matcher()._load_gofcards_variant_hgvsp(
            gofcards_exact_hgvsp_tsv=gofcards_exact_hgvsp_tsv
        )

    def _load_gofcards_variant_hgvsc(
        self,
    ) -> dict[tuple[str, str], list[dict[str, str]]]:
        return self._exact_matcher()._load_gofcards_variant_hgvsc()

    def _load_gofcards_variant_genomic(
        self,
        *,
        gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    ) -> dict[tuple[str, str, str, str], list[dict[str, str]]]:
        return self._exact_matcher()._load_gofcards_variant_genomic(
            gofcards_exact_hgvsp_tsv=gofcards_exact_hgvsp_tsv
        )

    @staticmethod
    def _gofcards_match_entry(
        variant: dict[str, Any],
        *,
        assembly: str = "",
        transcript: str = "",
        hgvsp: str = "",
        hgvsc: str = "",
    ) -> dict[str, str]:
        return ExactVariantMatcher._gofcards_match_entry(
            variant,
            assembly=assembly,
            transcript=transcript,
            hgvsp=hgvsp,
            hgvsc=hgvsc,
        )

    @staticmethod
    def _deduplicate_gofcards_matches(
        rows: list[dict[str, str]],
    ) -> list[dict[str, str]]:
        return ExactVariantMatcher._deduplicate_gofcards_matches(rows)

    def match_curated_nonlof_variant(
        self,
        gene_symbol: Any,
        *,
        hgvsp: Any = "",
        hgvsc: Any = "",
        chrom: Any = "",
        pos: Any = "",
        ref: Any = "",
        alt: Any = "",
        assembly: str = "hg38",
        key_type: str = "auto",
    ) -> dict[str, Any]:
        return self._exact_matcher().match_curated_nonlof_variant(
            gene_symbol,
            hgvsp=hgvsp,
            hgvsc=hgvsc,
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            assembly=assembly,
            key_type=key_type,
        )

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
        return self._exact_matcher().match_gofcards_variant_gof(
            gene_symbol,
            hgvsp,
            gofcards_exact_hgvsp_tsv=gofcards_exact_hgvsp_tsv,
            gofcards_step1_tsv=gofcards_step1_tsv,
            gofcards_active_tsv=gofcards_active_tsv,
            gofcards_raw_xlsx=gofcards_raw_xlsx,
            gofcards_conversion_audit_tsv=gofcards_conversion_audit_tsv,
        )

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
        return self._exact_matcher().match_gofcards_variant_gof_by_genomic(
            gene_symbol,
            chrom,
            pos,
            ref,
            alt,
            assembly=assembly,
            key_type=key_type,
            gofcards_exact_hgvsp_tsv=gofcards_exact_hgvsp_tsv,
        )

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
                            "source_record_id": _clean(data.get("source_record_id")),
                            "source_condition_id": _clean(
                                data.get("source_condition_id")
                            ),
                            "mondo_id": _clean(data.get("mondo_id")),
                            "disease": _clean(data.get("disease")),
                            "disease_scope": _clean(data.get("disease_scope")),
                            "priva_scope": _clean(data.get("priva_scope")),
                            "scope_review_status": _clean(
                                data.get("scope_review_status")
                            ),
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
        """Return condition-resolved and unresolved integrated-cache history."""
        symbol = self._resolved_symbol_key(gene_symbol)
        gene = self._load_condition_cache().get(symbol, {})
        entries = condition_cache_mechanism_entries(gene)
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
        entries = [
            entry
            for entry in condition_cache_mechanism_entries(
                self._load_condition_cache().get(symbol, {})
            )
            if entry.get("level") == "gene_level"
            and entry.get("source") == "G2P_DDG2P"
            and entry.get("mechanism") == "LOF"
            and _clean(entry.get("mechanism_confidence")).lower()
            in DDG2P_USABLE_MECHANISM_CONFIDENCE
            and _clean(entry.get("disease_confidence")).lower()
            in DDG2P_USABLE_DISEASE_CONFIDENCE
        ]
        inheritance_counts = Counter(
            _clean(entry.get("allelic_requirement")) for entry in entries
        )
        recessive = any(
            inheritance in DDG2P_RECESSIVE_LOF_INHERITANCE
            for inheritance in inheritance_counts
        )
        dominant = any(
            inheritance in DDG2P_DOMINANT_LOF_INHERITANCE
            for inheritance in inheritance_counts
        )
        x_linked_dominant = any(
            inheritance in DDG2P_X_LINKED_DOMINANT_LOF_INHERITANCE
            for inheritance in inheritance_counts
        )
        x_linked_unspecified = any(
            inheritance in DDG2P_X_LINKED_UNSPECIFIED_LOF_INHERITANCE
            for inheritance in inheritance_counts
        )
        y_linked = any(
            inheritance in DDG2P_Y_LINKED_LOF_INHERITANCE
            for inheritance in inheritance_counts
        )
        out: dict[str, Any] = {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "has_ddg2p_lof_history": bool(entries),
            "ddg2p_lof_recessive": recessive,
            "ddg2p_lof_dominant": dominant,
            "ddg2p_lof_x_linked_dominant": x_linked_dominant,
            "ddg2p_lof_x_linked_unspecified": x_linked_unspecified,
            "ddg2p_lof_y_linked": y_linked,
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
        """Return coarse inheritance/HI calls plus condition-cache histories."""
        symbol = self._resolved_symbol_key(gene_symbol)
        hpo_record = self._load_hpo().get(symbol, {})
        clingen_record = self._load_clingen().get(symbol, {})
        loeuf = self._load_loeuf().get(symbol, np.nan)
        clingen_hi_score = clingen_record.get("haploinsufficiency_score")
        ddg2p_lof = self.ddg2p_lof_history(symbol)
        condition_assertions = condition_cache_mechanism_assertions(
            self._load_condition_cache().get(symbol, {})
        )
        condition_inheritance: set[str] = set()
        condition_lof_inheritance: set[str] = set()
        for assertion in condition_assertions:
            normalized = {
                inheritance
                for inheritance, _ in normalize_inheritance_histories(
                    assertion.get("allelic_requirement", ""),
                    assertion.get("hpo_inheritance_modes", ()),
                )
                if inheritance
            }
            condition_inheritance.update(normalized)
            if _clean(assertion.get("mechanism")).upper() == "LOF":
                condition_lof_inheritance.update(normalized)

        row_dict = {
            "Gene": symbol,
            "SYMBOL": symbol,
            "LOEUF": loeuf,
            "HPO_IDs": hpo_record.get("hpo_id", ""),
            "HPO_gene_inheritance": hpo_record.get("inheritance_modes") or None,
        }
        hpo_inheritance = parse_hpo_inheritance(row_dict)
        result = identify_inheritance_mode_per_row(
            row_dict,
            gene_mean_am_score,
            clingen_hi_score,
            chromosome=clingen_record.get("genomic_location", ""),
        )
        recessive, dominant, non_monogenic, non_mendelian, haplo_insufficient, incomplete = result
        hpo_inheritance = hpo_inheritance if isinstance(hpo_inheritance, dict) else {}
        effective_recessive = bool(recessive) or "recessive" in condition_inheritance
        effective_dominant = bool(dominant) or "dominant" in condition_inheritance
        effective_haplo = bool(haplo_insufficient) or bool(
            condition_lof_inheritance & {"dominant", "x_linked_dominant"}
        )
        return {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "recessive": effective_recessive,
            "dominant": effective_dominant,
            "x_linked_recessive": bool(
                hpo_inheritance.get("hpo_x_linked_recessive")
            ) or "x_linked_recessive" in condition_inheritance,
            "x_linked_dominant": bool(
                hpo_inheritance.get("hpo_x_linked_dominant")
            ) or "x_linked_dominant" in condition_inheritance,
            "x_linked_unspecified": bool(
                hpo_inheritance.get("hpo_x_linked_unspecified")
            ) or "x_linked_unspecified" in condition_inheritance,
            "y_linked": "y_linked" in condition_inheritance,
            "non_monogenic": bool(non_monogenic) or bool(
                condition_inheritance & {"digenic", "oligogenic", "polygenic"}
            ),
            "non_mendelian": bool(non_mendelian)
            or "non_mendelian" in condition_inheritance,
            "haplo_insufficient": effective_haplo,
            "incomplete_penetrance": bool(incomplete),
            "hpo_inheritance_modes": hpo_record.get("inheritance_modes", ""),
            "hpo_inheritance_disease_ids": hpo_record.get("inheritance_disease_ids", ""),
            "hpo_scope_review_required": bool(
                hpo_record.get("scope_review_required", False)
            ),
            "hpo_scope_review_disease_ids": hpo_record.get(
                "scope_review_disease_ids", ""
            ),
            "hpo_scope_excluded_disease_ids": hpo_record.get(
                "scope_excluded_disease_ids", ""
            ),
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
            "decision_function": "acmg_inheritance.identify_inheritance_mode_per_row",
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
        """Return automatic condition histories from the integrated cache."""
        symbol = self._resolved_symbol_key(gene_symbol)
        gene = self._load_condition_cache().get(symbol, {})
        return condition_cache_mechanism_assertions(gene)
