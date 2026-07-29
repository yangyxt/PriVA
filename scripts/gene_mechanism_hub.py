#!/usr/bin/env python3
"""Gene-condition and query-variant mechanism hub for PriVA.

This module normalizes a gene query to one current HGNC symbol, then reports:

1. Condition-resolved mechanism history from the integrated HPO cache.
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
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from clinvar_vcv import open_text

SCRIPT_DIR = Path(__file__).resolve().parent
PRIVA_ROOT = SCRIPT_DIR.parent
DATA_DIR = PRIVA_ROOT / "data"

DEFAULT_HGNC_TABLE = DATA_DIR / "hgnc" / "non_alt_loci_set.tsv"
DEFAULT_HPO_ASSERTIONS = DATA_DIR / "hpo" / "genes_to_phenotype.assertions.tsv.gz"
# Keep the existing constructor argument name while migrating its file schema.
DEFAULT_HPO_COLLAPSED = DEFAULT_HPO_ASSERTIONS
DEFAULT_CLINGEN_DOSAGE = DATA_DIR / "clingen" / "gene_dosage_sensitivity.hg19.tsv"
DEFAULT_LOEUF_TABLE = DATA_DIR / "loeuf" / "loeuf_dataset.tsv.gz"
DEFAULT_HPO_CONDITION_MECHANISM_CACHE = (
    DATA_DIR
    / "gene_pathogenic_mechanism"
    / "prepared"
    / "hpo_condition_mechanism_cache.json.gz"
)
HPO_CONDITION_MECHANISM_SCHEMA_VERSION = "1.0"
PACKAGED_MECHANISM_JSON = (
    DATA_DIR / "gene_pathogenic_mechanism" / "prepared" / "gene_mechanism_curated_assertions.json"
)
CANONICAL_NONLOF_MECHANISM_JSON = (
    DATA_DIR
    / "gene_pathogenic_mechanism"
    / "prepared"
    / "gene_nonlof_mechanism_curated_assertions.json.gz"
)
DEFAULT_MECHANISM_JSON = (
    CANONICAL_NONLOF_MECHANISM_JSON
    if CANONICAL_NONLOF_MECHANISM_JSON.exists()
    else PACKAGED_MECHANISM_JSON
)
DEFAULT_DDG2P_MECHANISM_EVIDENCE = (
    DATA_DIR / "gene_pathogenic_mechanism" / "prepared" / "gene_pathogenic_mechanism_evidence.tsv"
)
DEFAULT_GOFCARDS_EXACT_GOF_HGVSP = DATA_DIR / "gofcards" / "gofcards_exact_gof.json.gz"
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
EXACT_SEQUENCE_MECHANISMS = {"GOF", "DOMINANT_NEGATIVE"}
VARIANT_MECHANISM_SCORE_KEYS = ("LOF", "GOF", "DOMINANT_NEGATIVE")
CONDITION_MECHANISM_SOURCES = {"G2P_DDG2P", "Orphadata"}
CONDITION_MECHANISM_EVIDENCE_COLUMNS = {
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
    "patho_mode_raw",
    "normalized_mechanisms",
    "mechanism_confidence",
    "disease_confidence",
    "pmids",
}
HPO_INHERITANCE_TERMS = {
    "HP:0000006": "Autosomal dominant inheritance",
    "HP:0000007": "Autosomal recessive inheritance",
    "HP:0001417": "X-linked inheritance",
    "HP:0001419": "X-linked recessive inheritance",
    "HP:0001423": "X-linked dominant inheritance",
    "HP:0001426": "Non-Mendelian inheritance",
    "HP:0001427": "Mitochondrial inheritance",
    "HP:0001450": "Y-linked inheritance",
    "HP:0010982": "Polygenic inheritance",
    "HP:0010983": "Oligogenic inheritance",
    "HP:0010984": "Digenic inheritance",
    "HP:0012275": "Autosomal dominant inheritance with maternal imprinting",
    "HP:0034341": "Pseudoautosomal recessive inheritance",
}
HPO_CACHE_INHERITANCE_LABELS = {
    "autosomal_dominant": "Autosomal dominant inheritance",
    "autosomal_recessive": "Autosomal recessive inheritance",
    "x_linked": "X-linked inheritance",
    "x_linked_recessive": "X-linked recessive inheritance",
    "x_linked_dominant": "X-linked dominant inheritance",
    "non_mendelian": "Non-Mendelian inheritance",
    "mitochondrial": "Mitochondrial inheritance",
    "y_linked": "Y-linked inheritance",
    "polygenic": "Polygenic inheritance",
    "oligogenic": "Oligogenic inheritance",
    "digenic": "Digenic inheritance",
    "autosomal_dominant_maternal_imprinting": (
        "Autosomal dominant inheritance with maternal imprinting"
    ),
    "pseudoautosomal_recessive": "Pseudoautosomal recessive inheritance",
}
HPO_PENETRANCE_TERMS = {
    "HP:0003829",  # Incomplete penetrance
    "HP:0034950",  # Complete penetrance
    "HP:4000159",  # Moderate penetrance
    "HP:4000160",  # Low penetrance
}
HPO_ONSET_TERMS = {
    "HP:0030674",  # Antenatal onset
    "HP:0011460",  # Embryonal onset
    "HP:0011461",  # Fetal onset
    "HP:0003577",  # Congenital onset
    "HP:0003623",  # Neonatal onset
    "HP:0003593",  # Infantile onset
    "HP:0011463",  # Childhood onset
    "HP:0003621",  # Juvenile onset
    "HP:0410280",  # Pediatric onset
    "HP:0003581",  # Adult onset
    "HP:0011462",  # Young adult onset
    "HP:0003596",  # Middle age onset
    "HP:0003584",  # Late onset
}
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
    "variant_effect_suppressed_evidence",
    "variant_effect_conflict",
    "variant_lof_score",
    "variant_gof_score",
    "variant_dn_score",
    "variant_mechanism_exclusive",
    "variant_exact_mechanisms",
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


def _bounded_mechanism_score(value: Any) -> int:
    try:
        score = int(float(_clean(value)))
    except ValueError:
        return 0
    return score if score in {0, 1, 2} else 0


def _exact_mechanisms_from_effect(value: Any) -> set[str]:
    effect = _clean(value)
    prefix = "exact_known_"
    if not effect.startswith(prefix):
        return set()
    return {
        _normalize_sequence_mechanism(token)
        for token in effect[len(prefix):].split("+")
        if _normalize_sequence_mechanism(token) in EXACT_SEQUENCE_MECHANISMS
    }


def infer_query_variant_effect(row: dict[str, Any] | pd.Series) -> dict[str, Any]:
    """Infer one exclusive variant effect while retaining suppressed predictions.

    Exact curated GoF or dominant-negative evidence has score ``2`` and defines
    the allowed mechanism set for this allele. Generic consequence, NMD, or
    LOFTEE signals remain visible as audit evidence but cannot create a competing
    LoF call. Without an exact assertion, prediction-based LoF evidence receives
    score ``1``. LOFTEE ``OS`` is its other-splice category.
    """
    predicted_lof_evidence: list[str] = []
    loftee_tokens = _split_annotation_tokens(row.get("LoF"))
    if "HC" in loftee_tokens:
        predicted_lof_evidence.append("LOFTEE_HC")
    if "OS" in loftee_tokens:
        predicted_lof_evidence.append("LOFTEE_OS")

    consequence = _clean(row.get("Consequence")).lower()
    nmd = _clean(row.get("NMD")).lower()
    lof_filter = _clean(row.get("LoF_filter")).upper()
    truncating = "stop_gained" in consequence or "frameshift" in consequence
    nmd_lof = truncating and "escaping" not in nmd and "END_TRUNC" not in lof_filter
    if nmd_lof:
        predicted_lof_evidence.append("NMD_PREDICTED_LOF")
    if _bool_value(row.get("vep_consq_lof")):
        predicted_lof_evidence.append("VEP_LOF")
    if _bool_value(row.get("splicing_lof")):
        predicted_lof_evidence.append("PRIVA_SPLICE_LOF")
    if _bool_value(row.get("5UTR_lof")):
        predicted_lof_evidence.append("PRIVA_5UTR_LOF")

    scores = {
        "LOF": _bounded_mechanism_score(row.get("variant_lof_score")),
        "GOF": _bounded_mechanism_score(row.get("variant_gof_score")),
        "DOMINANT_NEGATIVE": _bounded_mechanism_score(
            row.get("variant_dn_score")
        ),
    }
    gof_tokens = _split_annotation_tokens(row.get("variant_gof_tag"))
    if "GOF" in gof_tokens:
        scores["GOF"] = 2

    exact_mechanisms = {
        mechanism for mechanism, score in scores.items() if score == 2
    } & EXACT_SEQUENCE_MECHANISMS
    predicted_lof_evidence = list(dict.fromkeys(predicted_lof_evidence))
    suppressed_evidence: list[str] = []
    evidence: list[str] = []
    if exact_mechanisms:
        for mechanism in scores:
            scores[mechanism] = 2 if mechanism in exact_mechanisms else 0
        evidence.extend(
            f"CANONICAL_EXACT_{mechanism}" for mechanism in sorted(exact_mechanisms)
        )
        suppressed_evidence = predicted_lof_evidence
        effect = "exact_known_" + "+".join(sorted(exact_mechanisms))
    elif predicted_lof_evidence:
        scores["LOF"] = max(scores["LOF"], 1)
        evidence.extend(predicted_lof_evidence)
        effect = "predicted_LOF_high_confidence"
    else:
        effect = "uncertain"
    return {
        "variant_effect": effect,
        "variant_effect_evidence": ";".join(evidence),
        "variant_effect_suppressed_evidence": ";".join(suppressed_evidence),
        "variant_effect_conflict": "",
        "variant_lof_score": scores["LOF"],
        "variant_gof_score": scores["GOF"],
        "variant_dn_score": scores["DOMINANT_NEGATIVE"],
        "variant_mechanism_exclusive": bool(exact_mechanisms),
        "variant_exact_mechanisms": ";".join(sorted(exact_mechanisms)),
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


def _normalize_assembly(value: Any) -> str:
    aliases = {
        "grch37": "hg19",
        "hg19": "hg19",
        "b37": "hg19",
        "grch38": "hg38",
        "hg38": "hg38",
        "b38": "hg38",
    }
    cleaned = _clean(value).lower()
    return aliases.get(cleaned, cleaned)


def _requested_genomic_key_types(value: Any) -> tuple[str, ...]:
    cleaned = _clean(value).lower() or "auto"
    if cleaned == "auto":
        return ("vcf", "genomic")
    if cleaned in {"vcf", "genomic"}:
        return (cleaned,)
    return ()


def _norm_hgvsp(value: Any) -> str:
    text = _clean(value)
    if not text:
        return ""
    protein = text.split(":")[-1].strip()
    if protein.startswith("p."):
        protein = protein[2:]
    return protein.upper()


def _normalize_sequence_mechanism(value: Any) -> str:
    """Normalize an explicitly curated sequence-variant mechanism label."""
    token = re.sub(r"[^A-Z0-9]+", "_", _clean(value).upper()).strip("_")
    aliases = {
        "GAIN_OF_FUNCTION": "GOF",
        "ACTIVATING": "GOF",
        "DN": "DOMINANT_NEGATIVE",
        "DOMINANTNEGATIVE": "DOMINANT_NEGATIVE",
    }
    return aliases.get(token, token)


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


def build_hpo_condition_index(
    assertions: pd.DataFrame,
) -> dict[tuple[str, str], dict[str, Any]]:
    """Index HPO evidence by the inseparable gene-plus-disease identity.

    The index deliberately requires both a gene symbol and a stable disease ID.
    A MONDO ID alone is insufficient because one disease can involve multiple
    genes with different mechanisms. A gene symbol alone is also insufficient
    because the same gene can cause several disorders with different
    inheritance, penetrance, or onset.

    Every original HPO assertion is retained with its frequency, evidence code,
    and reference. Compact inheritance, penetrance, and onset lists are derived
    only for convenient downstream decisions; they never replace the auditable
    source rows. Each condition can be looked up by its native OMIM/ORPHA ID or
    by MONDO, allowing G2P and Orphadata to use the same condition record.
    """
    required = {
        "gene_symbol",
        "disease_id",
        "hpo_id",
        "frequency",
        "evidence",
        "reference",
    }
    missing = sorted(required - set(assertions.columns))
    if missing:
        raise ValueError(f"HPO assertion table missing columns: {', '.join(missing)}")

    frame = assertions.copy().fillna("")
    for column in (
        "mondo_id",
        "disease_scope",
        "priva_scope",
        "scope_review_status",
    ):
        if column not in frame.columns:
            frame[column] = ""
    frame = frame.loc[
        frame["gene_symbol"].astype(str).str.strip().ne("")
        & frame["disease_id"].astype(str).str.strip().ne("")
    ]

    index: dict[tuple[str, str], dict[str, Any]] = {}
    group_columns = [
        "gene_symbol",
        "disease_id",
        "mondo_id",
        "disease_scope",
        "priva_scope",
        "scope_review_status",
    ]
    for group_key, group in frame.groupby(group_columns, sort=False, dropna=False):
        (
            gene_symbol,
            disease_id,
            mondo_id,
            disease_scope,
            priva_scope,
            scope_review_status,
        ) = map(_clean, group_key)
        hpo_rows = [
            {
                "hpo_id": _clean(row.get("hpo_id")),
                "frequency": _clean(row.get("frequency")),
                "evidence": _clean(row.get("evidence")),
                "reference": _clean(row.get("reference")),
            }
            for row in group.to_dict(orient="records")
        ]
        hpo_ids = list(
            dict.fromkeys(row["hpo_id"] for row in hpo_rows if row["hpo_id"])
        )
        record = {
            "disease_id": disease_id,
            "mondo_id": mondo_id,
            "disease_scope": disease_scope,
            "priva_scope": priva_scope,
            "scope_review_status": scope_review_status,
            "hpo_ids": hpo_ids,
            "inheritance_modes": [
                HPO_INHERITANCE_TERMS[hpo_id]
                for hpo_id in hpo_ids
                if hpo_id in HPO_INHERITANCE_TERMS
            ],
            "penetrance_hpo_ids": [
                hpo_id for hpo_id in hpo_ids if hpo_id in HPO_PENETRANCE_TERMS
            ],
            "onset_hpo_ids": [
                hpo_id for hpo_id in hpo_ids if hpo_id in HPO_ONSET_TERMS
            ],
            "hpo_assertions": hpo_rows,
        }
        for condition_id in (disease_id, mondo_id):
            if condition_id:
                index[(gene_symbol.upper(), condition_id.upper())] = record
    return index


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
        self._canonical_exact_nonlof_rows: list[dict[str, Any]] | None = None
        self._gofcards_by_symbol_hgvsp: dict[tuple[str, str], list[dict[str, str]]] | None = None
        self._gofcards_by_symbol_hgvsc: dict[tuple[str, str], list[dict[str, str]]] | None = None
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

    def _load_canonical_exact_nonlof_rows(self) -> list[dict[str, Any]]:
        """Return runtime-eligible GoFCards variant blocks from the canonical JSON.

        ``build_gene_nonlof_mechanism_cache.py`` embeds each variant under
        ``exact_normalized_variants`` exactly as the GoFCards cache stores it,
        so the block still has its ``record`` and its per-assembly ``genomic``
        and ``transcripts``. Runtime reads that structure and never reinterprets
        the upstream cache independently.

        ClinVar VCV entries link an already curated variant to clinical
        assertions; ClinVar pathogenicity alone establishes neither gain of
        function nor a dominant-negative effect, so those links are excluded
        here on purpose.
        """
        if self._canonical_exact_nonlof_rows is not None:
            return self._canonical_exact_nonlof_rows

        entries: list[dict[str, Any]] = []
        seen: set[str] = set()
        for gene in self._load_mechanisms().values():
            symbol = self._try_resolve_symbol(gene.get("symbol"))
            if not symbol:
                continue
            for wrapped in gene.get("variant_level", []) or []:
                if not isinstance(wrapped, dict):
                    continue
                for source, assertion in wrapped.items():
                    if source == "ClinVar_VCV" or not isinstance(assertion, dict):
                        continue
                    status = assertion.get("exact_normalization_status")
                    if status and status != "matched_gene_concordant":
                        continue
                    for variant in assertion.get("exact_normalized_variants", []) or []:
                        if not isinstance(variant, dict):
                            continue
                        record = variant.get("record") or {}
                        if (record.get("eligibility") or {}).get("status") != "eligible":
                            continue
                        variant_id = _clean(variant.get("variant_id"))
                        # One allele can be asserted by two curated sources with
                        # different mechanisms, so the identity that de-duplicates
                        # must include the source. Keying on the allele alone
                        # would silently drop the second mechanism.
                        marker = f"{source}|{symbol}|{variant_id}"
                        if not variant_id or marker in seen:
                            continue
                        seen.add(marker)
                        entries.append({
                            **variant,
                            "symbol": symbol,
                            # Carried from the wrapping assertion so a curated
                            # dominant-negative source is not reported as GoFCards
                            # gain of function.
                            "assertion_source": source,
                            "assertion_mechanism": _clean(assertion.get("mechanism")),
                            # Which canonical file this evidence came from, so an
                            # audit row identifies its own provenance.
                            "canonical_mechanism_json": str(self.mechanism_json),
                        })

        self._canonical_exact_nonlof_rows = entries
        return entries

    def _load_gofcards_variant_hgvsp(
        self,
        *,
        gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    ) -> dict[tuple[str, str], list[dict[str, str]]]:
        """Index eligible variants by (HGNC symbol, normalized protein change).

        A variant registers under every protein change its transcripts present,
        because one variant is numbered differently on different isoforms --
        JAK2 V617F is also Val468Phe and Val212Phe. Indexing only one of them
        would miss a query annotated on any other transcript.
        """
        if self._gofcards_by_symbol_hgvsp is not None:
            return self._gofcards_by_symbol_hgvsp

        index: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for variant in self._load_canonical_exact_nonlof_rows():
            symbol = variant["symbol"]
            for assembly, block in (variant.get("assemblies") or {}).items():
                for transcript, view in (block.get("transcripts") or {}).items():
                    for hgvsp, coding in (view.get("by_hgvsp") or {}).items():
                        key = _norm_hgvsp(hgvsp)
                        if not key:
                            continue
                        index[(symbol, key)].append(
                            self._gofcards_match_entry(
                                variant, assembly=assembly, transcript=transcript,
                                hgvsp=hgvsp, hgvsc=coding[0] if coding else "",
                            )
                        )

        self._gofcards_by_symbol_hgvsp = {
            key: self._deduplicate_gofcards_matches(rows) for key, rows in index.items()
        }
        return self._gofcards_by_symbol_hgvsp

    def _load_gofcards_variant_hgvsc(
        self,
    ) -> dict[tuple[str, str], list[dict[str, str]]]:
        """Index eligible variants by (HGNC symbol, coding change).

        Measured unique across the catalogue at 25,106 keys with no collisions,
        so a coding change identifies a variant within its gene where a protein
        change may not.
        """
        if self._gofcards_by_symbol_hgvsc is not None:
            return self._gofcards_by_symbol_hgvsc

        index: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for variant in self._load_canonical_exact_nonlof_rows():
            symbol = variant["symbol"]
            for assembly, block in (variant.get("assemblies") or {}).items():
                for transcript, view in (block.get("transcripts") or {}).items():
                    for hgvsc, detail in (view.get("by_hgvsc") or {}).items():
                        key = _clean(hgvsc).upper()
                        if not key:
                            continue
                        index[(symbol, key)].append(
                            self._gofcards_match_entry(
                                variant, assembly=assembly, transcript=transcript,
                                hgvsp=detail.get("hgvsp") or "", hgvsc=hgvsc,
                            )
                        )

        self._gofcards_by_symbol_hgvsc = {
            key: self._deduplicate_gofcards_matches(rows) for key, rows in index.items()
        }
        return self._gofcards_by_symbol_hgvsc

    @staticmethod
    def _gofcards_match_entry(
        variant: dict[str, Any], *, assembly: str = "", transcript: str = "",
        hgvsp: str = "", hgvsc: str = "",
    ) -> dict[str, str]:
        """Flatten one variant into the shape the matcher reports back."""
        record = variant.get("record") or {}
        source = record.get("source") or {}
        genomic = ((variant.get("assemblies") or {}).get(assembly) or {}).get("genomic") or {}
        # GoFCards is the only variant-level source in the shipped cache, so it
        # is the default; a curated source that declares its own mechanism keeps
        # that mechanism rather than being reported as gain of function.
        assertion_source = _clean(variant.get("assertion_source")) or "GoFCards"
        assertion_mechanism = _clean(variant.get("assertion_mechanism")) or "GOF"
        return {
            "source": assertion_source,
            "mechanism": assertion_mechanism,
            "canonical_assertion_source": assertion_source,
            "canonical_mechanism_json": _clean(variant.get("canonical_mechanism_json")),
            "symbol": variant.get("symbol", ""),
            "gofcards_variant_id": _clean(variant.get("variant_id")),
            "gofcards_allele_key": _clean(source.get("gofcards_allele_key")),
            "gofcards_accession_id": _clean(
                (record.get("annotations") or {}).get("clinvar_variation_id")),
            "assembly": assembly,
            "vep_transcript": transcript,
            "HGVSp": hgvsp,
            "HGVSc": hgvsc,
            "chrom": _clean(genomic.get("chrom")),
            "pos": _clean(genomic.get("pos")),
            "ref": _clean(genomic.get("ref")),
            "alt": _clean(genomic.get("alt")),
            "match_eligibility": _clean((record.get("eligibility") or {}).get("status")),
            "pmids": ";".join(
                e.get("pmid", "") for e in (record.get("evidence") or []) if e.get("pmid")
            ),
        }

    @staticmethod
    def _deduplicate_gofcards_matches(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        seen: set[tuple[str, ...]] = set()
        unique_rows: list[dict[str, str]] = []
        for row in sorted(
            rows,
            key=lambda item: (
                item.get("gofcards_variant_id", ""),
                item.get("assembly", ""),
                item.get("vep_transcript", ""),
                item.get("HGVSp", ""),
            ),
        ):
            # The variant identifier already distinguishes variants, so a match
            # reported twice through different transcripts of the same variant
            # is one match, not two.
            dedup_key = (
                row.get("mechanism", ""),
                row.get("canonical_assertion_source", ""),
                row.get("gofcards_variant_id", ""),
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
        """Index eligible variants by (symbol, assembly, key type, chrom|pos|ref|alt).

        The coordinate is read from the assembly block it belongs to, so a
        variant that failed liftover simply contributes no hg38 key rather than
        an empty one. ``key_type`` is retained for the existing call contract;
        the cache stores a single normalized representation per build, so both
        types resolve to the same key.
        """
        if self._gofcards_by_symbol_genomic is not None:
            return self._gofcards_by_symbol_genomic

        index: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
        for variant in self._load_canonical_exact_nonlof_rows():
            symbol = variant["symbol"]
            for assembly, block in (variant.get("assemblies") or {}).items():
                coords = block.get("genomic") or {}
                if not coords.get("chrom") or coords.get("pos") is None:
                    continue
                key = _genomic_match_key(
                    coords.get("chrom"), coords.get("pos"),
                    coords.get("ref"), coords.get("alt"),
                )
                if not key:
                    continue
                entry = self._gofcards_match_entry(variant, assembly=assembly)
                for key_type in ("vcf", "genomic"):
                    index[(symbol, assembly, key_type, key)].append(entry)

        self._gofcards_by_symbol_genomic = {
            key: self._deduplicate_gofcards_matches(rows) for key, rows in index.items()
        }
        return self._gofcards_by_symbol_genomic

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
        """Score exact curated GoF/DN evidence for one query allele.

        Scores use the agreed variant-level contract: ``0`` means no exact
        support, and ``2`` means an exact curated assertion that excludes
        unasserted mechanisms. Score ``1`` is reserved for prediction-based
        plausibility and is added later by :func:`infer_query_variant_effect`;
        this exact matcher does not infer loss of function from ClinVar
        pathogenicity or consequence alone.

        ``current HGNC gene + normalized HGVSp`` is an exact match route.
        Normalized genomic alleles are the fallback when HGVSp is absent or not
        found. More than one mechanism receives score ``2`` only if separate
        exact records explicitly assert those mechanisms for the same allele.
        """
        symbol = self._resolved_symbol_key(gene_symbol)
        hgvsp_key = _norm_hgvsp(hgvsp)
        hgvsc_key = _clean(hgvsc).split(":")[-1].upper()
        assembly_key = _normalize_assembly(assembly)
        genomic_key = _genomic_match_key(chrom, pos, ref, alt)
        scores = {mechanism: 0 for mechanism in VARIANT_MECHANISM_SCORE_KEYS}
        matches: list[dict[str, str]] = []
        route = ""
        matched_key_type = ""

        if symbol and hgvsp_key:
            matches = list(
                self._load_gofcards_variant_hgvsp().get(
                    (symbol, hgvsp_key),
                    [],
                )
            )
            if matches:
                route = "hgvsp"

        # Coding change is the second route: it separates two different coding
        # changes that produce the same protein change.
        if symbol and not matches and hgvsc_key:
            matches = list(
                self._load_gofcards_variant_hgvsc().get((symbol, hgvsc_key), [])
            )
            if matches:
                route = "hgvsc"

        if symbol and not matches and genomic_key and assembly_key:
            genomic_lookup = self._load_gofcards_variant_genomic()
            for candidate_key_type in _requested_genomic_key_types(key_type):
                candidate = genomic_lookup.get(
                    (
                        symbol,
                        assembly_key,
                        candidate_key_type,
                        genomic_key,
                    ),
                    [],
                )
                if candidate:
                    matches.extend(candidate)
                    route = "genomic"
                    matched_key_type = candidate_key_type
            matches = self._deduplicate_gofcards_matches(matches)

        matches = [
            match
            for match in matches
            if _normalize_sequence_mechanism(match.get("mechanism"))
            in EXACT_SEQUENCE_MECHANISMS
        ]
        mechanisms = sorted(
            {
                _normalize_sequence_mechanism(match.get("mechanism"))
                for match in matches
            }
        )
        for mechanism in mechanisms:
            scores[mechanism] = 2

        accession_ids = sorted(
            {
                _clean(match.get("gofcards_accession_id"))
                for match in matches
                if _clean(match.get("gofcards_accession_id"))
            }
        )
        variant_ids = sorted(
            {
                _clean(match.get("gofcards_variant_id"))
                for match in matches
                if _clean(match.get("gofcards_variant_id"))
            }
        )
        return {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "input_hgvsp": _clean(hgvsp),
            "matched_hgvsp_key": hgvsp_key if route == "hgvsp" else "",
            "matched_hgvsc_key": hgvsc_key if route == "hgvsc" else "",
            "input_assembly": _clean(assembly),
            "assembly": assembly_key,
            "input_genomic_key": genomic_key,
            "matched_key_type": matched_key_type,
            "match_route": route,
            "mechanism_scores": scores,
            "lof_score": scores["LOF"],
            "gof_score": scores["GOF"],
            "dn_score": scores["DOMINANT_NEGATIVE"],
            "mechanism_exclusive": bool(mechanisms),
            "exclusive_mechanisms": mechanisms,
            "matched_mechanisms": mechanisms,
            "gofcards_accession_id": ";".join(accession_ids),
            "gofcards_variant_id": ";".join(variant_ids),
            "matches": matches,
        }

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
        """Compatibility view of canonical exact GoF HGVSp evidence."""
        result = self.match_curated_nonlof_variant(gene_symbol, hgvsp=hgvsp)
        gof_matches = [
            match
            for match in result["matches"]
            if _normalize_sequence_mechanism(match.get("mechanism")) == "GOF"
        ]
        result["matches"] = gof_matches
        result["variant_gof_tag"] = "GOF" if gof_matches else ""
        return result

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
        """Compatibility view of canonical exact GoF genomic evidence."""
        result = self.match_curated_nonlof_variant(
            gene_symbol,
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            assembly=assembly,
            key_type=key_type,
        )
        gof_matches = [
            match
            for match in result["matches"]
            if _normalize_sequence_mechanism(match.get("mechanism")) == "GOF"
        ]
        result.update(
            {
                "input_chrom": _clean(chrom),
                "input_pos": _clean(pos),
                "input_ref": _clean(ref),
                "input_alt": _clean(alt),
                "matched_genomic_key": result["input_genomic_key"],
                "variant_gof_tag": "GOF" if gof_matches else "",
                "matches": gof_matches,
            }
        )
        return result

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
        """Return automatic condition histories from the integrated cache."""
        symbol = self._resolved_symbol_key(gene_symbol)
        gene = self._load_condition_cache().get(symbol, {})
        return condition_cache_mechanism_assertions(gene)

    def matched_clinvar_vcv_for_gofcards(
        self,
        gene_symbol: Any,
        *,
        gofcards_variant_ids: Any = "",
        gofcards_accession_ids: Any = "",
    ) -> list[dict[str, Any]]:
        """Return integrated-cache variants linked to exact GoFCards IDs.

        The historical method name is retained for callers, but the full
        curated mechanism JSON is no longer read. Each result identifies
        whether the exact variant is nested under a condition or retained in
        ``unmapped_evidence`` for audit only.
        """
        variant_ids = {
            _norm(token)
            for token in re.split(r"[;,]", _clean(gofcards_variant_ids))
            if _clean(token)
        }
        accession_ids = {
            _norm(token)
            for token in re.split(r"[;,]", _clean(gofcards_accession_ids))
            if _clean(token)
        }
        if not variant_ids and not accession_ids:
            return []

        symbol = self._resolved_symbol_key(gene_symbol)
        gene = self._load_condition_cache().get(symbol, {})
        return match_condition_cache_gofcards_variants(
            gene,
            variant_ids=variant_ids,
            accession_ids=accession_ids,
        )


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


def condition_cache_context(
    condition_key: Any,
    condition: dict[str, Any],
) -> dict[str, Any]:
    """Return the shared PriVA context for one included cache condition.

    An empty result means the condition cannot influence automatic ACMG logic.
    The cache retains only the HPO axes used here; ``hpo_assertion_count`` keeps
    the size of the complete phenotype record visible without loading it.
    """
    scope = condition.get("priva_scope", {})
    if not isinstance(scope, dict) or _clean(scope.get("decision")).lower() != "include":
        return {}

    identifiers = condition.get("identifiers", {})
    if not isinstance(identifiers, dict):
        identifiers = {}
    mondo_ids = identifiers.get("MONDO", [])
    mondo_id = _clean(mondo_ids[0]) if isinstance(mondo_ids, list) and mondo_ids else ""

    inheritance = condition.get("inheritance", {})
    inheritance = inheritance if isinstance(inheritance, dict) else {}
    hpo_inheritance_modes = [
        HPO_CACHE_INHERITANCE_LABELS.get(_clean(mode), _clean(mode))
        for mode in inheritance.get("modes", [])
        if _clean(mode)
    ]
    penetrance = condition.get("penetrance", {})
    penetrance = penetrance if isinstance(penetrance, dict) else {}
    onset = condition.get("onset", {})
    onset = onset if isinstance(onset, dict) else {}

    axis_assertions: list[dict[str, str]] = []
    axis_seen: set[str] = set()
    for axis in (inheritance, penetrance, onset):
        for raw_assertion in axis.get("assertions", []) or []:
            if not isinstance(raw_assertion, dict):
                continue
            hpo_assertion = {
                "hpo_id": _clean(raw_assertion.get("hpo_id")),
                "frequency": _clean(raw_assertion.get("frequency")),
                "evidence": _clean(raw_assertion.get("evidence")),
                "reference": _clean(raw_assertion.get("reference")),
            }
            key = json.dumps(hpo_assertion, sort_keys=True, separators=(",", ":"))
            if key not in axis_seen:
                axis_seen.add(key)
                axis_assertions.append(hpo_assertion)

    return {
        "hpo_match_status": "matched_gene_and_condition",
        "hpo_disease_id": _clean(condition_key),
        "mondo_id": mondo_id,
        "disease_scope": _clean(scope.get("category")),
        "priva_scope": "include",
        "scope_review_status": _clean(scope.get("review_status")),
        "hpo_inheritance_modes": hpo_inheritance_modes,
        "penetrance_hpo_ids": [
            _clean(item.get("hpo_id"))
            for item in penetrance.get("assertions", []) or []
            if isinstance(item, dict) and _clean(item.get("hpo_id"))
        ],
        "onset_hpo_ids": [
            _clean(item.get("hpo_id"))
            for item in onset.get("assertions", []) or []
            if isinstance(item, dict) and _clean(item.get("hpo_id"))
        ],
        "hpo_assertions": axis_assertions,
        "hpo_assertion_count": condition.get("hpo_assertion_count", 0),
    }


def condition_cache_mechanism_entries(
    gene_record: dict[str, Any],
) -> list[dict[str, Any]]:
    """Flatten one integrated gene record for hub summaries and audits.

    Unlike automatic assertion selection, this audit view retains mechanisms
    from review/excluded conditions and evidence that could not be joined to an
    HPO condition. Scope and condition-link status remain explicit on every
    row. Exact GoFCards alleles are emitted at ``variant_level`` and never
    converted into gene-level mechanism evidence.
    """
    entries: list[dict[str, Any]] = []
    seen: set[str] = set()

    def append_entry(entry: dict[str, Any]) -> None:
        normalized = {
            **entry,
            "level": _clean(entry.get("level")),
            "source": _clean(entry.get("source")),
            "mechanism": _clean(entry.get("mechanism")).upper(),
            "pmids": list(entry.get("pmids", []) or []),
        }
        if normalized["mechanism"] not in CANONICAL_MECHANISMS:
            return
        identity = json.dumps(normalized, sort_keys=True, separators=(",", ":"))
        if identity not in seen:
            seen.add(identity)
            entries.append(normalized)

    conditions = gene_record.get("conditions", {})
    if isinstance(conditions, dict):
        for condition_key, condition in sorted(conditions.items()):
            if not isinstance(condition, dict):
                continue
            scope = condition.get("priva_scope", {})
            scope = scope if isinstance(scope, dict) else {}
            condition_common = {
                "source_condition_id": _clean(condition_key),
                "disease": _clean(condition.get("label")),
                "disease_scope": _clean(scope.get("category")),
                "priva_scope": _clean(scope.get("decision")),
                "scope_review_status": _clean(scope.get("review_status")),
                "condition_link_status": "exact",
            }
            mechanisms = condition.get("pathogenic_mechanisms", {})
            if not isinstance(mechanisms, dict):
                continue
            for mechanism, block in sorted(mechanisms.items()):
                mechanism = _clean(mechanism).upper()
                if mechanism not in CANONICAL_MECHANISMS or not isinstance(block, dict):
                    continue
                for evidence in block.get("evidence", []) or []:
                    if not isinstance(evidence, dict):
                        continue
                    source = _clean(evidence.get("source"))
                    if source == "GoFCards_exact+ClinVar_VCV":
                        continue
                    append_entry(
                        {
                            **condition_common,
                            "level": "gene_level",
                            "source": source,
                            "source_record_id": _clean(
                                evidence.get("source_record_id")
                            ),
                            "mechanism": mechanism,
                            "mechanism_raw": _clean(evidence.get("mechanism_raw")),
                            "allelic_requirement": _clean(
                                evidence.get("allelic_requirement")
                            ),
                            "confidence": _clean(
                                evidence.get("mechanism_confidence")
                            ),
                            "mechanism_confidence": _clean(
                                evidence.get("mechanism_confidence")
                            ),
                            "disease_confidence": _clean(
                                evidence.get("disease_confidence")
                            ),
                            "pmids": list(evidence.get("pmids", []) or []),
                        }
                    )
                for variant_key, variant in sorted(
                    (block.get("variants", {}) or {}).items()
                ):
                    if not isinstance(variant, dict):
                        continue
                    links = [
                        link
                        for link in variant.get("clinvar_links", []) or []
                        if isinstance(link, dict)
                    ]
                    append_entry(
                        {
                            **condition_common,
                            "level": "variant_level",
                            "source": "GoFCards_exact+ClinVar_VCV",
                            "source_record_id": _clean(variant_key),
                            "mechanism": mechanism,
                            "allelic_requirement": "",
                            "confidence": "exact_variant",
                            "mechanism_confidence": "exact_variant",
                            "disease_confidence": "ClinVar_germline_assertion",
                            "pmids": list(variant.get("pmids", []) or []),
                            "vcv_accessions": [
                                _clean(link.get("vcv_accession"))
                                for link in links
                                if _clean(link.get("vcv_accession"))
                            ],
                        }
                    )

    unmapped = gene_record.get("unmapped_evidence", {})
    unmapped = unmapped if isinstance(unmapped, dict) else {}
    for evidence in unmapped.get("mechanisms", []) or []:
        if not isinstance(evidence, dict):
            continue
        source_scope = evidence.get("source_scope", {})
        source_scope = source_scope if isinstance(source_scope, dict) else {}
        append_entry(
            {
                "level": "gene_level",
                "source": _clean(evidence.get("source")),
                "source_record_id": _clean(evidence.get("source_record_id")),
                "source_condition_id": next(
                    iter(evidence.get("condition_identifiers", []) or []), ""
                ),
                "disease": _clean(evidence.get("condition_label")),
                "mechanism": _clean(evidence.get("mechanism")),
                "mechanism_raw": _clean(evidence.get("mechanism_raw")),
                "allelic_requirement": _clean(evidence.get("allelic_requirement")),
                "confidence": _clean(evidence.get("mechanism_confidence")),
                "mechanism_confidence": _clean(
                    evidence.get("mechanism_confidence")
                ),
                "disease_confidence": _clean(evidence.get("disease_confidence")),
                "disease_scope": _clean(source_scope.get("category")),
                "priva_scope": _clean(source_scope.get("decision")),
                "scope_review_status": _clean(source_scope.get("review_status")),
                "condition_link_status": "unresolved",
                "unmapped_reason": _clean(evidence.get("reason")),
                "pmids": list(evidence.get("pmids", []) or []),
            }
        )
    variants = unmapped.get("variants", {})
    variants = variants if isinstance(variants, dict) else {}
    for variant_key, variant in sorted(variants.items()):
        if not isinstance(variant, dict):
            continue
        append_entry(
            {
                "level": "variant_level",
                "source": "GoFCards",
                "source_record_id": _clean(variant_key),
                "source_condition_id": "",
                # The disease name comes from the ClinVar records linked to
                # this variant. GoFCards' own free-text disease field is not
                # carried into the condition cache: it has no identifier, so
                # it can name a condition but never join to one.
                "disease": ";".join(
                    dict.fromkeys(
                        name
                        for link in variant.get("clinvar_links", []) or []
                        if isinstance(link, dict)
                        for name in link.get("condition_names", []) or []
                        if _clean(name)
                    )
                ),
                "mechanism": _clean(variant.get("mechanism")),
                "allelic_requirement": "",
                "confidence": "exact_variant_condition_unresolved",
                "mechanism_confidence": "exact_variant",
                "disease_confidence": "",
                "condition_link_status": "unresolved",
                "unmapped_reason": _clean(
                    (variant.get("condition_link", {}) or {}).get("reason")
                ),
                "pmids": list(variant.get("pmids", []) or []),
            }
        )
    return entries


def match_condition_cache_gofcards_variants(
    gene_record: dict[str, Any],
    *,
    variant_ids: set[str],
    accession_ids: set[str],
) -> list[dict[str, Any]]:
    """Find exact GoFCards query tokens in one integrated gene record.

    Variant identifiers look like ``loc_12:21995260:C->T_grch37`` and contain
    the separators used between records only by accident, so callers split on
    semicolon and comma alone. Tokens are compared after case normalization.
    Condition-linked and unresolved cache locations are returned in one shape,
    with an empty condition for unresolved variants.

    ``accession_ids`` is accepted and ignored. The condition cache stores one
    variant identifier per variant, minted the same way the matcher mints it,
    so that identifier settles the match on its own; the ClinVar accession was
    a second name for the same variant.
    """
    wanted_variant_ids = {_norm(value) for value in variant_ids if _clean(value)}
    if not wanted_variant_ids:
        return []

    matches: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()

    def add_matching_variants(
        variants: Any,
        *,
        condition_key: str = "",
        condition: dict[str, Any] | None = None,
    ) -> None:
        if not isinstance(variants, dict):
            return
        for variant_key, variant in variants.items():
            if not isinstance(variant, dict):
                continue
            # One identifier per variant, and it is the same identifier the
            # matcher reports, so it alone decides the match. The accession
            # route it replaced was a secondary key for the same variant.
            cached_variant_id = _norm(variant.get("gofcards_variant_id"))
            if not cached_variant_id or cached_variant_id not in wanted_variant_ids:
                continue
            identity = (_clean(condition_key), _clean(variant_key))
            if identity in seen:
                continue
            seen.add(identity)
            matches.append(
                {
                    "condition_key": _clean(condition_key),
                    "condition": condition if isinstance(condition, dict) else {},
                    "variant_key": _clean(variant_key),
                    "variant": variant,
                }
            )

    conditions = gene_record.get("conditions", {})
    if isinstance(conditions, dict):
        for condition_key, condition in sorted(conditions.items()):
            if not isinstance(condition, dict):
                continue
            mechanisms = condition.get("pathogenic_mechanisms", {})
            mechanisms = mechanisms if isinstance(mechanisms, dict) else {}
            gof = mechanisms.get("GOF", {})
            gof = gof if isinstance(gof, dict) else {}
            add_matching_variants(
                gof.get("variants", {}),
                condition_key=_clean(condition_key),
                condition=condition,
            )

    unmapped = gene_record.get("unmapped_evidence", {})
    unmapped = unmapped if isinstance(unmapped, dict) else {}
    add_matching_variants(unmapped.get("variants", {}))
    return matches


def condition_cache_mechanism_assertions(
    gene_record: dict[str, Any],
) -> list[dict[str, Any]]:
    """Translate one integrated gene record to PriVA's assertion contract.

    The integrated cache has already performed the only permitted join:
    gene + exact condition identifier. This adapter therefore does not repeat
    disease matching. It exposes only conditions whose effective PriVA scope is
    ``include`` and retains the condition's inheritance, penetrance, onset, and
    assertion provenance beside every mechanism record.

    The cache also stores synthetic ``GoFCards_exact+ClinVar_VCV`` evidence in a
    condition block when a particular allele has an exact ClinVar condition
    link. That source is intentionally excluded here: it becomes an assertion
    only when the query allele matches the nested cache variant. Otherwise one
    exact GOF allele would incorrectly create gene-wide GOF history.
    """
    assertions: list[dict[str, Any]] = []
    conditions = gene_record.get("conditions", {})
    if not isinstance(conditions, dict):
        return assertions

    for condition_key, condition in sorted(conditions.items()):
        if not isinstance(condition, dict):
            continue
        common = condition_cache_context(condition_key, condition)
        if not common:
            continue

        mechanism_blocks = condition.get("pathogenic_mechanisms", {})
        if not isinstance(mechanism_blocks, dict):
            continue
        for mechanism, block in sorted(mechanism_blocks.items()):
            mechanism = _clean(mechanism).upper()
            if mechanism not in CANONICAL_MECHANISMS or not isinstance(block, dict):
                continue
            for evidence in block.get("evidence", []) or []:
                if not isinstance(evidence, dict):
                    continue
                source = _clean(evidence.get("source"))
                if source not in CONDITION_MECHANISM_SOURCES:
                    continue
                evidence_identifiers = evidence.get("condition_identifiers", [])
                evidence_identifiers = (
                    evidence_identifiers if isinstance(evidence_identifiers, list) else []
                )
                source_condition_id = next(
                    (
                        _clean(identifier)
                        for identifier in evidence_identifiers
                        if _clean(identifier) and not _clean(identifier).upper().startswith("MONDO:")
                    ),
                    _clean(condition_key),
                )
                requirement = _clean(evidence.get("allelic_requirement"))
                requirements = [requirement] if requirement else sorted(
                    _hpo_allelic_requirements(
                        ";".join(common["hpo_inheritance_modes"])
                    )
                )
                if not requirements:
                    requirements = [""]

                for effective_requirement in requirements:
                    assertions.append(
                        {
                            **common,
                            "source": source,
                            "source_record_id": _clean(evidence.get("source_record_id")),
                            "source_condition_id": source_condition_id,
                            "disease": _clean(evidence.get("condition_label"))
                            or _clean(condition.get("label")),
                            "mechanism": mechanism,
                            "mechanism_raw": _clean(evidence.get("mechanism_raw")),
                            "allelic_requirement": effective_requirement,
                            "confidence": _clean(evidence.get("mechanism_confidence")),
                            "mechanism_confidence": _clean(
                                evidence.get("mechanism_confidence")
                            ),
                            "disease_confidence": _clean(
                                evidence.get("disease_confidence")
                            ),
                            "pmids": list(evidence.get("pmids", []) or []),
                        }
                    )
    return _deduplicate_assertions(assertions)


def enrich_condition_mechanism_assertion(
    assertion: dict[str, Any],
    *,
    gene_symbol: Any,
    hpo_condition_index: dict[tuple[str, str], dict[str, Any]],
) -> list[dict[str, Any]]:
    """Bind one mechanism assertion to HPO for the same gene and disease.

    Native disease identity is tried before MONDO identity. This avoids linking
    through names and prevents another disease of the same gene from donating
    inheritance, penetrance, or onset. Automatic transfer requires an effective
    ``priva_scope=include``; review, excluded, and still-unscoped diseases remain
    audit-only.

    G2P allelic requirements remain authoritative when present. Orphadata does
    not encode allelic requirements, so an unambiguous HPO inheritance term can
    supply a compact dominant, recessive, or mitochondrial requirement. If HPO
    legitimately records more than one mode, separate records are returned so
    no mode is hidden in a lossy combined string.
    """
    record = dict(assertion)
    source_scope = _clean(record.get("priva_scope")).lower()
    if source_scope in {"review", "exclude"}:
        return []

    symbol_key = _clean(gene_symbol).upper()
    hpo_record: dict[str, Any] | None = None
    for condition_id in (
        record.get("source_condition_id"),
        record.get("mondo_id"),
    ):
        condition_key = _clean(condition_id).upper()
        if condition_key:
            hpo_record = hpo_condition_index.get((symbol_key, condition_key))
        if hpo_record is not None:
            break

    if hpo_record is None:
        # A source assertion without a registry-supported inherited-disease
        # scope is useful for manual review, but it cannot safely influence an
        # automated germline classification. Explicitly included G2P diseases
        # remain usable even when HPO has no phenotype rows for that condition.
        if source_scope != "include":
            return []
        record.update(
            {
                "hpo_match_status": "no_matching_gene_condition_hpo_record",
                "hpo_inheritance_modes": [],
                "penetrance_hpo_ids": [],
                "onset_hpo_ids": [],
                "hpo_assertions": [],
            }
        )
        return [record]
    if _clean(hpo_record.get("priva_scope")).lower() in {"review", "exclude"}:
        return []

    record.update(
        {
            "hpo_match_status": "matched_gene_and_condition",
            "hpo_disease_id": _clean(hpo_record.get("disease_id")),
            "mondo_id": _clean(record.get("mondo_id"))
            or _clean(hpo_record.get("mondo_id")),
            "hpo_inheritance_modes": list(
                hpo_record.get("inheritance_modes", [])
            ),
            "penetrance_hpo_ids": list(
                hpo_record.get("penetrance_hpo_ids", [])
            ),
            "onset_hpo_ids": list(hpo_record.get("onset_hpo_ids", [])),
            "hpo_assertions": list(hpo_record.get("hpo_assertions", [])),
            "disease_scope": _clean(hpo_record.get("disease_scope"))
            or _clean(record.get("disease_scope")),
            "priva_scope": _clean(hpo_record.get("priva_scope"))
            or _clean(record.get("priva_scope")),
            "scope_review_status": _clean(
                hpo_record.get("scope_review_status")
            )
            or _clean(record.get("scope_review_status")),
        }
    )
    if _clean(record.get("priva_scope")).lower() != "include":
        return []
    if _clean(record.get("allelic_requirement")):
        return [record]

    requirements = sorted(
        _hpo_allelic_requirements(";".join(record["hpo_inheritance_modes"]))
    )
    if not requirements:
        return [record]
    return [
        {**record, "allelic_requirement": requirement}
        for requirement in requirements
    ]


def extract_exact_clinvar_condition_identities(
    condition_assertion: dict[str, Any],
) -> list[dict[str, str]]:
    """Extract only explicit disease identifiers from one matched ClinVar RCV.

    The aggregate condition commonly carries a MedGen identifier, while its
    contributing SCVs can carry exact OMIM, Orphanet, or MONDO cross-references
    in ``trait_mappings``. These database identifiers are safe for an identity
    join. Disease names, publications, and phenotype similarity are deliberately
    ignored because they do not prove that two condition records are identical.

    MedGen identifiers are retained for audit, but the current HPO condition
    table is keyed by OMIM/ORPHA/MONDO. Consequently, MedGen alone does not
    transfer HPO inheritance or penetrance.
    """
    raw_identities: list[tuple[str, str, str, str]] = []
    for condition in condition_assertion.get("conditions", []) or []:
        if not isinstance(condition, dict):
            continue
        raw_identities.append(
            (
                _clean(condition.get("database")),
                _clean(condition.get("id")),
                _clean(condition.get("name")),
                "ClinVar_RCV_condition",
            )
        )
    for scv in condition_assertion.get("matched_scvs", []) or []:
        if not isinstance(scv, dict):
            continue
        for mapping in scv.get("trait_mappings", []) or []:
            if not isinstance(mapping, dict):
                continue
            raw_identities.append(
                (
                    _clean(mapping.get("mapping_ref")),
                    _clean(mapping.get("mapping_value")),
                    _clean(mapping.get("medgen_name")),
                    "ClinVar_SCV_trait_mapping",
                )
            )

    prefixes = {
        "OMIM": "OMIM",
        "ORPHA": "ORPHA",
        "ORPHANET": "ORPHA",
        "MONDO": "MONDO",
        "MEDGEN": "MEDGEN",
    }
    identities: list[dict[str, str]] = []
    seen: set[str] = set()
    for database, identifier, name, provenance in raw_identities:
        prefix = prefixes.get(database.upper())
        if not prefix or not identifier:
            continue
        cleaned_id = identifier
        for known_prefix in ("OMIM:", "ORPHA:", "ORPHANET:", "MONDO:", "MEDGEN:"):
            if cleaned_id.upper().startswith(known_prefix):
                cleaned_id = cleaned_id.split(":", 1)[1]
                break
        condition_id = f"{prefix}:{cleaned_id}"
        if condition_id.upper() in seen:
            continue
        seen.add(condition_id.upper())
        identities.append(
            {
                "source_condition_id": condition_id,
                "mondo_id": condition_id if prefix == "MONDO" else "",
                "disease": name,
                "condition_id_provenance": provenance,
            }
        )
    return identities


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


def select_condition_histories_for_variant(
    assertions: list[dict[str, Any]],
    *,
    variant_effect: str,
    variant_effect_conflict: str = "",
) -> list[dict[str, Any]]:
    """Select only condition histories compatible with the query allele.

    A curated mechanism elsewhere in a gene is background history, not evidence
    that every variant acts through that mechanism. Selection therefore occurs
    before inheritance or penetrance can influence ACMG criteria:

    * an exact known GOF allele selects GOF histories;
    * a high-confidence predicted LOF allele selects LOF histories; and
    * an unresolved allele retains all histories, clearly marked as uncertain
      by the later applicability classifier.

    Exact curated mechanisms are exclusive. A consequence-based LoF prediction
    is retained in the row's suppressed-evidence field but cannot reintroduce a
    different condition history after an exact GoF or DN match.
    """
    effect = _clean(variant_effect)
    exact_mechanisms = _exact_mechanisms_from_effect(effect)
    if exact_mechanisms:
        allowed = exact_mechanisms
    elif effect == "predicted_LOF_high_confidence":
        allowed = {"LOF"}
    else:
        allowed = CANONICAL_MECHANISMS
    return [
        assertion
        for assertion in assertions
        if _clean(assertion.get("mechanism")).upper() in allowed
    ]


def _deduplicate_assertions(assertions: list[dict[str, Any]]) -> list[dict[str, Any]]:
    fields = (
        "source",
        "source_record_id",
        "source_condition_id",
        "mondo_id",
        "disease",
        "mechanism",
        "allelic_requirement",
        "confidence",
    )
    seen: set[str] = set()
    output: list[dict[str, Any]] = []
    for assertion in assertions:
        normalized = dict(assertion)
        normalized.update(
            {
                "source": _clean(assertion.get("source")),
                "source_record_id": _clean(assertion.get("source_record_id")),
                "source_condition_id": _clean(
                    assertion.get("source_condition_id")
                ),
                "mondo_id": _clean(assertion.get("mondo_id")),
                "disease": _clean(assertion.get("disease")),
                "mechanism": _clean(assertion.get("mechanism")).upper()
                or "UNRESOLVED",
                "allelic_requirement": _clean(
                    assertion.get("allelic_requirement")
                ),
                "confidence": _clean(assertion.get("confidence")),
            }
        )
        key = json.dumps(
            {field: normalized.get(field, "") for field in fields},
            sort_keys=True,
            separators=(",", ":"),
        )
        if key not in seen:
            seen.add(key)
            output.append(normalized)
    return output


def _classify_variant_applicability(
    assertions: list[dict[str, Any]],
    effect: str,
    effect_conflict: str = "",
) -> dict[str, Any]:
    exact_mechanisms = _exact_mechanisms_from_effect(effect)
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
                status, reason = (
                    "incompatible",
                    "exact_nonLOF_mechanism_excludes_predicted_LOF",
                )
        elif mechanism == "GOF":
            if "GOF" in exact_mechanisms:
                status, reason = "applicable", "exact_query_GOF_match"
            elif effect == "uncertain":
                status, reason = "uncertain", "query_GOF_effect_not_established"
            else:
                status, reason = "incompatible", "GOF_requires_exact_variant_match"
        elif mechanism == "DOMINANT_NEGATIVE":
            if "DOMINANT_NEGATIVE" in exact_mechanisms:
                status, reason = "applicable", "exact_query_DN_match"
            elif effect == "uncertain":
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


def summarize_condition_cache_exact_gof_matches(
    matches: list[dict[str, Any]],
) -> dict[str, Any]:
    """Build audit outputs and scoped assertions for exact cache variants.

    Audit fields include both condition-linked and unresolved matches. Automatic
    GOF assertions are more restrictive: the variant must be nested under a
    condition whose integrated PriVA scope is ``include``. This separation
    retains useful ClinVar provenance without letting unresolved, complex, or
    excluded disease links alter germline ACMG criteria.
    """
    output: dict[str, Any] = {
        "vcv_accessions": [],
        "condition_names": [],
        "review_stars": [],
        "hgvs": [],
        "assertions": [],
    }

    def append_unique(field: str, value: Any) -> None:
        cleaned = _clean(value)
        if cleaned and cleaned not in output[field]:
            output[field].append(cleaned)

    for match in matches:
        if not isinstance(match, dict):
            continue
        variant = match.get("variant", {})
        variant = variant if isinstance(variant, dict) else {}
        condition = match.get("condition", {})
        condition = condition if isinstance(condition, dict) else {}
        condition_key = _clean(match.get("condition_key"))
        condition_label = _clean(condition.get("label"))
        links = variant.get("clinvar_links", [])
        links = [link for link in links or [] if isinstance(link, dict)]
        has_link_condition_name = False
        has_link_hgvs = False

        for link in links:
            append_unique("vcv_accessions", link.get("vcv_accession"))
            for name in link.get("condition_names", []) or []:
                has_link_condition_name = has_link_condition_name or bool(_clean(name))
                append_unique("condition_names", name)
            for expression in link.get("hgvs", []) or []:
                has_link_hgvs = has_link_hgvs or bool(_clean(expression))
                append_unique("hgvs", expression)
            stars = link.get("review_stars")
            try:
                star_value = int(stars)
            except (TypeError, ValueError):
                pass
            else:
                if star_value not in output["review_stars"]:
                    output["review_stars"].append(star_value)

        if not has_link_hgvs:
            # The coding and protein changes travel per transcript view, each
            # carrying the transcript version they belong to.
            for view in variant.get("transcripts", []) or []:
                if not isinstance(view, dict):
                    continue
                append_unique("hgvs", view.get("hgvsc"))
                append_unique("hgvs", view.get("hgvsp"))
        if not has_link_condition_name and condition_label:
            append_unique("condition_names", condition_label)

        context = condition_cache_context(condition_key, condition)
        if not context:
            continue
        requirements = sorted(
            _hpo_allelic_requirements(
                ";".join(context.get("hpo_inheritance_modes", []))
            )
        ) or [""]
        assertion_links = links or [{}]
        for link in assertion_links:
            stars = link.get("review_stars")
            confidence = (
                f"ClinVar_{stars}_star"
                if stars is not None and _clean(stars)
                else "exact_variant_match"
            )
            for requirement in requirements:
                output["assertions"].append(
                    {
                        **context,
                        "source": "GoFCards_exact+ClinVar_VCV",
                        "source_record_id": _clean(link.get("vcv_accession"))
                        or _clean(match.get("variant_key")),
                        "source_condition_id": condition_key,
                        "disease": condition_label,
                        "mechanism": "GOF",
                        "mechanism_raw": "gain of function",
                        "allelic_requirement": requirement,
                        "confidence": confidence,
                        "mechanism_confidence": "exact_variant",
                        "disease_confidence": "ClinVar_germline_assertion",
                        "clinical_significance": _clean(
                            link.get("clinical_significance")
                        ),
                        "condition_identifiers": list(
                            link.get("condition_identifiers", []) or []
                        ),
                        "pmids": list(variant.get("pmids", []) or []),
                    }
                )

    output["assertions"] = _deduplicate_assertions(output["assertions"])
    return output


def annotate_gene_mechanism_categories(
    df: pd.DataFrame,
    *,
    clinvar_pathogenic_genes: set[str] | None = None,
    gene_to_am_score_map: dict[str, float] | None = None,
    condition_cache: str | Path = DEFAULT_HPO_CONDITION_MECHANISM_CACHE,
    mechanism_json: str | Path | None = None,
    ddg2p_evidence: str | Path | None = None,
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
    """Annotate condition-specific history and query-variant applicability.

    HPO inheritance, penetrance, onset, G2P/Orphadata mechanisms, and compact
    ClinVar links come from the integrated cache. Row-level HPO text and
    gene-wide ClinVar, LOEUF, AlphaMissense, or ClinGen dosage signals remain
    audit information and cannot create a condition-mechanism assertion.
    ``mechanism_json``, ``ddg2p_evidence``, and ``hpo_inheritance_col`` are
    retained only for call compatibility and are not runtime evidence sources.
    """
    if symbol_col not in df.columns:
        raise KeyError(f"missing symbol column: {symbol_col}")
    if gene_col not in df.columns:
        raise KeyError(f"missing gene column: {gene_col}")

    hub = GeneMechanismHub(
        condition_cache=condition_cache,
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
    assertion_cache: dict[str, list[dict[str, Any]]] = {}
    clinvar_pathogenic_genes = set(clinvar_pathogenic_genes or set())
    gene_to_am_score_map = gene_to_am_score_map or {}
    for _, row in out.iterrows():
        gene = row[symbol_col]
        symbol = hub.resolve_symbol(gene)
        if symbol not in assertion_cache:
            assertion_cache[symbol] = hub.condition_mechanism_assertions(gene)
        assertions = list(assertion_cache[symbol])

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
        effect_call = infer_query_variant_effect(row)
        assertions = select_condition_histories_for_variant(
            assertions,
            variant_effect=effect_call["variant_effect"],
            variant_effect_conflict=effect_call["variant_effect_conflict"],
        )
        vcv_accessions: list[str] = []
        vcv_conditions: list[str] = []
        vcv_hgvs: list[str] = []
        vcv_review_stars: list[int] = []
        if "GOF" in _exact_mechanisms_from_effect(effect_call["variant_effect"]):
            use_vcv_condition_history = not any(
                assertion.get("mechanism") == "GOF" for assertion in assertions
            )
            vcv_matches = hub.matched_clinvar_vcv_for_gofcards(
                symbol,
                gofcards_variant_ids=row.get("gofcards_variant_id", ""),
                gofcards_accession_ids=row.get("gofcards_accession_id", ""),
            )
            exact_gof = summarize_condition_cache_exact_gof_matches(vcv_matches)
            vcv_accessions.extend(exact_gof["vcv_accessions"])
            vcv_conditions.extend(exact_gof["condition_names"])
            vcv_hgvs.extend(exact_gof["hgvs"])
            vcv_review_stars.extend(exact_gof["review_stars"])
            if use_vcv_condition_history:
                assertions.extend(exact_gof["assertions"])
            if use_vcv_condition_history and not any(
                assertion.get("mechanism") == "GOF" for assertion in assertions
            ):
                assertions.append(
                    {
                        "source": "GoFCards",
                        "disease": "",
                        "mechanism": "GOF",
                        "allelic_requirement": "",
                        "confidence": "exact_variant_match_condition_unresolved",
                    }
                )
        assertions = _deduplicate_assertions(assertions)
        applicability = _classify_variant_applicability(
            assertions,
            effect_call["variant_effect"],
            effect_call["variant_effect_conflict"],
        )

        plausible_mechanism_values.append(applicability["plausible"])
        variant_outputs["gene_lof_evidence"].append(";".join(lof_evidence))
        variant_outputs["variant_effect"].append(effect_call["variant_effect"])
        variant_outputs["variant_effect_evidence"].append(effect_call["variant_effect_evidence"])
        variant_outputs["variant_effect_suppressed_evidence"].append(
            effect_call["variant_effect_suppressed_evidence"]
        )
        variant_outputs["variant_effect_conflict"].append(effect_call["variant_effect_conflict"])
        variant_outputs["variant_lof_score"].append(effect_call["variant_lof_score"])
        variant_outputs["variant_gof_score"].append(effect_call["variant_gof_score"])
        variant_outputs["variant_dn_score"].append(effect_call["variant_dn_score"])
        variant_outputs["variant_mechanism_exclusive"].append(
            effect_call["variant_mechanism_exclusive"]
        )
        variant_outputs["variant_exact_mechanisms"].append(
            effect_call["variant_exact_mechanisms"]
        )
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


def match_curated_nonlof_variant(
    gene_symbol: Any,
    *,
    hgvsp: Any = "",
    chrom: Any = "",
    pos: Any = "",
    ref: Any = "",
    alt: Any = "",
    assembly: str = "hg38",
    key_type: str = "auto",
    mechanism_json: str | Path = DEFAULT_MECHANISM_JSON,
    use_hgnc_package: bool = False,
) -> dict[str, Any]:
    """Convenience API for canonical exact GoF/dominant-negative matching."""
    return GeneMechanismHub(
        mechanism_json=mechanism_json,
        use_hgnc_package=use_hgnc_package,
    ).match_curated_nonlof_variant(
        gene_symbol,
        hgvsp=hgvsp,
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        assembly=assembly,
        key_type=key_type,
    )


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
    parser.add_argument(
        "--condition-cache",
        default=str(DEFAULT_HPO_CONDITION_MECHANISM_CACHE),
        help="Integrated HPO condition-mechanism JSON cache",
    )
    parser.add_argument(
        "--mechanism-json",
        default=str(DEFAULT_MECHANISM_JSON),
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--ddg2p-evidence",
        default=str(DEFAULT_DDG2P_MECHANISM_EVIDENCE),
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--hpo-collapsed", default=str(DEFAULT_HPO_COLLAPSED))
    parser.add_argument("--clingen-dosage", default=str(DEFAULT_CLINGEN_DOSAGE))
    parser.add_argument("--loeuf-table", default=str(DEFAULT_LOEUF_TABLE))
    parser.add_argument("--hgnc-table", default=str(DEFAULT_HGNC_TABLE))
    args = parser.parse_args()

    hub = GeneMechanismHub(
        condition_cache=args.condition_cache,
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
