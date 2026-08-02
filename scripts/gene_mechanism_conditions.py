"""Condition-linked inheritance, penetrance, and mechanism history.

All automatic paths in this module retain the gene-plus-condition identity and
apply the inherited-disease scope gate before returning criterion-facing data.
The module defines logic only; it is not a pipeline entry point.
"""

from __future__ import annotations

import json
import math
from typing import Any

import pandas as pd

from gene_mechanism_common import (
    CANONICAL_MECHANISMS,
    CONDITION_MECHANISM_SOURCES,
    HPO_CACHE_INHERITANCE_LABELS,
    HPO_INHERITANCE_TERMS,
    HPO_ONSET_TERMS,
    HPO_PENETRANCE_TERMS,
    UNRESOLVED_MECHANISM,
    _clean,
    _deduplicate_by,
)
from hpo_penetrance import normalize_penetrance_hpo_ids


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




ASSERTION_IDENTITY_FIELDS = (
    "source",
    "source_record_id",
    "source_condition_id",
    "mondo_id",
    "disease",
    "mechanism",
    "allelic_requirement",
    "confidence",
)




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


def condition_cache_mechanism_assertions(
    gene_record: dict[str, Any],
) -> list[dict[str, Any]]:
    """Translate one integrated gene record to PriVA's assertion contract.

    The integrated cache has already performed the only permitted join:
    gene + exact condition identifier. This adapter therefore does not repeat
    disease matching. It exposes only conditions whose effective PriVA scope is
    ``include`` and retains the condition's inheritance, penetrance, onset, and
    assertion provenance beside every mechanism record.

    An included recessive condition with no curated molecular mechanism may
    carry the cache builder's ``deduced_from_inheritance`` LOF assertion. It is
    exposed here for ordinary variant-to-condition matching, with its
    ``assertion_basis="deduced"`` provenance intact. PVS1 reads a separate
    curated-only view and therefore cannot mistake this inference for evidenced
    LOF history.

    An included condition with neither a curated nor deduced usable mechanism
    emits an ``UNRESOLVED`` assertion instead. This is the gene-history fallback:
    a query variant can still inherit the condition's inheritance and penetrance.
    If HPO states more than one Mendelian inheritance for such a condition, one
    unresolved assertion is emitted per mode so neither is lost.

    The cache also stores synthetic ``GoFCards_exact+ClinVar_VCV`` evidence in a
    condition block when a particular allele has an exact ClinVar condition
    link. It is retained as condition history, but an exact query-variant call
    is still established independently by matching the nested allele record.
    Thus the source says that GOF has occurred for this condition; it does not
    claim that every query variant is an exact GOF allele.
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

        emitted_mechanism = False
        mechanism_blocks = condition.get("pathogenic_mechanisms", {})
        if isinstance(mechanism_blocks, dict):
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
                        evidence_identifiers
                        if isinstance(evidence_identifiers, list)
                        else []
                    )
                    source_condition_id = next(
                        (
                            _clean(identifier)
                            for identifier in evidence_identifiers
                            if _clean(identifier)
                            and not _clean(identifier).upper().startswith("MONDO:")
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
                                "source_record_id": _clean(
                                    evidence.get("source_record_id")
                                ),
                                "source_condition_id": source_condition_id,
                                "disease": _clean(evidence.get("condition_label"))
                                or _clean(condition.get("label")),
                                "mechanism": mechanism,
                                "assertion_basis": _clean(
                                    evidence.get("assertion_basis")
                                ),
                                "mechanism_raw": _clean(
                                    evidence.get("mechanism_raw")
                                ),
                                "allelic_requirement": effective_requirement,
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
                        emitted_mechanism = True

        if emitted_mechanism:
            continue

        # No admitted curated or deduced mechanism was supplied. Keep the
        # condition as an unresolved history, with the full HPO context and
        # provenance in ``common``. Multiple Mendelian modes become separate
        # records; a non-Mendelian or unstated mode stays in the HPO list and is
        # normalized later without inventing an allelic requirement.
        requirements = sorted(
            _hpo_allelic_requirements(";".join(common["hpo_inheritance_modes"]))
        ) or [""]
        for requirement in requirements:
            assertions.append(
                {
                    **common,
                    "source": "HPO_condition_history",
                    "source_record_id": _clean(condition_key),
                    "source_condition_id": _clean(condition_key),
                    "disease": _clean(condition.get("label")),
                    "mechanism": UNRESOLVED_MECHANISM,
                    "assertion_basis": "condition_without_explicit_mechanism",
                    "mechanism_raw": "",
                    "allelic_requirement": requirement,
                    "confidence": "condition_without_explicit_mechanism",
                    "mechanism_confidence": "",
                    "disease_confidence": "",
                    "pmids": [],
                }
            )
    return _deduplicate_by(assertions, ASSERTION_IDENTITY_FIELDS)


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


# The reverse of HPO_CACHE_INHERITANCE_LABELS, so a mode is recognized whether
# it arrives as the cache's key or as the label substituted for display.
_INHERITANCE_LABEL_TO_KEY = {
    label.lower(): key for key, label in HPO_CACHE_INHERITANCE_LABELS.items()
}
# Downstream treats these three the same way: one allele in a carrier is
# enough for the condition to appear.
DOMINANT_LIKE_INHERITANCE = {"dominant", "y_linked", "mitochondrial"}
# Germline, but not single-gene. Delivered rather than discarded, so that
# benign-supporting criteria are not assigned easily against them.
NON_MENDELIAN_INHERITANCE = {
    "non_mendelian",
    "polygenic",
    "digenic",
    "oligogenic",
}


def normalize_inheritance(
    allelic_requirement: Any = "",
    hpo_inheritance_modes: Any = (),
) -> tuple[str, bool]:
    """Reduce the two source vocabularies to one value plus an X-linked flag.

    Two sources describe the same fact in different words: G2P states an
    allelic requirement (``biallelic_autosomal``, ``monoallelic_X_hemizygous``)
    and HPO states an inheritance mode (``autosomal_recessive``,
    ``x_linked_dominant``). Both fold to the same answer.

    The delivered value is ``recessive`` or ``dominant`` wherever that question
    has an answer, because that is what the criteria reason about. Being on the
    X chromosome is a separate fact and is returned separately: folding it into
    the value would erase it, and a hemizygous male affected by one allele is
    still the recessive pattern.

    ``y_linked`` and ``mitochondrial`` are delivered as themselves rather than
    forced into dominant; downstream treats them the same as dominant.
    ``non_mendelian``, ``polygenic``, ``digenic`` and ``oligogenic`` are also
    delivered as themselves. They are not discarded -- they exist so that
    benign-supporting criteria are not assigned easily against them.
    """
    requirement = _clean(allelic_requirement).lower()
    # Modes arrive in either form: the cache's own snake_case key, or the
    # human-readable HPO label that condition_cache_context substitutes. Both
    # name the same mode, so both are accepted here rather than requiring every
    # caller to know which one it is holding.
    modes = [
        _INHERITANCE_LABEL_TO_KEY.get(_clean(mode).lower(), _clean(mode).lower())
        for mode in (
            [hpo_inheritance_modes]
            if isinstance(hpo_inheritance_modes, str)
            else list(hpo_inheritance_modes or [])
        )
        if _clean(mode)
    ]
    x_linked = requirement.startswith("monoallelic_x") or any(
        mode.startswith("x_linked") for mode in modes
    )

    # This compatibility API returns one value, so a Mendelian mode wins when
    # both kinds are supplied together. The variant-condition path calls
    # normalize_inheritance_histories below and preserves every category; that
    # is what keeps the non-Mendelian gate from disappearing.
    mendelian_modes = [mode for mode in modes if mode not in NON_MENDELIAN_INHERITANCE]
    if not mendelian_modes:
        for mode in modes:
            if mode in NON_MENDELIAN_INHERITANCE:
                return mode, x_linked

    if requirement == "recessive" or requirement.startswith("biallelic_") or requirement in {
        "monoallelic_x",
        "monoallelic_x_hemizygous",
    }:
        return "recessive", x_linked
    if requirement == "dominant":
        return "dominant", x_linked
    if requirement == "mitochondrial":
        return "mitochondrial", False
    if requirement == "y_linked":
        return "y_linked", False
    if requirement.startswith("monoallelic_y"):
        return "y_linked", False
    if requirement.startswith("monoallelic_"):
        return "dominant", x_linked

    for mode in mendelian_modes:
        if mode in {"autosomal_recessive", "x_linked_recessive", "pseudoautosomal_recessive"}:
            return "recessive", x_linked
        # A bare "x_linked" with no recessive or dominant qualifier reads as
        # X-linked recessive.
        if mode == "x_linked":
            return "recessive", True
        if mode in {
            "autosomal_dominant",
            "x_linked_dominant",
            "autosomal_dominant_maternal_imprinting",
        }:
            return "dominant", x_linked
        if mode == "mitochondrial":
            return "mitochondrial", False
        if mode == "y_linked":
            return "y_linked", False
    return "", x_linked


def normalize_inheritance_histories(
    allelic_requirement: Any = "",
    hpo_inheritance_modes: Any = (),
) -> list[tuple[str, bool]]:
    """Return every condition-linked inheritance category without collapsing gates.

    ``normalize_inheritance`` remains the single-value compatibility API. The
    variant-condition path needs a stronger contract: if one condition is
    annotated as both autosomal dominant and polygenic, retaining only the
    Mendelian value would silently remove the non-monogenic gate used by the
    criteria. The explicit allelic requirement and each HPO mode are therefore
    normalized independently and de-duplicated.

    This does not invent a molecular mechanism or merge another condition's
    history. It preserves multiple curated assertions already attached to the
    same selected condition. A blank pair is returned only when neither source
    states a usable category, so the condition and its penetrance still travel.
    """
    modes = (
        [hpo_inheritance_modes]
        if isinstance(hpo_inheritance_modes, str)
        else list(hpo_inheritance_modes or [])
    )
    normalized: list[tuple[str, bool]] = []

    if _clean(allelic_requirement):
        normalized.append(normalize_inheritance(allelic_requirement, ()))
    normalized.extend(
        normalize_inheritance("", [mode])
        for mode in modes
        if _clean(mode)
    )

    values = list(
        dict.fromkeys(value for value in normalized if value[0])
    )
    if values:
        return values

    fallback = normalize_inheritance(allelic_requirement, modes)
    return [fallback] if fallback[0] or fallback[1] else [("", False)]


def gene_has_lof_mechanism_history(gene_record: dict[str, Any]) -> bool:
    """Does a curator say this gene causes germline disease by loss of function?

    A property of the GENE, deliberately independent of any one variant. PVS1
    asks whether loss of function is an established disease mechanism for the
    gene at all, and then decides on its own how strong the evidence is for the
    particular null variant in front of it. Filtering this by what the query
    variant does would answer a different question and pre-empt that decision.

    Two filters, and both are deliberate:

    ``priva_scope.decision == "include"`` -- only germline-inherited disease
    counts, the same gate every other reader of the cache applies.

    ``assertion_basis == "curated"`` -- only a curator's molecular-mechanism
    claim counts. A recessive condition may carry an LOF mechanism with
    ``assertion_basis == "deduced"`` for normal variant-condition matching, but
    that inference only restates the inheritance and cannot become circular
    PVS1 evidence. Only curated G2P/DDG2P, Orphadata and ClinGen
    haploinsufficiency assertions reach here.
    """
    for condition in (gene_record.get("conditions") or {}).values():
        if (condition.get("priva_scope") or {}).get("decision") != "include":
            continue
        lof = (condition.get("pathogenic_mechanisms") or {}).get("LOF")
        if not lof:
            continue
        for evidence in lof.get("evidence") or []:
            if (evidence.get("assertion_basis") or "") == "curated":
                return True
    return False


def gene_inheritance_consensus(gene_record: dict[str, Any]) -> tuple[str, bool, int]:
    """Return the one inheritance every germline condition of a gene agrees on.

    Current cache records normally do not need this fallback: every included
    condition emits a curated mechanism, a deduced recessive LOF mechanism, or
    an ``UNRESOLVED`` history upstream. It is retained defensively for malformed
    or legacy cache records that fail to produce any assertion.

    Unanimity is what makes it safe. A gene carrying both a dominant and a
    recessive condition is deliberately given nothing here, because there the
    question genuinely has two answers and only the per-history match can choose
    between them. Measured on the current cache: 2,772 genes are unanimously
    recessive, 1,453 unanimously dominant, and 501 are mixed.

    Conditions that are not germline-inherited disease never contribute, so a
    review-scoped or excluded condition cannot donate its inheritance here.

    The third return value is how many distinct inheritances the gene's germline
    conditions stated, which separates the two ways this can come back empty:
    ``0`` means nothing was stated at all and the caller may fall back further,
    while ``2`` or more means the gene genuinely disagrees with itself and no
    fallback is legitimate.
    """
    values: set[tuple[str, bool]] = set()
    for condition in (gene_record.get("conditions") or {}).values():
        if (condition.get("priva_scope") or {}).get("decision") != "include":
            continue
        inheritance, x_linked = normalize_inheritance(
            "", (condition.get("inheritance") or {}).get("modes") or []
        )
        if inheritance:
            values.add((inheritance, x_linked))
    distinct = len({inheritance for inheritance, _x in values})
    if distinct != 1:
        return "", False, distinct
    inheritance, x_linked = next(iter(values))
    return inheritance, x_linked, 1


def gene_inheritance_from_constraint(
    symbol: str,
    *,
    clingen: dict[str, dict[str, Any]],
    loeuf: dict[str, float],
) -> str:
    """Last resort for a gene HPO says nothing about: read the constraint data.

    Reached only when no germline condition of the gene states an inheritance at
    all -- 798 genes, for 755 of which HPO holds no inheritance annotation
    anywhere. Rather than deliver nothing, the inheritance is inferred from
    whether the gene tolerates losing one copy:

        dominant   ClinGen haploinsufficiency score 3, or LOEUF below 0.35.
                   Either is a statement that one lost copy already causes
                   disease, which is the dominant pattern.
        recessive  otherwise, as the default. Most disease genes are recessive,
                   and a gene with no haploinsufficiency signal is far more
                   likely to need both copies disabled.

    The two signals barely overlap -- of these genes 196 have only the LOEUF
    signal, 6 only the ClinGen one, 9 both -- so requiring both would reduce
    the rule to 9 genes. Either is therefore enough.

    No mechanism accompanies this. We do not know what a variant must do to
    cause disease here, only how many copies must be affected.
    """
    score = clingen.get(symbol, {}).get("haploinsufficiency_score")
    try:
        haploinsufficient = int(str(score)) == 3
    except (TypeError, ValueError):
        haploinsufficient = False
    constraint = loeuf.get(symbol, float("nan"))
    constrained = (
        isinstance(constraint, (int, float))
        and not math.isnan(constraint)
        and constraint < 0.35
    )
    return "dominant" if haploinsufficient or constrained else "recessive"


def normalize_penetrance(penetrance_hpo_ids: Any = ()) -> str:
    """Reduce condition-linked HPO assertions to PriVA's four categories.

    The shared vocabulary deliberately treats delayed or gradual onset and
    variable expressivity as ``incomplete`` for the downstream observation
    question, while preserving high penetrance as ``high``.  The cache retains
    every original HPO identifier, so this normalization does not replace the
    curator's assertion or detach it from its condition.
    """
    return normalize_penetrance_hpo_ids(penetrance_hpo_ids)
