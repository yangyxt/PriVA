#!/usr/bin/env python3
"""Build PriVA's HPO-framed condition and pathogenic-mechanism JSON cache.

The cache is intentionally organized for the production lookup performed by
PriVA: gene -> condition -> inheritance/penetrance/mechanism. HPO establishes
the gene-condition frame. G2P, Orphadata, ClinVar, and GoFCards may enrich a
condition only through stable identifiers; disease-name similarity is never a
join key.

The source HPO assertion table remains the complete phenotype resource. This
runtime cache retains assertion-level evidence for inheritance, penetrance,
and onset, which are the HPO axes used by ACMG logic, while recording the total
number of HPO assertions for each condition.
"""

from __future__ import annotations

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterator, TextIO


SCHEMA_VERSION = "1.0"

HPO_INHERITANCE_TERMS = {
    "HP:0000006": "autosomal_dominant",
    "HP:0000007": "autosomal_recessive",
    "HP:0001417": "x_linked",
    "HP:0001419": "x_linked_recessive",
    "HP:0001423": "x_linked_dominant",
    "HP:0001426": "non_mendelian",
    "HP:0001427": "mitochondrial",
    "HP:0001450": "y_linked",
    "HP:0010982": "polygenic",
    "HP:0010983": "oligogenic",
    "HP:0010984": "digenic",
    "HP:0012275": "autosomal_dominant_maternal_imprinting",
    "HP:0034341": "pseudoautosomal_recessive",
}

HPO_PENETRANCE_TERMS = {
    "HP:0003829": "incomplete",
    "HP:0034950": "complete",
    "HP:4000159": "moderate",
    "HP:4000160": "low",
}

HPO_ONSET_TERMS = {
    "HP:0030674": "antenatal",
    "HP:0011460": "embryonal",
    "HP:0011461": "fetal",
    "HP:0003577": "congenital",
    "HP:0003623": "neonatal",
    "HP:0003593": "infantile",
    "HP:0011463": "childhood",
    "HP:0003621": "juvenile",
    "HP:0410280": "pediatric",
    "HP:0003581": "adult",
    "HP:0011462": "young_adult",
    "HP:0003596": "middle_age",
    "HP:0003584": "late",
}

HPO_REQUIRED_COLUMNS = {
    "gene_symbol",
    "disease_id",
    "hpo_id",
    "frequency",
    "evidence",
    "reference",
    "mondo_id",
    "mondo_name",
    "disease_scope",
    "priva_scope",
    "scope_evidence",
    "scope_reference",
    "scope_review_status",
}

CONDITION_MECHANISM_SOURCES = {"G2P_DDG2P", "Orphadata"}
CANONICAL_MECHANISMS = {"LOF", "GOF", "DOMINANT_NEGATIVE"}
MECHANISM_REQUIRED_COLUMNS = {
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
    "evidence_url",
}


# ---------------------------------------------------------------------------
# Generic input and normalization helpers
# ---------------------------------------------------------------------------


def _clean(value: Any) -> str:
    """Return a stripped source value without converting missing data to text."""
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "null"} else text


def _split_multi(value: Any) -> list[str]:
    """Split PriVA's semicolon-delimited provenance fields deterministically."""
    return list(
        dict.fromkeys(
            part.strip()
            for part in _clean(value).split(";")
            if part.strip() and part.strip() != "-"
        )
    )


def _open_text(path: str | Path) -> TextIO:
    source = Path(path)
    if source.name.endswith(".gz"):
        return gzip.open(source, "rt", encoding="utf-8", newline="")
    return source.open(encoding="utf-8", newline="")


def _iter_tsv_rows(
    path: str | Path,
    required_columns: set[str],
) -> Iterator[dict[str, str]]:
    """Stream a TSV after validating its complete input contract."""
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required_columns - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
        yield from reader


def _add_identifier(target: dict[str, list[str]], identifier: Any) -> None:
    """Add one namespaced disease identifier without losing native HPO IDs."""
    cleaned = _clean(identifier).upper()
    if not cleaned or ":" not in cleaned:
        return
    namespace = cleaned.split(":", 1)[0]
    values = target.setdefault(namespace, [])
    if cleaned not in values:
        values.append(cleaned)


def _condition_key(row: dict[str, str]) -> str:
    """Prefer exact MONDO identity and otherwise retain HPO's native disease ID."""
    return _clean(row.get("mondo_id")).upper() or _clean(
        row.get("disease_id")
    ).upper()


def _evidence_record(row: dict[str, str]) -> dict[str, str]:
    """Preserve the inseparable HPO frequency/evidence/reference assertion."""
    return {
        "hpo_id": _clean(row.get("hpo_id")).upper(),
        "frequency": _clean(row.get("frequency")) or "-",
        "evidence": _clean(row.get("evidence")) or "-",
        "reference": _clean(row.get("reference")) or "-",
    }


# ---------------------------------------------------------------------------
# HPO frame construction
# ---------------------------------------------------------------------------


def _new_condition(row: dict[str, str]) -> dict[str, Any]:
    condition = {
        "label": _clean(row.get("mondo_name")),
        "identifiers": {},
        "priva_scope": {
            "decision": _clean(row.get("priva_scope")),
            "category": _clean(row.get("disease_scope")),
            "review_status": _clean(row.get("scope_review_status")),
            "evidence": [],
            "references": [],
        },
        "inheritance": {"modes": [], "assertions": []},
        "penetrance": {"statuses": [], "assertions": []},
        "onset": {"terms": [], "assertions": []},
        "pathogenic_mechanisms": {},
        "hpo_assertion_count": 0,
        "_assertion_keys": {
            "inheritance": set(),
            "penetrance": set(),
            "onset": set(),
        },
    }
    _add_identifier(condition["identifiers"], row.get("disease_id"))
    _add_identifier(condition["identifiers"], row.get("mondo_id"))
    return condition


def _merge_scope(condition: dict[str, Any], row: dict[str, str]) -> None:
    """Merge scope provenance and reject contradictory canonical conditions."""
    scope = condition["priva_scope"]
    for target, source in (
        ("decision", "priva_scope"),
        ("category", "disease_scope"),
        ("review_status", "scope_review_status"),
    ):
        incoming = _clean(row.get(source))
        current = _clean(scope.get(target))
        if current and incoming and current != incoming:
            raise ValueError(
                "Conflicting HPO condition scope after canonical identity merge: "
                f"{target}={current!r} versus {incoming!r}"
            )
        if incoming:
            scope[target] = incoming
    for value in _split_multi(row.get("scope_evidence")):
        if value not in scope["evidence"]:
            scope["evidence"].append(value)
    for value in _split_multi(row.get("scope_reference")):
        if value not in scope["references"]:
            scope["references"].append(value)


def _append_hpo_axis_assertion(
    condition: dict[str, Any],
    axis: str,
    normalized_value: str,
    row: dict[str, str],
) -> None:
    value_key = "terms" if axis == "onset" else (
        "statuses" if axis == "penetrance" else "modes"
    )
    if normalized_value not in condition[axis][value_key]:
        condition[axis][value_key].append(normalized_value)

    evidence = _evidence_record(row)
    evidence_key = tuple(evidence.values())
    if evidence_key not in condition["_assertion_keys"][axis]:
        condition["_assertion_keys"][axis].add(evidence_key)
        condition[axis]["assertions"].append(evidence)


def _summarize_hpo_gene(gene: dict[str, Any]) -> dict[str, Any]:
    conditions = list(gene["conditions"].values())
    return {
        "condition_count": len(conditions),
        "conditions_with_inheritance": sum(
            bool(condition["inheritance"]["modes"]) for condition in conditions
        ),
        "conditions_with_penetrance": sum(
            bool(condition["penetrance"]["statuses"]) for condition in conditions
        ),
        "conditions_with_onset": sum(
            bool(condition["onset"]["terms"]) for condition in conditions
        ),
        "pathogenic_mechanisms": [],
        "condition_resolved_mechanism_counts": {},
        "unmapped_mechanism_count": 0,
        "unmapped_variant_count": 0,
    }


def build_hpo_gene_condition_frame(
    hpo_assertions: str | Path,
) -> dict[str, dict[str, Any]]:
    """Build the complete HPO gene-condition frame used by later enrichers."""
    genes: dict[str, dict[str, Any]] = {}
    for row in _iter_tsv_rows(hpo_assertions, HPO_REQUIRED_COLUMNS):
        symbol = _clean(row.get("gene_symbol")).upper()
        condition_key = _condition_key(row)
        if not symbol or not condition_key:
            continue

        gene = genes.setdefault(
            symbol,
            {
                "summary": {},
                "conditions": {},
                "unmapped_evidence": {"mechanisms": [], "variants": {}},
            },
        )
        condition = gene["conditions"].setdefault(
            condition_key,
            _new_condition(row),
        )
        if not condition["label"] and _clean(row.get("mondo_name")):
            condition["label"] = _clean(row.get("mondo_name"))
        _add_identifier(condition["identifiers"], row.get("disease_id"))
        _add_identifier(condition["identifiers"], row.get("mondo_id"))
        _merge_scope(condition, row)
        condition["hpo_assertion_count"] += 1

        hpo_id = _clean(row.get("hpo_id")).upper()
        if hpo_id in HPO_INHERITANCE_TERMS:
            _append_hpo_axis_assertion(
                condition,
                "inheritance",
                HPO_INHERITANCE_TERMS[hpo_id],
                row,
            )
        if hpo_id in HPO_PENETRANCE_TERMS:
            _append_hpo_axis_assertion(
                condition,
                "penetrance",
                HPO_PENETRANCE_TERMS[hpo_id],
                row,
            )
        if hpo_id in HPO_ONSET_TERMS:
            _append_hpo_axis_assertion(
                condition,
                "onset",
                HPO_ONSET_TERMS[hpo_id],
                row,
            )

    ordered_genes: dict[str, dict[str, Any]] = {}
    for symbol in sorted(genes):
        gene = genes[symbol]
        gene["conditions"] = dict(sorted(gene["conditions"].items()))
        for condition in gene["conditions"].values():
            condition.pop("_assertion_keys", None)
            condition["identifiers"] = {
                namespace: sorted(values)
                for namespace, values in sorted(condition["identifiers"].items())
            }
            for axis, values_key in (
                ("inheritance", "modes"),
                ("penetrance", "statuses"),
                ("onset", "terms"),
            ):
                condition[axis][values_key] = sorted(
                    condition[axis][values_key]
                )
                condition[axis]["assertions"].sort(
                    key=lambda item: (
                        item["hpo_id"],
                        item["reference"],
                        item["evidence"],
                        item["frequency"],
                    )
                )
        gene["summary"] = _summarize_hpo_gene(gene)
        ordered_genes[symbol] = gene
    return ordered_genes


# ---------------------------------------------------------------------------
# Exact G2P/Orphadata condition-mechanism enrichment
# ---------------------------------------------------------------------------


def _empty_gene() -> dict[str, Any]:
    """Create a gene shell for curated histories absent from HPO's gene frame."""
    return {
        "summary": {},
        "conditions": {},
        "unmapped_evidence": {"mechanisms": [], "variants": {}},
    }


def _canonical_mechanisms(value: Any) -> list[str]:
    """Return only PriVA's sequence-variant pathogenic mechanism vocabulary."""
    tokens = {
        token.strip().upper()
        for token in re.split(r"[|;,]", _clean(value))
        if token.strip()
    }
    return sorted(tokens & CANONICAL_MECHANISMS)


def _condition_alias_index(
    genes: dict[str, dict[str, Any]],
) -> dict[str, dict[str, str]]:
    """Index only unambiguous stable identifiers within each gene.

    An identifier resolving to multiple canonical conditions is deliberately
    omitted. This makes a failed join visible in ``unmapped_evidence`` instead
    of selecting one condition according to input order.
    """
    aliases_by_gene: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    for symbol, gene in genes.items():
        for condition_key, condition in gene["conditions"].items():
            aliases_by_gene[symbol][condition_key.upper()].add(condition_key)
            for identifiers in condition["identifiers"].values():
                for identifier in identifiers:
                    aliases_by_gene[symbol][identifier.upper()].add(condition_key)
    return {
        symbol: {
            alias: next(iter(condition_keys))
            for alias, condition_keys in aliases.items()
            if len(condition_keys) == 1
        }
        for symbol, aliases in aliases_by_gene.items()
    }


def _mechanism_evidence_record(
    row: dict[str, str],
    mechanism: str,
) -> dict[str, Any]:
    """Retain mechanism provenance without copying unused source columns."""
    identifiers = list(
        dict.fromkeys(
            identifier
            for identifier in (
                _clean(row.get("source_condition_id")).upper(),
                _clean(row.get("mondo_id")).upper(),
            )
            if identifier
        )
    )
    return {
        "source": _clean(row.get("source")),
        "source_record_id": _clean(row.get("source_record_id")),
        "condition_identifiers": identifiers,
        "condition_label": _clean(row.get("disease_label")),
        "mechanism": mechanism,
        "mechanism_raw": _clean(row.get("patho_mode_raw")),
        "allelic_requirement": _clean(row.get("inheritance")),
        "mechanism_confidence": _clean(row.get("mechanism_confidence")),
        "disease_confidence": _clean(row.get("disease_confidence")),
        "source_scope": {
            "decision": _clean(row.get("priva_scope")),
            "category": _clean(row.get("disease_scope")),
            "review_status": _clean(row.get("scope_review_status")),
        },
        "pmids": _split_multi(row.get("pmids")),
        "evidence_url": _clean(row.get("evidence_url")),
    }


def _mechanism_evidence_key(evidence: dict[str, Any]) -> tuple[str, ...]:
    return (
        evidence["source"],
        evidence["source_record_id"],
        ";".join(evidence["condition_identifiers"]),
        evidence["mechanism"],
        evidence["allelic_requirement"],
    )


def _attach_mechanism_to_condition(
    condition: dict[str, Any],
    evidence: dict[str, Any],
) -> None:
    mechanism = evidence["mechanism"]
    block = condition["pathogenic_mechanisms"].setdefault(
        mechanism,
        {"allelic_requirements": [], "evidence": [], "variants": {}},
    )
    requirement = evidence["allelic_requirement"]
    if requirement and requirement not in block["allelic_requirements"]:
        block["allelic_requirements"].append(requirement)
    existing = {_mechanism_evidence_key(item) for item in block["evidence"]}
    if _mechanism_evidence_key(evidence) not in existing:
        block["evidence"].append(evidence)


def _refresh_gene_summary(gene: dict[str, Any]) -> None:
    """Derive the disposable gene summary from condition and unmapped records."""
    summary = _summarize_hpo_gene(gene)
    condition_counts: dict[str, int] = defaultdict(int)
    mechanisms: set[str] = set()
    for condition in gene["conditions"].values():
        for mechanism in condition["pathogenic_mechanisms"]:
            mechanisms.add(mechanism)
            condition_counts[mechanism] += 1
    unmapped_counts: dict[str, int] = defaultdict(int)
    for record in gene["unmapped_evidence"]["mechanisms"]:
        mechanism = _clean(record.get("mechanism")).upper()
        if mechanism:
            mechanisms.add(mechanism)
            unmapped_counts[mechanism] += 1
    summary.update(
        {
            "pathogenic_mechanisms": sorted(mechanisms),
            "condition_resolved_mechanism_counts": dict(
                sorted(condition_counts.items())
            ),
            "unmapped_mechanism_count": sum(unmapped_counts.values()),
            "unmapped_mechanism_counts": dict(sorted(unmapped_counts.items())),
            "unmapped_variant_count": len(
                gene["unmapped_evidence"]["variants"]
            ),
        }
    )
    gene["summary"] = summary


def attach_condition_mechanisms(
    genes: dict[str, dict[str, Any]],
    mechanism_evidence: str | Path,
) -> dict[str, int]:
    """Attach G2P/Orphadata mechanisms through exact condition identifiers.

    Returns build statistics so installation validation can distinguish an
    unexpectedly empty source from a legitimate set of unmatched histories.
    """
    aliases_by_gene = _condition_alias_index(genes)
    stats = {"source_rows": 0, "mechanism_records": 0, "matched": 0, "unmapped": 0}
    seen_unmapped: dict[str, set[tuple[str, ...]]] = defaultdict(set)

    for row in _iter_tsv_rows(mechanism_evidence, MECHANISM_REQUIRED_COLUMNS):
        if _clean(row.get("source")) not in CONDITION_MECHANISM_SOURCES:
            continue
        stats["source_rows"] += 1
        symbol = _clean(row.get("gene_symbol")).upper()
        if not symbol:
            continue
        gene = genes.setdefault(symbol, _empty_gene())
        mechanisms = _canonical_mechanisms(row.get("normalized_mechanisms"))
        for mechanism in mechanisms:
            stats["mechanism_records"] += 1
            evidence = _mechanism_evidence_record(row, mechanism)
            matched_keys = {
                aliases_by_gene.get(symbol, {}).get(identifier)
                for identifier in evidence["condition_identifiers"]
                if aliases_by_gene.get(symbol, {}).get(identifier)
            }
            if len(matched_keys) == 1:
                condition_key = next(iter(matched_keys))
                condition = gene["conditions"][condition_key]
                _attach_mechanism_to_condition(condition, evidence)
                if not condition["label"] and evidence["condition_label"]:
                    condition["label"] = evidence["condition_label"]
                stats["matched"] += 1
                continue

            reason = (
                "conflicting_condition_identifiers"
                if len(matched_keys) > 1
                else "no_exact_hpo_condition_identifier"
            )
            unmapped = {**evidence, "reason": reason}
            identity = _mechanism_evidence_key(unmapped)
            if identity not in seen_unmapped[symbol]:
                seen_unmapped[symbol].add(identity)
                gene["unmapped_evidence"]["mechanisms"].append(unmapped)
                stats["unmapped"] += 1

    for gene in genes.values():
        for condition in gene["conditions"].values():
            condition["pathogenic_mechanisms"] = dict(
                sorted(condition["pathogenic_mechanisms"].items())
            )
            for block in condition["pathogenic_mechanisms"].values():
                block["allelic_requirements"].sort()
                block["evidence"].sort(key=_mechanism_evidence_key)
        gene["unmapped_evidence"]["mechanisms"].sort(
            key=_mechanism_evidence_key
        )
        _refresh_gene_summary(gene)
    return stats
