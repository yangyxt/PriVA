#!/usr/bin/env python3
"""Build PriVA's HPO-framed condition and pathogenic-mechanism JSON cache.

The cache is organized for the production lookup PriVA performs:

    gene -> condition -> inheritance / penetrance / mechanism

HPO establishes the gene-condition frame. Every other source may enrich a
condition only through a stable identifier; disease-name similarity is never a
join key. The source HPO assertion table remains the complete phenotype
resource -- this cache keeps only the three axes ACMG logic reads, plus a count
of how many HPO assertions the full record holds.

WHAT THIS FILE FEEDS
====================

It produces the condition half of the chain documented in
``gene_mechanism_hub.py``. That chain asks, for one variant on one transcript:
which of this gene's germline condition histories does the variant's mechanism
reach, and what inheritance and penetrance do those conditions carry. This
build is what makes those three facts available and complete.

THE BUILD, IN ORDER
===================

Order is not incidental. A curator's assertion must never be overwritten by an
inference, so everything sourced runs before anything deduced.

    1. build_hpo_gene_condition_frame(hpo_assertions)
       ------------------------------------------------------------------
       One condition per (gene, disease) from HPO. The condition key is
       HPO's native disease identifier, preferring OMIM then ORPHA.

       Three HPO axes are retained, each as both a value list and the
       assertion rows behind it:

         inheritance.modes      from HPO_INHERITANCE_TERMS (13 HP: terms)
         penetrance.statuses    from the shared HPO penetrance vocabulary
         onset.terms            from HPO_ONSET_TERMS

    2. _merge_conditions_by_disease(...)   <- inside step 1, before return
       ------------------------------------------------------------------
       Reunites entries that describe ONE disease under different
       identifiers. This is the fix for the split that otherwise hides
       half the facts:

         OMIM:219700  "cystic fibrosis"  MONDO:0009061
                       inheritance: autosomal_recessive    mechanism: --
         ORPHA:586    "cystic fibrosis"  MONDO:0009061
                       inheritance: --                     mechanism: LOF

       HPO annotates inheritance against OMIM; G2P and Orphadata annotate
       mechanisms against their own identifiers. Both already state the
       same MONDO term, so MONDO is the identity that reunites them.

         group by      identifiers.MONDO[0], else the condition key
         surviving key OMIM > ORPHA > MONDO > existing key
         identifiers   union, per namespace
         inheritance / penetrance / onset   union of values and assertions
         mechanisms    union; per mechanism, union allelic_requirements,
                       concatenate evidence, merge variants
         priva_scope   MOST RESTRICTIVE wins, see _SCOPE_PRECEDENCE:
                          exclude  >  review  >  include  >  ""

       Measured on the current data: 12,645 entries -> 11,046 diseases,
       and NO merged group disagrees about recessive versus dominant, so
       the inheritance union needs no tie-break. 15 of 1,595 merging
       groups disagree on priva_scope, which is why that one field has a
       stated precedence rather than a union.

    3. attach_condition_mechanisms(genes, mechanism_evidence)
       ------------------------------------------------------------------
       G2P_DDG2P and Orphadata, joined by condition identifier.
       Writes assertion_basis="curated".

    4. attach_clingen_haploinsufficiency(genes, clingen_dosage)
       ------------------------------------------------------------------
       ClinGen's dosage curation. Reads three columns of the configured
       gene_dosage_sensitivity TSV (the assembly coordinate is not read):

         "#Gene Symbol"
         "Haploinsufficiency Score"
         "Haploinsufficiency Disease ID"     <- a MONDO term, 1,255 of 1,257

       Because that disease identifier is MONDO, this attaches at the
       CONDITION, never at the gene. A haploinsufficient gene is one where
       losing a single copy causes disease, which is a loss-of-function
       statement about that disease.

         score  admitted  ClinGen's own description
         -----  --------  --------------------------------------------
           3      yes     Sufficient evidence for dosage pathogenicity
           2      yes     Some evidence for dosage pathogenicity
           1      yes     Little evidence for dosage pathogenicity
           0      NO      No evidence available
          30      NO      Gene associated with autosomal recessive
                          phenotype  -- ClinGen saying haploinsufficiency
                          does NOT apply; it corroborates step 5, it is
                          not itself a mechanism
          40      NO      Dosage sensitivity unlikely

       0 and 40 are excluded because they assert the ABSENCE of the
       mechanism rather than any degree of it. Writes mechanism="LOF",
       assertion_basis="curated", mechanism_confidence one of
       sufficient_evidence | some_evidence | little_evidence.

    5. attach_gofcards_variants(genes, gofcards_variants)
       ------------------------------------------------------------------
       Exact curated gain-of-function alleles, nested under the condition
       their ClinVar record names. Writes assertion_basis="curated".

    6. deduce_mechanisms_from_inheritance(genes)             <- LAST
       ------------------------------------------------------------------
       Only for a condition with priva_scope.decision == "include" that
       still has NO mechanism from any curated source above.

         inheritance == "recessive"  ->  mechanism = "LOF"
                                         assertion_basis = "deduced"
         anything else               ->  left unresolved

       This is a high-likelihood working inference for variant-to-condition
       matching, not an established molecular-mechanism assertion. In
       particular, PVS1 admits only assertion_basis="curated", so this
       deduction can never open PVS1's gene-level LOF-history gate.

EVERY MECHANISM SAYS HOW IT WAS ESTABLISHED
===========================================

One field on every mechanism evidence entry distinguishes a curator's
molecular-mechanism assertion from PriVA's recessive-to-LOF inference:

    assertion_basis   "curated" | "deduced"

The deduced value is useful for matching a likely LOF variant to a recessive
condition, but it is not promoted to evidence that LOF is an established
disease mechanism. PVS1 therefore ignores it. A condition without either a
curated mechanism or a recessive inheritance remains unresolved for the
runtime hub to carry as condition history.

    source                        basis      count on current data
    ---------------------------   --------   ---------------------
    G2P_DDG2P                     curated    2278
    Orphadata                     curated     933
    ClinGen_haploinsufficiency    curated     330
    GoFCards_exact+ClinVar_VCV    curated     248
    deduced_from_inheritance      deduced    2342

Every one of these is admitted by the runtime reader --
CONDITION_MECHANISM_SOURCES in gene_mechanism_hub.py. Adding a source here
without adding it there writes a mechanism the reader silently ignores, which
is how ClinGen and the deduced entries were invisible when first added. The
two lists must be changed together.

WHAT THE BUILD PRODUCES
=======================

    conditions                            11,046
    conditions with a mechanism            5,433
    conditions with an inheritance          7,156
    condition-variant links                  248
    mechanisms                    LOF | GOF | DOMINANT_NEGATIVE
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator, TextIO

from clinvar_vcv import open_text
from hpo_penetrance import HPO_PENETRANCE_STATUS_BY_TERM


SCHEMA_VERSION = "1.0"

PROJECT_ROOT = Path(__file__).resolve().parents[1]
# ClinGen's dosage curation. Deliberately the GRCh38 copy: the file is a
# per-gene curation and the two assemblies differ only in the coordinate
# column, which this build does not read.
DEFAULT_CLINGEN_DOSAGE = (
    PROJECT_ROOT / "data" / "clingen" / "gene_dosage_sensitivity.hg38.tsv"
)

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

HPO_PENETRANCE_TERMS = HPO_PENETRANCE_STATUS_BY_TERM

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
    "HP:0034857": "variable_age",
    "HP:0003587": "insidious",
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

GOFCARDS_REQUIRED_COLUMNS = {
    "mechanism",
    "HGNC_Symbol",
    "GoFCards_HGNC_Symbol",
    "VEP_HGNC_Symbol",
    "gene_match_status",
    "match_eligibility",
    "HGVSc",
    "HGVSp",
    "hgvsp_key",
    "match_status",
    "gofcards_accession_id",
    "gofcards_variant_id",
    "disease",
    "pmids",
    "pscore",
    "function",
    "pathway",
    "allele_key",
    "hg19_genomic_key",
    "hg19_vcf_key",
    "hg38_genomic_key",
    "hg38_vcf_key",
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
    """Use HPO's native disease ID and retain MONDO as a secondary alias.

    MONDO occasionally maps multiple HPO disease records that have different
    PriVA scope decisions to one ontology term. Native OMIM/ORPHA identity is
    therefore the lossless condition key. MONDO remains usable by the alias
    index only when it resolves to one condition within the queried gene.
    """
    return _clean(row.get("disease_id")).upper() or _clean(
        row.get("mondo_id")
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


# Most restrictive first. When two entries for one disease disagree, the
# merged entry takes the most cautious decision, matching the runtime gate
# which already treats anything but "include" as unusable.
_SCOPE_PRECEDENCE = ("exclude", "review", "include", "")


def _merge_conditions_by_disease(
    conditions: dict[str, dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    """Collapse entries that describe one disease under different identifiers.

    HPO annotates inheritance against OMIM identifiers; G2P and Orphadata
    annotate mechanisms against their own. The same disease therefore arrives
    twice, each copy holding half the facts: cystic fibrosis appears as
    OMIM:219700 carrying autosomal recessive inheritance and no mechanism, and
    as ORPHA:586 carrying a loss-of-function mechanism and no inheritance. Both
    already state MONDO:0009061, so MONDO is the identity that reunites them --
    an exact identifier join, never a match on disease names.

    The surviving key prefers OMIM, then Orphanet, then MONDO, so identifiers
    that already match elsewhere keep matching. Everything else is a union;
    measured across the current cache, no two merged entries disagree about
    recessive versus dominant, so no tie-break is needed for inheritance.
    """
    groups: dict[str, list[tuple[str, dict[str, Any]]]] = defaultdict(list)
    for key, condition in conditions.items():
        mondo = (condition.get("identifiers") or {}).get("MONDO") or []
        groups[_clean(mondo[0]) if mondo else key].append((key, condition))

    merged: dict[str, dict[str, Any]] = {}
    for members in groups.values():
        if len(members) == 1:
            key, condition = members[0]
            merged[key] = condition
            continue

        def rank(item: tuple[str, dict[str, Any]]) -> tuple[int, str]:
            prefix = item[0].split(":", 1)[0].upper()
            return ({"OMIM": 0, "ORPHA": 1, "MONDO": 2}.get(prefix, 3), item[0])

        members.sort(key=rank)
        primary_key, primary = members[0]
        for _key, other in members[1:]:
            for namespace, values in (other.get("identifiers") or {}).items():
                target = primary.setdefault("identifiers", {}).setdefault(namespace, [])
                for value in values:
                    if value not in target:
                        target.append(value)
            for axis, value_field in (
                ("inheritance", "modes"),
                ("penetrance", "statuses"),
                ("onset", "terms"),
            ):
                source = other.get(axis) or {}
                target = primary.setdefault(axis, {"assertions": [], value_field: []})
                for value in source.get(value_field, []) or []:
                    if value not in target.setdefault(value_field, []):
                        target[value_field].append(value)
                for assertion in source.get("assertions", []) or []:
                    if assertion not in target.setdefault("assertions", []):
                        target["assertions"].append(assertion)
            for mechanism, block in (other.get("pathogenic_mechanisms") or {}).items():
                existing = primary.setdefault("pathogenic_mechanisms", {}).get(mechanism)
                if existing is None:
                    primary["pathogenic_mechanisms"][mechanism] = block
                    continue
                for requirement in block.get("allelic_requirements", []) or []:
                    if requirement not in existing.setdefault("allelic_requirements", []):
                        existing["allelic_requirements"].append(requirement)
                existing.setdefault("evidence", []).extend(block.get("evidence", []) or [])
                existing.setdefault("variants", {}).update(block.get("variants", {}) or {})
            primary["hpo_assertion_count"] = (
                primary.get("hpo_assertion_count", 0)
                + other.get("hpo_assertion_count", 0)
            )
            if not _clean(primary.get("label")):
                primary["label"] = other.get("label", "")
            primary_scope = primary.setdefault("priva_scope", {})
            other_scope = other.get("priva_scope") or {}
            if _SCOPE_PRECEDENCE.index(
                _clean(other_scope.get("decision"))
            ) < _SCOPE_PRECEDENCE.index(_clean(primary_scope.get("decision"))):
                primary_scope["decision"] = _clean(other_scope.get("decision"))
                primary_scope["category"] = _clean(other_scope.get("category"))
                primary_scope["review_status"] = _clean(other_scope.get("review_status"))
        merged[primary_key] = primary
    return dict(sorted(merged.items()))


def _new_condition(row: dict[str, str]) -> dict[str, Any]:
    condition = {
        "label": _clean(row.get("mondo_name")),
        "identifiers": {},
        "priva_scope": {
            "decision": _clean(row.get("priva_scope")),
            "category": _clean(row.get("disease_scope")),
            "review_status": _clean(row.get("scope_review_status")),
            "review_statuses": [],
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
    review_status = _clean(row.get("scope_review_status"))
    if review_status and review_status not in scope["review_statuses"]:
        scope["review_statuses"].append(review_status)
    if "manually_confirmed" in scope["review_statuses"]:
        scope["review_status"] = "manually_confirmed"
    elif len(scope["review_statuses"]) == 1:
        scope["review_status"] = scope["review_statuses"][0]
    elif len(scope["review_statuses"]) > 1:
        scope["review_status"] = "mixed"
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
        # Reunite the entries that describe one disease before anything else
        # reads them, so mechanism attachment and every later step see one
        # condition carrying both its inheritance and its mechanism.
        gene["conditions"] = _merge_conditions_by_disease(gene["conditions"])
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
        # Every mechanism records whether a curator asserted it or the build
        # inferred it, so the two can never be mistaken for each other.
        "assertion_basis": "curated",
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


# ClinGen scores that state some evidence for dosage pathogenicity, with the
# wording ClinGen itself uses. 0 ("No evidence available") and 40 ("Dosage
# sensitivity unlikely") are excluded because they assert the absence of the
# mechanism rather than any degree of it, and 30 ("Gene associated with
# autosomal recessive phenotype") is ClinGen saying haploinsufficiency does not
# apply -- it corroborates the recessive reading below, it is not a mechanism.
CLINGEN_HAPLOINSUFFICIENCY_SCORES = {
    "3": "sufficient_evidence",
    "2": "some_evidence",
    "1": "little_evidence",
}


def attach_clingen_haploinsufficiency(
    genes: dict[str, dict[str, Any]],
    clingen_dosage: str | Path,
) -> dict[str, int]:
    """Add ClinGen's haploinsufficiency curation as a loss-of-function mechanism.

    ClinGen states its haploinsufficiency verdict against a MONDO disease, which
    is the same identity the conditions were merged on, so this attaches at the
    condition and never at the gene. A haploinsufficient gene is one where losing
    a single copy causes disease, which is a loss-of-function statement about
    that disease specifically.

    This runs before any mechanism is deduced, so a curator's assertion is never
    overwritten by an inference.
    """
    stats: Counter[str] = Counter()
    required = {
        "#Gene Symbol",
        "Haploinsufficiency Score",
        "Haploinsufficiency Description",
        "Haploinsufficiency Disease ID",
    }
    for row in _iter_tsv_rows(clingen_dosage, required):
        score = _clean(row.get("Haploinsufficiency Score"))
        strength = CLINGEN_HAPLOINSUFFICIENCY_SCORES.get(score)
        disease = _clean(row.get("Haploinsufficiency Disease ID")).upper()
        symbol = _clean(row.get("#Gene Symbol")).upper()
        if not strength or not disease or not symbol:
            stats["rows_without_a_usable_score_or_disease"] += 1
            continue
        gene = genes.get(symbol)
        if gene is None:
            stats["gene_absent_from_the_hpo_frame"] += 1
            continue
        matched = [
            condition
            for condition in (gene.get("conditions") or {}).values()
            if disease in ((condition.get("identifiers") or {}).get("MONDO") or [])
        ]
        if not matched:
            stats["disease_absent_from_the_gene"] += 1
            continue
        for condition in matched:
            had_lof = "LOF" in (condition.get("pathogenic_mechanisms") or {})
            _attach_mechanism_to_condition(
                condition,
                {
                    "source": "ClinGen_haploinsufficiency",
                    "source_record_id": disease,
                    "condition_identifiers": [disease],
                    "condition_label": _clean(row.get("Haploinsufficiency Description")),
                    "mechanism": "LOF",
                    "mechanism_raw": _clean(row.get("Haploinsufficiency Description")),
                    "allelic_requirement": "",
                    "mechanism_confidence": strength,
                    "disease_confidence": "",
                    "assertion_basis": "curated",
                    "source_scope": {"decision": "", "category": "", "review_status": ""},
                    "pmids": [
                        _clean(row.get(f"Haploinsufficiency PMID{index}"))
                        for index in range(1, 7)
                        if _clean(row.get(f"Haploinsufficiency PMID{index}"))
                    ],
                    "evidence_url": "https://search.clinicalgenome.org/kb/gene-dosage",
                },
            )
            stats["conditions_matched"] += 1
            stats["conditions_that_already_had_lof" if had_lof
                  else "conditions_that_gained_lof"] += 1
    return dict(stats)


def deduce_mechanisms_from_inheritance(
    genes: dict[str, dict[str, Any]],
) -> dict[str, int]:
    """Infer LOF for an included recessive condition with no stated mechanism.

    Recessive inheritance makes loss of function a strong working hypothesis
    for matching a query variant to the condition. The record is nevertheless
    marked ``assertion_basis: deduced`` rather than ``curated``. This distinction
    is functional: PVS1's gene-level LOF-history gate accepts only curated
    assertions, while the variant-condition mechanism selector may use this
    inference.

    This runs after every curated source and never overwrites one. Dominant,
    mitochondrial, Y-linked, non-Mendelian, and inheritance-unknown conditions
    remain unresolved rather than receiving a guessed mechanism.
    """
    from gene_mechanism_hub import normalize_inheritance

    stats: Counter[str] = Counter()
    for symbol, gene in genes.items():
        for condition_key, condition in (gene.get("conditions") or {}).items():
            if (condition.get("priva_scope") or {}).get("decision") != "include":
                continue
            if condition.get("pathogenic_mechanisms"):
                stats["already_had_a_mechanism"] += 1
                continue
            inheritance, _x_linked = normalize_inheritance(
                "", (condition.get("inheritance") or {}).get("modes") or []
            )
            if inheritance != "recessive":
                stats[f"left_unresolved_{inheritance or 'no_inheritance'}"] += 1
                continue
            _attach_mechanism_to_condition(
                condition,
                {
                    "source": "deduced_from_inheritance",
                    "source_record_id": condition_key,
                    "condition_identifiers": [condition_key],
                    "condition_label": _clean(condition.get("label")),
                    "mechanism": "LOF",
                    "mechanism_raw": "recessive inheritance implies loss of function",
                    "allelic_requirement": "",
                    "mechanism_confidence": "",
                    "disease_confidence": "",
                    "assertion_basis": "deduced",
                    "source_scope": {"decision": "", "category": "", "review_status": ""},
                    "pmids": [],
                    "evidence_url": "",
                },
            )
            stats["recessive_conditions_given_lof"] += 1
    return dict(stats)


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
    unmapped_variant_counts: dict[str, int] = defaultdict(int)
    for variant in gene["unmapped_evidence"]["variants"].values():
        mechanism = _clean(variant.get("mechanism")).upper()
        if mechanism:
            mechanisms.add(mechanism)
            unmapped_variant_counts[mechanism] += 1
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
            "unmapped_variant_mechanism_counts": dict(
                sorted(unmapped_variant_counts.items())
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


# ---------------------------------------------------------------------------
# Exact GoFCards variant and ClinVar condition enrichment
# ---------------------------------------------------------------------------


def _normalize_clinvar_condition_identifier(
    namespace: Any,
    value: Any,
) -> str:
    """Normalize only registry identifiers supported by exact identity joins."""
    database = re.sub(r"[^A-Z0-9]", "", _clean(namespace).upper())
    identifier = _clean(value).upper()
    if not database or not identifier:
        return ""
    # Only the three registries the HPO table is keyed by. It carries 12,458
    # MONDO, 7,360 OMIM and 5,285 ORPHA identifiers and no MedGen at all, so a
    # MedGen identifier has nothing to join to and admitting it would only
    # produce identifiers that can never match. ClinVar's `Preferred` and `HP`
    # values sit in the same place in the structure and are excluded for a
    # different reason: one is a name qualifier and the other is a phenotype
    # term, neither is a disease identifier.
    prefixes = {
        "OMIM": "OMIM",
        "MONDO": "MONDO",
        "ORPHA": "ORPHA",
        "ORPHANET": "ORPHA",
    }
    prefix = prefixes.get(database)
    if not prefix:
        return ""
    if ":" in identifier:
        incoming_prefix, identifier = identifier.split(":", 1)
        if incoming_prefix in {"ORPHANET", "ORPHA"}:
            incoming_prefix = "ORPHA"
        if incoming_prefix != prefix:
            return ""
    return f"{prefix}:{identifier}"


def _clinvar_condition_identities(
    condition_assertion: dict[str, Any],
) -> list[str]:
    """Extract database identifiers without disease-name or PMID matching."""
    identities: list[str] = []
    for condition in condition_assertion.get("conditions", []) or []:
        if not isinstance(condition, dict):
            continue
        identifier = _normalize_clinvar_condition_identifier(
            condition.get("database"),
            condition.get("id"),
        )
        if identifier and identifier not in identities:
            identities.append(identifier)
    for scv in condition_assertion.get("matched_scvs", []) or []:
        if not isinstance(scv, dict):
            continue
        for mapping in scv.get("trait_mappings", []) or []:
            if not isinstance(mapping, dict):
                continue
            identifier = _normalize_clinvar_condition_identifier(
                mapping.get("mapping_ref"),
                mapping.get("mapping_value"),
            )
            if identifier and identifier not in identities:
                identities.append(identifier)
    return identities


def _append_unique(target: list[str], value: Any) -> None:
    cleaned = _clean(value)
    if cleaned and cleaned != "." and cleaned not in target:
        target.append(cleaned)


def gofcards_variant_key(variant_id: str) -> str:
    """The key a GoFCards variant is published under in this cache."""
    cleaned = _clean(variant_id)
    if not cleaned:
        raise ValueError("GoFCards variant has no identifier")
    return f"GOFCARDS:{cleaned}"


def _gofcards_handover(symbol: str, variant_id: str, variant: dict[str, Any]) -> dict[str, Any]:
    """Everything the GoFCards cache hands to the condition cache, and no more.

    Gene symbol, transcript with its version, and the HGVS notations. Nothing
    else: the curated disease text, pathway, functional description and score
    are evidence about the variant that this cache does not key, join or filter
    on, and restating them here would be a second copy of the GoFCards cache
    free to drift from it.

    Every transcript view on both builds is handed over, not a representative
    one. The transcript *version* differs between GRCh37 and GRCh38 for 1,982 of
    2,028 variants, so a single entry cannot describe both builds, and one
    variant is numbered differently on different isoforms -- JAK2 Val617Phe is
    also Val468Phe -- so a query annotated on any isoform must still be found.
    """
    views: list[dict[str, str]] = []
    for assembly, block in sorted((variant.get("assemblies") or {}).items()):
        for transcript, view in sorted((block.get("transcripts") or {}).items()):
            for hgvsc, detail in sorted((view.get("by_hgvsc") or {}).items()):
                views.append({
                    "assembly": assembly,
                    "transcript": transcript,
                    "hgvsc": hgvsc,
                    "hgvsp": _clean((detail or {}).get("hgvsp")),
                })
    return {
        "mechanism": "GOF",
        "symbol": symbol,
        "gofcards_variant_id": variant_id,
        "transcripts": views,
        "clinvar_links": [],
    }


def _condition_link_mechanism_evidence(
    variant_key: str,
    links: list[dict[str, Any]],
) -> dict[str, Any]:
    identifiers = sorted(
        {
            identifier
            for link in links
            for identifier in link["condition_identifiers"]
        }
    )
    vcv_accessions = sorted(
        {_clean(link.get("vcv_accession")) for link in links}
        - {""}
    )
    return {
        "source": "GoFCards_exact+ClinVar_VCV",
        "source_record_id": ";".join(vcv_accessions) or variant_key,
        "condition_identifiers": identifiers,
        "condition_label": "",
        "mechanism": "GOF",
        "mechanism_raw": "gain of function",
        "allelic_requirement": "",
        # A curator recorded this allele as gain of function in GoFCards; the
        # ClinVar record only supplies the condition it belongs to.
        "assertion_basis": "curated",
        "mechanism_confidence": "exact_variant",
        "disease_confidence": "ClinVar_germline_assertion",
        "source_scope": {"decision": "", "category": "", "review_status": ""},
        "pmids": [],
        "evidence_url": "https://www.ncbi.nlm.nih.gov/clinvar/",
    }


def _public_variant_record(
    variant: dict[str, Any],
    *,
    status: str,
    condition_key: str = "",
    clinvar_links: list[dict[str, Any]] | None = None,
    reason: str = "",
) -> dict[str, Any]:
    record = dict(variant)
    record["condition_link"] = {
        "status": status,
        "condition_key": condition_key,
    }
    if reason:
        record["condition_link"]["reason"] = reason
    record["clinvar_links"] = list(clinvar_links or [])
    return record


def _clinvar_links_by_variant(
    cache: dict[str, Any],
) -> dict[tuple[str, str], list[dict[str, Any]]]:
    """Index ClinVar condition assertions by (gene symbol, variant identifier).

    The injection step nests each ClinVar record under the variant it matched,
    so the link is already explicit. Nothing has to be re-resolved from tokens.

    Takes the parsed cache rather than a path, so the file is read once for this
    index and the eligibility partition together.
    """
    from clinvar_vcv import iter_gofcards_variants

    index: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for symbol, variant_id, variant in iter_gofcards_variants(cache):
        blocks = [variant["clinvar"]] if variant.get("clinvar") else []
        blocks += variant.get("clinvar_additional") or []
        for block in blocks:
            vcv = _clean(block.get("vcv_accession"))
            hgvs = [h for h in (block.get("hgvs") or []) if _clean(h)]
            for assertion in block.get("condition_assertions") or []:
                conditions = assertion.get("conditions") or []
                index[(symbol.strip().upper(), variant_id)].append({
                    "vcv_accession": vcv,
                    "hgvs": hgvs,
                    "rcv_accession": _clean(assertion.get("rcv_accession")),
                    "clinical_significance": _clean(
                        (assertion.get("germline_classification") or {}).get(
                            "clinical_significance"
                        )
                    ),
                    "review_stars": (assertion.get("germline_classification") or {}).get(
                        "review_stars"
                    ),
                    # Reuse the existing extractor: it also reads the trait
                    # mappings on each submitted record, which is where the OMIM
                    # and Orphanet identifiers live. The condition list alone
                    # carries only MedGen, which joins to nothing.
                    "condition_identifiers": _clinvar_condition_identities(assertion),
                    "condition_names": [
                        _clean(c.get("name")) for c in conditions if _clean(c.get("name"))
                    ],
                })
    return dict(index)


def attach_gofcards_variants(
    genes: dict[str, dict[str, Any]],
    gofcards_variants: str | Path,
) -> dict[str, int]:
    """Nest eligible condition-linked GOF alleles and retain unresolved ones.

    Quarantined variants remain auditable in the canonical non-LOF JSON; they
    must not be copied into the condition cache or exposed as candidate exact
    alleles.

    The injected GoFCards cache is the only input. ClinVar conditions were
    attached to each variant upstream by inject_clinvar_into_gofcards.py, so
    they are read straight off the variant rather than re-resolved through the
    non-LOF JSON.
    """
    from clinvar_vcv import load_gofcards_cache, partition_gofcards_variants

    # Read once; both the ClinVar link index and the eligibility partition are
    # built over that same parsed object.
    cache = load_gofcards_cache(Path(gofcards_variants))
    clinvar_links = _clinvar_links_by_variant(cache)
    aliases_by_gene = _condition_alias_index(genes)
    eligible, quarantined = partition_gofcards_variants(cache)

    stats = {
        "source_variants": len(eligible) + len(quarantined),
        "eligible_variants": len(eligible),
        "quarantined_variants": len(quarantined),
        "condition_linked_variants": 0,
        "condition_variant_links": 0,
        "unmapped_variants": 0,
    }
    for source_symbol, variant_id, source_variant in sorted(
        eligible, key=lambda item: (item[0].upper(), item[1])
    ):
        symbol = source_symbol.upper()
        variant_key = gofcards_variant_key(variant_id)
        variant = _gofcards_handover(source_symbol, variant_id, source_variant)
        gene = genes.setdefault(symbol, _empty_gene())
        links = list(
            {
                (
                    link["vcv_accession"],
                    link["rcv_accession"],
                    ";".join(link["condition_identifiers"]),
                ): link
                for link in clinvar_links.get((symbol, variant_id), [])
            }.values()
        )
        by_condition: dict[str, list[dict[str, Any]]] = defaultdict(list)
        ambiguous_link = False
        for link in links:
            matched_keys = {
                aliases_by_gene.get(symbol, {}).get(identifier)
                for identifier in link["condition_identifiers"]
                if aliases_by_gene.get(symbol, {}).get(identifier)
            }
            if len(matched_keys) == 1:
                by_condition[next(iter(matched_keys))].append(link)
            elif len(matched_keys) > 1:
                ambiguous_link = True

        if by_condition:
            stats["condition_linked_variants"] += 1
            for condition_key, condition_links in sorted(by_condition.items()):
                condition = gene["conditions"][condition_key]
                evidence = _condition_link_mechanism_evidence(
                    variant_key,
                    condition_links,
                )
                _attach_mechanism_to_condition(condition, evidence)
                condition["pathogenic_mechanisms"]["GOF"]["variants"][
                    variant_key
                ] = _public_variant_record(
                    variant,
                    status="exact",
                    condition_key=condition_key,
                    clinvar_links=condition_links,
                )
                stats["condition_variant_links"] += 1
            continue

        reason = (
            "ambiguous_clinvar_condition_identifiers"
            if ambiguous_link
            else (
                "no_exact_hpo_condition_identifier"
                if links
                else "no_exact_clinvar_condition_link"
            )
        )
        gene["unmapped_evidence"]["variants"][variant_key] = (
            _public_variant_record(
                variant,
                status="unresolved",
                clinvar_links=links,
                reason=reason,
            )
        )
        stats["unmapped_variants"] += 1

    for gene in genes.values():
        gene["unmapped_evidence"]["variants"] = dict(
            sorted(gene["unmapped_evidence"]["variants"].items())
        )
        for condition in gene["conditions"].values():
            condition["pathogenic_mechanisms"] = dict(
                sorted(condition["pathogenic_mechanisms"].items())
            )
            for block in condition["pathogenic_mechanisms"].values():
                block["variants"] = dict(sorted(block["variants"].items()))
                block["evidence"].sort(key=_mechanism_evidence_key)
        _refresh_gene_summary(gene)
    return stats


# ---------------------------------------------------------------------------
# Cache assembly, validation, and atomic publication
# ---------------------------------------------------------------------------


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_metadata(path: str | Path) -> dict[str, Any]:
    source = Path(path).resolve()
    return {
        "path": str(source),
        "size_bytes": source.stat().st_size,
        "sha256": _sha256_file(source),
    }


def _cache_counts(genes: dict[str, dict[str, Any]]) -> dict[str, int]:
    conditions = [
        condition
        for gene in genes.values()
        for condition in gene["conditions"].values()
    ]
    resolved_mechanism_conditions = sum(
        bool(condition["pathogenic_mechanisms"]) for condition in conditions
    )
    condition_variants = sum(
        len(block["variants"])
        for condition in conditions
        for block in condition["pathogenic_mechanisms"].values()
    )
    return {
        "genes": len(genes),
        "genes_with_hpo_conditions": sum(
            bool(gene["conditions"]) for gene in genes.values()
        ),
        "conditions": len(conditions),
        "conditions_with_inheritance": sum(
            bool(condition["inheritance"]["modes"]) for condition in conditions
        ),
        "conditions_with_penetrance": sum(
            bool(condition["penetrance"]["statuses"]) for condition in conditions
        ),
        "conditions_with_pathogenic_mechanisms": resolved_mechanism_conditions,
        "condition_variant_links": condition_variants,
        "unmapped_mechanisms": sum(
            len(gene["unmapped_evidence"]["mechanisms"])
            for gene in genes.values()
        ),
        "unmapped_variants": sum(
            len(gene["unmapped_evidence"]["variants"])
            for gene in genes.values()
        ),
    }


def build_cache_payload(
    *,
    hpo_assertions: str | Path,
    mechanism_evidence: str | Path,
    gofcards_variants: str | Path,
    clingen_dosage: str | Path = DEFAULT_CLINGEN_DOSAGE,
    hpo_release: str = "",
    mondo_release: str = "",
) -> dict[str, Any]:
    """Build the complete single-file cache from prepared PriVA resources."""
    genes = build_hpo_gene_condition_frame(hpo_assertions)
    mechanism_stats = attach_condition_mechanisms(genes, mechanism_evidence)
    # Curated sources run first. Only after they are exhausted may recessive
    # inheritance supply a clearly marked, deduced LOF working hypothesis.
    # PVS1 later admits only curated LOF history, so this inference guides
    # variant-condition matching without becoming PVS1 evidence.
    clingen_stats = attach_clingen_haploinsufficiency(genes, clingen_dosage)
    variant_stats = attach_gofcards_variants(genes, gofcards_variants)
    deduced_stats = deduce_mechanisms_from_inheritance(genes)
    genes = dict(sorted(genes.items()))
    payload = {
        "_meta": {
            "schema_version": SCHEMA_VERSION,
            "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "releases": {
                "HPO": _clean(hpo_release),
                "MONDO": _clean(mondo_release),
            },
            "sources": {
                "builder_script": _source_metadata(Path(__file__).resolve()),
                "hpo_assertions": _source_metadata(hpo_assertions),
                "mechanism_evidence": _source_metadata(mechanism_evidence),
                "clingen_dosage": _source_metadata(clingen_dosage),
                "gofcards_variants": _source_metadata(gofcards_variants),
            },
            "build_statistics": {
                "mechanisms": mechanism_stats,
                "clingen_haploinsufficiency": clingen_stats,
                "variants": variant_stats,
                "deduced_mechanisms": deduced_stats,
            },
            "counts": _cache_counts(genes),
        },
        "genes": genes,
    }
    validate_cache_payload(payload)
    return payload


def validate_cache_payload(payload: dict[str, Any]) -> dict[str, int]:
    """Validate structural and identity invariants before cache publication."""
    if not isinstance(payload, dict):
        raise ValueError("Cache payload must be a JSON object")
    meta = payload.get("_meta")
    genes = payload.get("genes")
    if not isinstance(meta, dict) or meta.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"Cache schema_version must be {SCHEMA_VERSION}")
    if not isinstance(genes, dict) or not genes:
        raise ValueError("Cache must contain at least one gene")

    for symbol, gene in genes.items():
        if not symbol or not isinstance(gene, dict):
            raise ValueError("Cache contains an invalid gene record")
        conditions = gene.get("conditions")
        unmapped = gene.get("unmapped_evidence")
        if not isinstance(conditions, dict) or not isinstance(unmapped, dict):
            raise ValueError(f"{symbol} is missing condition or unmapped sections")
        for condition_key, condition in conditions.items():
            identifiers = condition.get("identifiers", {})
            all_identifiers = {
                identifier
                for values in identifiers.values()
                for identifier in values
            }
            if condition_key not in all_identifiers:
                raise ValueError(
                    f"{symbol}/{condition_key} canonical key is absent from identifiers"
                )
            for mechanism, block in condition.get(
                "pathogenic_mechanisms", {}
            ).items():
                if mechanism not in CANONICAL_MECHANISMS:
                    raise ValueError(
                        f"{symbol}/{condition_key} has unsupported mechanism {mechanism}"
                    )
                for variant_key, variant in block.get("variants", {}).items():
                    link = variant.get("condition_link", {})
                    if link.get("status") != "exact" or link.get(
                        "condition_key"
                    ) != condition_key:
                        raise ValueError(
                            f"{symbol}/{condition_key}/{variant_key} lacks an exact link"
                        )
        for variant_key, variant in unmapped.get("variants", {}).items():
            if variant.get("condition_link", {}).get("status") == "exact":
                raise ValueError(
                    f"{symbol}/{variant_key} is exact but stored as unmapped"
                )

    actual_counts = _cache_counts(genes)
    recorded_counts = meta.get("counts")
    if recorded_counts != actual_counts:
        raise ValueError(
            f"Cache count mismatch: recorded={recorded_counts}, actual={actual_counts}"
        )
    return actual_counts


def load_and_validate_cache(path: str | Path) -> dict[str, int]:
    with open_text(path) as handle:
        payload = json.load(handle)
    return validate_cache_payload(payload)


def validate_cache_provenance(
    path: str | Path,
    *,
    hpo_assertions: str | Path,
    mechanism_evidence: str | Path,
    clingen_dosage: str | Path,
    gofcards_variants: str | Path,
    hpo_release: str,
    mondo_release: str,
) -> dict[str, int]:
    """Require a valid cache built from the exact files configured now."""
    with open_text(path) as handle:
        payload = json.load(handle)
    counts = validate_cache_payload(payload)
    expected_sources = {
        "builder_script": Path(__file__).resolve(),
        "hpo_assertions": hpo_assertions,
        "mechanism_evidence": mechanism_evidence,
        "clingen_dosage": clingen_dosage,
        "gofcards_variants": gofcards_variants,
    }
    recorded_sources = payload["_meta"].get("sources", {})
    mismatches = []
    for name, source in expected_sources.items():
        actual = _source_metadata(source)
        recorded = recorded_sources.get(name)
        if not isinstance(recorded, dict):
            mismatches.append(f"{name}: missing provenance")
            continue
        changed_fields = [
            field
            for field in ("path", "size_bytes", "sha256")
            if recorded.get(field) != actual[field]
        ]
        if changed_fields:
            mismatches.append(f"{name}: changed {','.join(changed_fields)}")

    recorded_releases = payload["_meta"].get("releases", {})
    for name, expected in (
        ("HPO", _clean(hpo_release)),
        ("MONDO", _clean(mondo_release)),
    ):
        if recorded_releases.get(name) != expected:
            mismatches.append(
                f"{name} release: recorded={recorded_releases.get(name)!r}, "
                f"configured={expected!r}"
            )
    if mismatches:
        raise ValueError("Cache provenance mismatch: " + "; ".join(mismatches))
    return counts


def write_json_atomic(
    payload: dict[str, Any],
    output: str | Path,
    *,
    pretty: bool = False,
) -> None:
    """Publish a validated JSON snapshot with same-filesystem atomic replace.

    The temporary file is fully written, flushed, reparsed, and validated
    before ``os.replace`` changes the live path. Existing readers therefore see
    either the previous complete cache or the next complete cache.
    """
    validate_cache_payload(payload)
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    # The temporary name ends the same way the destination does, because that
    # suffix is what decides whether the bytes are compressed -- and the
    # temporary file is reparsed before it is published, so it has to be
    # readable in the same form the destination will be.
    suffix = ".tmp.gz" if str(destination).endswith(".gz") else ".tmp"
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=suffix,
        dir=destination.parent,
        text=False,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with open_text(temporary, "wt") as handle:
            json.dump(
                payload,
                handle,
                ensure_ascii=True,
                indent=2 if pretty else None,
                separators=None if pretty else (",", ":"),
                sort_keys=True,
            )
            handle.write("\n")
        load_and_validate_cache(temporary)
        os.chmod(temporary, 0o644)
        os.replace(temporary, destination)
        try:
            directory_fd = os.open(destination.parent, os.O_RDONLY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
        except OSError:
            # Some network filesystems do not support directory fsync. The
            # same-directory os.replace above still prevents partial content.
            pass
    finally:
        temporary.unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hpo-assertions")
    parser.add_argument("--mechanism-evidence")
    parser.add_argument("--clingen-dosage")
    parser.add_argument("--gofcards-variants")
    parser.add_argument("--hpo-release", default=None)
    parser.add_argument("--mondo-release", default=None)
    parser.add_argument("--output")
    parser.add_argument(
        "--validate-only",
        metavar="CACHE_JSON",
        help="Validate an existing cache and exit without rebuilding",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Write indented JSON instead of the compact production form",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.validate_only:
        provenance_inputs = (
            args.hpo_assertions,
            args.mechanism_evidence,
            args.clingen_dosage,
            args.gofcards_variants,
            args.hpo_release,
            args.mondo_release,
        )
        if any(value is not None for value in provenance_inputs):
            if any(value is None for value in provenance_inputs):
                raise SystemExit(
                    "Provenance validation requires all build inputs and releases"
                )
            counts = validate_cache_provenance(
                args.validate_only,
                hpo_assertions=args.hpo_assertions,
                mechanism_evidence=args.mechanism_evidence,
                clingen_dosage=args.clingen_dosage,
                gofcards_variants=args.gofcards_variants,
                hpo_release=args.hpo_release,
                mondo_release=args.mondo_release,
            )
        else:
            counts = load_and_validate_cache(args.validate_only)
        print(json.dumps(counts, sort_keys=True))
        return
    clingen_dosage = args.clingen_dosage or str(DEFAULT_CLINGEN_DOSAGE)
    required = {
        "--hpo-assertions": args.hpo_assertions,
        "--mechanism-evidence": args.mechanism_evidence,
        "--clingen-dosage": clingen_dosage,
        "--gofcards-variants": args.gofcards_variants,
        "--output": args.output,
    }
    missing = [option for option, value in required.items() if not value]
    if missing:
        raise SystemExit(f"Missing required build arguments: {', '.join(missing)}")
    payload = build_cache_payload(
        hpo_assertions=args.hpo_assertions,
        mechanism_evidence=args.mechanism_evidence,
        clingen_dosage=clingen_dosage,
        gofcards_variants=args.gofcards_variants,
        hpo_release=args.hpo_release or "",
        mondo_release=args.mondo_release or "",
    )
    write_json_atomic(payload, args.output, pretty=args.pretty)
    print(json.dumps(payload["_meta"]["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
