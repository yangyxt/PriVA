#!/usr/bin/env python3
"""Build PriVA disease-scope annotations from MONDO and HPO evidence.

MONDO supplies disease cross-references and parent-child relationships. HPO
supplies disease-specific inheritance assertions. The resulting registry keeps
automatic classifications separate from manual overrides and never infers a
somatic-only mechanism merely from the absence of hereditary evidence.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


logger = logging.getLogger(__name__)

MONDO_HEREDITARY_ROOT = "MONDO:0003847"
MONDO_NEOPLASM_ROOT = "MONDO:0005070"

HPO_MENDELIAN_INHERITANCE = {
    "HP:0000006",  # Autosomal dominant inheritance
    "HP:0000007",  # Autosomal recessive inheritance
    "HP:0001417",  # X-linked inheritance
    "HP:0001419",  # X-linked recessive inheritance
    "HP:0001423",  # X-linked dominant inheritance
    "HP:0001427",  # Mitochondrial inheritance
    "HP:0001450",  # Y-linked inheritance
    "HP:0012275",  # Autosomal dominant inheritance with maternal imprinting
    "HP:0034341",  # Pseudoautosomal recessive inheritance
}
HPO_NON_MONOGENIC_INHERITANCE = {
    "HP:0001426",  # Non-Mendelian inheritance
    "HP:0010982",  # Polygenic inheritance
    "HP:0010983",  # Oligogenic inheritance
    "HP:0010984",  # Digenic inheritance
}

BASE_ASSERTION_COLUMNS = {
    "gene_symbol",
    "disease_id",
    "hpo_id",
    "frequency",
    "evidence",
    "reference",
}
REGISTRY_COLUMNS = [
    "disease_id",
    "mondo_id",
    "mondo_name",
    "mondo_version",
    "hpo_mendelian_inheritance",
    "hpo_non_monogenic_inheritance",
    "mondo_hereditary",
    "mondo_neoplasm",
    "mondo_predisposition",
    "mondo_omim_susceptibility",
    "disease_scope",
    "priva_scope",
    "scope_evidence",
    "scope_reference",
    "scope_review_status",
    "scope_note",
]
ASSERTION_SCOPE_COLUMNS = [
    "mondo_id",
    "mondo_name",
    "disease_scope",
    "priva_scope",
    "scope_evidence",
    "scope_reference",
    "scope_review_status",
]
OVERRIDE_REQUIRED_COLUMNS = {
    "disease_id",
    "disease_scope",
    "priva_scope",
    "scope_review_status",
}
VALID_DISEASE_SCOPES = {
    "mendelian_non_neoplastic",
    "mendelian_cancer_predisposition",
    "mendelian_susceptibility_uncertain",
    "complex_or_non_monogenic",
    "complex_susceptibility_candidate",
    "somatic_oncogenic",
    "neoplastic_uncertain",
    "conflicting_inheritance",
    "uncertain",
}
VALID_PRIVA_SCOPES = {"include", "exclude", "review"}


@dataclass(frozen=True)
class MondoTerm:
    mondo_id: str
    name: str
    parents: tuple[str, ...]
    xrefs: tuple[str, ...]
    subsets: tuple[str, ...]
    obsolete: bool


def _read_tsv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
        comment="#",
    )


def _normalize_disease_xref(value: str) -> str:
    if value.startswith("OMIM:"):
        return value
    if value.startswith("Orphanet:"):
        return f"ORPHA:{value.split(':', 1)[1]}"
    if value.startswith("ORPHA:"):
        return value
    return ""


def _finish_term(raw: dict[str, object]) -> MondoTerm | None:
    mondo_id = str(raw.get("id", ""))
    if not mondo_id.startswith("MONDO:"):
        return None
    return MondoTerm(
        mondo_id=mondo_id,
        name=str(raw.get("name", "")),
        parents=tuple(raw.get("parents", [])),
        xrefs=tuple(raw.get("xrefs", [])),
        subsets=tuple(raw.get("subsets", [])),
        obsolete=bool(raw.get("obsolete", False)),
    )


def parse_mondo_obo(path: str | Path) -> tuple[str, dict[str, MondoTerm]]:
    """Parse the MONDO fields needed for xref and ancestry classification."""
    version = ""
    terms: dict[str, MondoTerm] = {}
    raw: dict[str, object] | None = None

    with open(path, encoding="utf-8") as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if raw is None and stripped.startswith("data-version:"):
                version = stripped.split(":", 1)[1].strip()
            if stripped == "[Term]":
                if raw is not None:
                    term = _finish_term(raw)
                    if term is not None:
                        terms[term.mondo_id] = term
                raw = {"parents": [], "xrefs": [], "subsets": []}
                continue
            if stripped.startswith("[") and stripped != "[Term]":
                if raw is not None:
                    term = _finish_term(raw)
                    if term is not None:
                        terms[term.mondo_id] = term
                raw = None
                continue
            if raw is None:
                continue
            if stripped.startswith("id: "):
                raw["id"] = stripped[4:].strip()
            elif stripped.startswith("name: "):
                raw["name"] = stripped[6:].strip()
            elif stripped.startswith("is_a: "):
                raw["parents"].append(stripped[6:].split()[0])
            elif stripped.startswith("subset: "):
                raw["subsets"].append(stripped[8:].split()[0])
            elif stripped == "is_obsolete: true":
                raw["obsolete"] = True
            elif stripped.startswith("xref: "):
                xref = stripped[6:].split()[0]
                # Qualified related/broad/narrow mappings are not identity
                # mappings. Bare xrefs and MONDO equivalent mappings are safe.
                if "{" not in stripped or "MONDO:equivalentTo" in stripped:
                    normalized = _normalize_disease_xref(xref)
                    if normalized:
                        raw["xrefs"].append(normalized)

    if raw is not None:
        term = _finish_term(raw)
        if term is not None:
            terms[term.mondo_id] = term

    if not version:
        raise ValueError(f"MONDO OBO has no data-version header: {path}")
    for required_root in (MONDO_HEREDITARY_ROOT, MONDO_NEOPLASM_ROOT):
        if required_root not in terms:
            raise ValueError(f"MONDO OBO is missing required term {required_root}: {path}")
    return version, terms


def _ancestor_resolver(terms: dict[str, MondoTerm]):
    cache: dict[str, frozenset[str]] = {}

    def ancestors(mondo_id: str) -> frozenset[str]:
        if mondo_id in cache:
            return cache[mondo_id]
        seen: set[str] = set()
        stack = list(terms.get(mondo_id, MondoTerm("", "", (), (), (), False)).parents)
        while stack:
            parent = stack.pop()
            if parent in seen:
                continue
            seen.add(parent)
            if parent in cache:
                seen.update(cache[parent])
            else:
                stack.extend(terms.get(parent, MondoTerm("", "", (), (), (), False)).parents)
        cache[mondo_id] = frozenset(seen)
        return cache[mondo_id]

    return ancestors


def _ordered_join(values: Iterable[str]) -> str:
    return ";".join(dict.fromkeys(value for value in values if value))


def _load_overrides(path: str | Path | None) -> dict[str, dict[str, str]]:
    if path is None or not Path(path).exists() or Path(path).stat().st_size == 0:
        return {}
    overrides = _read_tsv(path)
    missing = OVERRIDE_REQUIRED_COLUMNS.difference(overrides.columns)
    if missing:
        raise ValueError(f"Disease-scope overrides are missing columns: {sorted(missing)}")
    if overrides["disease_id"].duplicated().any():
        duplicated = overrides.loc[overrides["disease_id"].duplicated(), "disease_id"].tolist()
        raise ValueError(f"Duplicate disease-scope overrides: {duplicated[:5]}")

    result: dict[str, dict[str, str]] = {}
    for row in overrides.to_dict(orient="records"):
        disease_id = row["disease_id"].strip()
        if not disease_id:
            continue
        disease_scope = row["disease_scope"].strip()
        priva_scope = row["priva_scope"].strip()
        if disease_scope not in VALID_DISEASE_SCOPES:
            raise ValueError(f"Invalid disease_scope for {disease_id}: {disease_scope}")
        if priva_scope not in VALID_PRIVA_SCOPES:
            raise ValueError(f"Invalid priva_scope for {disease_id}: {priva_scope}")
        result[disease_id] = {key: str(value).strip() for key, value in row.items()}
    return result


def _automatic_scope(
    *,
    hpo_mendelian: bool,
    hpo_non_monogenic: bool,
    mondo_hereditary: bool,
    mondo_neoplasm: bool,
    mondo_predisposition: bool,
    mondo_susceptibility: bool,
) -> tuple[str, str, str, str]:
    inherited = hpo_mendelian or mondo_hereditary
    if hpo_mendelian and hpo_non_monogenic:
        return (
            "conflicting_inheritance",
            "review",
            "needs_review",
            "HPO contains both Mendelian and non-monogenic inheritance assertions",
        )
    if hpo_non_monogenic:
        return (
            "complex_or_non_monogenic",
            "exclude",
            "auto_supported",
            "Explicit HPO non-monogenic inheritance",
        )
    if inherited and mondo_neoplasm:
        return (
            "mendelian_cancer_predisposition",
            "include",
            "auto_supported",
            "Inherited evidence and MONDO neoplasm ancestry",
        )
    if inherited and (mondo_predisposition or mondo_susceptibility):
        return (
            "mendelian_susceptibility_uncertain",
            "review",
            "needs_review",
            "Inherited susceptibility/predisposition lacks explicit MONDO neoplasm ancestry",
        )
    if inherited:
        return (
            "mendelian_non_neoplastic",
            "include",
            "auto_supported",
            "Inherited evidence without MONDO neoplasm ancestry",
        )
    if mondo_neoplasm:
        return (
            "neoplastic_uncertain",
            "review",
            "needs_review",
            "MONDO neoplasm ancestry does not establish somatic-only etiology",
        )
    if mondo_predisposition or mondo_susceptibility:
        return (
            "complex_susceptibility_candidate",
            "review",
            "needs_review",
            "Susceptibility/predisposition lacks explicit inherited evidence",
        )
    return (
        "uncertain",
        "review",
        "needs_review",
        "No explicit inherited, non-monogenic, or neoplastic scope evidence",
    )


def build_disease_scope_registry(
    mondo_obo: str | Path,
    hpo_assertions: str | Path,
    *,
    manual_overrides: str | Path | None = None,
) -> pd.DataFrame:
    """Build one scope record per OMIM/ORPHA disease in an HPO assertion table."""
    mondo_version, terms = parse_mondo_obo(mondo_obo)
    assertions = _read_tsv(hpo_assertions)
    missing = BASE_ASSERTION_COLUMNS.difference(assertions.columns)
    if missing:
        raise ValueError(f"HPO assertions are missing columns: {sorted(missing)}")

    disease_hpos = assertions.groupby("disease_id", sort=True)["hpo_id"].agg(set)
    inheritance_refs = (
        assertions.loc[
            assertions["hpo_id"].isin(
                HPO_MENDELIAN_INHERITANCE | HPO_NON_MONOGENIC_INHERITANCE
            )
        ]
        .groupby("disease_id", sort=False)["reference"]
        .agg(lambda values: _ordered_join(value for value in values if value != "-"))
        .to_dict()
    )

    xref_to_mondo: dict[str, set[str]] = defaultdict(set)
    for term in terms.values():
        if term.obsolete:
            continue
        for xref in term.xrefs:
            xref_to_mondo[xref].add(term.mondo_id)
    ancestors = _ancestor_resolver(terms)
    overrides = _load_overrides(manual_overrides)

    records: list[dict[str, object]] = []
    for disease_id, hpo_ids in disease_hpos.items():
        mondo_ids = sorted(xref_to_mondo.get(disease_id, set()))
        mapped_terms = [terms[mondo_id] for mondo_id in mondo_ids]
        mapped_ancestors = {
            ancestor
            for mondo_id in mondo_ids
            for ancestor in ({mondo_id} | set(ancestors(mondo_id)))
        }
        subsets = {subset for term in mapped_terms for subset in term.subsets}
        hpo_mendelian = bool(hpo_ids & HPO_MENDELIAN_INHERITANCE)
        hpo_non_monogenic = bool(hpo_ids & HPO_NON_MONOGENIC_INHERITANCE)
        mondo_hereditary = MONDO_HEREDITARY_ROOT in mapped_ancestors
        mondo_neoplasm = MONDO_NEOPLASM_ROOT in mapped_ancestors
        mondo_predisposition = "predisposition" in subsets
        mondo_susceptibility = "omim_susceptibility" in subsets

        disease_scope, priva_scope, review_status, note = _automatic_scope(
            hpo_mendelian=hpo_mendelian,
            hpo_non_monogenic=hpo_non_monogenic,
            mondo_hereditary=mondo_hereditary,
            mondo_neoplasm=mondo_neoplasm,
            mondo_predisposition=mondo_predisposition,
            mondo_susceptibility=mondo_susceptibility,
        )
        evidence_parts = []
        if hpo_mendelian:
            evidence_parts.append("HPO_mendelian_inheritance")
        if hpo_non_monogenic:
            evidence_parts.append("HPO_non_monogenic_inheritance")
        if mondo_hereditary:
            evidence_parts.append(f"MONDO_ancestor:{MONDO_HEREDITARY_ROOT}")
        if mondo_neoplasm:
            evidence_parts.append(f"MONDO_ancestor:{MONDO_NEOPLASM_ROOT}")
        if mondo_predisposition:
            evidence_parts.append("MONDO_subset:predisposition")
        if mondo_susceptibility:
            evidence_parts.append("MONDO_subset:omim_susceptibility")

        record: dict[str, object] = {
            "disease_id": disease_id,
            "mondo_id": _ordered_join(mondo_ids),
            "mondo_name": _ordered_join(term.name for term in mapped_terms),
            "mondo_version": mondo_version,
            "hpo_mendelian_inheritance": int(hpo_mendelian),
            "hpo_non_monogenic_inheritance": int(hpo_non_monogenic),
            "mondo_hereditary": int(mondo_hereditary),
            "mondo_neoplasm": int(mondo_neoplasm),
            "mondo_predisposition": int(mondo_predisposition),
            "mondo_omim_susceptibility": int(mondo_susceptibility),
            "disease_scope": disease_scope,
            "priva_scope": priva_scope,
            "scope_evidence": _ordered_join(evidence_parts) or "none",
            "scope_reference": _ordered_join(
                [f"MONDO_RELEASE:{mondo_version}", inheritance_refs.get(disease_id, "")]
            ),
            "scope_review_status": review_status,
            "scope_note": note,
        }

        override = overrides.get(disease_id)
        if override is not None:
            for column in (
                "disease_scope",
                "priva_scope",
                "scope_review_status",
                "scope_evidence",
                "scope_reference",
                "scope_note",
            ):
                if override.get(column, ""):
                    record[column] = override[column]
            if not override.get("scope_evidence", ""):
                record["scope_evidence"] = _ordered_join(
                    [str(record["scope_evidence"]), "manual_override"]
                )
            if not override.get("scope_review_status", ""):
                record["scope_review_status"] = "manually_confirmed"

        records.append(record)

    registry = pd.DataFrame.from_records(records, columns=REGISTRY_COLUMNS)
    if registry["disease_id"].duplicated().any():
        raise ValueError("Disease-scope registry contains duplicate disease IDs")
    logger.info(
        "Built scope registry for %d diseases: %s",
        len(registry),
        registry["disease_scope"].value_counts().to_dict(),
    )
    return registry


def annotate_hpo_assertions(
    hpo_assertions: pd.DataFrame,
    registry: pd.DataFrame,
) -> pd.DataFrame:
    """Attach disease scope and provenance to every HPO assertion row."""
    base = hpo_assertions.drop(columns=ASSERTION_SCOPE_COLUMNS, errors="ignore")
    scoped = base.merge(
        registry.loc[:, ["disease_id", *ASSERTION_SCOPE_COLUMNS]],
        on="disease_id",
        how="left",
        validate="many_to_one",
    )
    if scoped["priva_scope"].eq("").any() or scoped["priva_scope"].isna().any():
        missing = scoped.loc[
            scoped["priva_scope"].fillna("").eq(""), "disease_id"
        ].drop_duplicates()
        raise ValueError(f"Disease scope missing for HPO assertions: {missing.head(5).tolist()}")
    return scoped


def _write_tsv_atomic(df: pd.DataFrame, path: str | Path) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    suffix = ".tmp.gz" if output.name.endswith(".gz") else ".tmp"
    temporary = output.with_name(f".{output.name}.{os.getpid()}{suffix}")
    try:
        df.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def build_and_write_resources(
    mondo_obo: str | Path,
    hpo_assertions: str | Path,
    registry_output: str | Path,
    *,
    annotated_assertions_output: str | Path | None = None,
    manual_overrides: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    registry = build_disease_scope_registry(
        mondo_obo,
        hpo_assertions,
        manual_overrides=manual_overrides,
    )
    _write_tsv_atomic(registry, registry_output)

    annotated = None
    if annotated_assertions_output is not None:
        assertions = _read_tsv(hpo_assertions)
        annotated = annotate_hpo_assertions(assertions, registry)
        _write_tsv_atomic(annotated, annotated_assertions_output)
    return registry, annotated


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mondo-obo", required=True)
    parser.add_argument("--hpo-assertions", required=True)
    parser.add_argument("--registry-output", required=True)
    parser.add_argument("--annotated-assertions-output")
    parser.add_argument("--manual-overrides")
    return parser.parse_args()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)d:%(message)s",
    )
    args = parse_args()
    build_and_write_resources(
        args.mondo_obo,
        args.hpo_assertions,
        args.registry_output,
        annotated_assertions_output=args.annotated_assertions_output,
        manual_overrides=args.manual_overrides,
    )
