#!/usr/bin/env python3
"""Inventory exact condition accessions and isolate non-Mendelian gate candidates.

This audit never joins conditions by disease name, gene, PMID, phenotype, or
approximate text. Conditions are connected only when one source record carries
two or more supported stable identifiers (MONDO, OMIM, or ORPHA). Records with
no supported disease identifier remain separate and retain their source keys.

The candidate output is deliberately conservative. Missing inheritance or
mechanism evidence is not evidence of a non-Mendelian condition. Automatic
exclusion is proposed only for an identity component with an explicit
``priva_scope=exclude`` decision and no conflicting ``include`` decision.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, TextIO


SCRIPT_DIR = Path(__file__).resolve().parent
PRIVA_ROOT = SCRIPT_DIR.parent
DEFAULT_CACHE = (
    PRIVA_ROOT / "data" / "patho_mechanism" / "hpo_condition_mechanism_cache.json.gz"
)
DEFAULT_HPO = PRIVA_ROOT / "data" / "hpo" / "genes_to_phenotype.assertions.tsv.gz"
DEFAULT_MECHANISM_EVIDENCE = (
    PRIVA_ROOT / "data" / "patho_mechanism" / "gene_pathogenic_mechanism_evidence.tsv"
)
DEFAULT_SCOPE_REGISTRY = PRIVA_ROOT / "data" / "mondo" / "disease_scope.tsv.gz"
DEFAULT_CLINGEN = PRIVA_ROOT / "data" / "clingen" / "gene_dosage_sensitivity.hg19.tsv"
DEFAULT_GOFCARDS = PRIVA_ROOT / "data" / "gofcards" / "gofcards_exact_gof.json.gz"

SUPPORTED_PREFIXES = ("MONDO", "OMIM", "ORPHA")
IDENTIFIER_PATTERN = re.compile(
    r"(?i)(?:^|[^A-Z0-9])(?P<prefix>MONDO|OMIM|ORPHA|ORPHANET)\s*:\s*"
    r"(?P<value>[A-Z0-9_.-]+)"
)
NON_MENDELIAN_MODES = {"non_mendelian", "polygenic", "digenic", "oligogenic"}
MENDELIAN_MODES = {
    "autosomal_dominant",
    "autosomal_recessive",
    "x_linked",
    "x_linked_dominant",
    "x_linked_recessive",
    "mitochondrial",
    "y_linked",
    "autosomal_dominant_maternal_imprinting",
    "pseudoautosomal_recessive",
}
HPO_INHERITANCE_BY_ID = {
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

INVENTORY_COLUMNS = [
    "canonical_condition_id",
    "mondo_ids",
    "omim_ids",
    "orpha_ids",
    "unresolved_source_keys",
    "condition_labels",
    "gene_symbols",
    "source_names",
    "hpo_disease_ids",
    "hpo_inheritance_term_ids",
    "g2p_condition_ids",
    "g2p_record_ids",
    "orphadata_condition_ids",
    "orphadata_record_ids",
    "clingen_hi_disease_ids",
    "clingen_hi_scores",
    "clingen_ts_disease_ids",
    "clingen_ts_scores",
    "clingen_gene_ids",
    "gofcards_variant_ids",
    "clinvar_vcv_accessions",
    "clinvar_rcv_accessions",
    "cache_condition_keys",
    "priva_scope_decisions",
    "disease_scope_categories",
    "scope_review_statuses",
    "scope_evidence",
    "scope_references",
    "inheritance_modes",
    "source_inheritance_values",
    "pathogenic_mechanisms",
    "source_pathogenic_modes",
    "penetrance_statuses",
    "cache_condition_record_count",
    "source_observation_count",
]

CANDIDATE_COLUMNS = [
    "candidate_class",
    "recommended_cache_gate",
    "positive_non_mendelian_evidence",
    "mendelian_counterevidence",
    *INVENTORY_COLUMNS,
]


def _clean(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "null"}:
        return ""
    return text


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _load_json(path: Path) -> Any:
    with _open_text(path) as handle:
        return json.load(handle)


def _iter_tsv(path: Path, required: set[str]) -> Iterator[dict[str, str]]:
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = sorted(required - fields)
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
        for row in reader:
            yield {key: _clean(value) for key, value in row.items()}


def _split_values(value: Any) -> list[str]:
    return [token.strip() for token in re.split(r"[;|]", _clean(value)) if token.strip()]


def _normalize_identifier(prefix: Any, value: Any) -> str:
    namespace = re.sub(r"[^A-Z0-9]", "", _clean(prefix).upper())
    if namespace == "ORPHANET":
        namespace = "ORPHA"
    identifier = _clean(value).upper()
    if namespace not in SUPPORTED_PREFIXES or not identifier:
        return ""
    if ":" in identifier:
        incoming, identifier = identifier.split(":", 1)
        incoming = "ORPHA" if incoming in {"ORPHA", "ORPHANET"} else incoming
        if incoming != namespace:
            return ""
    elif identifier.startswith(namespace) and identifier[len(namespace) :].isdigit():
        # ClinVar SCV mappings sometimes encode ORPHA:207 as ORPHA207 while the
        # mapping_ref already says Orphanet. Avoid emitting ORPHA:ORPHA207.
        identifier = identifier[len(namespace) :]
    identifier = identifier.strip()
    return f"{namespace}:{identifier}" if identifier else ""


def _extract_identifiers(*values: Any) -> set[str]:
    identifiers: set[str] = set()
    for raw in values:
        text = _clean(raw)
        for match in IDENTIFIER_PATTERN.finditer(f" {text}"):
            identifier = _normalize_identifier(match.group("prefix"), match.group("value"))
            if identifier:
                identifiers.add(identifier)
    return identifiers


def _add_values(target: defaultdict[str, set[str]], key: str, values: Any) -> None:
    if values is None:
        return
    incoming = values if isinstance(values, (set, list, tuple)) else [values]
    for value in incoming:
        cleaned = _clean(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")
        if cleaned:
            target[key].add(cleaned)


@dataclass
class Observation:
    identifiers: set[str]
    values: defaultdict[str, set[str]] = field(default_factory=lambda: defaultdict(set))


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def add(self, item: str) -> None:
        self.parent.setdefault(item, item)

    def find(self, item: str) -> str:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        keep, merge = sorted((left_root, right_root))
        self.parent[merge] = keep


class InventoryBuilder:
    def __init__(self) -> None:
        self.union_find = UnionFind()
        self.observations: list[Observation] = []
        self.source_counts: Counter[str] = Counter()
        self.skipped_records: Counter[str] = Counter()

    def add(
        self,
        identifiers: Iterable[str],
        *,
        synthetic_key: str,
        source_name: str,
        **values: Any,
    ) -> None:
        stable_ids = {identifier for identifier in identifiers if identifier}
        identity_nodes = set(stable_ids)
        if not identity_nodes:
            identity_nodes.add(f"UNRESOLVED::{synthetic_key}")
        for identifier in identity_nodes:
            self.union_find.add(identifier)
        first = min(identity_nodes)
        for identifier in sorted(identity_nodes - {first}):
            self.union_find.union(first, identifier)

        observation = Observation(identity_nodes)
        _add_values(observation.values, "source_names", source_name)
        for key, value in values.items():
            _add_values(observation.values, key, value)
        self.observations.append(observation)
        self.source_counts[source_name] += 1

    def collapse(self) -> list[dict[str, str]]:
        components: dict[str, defaultdict[str, set[str]]] = {}
        component_nodes: dict[str, set[str]] = defaultdict(set)
        observation_counts: Counter[str] = Counter()
        for observation in self.observations:
            root = self.union_find.find(min(observation.identifiers))
            values = components.setdefault(root, defaultdict(set))
            component_nodes[root].update(observation.identifiers)
            observation_counts[root] += 1
            for key, members in observation.values.items():
                values[key].update(members)

        rows: list[dict[str, str]] = []
        for root, values in components.items():
            nodes = component_nodes[root]
            stable = sorted(node for node in nodes if not node.startswith("UNRESOLVED::"))
            unresolved = sorted(
                node.removeprefix("UNRESOLVED::")
                for node in nodes
                if node.startswith("UNRESOLVED::")
            )
            by_prefix = {
                prefix: [value for value in stable if value.startswith(f"{prefix}:")]
                for prefix in SUPPORTED_PREFIXES
            }
            canonical = next(
                (
                    values_for_prefix[0]
                    for prefix in SUPPORTED_PREFIXES
                    if (values_for_prefix := by_prefix[prefix])
                ),
                unresolved[0] if unresolved else root,
            )
            row = {column: "" for column in INVENTORY_COLUMNS}
            row.update(
                {
                    "canonical_condition_id": canonical,
                    "mondo_ids": ";".join(by_prefix["MONDO"]),
                    "omim_ids": ";".join(by_prefix["OMIM"]),
                    "orpha_ids": ";".join(by_prefix["ORPHA"]),
                    "unresolved_source_keys": ";".join(unresolved),
                }
            )
            for column in INVENTORY_COLUMNS:
                if column in {
                    "canonical_condition_id",
                    "mondo_ids",
                    "omim_ids",
                    "orpha_ids",
                    "unresolved_source_keys",
                    "cache_condition_record_count",
                    "source_observation_count",
                }:
                    continue
                row[column] = ";".join(sorted(values.get(column, set())))
            row["cache_condition_record_count"] = str(
                len(values.get("cache_condition_keys", set()))
            )
            row["source_observation_count"] = str(observation_counts[root])
            rows.append(row)
        return sorted(rows, key=lambda row: row["canonical_condition_id"])


def _scope_values(row: dict[str, str]) -> dict[str, Any]:
    return {
        "priva_scope_decisions": row.get("priva_scope", ""),
        "disease_scope_categories": row.get("disease_scope", ""),
        "scope_review_statuses": row.get("scope_review_status", ""),
        "scope_evidence": _split_values(row.get("scope_evidence", "")),
        "scope_references": _split_values(row.get("scope_reference", "")),
    }


def load_scope_registry(builder: InventoryBuilder, path: Path) -> dict[str, int]:
    required = {
        "disease_id",
        "mondo_id",
        "mondo_name",
        "disease_scope",
        "priva_scope",
        "scope_evidence",
        "scope_reference",
        "scope_review_status",
    }
    rows = 0
    for row in _iter_tsv(path, required):
        rows += 1
        identifiers = _extract_identifiers(row["disease_id"], row["mondo_id"])
        builder.add(
            identifiers,
            synthetic_key=f"MONDO_SCOPE:{row['disease_id'] or rows}",
            source_name="MONDO_scope_registry",
            condition_labels=row["mondo_name"],
            **_scope_values(row),
        )
    return {"rows": rows}


def load_hpo_assertions(builder: InventoryBuilder, path: Path) -> dict[str, int]:
    required = {
        "gene_symbol",
        "disease_id",
        "hpo_id",
        "mondo_id",
        "mondo_name",
        "disease_scope",
        "priva_scope",
        "scope_evidence",
        "scope_reference",
        "scope_review_status",
    }
    grouped: dict[tuple[str, str, str], defaultdict[str, set[str]]] = {}
    rows = 0
    for row in _iter_tsv(path, required):
        rows += 1
        key = (row["gene_symbol"], row["disease_id"], row["mondo_id"])
        values = grouped.setdefault(key, defaultdict(set))
        _add_values(values, "condition_labels", row["mondo_name"])
        _add_values(values, "hpo_disease_ids", row["disease_id"])
        if row["hpo_id"] in HPO_INHERITANCE_BY_ID:
            _add_values(values, "hpo_inheritance_term_ids", row["hpo_id"])
            _add_values(values, "inheritance_modes", HPO_INHERITANCE_BY_ID[row["hpo_id"]])
        for field, field_values in _scope_values(row).items():
            _add_values(values, field, field_values)

    for index, ((gene, disease_id, mondo_id), values) in enumerate(grouped.items(), start=1):
        identifiers = _extract_identifiers(disease_id, mondo_id)
        builder.add(
            identifiers,
            synthetic_key=f"HPO:{gene}:{disease_id or index}",
            source_name="HPO",
            gene_symbols=gene,
            **values,
        )
    return {"rows": rows, "gene_condition_groups": len(grouped)}


def load_mechanism_evidence(builder: InventoryBuilder, path: Path) -> dict[str, int]:
    required = {
        "gene_symbol",
        "source",
        "source_record_id",
        "source_condition_id",
        "mondo_id",
        "disease_label",
        "inheritance",
        "normalized_penetrance",
        "patho_mode_raw",
        "normalized_mechanisms",
        "disease_scope",
        "priva_scope",
        "scope_review_status",
    }
    counts: Counter[str] = Counter()
    for index, row in enumerate(_iter_tsv(path, required), start=1):
        source = row["source"] or "condition_mechanism_unknown_source"
        counts[source] += 1
        identifiers = _extract_identifiers(row["source_condition_id"], row["mondo_id"])
        values: dict[str, Any] = {
            "gene_symbols": row["gene_symbol"],
            "condition_labels": row["disease_label"],
            "source_inheritance_values": row["inheritance"],
            "source_pathogenic_modes": row["patho_mode_raw"],
            "pathogenic_mechanisms": _split_values(row["normalized_mechanisms"]),
            "penetrance_statuses": row["normalized_penetrance"],
            **_scope_values(row),
        }
        if source == "G2P_DDG2P":
            values["g2p_condition_ids"] = row["source_condition_id"]
            values["g2p_record_ids"] = row["source_record_id"]
        elif source == "Orphadata":
            values["orphadata_condition_ids"] = row["source_condition_id"]
            values["orphadata_record_ids"] = row["source_record_id"]
        builder.add(
            identifiers,
            synthetic_key=f"MECHANISM:{source}:{row['source_record_id'] or index}",
            source_name=source,
            **values,
        )
    return {"rows": sum(counts.values()), "by_source": dict(sorted(counts.items()))}


def load_clingen_dosage(builder: InventoryBuilder, path: Path) -> dict[str, int]:
    required = {
        "#Gene Symbol",
        "Gene ID",
        "Haploinsufficiency Score",
        "Haploinsufficiency Description",
        "Haploinsufficiency Disease ID",
        "Triplosensitivity Score",
        "Triplosensitivity Description",
        "Triplosensitivity Disease ID",
    }
    counts: Counter[str] = Counter()
    for row_number, row in enumerate(_iter_tsv(path, required), start=2):
        for axis, disease_field, score_field, description_field in (
            (
                "HI",
                "Haploinsufficiency Disease ID",
                "Haploinsufficiency Score",
                "Haploinsufficiency Description",
            ),
            (
                "TS",
                "Triplosensitivity Disease ID",
                "Triplosensitivity Score",
                "Triplosensitivity Description",
            ),
        ):
            disease_id = row[disease_field]
            if not disease_id:
                counts[f"{axis.lower()}_rows_without_disease_id"] += 1
                continue
            counts[f"{axis.lower()}_condition_rows"] += 1
            values: dict[str, Any] = {
                "gene_symbols": row["#Gene Symbol"],
                "clingen_gene_ids": row["Gene ID"],
                "condition_labels": row[description_field],
            }
            values[f"clingen_{axis.lower()}_disease_ids"] = disease_id
            values[f"clingen_{axis.lower()}_scores"] = row[score_field]
            builder.add(
                _extract_identifiers(disease_id),
                synthetic_key=f"CLINGEN:{axis}:{row['#Gene Symbol']}:{row_number}",
                source_name="ClinGen_dosage",
                **values,
            )
    counts["rows"] = row_number - 1 if "row_number" in locals() else 0
    return dict(sorted(counts.items()))


def _clinvar_condition_identifiers(
    assertion: dict[str, Any],
) -> tuple[set[str], set[str]]:
    """Return primary condition IDs separately from SCV trait mappings.

    Multiple IDs in ClinVar's primary condition block describe that RCV's
    condition. In contrast, the matched SCV trait mappings are candidate lookup
    keys and can span several distinct diseases for a pleiotropic variant. They
    must not be unioned into one disease identity component.
    """
    primary: set[str] = set()
    mapped: set[str] = set()
    for condition in assertion.get("conditions", []) or []:
        if not isinstance(condition, dict):
            continue
        identifier = _normalize_identifier(condition.get("database"), condition.get("id"))
        if identifier:
            primary.add(identifier)
    for scv in assertion.get("matched_scvs", []) or []:
        if not isinstance(scv, dict):
            continue
        for mapping in scv.get("trait_mappings", []) or []:
            if not isinstance(mapping, dict):
                continue
            identifier = _normalize_identifier(
                mapping.get("mapping_ref"), mapping.get("mapping_value")
            )
            if identifier:
                mapped.add(identifier)
    return primary, mapped


def load_gofcards_clinvar(builder: InventoryBuilder, path: Path) -> dict[str, int]:
    from clinvar_vcv import iter_gofcards_variants, load_gofcards_cache

    counts: Counter[str] = Counter()
    cache = load_gofcards_cache(path)
    for gene, variant_id, variant in iter_gofcards_variants(cache):
        counts["gofcards_variants"] += 1
        blocks = [variant["clinvar"]] if variant.get("clinvar") else []
        blocks.extend(variant.get("clinvar_additional") or [])
        for block_index, block in enumerate(blocks):
            if not isinstance(block, dict):
                continue
            vcv = _clean(block.get("vcv_accession"))
            for assertion_index, assertion in enumerate(
                block.get("condition_assertions") or []
            ):
                if not isinstance(assertion, dict):
                    continue
                counts["clinvar_condition_assertions"] += 1
                rcv = _clean(assertion.get("rcv_accession"))
                labels = [
                    _clean(condition.get("name"))
                    for condition in assertion.get("conditions", []) or []
                    if isinstance(condition, dict) and _clean(condition.get("name"))
                ]
                primary_ids, mapped_ids = _clinvar_condition_identifiers(assertion)
                base_key = (
                    f"CLINVAR:{vcv or variant_id}:{rcv or block_index}:{assertion_index}"
                )
                common_values = {
                    "gene_symbols": gene,
                    "condition_labels": labels,
                    "gofcards_variant_ids": variant_id,
                    "clinvar_vcv_accessions": vcv,
                    "clinvar_rcv_accessions": rcv,
                }
                if primary_ids:
                    builder.add(
                        primary_ids,
                        synthetic_key=f"{base_key}:PRIMARY",
                        source_name="GoFCards+ClinVar_primary_condition",
                        **common_values,
                    )
                    counts["primary_condition_identity_observations"] += 1

                additional_mappings = sorted(mapped_ids - primary_ids)
                for mapping_index, identifier in enumerate(additional_mappings):
                    builder.add(
                        {identifier},
                        synthetic_key=f"{base_key}:MAPPING:{mapping_index}",
                        source_name="GoFCards+ClinVar_trait_mapping",
                        **common_values,
                    )
                    counts["trait_mapping_identity_observations"] += 1
                if len(additional_mappings) > 1:
                    counts["assertions_with_multiple_additional_trait_mappings"] += 1

                if not primary_ids and not mapped_ids:
                    builder.add(
                        set(),
                        synthetic_key=base_key,
                        source_name="GoFCards+ClinVar_unresolved_condition",
                        **common_values,
                    )
                    counts["assertions_without_supported_condition_id"] += 1
    return dict(sorted(counts.items()))


def load_condition_cache(builder: InventoryBuilder, path: Path) -> tuple[dict[str, int], dict[str, Any]]:
    payload = _load_json(path)
    if not isinstance(payload, dict) or not isinstance(payload.get("genes"), dict):
        raise ValueError(f"{path} is not a PriVA condition cache")
    meta = payload.get("_meta") or {}
    counts: Counter[str] = Counter()
    for gene_symbol, gene in payload["genes"].items():
        if not isinstance(gene, dict):
            continue
        for condition_key, condition in (gene.get("conditions") or {}).items():
            if not isinstance(condition, dict):
                continue
            counts["condition_records"] += 1
            identifiers = _extract_identifiers(condition_key)
            for namespace, values in (condition.get("identifiers") or {}).items():
                for value in values or []:
                    identifier = _normalize_identifier(namespace, value)
                    if identifier:
                        identifiers.add(identifier)
            scope = condition.get("priva_scope") or {}
            values: dict[str, Any] = {
                "gene_symbols": gene_symbol,
                "condition_labels": condition.get("label"),
                "cache_condition_keys": f"{gene_symbol}|{condition_key}",
                "cache_condition_record_count": "1",
                "priva_scope_decisions": scope.get("decision"),
                "disease_scope_categories": scope.get("category"),
                "scope_review_statuses": scope.get("review_statuses")
                or scope.get("review_status"),
                "scope_evidence": scope.get("evidence"),
                "scope_references": scope.get("references"),
                "inheritance_modes": (condition.get("inheritance") or {}).get("modes"),
                "penetrance_statuses": (condition.get("penetrance") or {}).get("statuses"),
                "pathogenic_mechanisms": list(
                    (condition.get("pathogenic_mechanisms") or {}).keys()
                ),
            }
            for mechanism, block in (condition.get("pathogenic_mechanisms") or {}).items():
                for evidence in (block or {}).get("evidence", []) or []:
                    if not isinstance(evidence, dict):
                        continue
                    source = _clean(evidence.get("source"))
                    record_id = _clean(evidence.get("source_record_id"))
                    if source:
                        values.setdefault("source_names", []).append(source)
                    if source == "G2P_DDG2P":
                        values.setdefault("g2p_record_ids", []).append(record_id)
                    elif source == "Orphadata":
                        values.setdefault("orphadata_record_ids", []).append(record_id)
                    elif source == "ClinGen_haploinsufficiency":
                        values.setdefault("clingen_hi_disease_ids", []).append(record_id)
                    elif source == "GoFCards_exact+ClinVar_VCV":
                        values.setdefault("clinvar_vcv_accessions", []).extend(
                            _split_values(record_id)
                        )
                for variant_key, variant in (block or {}).get("variants", {}).items():
                    if not isinstance(variant, dict):
                        continue
                    values.setdefault("gofcards_variant_ids", []).append(
                        _clean(variant.get("gofcards_variant_id")) or variant_key
                    )
                    for link in variant.get("clinvar_links", []) or []:
                        if not isinstance(link, dict):
                            continue
                        values.setdefault("clinvar_vcv_accessions", []).append(
                            link.get("vcv_accession")
                        )
                        values.setdefault("clinvar_rcv_accessions", []).append(
                            link.get("rcv_accession")
                        )
            builder.add(
                identifiers,
                synthetic_key=f"CACHE:{gene_symbol}:{condition_key}",
                source_name="PriVA_condition_cache",
                **values,
            )
    expected = int((meta.get("counts") or {}).get("conditions") or 0)
    if expected and counts["condition_records"] != expected:
        raise ValueError(
            f"cache count mismatch: observed {counts['condition_records']}, metadata {expected}"
        )
    return dict(counts), meta


def _tokens(value: str) -> set[str]:
    return {token for token in value.split(";") if token}


def build_candidates(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    candidates: list[dict[str, str]] = []
    for row in rows:
        decisions = {value.lower() for value in _tokens(row["priva_scope_decisions"])}
        categories = {value.lower() for value in _tokens(row["disease_scope_categories"])}
        inheritance = {value.lower() for value in _tokens(row["inheritance_modes"])}
        scope_evidence = _tokens(row["scope_evidence"])
        raw_modes = {value.lower() for value in _tokens(row["source_pathogenic_modes"])}

        positive: set[str] = set()
        if "complex_or_non_monogenic" in categories:
            positive.add("scope_category:complex_or_non_monogenic")
        if "somatic_oncogenic" in categories:
            positive.add("scope_category:somatic_oncogenic")
        positive.update(
            f"scope_evidence:{value}"
            for value in scope_evidence
            if value == "HPO_non_monogenic_inheritance"
            or value.startswith("ORPHADATA_exact_gene_association:Disease_causing_somatic")
            or value.startswith("ORPHADATA_exact_gene_association:Part_of_a_fusion_gene")
        )
        positive.update(
            f"inheritance:{value}" for value in inheritance & NON_MENDELIAN_MODES
        )
        positive.update(
            f"source_mode:{value}"
            for value in raw_modes
            if "somatic" in value or "fusion gene" in value
        )

        counterevidence: set[str] = set()
        counterevidence.update(
            f"scope_decision:{value}" for value in decisions if value == "include"
        )
        counterevidence.update(
            f"inheritance:{value}" for value in inheritance & MENDELIAN_MODES
        )
        counterevidence.update(
            f"scope_evidence:{value}"
            for value in scope_evidence
            if value == "HPO_mendelian_inheritance"
            or value.startswith("ORPHADATA_exact_gene_association:Disease_causing_germline")
        )

        candidate_class = ""
        recommended_gate = ""
        if "exclude" in decisions:
            candidate_class = "explicit_scope_exclusion"
            recommended_gate = "exclude" if "include" not in decisions and positive else "manual_review"
        elif inheritance & NON_MENDELIAN_MODES:
            candidate_class = (
                "mixed_mendelian_non_mendelian_inheritance"
                if inheritance & MENDELIAN_MODES
                else "non_mendelian_inheritance_annotation"
            )
            recommended_gate = "manual_review"
        elif "review" in decisions:
            candidate_class = "unresolved_scope_review"
            recommended_gate = "manual_review"
        if not candidate_class:
            continue
        candidate = {
            "candidate_class": candidate_class,
            "recommended_cache_gate": recommended_gate,
            "positive_non_mendelian_evidence": ";".join(sorted(positive)),
            "mendelian_counterevidence": ";".join(sorted(counterevidence)),
            **row,
        }
        candidates.append(candidate)
    return sorted(
        candidates,
        key=lambda row: (row["recommended_cache_gate"], row["candidate_class"], row["canonical_condition_id"]),
    )


def _write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--condition-cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--hpo-assertions", type=Path, default=DEFAULT_HPO)
    parser.add_argument("--mechanism-evidence", type=Path, default=DEFAULT_MECHANISM_EVIDENCE)
    parser.add_argument("--scope-registry", type=Path, default=DEFAULT_SCOPE_REGISTRY)
    parser.add_argument("--clingen-dosage", type=Path, default=DEFAULT_CLINGEN)
    parser.add_argument("--gofcards", type=Path, default=DEFAULT_GOFCARDS)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    inputs = {
        "condition_cache": args.condition_cache,
        "hpo_assertions": args.hpo_assertions,
        "mechanism_evidence": args.mechanism_evidence,
        "scope_registry": args.scope_registry,
        "clingen_dosage": args.clingen_dosage,
        "gofcards": args.gofcards,
    }
    missing = [str(path) for path in inputs.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing audit inputs: " + ", ".join(missing))

    builder = InventoryBuilder()
    source_statistics = {
        "scope_registry": load_scope_registry(builder, args.scope_registry),
        "hpo_assertions": load_hpo_assertions(builder, args.hpo_assertions),
        "mechanism_evidence": load_mechanism_evidence(builder, args.mechanism_evidence),
        "clingen_dosage": load_clingen_dosage(builder, args.clingen_dosage),
        "gofcards_clinvar": load_gofcards_clinvar(builder, args.gofcards),
    }
    cache_statistics, cache_meta = load_condition_cache(builder, args.condition_cache)
    source_statistics["condition_cache"] = cache_statistics

    inventory = builder.collapse()
    candidates = build_candidates(inventory)
    inventory_path = args.output_dir / "condition_accession_inventory.tsv"
    candidate_path = args.output_dir / "non_mendelian_condition_review_candidates.tsv"
    summary_path = args.output_dir / "condition_scope_inventory_summary.json"
    _write_tsv(inventory_path, INVENTORY_COLUMNS, inventory)
    _write_tsv(candidate_path, CANDIDATE_COLUMNS, candidates)

    candidate_counts = Counter(row["candidate_class"] for row in candidates)
    gate_counts = Counter(row["recommended_cache_gate"] for row in candidates)
    multi_mondo = [row["canonical_condition_id"] for row in inventory if ";" in row["mondo_ids"]]
    unresolved = [row["canonical_condition_id"] for row in inventory if row["unresolved_source_keys"]]
    source_component_counts = Counter()
    for row in inventory:
        for source in _tokens(row["source_names"]):
            source_component_counts[source] += 1

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            name: {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for name, path in inputs.items()
        },
        "cache_meta": cache_meta,
        "source_statistics": source_statistics,
        "inventory": {
            "unique_identity_components": len(inventory),
            "components_by_source": dict(sorted(source_component_counts.items())),
            "components_with_multiple_mondo_ids": len(multi_mondo),
            "multiple_mondo_examples": multi_mondo[:25],
            "components_with_unresolved_source_keys": len(unresolved),
        },
        "candidates": {
            "rows": len(candidates),
            "by_class": dict(sorted(candidate_counts.items())),
            "by_recommended_gate": dict(sorted(gate_counts.items())),
        },
        "outputs": {
            "inventory": str(inventory_path.resolve()),
            "candidates": str(candidate_path.resolve()),
        },
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"inventory_rows\t{len(inventory)}")
    print(f"candidate_rows\t{len(candidates)}")
    print(f"recommended_exclusions\t{gate_counts.get('exclude', 0)}")
    print(f"manual_review\t{gate_counts.get('manual_review', 0)}")
    print(f"inventory_path\t{inventory_path.resolve()}")
    print(f"candidate_path\t{candidate_path.resolve()}")
    print(f"summary_path\t{summary_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
