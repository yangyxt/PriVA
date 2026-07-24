#!/usr/bin/env python3
"""Audit which exact GoFCards genes lack independent mechanism coverage.

This audit is intentionally identifier- and source-based. It never links
diseases by name, publication, phenotype similarity, or approximate variants.
Gene symbols are normalized through PriVA's pinned HGNC table before sets are
compared, so aliases do not create artificial source-specific genes.
"""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
from typing import Callable

import pandas as pd


GeneResolver = Callable[[str], str]
CONDITION_SOURCES = ("G2P_DDG2P", "Orphadata")
CANONICAL_MECHANISMS = {"LOF", "GOF", "DOMINANT_NEGATIVE", "TRIPLOSENSITIVITY"}


def load_gofcards_gene_counts(
    path: str | Path,
    *,
    resolve_symbol: GeneResolver,
) -> Counter[str]:
    """Count exact-cache GoFCards rows per HGNC-normalized gene.

    The exact GoFCards cache is the relevant starting population because these
    are the variant records PriVA can actually match at runtime. Counting rows,
    rather than reducing immediately to a set, preserves a useful audit signal:
    a remaining gene can represent either one isolated allele or many curated
    GoFCards alleles.

    Blank or unresolvable symbols are excluded. ``resolve_symbol`` is injected
    so production uses PriVA's HGNC resolver while unit tests remain independent
    of the large local HGNC table.
    """
    frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False).fillna("")
    if "HGNC_Symbol" not in frame.columns:
        raise ValueError(f"{path} is missing required column: HGNC_Symbol")

    counts: Counter[str] = Counter()
    for raw_symbol in frame["HGNC_Symbol"]:
        raw_symbol = str(raw_symbol).strip()
        if not raw_symbol:
            continue
        symbol = str(resolve_symbol(raw_symbol)).strip()
        if symbol:
            counts[symbol] += 1
    return counts


def load_condition_mechanism_gene_sets(
    path: str | Path,
    *,
    resolve_symbol: GeneResolver,
) -> dict[str, dict[str, set[str]]]:
    """Return broad, explicit-mechanism, and PriVA-usable genes by source.

    Three nested sets answer different audit questions without conflation:

    * ``all_rows`` means the gene occurs anywhere in that source;
    * ``canonical_mechanism`` requires an explicit normalized GOF/LOF/DN/TS
      assertion, regardless of disease-scope disposition; and
    * ``priva_included_mechanism`` additionally requires ``priva_scope=include``
      and is therefore eligible for automatic condition-history transfer.

    This distinction matters for review or excluded diseases. Their existence
    demonstrates source coverage, but it must not be presented as an inherited
    condition mechanism that PriVA can automatically apply.
    """
    frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False).fillna("")
    required = {
        "gene_symbol",
        "source",
        "normalized_mechanisms",
        "priva_scope",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")

    result = {
        source: {
            "all_rows": set(),
            "canonical_mechanism": set(),
            "priva_included_mechanism": set(),
        }
        for source in CONDITION_SOURCES
    }
    for row in frame.loc[frame["source"].isin(CONDITION_SOURCES)].to_dict(
        orient="records"
    ):
        raw_symbol = str(row["gene_symbol"]).strip()
        if not raw_symbol:
            continue
        symbol = str(resolve_symbol(raw_symbol)).strip()
        if not symbol:
            continue

        source = row["source"]
        result[source]["all_rows"].add(symbol)
        mechanisms = {
            token.strip().upper()
            for token in str(row["normalized_mechanisms"]).replace("|", ";").split(";")
            if token.strip()
        }
        if not mechanisms & CANONICAL_MECHANISMS:
            continue
        result[source]["canonical_mechanism"].add(symbol)
        if str(row["priva_scope"]).strip().lower() == "include":
            result[source]["priva_included_mechanism"].add(symbol)
    return result


def load_exact_clinvar_linked_gene_counts(
    path: str | Path,
    *,
    resolve_symbol: GeneResolver,
) -> Counter[str]:
    """Count ClinVar VCV records linked to GoFCards by an exact allele key.

    The shared mechanism JSON records the match method and the exact GoFCards
    rows supporting each link. Both must be present. Merely having a ClinVar
    record in the same gene, sharing a disease name, or sharing a publication is
    not counted as coverage of a GoFCards allele.

    Counts are VCV records, not unique variants or conditions. The final audit
    uses only whether the normalized gene has at least one exact link, while the
    count remains available to inspect the strength of that coverage.
    """
    with open(path, encoding="utf-8") as handle:
        payload = json.load(handle)

    counts: Counter[str] = Counter()
    for key, gene in payload.items():
        if key == "_meta" or not isinstance(gene, dict):
            continue
        raw_symbol = str(gene.get("symbol", "")).strip()
        if not raw_symbol:
            continue
        symbol = str(resolve_symbol(raw_symbol)).strip()
        if not symbol:
            continue

        for wrapper in gene.get("variant_level", []) or []:
            if not isinstance(wrapper, dict):
                continue
            record = wrapper.get("ClinVar_VCV")
            if not isinstance(record, dict):
                continue
            match = record.get("match") or {}
            if match.get("method") != "exact_normalized_vcf_allele":
                continue
            if not match.get("matched_gofcards_records"):
                continue
            counts[symbol] += 1
    return counts
