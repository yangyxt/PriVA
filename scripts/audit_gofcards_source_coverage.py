#!/usr/bin/env python3
"""Audit which exact GoFCards genes lack independent mechanism coverage.

This audit is intentionally identifier- and source-based. It never links
diseases by name, publication, phenotype similarity, or approximate variants.
Gene symbols are normalized through PriVA's pinned HGNC table before sets are
compared, so aliases do not create artificial source-specific genes.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Callable

import pandas as pd


GeneResolver = Callable[[str], str]


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
