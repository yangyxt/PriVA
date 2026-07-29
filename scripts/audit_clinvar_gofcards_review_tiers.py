#!/usr/bin/env python3
"""Audit ClinVar RCV review tiers among exact GoFCards allele matches."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from clinvar_vcv import (
    gofcards_genomic_index_key,
    index_gofcards_variants,
    load_gofcards_cache,
    stream_parse_clinvar_vcv,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xml", type=Path, required=True)
    parser.add_argument("--gofcards-exact", type=Path, required=True)
    parser.add_argument("--summary-json", type=Path, required=True)
    parser.add_argument("--rcv-tsv", type=Path, required=True)
    return parser.parse_args()


def summarize_review_tiers(
    matches: list[dict[str, Any]], parser_stats: dict[str, int]
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rcv_rows: list[dict[str, Any]] = []
    vcv_stars: dict[str, list[int]] = defaultdict(list)
    vcv_genes: dict[str, set[str]] = defaultdict(set)
    review_status_counts: Counter[str] = Counter()
    star_counts: Counter[int] = Counter()

    seen_rcv_rows: set[tuple[str, str, str]] = set()
    for match in matches:
        payload = match.get("ClinVar_VCV", {})
        variation = payload.get("variation", {})
        vcv = str(variation.get("vcv_accession") or "")
        genes = sorted(str(gene) for gene in match.get("gene_symbols", []) if gene)
        vcv_genes[vcv].update(genes)
        for condition in payload.get("condition_assertions", []):
            classification = condition.get("germline_classification", {})
            rcv = str(condition.get("rcv_accession") or "")
            status = str(classification.get("review_status") or "")
            stars = int(classification.get("review_stars") or 0)
            dedup_key = (vcv, rcv, status)
            if dedup_key in seen_rcv_rows:
                continue
            seen_rcv_rows.add(dedup_key)
            vcv_stars[vcv].append(stars)
            review_status_counts[status or "<empty>"] += 1
            star_counts[stars] += 1
            rcv_rows.append(
                {
                    "vcv_accession": vcv,
                    "rcv_accession": rcv,
                    "gene_symbols": ";".join(genes),
                    "review_stars": stars,
                    "review_status": status,
                    "clinical_significance": str(
                        classification.get("clinical_significance") or ""
                    ),
                    "conditions": json.dumps(
                        condition.get("conditions", []),
                        sort_keys=True,
                    ),
                }
            )

    vcv_max_star_counts: Counter[int] = Counter(
        max(stars) for stars in vcv_stars.values() if stars
    )
    vcvs_with_two_plus = {
        vcv for vcv, stars in vcv_stars.items() if max(stars, default=0) >= 2
    }
    vcvs_below_two = set(vcv_stars) - vcvs_with_two_plus
    genes_with_two_plus = {
        gene for vcv in vcvs_with_two_plus for gene in vcv_genes.get(vcv, set())
    }
    genes_below_two = {
        gene for vcv in vcvs_below_two for gene in vcv_genes.get(vcv, set())
    }

    summary = {
        "parser_stats_minimum_review_stars_0": parser_stats,
        "gene_concordant_vcvs_with_any_germline_rcv": len(vcv_stars),
        "gene_concordant_vcvs_with_any_rcv_at_least_2_stars": len(vcvs_with_two_plus),
        "gene_concordant_vcvs_with_rcvs_but_none_at_least_2_stars": len(
            vcvs_below_two
        ),
        "genes_with_any_vcv_at_least_2_stars": len(genes_with_two_plus),
        "genes_with_only_lower_review_matched_vcvs": len(
            genes_below_two - genes_with_two_plus
        ),
        "genes_with_any_lower_review_matched_vcv": len(genes_below_two),
        "rcv_review_star_counts": {
            str(key): value for key, value in sorted(star_counts.items())
        },
        "vcv_maximum_rcv_review_star_counts": {
            str(key): value for key, value in sorted(vcv_max_star_counts.items())
        },
        "rcv_review_status_counts": dict(sorted(review_status_counts.items())),
    }
    return summary, sorted(
        rcv_rows,
        key=lambda row: (
            row["vcv_accession"],
            row["rcv_accession"],
            row["review_stars"],
        ),
    )


def main() -> None:
    args = parse_args()
    lookup = index_gofcards_variants(
        load_gofcards_cache(args.gofcards_exact), gofcards_genomic_index_key
    )
    matches, parser_stats = stream_parse_clinvar_vcv(
        args.xml,
        lookup,
        min_review_stars=0,
    )
    summary, rcv_rows = summarize_review_tiers(matches, parser_stats)
    summary["clinvar_xml"] = str(args.xml.resolve())
    summary["gofcards_exact"] = str(args.gofcards_exact.resolve())

    args.summary_json.parent.mkdir(parents=True, exist_ok=True)
    args.summary_json.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    args.rcv_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.rcv_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rcv_rows[0]) if rcv_rows else [],
            delimiter="\t",
        )
        if rcv_rows:
            writer.writeheader()
            writer.writerows(rcv_rows)

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"rcv_tsv={args.rcv_tsv.resolve()}")


if __name__ == "__main__":
    main()
