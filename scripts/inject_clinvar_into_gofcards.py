#!/usr/bin/env python3
"""Attach ClinVar VCV condition evidence to each GoFCards variant.

This runs immediately after normalization and before anything downstream. It
reads the normalized GoFCards cache, matches each eligible variant against
ClinVar by exact genomic allele, and writes the same structure back with a
``clinvar`` block nested under the variants that matched.

Why the annotation belongs on the variant rather than beside it: a ClinVar
record is evidence *about one variant*. Storing it as a sibling entry forces
every reader to re-resolve which variant it belongs to, which is what the
previous token-matching indirection did.

What ClinVar does and does not contribute:

  * conditions  -- yes. This is what the HPO condition cache is built from.
  * mechanism   -- NO. ClinVar pathogenicity establishes neither gain of
                   function nor a dominant-negative effect, so the variant
                   matcher ignores the ``clinvar`` block entirely.

Matching is by exact normalized allele on either build. A variant that failed
liftover is still matchable through its GRCh37 coordinates.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

from clinvar_vcv import (  # noqa: E402
    flatten_match_audit,
    gofcards_genomic_index_key,
    gofcards_variation_id,
    gofcards_variation_id_index_key,
    index_gofcards_variants,
    iter_gofcards_variants,
    load_gofcards_cache,
    open_text,
    stream_parse_clinvar_vcv,
)


def log(message: str) -> None:
    print(f"[clinvar-inject] {message}", file=sys.stderr, flush=True)


def _clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def clinvar_block(payload: dict[str, Any]) -> dict[str, Any]:
    """Reduce one ClinVar VCV payload to what a condition consumer needs."""
    variation = payload.get("variation") or {}
    hgvs = [
        _clean(entry.get("expression"))
        for entry in variation.get("hgvs", []) or []
        if isinstance(entry, dict) and _clean(entry.get("expression"))
    ]
    return {
        "vcv_accession": _clean(variation.get("vcv_accession")),
        "variation_id": _clean(variation.get("variation_id")),
        "variation_name": _clean(variation.get("name")),
        "hgvs": list(dict.fromkeys(hgvs)),
        "condition_assertions": payload.get("condition_assertions") or [],
    }


def inject(cache: dict[str, Any], matches: list[dict[str, Any]]) -> dict[str, int]:
    """Nest each ClinVar match under the GoFCards variant it matched.

    The link is by variant identifier, which the lookup already carries, so no
    token re-resolution is needed here or in any downstream reader.
    """
    variants_by_id = {
        variant_id: variant for _symbol, variant_id, variant in iter_gofcards_variants(cache)
    }
    counts: Counter[str] = Counter()
    for match in matches:
        payload = match.get("ClinVar_VCV") or {}
        block = clinvar_block(payload)
        if not block["condition_assertions"]:
            counts["clinvar_records_without_conditions"] += 1
            continue
        matched = (payload.get("match") or {}).get("matched_gofcards_records") or []
        for record in matched:
            if not isinstance(record, dict):
                continue
            variant_id = _clean(record.get("gofcards_variant_id"))
            variant = variants_by_id.get(variant_id)
            if variant is None:
                counts["clinvar_matches_without_a_variant"] += 1
                continue
            existing = variant.get("clinvar")
            if existing is None:
                # Record every route that confirms this link, not just the first.
                # Both routes usually agree, and saying so is more informative
                # than picking one and implying the other did not apply.
                routes = {record.get("link_route", "")} - {""}
                if gofcards_variation_id(variant) == block["variation_id"] != "":
                    routes.add("gofcards_variation_id")
                variant["clinvar"] = {**block, "link_routes": sorted(routes)}
                counts["variants_annotated"] += 1
                for route in routes:
                    counts[f"confirmed_by_{route}"] += 1
                if len(routes) > 1:
                    counts["confirmed_by_both_routes"] += 1
                counts["condition_assertions"] += len(block["condition_assertions"])
            elif existing.get("vcv_accession") != block["vcv_accession"]:
                # One variant can appear in more than one ClinVar record; keep
                # them all rather than letting the last one win.
                extra = variant.setdefault("clinvar_additional", [])
                if block["vcv_accession"] not in {e.get("vcv_accession") for e in extra}:
                    extra.append(block)
                    counts["additional_clinvar_records"] += 1
    return dict(counts)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--gofcards-cache", required=True, type=Path,
                        help="normalized GoFCards cache produced by build_gofcards_exact_gof_cache.py")
    parser.add_argument("--clinvar-vcv", required=True, type=Path,
                        help="ClinVar VariationArchive XML (optionally gzipped)")
    parser.add_argument("--out-json", required=True, type=Path,
                        help="gzipped JSON with a clinvar block nested per matched variant")
    parser.add_argument("--min-review-stars", type=int, default=2, choices=(2, 3, 4))
    parser.add_argument("--audit-tsv", type=Path, default=None,
                        help="per-condition audit of every ClinVar match")
    parser.add_argument("--stats-json", type=Path, default=None)
    args = parser.parse_args(argv)

    if not args.clinvar_vcv.exists() or args.clinvar_vcv.stat().st_size == 0:
        raise SystemExit(f"ClinVar VCV is required and was not found: {args.clinvar_vcv}")

    # Parsed once. Both indexes are built over that same parsed object, so the
    # file is not read three times and the two indexes point at the very same
    # variant objects the injection will write onto.
    cache = load_gofcards_cache(args.gofcards_cache)
    lookup = index_gofcards_variants(cache, gofcards_genomic_index_key)
    by_variation_id = index_gofcards_variants(cache, gofcards_variation_id_index_key)
    log(f"indexed {len(by_variation_id)} GoFCards ClinVar variation ids "
        f"and {len(lookup)} genomic keys from eligible variants")

    matches, parse_stats = stream_parse_clinvar_vcv(
        args.clinvar_vcv, lookup, min_review_stars=args.min_review_stars,
        variation_id_lookup=by_variation_id,
    )
    log(f"ClinVar returned {len(matches)} matched records")

    counts = inject(cache, matches)
    cache.setdefault("metadata", {})["clinvar"] = {
        "source": str(args.clinvar_vcv),
        "min_review_stars": args.min_review_stars,
        **{k: int(v) for k, v in parse_stats.items()},
        **counts,
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.out_json, "wt") as handle:
        json.dump(cache, handle, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    log("injected: " + "; ".join(f"{k}={v}" for k, v in counts.items()))

    if args.audit_tsv:
        # Written here rather than downstream: this step holds the full ClinVar
        # match objects, so the audit does not have to be reconstructed from a
        # reduced copy.
        columns = [
            "gene_symbol", "vcv_accession", "variation_id", "variation_name",
            "classification_scope", "matched_component_context", "rcv_accession",
            "condition", "clinical_significance", "review_status", "review_stars",
            "matched_scv_count",
        ]
        args.audit_tsv.parent.mkdir(parents=True, exist_ok=True)
        rows = flatten_match_audit(matches)
        with args.audit_tsv.open("w", encoding="utf-8", newline="") as handle:
            handle.write("\t".join(columns) + "\n")
            for row in rows:
                handle.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")
        log(f"wrote {len(rows)} audit rows to {args.audit_tsv}")

    if args.stats_json:
        args.stats_json.write_text(
            json.dumps({**parse_stats, **counts}, indent=2), encoding="utf-8"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
