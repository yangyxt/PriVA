#!/usr/bin/env python3
"""Audit which exact GoFCards genes lack independent mechanism coverage.

This audit is intentionally identifier- and source-based. It never links
diseases by name, publication, phenotype similarity, or approximate variants.
Gene symbols are normalized through PriVA's pinned HGNC table before sets are
compared, so aliases do not create artificial source-specific genes.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Callable

import pandas as pd


GeneResolver = Callable[[str], str]
CONDITION_SOURCES = ("G2P_DDG2P", "Orphadata")
CANONICAL_MECHANISMS = {"LOF", "GOF", "DOMINANT_NEGATIVE", "TRIPLOSENSITIVITY"}
SCRIPT_DIR = Path(__file__).resolve().parent
PRIVA_ROOT = SCRIPT_DIR.parent
DEFAULT_GOFCARDS = PRIVA_ROOT / "data" / "gofcards" / "gofcards_exact_gof_hgvsp.tsv.gz"
DEFAULT_EVIDENCE = (
    PRIVA_ROOT
    / "data"
    / "gene_pathogenic_mechanism"
    / "prepared"
    / "gene_pathogenic_mechanism_evidence.tsv"
)
DEFAULT_NONLOF_JSON = (
    PRIVA_ROOT
    / "data"
    / "gene_pathogenic_mechanism"
    / "prepared"
    / "gene_nonlof_mechanism_curated_assertions.json"
)
DEFAULT_HGNC = PRIVA_ROOT / "data" / "hgnc" / "non_alt_loci_set.tsv"
DEFAULT_OUTPUT = (
    PRIVA_ROOT / "data" / "gofcards" / "gofcards_gene_source_coverage.audit.tsv"
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse explicit, overridable resource paths for reproducible audits.

    Defaults point to the same deployed resources used by PriVA. Every path is
    exposed because archived releases should be auditable without editing the
    script or silently mixing current and historical source versions.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gofcards", type=Path, default=DEFAULT_GOFCARDS)
    parser.add_argument("--condition-evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument(
        "--nonlof-mechanism-json",
        "--mechanism-json",
        dest="mechanism_json",
        type=Path,
        default=DEFAULT_NONLOF_JSON,
        help="Shared non-LOF mechanism JSON; --mechanism-json is a compatibility alias.",
    )
    parser.add_argument("--hgnc-table", type=Path, default=DEFAULT_HGNC)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=None,
        help="JSON summary path (default: <output stem>.summary.json).",
    )
    return parser.parse_args(argv)


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


def build_gofcards_coverage_audit(
    gofcards_counts: Counter[str],
    clinvar_counts: Counter[str],
    condition_sets: dict[str, dict[str, set[str]]],
) -> pd.DataFrame:
    """Build one transparent source-coverage row per exact GoFCards gene.

    ``only_gofcards_vs_explicit_sources`` is the primary strict answer: no exact
    ClinVar allele link and no canonical mechanism assertion from either G2P or
    Orphadata. Two companion flags make the policy boundary inspectable:

    * ``only_gofcards_vs_priva_included`` treats review/excluded diseases as not
      usable because PriVA cannot automatically transfer their histories; and
    * ``only_gofcards_vs_any_source_record`` asks the broader bibliographic
      question of whether the gene appears in G2P/Orphadata at all.

    All flags are integers in the written table for stable TSV interoperability.
    The underlying source-specific columns prevent a composite flag from hiding
    which database supplied coverage.
    """
    g2p = condition_sets["G2P_DDG2P"]
    orphadata = condition_sets["Orphadata"]
    rows: list[dict[str, int | str]] = []
    for symbol in sorted(gofcards_counts):
        clinvar_rows = clinvar_counts.get(symbol, 0)
        g2p_any = symbol in g2p["all_rows"]
        g2p_mechanism = symbol in g2p["canonical_mechanism"]
        g2p_included = symbol in g2p["priva_included_mechanism"]
        orphadata_any = symbol in orphadata["all_rows"]
        orphadata_mechanism = symbol in orphadata["canonical_mechanism"]
        orphadata_included = symbol in orphadata["priva_included_mechanism"]

        explicit_coverage = bool(
            clinvar_rows or g2p_mechanism or orphadata_mechanism
        )
        included_coverage = bool(
            clinvar_rows or g2p_included or orphadata_included
        )
        any_source_coverage = bool(clinvar_rows or g2p_any or orphadata_any)
        rows.append(
            {
                "gene_symbol": symbol,
                "gofcards_exact_rows": gofcards_counts[symbol],
                "exact_clinvar_vcv_rows": clinvar_rows,
                "g2p_any_source_record": int(g2p_any),
                "g2p_canonical_mechanism": int(g2p_mechanism),
                "g2p_priva_included_mechanism": int(g2p_included),
                "orphadata_any_source_record": int(orphadata_any),
                "orphadata_canonical_mechanism": int(orphadata_mechanism),
                "orphadata_priva_included_mechanism": int(orphadata_included),
                "only_gofcards_vs_explicit_sources": int(not explicit_coverage),
                "only_gofcards_vs_priva_included": int(not included_coverage),
                "only_gofcards_vs_any_source_record": int(
                    not any_source_coverage
                ),
            }
        )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    """Run the normalized cross-source audit and write TSV plus JSON outputs.

    Output creation is intentionally separate from each loader. This keeps the
    comparison reusable in tests and notebooks, while the command-line entry
    point records the exact input paths alongside aggregate counts. The HGNC
    resolver is the same local resolver used by the runtime mechanism hub.
    """
    args = parse_args(argv)

    # Import lazily so the pure loader/classifier functions can be reused without
    # initializing any of the runtime hub's optional resources.
    from gene_mechanism_hub import GeneMechanismHub

    hub = GeneMechanismHub(
        hgnc_table=args.hgnc_table,
        use_hgnc_package=False,
    )
    gofcards_counts = load_gofcards_gene_counts(
        args.gofcards,
        resolve_symbol=hub.resolve_symbol,
    )
    condition_sets = load_condition_mechanism_gene_sets(
        args.condition_evidence,
        resolve_symbol=hub.resolve_symbol,
    )
    clinvar_counts = load_exact_clinvar_linked_gene_counts(
        args.mechanism_json,
        resolve_symbol=hub.resolve_symbol,
    )
    audit = build_gofcards_coverage_audit(
        gofcards_counts,
        clinvar_counts,
        condition_sets,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.output, sep="\t", index=False)
    summary_output = args.summary_output or args.output.with_name(
        f"{args.output.stem}.summary.json"
    )
    summary_output.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "inputs": {
            "gofcards": str(args.gofcards.resolve()),
            "condition_evidence": str(args.condition_evidence.resolve()),
            "nonlof_mechanism_json": str(args.mechanism_json.resolve()),
            "hgnc_table": str(args.hgnc_table.resolve()),
        },
        "outputs": {
            "gene_audit_tsv": str(args.output.resolve()),
            "summary_json": str(summary_output.resolve()),
        },
        "gofcards_exact_rows": int(sum(gofcards_counts.values())),
        "gofcards_exact_genes": len(gofcards_counts),
        "genes_with_exact_clinvar_link": int(
            audit["exact_clinvar_vcv_rows"].astype(bool).sum()
        ),
        "genes_with_g2p_canonical_mechanism": int(
            audit["g2p_canonical_mechanism"].astype(bool).sum()
        ),
        "genes_with_orphadata_canonical_mechanism": int(
            audit["orphadata_canonical_mechanism"].astype(bool).sum()
        ),
        "genes_only_gofcards_vs_explicit_sources": int(
            audit["only_gofcards_vs_explicit_sources"].astype(bool).sum()
        ),
        "genes_only_gofcards_vs_priva_included": int(
            audit["only_gofcards_vs_priva_included"].astype(bool).sum()
        ),
        "genes_only_gofcards_vs_any_source_record": int(
            audit["only_gofcards_vs_any_source_record"].astype(bool).sum()
        ),
    }
    with open(summary_output, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
