#!/usr/bin/env python3
"""Build disease-specific gene-HPO assertions without gene-wide collapsing.

HPO's ``genes_to_phenotype.txt`` supplies the gene-to-disease relationship,
while ``phenotype.hpoa`` supplies assertion-level frequency and provenance.
Both files must come from the same HPO release.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from build_mondo_disease_scope import annotate_hpo_assertions


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
if not logger.handlers:
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s"
    )
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)


GENE_INPUT_COLUMNS = {
    "gene_symbol",
    "disease_id",
    "hpo_id",
    "frequency",
}
HPOA_INPUT_COLUMNS = {
    "database_id",
    "qualifier",
    "hpo_id",
    "frequency",
    "evidence",
    "reference",
}
ASSERTION_COLUMNS = [
    "gene_symbol",
    "disease_id",
    "hpo_id",
    "frequency",
    "evidence",
    "reference",
]


def _read_tsv(path: str | Path, *, comments: bool = False) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
        comment="#" if comments else None,
    )


def _require_columns(
    df: pd.DataFrame,
    required: set[str],
    source_name: str,
) -> None:
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{source_name} is missing columns: {sorted(missing)}")


def _normalize_optional(series: pd.Series) -> pd.Series:
    cleaned = series.astype(str).str.strip()
    return cleaned.mask(cleaned.eq(""), "-")


def build_gene_disease_hpo_assertions(
    genes_to_phenotype: str | Path,
    phenotype_hpoa: str | Path,
    *,
    minimum_match_rate: float = 0.999,
) -> pd.DataFrame:
    """Return one row per gene, disease, HPO, frequency, and provenance assertion."""
    gene_df = _read_tsv(genes_to_phenotype)
    hpoa_df = _read_tsv(phenotype_hpoa, comments=True)
    _require_columns(gene_df, GENE_INPUT_COLUMNS, "genes_to_phenotype")
    _require_columns(hpoa_df, HPOA_INPUT_COLUMNS, "phenotype.hpoa")

    gene_df = gene_df.loc[:, sorted(GENE_INPUT_COLUMNS)].copy()
    hpoa_df = hpoa_df.loc[:, sorted(HPOA_INPUT_COLUMNS)].copy()
    for column in gene_df.columns:
        gene_df[column] = gene_df[column].astype(str).str.strip()
    for column in hpoa_df.columns:
        hpoa_df[column] = hpoa_df[column].astype(str).str.strip()

    invalid_identifier = gene_df[["disease_id", "hpo_id"]].isin(["", "-"]).any(axis=1)
    if invalid_identifier.any():
        raise ValueError(
            f"genes_to_phenotype contains {int(invalid_identifier.sum())} rows with "
            "an empty disease_id or hpo_id"
        )
    unmapped_gene = gene_df["gene_symbol"].isin(["", "-"])
    if unmapped_gene.any():
        logger.info(
            "Excluded %d HPO rows without a mapped gene symbol",
            int(unmapped_gene.sum()),
        )
        gene_df = gene_df.loc[~unmapped_gene].copy()

    # genes_to_phenotype does not carry HPOA's qualifier. Identify NOT-only
    # disease/HPO links so they can be excluded rather than turned into false
    # positive assertions with blank provenance.
    negative_pairs = set(
        map(
            tuple,
            hpoa_df.loc[
                hpoa_df["qualifier"].eq("NOT"), ["database_id", "hpo_id"]
            ].itertuples(index=False, name=None),
        )
    )
    positive_hpoa = hpoa_df.loc[hpoa_df["qualifier"].eq("")].copy()
    positive_hpoa.rename(columns={"database_id": "disease_id"}, inplace=True)
    positive_hpoa["frequency"] = _normalize_optional(positive_hpoa["frequency"])
    positive_hpoa["evidence"] = _normalize_optional(positive_hpoa["evidence"])
    positive_hpoa["reference"] = _normalize_optional(positive_hpoa["reference"])

    hpoa_assertions = positive_hpoa.loc[
        :, ["disease_id", "hpo_id", "frequency", "evidence", "reference"]
    ].drop_duplicates()
    gene_links = gene_df.loc[
        :, ["gene_symbol", "disease_id", "hpo_id", "frequency"]
    ].rename(columns={"frequency": "gene_summary_frequency"})

    joined = gene_links.merge(
        hpoa_assertions,
        on=["disease_id", "hpo_id"],
        how="left",
        indicator=True,
        validate="many_to_many",
    )

    matched = joined["_merge"].eq("both")
    joined_pairs = list(zip(joined["disease_id"], joined["hpo_id"]))
    negative_only = ~matched & pd.Series(
        (pair in negative_pairs for pair in joined_pairs),
        index=joined.index,
    )
    unknown = ~matched & ~negative_only
    match_rate = 1.0 - (float(unknown.sum()) / len(joined) if len(joined) else 0.0)
    if match_rate < minimum_match_rate:
        raise ValueError(
            "Only "
            f"{match_rate:.3%} of gene-disease-HPO rows matched phenotype.hpoa; "
            "the source files are probably from different HPO releases"
        )

    if negative_only.any():
        logger.info(
            "Excluded %d gene-disease-HPO rows supported only by qualifier=NOT",
            int(negative_only.sum()),
        )

    unmatched = joined.loc[unknown, ["gene_symbol", "disease_id", "hpo_id"]]
    if not unmatched.empty:
        examples = unmatched.drop_duplicates().head(5).to_dict(orient="records")
        logger.warning(
            "%d gene-disease-HPO rows have no positive phenotype.hpoa assertion; "
            "retaining them with '-' evidence/reference. Examples: %s",
            len(unmatched),
            examples,
        )

    joined["frequency"] = joined["frequency"].where(
        matched,
        _normalize_optional(joined["gene_summary_frequency"]),
    )
    joined["evidence"] = _normalize_optional(joined["evidence"])
    joined["reference"] = _normalize_optional(joined["reference"])

    assertions = (
        joined.loc[~negative_only, ASSERTION_COLUMNS]
        .drop_duplicates()
        .sort_values(ASSERTION_COLUMNS, kind="stable")
        .reset_index(drop=True)
    )
    logger.info(
        "Built %d unique assertions for %d genes and %d diseases (match rate %.3f%%)",
        len(assertions),
        assertions["gene_symbol"].nunique(),
        assertions["disease_id"].nunique(),
        match_rate * 100,
    )
    return assertions


def write_gene_disease_hpo_assertions(
    genes_to_phenotype: str | Path,
    phenotype_hpoa: str | Path,
    output_file: str | Path,
    *,
    minimum_match_rate: float = 0.999,
    disease_scope_registry: str | Path | None = None,
) -> pd.DataFrame:
    assertions = build_gene_disease_hpo_assertions(
        genes_to_phenotype,
        phenotype_hpoa,
        minimum_match_rate=minimum_match_rate,
    )
    if disease_scope_registry is not None:
        registry = _read_tsv(disease_scope_registry)
        assertions = annotate_hpo_assertions(assertions, registry)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    assertions.to_csv(output_path, sep="\t", index=False)
    logger.info("Wrote assertion table to %s", output_path)
    return assertions


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Join release-matched HPO gene links and HPOA provenance."
    )
    parser.add_argument("genes_to_phenotype", help="HPO genes_to_phenotype.txt")
    parser.add_argument("phenotype_hpoa", help="Matching HPO phenotype.hpoa")
    parser.add_argument("output", help="Output assertion TSV; .gz is supported")
    parser.add_argument(
        "--minimum-match-rate",
        type=float,
        default=0.999,
        help="Fail below this gene/HPO source match rate (default: 0.999)",
    )
    parser.add_argument(
        "--disease-scope-registry",
        help="Optional MONDO/HPO disease-scope registry to attach to each assertion",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    write_gene_disease_hpo_assertions(
        args.genes_to_phenotype,
        args.phenotype_hpoa,
        args.output,
        minimum_match_rate=args.minimum_match_rate,
        disease_scope_registry=args.disease_scope_registry,
    )
