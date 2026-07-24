#!/usr/bin/env python3
"""Re-sort existing .filtered.tsv files using sort_and_rank_variants only (skip ACMG)."""

import argparse
import os
import sys
import logging

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from acmg_criteria_assign import sort_and_rank_variants

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def resort_family(filtered_tsv: str, ped_table: str,
                  fam_name: str,
                  pext_tissues: str, relevant_gene_list: str,
                  dispensable_gene_list: str,
                  pext_low_expression_cutoff: float,
                  pext_penalty_floor: float,
                  pext_penalty_shape: float):
    if not os.path.isfile(filtered_tsv):
        logger.warning(f"SKIP: {filtered_tsv} not found")
        return
    anno_df = pd.read_table(filtered_tsv, low_memory=False)
    if anno_df.empty or "ACMG_quant_score" not in anno_df.columns:
        logger.warning(f"SKIP: {filtered_tsv} has no ACMG_quant_score")
        return

    ped_df = pd.read_table(ped_table, low_memory=False) if ped_table else None

    # Drop columns that sort_and_rank_variants will recreate
    drop_cols = [c for c in ("sort_index", "pext_sort_index", "variant_rank",
                             "control_common_index", "haplo_insuf_index",
                             "zygosity_inheritance_mechanism_index",
                             "zygosity_inheritance_mechanism_compatibility",
                             "dispensable_gene_index", "relevant_gene_index",
                             "max_variant_score") if c in anno_df.columns]
    anno_df = anno_df.drop(columns=drop_cols)

    sorted_df = sort_and_rank_variants(
        anno_df, ped_df, fam_name,
        pext_tissues=pext_tissues,
        pext_low_expression_cutoff=pext_low_expression_cutoff,
        pext_penalty_floor=pext_penalty_floor,
        pext_penalty_shape=pext_penalty_shape,
        relevant_gene_list=relevant_gene_list,
        dispensable_gene_list=dispensable_gene_list,
    )
    sorted_df.to_csv(filtered_tsv, sep="\t", index=False)
    logger.info(f"DONE: {fam_name} → {len(sorted_df)} rows, ranks 1-{sorted_df['variant_rank'].max()}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Re-sort .filtered.tsv (step 3 only)")
    ap.add_argument("--filtered_tsv", required=True)
    ap.add_argument("--ped_table", default=None)
    ap.add_argument("--fam_name", default=None)
    ap.add_argument("--pext_tissues", default="")
    ap.add_argument("--pext_low_expression_cutoff", type=float, default=0.1)
    ap.add_argument("--pext_penalty_floor", type=float, default=0.8)
    ap.add_argument("--pext_penalty_shape", type=float, default=0.5)
    ap.add_argument("--relevant_gene_list", default=None)
    ap.add_argument("--dispensable_gene_list", default=None)
    args = ap.parse_args()
    resort_family(args.filtered_tsv, args.ped_table,
                  args.fam_name,
                  args.pext_tissues, args.relevant_gene_list,
                  args.dispensable_gene_list,
                  args.pext_low_expression_cutoff,
                  args.pext_penalty_floor,
                  args.pext_penalty_shape)
