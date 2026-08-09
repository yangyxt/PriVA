#!/usr/bin/env python3
"""PM1, PP2 and BP1 -- is a missense change tolerated here?

Three scales of the same question.

    PM1  the residue sits in a mutational hotspot or a well-studied
         functional domain without benign variation. Three arms: an
         AlphaMissense-derived intolerant domain, a HCSeeker hotspot, and a
         gnomAD regional missense constraint region.
    PP2  missense is a common mechanism of disease in this gene and benign
         missense variation is rare -- so a missense change is itself weak
         evidence for pathogenicity.
    BP1  the mirror image: truncating variants are the known mechanism, so a
         missense change is weak evidence against.

analyze_bp1_pp2 derives the per-gene pathogenic-to-benign missense ratio that
PP2 and BP1 both threshold.
"""

import logging
import pickle
import gzip
import multiprocessing as mp
import pandas as pd
import numpy as np
from typing import Tuple, Dict
from scipy import stats


logger = logging.getLogger(__name__)


def locate_pm1_region(row: dict,
                      hcseeker_hotspots: dict,
                      rmc_constrained: dict) -> Tuple[bool, bool]:
    """Check if the variant falls in an HCSeeker hotspot or gnomAD RMC-constrained region.

    Returns (in_hcseeker, in_rmc).
    Both lookups are keyed by Ensembl transcript (Feature) + protein position.
    """
    enst = row.get('Feature', '')
    raw_protein_pos = row.get('Protein_position', '')
    pos_str = str(raw_protein_pos).split("/")[0] if raw_protein_pos and raw_protein_pos not in [np.nan, np.inf] else ''

    if not pos_str or not pos_str.isdigit():
        return False, False

    pos_int = int(pos_str)
    in_hcseeker = enst in hcseeker_hotspots and pos_int in hcseeker_hotspots[enst]
    in_rmc = enst in rmc_constrained and pos_int in rmc_constrained[enst]
    return in_hcseeker, in_rmc


def PM1_criteria(df: pd.DataFrame,
                 pvs1_criteria: np.ndarray,
                 loc_intol_domain: np.ndarray,
                 clinvar_patho: np.ndarray = None,
                 pm1_regions_pkl: str=None,
                 threads: int = 10) -> np.ndarray:
    """Assign PM1 criteria using the 3-arm PM1 regions pickle.

    Arms:
      1. DAS intolerant domains — via loc_intol_domain (pre-computed in PVS1)
      2. HCSeeker hotspots — transcript + protein position lookup
      3. gnomAD RMC constrained — transcript + protein position lookup
    """
    row_dicts = df.to_dict('records')

    logger.info(f"There are {np.sum(loc_intol_domain)} variants located in a protein domain that is seemingly intolerant to AA changes according to AM scores")

    missense = df['Consequence'].str.contains('missense_variant')
    logger.info(f"There are {missense.sum()} missense variants")
    missense_damaging = df["am_class"].fillna("").str.contains('athogenic') | (df['PrimateAI'].fillna(0) > 0.9)
    logger.info(f"There are {missense_damaging.sum()} missense variants that are considered damaging by AlphaMissense and PrimateAI")

    # AM-benign suppression of PM1.
    #
    # Two independent ways for AlphaMissense to call a substitution benign:
    #   1. the categorical label (am_class contains 'benign'), and
    #   2. the numeric score below the same 0.564 cutoff BP4 uses for
    #      missense_benign in acmg_pp3_bp4_bp7_insilico.py.
    # The label alone misses 'ambiguous' rows whose score is clearly benign, and
    # misses rows where am_class is absent but am_pathogenicity is populated.
    #
    # PrimateAI acts only as a veto when it positively disagrees (>= 0.8, matching
    # BP4's threshold). A MISSING PrimateAI score must not veto an AM-benign call:
    # the previous `fillna(1.0) < 0.5` turned absent PrimateAI evidence into an
    # assertion of damage, so any AM-benign variant lacking PrimateAI still fired PM1.
    am_pathogenicity_num = pd.to_numeric(df["am_pathogenicity"], errors='coerce')
    primateai_num = pd.to_numeric(df["PrimateAI"], errors='coerce')
    am_says_benign = (
        df["am_class"].fillna("").str.contains('benign')
        | (am_pathogenicity_num < 0.564)
    )
    primateai_not_damaging = ~(primateai_num >= 0.8).fillna(False)
    missense_benign = am_says_benign & primateai_not_damaging
    logger.info(f"There are {missense_benign.sum()} missense variants that are considered benign by AlphaMissense and PrimateAI")
    if clinvar_patho is not None:
        missense_damaging = missense_damaging | clinvar_patho
        logger.info(f"There are {missense_damaging.sum()} missense variants that are considered damaging after considering ClinVar high-confidence pathogenic status")
    missense_damaging = missense_damaging & missense

    # Load the 3-arm PM1 regions pickle
    pm1_regions = pickle.load(gzip.open(pm1_regions_pkl, 'rb')) if pm1_regions_pkl.endswith(".gz") else pickle.load(open(pm1_regions_pkl, 'rb'))
    hcseeker_hotspots = pm1_regions.get("hcseeker_hotspots", {})
    rmc_constrained = pm1_regions.get("rmc_constrained", {})
    meta = pm1_regions.get("metadata", {})
    logger.info(f"PM1 regions loaded: {meta.get('n_das_intolerant_domains', '?')} DAS domains, "
                f"{meta.get('n_hcseeker_hotspot_residues', '?')} HCSeeker residues, "
                f"{meta.get('n_rmc_constrained_residues', '?')} RMC residues")

    # Check each variant against HCSeeker hotspots and RMC constrained regions
    args = [(row, hcseeker_hotspots, rmc_constrained) for row in row_dicts]
    with mp.Pool(threads) as pool:
        results = pool.starmap(locate_pm1_region, args)
        hcseeker_hit = np.array([r[0] for r in results])
        rmc_hit = np.array([r[1] for r in results])

    logger.info(f"There are {np.sum(hcseeker_hit)} variants located in HCSeeker mutational hotspots")
    logger.info(f"There are {np.sum(rmc_hit)} variants located in gnomAD RMC-constrained regions")


    consq = df["Consequence"].fillna("").astype(str)
    nmd = df["NMD"].fillna(".").astype(str)
    lof_filter = df.get("LoF_filter", pd.Series(".", index=df.index)).fillna(".").astype(str)

    truncating = (
        consq.str.contains("stop_gained", regex=False)
        | consq.str.contains("frameshift", regex=False)
    )

    non_nmd_truncating = truncating & (
        nmd.str.contains("escaping", regex=False)
        | lof_filter.str.contains("END_TRUNC", regex=False)
    )

    inframe_indels = (
        consq.str.contains("inframe_deletion", regex=False)
        | consq.str.contains("inframe_insertion", regex=False)
    )

    indels = inframe_indels | non_nmd_truncating

    pvs1_double_count = pvs1_criteria >= 2
    pm1_criteria = ( hcseeker_hit & ((missense & ~missense_benign) | indels) ) | \
                   ( rmc_hit & ((missense & ~missense_benign) | indels) ) | \
                   ( loc_intol_domain & ((missense & ~missense_benign) | indels) )
    pm1_criteria = pm1_criteria & np.logical_not(pvs1_double_count)
    pm1_array = np.zeros(len(df), dtype=int)
    pm1_array[pm1_criteria] = 2
    return pm1_array, loc_intol_domain


def analyze_bp1_pp2(gene_stat_dict: Dict):
    """
    Analyze BP1 and PP2 criteria for a gene based on variant statistics.

    Args:
        gene_stat_dict: Dictionary containing gene-level variant statistics

    Returns:
        Dictionary with BP1 and PP2 analysis results
    """
    gene = gene_stat_dict['ensembl_id']
    if len(gene_stat_dict) <= 1:
        return (gene, (False, False))

    # Get pathogenic variants
    patho_variants = gene_stat_dict['clinvar_pathogenic']
    total_patho = len(patho_variants)

    if total_patho == 0:
        return (gene, (False, False))

    # BP1 Analysis
    # Combine large AA change and NMD variants
    large_aa_nmd_variants = (gene_stat_dict['large_aachange'] | gene_stat_dict['putative_nmd_variants'])
    patho_large_aa_nmd = len(patho_variants & large_aa_nmd_variants)
    bp1_fraction = patho_large_aa_nmd / total_patho
    bp1_granted = bp1_fraction >= 0.7

    # PP2 Analysis
    # Create contingency table for small vs large AA changes
    small_aa_variants = gene_stat_dict['small_aachange']
    large_aa_variants = gene_stat_dict['large_aachange']

    # Count variants in each category
    patho_small_aa = len(patho_variants & small_aa_variants)
    patho_large_aa = len(patho_variants & large_aa_variants)
    non_patho_small_aa = len(small_aa_variants - patho_variants)
    non_patho_large_aa = len(large_aa_variants - patho_variants)

    # Create contingency table
    table = [[patho_small_aa, patho_large_aa],
             [non_patho_small_aa, non_patho_large_aa]]

    # Initialize PP2 results
    pp2_granted = False
    pp2_test_used = 'none'
    pp2_pvalue = np.nan
    pp2_odds_ratio = np.nan

    # Check if we can perform statistical tests
    total_count = sum(sum(row) for row in table)
    min_cell_count = min(min(row) for row in table)

    if (total_count >= 5 and
        all(sum(row) > 0 for row in table) and
        all(sum(col) > 0 for col in zip(*table))):

        if total_count >= 20 and min_cell_count >= 5:
            # Use Chi-square test
            chi2, pp2_pvalue = stats.chi2_contingency(table)[0:2]
            pp2_odds_ratio = ((table[0][1] * table[1][0]) /
                            (table[0][0] * table[1][1]))  # Note: reversed for large AA changes
            pp2_test_used = 'chi_square'
        else:
            # Use Fisher's exact test
            pp2_odds_ratio, pp2_pvalue = stats.fisher_exact(table, alternative='greater')
            pp2_test_used = 'fisher'

        # Adjust p-value of 1 to 0.99999 to avoid log10(1) = 0
        pp2_pvalue = min(pp2_pvalue, 0.99999)

        # Grant PP2 if there is NO significant enrichment of large AA changes
        # (p-value > 0.05) and missense variants are ≥50% of pathogenic variants
        pp2_fraction = patho_small_aa / total_patho
        pp2_granted = (pp2_pvalue > 0.05) and (pp2_fraction >= 0.4)
    else:
        # Fall back to fraction-based approach
        pp2_fraction = patho_small_aa / total_patho
        pp2_granted = pp2_fraction >= 0.5
        pp2_test_used = 'fraction'

    result_dict = {
        'bp1_granted': bp1_granted,
        'bp1_fraction': bp1_fraction,
        'pp2_granted': pp2_granted,
        'pp2_fraction': patho_small_aa / total_patho if total_patho > 0 else 0.0,
        'pp2_test_used': pp2_test_used,
        'pp2_pvalue': pp2_pvalue,
        'pp2_odds_ratio': pp2_odds_ratio
    }

    logger.debug(f"The result_dict looks like {result_dict}")

    return (gene, (result_dict['pp2_granted'], result_dict['bp1_granted']))


def PP2_BP1_criteria(df: pd.DataFrame,
                     clinvar_stats_dict: dict,
                     am_intol_domains_tsv: str,
                     threads: int = 10) -> Tuple[pd.Series, pd.Series]:
    """
    Efficient implementation of PP2/BP1 criteria evaluation using
    ClinVar statistics dictionary.

    Args:
        df: Variant annotation DataFrame
        clinvar_stats_dict: Dictionary containing ClinVar statistics per gene
        threads: Number of threads for multiprocessing

    Returns:
        Tuple[pd.Series, pd.Series]: PP2 and BP1 criteria boolean series
    """

    # Get unique genes
    unique_genes = df['Gene'].dropna().unique()
    logger.info(f"Processing BP1/PP2 criteria for {len(unique_genes)} unique genes")

    # Process genes in parallel
    with mp.Pool(threads) as pool:
        gene_results = pool.imap_unordered(analyze_bp1_pp2, (clinvar_stats_dict.get(gene, {"ensembl_id": gene}) for gene in unique_genes))
        # Convert to dictionary for faster lookups
        gene_criteria_dict = dict(gene_results)

    # Create result Series using map operation (vectorized)
    gene_pp2 = df['Gene'].map(lambda g: gene_criteria_dict.get(g, (False, False))[0] if pd.notna(g) else False)
    gene_bp1 = df['Gene'].map(lambda g: gene_criteria_dict.get(g, (False, False))[1] if pd.notna(g) else False)

    # Apply PP2 only to missense variants
    is_missense = df['Consequence'].str.contains('missense_variant', na=False)
    pp2_criteria = gene_pp2 & is_missense

    # BP1 applies to all variants
    bp1_criteria = gene_bp1 & is_missense

    # Log stats
    logger.info(f"Granted PP2 to {pp2_criteria.sum()} variants based on ClinVar stats")
    logger.info(f"Granted BP1 to {bp1_criteria.sum()} variants based on ClinVar stats")

    # Load the AlphaMissense intolerant domains
    am_intol_domains_df = pd.read_table(am_intol_domains_tsv, low_memory=False)
    am_intol_domains_df["Ensembl_Gene_ID"] = am_intol_domains_df["domain"].str.split(":").str[0]
    by_gene = am_intol_domains_df.groupby("Ensembl_Gene_ID", as_index=False)
    gene_all_intol_domains = {}
    for gene, gene_df in by_gene:
        if gene_df["is_more_tolerant"].any():
            gene_all_intol_domains[gene] = False
        else:
            gene_all_intol_domains[gene] = True

    # Apply PP2 to all variants
    pp2_am = df['Gene'].map(lambda g: gene_all_intol_domains.get(g, False))
    pp2_gene_am = df["Gene_avg_AM_score"].fillna(0) > 0.564
    pp2_criteria = pp2_criteria | pp2_am | pp2_gene_am
    pp2_criteria = pp2_criteria & is_missense
    logger.info(f"Granted PP2 to {pp2_criteria.sum()} variants based on updates from AlphaMissense intolerant domains")

    pp2_array = np.zeros(len(df), dtype=int)
    pp2_array[pp2_criteria] = 1
    bp1_array = np.zeros(len(df), dtype=int)
    bp1_array[bp1_criteria] = 1
    return pp2_array, bp1_array
