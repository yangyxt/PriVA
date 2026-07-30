#!/usr/bin/env python3
"""Combine the criteria into a classification, a probability, and a ranking.

The end of the chain. Every criterion has been assigned by now; this module
decides what they add up to and which variants a clinician sees first.

    summarize_acmg_criteria         collapse the per-criterion arrays into one
                                    readable summary column and a boolean
                                    matrix
    calculate_posterior_probability convert the combined evidence into a
                                    posterior probability of pathogenicity,
                                    on the Tavtigian Bayesian framework
    compute_control_common          flag alleles common in the local control
                                    cohort, which are ranked down regardless
                                    of what the criteria say
    sort_and_rank_variants          order candidates, applying the expression
                                    penalty for poorly expressed transcripts
                                    and the relevant- and dispensable-gene
                                    lists

sort_and_rank_variants is also called on its own by
/paedyl01/disk1/yangyxt/PriVA/scripts/resort_step3_only.py to re-rank an
existing result table without recomputing the criteria.
"""

import logging
import os
import pandas as pd
import numpy as np
from typing import Tuple, Dict
from scipy.stats import binom, beta as beta_dist

from acmg_consequence import coerce_pext_value
from acmg_variant_mechanism import _variant_mechanism_masks


self_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(os.path.dirname(self_dir), "data")
DEFAULT_DISPENSABLE_GENE_LIST = os.path.join(
    data_dir, "dispensable_genes", "dispensable_gene_list.txt"
)

logger = logging.getLogger(__name__)


def create_criteria_summary(row, criteria_order=None, clingen_evidence=None):
    """
    Create a summary string of active criteria for a given row.

    Args:
        row: A row from the criteria matrix
        criteria_order: List of criteria in order

    Returns:
        A semicolon-separated string of active criteria
    """

    strength_level_suffix = {1: "Supporting", 2: "Moderate", 3: "Strong", 4: "Very Strong", 5: "Standalone"}
    label_strength_level = {"P": "Supporting", "M": "Moderate", "S": "Strong", "V": "Very Strong", "A": "Standalone"}

    active_criteria = [ col if strength_level_suffix.get(int(row[col]), None) == label_strength_level.get(col[1], None) else col + "_" + strength_level_suffix.get(int(row[col]), None) for col in criteria_order if int(row[col]) > 0 ]

    # TEMPORARILY COMMENTED OUT: ClinGen evidence inheritance
    # ClinGen evidence lookup is deliberately not applied here. clingen_evidence
    # is accepted so the parameter stays wired through the pipeline, but a
    # ClinGen assertion does not override the criteria PriVA derived itself.

    return ";".join(active_criteria) if active_criteria else np.nan


def summarize_acmg_criteria(df: pd.DataFrame, criteria_dict: Dict[str, np.ndarray], clingen_map_pkl: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create a summary column of assigned criteria and a boolean matrix of all criteria.

    Args:
        df: Input DataFrame
        criteria_dict: Dictionary mapping criteria names to boolean arrays

    Returns:
        Tuple of:
        - Original DataFrame with added criteria summary column
        - Boolean DataFrame of all criteria (rows=variants, columns=criteria)
    """
    # Define criteria order
    criteria_order = [
        # Very Strong
        'PVS1',
        # Strong
        'PS1', 'PS2', 'PS3', 'PS4',
        # Moderate
        'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6',
        # Supporting
        'PP1', 'PP2', 'PP3', 'PP4', 'PP5',
        # Stand-alone
        'BA1',
        # Strong Benign
        'BS1', 'BS2', 'BS3', 'BS4',
        # Supporting Benign
        'BP1', 'BP2', 'BP3', 'BP4', 'BP5', 'BP6', 'BP7'
    ]

    # Create criteria matrix
    criteria_matrix = pd.DataFrame(
        {name: criteria_dict.get(name, np.zeros(len(df), dtype=int))
         for name in criteria_order},
        index=df.index
    )

    criteria_matrix["variant_id"] = df["variant_id"]
    logger.info(f"Criteria matrix: \n{criteria_matrix[:20].to_string(index=False)}")

    # Check which column has no integer values and return error
    for col in criteria_order:
        if not criteria_matrix[col].isin([0, 1, 2, 3, 4, 5]).all():
            raise ValueError(f"Column {col} has non-integer values")

    # TEMPORARILY COMMENTED OUT: ClinGen evidence inheritance
    # Uncomment the following lines to re-enable ClinGen evidence lookup
    # Load ClinGen evidence map
    # clingen_map = pickle.load(gzip.open(clingen_map_pkl, "rb")) if clingen_map_pkl.endswith(".gz") else pickle.load(open(clingen_map_pkl, "rb"))
    clingen_map = {}  # Empty map - PriVA will use its own criteria assignments

    # Add summary column to original DataFrame
    df['ACMG_criteria'] = criteria_matrix.apply(create_criteria_summary, axis=1, criteria_order=criteria_order, clingen_evidence=clingen_map)

    # The criteria matrix is deliberately not overridden with ClinGen evidence.
    # clingen_map_pkl is accepted so the parameter stays wired through, but
    # PriVA reports the criteria it derived rather than a curated verdict.

    return df, criteria_matrix


def calculate_posterior_probability(row, prior_probability=0.1, exp_base=2, odds_pvst=350):
    """
    Calculates the posterior probability of pathogenicity for a variant based on its criteria assignment using a Bayesian framework.

    Args:
        criteria_assignment (dict): A dictionary representing the strength of each evidence criterion.
                                     Keys are criteria names (e.g., "PVS1", "PS1", "PM2", "PP3", "BS1", "BP4"),
                                     and values are the corresponding strength scores.
        prior_probability (float): The prior probability of pathogenicity (default: 0.1).
        evidence_weights (dict): A dictionary mapping criteria names to their numerical weights,
                                 according to the latest ClinGen SVI recommendations (default: illustrative values).
                                 PM2 is downgraded to a PP criterion

    Returns:
        float: The posterior probability of pathogenicity.
    """
    if isinstance(row, pd.Series):
        # Convert the Series to a dictionary
        criteria_assignment = row.to_dict()

    # evidence_weights = {"PVS1": 1,
    #                     "PS1": 1/exp_base, "PS2": 1/exp_base, "PS3": 1/exp_base, "PS4": 1/exp_base,
    #                     "PM1": 1/exp_base**2, "PM2": 1/exp_base**3, "PM3": 1/exp_base**2, "PM4": 1/exp_base**2, "PM5": 1/exp_base**2, "PM6": 1/exp_base**2,
    #                     "PP1": 1/exp_base**3, "PP2": 1/exp_base**3, "PP3": 1/exp_base**3, "PP4": 1/exp_base**3, "PP5": 1/exp_base**3,
    #                     "BA1": -5,
    #                     "BS1": -1/exp_base, "BS2": -1/exp_base, "BS3": -1/exp_base, "BS4": -1/exp_base,
    #                     "BP1": -1/exp_base**3, "BP2": -1/exp_base**3, "BP3": -1/exp_base**3, "BP4": -1/exp_base**3, "BP5": -1/exp_base**3, "BP6": -1/exp_base**3, "BP7": -1/exp_base**3}

    evidence_weights = {4: 1.0, 3: 1/exp_base, 2: 1/exp_base**2, 1: 1/exp_base**3, 5: 2.0}

    # Default odds function (Illustrative - replace with a calibrated function based on SVI guidelines)
    odds_function = lambda total_weight: np.power(odds_pvst, total_weight)

    # There are 2 internal inconsistencies in the SVI guidelines.
    # Pathogenic(ii) and Likely Pathogenic(i) do not generate a posterior probability of the corresponding class.
    # We need to handle this by setting the odds of pathogenic to 0 if the total weight is negative. If we run into 2 PS evidences, we adjust the total_weight_pathogenic to 1514
    # If there is a PVS + 1 PM, we need to adjust their sum to 350

    patho_matches = [n for c,n in criteria_assignment.items() if c.startswith("P") and int(n) > 0]
    ps_matches = [n for n in patho_matches if n == 3]
    pm_matches = [n for n in patho_matches if n == 2]
    pp_matches = [n for n in patho_matches if n == 1]
    pvs_matches = [n for n in patho_matches if n == 4]
    benign_matches = [n for c,n in criteria_assignment.items() if c.startswith("B") and int(n) > 0]
    bs_matches = [n for n in benign_matches if n == 3]
    bp_matches = [n for n in benign_matches if n == 1]
    ba_matches = [n for n in benign_matches if n == 5]
    bm_matches = [n for n in benign_matches if n == 2]

    odds_benign = odds_function(-float(sum([evidence_weights[n] for n in bs_matches + bp_matches + ba_matches + bm_matches])))
    if len(ps_matches) >= 2:
        ps_odds = 1514 * odds_function(float(sum([evidence_weights[int(n)] for n in ps_matches[:-2]])))
    else:
        ps_odds = odds_function(float(sum([evidence_weights[int(n)] for n in ps_matches])))

    if len(pvs_matches) > 0 and len(pm_matches) > 0 and len(pp_matches) == 0:
        pv_pm_odds = odds_pvst * odds_function(float(sum([evidence_weights[int(n)] for n in pm_matches[:-1]])))
    else:
        pv_pm_odds = odds_function(float(sum([evidence_weights[int(n)] for n in pvs_matches + pm_matches])))

    pp_odds = odds_function(float(sum([evidence_weights[int(n)] for n in pp_matches])))

    odds_path = ps_odds * pp_odds * pv_pm_odds

    # Calculate odds of pathogenicity
    odds_path = odds_path * odds_benign

    # Calculate posterior probability
    posterior_probability = (odds_path * prior_probability) / ((odds_path - 1) * prior_probability + 1)
    # Round posterior probability to 3 decimal places
    posterior_probability = round(posterior_probability, 4)

    if len(ba_matches) > 0:
        posterior_probability = 0

    # Classify the variant based on the posterior probability
    if posterior_probability > 0.99:
        acmg_class = "Pathogenic"
    elif posterior_probability > 0.899:
        acmg_class = "Likely Pathogenic"
    elif posterior_probability < 0.0005:
        acmg_class = "Benign"
    elif posterior_probability < 0.05:
        acmg_class = "Likely Benign"
    else:
        acmg_class = "Uncertain Significance"

    return posterior_probability, acmg_class


def compute_control_common(df: pd.DataFrame, proband: str,
                           cutoff: float = 0.01, conf_level: float = 0.999) -> pd.DataFrame:
    """Flag variants common in the control cohort via binomial CDF test.

    Uses control_AC/AN/AF/nhomalt fields. For homozygous proband calls,
    tests nhomalt count; for heterozygous, tests AC count. Falls back to
    AF > cutoff when counts are unavailable.

    Flagged variants receive control_common_index = 0.5 (rank penalty).
    """
    control_cols = {"control_AC", "control_AN", "control_AF", "control_nhomalt"}
    if not control_cols.intersection(set(df.columns)):
        df["control_common"] = False
        df["control_common_index"] = 1.0
        return df

    an = pd.to_numeric(df.get("control_AN"), errors='coerce')
    ac = pd.to_numeric(df.get("control_AC"), errors='coerce')
    af = pd.to_numeric(df.get("control_AF"), errors='coerce')
    nhomo = pd.to_numeric(df.get("control_nhomalt"), errors='coerce')

    valid_an = an.fillna(0) > 0

    # Determine proband homozygosity from GT column
    if proband and proband in df.columns:
        gt = df[proband].astype(str).str.split(":").str[0]
        is_hom = gt.isin(["1/1", "1|1", "1"])
    else:
        is_hom = pd.Series(False, index=df.index)

    # Route 1: hom proband + nhomo available → binomial test on nhomo vs AN/2
    use_nhomo = valid_an & is_hom & nhomo.notna()
    # Route 2: AC available (and not using nhomo) → binomial test on AC vs AN
    use_ac = valid_an & ~use_nhomo & ac.notna()
    # Route 3: only AF available → simple threshold
    use_af = valid_an & ~use_nhomo & ~use_ac & af.notna()

    common = pd.Series(False, index=df.index)

    if use_nhomo.any():
        m = use_nhomo
        common[m] = binom.cdf(
            nhomo[m].clip(lower=0).astype(int).values,
            (an[m] / 2).round().clip(lower=1).astype(int).values,
            cutoff
        ) > conf_level

    if use_ac.any():
        m = use_ac
        common[m] = binom.cdf(
            ac[m].clip(lower=0).astype(int).values,
            an[m].round().clip(lower=1).astype(int).values,
            cutoff
        ) > conf_level

    if use_af.any():
        m = use_af
        common[m] = af[m] > cutoff

    df["control_common"] = common

    # Graded penalty via Beta posterior lower bound: penalize more when
    # both AF and sample size are large (high statistical confidence).
    # Floor at 0.01, ceiling at 0.5 (minimum penalty for any common variant).
    ac_safe = ac.fillna(0).values.astype(float)
    an_safe = an.fillna(0).values.astype(float)
    confident_af = np.zeros(len(df))
    valid_for_beta = common.values & (an_safe > 0)
    if valid_for_beta.any():
        confident_af[valid_for_beta] = beta_dist.ppf(
            0.01, ac_safe[valid_for_beta] + 1,
            an_safe[valid_for_beta] - ac_safe[valid_for_beta] + 1
        )
    df["control_common_index"] = np.where(
        common, np.clip(1.0 - confident_af, 0.01, 0.5), 1.0
    )

    n_flagged = common.sum()
    if n_flagged > 0:
        logger.info(f"Control-common penalty: {n_flagged} variants flagged as common in control cohort (graded Beta penalty, ceiling 0.5)")

    return df


def sort_and_rank_variants(df: pd.DataFrame,
                           ped_df: pd.DataFrame,
                           fam_name: str,
                           pext_tissues: str = "",
                           pext_low_expression_cutoff: float = 0.1,
                           pext_penalty_floor: float = 0.8,
                           pext_penalty_shape: float = 0.5,
                           relevant_gene_list: str = None,
                           dispensable_gene_list: str = os.path.join(data_dir, "dispensable_genes", "dispensable_gene_list.txt")) -> pd.DataFrame:
    """
    Sort variants by their maximum ACMG quantitative score and add ranking.

    Args:
        df: DataFrame with ACMG_quant_score column and variant information

    Returns:
        DataFrame sorted by variant's max ACMG score with added variant rank column
    """
    df = df.copy()
    for stale_col in ("haplo_insuf_index", "zygosity_inheritance_mechanism_index"):
        if stale_col in df.columns:
            df = df.drop(columns=[stale_col])

    df["sort_index"] = df.loc[:, "ACMG_quant_score"]
    df = df.loc[df["BIOTYPE"] == "protein_coding", :]

    # NOTE: proband_het must be created AFTER filtering to match df's shape
    if ped_df is not None and fam_name is not None:
        proband = ped_df.loc[(ped_df['#FamilyID'] == fam_name) & (ped_df['Phenotype'].isin(["2", 2])), 'IndividualID'].values[0]
        proband_het = (df.loc[:, proband].str.count("1") == 1)
    else:
        logger.warning("No ped_table provided, skip the zygosity/inheritance/mechanism compatibility penalty")
        proband = None
        proband_het = pd.Series([False] * len(df), index=df.index)  # Use Series with matching index

    # Control-common penalty: halve sort_index for variants common in control cohort
    df = compute_control_common(df, proband)
    df.loc[:, "sort_index"] = df.loc[:, "sort_index"] * df.loc[:, "control_common_index"]

    # ---------------------------------------------------------------------
    # Zygosity / inheritance / mechanism compatibility penalty.
    #
    # Down-weight a single-allele predicted-LOF call when the upstream hub does
    # not provide a dominant-compatible assertion. No raw annotation or
    # gene-level fallback is allowed at this stage.
    # ---------------------------------------------------------------------
    mechanism_masks = _variant_mechanism_masks(df)
    has_dom_lof_compatible = mechanism_masks["has_dominant_compatible"]
    exact_gof = mechanism_masks["is_exact_gof"]
    asserted_lof_effect = mechanism_masks["is_predicted_lof"]
    dn_without_hi_or_ambiguous = (
        mechanism_masks["has_dom_dn_history"]
        & ~mechanism_masks["has_dom_lof_history"]
        & ~mechanism_masks["has_dom_unresolved_history"]
    )
    heterozygous_lof_incompatible = (
        asserted_lof_effect
        & ~exact_gof
        & ~has_dom_lof_compatible
    )

    high_acmg_lof = df["ACMG_quant_score"] > 0.89
    compatibility_penalty_rows = high_acmg_lof & proband_het & heterozygous_lof_incompatible
    pext_modulation_rows = high_acmg_lof

    logger.info(
        "Zygosity/inheritance/mechanism compatibility signals for ranking: "
        f"hub_profile_rows={mechanism_masks['modern_profile'].sum()}, "
        f"asserted_lof_effect={pd.Series(asserted_lof_effect, index=df.index).sum()}, "
        f"dominant_DN_without_HI_or_ambiguous={dn_without_hi_or_ambiguous.sum()}, "
        f"heterozygous_lof_incompatible={heterozygous_lof_incompatible.sum()}, "
        f"compatibility_penalty_rows={compatibility_penalty_rows.sum()}"
    )
    df["zygosity_inheritance_mechanism_compatibility"] = 1.0
    df.loc[compatibility_penalty_rows, "zygosity_inheritance_mechanism_compatibility"] = 0.9

    if dispensable_gene_list is not None:
        dispensable_genes = pd.read_table(dispensable_gene_list, header=None, names=["SYMBOL"], comment="#") # Skip the header row starting with #
        dispensable_genes = set(dispensable_genes["SYMBOL"].tolist())
        df["dispensable_gene_index"] = 1.0
        df.loc[df["SYMBOL"].isin(dispensable_genes), "dispensable_gene_index"] = 0.9
        df.loc[:, "sort_index"] = df.loc[:, "sort_index"] * df.loc[:, "dispensable_gene_index"]

    if relevant_gene_list is not None:
        relevant_genes = pd.read_table(relevant_gene_list, header=None, names=["SYMBOL"], comment="#") # Skip the header row starting with #
        relevant_genes = set(relevant_genes["SYMBOL"].tolist())
        df["relevant_gene_index"] = 1.0
        df.loc[df["SYMBOL"].isin(relevant_genes), "relevant_gene_index"] = 1.2  # Prioritize the variants in relevant genes
        df.loc[:, "sort_index"] = df.loc[:, "sort_index"] * df.loc[:, "relevant_gene_index"]

    df.loc[compatibility_penalty_rows, "sort_index"] = (
        df.loc[compatibility_penalty_rows, "sort_index"]
        * df.loc[compatibility_penalty_rows, "zygosity_inheritance_mechanism_compatibility"]
    )

    # Apply pext-based expression modulation if pext columns are available
    # Low pext = variant in low-expression region = deprioritize for LoF interpretation
    # Use the first tissue-specific pext column if provided, otherwise use PEXT_MEAN
    pext_col_to_use = None
    if pext_tissues:
        pext_tissue_list = [t.strip() for t in pext_tissues.split(",") if t.strip()]
        if pext_tissue_list:
            # Try the first tissue-specific column
            first_tissue_col = f"PEXT_{pext_tissue_list[0].replace('-', '_')}"
            if first_tissue_col in df.columns:
                pext_col_to_use = first_tissue_col
                logger.info(f"Using tissue-specific pext column for sorting: {first_tissue_col}")

    # Fallback to PEXT_MEAN if no tissue-specific column found
    if pext_col_to_use is None and "PEXT_MEAN" in df.columns:
        pext_col_to_use = "PEXT_MEAN"
        logger.info("Using PEXT_MEAN column for expression-aware sorting")

    if pext_col_to_use is not None:
        if pext_low_expression_cutoff <= 0:
            raise ValueError("pext_low_expression_cutoff must be > 0")
        if not 0 <= pext_penalty_floor <= 1:
            raise ValueError("pext_penalty_floor must be between 0 and 1")
        if pext_penalty_shape <= 0:
            raise ValueError("pext_penalty_shape must be > 0")

        # Cummings 2020 used PEXT < 0.1 as the low-expression bin.
        # Bound the ranking penalty so PEXT cannot overwhelm ACMG evidence.
        pext_numeric = df[pext_col_to_use].map(coerce_pext_value).fillna(1.0).clip(0, 1)
        scaled_pext = (pext_numeric / pext_low_expression_cutoff).clip(0, 1)
        df["pext_sort_index"] = np.where(
            pext_numeric >= pext_low_expression_cutoff,
            1.0,
            pext_penalty_floor + (1.0 - pext_penalty_floor) * np.power(scaled_pext, pext_penalty_shape)
        )
        # Historical PriVA behavior applied PEXT modulation to high ACMG-score rows.
        df.loc[pext_modulation_rows, "sort_index"] = (
            df.loc[pext_modulation_rows, "sort_index"]
            * df.loc[pext_modulation_rows, "pext_sort_index"]
        )
        logger.info(
            "Applied pext-based sort modulation using %s "
            "(cutoff=%s, floor=%s, shape=%s)",
            pext_col_to_use,
            pext_low_expression_cutoff,
            pext_penalty_floor,
            pext_penalty_shape,
        )

    # Group by variant coordinates to get max score per variant
    variant_groups = df.groupby(['chrom', 'pos', 'ref', 'alt'])
    max_scores = variant_groups['sort_index'].transform('max')

    # Sort the DataFrame by max scores (descending) while maintaining variant grouping
    df['max_variant_score'] = max_scores
    df_sorted = df.sort_values(
        by=['max_variant_score', 'ACMG_quant_score', 'chrom', 'pos', 'ref', 'alt'],
        ascending=[False, False, True, True, True, True]
    )

    # Add variant rank based on actual sorted order (not lexicographic ngroup)
    # Use factorize() which assigns codes based on order of first appearance (0-indexed)
    # Adding 1 gives us 1-indexed ranks where rank 1 = first (highest-scored) variant
    if 'variant_id' in df_sorted.columns:
        # factorize returns (codes, uniques) where codes are assigned in order of first appearance
        codes, _ = pd.factorize(df_sorted['variant_id'], sort=False)
        df_sorted['variant_rank'] = codes + 1  # Convert 0-indexed to 1-indexed
    else:
        # Create variant key for ranking
        var_key = (df_sorted['chrom'].astype(str) + ':' +
                   df_sorted['pos'].astype(str) + ':' +
                   df_sorted['ref'].astype(str) + ':' +
                   df_sorted['alt'].astype(str))
        codes, _ = pd.factorize(var_key, sort=False)
        df_sorted['variant_rank'] = codes + 1  # Convert 0-indexed to 1-indexed

    # Clean up temporary column
    df_sorted = df_sorted.drop('max_variant_score', axis=1)

    return df_sorted
