#!/usr/bin/env python3
"""PM2, BS1 and BA1 -- how common the allele is in the population.

One axis, three thresholds, read from gnomAD.

    PM2  absent, or extremely rare, in population databases
    BS1  more frequent than the disorder's expected incidence allows
    BA1  above 5% in any reasonably sized population -- stand-alone benign

control_false_neg_rate is shared by BS1 and BA1. It asks whether the observed
allele number is large enough for an observed frequency to be trusted: at low
coverage a variant can look rare simply because few alleles were called, and
calling that benign-supporting would be an error. BS1 is mechanism-aware
through the variant-level assertions; BA1 is not, because at 5% no disease
mechanism makes a pathogenic call defensible.
"""

import logging
import pickle
import gzip
import pandas as pd
import numpy as np
from typing import Tuple
from scipy.stats import binom

from acmg_variant_mechanism import _variant_condition_masks


logger = logging.getLogger(__name__)


def control_false_neg_rate(
    allele_frequencies: pd.Series,
    allele_numbers: pd.Series,
    af_threshold = 0.001,
    alpha: float = 0.01
) -> Tuple[pd.Series, pd.Series]:
    """
    Performs a one-sided binomial test for each variant with potentially
    different allele frequency thresholds for each variant.

    Args:
        allele_frequencies (pd.Series): Series of observed allele frequencies.
        allele_numbers (pd.Series): Series of total observed allele numbers.
        af_threshold (float, pd.Series, or np.ndarray): The allele frequency
            threshold(s) for the null hypothesis (H₀: p_true < af_threshold).
            If array-like, must match the length of allele_frequencies.
        alpha (float): The significance level for rejecting H₀. Default is 0.01.

    Returns:
        Tuple[pd.Series, pd.Series]: p_values and reject_h0 decisions.
    """
    if not isinstance(allele_frequencies, pd.Series):
        allele_frequencies = pd.Series(allele_frequencies)
    if not isinstance(allele_numbers, pd.Series):
        allele_numbers = pd.Series(allele_numbers)

    if not len(allele_frequencies) == len(allele_numbers):
        raise ValueError("Input Series must have the same length.")

    # Handle af_threshold input
    if isinstance(af_threshold, (float, int)):
        # If scalar, convert to Series with same index as allele_frequencies
        af_thresh_series = pd.Series(af_threshold, index=allele_frequencies.index)
    elif isinstance(af_threshold, np.ndarray):
        # Convert numpy array to Series, aligning index
        if len(af_threshold) != len(allele_frequencies):
            raise ValueError("If af_threshold is an array, it must have the same length as allele_frequencies.")
        af_thresh_series = pd.Series(af_threshold, index=allele_frequencies.index)
    elif isinstance(af_threshold, pd.Series):
        af_thresh_series = af_threshold
        if len(af_thresh_series) != len(allele_frequencies):
            raise ValueError("If af_threshold is a Series, it must have the same length as allele_frequencies.")
    else:
        raise TypeError("af_threshold must be a float, int, pandas Series, or numpy array.")

    # Ensure threshold is numeric and within [0,1]
    af_thresh_series = pd.to_numeric(af_thresh_series, errors='coerce').clip(0, 1)
    if af_thresh_series.isnull().any():
        raise ValueError("Non-numeric values found in af_threshold Series/array after coercion.")

    # --- Input Cleaning and Preparation (same as before) ---
    an = pd.to_numeric(allele_numbers, errors='coerce')
    an = an.fillna(0).astype(int)
    an[an < 0] = 0

    af = pd.to_numeric(allele_frequencies, errors='coerce')
    af = af.clip(0, 1)

    ac_float = np.where(np.isnan(af), 0, af * an)
    ac_observed = np.where(np.isnan(ac_float), 0, np.round(ac_float).astype(int))
    ac_observed = np.minimum(ac_observed, an)
    ac_observed[an == 0] = 0

    valid_mask = (an > 0) & (~af.isna()) & (~allele_numbers.isna())

    # --- P-value Calculation ---
    p_values = pd.Series(np.nan, index=allele_frequencies.index)

    # Handle edge case k=0: P(AC >= 0) is always 1
    p_values.loc[valid_mask & (ac_observed == 0)] = 1.0

    # Calculate for k > 0
    mask_k_pos = valid_mask & (ac_observed > 0)
    if mask_k_pos.any():
        # Get the threshold values for each valid variant
        thresholds_to_use = af_thresh_series[mask_k_pos].values

        p_values.loc[mask_k_pos] = binom.sf(
            k=ac_observed[mask_k_pos] - 1,
            n=an[mask_k_pos],
            p=thresholds_to_use  # Now using the array of thresholds
        )

    # --- Decision ---
    reject_h0 = p_values <= alpha
    logger.info(f"There are {reject_h0.sum()} variants having their p-value less than {alpha}, {reject_h0.isna().sum()} datapoints are NA")

    return p_values, reject_h0


def gnomAD_rare_AF(df: pd.DataFrame, cutoff: float) -> np.ndarray:
    if isinstance(cutoff, float):
        return np.where(df['gnomAD_joint_AN_max'].fillna(1000000) >= 1/cutoff, \
                        df['gnomAD_joint_AF_max'].fillna(0) <= cutoff, \
                        df['gnomAD_joint_AF'].fillna(0) <= cutoff)
    elif isinstance(cutoff, str):
        return np.where(df['gnomAD_joint_AN_max'].fillna(1000000) >= 1/(df[cutoff].fillna(10e-6)), \
                        df['gnomAD_joint_AF_max'].fillna(0) <= df[cutoff].fillna(0), \
                        df['gnomAD_joint_AF'].fillna(0) <= df[cutoff].fillna(0))


def PM2_criteria(df: pd.DataFrame,
                 clinvar_patho_af_stat: str,
                 clinvar_patho_exon_af_stat: str,
                 gnomAD_extreme_rare_threshold: float = 0.0001) -> np.ndarray:

    logger.info(f"Loading the clinvar pathogenic AF stat from {clinvar_patho_af_stat}")
    clinvar_patho_af_dict = pickle.load(gzip.open(clinvar_patho_af_stat)) if clinvar_patho_af_stat.endswith(".gz") else pickle.load(open(clinvar_patho_af_stat, 'rb'))
    clinvar_patho_af_dict = {k: v for k, v in clinvar_patho_af_dict.items() if v is not None}
    logger.info(f"The clinvar pathogenic AF stat type is {type(clinvar_patho_af_dict)}")
    df.loc[:, "Gene"] = df["Gene"].fillna(np.nan)
    gene_max_patho_af = df['Gene'].map(lambda gene: clinvar_patho_af_dict.get(gene, {"af": 0}).get('af', 0))
    df["clinvar_patho_gene_max_af"] = gene_max_patho_af

    clinvar_patho_exon_af_stat_dict = pickle.load(gzip.open(clinvar_patho_exon_af_stat)) if clinvar_patho_exon_af_stat.endswith(".gz") else pickle.load(open(clinvar_patho_exon_af_stat, 'rb'))
    df["exon_patho_median_af"] = df.apply(lambda row: clinvar_patho_exon_af_stat_dict.get(row['Feature'], {}).get(row["EXON"], (np.nan,))[0], axis=1)
    df["exon_patho_max_af"] = df.apply(lambda row: clinvar_patho_exon_af_stat_dict.get(row['Feature'], {}).get(row["EXON"], (np.nan,np.nan,np.nan))[2], axis=1)

    pm2_supporting_exon = gnomAD_rare_AF(df, "exon_patho_max_af")
    pm2_supporting_gene = gnomAD_rare_AF(df, "clinvar_patho_gene_max_af")

    pm2_supporting = np.where((df["exon_patho_max_af"] == 0) | df["exon_patho_max_af"].isna(), pm2_supporting_gene, pm2_supporting_exon)
    gnomAD_AF_rare = gnomAD_rare_AF(df, gnomAD_extreme_rare_threshold)
    pm2_supporting = np.where((df["clinvar_patho_gene_max_af"] == 0) | df["clinvar_patho_gene_max_af"].isna(), gnomAD_AF_rare, pm2_supporting)

    pm2_moderate = gnomAD_rare_AF(df, "exon_patho_median_af") & gnomAD_AF_rare
    gnomAD_absent = gnomAD_rare_AF(df, 1e-7)
    pm2_moderate = np.where((df["exon_patho_median_af"] == 0) | df["exon_patho_median_af"].isna(), gnomAD_absent, pm2_moderate)

    pm2_array = np.zeros(len(df), dtype=int)
    pm2_array[pm2_supporting] = 1
    pm2_array[pm2_moderate] = 2
    return pm2_array


def BS1_criteria(df: pd.DataFrame,
                 expected_incidence: float = 0.0001,
                 pm2_criteria: np.ndarray = None,
                 return_frequency_components: bool = False):
    '''
    BS1: allele frequency is greater than expected for the disorder.

    The upstream hub has already selected condition histories compatible with
    this variant's mechanism and added every included history whose mechanism
    is unresolved. BS1 therefore reads ``variant_inheritance`` and
    ``variant_penetrance`` directly; it must not recreate a gene-wide answer.

    Autosomal dominant frequency model:
      - dominant only: any variant state remains uncertain-compatible.
      - dominant_LOF only: predicted LOF or uncertain; exact GOF is excluded.
      - dominant_GOF only: exact GOF only.
      - dominant_DN only: uncertain only; PriVA has no exact DN database.
      - dominant_GOF + dominant_DN: exact GOF or uncertain; predicted LOF excluded.
      - dominant_LOF + dominant_GOF/DN: any variant state can be interpreted
        under at least one dominant history.
      - Any compatible biallelic requirement plus dominant history: do not use carrier AF for
        BS1, because heterozygous population observations may simply be
        recessive carriers. Use the recessive homozygous/hemizygous frequency
        model instead.

    Recessive frequency model:
      - A compatible biallelic requirement uses homozygous/hemizygous population frequency, not
        carrier allele count.

    X-linked frequency models:
      - ``x_linked_recessive`` uses the XY allele frequency. XX carrier
        frequency is not evidence against that model.
      - ``x_linked_dominant`` uses the XX allele frequency.
      - ``x_linked_unspecified`` supplies no disease-incidence BS1 model because
        the affected allele state is unknown.

    Gene pathogenic-frequency comparator:
      - For any valid Mendelian model, BS1 can also be assigned when the
        variant AF is significantly greater than the non-zero maximum AF of
        known pathogenic ClinVar variants in that gene.

    Global gates:
      - no BS1 for non-monogenic/polygenic, non-Mendelian, or any selected
        incomplete-penetrance-equivalent history.
      - no BS1 while a relevant HPO disease context still requires scope review.
      - PM2 blocks only the gene pathogenic-frequency comparator.
      - A disease-incidence BS1 assignment takes precedence over PM2.
    '''
    if pm2_criteria is None:
        pm2_present = np.zeros(len(df), dtype=bool)
    else:
        if len(pm2_criteria) != len(df):
            raise ValueError("pm2_criteria must have the same number of rows as df")
        pm2_present = np.asarray(pm2_criteria) > 0

    def _series_or_empty(column: str) -> pd.Series:
        return df.get(column, pd.Series("", index=df.index)).fillna("").astype(str)

    def _numeric_series(column: str) -> pd.Series:
        return pd.to_numeric(
            df.get(column, pd.Series(0, index=df.index)), errors="coerce"
        ).fillna(0)

    def _frequency_above(column_af: str, column_an: str, threshold: float) -> pd.Series:
        allele_frequency = _numeric_series(column_af)
        allele_number = _numeric_series(column_an)
        _, significant = control_false_neg_rate(
            allele_frequency,
            allele_number,
            af_threshold=threshold,
            alpha=0.01,
        )
        return pd.Series(
            np.where(significant.isna(), allele_frequency > threshold, significant),
            index=df.index,
            dtype=bool,
        )

    chrom = _series_or_empty("chrom").str.upper().str.removeprefix("CHR")
    x_linked = chrom.eq("X")
    y_locus = chrom.eq("Y")
    mitochondrial_locus = chrom.isin({"M", "MT"})
    autosomal = ~(x_linked | y_locus | mitochondrial_locus)

    condition_masks = _variant_condition_masks(df)
    scope_review_required = _series_or_empty("HPO_scope_review_required").str.lower().isin(
        {"1", "true", "yes"}
    )

    _, common_vars = control_false_neg_rate(df['gnomAD_joint_AF_max'], df['gnomAD_joint_AN_max'], af_threshold=expected_incidence, alpha=0.01)

    max_ind_incidence = np.where(df['gnomAD_joint_AN_max']/2 > 10/expected_incidence, df['gnomAD_nhomalt_max']/(df['gnomAD_joint_AN_max']/2), (df['gnomAD_nhomalt_XX'] + df['gnomAD_nhomalt_XY'])/(df['gnomAD_joint_AN']/2))
    max_af_larger_incidence = np.where(common_vars.isna(), df['gnomAD_joint_AF'] > expected_incidence, common_vars)
    logger.info(f"There are {max_af_larger_incidence.sum()} variants having their PAF greater than the expected incidence of the disease")

    has_recessive_requirement = condition_masks["has_recessive"]
    autosomal_dominant_without_recessive = (
        condition_masks["has_dominant"] & ~has_recessive_requirement
    )
    x_recessive_requirement = condition_masks["has_x_linked_recessive"]
    x_dominant_requirement = (
        condition_masks["has_x_linked_dominant"]
        & ~has_recessive_requirement
    )

    valid_model = (
        condition_masks["has_mendelian"]
        & np.logical_not(condition_masks["has_non_monogenic"])
        & np.logical_not(condition_masks["has_non_mendelian"])
        & np.logical_not(condition_masks["has_incomplete_penetrance"])
        & np.logical_not(scope_review_required)
    )

    x_xy_af_larger_incidence = _frequency_above(
        "gnomAD_joint_AF_XY", "gnomAD_joint_AN_XY", expected_incidence
    )
    x_xx_af_larger_incidence = _frequency_above(
        "gnomAD_joint_AF_XX", "gnomAD_joint_AN_XX", expected_incidence
    )

    # The two explicit X-linked inheritance models use different population
    # strata. A generic X-linked assertion supplies neither model.
    autosomal_dominant = (
        autosomal
        & autosomal_dominant_without_recessive
        & max_af_larger_incidence
        & valid_model
    )
    x_linked_recessive = (
        x_linked
        & x_recessive_requirement
        & x_xy_af_larger_incidence
        & valid_model
    )
    x_linked_dominant = (
        x_linked
        & x_dominant_requirement
        & x_xx_af_larger_incidence
        & valid_model
    )

    # Recessive BS1 uses homozygous/hemizygous frequency, not population carrier
    # count, because heterozygous carriers are expected for AR LoF disease.
    autosomal_recessive = autosomal & has_recessive_requirement & (max_ind_incidence > expected_incidence) & valid_model
    y_linked = (
        y_locus
        & condition_masks["has_y_linked"]
        & max_af_larger_incidence
        & valid_model
    )
    mitochondrial = (
        mitochondrial_locus
        & condition_masks["has_mitochondrial"]
        & max_af_larger_incidence
        & valid_model
    )
    greater_than_disease_incidence = (
        autosomal_dominant
        | autosomal_recessive
        | x_linked_recessive
        | x_linked_dominant
        | y_linked
        | mitochondrial
    )
    logger.info(f"There are {greater_than_disease_incidence.sum()} variants having their PAF greater than the expected incidence of the disease")
    gene_max_patho_af = df["clinvar_patho_gene_max_af"].fillna(0)

    _, greater_than_clinvar_patho_af = control_false_neg_rate(df['gnomAD_joint_AF_max'], df['gnomAD_joint_AN_max'], af_threshold=gene_max_patho_af, alpha=0.01)
    greater_than_clinvar_patho_af = np.where(greater_than_clinvar_patho_af.isna() & (gene_max_patho_af > 0), df['gnomAD_joint_AF'] > gene_max_patho_af, greater_than_clinvar_patho_af)
    greater_than_clinvar_patho_af = np.where(np.isnan(greater_than_clinvar_patho_af), False, greater_than_clinvar_patho_af)
    greater_than_clinvar_patho_af = greater_than_clinvar_patho_af & (gene_max_patho_af > 0)
    _, greater_than_basic_af = control_false_neg_rate(df['gnomAD_joint_AF_max'], df['gnomAD_joint_AN_max'], af_threshold=0.0001, alpha=0.01)
    greater_than_basic_af = np.where(greater_than_basic_af.isna(), df['gnomAD_joint_AF'] > 0.0001, greater_than_basic_af)
    x_xy_greater_than_basic_af = _frequency_above(
        "gnomAD_joint_AF_XY", "gnomAD_joint_AN_XY", 0.0001
    )
    x_xx_greater_than_basic_af = _frequency_above(
        "gnomAD_joint_AF_XX", "gnomAD_joint_AN_XX", 0.0001
    )
    clinvar_af_model_compatible = valid_model
    disease_incidence_bs1 = np.asarray(
        (
            (autosomal_dominant | autosomal_recessive | y_linked | mitochondrial)
            & greater_than_basic_af
        )
        | (x_linked_recessive & x_xy_greater_than_basic_af)
        | (
            x_linked_dominant
            & x_xx_greater_than_basic_af
        ),
        dtype=bool,
    )
    gene_pathogenic_af_bs1_candidate = np.asarray(
        greater_than_clinvar_patho_af
        & clinvar_af_model_compatible
        & greater_than_basic_af,
        dtype=bool,
    )
    gene_pathogenic_af_blocked_by_pm2 = gene_pathogenic_af_bs1_candidate & pm2_present
    gene_pathogenic_af_bs1 = gene_pathogenic_af_bs1_candidate & ~pm2_present
    bs1_criteria = disease_incidence_bs1 | gene_pathogenic_af_bs1
    bs1_array = np.zeros(len(df), dtype=int)
    bs1_array[bs1_criteria] = 3
    frequency_components = {
        "disease_incidence_bs1": disease_incidence_bs1,
        "gene_pathogenic_af_bs1": gene_pathogenic_af_bs1,
        "gene_pathogenic_af_blocked_by_pm2": gene_pathogenic_af_blocked_by_pm2,
    }
    if return_frequency_components:
        return bs1_array, frequency_components
    return bs1_array


def BA1_criteria(anno_df: pd.DataFrame) -> pd.Series:
    """
    Apply BA1 criteria to the annotation DataFrame.
    """
    false_neg_rate, common_vars = control_false_neg_rate(anno_df['gnomAD_joint_AF_max'], anno_df['gnomAD_joint_AN_max'], af_threshold=0.05, alpha=0.01)
    ba1_criteria = np.where(common_vars.isna(), anno_df['gnomAD_joint_AF'] > 0.05, common_vars)
    logger.info(f"BA1 criteria applied, {ba1_criteria.sum()} variants are having the BA1 criteria")

    ba1_array = np.zeros(len(anno_df), dtype=int)
    ba1_array[ba1_criteria] = 5
    return ba1_array
