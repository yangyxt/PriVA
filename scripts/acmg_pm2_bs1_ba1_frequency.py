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
    # Three outcomes, not two. A plain ``p_values <= alpha`` would fold "the
    # test could not be run" into "the test said no", because a comparison
    # against NaN is False. Callers need to tell those apart so they can fall
    # back to a direct frequency comparison, so the undecidable rows carry NA.
    reject_h0 = (p_values <= alpha).astype("boolean").mask(p_values.isna())
    logger.info(
        f"There are {int(reject_h0.sum())} variants having their p-value less than {alpha}, "
        f"{int(reject_h0.isna().sum())} datapoints are NA (no allele number or no frequency, "
        "test undecidable)"
    )

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

    Returns an int array over ``df`` rows: 0 = no BS1, 3 = BS1. When
    ``return_frequency_components`` is true, returns that array plus a dict of
    the three per-route boolean masks described under "Returned components".

    Inheritance and penetrance are read from ``variant_inheritance`` and
    ``variant_penetrance`` through ``_variant_condition_masks``, which raises
    KeyError if either column is absent. Those columns already hold the
    condition histories the upstream hub attached to this variant, so BS1 does
    not recompute a gene-wide answer and consults no HPO array.

    Pathogenic mechanism is not an input, and this function never calls
    ``_variant_mechanism_masks``. What BS1 tests is a maximum credible allele
    frequency, whose parameters are disease incidence, heterogeneity and
    penetrance. Whether an allele acts by loss of function, gain of function or
    a dominant-negative effect does not change how common it may be among
    unaffected people.

    HOW SIGNIFICANCE IS DECIDED
    ===========================

    Every frequency comparison except the autosomal-recessive one goes through
    ``control_false_neg_rate``: a one-sided binomial test of the observed
    allele count against the route's threshold, at alpha = 0.01. The point
    estimate being above the threshold is not sufficient; the excess must be
    significant given the allele number, so small-cohort strata cannot produce
    BS1 on their own.

    That test returns three outcomes, not two: yes, no, and undecidable. A row
    whose allele number is zero, or whose frequency is missing, cannot be
    tested at all and comes back undecidable rather than counting as a "no".
    Each route then falls back to comparing ``gnomAD_joint_AF`` against its
    threshold directly, unweighted by any sample size. A row with neither a
    testable stratum nor a joint frequency gets no BS1.

    The one exception to the fallback is the gene pathogenic-frequency route:
    it falls back only where the gene actually has a pathogenic frequency above
    zero to be compared against, since a gene with none is not a candidate for
    that route in the first place.

    The autosomal-recessive route uses a one-sided binomial test on the
    observed homozygote count against ``expected_incidence``, like the other
    incidence-based routes.  The count is read from the max-population stratum
    when that stratum is large enough to be trusted, otherwise from the joint
    sample.

    TWO INDEPENDENT ROUTES TO BS1
    =============================

    BS1 fires if either route fires. Both additionally require the variant to
    clear an absolute floor of AF > 0.0001, tested in the same stratum the
    route uses. With the default ``expected_incidence`` of 0.0001 that floor
    coincides with the incidence threshold; passing a rarer incidence leaves
    the floor in force as the binding constraint.

    Route 1, disease incidence. The locus and the inherited model must agree,
    and each pairing reads its own population stratum:

        autosomal + dominant history, no biallelic requirement
            carrier allele frequency (``gnomAD_joint_AF``/``AN``). A dominant
            history that also carries any recessive requirement is excluded
            here, because population heterozygotes may simply be recessive
            carriers; such genes are served by the recessive route instead.
        autosomal + recessive history
            observed homozygote count, not carrier count. Uses
            ``gnomAD_nhomalt_max`` over the max-population sample when that
            sample has at least 10,000 called individuals (AN_max/2 > 10,000),
            otherwise ``nhomalt_XX + nhomalt_XY`` over the joint sample.  The
            stratum-selection gate is fixed and does not change with
            ``expected_incidence``.  BS1 fires when the count is significantly
            too high under the one-sided binomial test
            P(X >= x | n, p=expected_incidence) at alpha = 0.01 and the
            observed homozygote frequency exceeds the 0.0001 floor.
        chrX + x_linked_recessive
            the XY stratum. XX carrier frequency is not evidence against this
            model.
        chrX + x_linked_dominant, no biallelic requirement
            the XX stratum.
        chrY + y_linked, chrM + mitochondrial
            carrier allele frequency.

    ``x_linked_unspecified`` supplies no route, because the affected allele
    state is unknown. A history whose locus does not match, such as a dominant
    history on chrX, likewise supplies no route.

    Route 2, gene pathogenic-frequency comparator. For any valid Mendelian
    model, BS1 also fires when the variant frequency significantly exceeds
    ``clinvar_patho_gene_max_af``, the maximum frequency among known pathogenic
    ClinVar variants in that gene. Genes with no such variant have the value 0
    and are excluded. This route uses the joint max stratum regardless of
    inheritance, and it is the only route PM2 can block.

    GATES APPLIED TO BOTH ROUTES
    ============================

    A row must carry a Mendelian history and must not carry a non-monogenic
    (polygenic, digenic, oligogenic) history, a literal non-Mendelian history,
    an incomplete-penetrance history, or a set ``HPO_scope_review_required``
    flag. Upstream folds moderate penetrance, low penetrance and late onset
    into the ``incomplete`` token, so those block here too.

    RELATIONSHIP WITH PM2
    =====================

    ``pm2_criteria`` is optional; omitting it treats PM2 as absent everywhere.
    Where PM2 is present it suppresses route 2 only, leaving route 1 intact.
    BS1 does not modify PM2 itself. The reverse precedence, in which a
    disease-incidence BS1 backs PM2 down to 0, is applied by the caller in
    ``acmg_criteria_assign`` using the components below.

    Returned components
    -------------------
        disease_incidence_bs1               route 1 fired
        gene_pathogenic_af_bs1              route 2 fired and PM2 was absent
        gene_pathogenic_af_blocked_by_pm2   route 2 would have fired but PM2
                                            was present

    Required columns
    ----------------
    ``variant_inheritance``, ``variant_penetrance``, ``gnomAD_joint_AF``,
    ``gnomAD_joint_AF_max``, ``gnomAD_joint_AN``, ``gnomAD_joint_AN_max``,
    ``gnomAD_nhomalt_max``, ``gnomAD_nhomalt_XX``, ``gnomAD_nhomalt_XY`` and
    ``clinvar_patho_gene_max_af`` raise if absent. ``chrom``,
    ``HPO_scope_review_required`` and the four XX/XY frequency columns default
    to empty or zero when absent, which silently removes the X-linked routes
    and treats every locus as autosomal.
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

    # The max-population homozygote estimate is only used when that stratum
    # has enough called individuals to be trusted.  This is a sample-size gate
    # for the estimator, not part of the incidence threshold, so it is fixed.
    hom_max_pop_min_individuals = 10_000  # AN_max/2, i.e. called individuals
    use_max_pop_hom = (df['gnomAD_joint_AN_max'] / 2) > hom_max_pop_min_individuals
    hom_n_individuals = pd.Series(
        np.where(use_max_pop_hom, df['gnomAD_joint_AN_max'] / 2, df['gnomAD_joint_AN'] / 2),
        index=df.index,
    )
    hom_observed = pd.Series(
        np.where(use_max_pop_hom, df['gnomAD_nhomalt_max'], df['gnomAD_nhomalt_XX'] + df['gnomAD_nhomalt_XY']),
        index=df.index,
    )
    hom_n_individuals = pd.to_numeric(hom_n_individuals, errors='coerce').fillna(0).clip(lower=0).astype(int)
    hom_observed = pd.to_numeric(hom_observed, errors='coerce').fillna(0).clip(lower=0).astype(int)

    # One-sided binomial test: P(X >= x | n, p=expected_incidence).  BS1 fires
    # only when the observed homozygote count is too high to be a chance
    # observation under the rare-disease model.  This controls the Type I
    # error -- falsely calling a truly rare variant common.
    hom_freq = hom_observed / hom_n_individuals.replace(0, np.nan)
    n_arr = hom_n_individuals.to_numpy(dtype=int)
    x_arr = hom_observed.to_numpy(dtype=int)
    pval_arr = np.full(len(df), np.nan)
    testable = n_arr > 0
    if testable.any():
        k = np.maximum(x_arr[testable] - 1, 0)
        pval_arr[testable] = binom.sf(k, n_arr[testable], expected_incidence)
        pval_arr[testable & (x_arr == 0)] = 1.0
    recessive_hom_bs1 = (hom_freq > 0.0001) & (pval_arr <= 0.01)
    # Where the max-population stratum had no allele number to test, fall back to
    # the joint frequency. fillna(False) then covers a row that has neither.
    max_af_larger_incidence = (
        common_vars
        .mask(common_vars.isna(), df['gnomAD_joint_AF'] > expected_incidence)
        .fillna(False)
        .astype(bool)
    )
    logger.info(f"There are {int(max_af_larger_incidence.sum())} variants having their PAF greater than the expected incidence of the disease")

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

    # Recessive BS1 uses the homozygous/hemizygous count, not population carrier
    # count, because heterozygous carriers are expected for AR LoF disease.
    autosomal_recessive = autosomal & has_recessive_requirement & recessive_hom_bs1 & valid_model
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
    # Fall back only where the gene actually has a pathogenic frequency to beat;
    # a gene with none is not a candidate for this route, so an undecidable test
    # there stays False rather than being compared against zero.
    greater_than_clinvar_patho_af = (
        greater_than_clinvar_patho_af
        .mask(
            greater_than_clinvar_patho_af.isna() & (gene_max_patho_af > 0),
            df['gnomAD_joint_AF'] > gene_max_patho_af,
        )
        .fillna(False)
        .astype(bool)
    )
    greater_than_clinvar_patho_af = greater_than_clinvar_patho_af & (gene_max_patho_af > 0)
    _, greater_than_basic_af = control_false_neg_rate(df['gnomAD_joint_AF_max'], df['gnomAD_joint_AN_max'], af_threshold=0.0001, alpha=0.01)
    greater_than_basic_af = (
        greater_than_basic_af
        .mask(greater_than_basic_af.isna(), df['gnomAD_joint_AF'] > 0.0001)
        .fillna(False)
        .astype(bool)
    )
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


def BA1_criteria(anno_df: pd.DataFrame) -> np.ndarray:
    """
    Apply BA1 criteria to the annotation DataFrame.

    Returns an int array over ``anno_df`` rows: 0 = no BA1, 5 = BA1. A variant
    is BA1 when its frequency is significantly above 5%. Where the max-population
    stratum has no allele number to test, the joint frequency is compared to 5%
    directly; where neither is available the variant is not BA1.
    """
    false_neg_rate, common_vars = control_false_neg_rate(anno_df['gnomAD_joint_AF_max'], anno_df['gnomAD_joint_AN_max'], af_threshold=0.05, alpha=0.01)
    ba1_criteria = (
        common_vars
        .mask(common_vars.isna(), anno_df['gnomAD_joint_AF'] > 0.05)
        .fillna(False)
        .astype(bool)
    )
    logger.info(f"BA1 criteria applied, {int(ba1_criteria.sum())} variants are having the BA1 criteria")

    ba1_array = np.zeros(len(anno_df), dtype=int)
    ba1_array[ba1_criteria] = 5
    return ba1_array
