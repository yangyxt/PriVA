#!/usr/bin/env python3
"""BS2 and BP5 -- observations about the people carrying the allele.

Neither is about the molecule. Both ask whether the carriers are consistent
with the allele causing the disease in question.

    BS2  observed in an apparently healthy adult, in a genotype the disease
         model would require to be affected. Late onset is already normalized
         to the linked condition's ``incomplete`` penetrance value,
         so a healthy thirty-year-old cannot contradict a late-onset disorder.
    BP5  found in a patient whose disease already has an alternative
         molecular explanation. extract_and_summarize_alt_disease_variants
         builds that cohort from the alternative-disease VCF.

Both are withheld when a selected history is non-monogenic, non-Mendelian or
incomplete-penetrance, because a healthy carrier or second diagnosis is then
not contradictory.
"""

import logging
import subprocess
import functools
import multiprocessing as mp
import pandas as pd
import numpy as np
from io import StringIO

from acmg_variant_mechanism import _variant_condition_masks


logger = logging.getLogger(__name__)


def BS2_criteria(
    df: pd.DataFrame,
    pm2_criteria: np.ndarray,
) -> np.ndarray:
    '''
    BS2: variant observed in apparently healthy individuals at a genotype dose
    incompatible with the inherited model.

    Returns an int array over ``df`` rows: 0 = no BS2, 3 = BS2. ``pm2_criteria``
    is required and must be an array over the same rows.

    Inheritance and penetrance are read from ``variant_inheritance`` and
    ``variant_penetrance`` through ``_variant_condition_masks``, which raises
    KeyError if either column is absent. Those columns already hold the
    condition histories the upstream hub attached to this variant, so BS2 does
    not recompute a gene-wide answer and consults no HPO array.

    Pathogenic mechanism is not an input. This module imports only
    ``_variant_condition_masks``, and no branch below distinguishes a
    loss-of-function variant from a gain-of-function or dominant-negative one.
    The question BS2 asks is whether a healthy person carries more copies of
    the allele than the inherited model tolerates, which is a statement about
    dose and inheritance alone.

    WHAT COUNTS AS A HEALTHY OBSERVATION
    ====================================

    Two sources contribute, and either is sufficient:

        gnomAD          for carrier counts, > 10 normally and > 3 when the
                        history asserts complete penetrance. For homozygote
                        counts the threshold is not a fixed number: it is the
                        count Hardy-Weinberg predicts at a carrier allele
                        frequency of 1e-3, that is cohort_individuals x 1e-6,
                        relaxed to 3e-4 (cohort_individuals x 9e-8) under
                        complete penetrance. At least one observed homozygote
                        is required in either case.
        control cohort  ``control_AC`` or ``control_nhomalt`` above zero, when
                        the user supplied a control VCF. Any observation counts;
                        no threshold and no penetrance relaxation apply.

    Why the homozygote threshold is stated as a frequency: homozygote count is
    roughly cohort size times allele frequency squared, so a flat "> 10" is a
    much stricter test than it looks. In a cohort of this size it demands a
    carrier frequency near 3.5e-3, about 35 times the 1e-4 the single-copy
    models use, and it moves with every gnomAD release. Naming the carrier
    frequency makes the requirement explicit and keeps it stable.

    The observation is read from a different stratum per model:

        carrier         ``gnomAD_joint_AF`` x ``gnomAD_joint_AN``, or
                        ``control_AC`` > 0, or a chrY variant with a non-zero
                        XY frequency.
        hom/hemi        ``gnomAD_nhomalt_XX`` + ``gnomAD_nhomalt_XY``, or
                        ``control_nhomalt`` > 0.
        XY hemizygous   ``gnomAD_nhomalt_XY`` alone.
        XX carrier      ``gnomAD_joint_AC_XX`` alone.

    WHICH OBSERVATION EACH MODEL ACCEPTS
    ====================================

        recessive history
            hom/hemi evidence only. A healthy heterozygote is expected under a
            biallelic model and is never BS2. This route is not gated on locus,
            so an autosomal-recessive history also admits hom/hemi evidence for
            a variant on chrX.
        autosomal + dominant history, no recessive requirement
            carrier or hom/hemi evidence. One copy in a healthy person already
            contradicts a dominant model.
        chrY + y_linked, chrM + mitochondrial
            same single-copy treatment as autosomal dominant.
        chrX + x_linked_recessive
            healthy XY hemizygotes only. XX carriers are expected.
        chrX + x_linked_dominant, no recessive requirement
            healthy XX carriers only. XY observations cannot substitute.

    A history whose locus does not match supplies no route, and
    ``x_linked_unspecified`` supplies none at all because the contradicting
    dose is unknown. When a gene carries both recessive and dominant history
    the recessive reading wins: carrier-only evidence is not BS2, while
    hom/hemi evidence still is.

    GATES
    =====

    A row must carry a Mendelian history and must not carry a non-monogenic
    (polygenic, digenic, oligogenic) history, a literal non-Mendelian history,
    an incomplete-penetrance history, or a set ``HPO_scope_review_required``
    flag. Upstream folds moderate penetrance, low penetrance and late onset
    into the ``incomplete`` token, so those block here too; other onset terms
    and variable expressivity do not.

    Penetrance therefore acts twice, in opposite directions. An ``incomplete``
    token blocks BS2 outright, because a healthy carrier of a partly penetrant
    allele is not contradictory. A ``complete`` token, which upstream also uses
    for high penetrance, relaxes the gnomAD thresholds: carrier counts from > 10
    to > 3, and the homozygote carrier-frequency basis from 1e-3 to 3e-4.

    Finally, BS2 is withdrawn wherever PM2 is set. The caller in
    ``acmg_criteria_assign`` passes the PM2 array after a disease-incidence BS1
    has already backed PM2 down, so a variant common enough for BS1 does not
    lose BS2 to the PM2 it no longer holds.

    Required columns
    ----------------
    Only ``variant_inheritance`` and ``variant_penetrance`` raise if absent.
    ``chrom``, ``HPO_scope_review_required``, every gnomAD column and both
    control columns default to empty or zero, so a missing frequency column
    silently yields no evidence and a missing ``chrom`` treats every locus as
    autosomal.
    '''
    def _series_or_empty(column: str) -> pd.Series:
        return df.get(column, pd.Series("", index=df.index)).fillna("").astype(str)

    def _numeric_column(column: str) -> pd.Series:
        return pd.to_numeric(df.get(column, pd.Series(0, index=df.index)), errors="coerce").fillna(0)

    chrom = _series_or_empty("chrom").str.upper().str.replace(
        r"^CHR", "", regex=True
    )
    x_locus = chrom.eq("X")
    y_locus = chrom.eq("Y")
    mitochondrial_locus = chrom.isin({"M", "MT"})
    autosomal_locus = ~(x_locus | y_locus | mitochondrial_locus)

    scope_review_required = _series_or_empty("HPO_scope_review_required").str.lower().isin(
        {"1", "true", "yes"}
    )
    condition_masks = _variant_condition_masks(df)
    complete_penetrance = condition_masks["has_complete_penetrance"]
    valid_model = (
        condition_masks["has_mendelian"]
        & np.logical_not(condition_masks["has_non_monogenic"])
        & np.logical_not(condition_masks["has_non_mendelian"])
    )
    baseline_eligible = (
        valid_model
        & np.logical_not(condition_masks["has_incomplete_penetrance"])
        & np.logical_not(scope_review_required)
    )
    complete_eligible = baseline_eligible & complete_penetrance

    has_recessive_requirement = condition_masks["has_recessive"]
    non_x_single_copy_model = (
        (
            autosomal_locus
            & condition_masks["has_dominant"]
            & ~has_recessive_requirement
        )
        | (y_locus & condition_masks["has_y_linked"])
        | (mitochondrial_locus & condition_masks["has_mitochondrial"])
    )
    x_recessive_model = x_locus & condition_masks["has_x_linked_recessive"]
    x_dominant_model = (
        x_locus
        & condition_masks["has_x_linked_dominant"]
        & ~has_recessive_requirement
    )

    gnomad_ac = (_numeric_column("gnomAD_joint_AF") * _numeric_column("gnomAD_joint_AN"))
    gnomad_hom_hemi = _numeric_column("gnomAD_nhomalt_XX") + _numeric_column("gnomAD_nhomalt_XY")
    gnomad_x_xy = _numeric_column("gnomAD_nhomalt_XY")
    gnomad_xx_ac = _numeric_column("gnomAD_joint_AC_XX")
    y_allele_observed = y_locus & (_numeric_column("gnomAD_joint_AF_XY") > 0)

    # The homozygote threshold is expressed as an equivalent carrier allele
    # frequency rather than a fixed count. A flat "> 10 homozygotes" is a far
    # stricter test than it appears: homozygote count is roughly cohort size
    # times AF squared, so > 10 in this cohort needs a carrier AF near 3.5e-3,
    # about 35 times the 1e-4 the single-copy models use. Deriving the count
    # from a stated carrier AF keeps the requirement explicit and makes it
    # scale with cohort size instead of with gnomAD's release.
    BS2_HOM_CARRIER_AF = 1e-3
    BS2_HOM_CARRIER_AF_COMPLETE = 3e-4
    gnomad_individuals = _numeric_column("gnomAD_joint_AN") / 2.0
    hom_hemi_min = gnomad_individuals * (BS2_HOM_CARRIER_AF ** 2)
    hom_hemi_min_complete = gnomad_individuals * (BS2_HOM_CARRIER_AF_COMPLETE ** 2)

    if 'control_AC' in df.columns:
        logger.warning("Seems user has provided a control cohort VCF, using control allele counts to assign BS2 criteria")
    if 'control_nhomalt' in df.columns:
        logger.warning("Seems user has provided a control cohort VCF, using control homozygous counts to assign BS2 criteria")
    control_ac_observed = _numeric_column("control_AC") > 0
    control_hom_observed = _numeric_column("control_nhomalt") > 0

    # A cohort with no called alleles yields a zero threshold, which any count
    # would clear, so require at least one observed homozygote as well.
    hom_hemi_evidence = (
        (gnomad_hom_hemi > hom_hemi_min) & (gnomad_hom_hemi > 0)
    ) | control_hom_observed
    carrier_evidence = (
        (gnomad_ac > 10) | control_ac_observed | y_allele_observed
    )
    complete_hom_hemi_evidence = (
        (gnomad_hom_hemi > hom_hemi_min_complete) & (gnomad_hom_hemi > 0)
    ) | control_hom_observed
    complete_carrier_evidence = (
        (gnomad_ac > 3) | control_ac_observed | y_allele_observed
    )
    logger.info(
        "BS2 homozygote thresholds from carrier AF %g (complete: %g): median required "
        "count %.2f (complete %.2f); rows clearing it: %d standard, %d complete",
        BS2_HOM_CARRIER_AF, BS2_HOM_CARRIER_AF_COMPLETE,
        float(np.nanmedian(hom_hemi_min)) if len(hom_hemi_min) else float("nan"),
        float(np.nanmedian(hom_hemi_min_complete)) if len(hom_hemi_min_complete) else float("nan"),
        int(hom_hemi_evidence.sum()), int(complete_hom_hemi_evidence.sum()),
    )
    x_xy_evidence = gnomad_x_xy > 10
    complete_x_xy_evidence = gnomad_x_xy > 3
    x_xx_evidence = gnomad_xx_ac > 10
    complete_x_xx_evidence = gnomad_xx_ac > 3
    bs2_standard = baseline_eligible & (
        (has_recessive_requirement & hom_hemi_evidence)
        | (non_x_single_copy_model & (carrier_evidence | hom_hemi_evidence))
        | (x_recessive_model & x_xy_evidence)
        | (x_dominant_model & x_xx_evidence)
    )
    bs2_complete_penetrance = complete_eligible & (
        (has_recessive_requirement & complete_hom_hemi_evidence)
        | (x_recessive_model & complete_x_xy_evidence)
        | (
            x_dominant_model
            & complete_x_xx_evidence
        )
        | (
            non_x_single_copy_model
            & (complete_carrier_evidence | complete_hom_hemi_evidence)
        )
    )
    bs2_criteria = bs2_standard | bs2_complete_penetrance
    bs2_criteria = bs2_criteria & (pm2_criteria == 0)
    bs2_array = np.zeros(len(df), dtype=int)
    bs2_array[bs2_criteria] = 3
    return bs2_array


def extract_and_summarize_alt_disease_variants(alt_disease_vcf: str, chrom = None) -> pd.DataFrame:
    '''
    Extract variants from alt_disease_vcf using bcftools and summarize genotype information.

    Args:
        alt_disease_vcf: Path to VCF file with known alternative molecular basis for disease

    Returns:
        DataFrame with columns: chrom, pos, ref, alt, alt_disease_hets, alt_disease_homs
    '''
    # Use bcftools to query genotype
    if chrom is None:
        cmd = f'bcftools +fill-tags {alt_disease_vcf} -Ou -- -t AC_Het,AC_Hom,AC_Hemi | bcftools filter -e \'(INFO/AC_Het == 0) && (INFO/AC_Hom == 0) && (INFO/AC_Hemi == 0)\' -Ou | bcftools query -f "%CHROM:%POS:%REF-%ALT\\t%INFO/AC_Het\\t%INFO/AC_Hom\\t%INFO/AC_Hemi\\n"'
    else:
        cmd = f'bcftools view -Ou {alt_disease_vcf} {chrom} | bcftools +fill-tags - -Ou -- -t AC_Het,AC_Hom,AC_Hemi | bcftools filter -e \'(INFO/AC_Het == 0) && (INFO/AC_Hom == 0) && (INFO/AC_Hemi == 0)\' -Ou | bcftools query -f "%CHROM:%POS:%REF-%ALT\\t%INFO/AC_Het\\t%INFO/AC_Hom\\t%INFO/AC_Hemi\\n"'

    logger.info(f"Extracting and normalizing variants from {alt_disease_vcf} using bcftools...")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)

    if not result.stdout.strip():
        logger.warning(f"No variants found in {alt_disease_vcf}")
        return {}, {}

    # Read bcftools output directly into DataFrame
    df = pd.read_csv(StringIO(result.stdout), sep='\t', header=None)

    # Set column names - first 4 are variant info, rest are genotypes
    col_names = ['variant_id'] + ['AC_Het', 'AC_Hom', 'AC_Hemi']
    df.columns = col_names

    # Create boolean masks across all samples (vectorized)
    het_mask = df['AC_Het'] > 0
    hom_mask = (df['AC_Hom'] > 0) | (df['AC_Hemi'] > 0)

    df['alt_disease_hets'] = het_mask
    df['alt_disease_homs'] = hom_mask

    # Return only the required columns
    result_df = df[['variant_id', 'alt_disease_hets', 'alt_disease_homs']]
    # Convert the df to 2 dicts efficiently, one is het, one is hom
    result_df.set_index('variant_id', inplace=True)
    het_dict = result_df['alt_disease_hets'].to_dict()
    hom_dict = result_df['alt_disease_homs'].to_dict()

    return het_dict, hom_dict


def BP5_criteria(df: pd.DataFrame,
                 alt_disease_vcf: str,
                 n_processes: int = 1) -> np.ndarray:
    '''
    BP5: variant found in a sample with known alternative molecular basis for disease

    Logic: If a variant is found in patients with alternative diseases, it suggests benignity
    for the query disease, UNLESS the gene has non-monogenic inheritance or incomplete penetrance.

    For selected autosomal dominant or X-linked-dominant histories:
    - Variants found in alt disease patients (het or hom) get BP5

    For selected recessive or X-linked-recessive histories:
    - Only variants found as homozygous or hemizygous get BP5

    The inheritance and penetrance inputs are the hub's condition-linked
    categorical columns. Known incompatible mechanisms were already removed;
    unresolved condition histories remain, so this is the highest available
    resolution without pretending that an unknown mechanism is known.

    THE DECISION TREE
    =================

    BP5 reads two variant-condition columns:

        variant_inheritance / variant_penetrance
              |
              v
        eligible_genes = NOT non_monogenic
                       & NOT non_mendelian
                       & NOT incomplete_penetrance
              |
              |   a genotype contradiction in a gene with incomplete
              |   penetrance or a non-Mendelian model is not interpretable,
              |   which is why all three are hard gates rather than weights
              v
        variant seen in a patient whose disease has another established
        molecular basis
              |
              +-- dominant gene  -> heterozygous OR homozygous carrier -> BP5
              +-- recessive gene -> homozygous carrier ONLY            -> BP5
                                    (a heterozygous carrier of a recessive
                                     gene is unremarkable and says nothing)

    ``variant_inheritance_basis`` records whether the values came from matched,
    unresolved, or mixed histories. It is audit provenance, not a reason to
    change this genotype decision tree.
    '''
    if {"alt_disease_hets", "alt_disease_homs"}.issubset(df.columns):
        logger.info(
            "Reusing precomputed alt_disease_hets/alt_disease_homs columns for BP5; "
            "skip re-scanning %s",
            alt_disease_vcf,
        )
        het_col = df["alt_disease_hets"]
        hom_col = df["alt_disease_homs"]
        if het_col.dtype == object:
            df["alt_disease_hets"] = het_col.astype(str).str.lower().isin({"true", "1", "yes"})
        if hom_col.dtype == object:
            df["alt_disease_homs"] = hom_col.astype(str).str.lower().isin({"true", "1", "yes"})
    else:
        # First, efficiently extract unique chromosomes from alt_disease_vcf
        logger.info(f"Extracting unique chromosomes from {alt_disease_vcf}")
        chrom_cmd = f'bcftools query -f "%CHROM\\n" {alt_disease_vcf} | sort -u'
        chrom_result = subprocess.run(chrom_cmd, shell=True, capture_output=True, text=True, check=True)

        if not chrom_result.stdout.strip():
            logger.warning(f"No chromosomes found in {alt_disease_vcf}")
            return np.zeros(len(df), dtype=int)

        unique_chroms = [chrom.strip() for chrom in chrom_result.stdout.strip().split('\n') if chrom.strip()]
        logger.info(f"Found {len(unique_chroms)} unique chromosomes: {unique_chroms}")

        # Create partial function with fixed alt_disease_vcf parameter
        extract_func = functools.partial(extract_and_summarize_alt_disease_variants, alt_disease_vcf)

        logger.info(f"Using {n_processes} processes for parallel chromosome processing")

        # Execute in parallel
        with mp.Pool(processes=n_processes) as pool:
            results = pool.map(extract_func, unique_chroms)

        # Combine results from all chromosomes
        combined_het_dict = {}
        combined_hom_dict = {}

        for het_dict, hom_dict in results:
            combined_het_dict.update(het_dict)
            combined_hom_dict.update(hom_dict)

        logger.info(f"Combined results: {len(combined_het_dict)} het variants, {len(combined_hom_dict)} hom variants")

        # Create variant_id column for merging (matching the format used in extract function)
        df['variant_id'] = df['chrom'].astype(str) + ':' + df['pos'].astype(str) + ':' + df['ref'] + '-' + df['alt']

        # Map the dict results back to the dataframe
        df['alt_disease_hets'] = df['variant_id'].map(combined_het_dict).fillna(False)
        df['alt_disease_homs'] = df['variant_id'].map(combined_hom_dict).fillna(False)

    # Only an included monogenic Mendelian history without an incomplete-
    # equivalent condition makes an alternative diagnosis contradictory.
    condition_masks = _variant_condition_masks(df)
    eligible_genes = (
        condition_masks["has_mendelian"]
        & ~condition_masks["has_non_monogenic"]
        & ~condition_masks["has_non_mendelian"]
        & ~condition_masks["has_incomplete_penetrance"]
    )

    # A heterozygous observation is informative for autosomal dominant or
    # explicitly X-linked-dominant disease. An X-linked-recessive history does
    # not erase a separate dominant assertion for the same condition.
    heterozygous_disease = (
        condition_masks["has_dominant_like"]
        & ~condition_masks["has_recessive"]
    )
    dominant_bp5 = (heterozygous_disease & \
                   eligible_genes & \
                   (df['alt_disease_hets'] | df['alt_disease_homs']))

    # Biallelic and X-linked-recessive models require a homozygous/hemizygous
    # observation. ``alt_disease_homs`` includes AC_Hemi from the input VCF.
    hom_hemi_disease = (
        condition_masks["has_recessive"]
        | condition_masks["has_x_linked_recessive"]
    )
    recessive_bp5 = hom_hemi_disease & \
                    eligible_genes & \
                    df['alt_disease_homs']

    # Combine conditions
    bp5_criteria = dominant_bp5 | recessive_bp5

    bp5_array = np.zeros(len(df), dtype=int)
    bp5_array[bp5_criteria] = 1

    logger.info(
        "BP5 criteria evaluated for %s variants with an interpretable "
        "monogenic Mendelian history",
        int(eligible_genes.sum()),
    )
    logger.info(f"BP5 assigned to {bp5_criteria.sum()} variants found in alternative disease patients")

    return bp5_array
