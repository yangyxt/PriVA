#!/usr/bin/env python3
"""BS2 and BP5 -- observations about the people carrying the allele.

Neither is about the molecule. Both ask whether the carriers are consistent
with the allele causing the disease in question.

    BS2  observed in an apparently healthy adult, in a genotype the disease
         model would require to be affected. Whether a carrier counts as
         informative depends on the disorder's age of onset, which is why
         this reads hpo_onset_modes from acmg_inheritance -- a healthy
         thirty-year-old says nothing about a late-onset disorder.
    BP5  found in a patient whose disease already has an alternative
         molecular explanation. extract_and_summarize_alt_disease_variants
         builds that cohort from the alternative-disease VCF.

Both are withheld for non-monogenic, non-Mendelian and incomplete-penetrance
genes, where a healthy carrier or a second diagnosis is not contradictory.
"""

import logging
import subprocess
import functools
import multiprocessing as mp
import pandas as pd
import numpy as np
from io import StringIO

from acmg_inheritance import hpo_onset_modes
from acmg_variant_mechanism import _variant_mechanism_masks


logger = logging.getLogger(__name__)


def BS2_criteria(df: pd.DataFrame,
                 non_monogenic: np.ndarray,
                 non_mendelian: np.ndarray,
                 incomplete_penetrance: np.ndarray,
                 pm2_criteria: np.ndarray) -> pd.Series:
    '''
    BS2: variant observed in apparently healthy adults for a disease expected
    to be penetrant at the observed age. Output keeps the existing PriVA
    strength encoding:
      0 = no BS2, 3 = BS2.

    The decision requires the upstream variant-effect and mechanism-
    applicability fields. There is no legacy mechanism fallback.

    Global gates:
      - no BS2 for non-monogenic/polygenic, non-Mendelian, incomplete-
        penetrance, moderate/low-penetrance, or late-onset-compatible rows.
      - no BS2 while a relevant HPO disease context still requires scope review.
      - if HPO explicitly marks complete penetrance, use lower gnomAD count
        thresholds (>3 instead of >10), still requiring no incomplete
        penetrance and no late-onset confounder.
      - final BS2 is removed when PM2 is assigned.

    Healthy observation classes:
      - carrier evidence = gnomAD allele count above threshold, internal
        control_AC > 0, or healthy Y-linked allele observation.
      - hom/hemi evidence = gnomAD_nhomalt_XX + gnomAD_nhomalt_XY above
        threshold, or internal control_nhomalt > 0.

    BS2 decision matrix:
      Variant state:
        predicted_LOF = LOFTEE HC/OS, NMD pLoF, or PriVA splice/UTR LOF.
        exact_GOF = variant-level exact GoFCards GOF match.
        uncertain = neither predicted_LOF nor exact_GOF.

      Important DN limitation:
        dominant_DN is currently curated condition history only. PriVA has no
        variant-level DN database, so there is no exact_DN variant state. DN
        history is treated as compatible with ambiguous variants only; NMD_LOF
        and exact_GOF are not treated as DN-compatible.

      Matrix values are PriVA strengths: 0 = no BS2, 3 = BS2.
      Every "no healthy evidence" combination is 0.

      mechanism profile                   carrier NMD  carrier GOF  carrier amb  hom NMD  hom GOF  hom amb
      no usable history                             0            0            0        0        0        0
      recessive only                                0            0            0        3        3        3
      dominant only                                 3            3            3        3        3        3
      dominant_LOF only                             3            0            3        3        0        3
      dominant_GOF only                             0            3            0        0        3        0
      dominant_DN only                              0            0            3        0        0        3
      dominant_GOF + dominant_DN                    0            3            3        0        3        3
      dominant_LOF + dominant_GOF/DN                3            3            3        3        3        3
      recessive + any dominant history              0            0            0        3        3        3

      Consequences:
        - GOF_only_NMD_control_AC and DN_only_NMD_control_AC are 0.
        - GOF_only_ambiguous_control_AC is 0.
        - DN_only_ambiguous_control_AC is 3.
        - DN-only ambiguous carrier/hom evidence is BS2-compatible; NMD_LOF is
          not DN-compatible, and exact_GOF is not treated as DN evidence.
        - dominant only means dominant inheritance without a curated
          LOF/GOF/DN mechanism and without recessive history; healthy carrier or
          hom/hemi evidence is BS2-compatible regardless of query consequence.
        - dominant_LOF healthy carriers are BS2 unless the exact variant is already
          ascertained as GOF.
        - If recessive history is present, carrier-only evidence is not BS2 even
          when dominant history is also present; hom/hemi evidence remains BS2.
    '''
    def _series_or_empty(column: str) -> pd.Series:
        return df.get(column, pd.Series("", index=df.index)).fillna("").astype(str)

    def _numeric_column(column: str) -> pd.Series:
        return pd.to_numeric(df.get(column, pd.Series(0, index=df.index)), errors="coerce").fillna(0)

    chrom = _series_or_empty("chrom")
    y_locus = chrom.str.contains("Y", regex=False)

    hpo_ids = _series_or_empty("HPO_IDs")
    scope_review_required = _series_or_empty("HPO_scope_review_required").str.lower().isin(
        {"1", "true", "yes"}
    )
    not_late_onsets = hpo_ids.map(hpo_onset_modes)
    complete_penetrance = hpo_ids.str.contains("HP:0034950", regex=False)  # Complete penetrance
    moderate_to_low_penetrance = (
        hpo_ids.str.contains("HP:4000159", regex=False)
        | hpo_ids.str.contains("HP:4000160", regex=False)
    )
    valid_model = pd.Series(
        np.logical_not(non_monogenic)
        & np.logical_not(non_mendelian),
        index=df.index,
    )
    baseline_eligible = (
        valid_model
        & np.logical_not(incomplete_penetrance)
        & np.logical_not(moderate_to_low_penetrance)
        & not_late_onsets
        & np.logical_not(scope_review_required)
    )
    complete_eligible = (
        valid_model
        & np.logical_not(incomplete_penetrance)
        & complete_penetrance
        & not_late_onsets
        & np.logical_not(scope_review_required)
    )

    mechanism_masks = _variant_mechanism_masks(df)
    has_recessive_requirement = mechanism_masks["has_recessive_compatible"]
    dominant_mechanism_compatible = mechanism_masks["has_dominant_compatible"]

    gnomad_ac = (_numeric_column("gnomAD_joint_AF") * _numeric_column("gnomAD_joint_AN"))
    gnomad_hom_hemi = _numeric_column("gnomAD_nhomalt_XX") + _numeric_column("gnomAD_nhomalt_XY")
    y_allele_observed = y_locus & (_numeric_column("gnomAD_joint_AF_XY") > 0)

    if 'control_AC' in df.columns:
        logger.warning("Seems user has provided a control cohort VCF, using mechanism-aware control allele counts to assign BS2 criteria")
    if 'control_nhomalt' in df.columns:
        logger.warning("Seems user has provided a control cohort VCF, using control homozygous counts to assign BS2 criteria")
    control_ac_observed = _numeric_column("control_AC") > 0
    control_hom_observed = _numeric_column("control_nhomalt") > 0

    hom_hemi_evidence = (
        ((gnomad_hom_hemi > 10) | control_hom_observed)
        | (complete_penetrance & (gnomad_hom_hemi > 3))
    )
    carrier_evidence = (
        ((gnomad_ac > 10) | control_ac_observed | y_allele_observed)
        | (complete_penetrance & ((gnomad_ac > 3) | y_allele_observed))
    )
    bs2_standard = baseline_eligible & (
        (has_recessive_requirement & hom_hemi_evidence)
        | (dominant_mechanism_compatible & (carrier_evidence | hom_hemi_evidence))
    )
    bs2_complete_penetrance = complete_eligible & (
        (has_recessive_requirement & hom_hemi_evidence)
        | (dominant_mechanism_compatible & (carrier_evidence | hom_hemi_evidence))
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
                 recessive: np.ndarray,
                 dominant: np.ndarray,
                 non_monogenic: np.ndarray,
                 non_mendelian: np.ndarray,
                 incomplete_penetrance: np.ndarray,
                 n_processes: int = 1) -> np.ndarray:
    '''
    BP5: variant found in a sample with known alternative molecular basis for disease

    Logic: If a variant is found in patients with alternative diseases, it suggests benignity
    for the query disease, UNLESS the gene has non-monogenic inheritance or incomplete penetrance.

    For genes with dominant inheritance + monogenic + complete penetrance:
    - Variants found in alt disease patients (het or hom) get BP5

    For genes with recessive inheritance + monogenic + complete penetrance:
    - Only variants found as homozygous in alt disease patients get BP5

    Optimized version using pre-computed inheritance arrays and vectorized operations

    THE DECISION TREE
    =================

    BP5 reads five gene-level arrays and no variant-level input. Like PP1, the
    inheritance it reasons about belongs to the gene, not to this variant.

        recessive / dominant / non_monogenic / non_mendelian /
        incomplete_penetrance          <- from identify_inheritance_mode
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

    WHAT IT WOULD READ INSTEAD
    ==========================

        variant_inheritance          replaces recessive / dominant
        variant_inheritance_basis    says whether that came from this
                                     variant's own matched history, the
                                     gene's unanimous consensus, or the
                                     gene's constraint data
        variant_penetrance           replaces incomplete_penetrance
        the non-Mendelian inheritance values replace non_monogenic and
        non_mendelian
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

    # Vectorized BP5 logic
    # Only apply BP5 to genes with monogenic inheritance and complete penetrance
    eligible_genes = ~non_monogenic & ~incomplete_penetrance & ~non_mendelian

    # For dominant genes: any presence (het or hom) in alt disease patients = BP5
    only_dominant = (dominant | incomplete_penetrance) & ~recessive
    dominant_bp5 = (only_dominant & \
                   eligible_genes & \
                   (df['alt_disease_hets'] | df['alt_disease_homs']))

    # For recessive genes: only homozygous presence in alt disease patients = BP5
    recessive_bp5 = (recessive | ~incomplete_penetrance) & \
                    eligible_genes & \
                    df['alt_disease_homs']

    # Combine conditions
    bp5_criteria = dominant_bp5 | recessive_bp5

    bp5_array = np.zeros(len(df), dtype=int)
    bp5_array[bp5_criteria] = 1

    logger.info(f"BP5 criteria applied to {eligible_genes.sum()} variants with monogenic + complete penetrance")
    logger.info(f"BP5 assigned to {bp5_criteria.sum()} variants found in alternative disease patients")

    return bp5_array
