#!/usr/bin/env python3
"""PP3, BP4 and BP7 -- what the prediction tools say.

    PP3  multiple lines of computational evidence support a damaging effect
    BP4  multiple lines support no impact
    BP7  a synonymous change with no predicted splicing effect at a
         position that is not conserved

Predictors are read per consequence class, because none of them is calibrated
across all of them: AlphaMissense and PrimateAI for missense, SpliceAI and
SpliceVault for splicing, UTRAnnotator for 5' untranslated changes, and CADD
plus GERP++ conservation throughout.
"""

import logging
import pandas as pd
import numpy as np


logger = logging.getLogger(__name__)


def PP3_BP4_criteria(df: pd.DataFrame, pvs1_criteria: np.ndarray = None, high_confidence_status = {
        # Higher confidence (2+ stars)
        'practice_guideline': 4,                                   # 4 stars
        'reviewed_by_expert_panel': 3,                             # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts': 2,  # 2 stars
    }) -> pd.Series:
    # PP3: predicted to be deleterious by in-silico tools
    # Including PrimateAI (Missense), CADD, AlphaMissense (Missense), VEP, SpliceAI, SpliceVault, UTRAnnotator

    missense_variant = df['Consequence'].str.contains('missense_variant') & np.logical_not(df['Consequence'].str.contains('stop')) & np.logical_not(df['Consequence'].str.contains('frameshift'))
    splice_variant = df['Consequence'].str.contains('splice') & np.logical_not(df['Consequence'].str.contains('stop')) & np.logical_not(df['Consequence'].str.contains('frameshift'))
    five_utr_variant = df['Consequence'].str.contains('5_prime_UTR_variant').fillna(False)
    missense_variant = missense_variant.fillna(False)
    splice_variant = splice_variant.fillna(False)
    loftee_high_confidence = df.get(
        "LoF", pd.Series("", index=df.index)
    ).fillna("").astype(str).str.contains(
        r"(?:^|[&,;|])(?:HC|OS)(?:[&,;|]|$)", regex=True
    )
    loftee_splice_lof = splice_variant & loftee_high_confidence

    def _coerce_score_column(column_name: str, missing_default: float | None = None) -> pd.Series:
        raw_values = df[column_name]
        numeric_values = pd.to_numeric(raw_values, errors='coerce')
        non_empty_mask = raw_values.notna() & (raw_values.astype(str).str.strip() != "")
        non_numeric_mask = non_empty_mask & numeric_values.isna()
        if non_numeric_mask.any():
            examples = sorted(raw_values[non_numeric_mask].astype(str).str.strip().unique())[:5]
            fill_msg = (
                f"filling default={missing_default}"
                if missing_default is not None else
                "leaving as missing"
            )
            logger.warning(
                f"PP3/BP4 score coercion: column '{column_name}' has {non_numeric_mask.sum()} non-numeric values "
                f"(examples: {examples}); treating as missing and {fill_msg}"
            )
        if missing_default is None:
            return numeric_values
        return numeric_values.fillna(missing_default)

    primateai = _coerce_score_column('PrimateAI', 0)
    am_pathogenicity = _coerce_score_column('am_pathogenicity', 0)
    cadd_phred = _coerce_score_column('CADD_phred')
    cadd_reg_phred = _coerce_score_column('CADD_reg_phred')
    cadd_phred_high = cadd_phred >= 20
    # PriVA heuristic: use CADD < 15 only as computational benign support.
    # Missing CADD remains no evidence, and this is not a universal benign cutoff.
    cadd_bp4_phred_cutoff = 15.0
    cadd_phred_low = cadd_phred.notna() & (cadd_phred < cadd_bp4_phred_cutoff)
    cadd_reg_phred_not_high = ~(cadd_reg_phred >= cadd_bp4_phred_cutoff).fillna(False)
    splice_computational_lof = df['splicing_lof'].fillna(False) | loftee_splice_lof

    # BP4: variant is reported benign
    pp3_criteria = ((primateai > 0.8) & missense_variant) | \
                    (cadd_phred_high & np.logical_not(splice_variant) & np.logical_not(missense_variant) & np.logical_not(five_utr_variant)) | \
                    (df['am_class'].fillna("").str.contains('pathogenic') & missense_variant) | \
                    (df['vep_consq_lof'] & np.logical_not(splice_variant) & np.logical_not(missense_variant) & np.logical_not(five_utr_variant)) | \
                    (((splice_computational_lof | cadd_phred_high) & splice_variant) | (df['5UTR_lof'] & five_utr_variant))
    clinvar_benign = df['CLNSIG'].fillna("").str.contains('enign') & (df['CLNREVSTAT'].map(high_confidence_status, na_action="ignore") >= 2)
    pp3_criteria = pp3_criteria & ~clinvar_benign

    # Per ClinGen SVI guidelines (svi_questions_and_updates_09232021.pdf):
    # PP3 should NOT be applied when PVS1 is assigned to avoid double-counting evidence.
    # Both criteria assess similar aspects of variant impact (LoF prediction), and using
    # them together leads to overestimation of pathogenicity.
    if pvs1_criteria is not None:
        # Extended to block PP3 also when PVS1 is Strong (>=3), not only Very Strong (>=4).
        # PVS1_Strong on NMD-escaping truncating variants and PVS1_Very_Strong both already
        # capture the LOFTEE HC / vep_consq_lof signal that PP3 re-uses, so co-applying PP3
        # at PVS1_Strong still double-counts the same underlying LOF evidence.
        pp3_criteria = pp3_criteria & (pvs1_criteria < 3)
        logger.info(f"PP3 double-counting prevention: blocked {(pvs1_criteria >= 3).sum()} variants with PVS1_Strong or stronger")

    missense_benign = (primateai < 0.8).fillna(True) & (am_pathogenicity < 0.564).fillna(True) & missense_variant
    splice_benign = np.logical_not(
        splice_computational_lof
    ) & splice_variant & cadd_phred_low
    utr_benign = np.logical_not(df['5UTR_lof'].fillna(False)) & five_utr_variant
    other_benign = np.logical_not(df['vep_consq_lof'].fillna(False)) & cadd_phred_low & cadd_reg_phred_not_high & np.logical_not(splice_variant) & np.logical_not(missense_variant) & np.logical_not(five_utr_variant)
    bp4_criteria = missense_benign | splice_benign | utr_benign | other_benign
    clinvar_patho = df['CLNSIG'].fillna("").str.contains('athogenic') & (df['CLNREVSTAT'].map(high_confidence_status, na_action="ignore") >= 2)
    bp4_criteria = bp4_criteria & ~clinvar_patho
    bp4_criteria = bp4_criteria & ~pp3_criteria
    if pvs1_criteria is not None:
        # Strong-or-higher PVS1 reflects an explicit LoF model; do not also emit BP4
        # from generic computational scores for the same variant.
        strong_pvs1 = np.asarray(pvs1_criteria) >= 3
        suppressed_bp4 = int((bp4_criteria & strong_pvs1).sum())
        bp4_criteria = bp4_criteria & ~strong_pvs1
        logger.info(f"BP4 suppression: blocked {suppressed_bp4} variants with PVS1_Strong or stronger")

    pp3_array = np.zeros(len(df), dtype=int)
    bp4_array = np.zeros(len(df), dtype=int)
    pp3_array[pp3_criteria] = 1
    bp4_array[bp4_criteria] = 1
    return pp3_array, bp4_array


def BP7_criteria(df: pd.DataFrame) -> pd.Series:
    # BP7: synonymous variant, no splicing-altering consequence, not conserved.
    synonymous = df['Consequence'] == 'synonymous_variant'
    no_splicing_altering = (df['splicing_len_changing'] == False) | (df['splicing_lof'].isna()) | (df['splicing_lof'] == False) | (df['splicing_len_changing'].isna())
    not_conserved = df['Conservation'] <= 5 # 5 is the cutoff for highly conserved of a GERP++ score.
    bp7_array = np.zeros(len(df), dtype=int)
    bp7_array[synonymous & no_splicing_altering & not_conserved] = 1
    return bp7_array
