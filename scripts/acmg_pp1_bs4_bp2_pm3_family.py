#!/usr/bin/env python3
"""PP1, BS4, BP2 and PM3 -- evidence from the genotypes of relatives.

All four need somebody other than the proband genotyped, which is why they
sit together.

    PP1  the variant co-segregates with disease across affected relatives
    BS4  the variant fails to segregate -- an affected relative lacks it, or
         an unaffected relative carries it in a pattern the disease model
         forbids
    PM3  seen in trans with a pathogenic variant in an autosomal recessive disorder,
         which is what a second hit looks like
    BP2  seen in trans with a pathogenic variant in a dominant disorder, or
         in cis under an autosomal recessive model; single-copy X phase is not
         treated as biallelic evidence

All four read the mechanism hub's selected condition histories and are gated
off for non-monogenic, non-Mendelian and incomplete-equivalent penetrance,
where a genotype contradiction is not interpretable.
"""

import logging
import pandas as pd
import numpy as np
from typing import Tuple
from scipy import stats

from determine_phase import determine_cis_trans_relationships
from find_cosegregation_vars import find_cosegregating_variants

from acmg_variant_mechanism import (
    _variant_condition_masks,
    _variant_mechanism_masks,
)


logger = logging.getLogger(__name__)


def PP1_criteria(df: pd.DataFrame,
                 multi_fam_vcf: str = "",
                 multi_fam_ped: str = "",
                 mode: str = "both",) -> np.ndarray:
    '''
    PP1: The variant is cosegregating with a pathogenic variant in one or more families

    THE DECISION TREE
    =================

    PP1 reads the inheritance and penetrance of the condition histories the
    hub selected for THIS variant. Compatible curated mechanisms contribute at
    high resolution; mechanism-free included conditions contribute unresolved
    gene history instead of being discarded or guessed.

        variant_inheritance / variant_penetrance
              |
              v
        Every state first passes the non-monogenic, non-Mendelian, and
        incomplete-penetrance gates. Evidence is then counted by affected
        allele state:

        recessive                  autosomal segregation counts
        dominant                   autosomal segregation counts
        x_linked_recessive         male segregation counts
        x_linked_dominant          male and female segregation counts
        x_linked_unspecified       no PP1 genotype model
              |
              v
        segregation counts, per family, from find_cosegregating_variants:
              Affected_segregated_inds     Unaffected_segregated_inds
              Male_affected_segregated_inds
              Male_unaffected_segregated_inds     <- used for chrX and chrY
              Female_affected_segregated_inds
              Female_unaffected_segregated_inds   <- used for chrX
              |
              v
        autosomal = NOT chrX & NOT chrY & NOT chrM
              |
              +-- each explicit state -> non-overlapping points -> strength

    Explicit incomplete, low/moderate penetrance, and late onset have already
    been linked to their condition and normalized to ``incomplete`` upstream;
    high penetrance is normalized to ``complete``.
    '''
    pp1_array = np.zeros(len(df), dtype=int)
    if multi_fam_vcf and multi_fam_ped:
        cosegregating_variants = find_cosegregating_variants(multi_fam_vcf, multi_fam_ped, mode)
        affected_segs = {}
        affected_homs = {}
        control_segs = {}
        control_nonhoms = {}
        male_affected_segs = {}
        male_control_segs = {}
        female_affected_segs = {}
        female_control_segs = {}
        for variant, stats in cosegregating_variants.items():
            varid_key = str(variant[0]) + ":" + str(variant[1]) + ":" + str(variant[2]) + "-" + str(variant[3])
            affected_segs[varid_key] = stats["Affected_segregated_inds"]
            affected_homs[varid_key] = stats["Affected_homozygous_inds"]
            control_segs[varid_key] = stats["Unaffected_segregated_inds"]
            control_nonhoms[varid_key] = stats["Unaffected_nonhomo_inds"]
            male_affected_segs[varid_key] = stats["Male_affected_segregated_inds"]
            male_control_segs[varid_key] = stats["Male_unaffected_segregated_inds"]
            female_affected_segs[varid_key] = stats[
                "Female_affected_segregated_inds"
            ]
            female_control_segs[varid_key] = stats[
                "Female_unaffected_segregated_inds"
            ]

        affected_seg_count = df['variant_id'].map(affected_segs).fillna(0).astype(int)
        affected_hom_count = df['variant_id'].map(affected_homs).fillna(0).astype(int)
        control_seg_count = df['variant_id'].map(control_segs).fillna(0).astype(int)
        control_nonhom_count = df['variant_id'].map(control_nonhoms).fillna(0).astype(int)
        male_affected_seg_count = df['variant_id'].map(male_affected_segs).fillna(0).astype(int)
        male_control_seg_count = df['variant_id'].map(male_control_segs).fillna(0).astype(int)
        female_affected_seg_count = df['variant_id'].map(female_affected_segs).fillna(0).astype(int)
        female_control_seg_count = df['variant_id'].map(female_control_segs).fillna(0).astype(int)

        condition_masks = _variant_condition_masks(df)
        valid_model = (
            condition_masks["has_mendelian"]
            & ~condition_masks["has_non_monogenic"]
            & ~condition_masks["has_non_mendelian"]
            & ~condition_masks["has_incomplete_penetrance"]
        )
        recessive_ih = valid_model & condition_masks["has_recessive"]
        x_recessive_ih = valid_model & condition_masks["has_x_linked_recessive"]
        x_dominant_ih = (
            valid_model
            & condition_masks["has_x_linked_dominant"]
            & ~condition_masks["has_recessive"]
        )
        dominant_ih = (
            valid_model
            & condition_masks["has_dominant"]
            & ~condition_masks["has_recessive"]
        )

        pp1_points = [0] * len(df)
        recessive_encode = recessive_ih.astype(int)
        x_recessive_encode = x_recessive_ih.astype(int)
        x_dominant_encode = x_dominant_ih.astype(int)
        dominant_encode = dominant_ih.astype(int)

        autosomal = np.logical_not(df["chrom"].str.contains("X") | df["chrom"].str.contains("Y") | df['chrom'].str.contains("M"))
        x_chr = df["chrom"].str.contains("X")
        autosomal_encode = autosomal.astype(int)
        x_encode = x_chr.astype(int)

        pp1_recessive_points = (
            autosomal_encode * recessive_encode * affected_seg_count * 2
            + autosomal_encode * recessive_encode * control_seg_count * 0.4
        )
        pp1_x_recessive_points = x_encode * x_recessive_encode * (
            male_affected_seg_count + male_control_seg_count
        )

        pp1_x_dominant_points = (
            x_encode
            * x_dominant_encode
            * (
                male_affected_seg_count
                + male_control_seg_count
                + female_affected_seg_count
                + female_control_seg_count
            )
        )

        pp1_dominant_points = autosomal_encode * dominant_encode * (
            affected_seg_count
            + control_seg_count
            + affected_hom_count
            + control_nonhom_count
        )

        pp1_points = (
            pp1_recessive_points
            + pp1_x_recessive_points
            + pp1_x_dominant_points
            + pp1_dominant_points
        )

        pp1_array[(pp1_points <= 1.9) & (pp1_points >= 1)] = 1  # Supporting
        pp1_array[(pp1_points >= 2) & (pp1_points <= 3.9)] = 2  # Moderate
        pp1_array[(pp1_points >= 4)] = 3   # Strong

        return pp1_array

    return pp1_array


def BS4_criteria(
    df: pd.DataFrame,
    ped_df: pd.DataFrame,
    fam_name: str,
) -> pd.Series:
    # BS4: lack of segregation / genotype incompatibility within the submitted
    # family. Output keeps the existing PriVA strength encoding:
    #   0 = no BS4, 1 = BS4_Supporting, 3 = BS4.
    #
    # The hub has already limited known histories to mechanisms compatible with
    # this variant and retained every mechanism-unresolved included condition.
    # Do not assign BS4 for non-monogenic/polygenic, non-Mendelian, or any
    # selected incomplete-equivalent penetrance history, because genotype
    # contradictions are not interpretable there.
    #
    # Biallelic, X-linked recessive, autosomal dominant, and X-linked dominant
    # compatibility come only from the upstream compact
    # variant-level assertions.
    #
    # Inputs per variant:
    #   patient_GTs = all affected family members with callable GT.
    #   control_GTs = all unaffected family members with callable GT.
    #   has_recessive / has_x_linked_recessive identify biallelic/XY models.
    #   has_dominant / has_x_linked_dominant identify dominant models.
    #   is_predicted_LOF / is_exact_GOF come from variant_effect.
    #
    # Variant compatibility under dominant-only history:
    #   dominant only           -> any variant state remains uncertain-compatible.
    #   dominant_LOF only       -> predicted LOF or uncertain; exact GOF excluded.
    #   dominant_GOF only       -> exact_GOF only.
    #   dominant_DN only        -> ambiguous only; PriVA has no exact DN DB.
    #   dominant_GOF + DN       -> exact_GOF or ambiguous; NMD_LOF excluded.
    #   dominant_LOF + GOF/DN   -> any variant state can be interpreted under
    #                              at least one dominant history.
    #
    # Sex-chromosome variants use the same tree after sex-aware allele-state
    # normalization: male chrX/chrY with any ALT allele is treated as hom/hemi
    # (2), female chrX remains diploid (0/1/2), female chrY is non-informative,
    # and non-sex chromosomes remain unchanged.
    #
    # Patient-patient comparison:
    #   1. hom patient vs WT patient -> BS4.
    #   2. hom patient vs het patient:
    #      - recessive only -> BS4_Supporting.
    #      - recessive + dominant_GOF/DN + predicted LOF, with no dominant_LOF
    #        or unresolved dominant history -> BS4_Supporting. The heterozygous
    #        affected genotype is incompatible with AR, and the NMD_LOF query
    #        variant is not compatible with GOF/DN.
    #      - recessive + dominant_GOF/DN + uncertain -> no BS4; uncertainty
    #        cannot rule out a dominant DN/GOF-compatible mechanism.
    #      - recessive + dominant_LOF or unresolved dominant -> no BS4.
    #      - dominant-only history -> no BS4, because both affected individuals
    #        carry ALT and dominant disease only requires one ALT allele.
    #   3. het patient vs WT patient:
    #      - no recessive history + dominant-compatible variant -> BS4.
    #      - recessive history present -> no BS4, because the heterozygous affected
    #        genotype may reflect a second allele not represented by this row or
    #        a dominant-compatible history.
    #   4. het patient vs het patient -> no BS4.
    #
    # Patient-control comparison:
    #   5. hom patient vs WT control -> no BS4.
    #   6. hom patient vs het control:
    #      - recessive history present -> no BS4; a healthy heterozygous carrier
    #        is segregation-compatible for a recessive LoF model.
    #      - no recessive history + dominant-compatible variant -> BS4_Supporting.
    #   7. carrier patient vs hom/hemi control -> BS4.
    #   8. het patient vs WT control -> no BS4.
    #      Unaffected WT relatives are segregation-compatible for a dominant
    #      heterozygous disease model and are not contradictory for AR logic.
    #   9. het patient vs het control:
    #      - recessive history present -> no BS4, even when dominant history also
    #        exists, because a healthy heterozygous carrier can fit the AR model.
    #      - no recessive history + dominant-compatible variant -> BS4_Supporting.
    fam_ped_df = ped_df.loc[ped_df['#FamilyID'] == fam_name, :].copy()
    if fam_ped_df.empty:
        bs4_array = np.zeros(len(df), dtype=int)
        return bs4_array
    fam_ped_df.loc[:, "Phenotype"] = pd.to_numeric(fam_ped_df["Phenotype"], errors="coerce").astype("Int64")
    fam_ped_df.loc[:, "IndividualID"] = fam_ped_df["IndividualID"].astype(str)
    fam_ped_df.loc[:, "Sex"] = pd.to_numeric(fam_ped_df.get("Sex", pd.Series(index=fam_ped_df.index)), errors="coerce")
    sex_by_sample = dict(zip(fam_ped_df["IndividualID"], fam_ped_df["Sex"]))

    patient_cols = [
        sample for sample in fam_ped_df.loc[fam_ped_df["Phenotype"] == 2, "IndividualID"].tolist()
        if sample in df.columns
    ]
    control_cols = [
        sample for sample in fam_ped_df.loc[fam_ped_df["Phenotype"] == 1, "IndividualID"].tolist()
        if sample in df.columns
    ]

    final_criteria = pd.Series(False, index=df.index)
    supporting_criteria = pd.Series(False, index=df.index)
    if not patient_cols:
        bs4_array = np.zeros(len(df), dtype=int)
        return bs4_array

    def _gt_alt_count(sample_col: str) -> pd.Series:
        return df[sample_col].fillna(".").astype(str).str.split(":", n=1).str[0].str.count("1")

    def _effective_alt_state(sample_col: str, raw_alt_count: pd.Series) -> pd.Series:
        if sex_by_sample.get(sample_col) == 1:
            return raw_alt_count.mask((x_linked | y_linked) & (raw_alt_count >= 1), 2)
        return raw_alt_count

    def _gt_called(sample_col: str) -> pd.Series:
        gt = df[sample_col].fillna(".").astype(str).str.split(":", n=1).str[0]
        called = (~gt.str.contains(".", regex=False)) & gt.str.contains(r"[01]", regex=True)
        if sex_by_sample.get(sample_col) == 2:
            called = called & ~y_linked
        return called

    def _series_or_empty(column: str) -> pd.Series:
        return df.get(column, pd.Series("", index=df.index)).fillna("").astype(str)

    condition_masks = _variant_condition_masks(df)
    valid_model = (
        condition_masks["has_mendelian"]
        & ~condition_masks["has_non_monogenic"]
        & ~condition_masks["has_non_mendelian"]
        & ~condition_masks["has_incomplete_penetrance"]
    )

    chrom = _series_or_empty("chrom")
    autosomal = ~chrom.str.contains("X|Y|M", regex=True)
    x_linked = chrom.str.contains("X", regex=False)
    y_linked = chrom.str.contains("Y", regex=False)
    mendelian_locus = autosomal | x_linked | y_linked

    def _sample_has_applicable_model(sample_col: str) -> pd.Series:
        sex = sex_by_sample.get(sample_col)
        x_model = pd.Series(False, index=df.index)
        y_model = pd.Series(False, index=df.index)
        if sex == 1:
            x_model = (
                condition_masks["has_x_linked_recessive"]
                | condition_masks["has_x_linked_dominant"]
            )
            y_model = condition_masks["has_y_linked"]
        elif sex == 2:
            x_model = condition_masks["has_x_linked_dominant"]
        return autosomal | (x_linked & x_model) | (y_linked & y_model)

    mechanism_masks = _variant_mechanism_masks(df)

    has_recessive_requirement = condition_masks["has_recessive"]
    has_hom_hemi_requirement = (
        has_recessive_requirement | condition_masks["has_x_linked_recessive"]
    )
    has_any_dominant_history = (
        mechanism_masks["has_dom_lof_history"]
        | mechanism_masks["has_dom_gof_history"]
        | mechanism_masks["has_dom_dn_history"]
        | mechanism_masks["has_dom_unresolved_history"]
        | mechanism_masks["has_x_dominant_lof_history"]
        | mechanism_masks["has_x_dominant_gof_history"]
        | mechanism_masks["has_x_dominant_dn_history"]
        | mechanism_masks["has_x_dominant_unresolved_history"]
    )
    hom_hemi_only = has_hom_hemi_requirement & ~has_any_dominant_history
    dominant_gof_dn_history = (
        mechanism_masks["has_dom_gof_history"]
        | mechanism_masks["has_dom_dn_history"]
        | mechanism_masks["has_x_dominant_gof_history"]
        | mechanism_masks["has_x_dominant_dn_history"]
    )
    dominant_lof_or_unresolved_history = (
        mechanism_masks["has_dom_lof_history"]
        | mechanism_masks["has_dom_unresolved_history"]
        | mechanism_masks["has_x_dominant_lof_history"]
        | mechanism_masks["has_x_dominant_unresolved_history"]
    )
    hom_hemi_plus_gof_dn_nmd_only = (
        has_hom_hemi_requirement
        & dominant_gof_dn_history
        & ~dominant_lof_or_unresolved_history
    )
    is_nmd_lof = mechanism_masks["is_predicted_lof"]
    dominant_variant_compatible_no_biallelic = (
        condition_masks["has_dominant_like"] & ~has_recessive_requirement
    )

    patient_alt_counts = {
        sample: _effective_alt_state(sample, _gt_alt_count(sample))
        for sample in patient_cols
    }
    patient_called = {sample: _gt_called(sample) for sample in patient_cols}
    control_alt_counts = {
        sample: _effective_alt_state(sample, _gt_alt_count(sample))
        for sample in control_cols
    }
    control_called = {sample: _gt_called(sample) for sample in control_cols}

    # Patient-patient comparisons.
    for i, sample_a in enumerate(patient_cols):
        for sample_b in patient_cols[i + 1:]:
            a_alt = patient_alt_counts[sample_a]
            b_alt = patient_alt_counts[sample_b]
            called = (
                patient_called[sample_a]
                & patient_called[sample_b]
                & valid_model
                & _sample_has_applicable_model(sample_a)
                & _sample_has_applicable_model(sample_b)
            )

            hom_wt = (
                ((a_alt >= 2) & (b_alt == 0))
                | ((b_alt >= 2) & (a_alt == 0))
            )
            hom_het = (
                ((a_alt >= 2) & (b_alt == 1))
                | ((b_alt >= 2) & (a_alt == 1))
            )
            het_wt = (
                ((a_alt == 1) & (b_alt == 0))
                | ((b_alt == 1) & (a_alt == 0))
            )

            final_criteria = final_criteria | (called & mendelian_locus & hom_wt)
            supporting_criteria = supporting_criteria | (
                called
                & mendelian_locus
                & hom_het
                & (
                    hom_hemi_only
                    | (hom_hemi_plus_gof_dn_nmd_only & is_nmd_lof)
                )
            )
            final_criteria = final_criteria | (
                called
                & mendelian_locus
                & het_wt
                & dominant_variant_compatible_no_biallelic
            )

    # Patient-control comparisons.
    for patient in patient_cols:
        p_alt = patient_alt_counts[patient]
        p_called = patient_called[patient]
        for control in control_cols:
            c_alt = control_alt_counts[control]
            called = (
                p_called
                & control_called[control]
                & valid_model
                & _sample_has_applicable_model(patient)
                & _sample_has_applicable_model(control)
            )

            patient_hom_control_het = (p_alt >= 2) & (c_alt == 1)
            patient_het_control_het = (p_alt == 1) & (c_alt == 1)
            patient_carrier_control_hom = (p_alt >= 1) & (c_alt >= 2)

            final_criteria = final_criteria | (called & mendelian_locus & patient_carrier_control_hom)
            supporting_criteria = supporting_criteria | (
                called
                & mendelian_locus
                & patient_hom_control_het
                & dominant_variant_compatible_no_biallelic
            )
            supporting_criteria = supporting_criteria | (
                called
                & mendelian_locus
                & patient_het_control_het
                & dominant_variant_compatible_no_biallelic
            )

    bs4_array = np.zeros(len(df), dtype=int)
    bs4_array[supporting_criteria.to_numpy()] = 1
    bs4_array[final_criteria.to_numpy()] = 3
    return bs4_array


def BP2_PM3_criteria(df: pd.DataFrame,
                     ped_df: pd.DataFrame,
                     pm2_criteria: np.ndarray,
                     pvs1_criteria: np.ndarray,
                     ps1_criteria: np.ndarray,
                     ps3_criteria: np.ndarray,
                     threads: int = 1) -> Tuple[pd.Series, pd.Series]:
    # BP2: observed in trans with a pathogenic variant in a dominant disease, or
    # in cis with a pathogenic variant in an autosomal recessive disease.
    # PM3: observed in trans with a pathogenic variant in recessive disease.
    #
    # THE DECISION TREE
    # =================
    #
    # Two criteria of opposite direction share one function because they read
    # the same phase evidence and split on inheritance alone. Both take the
    # condition-linked inheritance and penetrance selected upstream for this
    # variant. Known incompatible mechanisms have already been removed, while
    # included conditions with no explicit mechanism remain as unresolved
    # histories.
    #
    #   what counts as a pathogenic partner variant:
    #       vep_consq_lof | splicing_lof | 5UTR_lof
    #                     | (PS1 AND PS3) | PVS1
    #             |
    #             v
    #   determine_cis_trans_relationships  ->  in_cis / in_trans
    #             |
    #             +-- BP2, benign supporting
    #             |     in_trans AND (autosomal dominant OR
    #             |                           x_linked_dominant)
    #             |                    AND NOT autosomal recessive
    #             |         a second pathogenic allele on the other copy is
    #             |         incompatible with dominant disease
    #             |     OR in_cis AND autosomal recessive
    #             |         both hits on one copy leave the other intact, so
    #             |         a recessive disease is not explained
    #             |     x_linked_recessive is excluded from both branches: a person
    #             |         with one X copy has no second homolog on which phase
    #             |         can provide BP2 evidence
    #             |
    #             +-- PM3, pathogenic moderate
    #                   in_trans AND is_recessive AND PM2
    #                       the partner completes a biallelic genotype, and
    #                       PM2 requires the variant to be rare enough for
    #                       that to be meaningful
    #                   x_linked_recessive is deliberately excluded: a male with one
    #                       X copy does not need a second allele. Male X variants
    #                       are also normalized to in-cis by the phase module.
    #
    # Note the asymmetry: dominant is additionally guarded by NOT recessive.
    # This remains necessary when the selected evidence contains unresolved
    # histories under both models: the improved resolution is used wherever a
    # mechanism is known, but unknown mechanisms must retain the gene's full
    # included disease history.
    #
    # Any selected ``incomplete`` value gates both criteria because phase cannot
    # contradict a model when carriers may be unaffected. The upstream category
    # also covers low/moderate penetrance and late onset; other onset terms,
    # age-dependent penetrance and variable expressivity are non-blocking; high
    # penetrance has already normalized to complete.
    pathogenic = df["vep_consq_lof"] | df["splicing_lof"] | df["5UTR_lof"] | (ps1_criteria & ps3_criteria) | pvs1_criteria

    in_cis_pathogenic, in_trans_pathogenic, df = determine_cis_trans_relationships( df,
                                                                                    pathogenic,
                                                                                    ped_df,
                                                                                    threads )

    in_cis_pathogenic = in_cis_pathogenic > 0
    in_trans_pathogenic = in_trans_pathogenic > 0

    condition_masks = _variant_condition_masks(df)
    is_recessive = condition_masks["has_recessive"]
    is_x_linked_recessive = condition_masks["has_x_linked_recessive"]
    is_dominant = condition_masks["has_dominant_like"]
    valid_model = (
        condition_masks["has_mendelian"]
        & ~condition_masks["has_non_monogenic"]
        & ~condition_masks["has_non_mendelian"]
        & ~condition_masks["has_incomplete_penetrance"]
    )

    bp2_criteria = (
        (
            in_trans_pathogenic
            & is_dominant
            & np.logical_not(is_recessive | is_x_linked_recessive)
        )
        | (in_cis_pathogenic & is_recessive)
    ) & valid_model
    pm3_criteria = (
        in_trans_pathogenic
        & is_recessive
        & (np.asarray(pm2_criteria) > 0)
        & valid_model
    )

    bp2_array = np.zeros(len(df), dtype=int)
    pm3_array = np.zeros(len(df), dtype=int)
    bp2_array[bp2_criteria] = 1
    pm3_array[pm3_criteria] = 2

    return bp2_array, pm3_array
