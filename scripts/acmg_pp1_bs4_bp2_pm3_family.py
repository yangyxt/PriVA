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

from acmg_variant_mechanism import _variant_condition_masks


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
) -> np.ndarray:
    '''
    BS4: lack of segregation within the submitted family.

    Returns an int array over ``df`` rows: 0 = no BS4, 1 = BS4_Supporting,
    3 = BS4. Only family ``fam_name`` in ``ped_df`` is read; an absent family,
    or one with no affected sample present as a column of ``df``, yields all
    zeros. Where branches of both strengths fire for a row, 3 wins.

    BS4 asks whether a relative's genotype dose contradicts the inherited model
    given their affected status. Five of the six branches answer that from dose
    and inheritance alone. Mechanism enters one branch only, and only to decide
    whether a single contradiction counts as supporting evidence.

    GENOTYPE DOSE
    =============

    Dose is the count of ``1`` alleles in the GT field, the part of the sample
    column before the first colon. Only the first ALT allele is counted, so a
    ``2/2`` call scores 0. A sample is callable when its GT contains no ``.``
    and at least one ``0`` or ``1``, which also excludes ``2/2``.

    Sex from the pedigree, 1 = male and 2 = female, then rewrites dose by
    locus:

        male, chrX or chrY      any ALT is promoted to dose 2, because a
                                hemizygote carries a full dose.
        female, chrY            never callable, so non-informative.
        female, chrX            left diploid, 0/1/2.
        any sex, autosome       left diploid, 0/1/2.

    Unparseable or absent sex leaves a sample judgeable on autosomes only.

    WHICH PAIRS ARE COMPARED
    ========================

    Every affected-affected pair and every affected-unaffected pair is examined
    independently, and one qualifying pair anywhere in the family flags the
    variant. A pair is examined only when both members are callable, the row
    passes the history gate, and the locus is one both members can be judged
    on:

        autosome    always.
        chrX        for a male, requires x_linked_recessive or
                    x_linked_dominant history; for a female, requires
                    x_linked_dominant history only, since a female
                    heterozygote under an X-linked recessive model is an
                    expected carrier.
        chrY        males with y_linked history only; females excluded.

    Mitochondrial variants are out of scope entirely: the locus test admits
    only autosomes, chrX and chrY.

    The history gate requires a Mendelian history and rejects non-monogenic
    (polygenic, digenic, oligogenic), literal non-Mendelian, and
    incomplete-penetrance histories, the last because a healthy carrier of a
    partly penetrant allele contradicts nothing. There is no
    ``HPO_scope_review_required`` gate here, unlike BS1 and BS2.

    THE SIX BRANCHES THAT FIRE
    ==========================

    Affected versus affected, both members carrying this row:

        1. one hom, one WT                                          -> BS4
           No inheritance gate beyond the history and locus gates. No single
           Mendelian model explains two affected relatives at doses 2 and 0.
        2. one hom, one het                             -> BS4_Supporting
           Requires a hom/hemi requirement with no dominant history of any
           kind, so nothing explains the heterozygote.
        3. one het, one WT                                          -> BS4
           Requires dominant-like history, defined below.

    Affected versus unaffected:

        4. affected dose >= 1, unaffected dose >= 2                 -> BS4
           No inheritance gate. Note the sex normalization above makes a
           healthy male hemizygote on chrX dose 2, so he reaches this branch.
        5. affected hom, unaffected het                 -> BS4_Supporting
           Requires dominant-like history.
        6. affected het, unaffected het                 -> BS4_Supporting
           Requires dominant-like history.

    Every other dose combination yields nothing, by omission rather than by an
    explicit rule: affected het against affected het, affected hom against
    unaffected WT, and affected het against unaffected WT are each
    segregation-compatible with some model.

    "Dominant-like history" is read from ``variant_inheritance``: dominant,
    x_linked_dominant, y_linked or mitochondrial, minus any autosomal-recessive
    requirement. Two consequences follow. The mitochondrial member is
    unreachable, because chrM fails the locus test. And x_linked_recessive is
    not excluded, so a gene carrying both x_linked_recessive and
    x_linked_dominant history reaches branches 3, 5 and 6, while branch 2
    treats that same x_linked_recessive history as a biallelic requirement.

    Branch 2 also reads ``variant_inheritance`` for the "any dominant history"
    question. All other inputs to BS4 come from the same column — and from
    ``variant_penetrance``, the pedigree, the locus, and the genotype. BS4
    reads no mechanism column.

    Why mechanism does not enter here: upstream
    ``select_condition_histories_for_variant`` selects histories by variant
    effect and delivers the result through ``variant_inheritance``. By the time
    BS4 runs, the inheritance already reflects which conditions and mechanisms
    are compatible with this variant. Re-reading the mechanism column to
    confirm variant-history compatibility would duplicate that work, and would
    observe a different and narrower fact: what tags exist in
    ``var_plausible_patho_mechs``, rather than what inheritance was actually
    selected for this variant. The segregation question is purely about whether
    the observed genotype dose is consistent with the inherited model, which is
    an inheritance and dose question.

    Required columns
    ----------------
    Only ``variant_inheritance`` and ``variant_penetrance`` raise KeyError if
    absent. ``chrom`` defaults to empty when absent, treating every locus as
    autosomal. ``ped_df`` must carry ``#FamilyID``, ``IndividualID``,
    ``Phenotype`` (2 = affected, 1 = unaffected) and ``Sex``.
    '''
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

    has_recessive_requirement = condition_masks["has_recessive"]
    has_hom_hemi_requirement = (
        has_recessive_requirement | condition_masks["has_x_linked_recessive"]
    )
    has_any_dominant_history = (
        condition_masks["has_dominant"] | condition_masks["has_x_linked_dominant"]
    )
    hom_hemi_only = has_hom_hemi_requirement & ~has_any_dominant_history

    # For branch 3 (affected het vs affected WT), require dominant-like with no
    # biallelic requirement of any kind. X-linked recessive counts as a biallelic
    # requirement here — if the gene has both x_linked_dominant and x_linked_recessive,
    # an affected person with zero copies visible may be explained by the recessive route.
    dominant_compatible_affected_vs_affected = (
        condition_masks["has_dominant_like"] & ~has_hom_hemi_requirement
    )

    # For branches 5 & 6 (affected vs unaffected), require dominant-like with no
    # *autosomal* biallelic requirement. X-linked recessive is NOT excluded here,
    # because the contradiction is in the affected/unaffected status at the same dose,
    # not in missing genotype data. A female with one X-linked copy being affected
    # contradicts the recessive reading, and being unaffected contradicts the dominant
    # reading — no single model explains both.
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
                & hom_hemi_only
            )
            final_criteria = final_criteria | (
                called
                & mendelian_locus
                & het_wt
                & dominant_compatible_affected_vs_affected
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
                     threads: int = 1) -> Tuple[np.ndarray, np.ndarray]:
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
