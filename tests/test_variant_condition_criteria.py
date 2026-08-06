import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import acmg_pp1_bs4_bp2_pm3_family as family_criteria  # noqa: E402
from acmg_bs2_bp5_observation import BP5_criteria, BS2_criteria  # noqa: E402
from acmg_pm2_bs1_ba1_frequency import BS1_criteria  # noqa: E402
from acmg_variant_mechanism import (  # noqa: E402
    _variant_condition_masks,
    _variant_mechanism_masks,
)


def test_variant_condition_masks_are_token_exact() -> None:
    frame = pd.DataFrame(
        {
            "variant_inheritance": [
                "dominant;non_mendelian",
                "dominantly_inherited",
                "recessive;polygenic",
                "x_linked_recessive;x_linked_dominant",
                "x_linked_unspecified",
            ],
            "variant_penetrance": [
                "complete",
                "incomplete",
                "",
                "",
                "",
            ],
        }
    )

    masks = _variant_condition_masks(frame)

    assert masks["has_dominant"].tolist() == [True, False, False, False, False]
    assert masks["has_recessive"].tolist() == [False, False, True, False, False]
    assert masks["has_x_linked_recessive"].tolist() == [False, False, False, True, False]
    assert masks["has_x_linked_dominant"].tolist() == [False, False, False, True, False]
    assert masks["has_x_linked_unspecified"].tolist() == [False, False, False, False, True]
    assert masks["has_non_mendelian"].tolist() == [True, False, False, False, False]
    assert masks["has_non_monogenic"].tolist() == [False, False, True, False, False]
    assert masks["has_complete_penetrance"].tolist() == [True, False, False, False, False]
    assert masks["has_incomplete_penetrance"].tolist() == [False, True, False, False, False]


def test_variant_mechanism_masks_keep_x_linked_profiles_distinct() -> None:
    frame = pd.DataFrame(
        {
            "var_plausible_patho_mechs": [
                "x_linked_recessive_LOF",
                "x_linked_dominant_LOF",
                "x_linked_unspecified_LOF",
            ],
            "variant_effect": ["predicted_LOF_high_confidence"] * 3,
            "variant_lof_score": [1] * 3,
        }
    )

    masks = _variant_mechanism_masks(frame)

    assert masks["has_x_linked_recessive_compatible"].tolist() == [True, False, False]
    assert masks["has_x_linked_dominant_compatible"].tolist() == [False, True, False]
    assert masks["has_x_linked_unspecified_compatible"].tolist() == [False, False, True]
    assert masks["has_x_recessive_lof_history"].tolist() == [True, False, False]
    assert masks["has_x_dominant_lof_history"].tolist() == [False, True, False]


def _x_linked_frequency_frame() -> pd.DataFrame:
    inheritance = [
        "x_linked_recessive",
        "x_linked_recessive",
        "x_linked_dominant",
        "x_linked_dominant",
        "x_linked_unspecified",
        "x_linked_recessive;x_linked_dominant",
    ]
    size = len(inheritance)
    return pd.DataFrame(
        {
            "chrom": ["chrX"] * size,
            "variant_inheritance": inheritance,
            "variant_penetrance": [""] * size,
            "HPO_scope_review_required": [0] * size,
            "gnomAD_joint_AF_max": [0.1] * size,
            "gnomAD_joint_AN_max": [200] * size,
            "gnomAD_nhomalt_max": [0] * size,
            "gnomAD_nhomalt_XX": [0] * size,
            "gnomAD_nhomalt_XY": [11, 0, 0, 11, 11, 0],
            "gnomAD_joint_AF": [0.05] * size,
            "gnomAD_joint_AN": [400] * size,
            "gnomAD_joint_AC_XX": [0, 20, 11, 0, 20, 11],
            "gnomAD_joint_AN_XX": [200] * size,
            "gnomAD_joint_AF_XX": [0.0, 0.1, 0.055, 0.0, 0.1, 0.055],
            "gnomAD_joint_AC_XY": [11, 0, 0, 11, 11, 0],
            "gnomAD_joint_AN_XY": [100] * size,
            "gnomAD_joint_AF_XY": [0.11, 0.0, 0.0, 0.11, 0.11, 0.0],
            "clinvar_patho_gene_max_af": [0.0] * size,
        }
    )


def test_bs1_and_bs2_use_explicit_x_linked_population_strata() -> None:
    frame = _x_linked_frequency_frame()
    zeros = np.zeros(len(frame), dtype=int)

    bs1 = BS1_criteria(frame, pm2_criteria=zeros)
    bs2 = BS2_criteria(frame, pm2_criteria=zeros)

    # XX carrier observations do not count against X-linked recessive disease,
    # while X-linked dominant disease uses the XX population stratum.
    assert bs1.tolist() == [3, 0, 3, 0, 0, 3]
    assert bs2.tolist() == [3, 0, 3, 0, 0, 3]


def _bs1_frame() -> pd.DataFrame:
    inheritance = [
        "dominant",
        "dominant",
        "dominant",
        "dominant;polygenic",
        "dominant;recessive",
        "recessive",
    ]
    penetrance = ["", "complete", "incomplete", "", "", ""]
    homalts = [0, 0, 0, 0, 0, 12]
    size = len(inheritance)
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * size,
            "variant_inheritance": inheritance,
            "variant_penetrance": penetrance,
            "HPO_scope_review_required": [0] * size,
            "gnomAD_joint_AF_max": [0.1] * size,
            "gnomAD_joint_AN_max": [200] * size,
            "gnomAD_nhomalt_max": homalts,
            "gnomAD_nhomalt_XX": homalts,
            "gnomAD_nhomalt_XY": [0] * size,
            "gnomAD_joint_AF": [0.1] * size,
            "gnomAD_joint_AN": [200] * size,
            "gnomAD_joint_AF_XY": [0.0] * size,
            "clinvar_patho_gene_max_af": [0.0] * size,
        }
    )


def test_bs1_uses_selected_inheritance_and_penetrance_categories() -> None:
    frame = _bs1_frame()

    result = BS1_criteria(frame, pm2_criteria=np.zeros(len(frame), dtype=int))

    # Complete penetrance is eligible. Incomplete and polygenic histories block.
    # When dominant and recessive histories both remain, carrier AF cannot
    # contradict the recessive model; homozygous evidence is required.
    assert result.tolist() == [3, 3, 0, 0, 0, 3]


def test_bs2_uses_selected_categories_and_complete_threshold() -> None:
    inheritance = [
        "dominant",
        "dominant",
        "dominant",
        "dominant",
        "dominant",
        "dominant;polygenic",
        "non_mendelian",
        "recessive",
        "recessive",
        "dominant;recessive",
    ]
    penetrance = [
        "",
        "complete",
        "complete",
        "",
        "incomplete",
        "",
        "",
        "",
        "",
        "",
    ]
    allele_frequency = [0.1, 0.1, 0.02, 0.02, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    homalts = [0, 0, 0, 0, 0, 0, 0, 11, 0, 0]
    size = len(inheritance)
    frame = pd.DataFrame(
        {
            "chrom": ["chr1"] * size,
            "variant_inheritance": inheritance,
            "variant_penetrance": penetrance,
            "HPO_scope_review_required": [0] * size,
            "gnomAD_joint_AF": allele_frequency,
            "gnomAD_joint_AN": [200] * size,
            "gnomAD_nhomalt_XX": homalts,
            "gnomAD_nhomalt_XY": [0] * size,
            "gnomAD_joint_AF_XY": [0.0] * size,
        }
    )

    result = BS2_criteria(frame, pm2_criteria=np.zeros(size, dtype=int))

    assert result.tolist() == [3, 3, 3, 0, 0, 0, 0, 3, 0, 0]


def test_pp1_uses_selected_categories(monkeypatch) -> None:
    variant_ids = [f"1:{position}:A-G" for position in range(100, 105)]
    stats = {
        ("1", position, "A", "G"): {
            "Affected_segregated_inds": 1,
            "Affected_homozygous_inds": 0,
            "Unaffected_segregated_inds": 0,
            "Unaffected_nonhomo_inds": 0,
            "Male_affected_segregated_inds": 0,
            "Male_unaffected_segregated_inds": 0,
            "Female_affected_segregated_inds": 0,
            "Female_unaffected_segregated_inds": 0,
        }
        for position in range(100, 105)
    }
    monkeypatch.setattr(
        family_criteria,
        "find_cosegregating_variants",
        lambda *_args, **_kwargs: stats,
    )
    frame = pd.DataFrame(
        {
            "variant_id": variant_ids,
            "chrom": ["chr1"] * 5,
            "variant_inheritance": [
                "recessive",
                "dominant",
                "dominant",
                "dominant;polygenic",
                "dominant",
            ],
            "variant_penetrance": [
                "",
                "",
                "incomplete",
                "",
                "complete",
            ],
        }
    )

    result = family_criteria.PP1_criteria(
        frame,
        multi_fam_vcf="fixture.vcf",
        multi_fam_ped="fixture.ped",
    )

    assert result.tolist() == [2, 1, 0, 0, 1]


def test_pp1_counts_x_recessive_males_and_x_dominant_both_sexes(
    monkeypatch,
) -> None:
    variant_ids = [f"X:{position}:A-G" for position in range(100, 104)]
    stats = {
        ("X", position, "A", "G"): {
            "Affected_segregated_inds": 3,
            "Affected_homozygous_inds": 0,
            "Unaffected_segregated_inds": 0,
            "Unaffected_nonhomo_inds": 0,
            "Male_affected_segregated_inds": 1,
            "Male_unaffected_segregated_inds": 0,
            "Female_affected_segregated_inds": 1,
            "Female_unaffected_segregated_inds": 0,
        }
        for position in range(100, 104)
    }
    monkeypatch.setattr(
        family_criteria,
        "find_cosegregating_variants",
        lambda *_args, **_kwargs: stats,
    )
    frame = pd.DataFrame(
        {
            "variant_id": variant_ids,
            "chrom": ["chrX"] * 4,
            "variant_inheritance": [
                "x_linked_recessive",
                "x_linked_dominant",
                "x_linked_unspecified",
                "x_linked_recessive;x_linked_dominant",
            ],
            "variant_penetrance": [""] * 4,
        }
    )

    result = family_criteria.PP1_criteria(
        frame,
        multi_fam_vcf="fixture.vcf",
        multi_fam_ped="fixture.ped",
    )

    assert result.tolist() == [1, 2, 0, 2]


def test_bs4_uses_selected_categories() -> None:
    frame = pd.DataFrame(
        {
            "chrom": ["chr1"] * 5,
            "variant_inheritance": [
                "dominant",
                "dominant",
                "dominant",
                "dominant;polygenic",
                "recessive",
            ],
            "variant_penetrance": [
                "",
                "complete",
                "incomplete",
                "",
                "",
            ],
            "var_plausible_patho_mechs": [
                "dominant",
                "dominant",
                "dominant",
                "dominant",
                "recessive",
            ],
            "variant_effect": ["uncertain"] * 5,
            "PROBAND": ["0/1:20"] * 5,
            "CONTROL": ["0/1:20"] * 5,
        }
    )
    pedigree = pd.DataFrame(
        {
            "#FamilyID": ["F1", "F1"],
            "IndividualID": ["PROBAND", "CONTROL"],
            "Phenotype": ["2", "1"],
            "Sex": ["2", "2"],
        }
    )

    result = family_criteria.BS4_criteria(frame, pedigree, "F1")

    assert result.tolist() == [1, 1, 0, 0, 0]


def test_bs4_does_not_treat_female_x_recessive_carriers_as_dominant() -> None:
    frame = pd.DataFrame(
        {
            "chrom": ["chrX"] * 3,
            "variant_inheritance": [
                "x_linked_recessive",
                "x_linked_dominant",
                "x_linked_recessive;x_linked_dominant",
            ],
            "variant_penetrance": [""] * 3,
            "var_plausible_patho_mechs": [
                "x_linked_recessive_LOF",
                "x_linked_dominant_LOF",
                "x_linked_recessive_LOF;x_linked_dominant_LOF",
            ],
            "variant_effect": ["predicted_LOF_high_confidence"] * 3,
            "variant_lof_score": [1] * 3,
            "PROBAND": ["0/1:20"] * 3,
            "CONTROL": ["0/1:20"] * 3,
        }
    )
    pedigree = pd.DataFrame(
        {
            "#FamilyID": ["F1", "F1"],
            "IndividualID": ["PROBAND", "CONTROL"],
            "Phenotype": ["2", "1"],
            "Sex": ["2", "2"],
        }
    )

    result = family_criteria.BS4_criteria(frame, pedigree, "F1")

    assert result.tolist() == [0, 1, 1]


def test_bs4_uses_only_sex_compatible_explicit_x_models() -> None:
    frame = pd.DataFrame(
        {
            "chrom": ["chrX"] * 4,
            "variant_inheritance": [
                "x_linked_recessive",
                "x_linked_dominant",
                "x_linked_unspecified",
                "x_linked_recessive;x_linked_dominant",
            ],
            "variant_penetrance": [""] * 4,
            "var_plausible_patho_mechs": [
                "x_linked_recessive_LOF",
                "x_linked_dominant_LOF",
                "x_linked_unspecified_LOF",
                "x_linked_recessive_LOF;x_linked_dominant_LOF",
            ],
            "variant_effect": ["predicted_LOF_high_confidence"] * 4,
            "variant_lof_score": [1] * 4,
            "PROBAND": ["0/1:20"] * 4,
            "CONTROL": ["0/0:20"] * 4,
        }
    )

    def pedigree(sex: str) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "#FamilyID": ["F1", "F1"],
                "IndividualID": ["PROBAND", "CONTROL"],
                "Phenotype": ["2", "2"],
                "Sex": [sex, sex],
            }
        )

    assert family_criteria.BS4_criteria(frame, pedigree("1"), "F1").tolist() == [
        3,
        3,
        0,
        3,
    ]
    assert family_criteria.BS4_criteria(frame, pedigree("3"), "F1").tolist() == [
        0,
        0,
        0,
        0,
    ]


def test_bp5_uses_selected_categories_without_rescanning_vcf() -> None:
    frame = pd.DataFrame(
        {
            "variant_inheritance": [
                "dominant",
                "recessive",
                "recessive",
                "dominant",
                "dominant",
                "dominant;polygenic",
                "non_mendelian",
                "dominant;recessive",
                "dominant;recessive",
            ],
            "variant_penetrance": [
                "",
                "",
                "",
                "incomplete",
                "complete",
                "",
                "",
                "",
                "",
            ],
            "alt_disease_hets": [True, True, False, True, True, True, False, True, False],
            "alt_disease_homs": [False, False, True, False, False, False, True, False, True],
        }
    )

    result = BP5_criteria(frame, alt_disease_vcf="unused.vcf")

    assert result.tolist() == [1, 0, 1, 0, 1, 0, 0, 0, 1]


def test_bp5_keeps_x_recessive_and_x_dominant_models_separate() -> None:
    frame = pd.DataFrame(
        {
            "variant_inheritance": [
                "x_linked_recessive",
                "x_linked_recessive",
                "x_linked_dominant",
                "x_linked_unspecified",
                "x_linked_recessive;x_linked_dominant",
            ],
            "variant_penetrance": [""] * 5,
            "alt_disease_hets": [True, False, True, False, True],
            "alt_disease_homs": [False, True, False, True, False],
        }
    )

    result = BP5_criteria(frame, alt_disease_vcf="unused.vcf")

    assert result.tolist() == [0, 1, 1, 0, 1]


def test_bp2_pm3_use_selected_categories(monkeypatch) -> None:
    inheritance = [
        "dominant",
        "recessive",
        "recessive",
        "dominant;recessive",
        "dominant",
        "dominant;polygenic",
        "dominant",
    ]
    penetrance = [
        "",
        "",
        "",
        "",
        "incomplete",
        "",
        "complete",
    ]
    in_cis = np.array([0, 0, 1, 0, 0, 0, 0])
    in_trans = np.array([1, 1, 0, 1, 1, 1, 1])

    def fake_phase(frame, _pathogenic, _pedigree, _threads):
        return in_cis, in_trans, frame

    monkeypatch.setattr(
        family_criteria,
        "determine_cis_trans_relationships",
        fake_phase,
    )
    frame = pd.DataFrame(
        {
            "variant_inheritance": inheritance,
            "variant_penetrance": penetrance,
            "vep_consq_lof": [False] * 7,
            "splicing_lof": [False] * 7,
            "5UTR_lof": [False] * 7,
        }
    )
    zeros = np.zeros(7, dtype=int)
    pm2 = np.array([0, 1, 0, 1, 0, 0, 0])

    bp2, pm3 = family_criteria.BP2_PM3_criteria(
        frame,
        pd.DataFrame(),
        pm2,
        zeros,
        zeros,
        zeros,
    )

    assert bp2.tolist() == [1, 0, 1, 0, 0, 0, 1]
    assert pm3.tolist() == [0, 2, 0, 2, 0, 0, 0]


def test_bp2_pm3_do_not_treat_x_recessive_as_autosomal_biallelic(monkeypatch) -> None:
    inheritance = [
        "x_linked_recessive",
        "x_linked_recessive",
        "x_linked_dominant",
        "x_linked_unspecified",
        "x_linked_recessive;x_linked_dominant",
        "recessive",
    ]
    in_cis = np.array([0, 1, 0, 0, 0, 0])
    in_trans = np.array([1, 0, 1, 1, 1, 1])

    def fake_phase(frame, _pathogenic, _pedigree, _threads):
        return in_cis, in_trans, frame

    monkeypatch.setattr(
        family_criteria,
        "determine_cis_trans_relationships",
        fake_phase,
    )
    frame = pd.DataFrame(
        {
            "variant_inheritance": inheritance,
            "variant_penetrance": [""] * len(inheritance),
            "vep_consq_lof": [False] * len(inheritance),
            "splicing_lof": [False] * len(inheritance),
            "5UTR_lof": [False] * len(inheritance),
        }
    )
    zeros = np.zeros(len(frame), dtype=int)
    pm2 = np.ones(len(frame), dtype=int)

    bp2, pm3 = family_criteria.BP2_PM3_criteria(
        frame,
        pd.DataFrame(),
        pm2,
        zeros,
        zeros,
        zeros,
    )

    assert bp2.tolist() == [0, 0, 1, 0, 0, 0]
    assert pm3.tolist() == [0, 0, 0, 0, 0, 2]
