import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import acmg_pp1_bs4_bp2_pm3_family as family_criteria  # noqa: E402
from acmg_bs2_bp5_observation import BP5_criteria, BS2_criteria  # noqa: E402
from acmg_pm2_bs1_ba1_frequency import BS1_criteria  # noqa: E402
from acmg_variant_mechanism import _variant_condition_masks  # noqa: E402


def test_variant_condition_masks_are_token_exact() -> None:
    frame = pd.DataFrame(
        {
            "variant_inheritance": [
                "dominant;non_mendelian",
                "dominantly_inherited",
                "recessive;polygenic",
            ],
            "variant_penetrance": [
                "high;complete",
                "incomplete",
                "unknown",
            ],
        }
    )

    masks = _variant_condition_masks(frame)

    assert masks["has_dominant"].tolist() == [True, False, False]
    assert masks["has_recessive"].tolist() == [False, False, True]
    assert masks["has_non_mendelian"].tolist() == [True, False, False]
    assert masks["has_non_monogenic"].tolist() == [False, False, True]
    assert masks["has_high_penetrance"].tolist() == [True, False, False]
    assert masks["has_complete_penetrance"].tolist() == [True, False, False]
    assert masks["has_incomplete_penetrance"].tolist() == [False, True, False]


def _bs1_frame() -> pd.DataFrame:
    inheritance = [
        "dominant",
        "dominant",
        "dominant",
        "dominant;polygenic",
        "dominant;recessive",
        "recessive",
    ]
    penetrance = ["unknown", "high", "incomplete", "unknown", "unknown", "unknown"]
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

    # High penetrance is eligible. Incomplete and polygenic histories block.
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
        "unknown",
        "high",
        "complete",
        "unknown",
        "incomplete",
        "unknown",
        "unknown",
        "unknown",
        "unknown",
        "unknown",
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
                "unknown",
                "unknown",
                "incomplete",
                "unknown",
                "high",
            ],
        }
    )

    result = family_criteria.PP1_criteria(
        frame,
        multi_fam_vcf="fixture.vcf",
        multi_fam_ped="fixture.ped",
    )

    assert result.tolist() == [2, 1, 0, 0, 1]


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
                "unknown",
                "high",
                "incomplete",
                "unknown",
                "unknown",
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
                "unknown",
                "unknown",
                "unknown",
                "incomplete",
                "high",
                "unknown",
                "unknown",
                "unknown",
                "unknown",
            ],
            "alt_disease_hets": [True, True, False, True, True, True, False, True, False],
            "alt_disease_homs": [False, False, True, False, False, False, True, False, True],
        }
    )

    result = BP5_criteria(frame, alt_disease_vcf="unused.vcf")

    assert result.tolist() == [1, 0, 1, 0, 1, 0, 0, 0, 1]


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
        "unknown",
        "unknown",
        "unknown",
        "unknown",
        "incomplete",
        "unknown",
        "high",
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
