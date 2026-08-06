from pathlib import Path
import sys

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from acmg_scoring import sort_and_rank_variants  # noqa: E402


@pytest.mark.parametrize(
    ("sex", "mechanism_profile", "expected_compatibility"),
    [
        ("1", "x_linked_recessive_LOF", 1.0),
        ("2", "x_linked_recessive_LOF", 0.9),
        ("1", "x_linked_dominant_LOF", 1.0),
        ("2", "x_linked_dominant_LOF", 1.0),
        ("3", "x_linked_recessive_LOF", 0.9),
        ("3", "x_linked_dominant_LOF", 0.9),
        ("1", "x_linked_unspecified_LOF", 0.9),
    ],
)
def test_x_linked_ranking_requires_an_explicit_compatible_proband_sex(
    sex: str,
    mechanism_profile: str,
    expected_compatibility: float,
) -> None:
    variants = pd.DataFrame(
        {
            "ACMG_quant_score": [0.95],
            "BIOTYPE": ["protein_coding"],
            "PROBAND": ["0/1:20"],
            "var_plausible_patho_mechs": [mechanism_profile],
            "variant_effect": ["predicted_LOF_high_confidence"],
            "variant_lof_score": [2],
            "chrom": ["chrX"],
            "pos": [100],
            "ref": ["A"],
            "alt": ["G"],
            "variant_id": ["X:100:A-G"],
            "SYMBOL": ["TESTX"],
        }
    )
    pedigree = pd.DataFrame(
        {
            "#FamilyID": ["F1"],
            "IndividualID": ["PROBAND"],
            "Phenotype": ["2"],
            "Sex": [sex],
        }
    )

    ranked = sort_and_rank_variants(
        variants,
        pedigree,
        "F1",
        dispensable_gene_list=None,
    )

    assert ranked.iloc[0][
        "zygosity_inheritance_mechanism_compatibility"
    ] == expected_compatibility
