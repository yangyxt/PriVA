from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gene_mechanism_hub import (  # noqa: E402
    build_hpo_condition_index,
    enrich_condition_mechanism_assertion,
    select_condition_histories_for_variant,
)


def test_build_hpo_condition_index_preserves_assertion_provenance() -> None:
    assertions = pd.DataFrame(
        [
            {
                "gene_symbol": "GENE1",
                "disease_id": "OMIM:123",
                "mondo_id": "MONDO:0000123",
                "hpo_id": "HP:0000006",
                "frequency": "-",
                "evidence": "TAS",
                "reference": "OMIM:123",
                "disease_scope": "mendelian_non_neoplastic",
                "priva_scope": "include",
                "scope_review_status": "auto_supported",
            },
            {
                "gene_symbol": "GENE1",
                "disease_id": "OMIM:123",
                "mondo_id": "MONDO:0000123",
                "hpo_id": "HP:0003829",
                "frequency": "2/10",
                "evidence": "PCS",
                "reference": "PMID:111",
                "disease_scope": "mendelian_non_neoplastic",
                "priva_scope": "include",
                "scope_review_status": "auto_supported",
            },
            {
                "gene_symbol": "GENE1",
                "disease_id": "OMIM:123",
                "mondo_id": "MONDO:0000123",
                "hpo_id": "HP:0003584",
                "frequency": "1/10",
                "evidence": "PCS",
                "reference": "PMID:222",
                "disease_scope": "mendelian_non_neoplastic",
                "priva_scope": "include",
                "scope_review_status": "auto_supported",
            },
        ]
    )

    index = build_hpo_condition_index(assertions)
    by_omim = index[("GENE1", "OMIM:123")]
    by_mondo = index[("GENE1", "MONDO:0000123")]

    assert by_mondo is by_omim
    assert by_omim["inheritance_modes"] == ["Autosomal dominant inheritance"]
    assert by_omim["penetrance_hpo_ids"] == ["HP:0003829"]
    assert by_omim["onset_hpo_ids"] == ["HP:0003584"]
    assert by_omim["hpo_assertions"][1] == {
        "hpo_id": "HP:0003829",
        "frequency": "2/10",
        "evidence": "PCS",
        "reference": "PMID:111",
    }


def test_enrich_condition_mechanism_assertion_requires_gene_and_condition() -> None:
    hpo_index = {
        ("GENE1", "MONDO:1"): {
            "disease_id": "OMIM:1",
            "mondo_id": "MONDO:1",
            "disease_scope": "mendelian_non_neoplastic",
            "priva_scope": "include",
            "scope_review_status": "auto_supported",
            "inheritance_modes": ["Autosomal dominant inheritance"],
            "penetrance_hpo_ids": ["HP:0003829"],
            "onset_hpo_ids": ["HP:0003584"],
            "hpo_assertions": [
                {
                    "hpo_id": "HP:0003829",
                    "frequency": "2/10",
                    "evidence": "PCS",
                    "reference": "PMID:111",
                }
            ],
        }
    }
    assertion = {
        "source": "Orphadata",
        "source_condition_id": "ORPHA:1",
        "mondo_id": "MONDO:1",
        "disease": "Condition one",
        "mechanism": "GOF",
        "allelic_requirement": "",
        "priva_scope": "",
    }

    matched = enrich_condition_mechanism_assertion(
        assertion,
        gene_symbol="GENE1",
        hpo_condition_index=hpo_index,
    )
    wrong_gene = enrich_condition_mechanism_assertion(
        assertion,
        gene_symbol="GENE2",
        hpo_condition_index=hpo_index,
    )

    assert len(matched) == 1
    assert matched[0]["allelic_requirement"] == "dominant"
    assert matched[0]["penetrance_hpo_ids"] == ["HP:0003829"]
    assert matched[0]["hpo_assertions"][0]["frequency"] == "2/10"
    assert wrong_gene[0]["hpo_match_status"] == (
        "no_matching_gene_condition_hpo_record"
    )
    assert wrong_gene[0]["penetrance_hpo_ids"] == []


def test_enrich_condition_mechanism_assertion_blocks_review_scope() -> None:
    assertion = {
        "source": "G2P_DDG2P",
        "source_condition_id": "OMIM:2",
        "mondo_id": "MONDO:2",
        "mechanism": "LOF",
        "priva_scope": "review",
    }

    assert (
        enrich_condition_mechanism_assertion(
            assertion,
            gene_symbol="GENE1",
            hpo_condition_index={},
        )
        == []
    )


def test_select_condition_histories_for_variant_is_mechanism_specific() -> None:
    assertions = [
        {"disease": "LOF disorder", "mechanism": "LOF"},
        {"disease": "GOF disorder", "mechanism": "GOF"},
        {"disease": "DN disorder", "mechanism": "DOMINANT_NEGATIVE"},
    ]

    exact_gof = select_condition_histories_for_variant(
        assertions,
        variant_effect="exact_known_GOF",
    )
    predicted_lof = select_condition_histories_for_variant(
        assertions,
        variant_effect="predicted_LOF_high_confidence",
    )
    unresolved = select_condition_histories_for_variant(
        assertions,
        variant_effect="uncertain",
    )
    conflict = select_condition_histories_for_variant(
        assertions,
        variant_effect="exact_known_GOF",
        variant_effect_conflict="predicted_LOF_vs_exact_GOF",
    )

    assert [row["mechanism"] for row in exact_gof] == ["GOF"]
    assert [row["mechanism"] for row in predicted_lof] == ["LOF"]
    assert [row["mechanism"] for row in unresolved] == [
        "LOF",
        "GOF",
        "DOMINANT_NEGATIVE",
    ]
    assert [row["mechanism"] for row in conflict] == ["LOF", "GOF"]
