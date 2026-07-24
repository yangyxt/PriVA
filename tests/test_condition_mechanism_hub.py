import json
from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gene_mechanism_hub import (  # noqa: E402
    GeneMechanismHub,
    build_hpo_condition_index,
    enrich_condition_mechanism_assertion,
    extract_exact_clinvar_condition_identities,
    select_condition_histories_for_variant,
)


def test_load_condition_mechanism_evidence_preserves_identity_and_mechanism(
    tmp_path: Path,
) -> None:
    """The TSV remains authoritative even when the JSON has no G2P records."""
    hgnc = tmp_path / "hgnc.tsv"
    hgnc.write_text(
        "hgnc_id\tsymbol\tensembl_gene_id\tentrez_id\tprev_symbol\t"
        "alias_symbol\trefseq_accession\tuniprot_ids\tmane_select\n"
        "HGNC:1\tGENE1\tENSG1\t1\tOLD1\t\t\t\t\n",
        encoding="utf-8",
    )
    mechanism_json = tmp_path / "mechanisms.json"
    mechanism_json.write_text(json.dumps({"_meta": {}}), encoding="utf-8")
    hpo = tmp_path / "hpo.tsv"
    hpo.write_text(
        "gene_symbol\tdisease_id\tmondo_id\thpo_id\tfrequency\tevidence\t"
        "reference\tdisease_scope\tpriva_scope\tscope_review_status\n"
        "GENE1\tOMIM:1\tMONDO:0000001\tHP:0000006\t-\tTAS\tOMIM:1\t"
        "mendelian_non_neoplastic\tinclude\tauto_supported\n"
        "GENE1\tOMIM:1\tMONDO:0000001\tHP:0003829\t2/10\tPCS\tPMID:1\t"
        "mendelian_non_neoplastic\tinclude\tauto_supported\n",
        encoding="utf-8",
    )
    evidence = tmp_path / "evidence.tsv"
    columns = [
        "gene_symbol",
        "source",
        "source_record_id",
        "source_condition_id",
        "mondo_id",
        "disease_scope",
        "priva_scope",
        "scope_review_status",
        "disease_label",
        "inheritance",
        "patho_mode_raw",
        "normalized_mechanisms",
        "mechanism_confidence",
        "disease_confidence",
        "pmids",
    ]
    evidence.write_text(
        "\t".join(columns)
        + "\n"
        + "\t".join(
            [
                "OLD1",
                "G2P_DDG2P",
                "G2P0001",
                "OMIM:1",
                "MONDO:0000001",
                "mendelian_non_neoplastic",
                "include",
                "auto_supported",
                "Condition one",
                "monoallelic_autosomal",
                "loss or gain of function",
                "LOF;GOF",
                "high",
                "definitive",
                "11;12",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    hub = GeneMechanismHub(
        mechanism_json=mechanism_json,
        ddg2p_evidence=evidence,
        hpo_collapsed=hpo,
        hgnc_table=hgnc,
        use_hgnc_package=False,
    )

    loaded = hub._load_condition_mechanism_evidence()

    assert loaded["OLD1"] is not loaded["GENE1"]
    assert [row["mechanism"] for row in loaded["GENE1"]] == ["GOF", "LOF"]
    assert loaded["GENE1"][0]["source_condition_id"] == "OMIM:1"
    assert loaded["GENE1"][0]["mondo_id"] == "MONDO:0000001"
    assert loaded["GENE1"][0]["disease_scope"] == "mendelian_non_neoplastic"
    assert loaded["GENE1"][0]["pmids"] == ["11", "12"]

    condition_assertions = hub.condition_mechanism_assertions("GENE1")
    assert [row["mechanism"] for row in condition_assertions] == ["GOF", "LOF"]
    assert all(
        row["hpo_match_status"] == "matched_gene_and_condition"
        for row in condition_assertions
    )
    assert all(
        row["penetrance_hpo_ids"] == ["HP:0003829"]
        for row in condition_assertions
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


def test_extract_exact_clinvar_condition_identities_uses_database_ids_only() -> None:
    condition_assertion = {
        "conditions": [
            {
                "database": "MedGen",
                "id": "C3809233",
                "name": "Noonan syndrome 8",
            }
        ],
        "matched_scvs": [
            {
                "trait_mappings": [
                    {
                        "mapping_ref": "OMIM",
                        "mapping_value": "615355",
                        "medgen_name": "Noonan syndrome 8",
                    },
                    {
                        "mapping_ref": "Orphanet",
                        "mapping_value": "648",
                        "medgen_name": "Noonan syndrome",
                    },
                    {
                        "mapping_ref": "Preferred",
                        "mapping_value": "Noonan syndrome 8",
                        "medgen_name": "Noonan syndrome 8",
                    },
                ]
            }
        ],
    }

    identities = extract_exact_clinvar_condition_identities(condition_assertion)

    assert [row["source_condition_id"] for row in identities] == [
        "MEDGEN:C3809233",
        "OMIM:615355",
        "ORPHA:648",
    ]
    assert all("Preferred" not in row.values() for row in identities)
