import json
from pathlib import Path
import sys

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gene_mechanism_hub import (  # noqa: E402
    GeneMechanismHub,
    build_hpo_condition_index,
    condition_cache_context,
    condition_cache_mechanism_assertions,
    enrich_condition_mechanism_assertion,
    extract_exact_clinvar_condition_identities,
    select_condition_histories_for_variant,
)


def test_load_integrated_condition_cache_validates_schema_and_resolves_alias(
    tmp_path: Path,
) -> None:
    hgnc = tmp_path / "hgnc.tsv"
    hgnc.write_text(
        "hgnc_id\tsymbol\tensembl_gene_id\tentrez_id\tprev_symbol\t"
        "alias_symbol\trefseq_accession\tuniprot_ids\tmane_select\n"
        "HGNC:1\tGENE1\tENSG1\t1\tOLD1\t\t\t\t\n",
        encoding="utf-8",
    )
    cache = tmp_path / "condition-cache.json"
    gene_record = {
        "conditions": {},
        "summary": {"pathogenic_mechanisms": []},
        "unmapped_evidence": {"mechanisms": [], "variants": {}},
    }
    cache.write_text(
        json.dumps(
            {
                "_meta": {"schema_version": "1.0"},
                "genes": {"OLD1": gene_record},
            }
        ),
        encoding="utf-8",
    )

    hub = GeneMechanismHub(
        condition_cache=cache,
        hgnc_table=hgnc,
        use_hgnc_package=False,
    )
    loaded = hub._load_condition_cache()

    assert loaded["OLD1"] is loaded["GENE1"]
    assert hub._condition_cache_meta == {"schema_version": "1.0"}

    cache.write_text(
        json.dumps({"_meta": {"schema_version": "0.9"}, "genes": {}}),
        encoding="utf-8",
    )
    stale_hub = GeneMechanismHub(
        condition_cache=cache,
        hgnc_table=hgnc,
        use_hgnc_package=False,
    )
    with pytest.raises(ValueError, match="Unsupported condition cache schema"):
        stale_hub._load_condition_cache()


def test_condition_cache_assertions_preserve_context_and_exclude_variant_only_gof(
) -> None:
    hpo_assertion = {
        "hpo_id": "HP:0000006",
        "frequency": "-",
        "evidence": "TAS",
        "reference": "OMIM:1",
    }
    included_condition = {
        "label": "Condition one",
        "identifiers": {"OMIM": ["OMIM:1"], "MONDO": ["MONDO:0000001"]},
        "priva_scope": {
            "decision": "include",
            "category": "mendelian_non_neoplastic",
            "review_status": "auto_supported",
        },
        "inheritance": {
            "modes": ["autosomal_dominant"],
            "assertions": [hpo_assertion],
        },
        "penetrance": {
            "statuses": ["incomplete"],
            "assertions": [
                {
                    "hpo_id": "HP:0003829",
                    "frequency": "2/10",
                    "evidence": "PCS",
                    "reference": "PMID:1",
                }
            ],
        },
        "onset": {"terms": [], "assertions": []},
        "hpo_assertion_count": 12,
        "pathogenic_mechanisms": {
            "GOF": {
                "allelic_requirements": ["monoallelic_autosomal"],
                "evidence": [
                    {
                        "source": "G2P_DDG2P",
                        "source_record_id": "G2P0001",
                        "condition_identifiers": ["OMIM:1", "MONDO:0000001"],
                        "condition_label": "Condition one",
                        "mechanism": "GOF",
                        "mechanism_raw": "gain of function",
                        "allelic_requirement": "monoallelic_autosomal",
                        "mechanism_confidence": "high",
                        "disease_confidence": "definitive",
                        "pmids": ["1"],
                    },
                    {
                        "source": "GoFCards_exact+ClinVar_VCV",
                        "source_record_id": "VCV1",
                        "condition_identifiers": ["OMIM:1"],
                        "mechanism": "GOF",
                        "mechanism_confidence": "exact_variant",
                    },
                ],
                "variants": {},
            }
        },
    }
    review_condition = {
        **included_condition,
        "priva_scope": {
            **included_condition["priva_scope"],
            "decision": "review",
        },
    }
    gene = {
        "conditions": {
            "OMIM:1": included_condition,
            "OMIM:2": review_condition,
        }
    }

    assertions = condition_cache_mechanism_assertions(gene)

    assert len(assertions) == 1
    assertion = assertions[0]
    assert assertion["source"] == "G2P_DDG2P"
    assert assertion["source_condition_id"] == "OMIM:1"
    assert assertion["mondo_id"] == "MONDO:0000001"
    assert assertion["hpo_inheritance_modes"] == [
        "Autosomal dominant inheritance"
    ]
    assert assertion["penetrance_hpo_ids"] == ["HP:0003829"]
    assert assertion["hpo_assertion_count"] == 12
    assert assertion["hpo_assertions"][1]["frequency"] == "2/10"

    assert condition_cache_context("OMIM:2", review_condition) == {}


def test_legacy_condition_evidence_loader_remains_audit_only(
    tmp_path: Path,
) -> None:
    """The old TSV parser remains inspectable but no longer drives runtime."""
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

    assert hub.condition_mechanism_assertions("GENE1") == []


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
        "priva_scope": "include",
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


def test_enrich_condition_mechanism_assertion_blocks_unscoped_disease() -> None:
    assertion = {
        "source": "G2P_DDG2P",
        "source_condition_id": "OMIM:3",
        "mechanism": "GOF",
        "priva_scope": "",
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
