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
    condition_cache_mechanism_entries,
    condition_cache_mechanism_assertions,
    enrich_condition_mechanism_assertion,
    gene_has_lof_mechanism_history,
    infer_query_variant_effect,
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


def test_condition_cache_assertions_preserve_context_and_gof_history(
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

    # Both sources are admitted. A curated GoFCards allele is this gene's
    # curated history for the condition it was curated against; it cannot leak
    # onto an unrelated variant because select_condition_histories_for_variant
    # keeps only the histories the query variant's own mechanism reaches.
    assert {assertion["source"] for assertion in assertions} == {
        "G2P_DDG2P",
        "GoFCards_exact+ClinVar_VCV",
    }

    # The review-scoped condition contributes nothing from either source.
    assert {assertion["source_condition_id"] for assertion in assertions} == {"OMIM:1"}

    g2p = next(a for a in assertions if a["source"] == "G2P_DDG2P")
    assert g2p["mondo_id"] == "MONDO:0000001"
    assert g2p["hpo_inheritance_modes"] == ["Autosomal dominant inheritance"]
    assert g2p["penetrance_hpo_ids"] == ["HP:0003829"]
    assert g2p["hpo_assertion_count"] == 12
    assert g2p["hpo_assertions"][1]["frequency"] == "2/10"

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
    assert by_omim["penetrance_hpo_ids"] == ["HP:0003829", "HP:0003584"]
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


def test_mechanism_free_condition_is_unresolved_and_deduced_lof_is_admitted() -> None:
    base_condition = {
        "label": "Condition one",
        "identifiers": {"OMIM": ["OMIM:1"]},
        "priva_scope": {
            "decision": "include",
            "category": "mendelian_non_neoplastic",
            "review_status": "auto_supported",
        },
        "inheritance": {
            "modes": ["autosomal_recessive"],
            "assertions": [
                {
                    "hpo_id": "HP:0000007",
                    "frequency": "-",
                    "evidence": "TAS",
                    "reference": "OMIM:1",
                }
            ],
        },
        "penetrance": {
            "statuses": ["incomplete"],
            "assertions": [
                {
                    "hpo_id": "HP:0003584",
                    "frequency": "-",
                    "evidence": "PCS",
                    "reference": "PMID:1",
                }
            ],
        },
        "onset": {
            "terms": ["late"],
            "assertions": [
                {
                    "hpo_id": "HP:0003584",
                    "frequency": "-",
                    "evidence": "PCS",
                    "reference": "PMID:1",
                }
            ],
        },
        "hpo_assertion_count": 2,
        "pathogenic_mechanisms": {},
    }
    legacy_deduced = json.loads(json.dumps(base_condition))
    legacy_deduced["label"] = "Legacy condition"
    legacy_deduced["pathogenic_mechanisms"] = {
        "LOF": {
            "evidence": [
                {
                    "source": "deduced_from_inheritance",
                    "source_record_id": "OMIM:2",
                    "assertion_basis": "deduced",
                    "mechanism": "LOF",
                }
            ]
        }
    }
    review_only = json.loads(json.dumps(base_condition))
    review_only["priva_scope"]["decision"] = "review"

    assertions = condition_cache_mechanism_assertions(
        {
            "conditions": {
                "OMIM:1": base_condition,
                "OMIM:2": legacy_deduced,
                "OMIM:3": review_only,
            }
        }
    )

    by_condition = {row["hpo_disease_id"]: row for row in assertions}
    assert set(by_condition) == {"OMIM:1", "OMIM:2"}
    assert by_condition["OMIM:1"]["mechanism"] == "UNRESOLVED"
    assert by_condition["OMIM:1"]["source"] == "HPO_condition_history"
    assert by_condition["OMIM:1"]["assertion_basis"] == (
        "condition_without_explicit_mechanism"
    )
    assert by_condition["OMIM:2"]["mechanism"] == "LOF"
    assert by_condition["OMIM:2"]["source"] == "deduced_from_inheritance"
    assert by_condition["OMIM:2"]["assertion_basis"] == "deduced"
    assert all(row["penetrance_hpo_ids"] == ["HP:0003584"] for row in assertions)

    deduced_assertion = by_condition["OMIM:2"]
    assert select_condition_histories_for_variant(
        [deduced_assertion], variant_effect="predicted_LOF_high_confidence"
    ) == [
        {
            "mechanism": "LOF",
            "condition": "OMIM:2",
            "inheritance": "recessive",
            "x_linked": False,
            "penetrance": "incomplete",
        }
    ]
    assert select_condition_histories_for_variant(
        [deduced_assertion], variant_effect="exact_known_GOF"
    ) == []


def test_deduced_lof_does_not_count_as_pvs1_lof_history() -> None:
    condition = {
        "priva_scope": {"decision": "include"},
        "pathogenic_mechanisms": {
            "LOF": {
                "evidence": [
                    {
                        "source": "deduced_from_inheritance",
                        "assertion_basis": "deduced",
                    }
                ]
            }
        },
    }
    gene = {"conditions": {"OMIM:1": condition}}

    assert gene_has_lof_mechanism_history(gene) is False

    condition["pathogenic_mechanisms"]["LOF"]["evidence"].append(
        {"source": "G2P_DDG2P", "assertion_basis": "curated"}
    )
    assert gene_has_lof_mechanism_history(gene) is True


def test_select_condition_histories_for_variant_is_mechanism_specific() -> None:
    assertions = [
        {
            "mechanism": "LOF",
            "hpo_disease_id": "OMIM:1",
            "allelic_requirement": "biallelic_autosomal",
            "penetrance_hpo_ids": ["HP:0003829"],
        },
        {
            "mechanism": "GOF",
            "hpo_disease_id": "OMIM:2",
            "allelic_requirement": "monoallelic_autosomal",
            "penetrance_hpo_ids": [],
        },
        {
            "mechanism": "DOMINANT_NEGATIVE",
            "hpo_disease_id": "OMIM:3",
            "hpo_inheritance_modes": ["x_linked_dominant"],
            "penetrance_hpo_ids": [],
        },
        {
            "mechanism": "UNRESOLVED",
            "hpo_disease_id": "OMIM:4",
            "hpo_inheritance_modes": ["autosomal_recessive"],
            "penetrance_hpo_ids": ["HP:0003584"],
        },
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

    assert [row["mechanism"] for row in exact_gof] == ["GOF", "UNRESOLVED"]
    assert [row["mechanism"] for row in predicted_lof] == ["LOF", "UNRESOLVED"]
    assert [row["mechanism"] for row in unresolved] == [
        "LOF",
        "GOF",
        "DOMINANT_NEGATIVE",
        "UNRESOLVED",
    ]

    # Each selected history carries the three facts the chain delivers, and
    # nothing else travels with it.
    assert exact_gof[0] == {
        "mechanism": "GOF",
        "condition": "OMIM:2",
        "inheritance": "dominant",
        "x_linked": False,
        "penetrance": "unknown",
    }
    assert predicted_lof[0] == {
        "mechanism": "LOF",
        "condition": "OMIM:1",
        "inheritance": "recessive",
        "x_linked": False,
        "penetrance": "incomplete",
    }
    # X-linkage survives as its own fact rather than being folded into the value.
    assert unresolved[2]["inheritance"] == "dominant"
    assert unresolved[2]["x_linked"] is True
    assert exact_gof[1] == {
        "mechanism": "UNRESOLVED",
        "condition": "OMIM:4",
        "inheritance": "recessive",
        "x_linked": False,
        "penetrance": "incomplete",
    }


def test_mixed_mendelian_and_non_monogenic_condition_keeps_both_gates() -> None:
    histories = select_condition_histories_for_variant(
        [
            {
                "mechanism": "GOF",
                "hpo_disease_id": "OMIM:1",
                "allelic_requirement": "monoallelic_autosomal",
                "hpo_inheritance_modes": [
                    "Autosomal dominant inheritance",
                    "Polygenic inheritance",
                ],
                "penetrance_hpo_ids": [],
            }
        ],
        variant_effect="exact_known_GOF",
    )

    assert [(row["inheritance"], row["condition"]) for row in histories] == [
        ("dominant", "OMIM:1"),
        ("polygenic", "OMIM:1"),
    ]


def test_exact_gof_suppresses_but_retains_predicted_lof_evidence() -> None:
    # A curated gain-of-function allele outranks predicted loss of function,
    # including LOFTEE's high-confidence call: those are predictions, this is a
    # curator's verdict on this exact allele. The predictions stay visible.
    effect = infer_query_variant_effect(
        {
            "Consequence": "stop_gained",
            "LoF": "HC",
            "NMD": "NMD_escaping_variant",
            "vep_consq_lof": True,
            "variant_gof_score": 2,
            "variant_dn_score": 0,
            "variant_lof_score": 0,
        }
    )

    assert effect["variant_effect"] == "exact_known_GOF"
    assert effect["variant_gof_score"] == 2
    assert effect["variant_lof_score"] == 0
    assert effect["variant_mechanism_exclusive"] is True
    assert "LOFTEE_HC" in effect["variant_effect_suppressed_evidence"]


def test_nonsense_mediated_decay_outranks_even_a_curated_gof_allele() -> None:
    # Decay destroys the transcript, so there is no protein left to gain a
    # function. This is the one case where a prediction outranks curation, and
    # the overridden curated call is recorded rather than silently dropped.
    effect = infer_query_variant_effect(
        {
            "Consequence": "stop_gained",
            "NMD": "NMD",
            "variant_gof_score": 2,
            "variant_dn_score": 0,
            "variant_lof_score": 0,
        }
    )

    assert effect["variant_effect"] == "exact_known_LOF"
    assert effect["variant_lof_score"] == 2
    assert effect["variant_gof_score"] == 0
    assert effect["variant_mechanism_exclusive"] is True
    assert "CANONICAL_EXACT_GOF" in effect["variant_effect_suppressed_evidence"]
