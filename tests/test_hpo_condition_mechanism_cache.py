import json
import os
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_hpo_condition_mechanism_cache import (  # noqa: E402
    attach_clingen_haploinsufficiency,
    attach_condition_mechanisms,
    attach_gofcards_variants,
    apply_condition_scope_gate,
    build_cache_payload,
    build_hpo_gene_condition_frame,
    deduce_mechanisms_from_inheritance,
    load_and_validate_cache,
    validate_cache_provenance,
    write_json_atomic,
)
from clinvar_vcv import partition_gofcards_variants  # noqa: E402


HPO_HEADER = (
    "gene_symbol\tdisease_id\thpo_id\tfrequency\tevidence\treference\t"
    "mondo_id\tmondo_name\tdisease_scope\tpriva_scope\tscope_evidence\t"
    "scope_reference\tscope_review_status\n"
)

MECHANISM_HEADER = (
    "gene_symbol\tsource\tsource_record_id\tsource_condition_id\tmondo_id\t"
    "disease_scope\tpriva_scope\tscope_review_status\tdisease_label\t"
    "inheritance\tpenetrance_raw\tpenetrance_hpo_ids\tnormalized_penetrance\t"
    "patho_mode_raw\tnormalized_mechanisms\tmechanism_confidence\t"
    "disease_confidence\tpmids\tevidence_url\n"
)


def _gofcards_variant(
    chrom: str,
    pos: int,
    disease: str,
    pmid: str,
    *,
    eligibility: str = "eligible",
    vep_symbol: str = "GENE1",
    gene_match_status: str = "gene_concordant",
    clinvar: dict | None = None,
) -> dict:
    """One variant in the shape the injected GoFCards cache stores."""
    variant = {
        "record": {
            "source": {
                "gofcards_allele_key": f"SNV|{chrom}|{pos}|{pos}|A|G",
                "variant_type_label": "SNV",
                "assembly": "hg19", "chrom": chrom, "start": str(pos),
                "ref": "A", "alt": "G",
            },
            "eligibility": {
                "status": eligibility, "gene_match_status": gene_match_status,
                "vep_symbol": vep_symbol,
                "reason": None if eligibility == "eligible" else "gene_discordant",
            },
            "liftover_status": "mapped",
            "annotations": {"clinvar_variation_id": f"ID{pos}"},
            "evidence": [{"pmid": pmid, "disease": disease, "pscore": "5",
                          "function": "Activating", "pathway": "Pathway"}],
        },
        "assemblies": {
            "hg19": {
                "genomic": {"chrom": chrom, "pos": pos, "ref": "A", "alt": "G",
                            "status": "raw_ref_alt"},
                "transcripts": {"ENST0000000001.1": {
                    "by_hgvsc": {f"c.{pos}A>G": {
                        "hgvsp": f"p.Lys{pos}Arg", "consequence": "missense_variant",
                        "canonical": True, "mane_select": "NM_1.1"}},
                    "by_hgvsp": {f"p.Lys{pos}Arg": [f"c.{pos}A>G"]}}},
            },
            "hg38": {
                "genomic": {"chrom": chrom, "pos": pos + 10, "ref": "A", "alt": "G",
                            "status": "lifted_ref_match"},
                "transcripts": {"ENST0000000001.2": {
                    "by_hgvsc": {f"c.{pos}A>G": {
                        "hgvsp": f"p.Lys{pos}Arg", "consequence": "missense_variant",
                        "canonical": True, "mane_select": "NM_1.1"}},
                    "by_hgvsp": {f"p.Lys{pos}Arg": [f"c.{pos}A>G"]}}},
            },
        },
    }
    if clinvar:
        variant["clinvar"] = clinvar
    return variant


def _write_gofcards_cache(path: Path, variants: dict[str, dict]) -> Path:
    """Write a nested GoFCards cache keyed gene -> variant identifier."""
    path.write_text(
        json.dumps({
            "metadata": {"source": "GoFCards", "mechanism": "GOF",
                         "derived_on": "2026-07-28"},
            "genes": {"GENE1": {"hgnc_id": "HGNC:1", "variants": variants}},
        }),
        encoding="utf-8",
    )
    return path


def test_variant_partition_follows_the_upstream_eligibility_verdict() -> None:
    # A quarantined variant must never reach the condition cache, and its
    # quarantine is decided once by the normalizer and stored on the variant.
    cache = {
        "metadata": {"source": "GoFCards", "mechanism": "GOF"},
        "genes": {
            "GENE1": {"hgnc_id": "HGNC:1", "variants": {
                "loc_1:10:A->G_grch37": _gofcards_variant("1", 10, "d", "1"),
                "loc_1:11:A->G_grch37": _gofcards_variant(
                    "1", 11, "d", "1", eligibility="quarantined_gene_discordance"),
            }},
            "GENE2": {"hgnc_id": "HGNC:2", "variants": {
                "loc_2:10:A->G_grch37": _gofcards_variant("2", 10, "d", "1"),
            }},
        },
    }

    eligible, quarantined = partition_gofcards_variants(cache)

    assert [(s, v) for s, v, _ in eligible] == [
        ("GENE1", "loc_1:10:A->G_grch37"),
        ("GENE2", "loc_2:10:A->G_grch37"),
    ]
    assert [(s, v) for s, v, _ in quarantined] == [("GENE1", "loc_1:11:A->G_grch37")]


def test_variant_partition_rejects_a_variant_reviewed_as_loss_of_function() -> None:
    cache = {
        "metadata": {"source": "GoFCards", "mechanism": "GOF"},
        "genes": {"CFTR": {"hgnc_id": "HGNC:1", "variants": {
            "loc_7:10:A->G_grch37": _gofcards_variant(
                "7", 10, "d", "1", eligibility="quarantined_reviewed_lof"),
        }}},
    }

    eligible, quarantined = partition_gofcards_variants(cache)

    assert eligible == []
    assert [(s, v) for s, v, _ in quarantined] == [("CFTR", "loc_7:10:A->G_grch37")]


def test_recessive_condition_without_curated_mechanism_gets_deduced_lof() -> None:
    def condition(mode: str, decision: str = "include") -> dict:
        return {
            "label": f"{mode} condition",
            "priva_scope": {"decision": decision},
            "inheritance": {"modes": [mode]},
            "pathogenic_mechanisms": {},
        }

    genes = {
        "GENE1": {
            "conditions": {
                "OMIM:1": condition("autosomal_recessive"),
                "OMIM:2": condition("autosomal_dominant"),
                "OMIM:3": condition("autosomal_recessive", decision="review"),
            }
        }
    }

    stats = deduce_mechanisms_from_inheritance(genes)

    evidence = genes["GENE1"]["conditions"]["OMIM:1"][
        "pathogenic_mechanisms"
    ]["LOF"]["evidence"]
    assert evidence == [
        {
            "source": "deduced_from_inheritance",
            "source_record_id": "OMIM:1",
            "condition_identifiers": ["OMIM:1"],
            "condition_label": "autosomal_recessive condition",
            "mechanism": "LOF",
            "mechanism_raw": "recessive inheritance implies loss of function",
            "allelic_requirement": "",
            "mechanism_confidence": "",
            "disease_confidence": "",
            "assertion_basis": "deduced",
            "source_scope": {"decision": "", "category": "", "review_status": ""},
            "pmids": [],
            "evidence_url": "",
        }
    ]
    assert genes["GENE1"]["conditions"]["OMIM:2"]["pathogenic_mechanisms"] == {}
    assert genes["GENE1"]["conditions"]["OMIM:3"]["pathogenic_mechanisms"] == {}
    assert stats == {
        "recessive_conditions_given_lof": 1,
        "left_unresolved_dominant": 1,
    }


def test_x_recessive_but_not_unspecified_x_linked_condition_gets_deduced_lof(
) -> None:
    genes = {
        "GENEX": {
            "conditions": {
                "OMIM:1": {
                    "label": "X-linked recessive condition",
                    "priva_scope": {"decision": "include"},
                    "inheritance": {"modes": ["x_linked_recessive"]},
                    "pathogenic_mechanisms": {},
                },
                "OMIM:2": {
                    "label": "Unspecified X-linked condition",
                    "priva_scope": {"decision": "include"},
                    "inheritance": {"modes": ["x_linked"]},
                    "pathogenic_mechanisms": {},
                },
            }
        }
    }

    stats = deduce_mechanisms_from_inheritance(genes)

    assert "LOF" in genes["GENEX"]["conditions"]["OMIM:1"][
        "pathogenic_mechanisms"
    ]
    assert genes["GENEX"]["conditions"]["OMIM:2"][
        "pathogenic_mechanisms"
    ] == {}
    assert stats == {
        "x_linked_recessive_conditions_given_lof": 1,
        "left_unresolved_x_linked_unspecified": 1,
    }


def test_hpo_frame_groups_gene_conditions_and_preserves_axis_evidence(
    tmp_path: Path,
) -> None:
    hpo = tmp_path / "hpo.tsv"
    hpo.write_text(
        HPO_HEADER
        + "GENE1\tOMIM:1\tHP:0000006\t-\tTAS\tOMIM:1\tMONDO:1\t"
        "Condition one\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1;OMIM:1\tauto_supported\n"
        + "GENE1\tOMIM:1\tHP:0003829\t2/10\tPCS\tPMID:1\tMONDO:1\t"
        "Condition one\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1;OMIM:1\tauto_supported\n"
        + "GENE1\tOMIM:1\tHP:0003584\t1/10\tPCS\tPMID:2\tMONDO:1\t"
        "Condition one\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1;OMIM:1\tmanually_confirmed\n"
        + "GENE1\tOMIM:1\tHP:0001250\t5/10\tPCS\tPMID:3\tMONDO:1\t"
        "Condition one\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1;OMIM:1\tauto_supported\n",
        encoding="utf-8",
    )

    genes = build_hpo_gene_condition_frame(hpo)

    gene = genes["GENE1"]
    condition = gene["conditions"]["OMIM:1"]
    assert gene["summary"]["condition_count"] == 1
    assert gene["summary"]["conditions_with_inheritance"] == 1
    assert gene["summary"]["conditions_with_penetrance"] == 1
    assert condition["identifiers"] == {
        "MONDO": ["MONDO:1"],
        "OMIM": ["OMIM:1"],
    }
    assert condition["inheritance"]["modes"] == ["autosomal_dominant"]
    assert condition["penetrance"] == {
        "statuses": ["incomplete"],
        "assertions": [
            {
                "source": "HPO",
                "source_record_id": "OMIM:1|HP:0003584|PMID:2",
                "condition_identifiers": ["OMIM:1", "MONDO:1"],
                "raw_values": ["HP:0003584"],
                "normalized_value": "incomplete",
                "hpo_ids": ["HP:0003584"],
                "frequency": "1/10",
                "evidence": "PCS",
                "reference": "PMID:2",
                "evidence_url": "",
            },
            {
                "source": "HPO",
                "source_record_id": "OMIM:1|HP:0003829|PMID:1",
                "condition_identifiers": ["OMIM:1", "MONDO:1"],
                "raw_values": ["HP:0003829"],
                "normalized_value": "incomplete",
                "hpo_ids": ["HP:0003829"],
                "frequency": "2/10",
                "evidence": "PCS",
                "reference": "PMID:1",
                "evidence_url": "",
            }
        ],
    }
    assert condition["onset"]["terms"] == ["late"]
    assert condition["hpo_assertion_count"] == 4
    assert condition["priva_scope"]["references"] == [
        "MONDO:v1",
        "OMIM:1",
    ]
    assert condition["priva_scope"]["review_status"] == "manually_confirmed"
    assert condition["priva_scope"]["review_statuses"] == [
        "auto_supported",
        "manually_confirmed",
    ]


def test_final_scope_gate_prunes_only_explicit_exclusions() -> None:
    def condition(decision: str) -> dict:
        return {
            "label": decision,
            "identifiers": {"OMIM": [f"OMIM:{decision}"]},
            "priva_scope": {"decision": decision},
            "inheritance": {"modes": [], "assertions": []},
            "penetrance": {"statuses": [], "assertions": []},
            "onset": {"terms": [], "assertions": []},
            "pathogenic_mechanisms": {},
            "hpo_assertion_count": 1,
        }

    genes = {
        "GENE1": {
            "conditions": {
                "OMIM:include": condition("include"),
                "OMIM:review": condition("review"),
                "OMIM:exclude": condition("exclude"),
            },
            "summary": {},
            "unmapped_evidence": {
                "mechanisms": [],
                "penetrance": [],
                "variants": {},
            },
        },
        "GENE2": {
            "conditions": {"OMIM:exclude": condition("exclude")},
            "summary": {},
            "unmapped_evidence": {
                "mechanisms": [],
                "penetrance": [],
                "variants": {},
            },
        },
    }

    stats = apply_condition_scope_gate(genes)

    assert set(genes["GENE1"]["conditions"]) == {"OMIM:include", "OMIM:review"}
    assert genes["GENE1"]["scope_gate"] == {"excluded_condition_count": 1}
    assert genes["GENE2"]["conditions"] == {}
    assert genes["GENE2"]["scope_gate"] == {"excluded_condition_count": 1}
    assert stats == {
        "conditions_examined": 4,
        "conditions_excluded": 2,
        "conditions_retained_include": 1,
        "conditions_retained_review": 1,
        "conditions_retained_unscoped": 0,
        "genes_with_excluded_conditions": 2,
        "genes_left_without_conditions": 1,
    }


@pytest.mark.parametrize(
    ("hpo_id", "expected_status", "expected_onset"),
    [
        ("HP:0003829", "incomplete", None),
        ("HP:0003831", None, None),
        ("HP:4000159", "incomplete", None),
        ("HP:4000160", "incomplete", None),
        ("HP:0034857", None, "variable_age"),
        ("HP:0003581", None, "adult"),
        ("HP:0011462", None, "young_adult"),
        ("HP:0003596", None, "middle_age"),
        ("HP:0003584", "incomplete", "late"),
        ("HP:0003587", None, "insidious"),
        ("HP:0003828", None, None),
        ("HP:0034950", "complete", None),
        ("HP:4000158", "complete", None),
        ("HP:0001470", None, None),
        ("HP:0003677", None, None),
    ],
)
def test_penetrance_assertions_stay_linked_to_their_condition(
    tmp_path: Path,
    hpo_id: str,
    expected_status: str | None,
    expected_onset: str | None,
) -> None:
    hpo = tmp_path / "hpo.tsv"
    hpo.write_text(
        HPO_HEADER
        + f"GENE1\tOMIM:1\t{hpo_id}\t-\tTAS\tOMIM:1\tMONDO:1\t"
        "Condition one\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1\tauto_supported\n",
        encoding="utf-8",
    )

    condition = build_hpo_gene_condition_frame(hpo)["GENE1"]["conditions"]["OMIM:1"]

    if expected_status is None:
        assert condition["penetrance"] == {"statuses": [], "assertions": []}
    else:
        assert condition["penetrance"]["statuses"] == [expected_status]
        assert [
            assertion["hpo_ids"][0]
            for assertion in condition["penetrance"]["assertions"]
        ] == [hpo_id]

    if expected_onset is None:
        assert condition["onset"] == {"terms": [], "assertions": []}
    else:
        assert condition["onset"]["terms"] == [expected_onset]
        assert condition["onset"]["assertions"][0]["hpo_id"] == hpo_id


def test_mechanisms_attach_only_through_exact_condition_identifiers(
    tmp_path: Path,
) -> None:
    hpo = tmp_path / "hpo.tsv"
    hpo.write_text(
        HPO_HEADER
        + "GENE1\tOMIM:1\tHP:0000006\t-\tTAS\tOMIM:1\tMONDO:1\t"
        "Condition one\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1\tauto_supported\n",
        encoding="utf-8",
    )
    mechanism = tmp_path / "mechanism.tsv"
    mechanism.write_text(
        MECHANISM_HEADER
        + "GENE1\tG2P_DDG2P\tG2P1\tOMIM:1\tMONDO:1\t"
        "mendelian_non_neoplastic\tinclude\tauto_supported\tCondition one\t"
        "monoallelic_autosomal\t\t\t\tloss of function\tLOF\thigh\tdefinitive\t"
        "PMID:1;PMID:2\thttps://g2p.example/1\n"
        "GENE1\tOrphadata\tORPHA2\tORPHA:2\tMONDO:2\t"
        "mendelian_non_neoplastic\tinclude\tauto_supported\tCondition two\t"
        "\t\t\t\tgain of function\tGOF\thigh\tAssessed\tPMID:3\t"
        "https://orpha.example/2\n"
        "GENE2\tG2P_DDG2P\tG2P3\tOMIM:3\tMONDO:3\t"
        "mendelian_non_neoplastic\tinclude\tauto_supported\tCondition three\t"
        "biallelic_autosomal\t\t\t\tloss of function\tLOF\thigh\tdefinitive\t"
        "PMID:4\thttps://g2p.example/3\n",
        encoding="utf-8",
    )
    genes = build_hpo_gene_condition_frame(hpo)

    stats = attach_condition_mechanisms(genes, mechanism)

    lof = genes["GENE1"]["conditions"]["OMIM:1"][
        "pathogenic_mechanisms"
    ]["LOF"]
    assert lof["allelic_requirements"] == ["monoallelic_autosomal"]
    assert lof["evidence"][0]["source_record_id"] == "G2P1"
    assert lof["evidence"][0]["pmids"] == ["PMID:1", "PMID:2"]
    assert genes["GENE1"]["unmapped_evidence"]["mechanisms"][0][
        "reason"
    ] == "no_exact_hpo_condition_identifier"
    assert genes["GENE1"]["summary"]["pathogenic_mechanisms"] == ["GOF", "LOF"]
    assert genes["GENE1"]["summary"]["condition_resolved_mechanism_counts"] == {
        "LOF": 1
    }
    assert genes["GENE2"]["conditions"] == {}
    assert genes["GENE2"]["summary"]["unmapped_mechanism_count"] == 1
    assert stats == {
        "source_rows": 3,
        "mechanism_records": 3,
        "mechanisms_matched": 1,
        "mechanisms_unmapped": 2,
        "penetrance_records": 0,
        "penetrance_matched": 0,
        "penetrance_unmapped": 0,
    }


def test_g2p_penetrance_attaches_without_a_resolved_molecular_mechanism(
    tmp_path: Path,
) -> None:
    hpo = tmp_path / "hpo.tsv"
    hpo.write_text(
        HPO_HEADER
        + "GENE1\tOMIM:1\tHP:0000006\t-\tTAS\tOMIM:1\tMONDO:1\t"
        "Condition one\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1\tauto_supported\n",
        encoding="utf-8",
    )
    mechanism = tmp_path / "mechanism.tsv"
    mechanism.write_text(
        MECHANISM_HEADER
        + "GENE1\tG2P_DDG2P\tG2P03019\tOMIM:1\tMONDO:1\t"
        "mendelian_non_neoplastic\tinclude\tauto_supported\tCondition one\t"
        "monoallelic_autosomal\ttypified by incomplete penetrance\t"
        "HP:0003829\tincomplete\tundetermined\t\t\tdefinitive\tPMID:1\t"
        "https://g2p.example/1\n",
        encoding="utf-8",
    )
    genes = build_hpo_gene_condition_frame(hpo)

    stats = attach_condition_mechanisms(genes, mechanism)

    condition = genes["GENE1"]["conditions"]["OMIM:1"]
    assert condition["pathogenic_mechanisms"] == {}
    assert condition["penetrance"]["statuses"] == ["incomplete"]
    assert condition["penetrance"]["assertions"] == [
        {
            "source": "G2P_DDG2P",
            "source_record_id": "G2P03019",
            "condition_identifiers": ["OMIM:1", "MONDO:1"],
            "raw_values": [
                "typified by incomplete penetrance",
                "HP:0003829",
            ],
            "normalized_value": "incomplete",
            "hpo_ids": ["HP:0003829"],
            "frequency": "",
            "evidence": "",
            "reference": "",
            "evidence_url": "https://g2p.example/1",
        }
    ]
    assert stats["mechanism_records"] == 0
    assert stats["penetrance_records"] == 1
    assert stats["penetrance_matched"] == 1


def test_gofcards_variants_require_exact_clinvar_hpo_condition_identity(
    tmp_path: Path,
) -> None:
    hpo = tmp_path / "hpo.tsv"
    hpo.write_text(
        HPO_HEADER
        + "GENE1\tOMIM:1\tHP:0000006\t-\tTAS\tOMIM:1\tMONDO:1\t"
        "Condition one\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1\tauto_supported\n",
        encoding="utf-8",
    )
    # VAR1 carries the ClinVar block the injection step nests on the variant.
    # Its OMIM identity comes from the submitted record's trait mapping, which
    # is the only place an identifier that joins to HPO appears.
    linked = _gofcards_variant(
        "1", 10, "Condition one", "PMID:1",
        clinvar={
            "vcv_accession": "VCV0001",
            "variation_id": "ID10",
            "hgvs": ["NM_1.1:c.1A>G", "NP_1.1:p.Lys1Arg"],
            "condition_assertions": [
                {
                    "rcv_accession": "RCV0001",
                    "conditions": [
                        {"database": "MedGen", "id": "C1", "name": "Condition one"}
                    ],
                    "matched_scvs": [
                        {"trait_mappings": [
                            {"mapping_ref": "OMIM", "mapping_value": "1"}
                        ]}
                    ],
                    "germline_classification": {
                        "clinical_significance": "Pathogenic", "review_stars": 2
                    },
                }
            ],
        },
    )
    # VAR2 has no ClinVar block, so it cannot reach a condition.
    unlinked = _gofcards_variant("1", 11, "Condition two", "PMID:2")
    # VAR3 is quarantined upstream and must never enter the cache.
    discordant = _gofcards_variant(
        "1", 12, "Condition three", "PMID:3",
        eligibility="quarantined_gene_discordance",
        vep_symbol="OTHER1", gene_match_status="gene_discordant",
    )
    gofcards = _write_gofcards_cache(tmp_path / "gofcards.json", {
        "loc_1:10:A->G_grch37": linked,
        "loc_1:11:A->G_grch37": unlinked,
        "loc_1:12:A->G_grch37": discordant,
    })
    genes = build_hpo_gene_condition_frame(hpo)

    stats = attach_gofcards_variants(genes, gofcards)

    gof = genes["GENE1"]["conditions"]["OMIM:1"][
        "pathogenic_mechanisms"
    ]["GOF"]
    exact = gof["variants"]["GOFCARDS:loc_1:10:A->G_grch37"]
    assert exact["condition_link"] == {
        "status": "exact",
        "condition_key": "OMIM:1",
    }
    assert exact["clinvar_links"][0]["vcv_accession"] == "VCV0001"
    assert exact["clinvar_links"][0]["condition_names"] == ["Condition one"]
    assert exact["clinvar_links"][0]["hgvs"] == [
        "NM_1.1:c.1A>G",
        "NP_1.1:p.Lys1Arg",
    ]
    # The handover is symbol, transcript with its version, and the HGVS
    # notations -- every view on both builds, not a representative one.
    assert exact["symbol"] == "GENE1"
    assert exact["transcripts"] == [
        {"assembly": "hg19", "transcript": "ENST0000000001.1",
         "hgvsc": "c.10A>G", "hgvsp": "p.Lys10Arg"},
        {"assembly": "hg38", "transcript": "ENST0000000001.2",
         "hgvsc": "c.10A>G", "hgvsp": "p.Lys10Arg"},
    ]
    unresolved = genes["GENE1"]["unmapped_evidence"]["variants"][
        "GOFCARDS:loc_1:11:A->G_grch37"
    ]
    assert unresolved["condition_link"]["reason"] == (
        "no_exact_clinvar_condition_link"
    )
    assert genes["GENE1"]["summary"]["pathogenic_mechanisms"] == ["GOF"]
    assert stats == {
        "source_variants": 3,
        "eligible_variants": 2,
        "quarantined_variants": 1,
        "condition_linked_variants": 1,
        "condition_variant_links": 1,
        "unmapped_variants": 1,
    }


def test_curated_gof_precedes_recessive_lof_deduction(tmp_path: Path) -> None:
    """A curated condition mechanism must prevent the fallback inference."""
    hpo = tmp_path / "hpo.tsv"
    hpo.write_text(
        HPO_HEADER
        + "GENE1\tOMIM:1\tHP:0000007\t-\tTAS\tOMIM:1\tMONDO:1\t"
        "Condition one\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1;OMIM:1\tauto_supported\n",
        encoding="utf-8",
    )
    mechanism = tmp_path / "mechanism.tsv"
    mechanism.write_text(
        MECHANISM_HEADER,
        encoding="utf-8",
    )
    clingen = tmp_path / "clingen.tsv"
    clingen.write_text(
        "#Gene Symbol\tGenomic Location\tHaploinsufficiency Score\t"
        "Haploinsufficiency Description\tHaploinsufficiency Disease ID\n",
        encoding="utf-8",
    )
    linked = _gofcards_variant(
        "1",
        10,
        "Condition one",
        "PMID:1",
        clinvar={
            "vcv_accession": "VCV0001",
            "variation_id": "ID10",
            "condition_assertions": [
                {
                    "rcv_accession": "RCV0001",
                    "conditions": [{"name": "Condition one"}],
                    "matched_scvs": [
                        {
                            "trait_mappings": [
                                {"mapping_ref": "OMIM", "mapping_value": "1"}
                            ]
                        }
                    ],
                    "germline_classification": {
                        "clinical_significance": "Pathogenic",
                        "review_stars": 2,
                    },
                }
            ],
        },
    )
    gofcards = _write_gofcards_cache(
        tmp_path / "gofcards.json",
        {"loc_1:10:A->G_grch37": linked},
    )

    payload = build_cache_payload(
        hpo_assertions=hpo,
        mechanism_evidence=mechanism,
        gofcards_variants=gofcards,
        clingen_dosage=clingen,
    )

    condition = payload["genes"]["GENE1"]["conditions"]["OMIM:1"]
    assert set(condition["pathogenic_mechanisms"]) == {"GOF"}
    assert condition["pathogenic_mechanisms"]["GOF"]["evidence"][0][
        "source"
    ] == "GoFCards_exact+ClinVar_VCV"
    assert payload["_meta"]["build_statistics"]["deduced_mechanisms"] == {
        "already_had_a_mechanism": 1
    }


def test_clingen_hi3_uses_condition_linked_x_state_or_chromosome_default(
    tmp_path: Path,
) -> None:
    hpo = tmp_path / "hpo.tsv"
    hpo.write_text(
        HPO_HEADER
        + "GENEA\tOMIM:1\tHP:0000007\t-\tTAS\tOMIM:1\tMONDO:1\t"
        "Autosomal condition\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1\tauto_supported\n"
        + "GENEX\tOMIM:2\tHP:0001419\t-\tTAS\tOMIM:2\tMONDO:2\t"
        "X-linked condition\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1\tauto_supported\n"
        + "GENELOW\tOMIM:3\tHP:0000007\t-\tTAS\tOMIM:3\tMONDO:3\t"
        "Lower-evidence condition\tmendelian_non_neoplastic\tinclude\t"
        "MONDO_ancestor\tMONDO:v1\tauto_supported\n"
        + "GENEXNONE\tOMIM:4\tHP:0001250\t-\tTAS\tOMIM:4\tMONDO:4\t"
        "X-linked condition without inheritance\tmendelian_non_neoplastic\t"
        "include\tMONDO_ancestor\tMONDO:v1\tauto_supported\n"
        + "GENEY\tOMIM:5\tHP:0001250\t-\tTAS\tOMIM:5\tMONDO:5\t"
        "Y-linked condition without inheritance\tmendelian_non_neoplastic\t"
        "include\tMONDO_ancestor\tMONDO:v1\tauto_supported\n",
        encoding="utf-8",
    )
    genes = build_hpo_gene_condition_frame(hpo)
    clingen = tmp_path / "clingen.tsv"
    clingen.write_text(
        "#Gene Symbol\tGenomic Location\tHaploinsufficiency Score\t"
        "Haploinsufficiency Description\tHaploinsufficiency Disease ID\n"
        "GENEA\tchr1:10-20\t3\tSufficient evidence\tMONDO:1\n"
        "GENEX\tchrX:10-20\t3\tSufficient evidence\tMONDO:2\n"
        "GENELOW\tchr2:10-20\t2\tSome evidence\tMONDO:3\n"
        "GENEXNONE\tchrX:30-40\t3\tSufficient evidence\tMONDO:4\n"
        "GENEY\tchrY:30-40\t3\tSufficient evidence\tMONDO:5\n",
        encoding="utf-8",
    )

    stats = attach_clingen_haploinsufficiency(genes, clingen)

    autosomal = genes["GENEA"]["conditions"]["OMIM:1"]
    autosomal_evidence = autosomal["pathogenic_mechanisms"]["LOF"]["evidence"][0]
    assert autosomal["inheritance"]["modes"] == ["autosomal_recessive"]
    assert autosomal_evidence["allelic_requirement"] == "monoallelic_autosomal"

    x_linked = genes["GENEX"]["conditions"]["OMIM:2"]
    x_evidence = x_linked["pathogenic_mechanisms"]["LOF"]["evidence"][0]
    assert x_linked["inheritance"]["modes"] == ["x_linked_recessive"]
    assert x_evidence["allelic_requirement"] == ""

    x_without_inheritance = genes["GENEXNONE"]["conditions"]["OMIM:4"]
    x_without_inheritance_evidence = x_without_inheritance[
        "pathogenic_mechanisms"
    ]["LOF"]["evidence"][0]
    assert x_without_inheritance["inheritance"]["modes"] == []
    assert x_without_inheritance_evidence["allelic_requirement"] == (
        "monoallelic_X_heterozygous"
    )

    y_without_inheritance = genes["GENEY"]["conditions"]["OMIM:5"]
    y_without_inheritance_evidence = y_without_inheritance[
        "pathogenic_mechanisms"
    ]["LOF"]["evidence"][0]
    assert y_without_inheritance["inheritance"]["modes"] == []
    assert y_without_inheritance_evidence["allelic_requirement"] == (
        "monoallelic_Y_hemizygous"
    )

    lower_evidence = genes["GENELOW"]["conditions"]["OMIM:3"]
    assert lower_evidence["pathogenic_mechanisms"]["LOF"]["evidence"][0][
        "allelic_requirement"
    ] == ""
    assert stats["hi3_monoallelic_autosomal"] == 1
    assert stats["hi3_allelic_state_deferred_to_condition_evidence"] == 1
    assert stats["hi3_monoallelic_X_heterozygous"] == 1
    assert stats["hi3_monoallelic_Y_hemizygous"] == 1


def test_complete_cache_is_validated_and_published_atomically(
    tmp_path: Path,
) -> None:
    hpo = tmp_path / "hpo.tsv"
    hpo.write_text(
        HPO_HEADER
        + "GENE1\tOMIM:1\tHP:0000006\t-\tTAS\tOMIM:1\tMONDO:1\t"
        "Condition one\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1\tauto_supported\n"
        + "GENE2\tOMIM:2\tHP:0001426\t-\tTAS\tOMIM:2\tMONDO:2\t"
        "Complex condition\tcomplex_or_non_monogenic\texclude\t"
        "HPO_non_monogenic_inheritance\tMONDO:v1\tauto_supported\n",
        encoding="utf-8",
    )
    mechanism = tmp_path / "mechanism.tsv"
    mechanism.write_text(
        MECHANISM_HEADER,
        encoding="utf-8",
    )
    clingen = tmp_path / "clingen.tsv"
    clingen.write_text(
        "#Gene Symbol\tGenomic Location\tHaploinsufficiency Score\t"
        "Haploinsufficiency Description\tHaploinsufficiency Disease ID\n",
        encoding="utf-8",
    )
    gofcards = _write_gofcards_cache(tmp_path / "gofcards.json", {})
    output = tmp_path / "cache.json"
    payload = build_cache_payload(
        hpo_assertions=hpo,
        mechanism_evidence=mechanism,
        clingen_dosage=clingen,
        gofcards_variants=gofcards,
        hpo_release="v1",
        mondo_release="v2",
    )

    write_json_atomic(payload, output)

    counts = load_and_validate_cache(output)
    assert counts["genes"] == 2
    assert counts["conditions"] == 1
    assert counts["excluded_conditions"] == 1
    assert payload["genes"]["GENE2"]["conditions"] == {}
    assert payload["genes"]["GENE2"]["scope_gate"] == {
        "excluded_condition_count": 1
    }
    assert payload["_meta"]["build_statistics"]["scope_gate"][
        "conditions_excluded"
    ] == 1
    assert json.loads(output.read_text(encoding="utf-8"))["_meta"]["releases"] == {
        "HPO": "v1",
        "MONDO": "v2",
    }
    assert set(payload["_meta"]["sources"]) == {
        "builder_script",
        "hpo_assertions",
        "mechanism_evidence",
        "clingen_dosage",
        "gofcards_variants",
    }
    assert not list(tmp_path.glob(".cache.json.*.tmp"))

    original = output.read_bytes()
    with pytest.raises(ValueError, match="schema_version"):
        write_json_atomic({"_meta": {}, "genes": {}}, output)
    assert output.read_bytes() == original


def test_cache_provenance_detects_changed_inputs(tmp_path: Path) -> None:
    hpo = tmp_path / "hpo.tsv"
    hpo.write_text(
        HPO_HEADER
        + "GENE1\tOMIM:1\tHP:0000006\t-\tTAS\tOMIM:1\tMONDO:1\t"
        "Condition one\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1\tauto_supported\n",
        encoding="utf-8",
    )
    mechanism = tmp_path / "mechanism.tsv"
    mechanism.write_text(
        MECHANISM_HEADER,
        encoding="utf-8",
    )
    clingen = tmp_path / "clingen.tsv"
    clingen.write_text(
        "#Gene Symbol\tGenomic Location\tHaploinsufficiency Score\t"
        "Haploinsufficiency Description\tHaploinsufficiency Disease ID\n",
        encoding="utf-8",
    )
    gofcards = _write_gofcards_cache(tmp_path / "gofcards.json", {})
    output = tmp_path / "cache.json"
    payload = build_cache_payload(
        hpo_assertions=hpo,
        mechanism_evidence=mechanism,
        clingen_dosage=clingen,
        gofcards_variants=gofcards,
        hpo_release="v1",
        mondo_release="v2",
    )
    write_json_atomic(payload, output)

    expected = {
        "hpo_assertions": hpo,
        "mechanism_evidence": mechanism,
        "clingen_dosage": clingen,
        "gofcards_variants": gofcards,
        "hpo_release": "v1",
        "mondo_release": "v2",
    }
    assert validate_cache_provenance(output, **expected)["genes"] == 1

    touched_timestamp_ns = mechanism.stat().st_mtime_ns + 1_000_000_000
    os.utime(mechanism, ns=(touched_timestamp_ns, touched_timestamp_ns))
    assert mechanism.stat().st_mtime_ns == touched_timestamp_ns
    assert validate_cache_provenance(output, **expected)["genes"] == 1

    mechanism.write_text(mechanism.read_text() + "# changed\n", encoding="utf-8")
    with pytest.raises(ValueError, match="mechanism_evidence: changed"):
        validate_cache_provenance(output, **expected)
