import json
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_hpo_condition_mechanism_cache import (  # noqa: E402
    attach_condition_mechanisms,
    attach_gofcards_variants,
    build_cache_payload,
    build_hpo_gene_condition_frame,
    load_and_validate_cache,
    write_json_atomic,
)


HPO_HEADER = (
    "gene_symbol\tdisease_id\thpo_id\tfrequency\tevidence\treference\t"
    "mondo_id\tmondo_name\tdisease_scope\tpriva_scope\tscope_evidence\t"
    "scope_reference\tscope_review_status\n"
)


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
                "hpo_id": "HP:0003829",
                "frequency": "2/10",
                "evidence": "PCS",
                "reference": "PMID:1",
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
        "gene_symbol\tsource\tsource_record_id\tsource_condition_id\tmondo_id\t"
        "disease_scope\tpriva_scope\tscope_review_status\tdisease_label\t"
        "inheritance\tpatho_mode_raw\tnormalized_mechanisms\t"
        "mechanism_confidence\tdisease_confidence\tpmids\tevidence_url\n"
        "GENE1\tG2P_DDG2P\tG2P1\tOMIM:1\tMONDO:1\t"
        "mendelian_non_neoplastic\tinclude\tauto_supported\tCondition one\t"
        "monoallelic_autosomal\tloss of function\tLOF\thigh\tdefinitive\t"
        "PMID:1;PMID:2\thttps://g2p.example/1\n"
        "GENE1\tOrphadata\tORPHA2\tORPHA:2\tMONDO:2\t"
        "mendelian_non_neoplastic\tinclude\tauto_supported\tCondition two\t"
        "\tgain of function\tGOF\thigh\tAssessed\tPMID:3\t"
        "https://orpha.example/2\n"
        "GENE2\tG2P_DDG2P\tG2P3\tOMIM:3\tMONDO:3\t"
        "mendelian_non_neoplastic\tinclude\tauto_supported\tCondition three\t"
        "biallelic_autosomal\tloss of function\tLOF\thigh\tdefinitive\t"
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
        "matched": 1,
        "unmapped": 2,
    }


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
    gofcards = tmp_path / "gofcards.tsv"
    gofcards.write_text(
        "mechanism\tHGNC_Symbol\tGoFCards_HGNC_Symbol\tVEP_HGNC_Symbol\t"
        "gene_match_status\tmatch_eligibility\tHGVSc\tHGVSp\thgvsp_key\tmatch_status\t"
        "gofcards_accession_id\tgofcards_variant_id\tdisease\tpmids\tpscore\t"
        "function\tpathway\tallele_key\thg19_genomic_key\thg19_vcf_key\t"
        "hg38_genomic_key\thg38_vcf_key\n"
        "GOF\tGENE1\tGENE1\tGENE1\tgene_concordant\teligible\t"
        "NM_1:c.1A>G\tNP_1:p.Lys1Arg\tNP_1:p.Lys1Arg\tmatched\t"
        "rs1\tVAR1\tCondition one\tPMID:1\t5\tActivating\tPathway\tVAR1\t"
        "1|10|A|G\t1|10|A|G\t1|20|A|G\t1|20|A|G\n"
        "GOF\tGENE1\tGENE1\tGENE1\tgene_concordant\teligible\t"
        "NM_1:c.2A>G\tNP_1:p.Lys2Arg\tNP_1:p.Lys2Arg\tmatched\t"
        "rs2\tVAR2\tCondition two\tPMID:2\t3\tActivating\tPathway\tVAR2\t"
        "1|11|A|G\t1|11|A|G\t1|21|A|G\t1|21|A|G\n"
        "GOF\tGENE1\tGENE1\tOTHER1\tgene_discordant\t"
        "quarantined_gene_discordance\tNM_1:c.3A>G\tNP_1:p.Lys3Arg\t"
        "NP_1:p.Lys3Arg\tmatched\trs3\tVAR3\tCondition three\tPMID:3\t5\t"
        "Activating\tPathway\tVAR3\t1|12|A|G\t1|12|A|G\t1|22|A|G\t1|22|A|G\n"
        "GOF\tGENE1\tGENE1\tGENE1\tgene_concordant\teligible\t"
        "NM_1:c.3A>G\tNP_1:p.Lys3Arg\tNP_1:p.Lys3Arg\tmatched\t"
        "rs3\tVAR3\tCondition three\tPMID:3\t5\tActivating\tPathway\tVAR3\t"
        "1|12|A|G\t1|12|A|G\t1|22|A|G\t1|22|A|G\n",
        encoding="utf-8",
    )
    mechanism_json = tmp_path / "mechanisms.json"
    mechanism_json.write_text(
        json.dumps(
            {
                "HGNC:1": {
                    "symbol": "GENE1",
                    "variant_level": [
                        {
                            "ClinVar_VCV": {
                                "variation": {
                                    "vcv_accession": "VCV0001",
                                    "hgvs": [
                                        {"expression": "NM_1.1:c.1A>G"},
                                        {"expression": "NP_1.1:p.Lys1Arg"},
                                    ],
                                },
                                "match": {
                                    "matched_gofcards_records": [
                                        {
                                            "gofcards_variant_id": "VAR1",
                                            "gofcards_accession_id": "rs1",
                                        }
                                    ]
                                },
                                "condition_assertions": [
                                    {
                                        "rcv_accession": "RCV0001",
                                        "conditions": [
                                            {
                                                "database": "MedGen",
                                                "id": "C1",
                                                "name": "Condition one",
                                            }
                                        ],
                                        "matched_scvs": [
                                            {
                                                "trait_mappings": [
                                                    {
                                                        "mapping_ref": "OMIM",
                                                        "mapping_value": "1",
                                                    }
                                                ]
                                            }
                                        ],
                                        "germline_classification": {
                                            "clinical_significance": "Pathogenic",
                                            "review_stars": 2,
                                        },
                                    }
                                ],
                            }
                        }
                    ],
                }
            }
        ),
        encoding="utf-8",
    )
    genes = build_hpo_gene_condition_frame(hpo)

    stats = attach_gofcards_variants(genes, gofcards, mechanism_json)

    gof = genes["GENE1"]["conditions"]["OMIM:1"][
        "pathogenic_mechanisms"
    ]["GOF"]
    exact = gof["variants"]["GOFCARDS:VAR1"]
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
    assert exact["match_keys"]["GRCh38"] == ["1|20|A|G"]
    unresolved = genes["GENE1"]["unmapped_evidence"]["variants"][
        "GOFCARDS:VAR2"
    ]
    assert unresolved["condition_link"]["reason"] == (
        "no_exact_clinvar_condition_link"
    )
    assert genes["GENE1"]["summary"]["pathogenic_mechanisms"] == ["GOF"]
    assert stats == {
        "source_rows": 4,
        "eligible_source_rows": 2,
        "quarantined_source_rows": 2,
        "quarantined_unique_variants": 1,
        "unique_variants": 2,
        "condition_linked_variants": 1,
        "condition_variant_links": 1,
        "unmapped_variants": 1,
    }


def test_complete_cache_is_validated_and_published_atomically(
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
        "gene_symbol\tsource\tsource_record_id\tsource_condition_id\tmondo_id\t"
        "disease_scope\tpriva_scope\tscope_review_status\tdisease_label\t"
        "inheritance\tpatho_mode_raw\tnormalized_mechanisms\t"
        "mechanism_confidence\tdisease_confidence\tpmids\tevidence_url\n",
        encoding="utf-8",
    )
    gofcards = tmp_path / "gofcards.tsv"
    gofcards.write_text(
        "mechanism\tHGNC_Symbol\tGoFCards_HGNC_Symbol\tVEP_HGNC_Symbol\t"
        "gene_match_status\tmatch_eligibility\tHGVSc\tHGVSp\thgvsp_key\tmatch_status\t"
        "gofcards_accession_id\tgofcards_variant_id\tdisease\tpmids\tpscore\t"
        "function\tpathway\tallele_key\thg19_genomic_key\thg19_vcf_key\t"
        "hg38_genomic_key\thg38_vcf_key\n",
        encoding="utf-8",
    )
    mechanism_json = tmp_path / "mechanisms.json"
    mechanism_json.write_text(json.dumps({"_meta": {}}), encoding="utf-8")
    output = tmp_path / "cache.json"
    payload = build_cache_payload(
        hpo_assertions=hpo,
        mechanism_evidence=mechanism,
        mechanism_json=mechanism_json,
        gofcards_variants=gofcards,
        hpo_release="v1",
        mondo_release="v2",
    )

    write_json_atomic(payload, output)

    assert load_and_validate_cache(output)["genes"] == 1
    assert json.loads(output.read_text(encoding="utf-8"))["_meta"]["releases"] == {
        "HPO": "v1",
        "MONDO": "v2",
    }
    assert not list(tmp_path.glob(".cache.json.*.tmp"))

    original = output.read_bytes()
    with pytest.raises(ValueError, match="schema_version"):
        write_json_atomic({"_meta": {}, "genes": {}}, output)
    assert output.read_bytes() == original
