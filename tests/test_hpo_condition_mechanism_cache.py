from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_hpo_condition_mechanism_cache import (  # noqa: E402
    attach_condition_mechanisms,
    build_hpo_gene_condition_frame,
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
        + "GENE1\tORPHA:1\tHP:0003584\t1/10\tPCS\tPMID:2\tMONDO:1\t"
        "Condition one\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1;ORPHA:1\tauto_supported\n"
        + "GENE1\tOMIM:1\tHP:0001250\t5/10\tPCS\tPMID:3\tMONDO:1\t"
        "Condition one\tmendelian_non_neoplastic\tinclude\tMONDO_ancestor\t"
        "MONDO:v1;OMIM:1\tauto_supported\n",
        encoding="utf-8",
    )

    genes = build_hpo_gene_condition_frame(hpo)

    gene = genes["GENE1"]
    condition = gene["conditions"]["MONDO:1"]
    assert gene["summary"]["condition_count"] == 1
    assert gene["summary"]["conditions_with_inheritance"] == 1
    assert gene["summary"]["conditions_with_penetrance"] == 1
    assert condition["identifiers"] == {
        "MONDO": ["MONDO:1"],
        "OMIM": ["OMIM:1"],
        "ORPHA": ["ORPHA:1"],
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
        "ORPHA:1",
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

    lof = genes["GENE1"]["conditions"]["MONDO:1"][
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
