from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_hpo_condition_mechanism_cache import (  # noqa: E402
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
