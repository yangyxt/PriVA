import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from acmg_criteria_assign import (  # noqa: E402
    BS1_criteria,
    BS2_criteria,
    BS4_criteria,
    PP3_BP4_criteria,
    _variant_mechanism_masks,
    summarize_clinvar_gene_pathogenicity,
    vep_consq_interpret_per_row,
)
from gene_mechanism_hub import (  # noqa: E402
    GeneMechanismHub,
    annotate_gene_mechanism_categories,
)


def _write_fixture_sources(tmp_path: Path) -> dict[str, Path]:
    mechanism_json = tmp_path / "mechanisms.json"
    mechanism_json.write_text(
        json.dumps(
            {
                "HGNC:3": {
                    "symbol": "TESTGOF",
                    "gene_level": [
                        {
                            "G2P_DDG2P": {
                                "disease": "biallelic GOF disorder",
                                "inheritance": "biallelic_autosomal",
                                "mechanism": "GOF",
                                "confidence": "high",
                            }
                        }
                    ],
                    "variant_level": [
                        {
                            "ClinVar_VCV": {
                                "match": {
                                    "matched_gofcards_records": [
                                        {
                                            "gofcards_variant_id": "SNV|1|100|100|A|G",
                                            "gofcards_accession_id": "rs1",
                                        }
                                    ]
                                },
                                "variation": {
                                    "vcv_accession": "VCV000000001",
                                    "hgvs": [
                                        {"expression": "NM_000001.1:c.1A>G"},
                                        {"expression": "NP_000001.1:p.Met1Val"},
                                    ],
                                },
                                "condition_assertions": [
                                    {
                                        "conditions": [
                                            {"name": "biallelic GOF disorder"}
                                        ],
                                        "germline_classification": {
                                            "review_stars": 2
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

    evidence_tsv = tmp_path / "evidence.tsv"
    evidence_tsv.write_text(
        "\t".join(
            [
                "gene_symbol",
                "source",
                "source_record_id",
                "disease_label",
                "inheritance",
                "patho_mode_raw",
                "normalized_mechanisms",
                "mechanism_confidence",
                "disease_confidence",
                "pmids",
            ]
        )
        + "\n"
        + "TESTLOF\tG2P_DDG2P\t1\tbiallelic LOF disorder\t"
        "biallelic_autosomal\tloss of function\tLOF\thigh\tdefinitive\t1\n"
        + "TESTMONO\tG2P_DDG2P\t2\tmonoallelic LOF disorder\t"
        "monoallelic_autosomal\tloss of function\tLOF\thigh\tdefinitive\t2\n",
        encoding="utf-8",
    )

    hpo_tsv = tmp_path / "hpo.tsv"
    hpo_tsv.write_text(
        "gene_symbol\tdisease_id\thpo_id\tfrequency\tevidence\treference\n",
        encoding="utf-8",
    )
    clingen_tsv = tmp_path / "clingen.tsv"
    clingen_tsv.write_text(
        "#Gene Symbol\tHaploinsufficiency Score\tHaploinsufficiency Description\t"
        "Triplosensitivity Score\tTriplosensitivity Description\n",
        encoding="utf-8",
    )
    loeuf_tsv = tmp_path / "loeuf.tsv"
    loeuf_tsv.write_text(
        "#gene\tcanonical\toe_lof_upper\n"
        "TESTREC\ttrue\t1.0\n"
        "TESTLOF\ttrue\t1.0\n"
        "TESTGOF\ttrue\t1.0\n"
        "TESTMONO\ttrue\t1.0\n",
        encoding="utf-8",
    )
    hgnc_tsv = tmp_path / "hgnc.tsv"
    hgnc_tsv.write_text(
        "hgnc_id\tsymbol\tensembl_gene_id\tentrez_id\tprev_symbol\t"
        "alias_symbol\trefseq_accession\tuniprot_ids\tmane_select\n"
        "HGNC:1\tTESTREC\tENSG1\t1\t\t\t\t\t\n"
        "HGNC:2\tTESTLOF\tENSG2\t2\t\t\t\t\t\n"
        "HGNC:3\tTESTGOF\tENSG3\t3\t\t\t\t\t\n"
        "HGNC:4\tTESTMONO\tENSG4\t4\t\t\t\t\t\n",
        encoding="utf-8",
    )
    return {
        "mechanism_json": mechanism_json,
        "ddg2p_evidence": evidence_tsv,
        "hpo_collapsed": hpo_tsv,
        "clingen_dosage": clingen_tsv,
        "loeuf_table": loeuf_tsv,
        "hgnc_table": hgnc_tsv,
    }


def test_loftee_hc_and_os_are_high_confidence_lof() -> None:
    for row in (
        {
            "Consequence": "splice_donor_5th_base_variant&intron_variant",
            "LoF": "OS",
            "NMD": "",
        },
        {
            "Consequence": "stop_gained",
            "LoF": "HC",
            "NMD": "NMD",
        },
    ):
        is_lof, is_length_changing = vep_consq_interpret_per_row(row)
        assert is_lof is True
        assert is_length_changing is True

    is_lof, _ = vep_consq_interpret_per_row(
        {"Consequence": "stop_gained", "LoF": "LC", "NMD": "NMD"}
    )
    assert is_lof is False


def test_loftee_os_enters_splice_pp3_and_blocks_bp4() -> None:
    frame = pd.DataFrame(
        {
            "Consequence": [
                "splice_donor_5th_base_variant&intron_variant",
                "splice_donor_5th_base_variant&intron_variant",
            ],
            "LoF": ["OS", "LC"],
            "PrimateAI": [0.0, 0.0],
            "am_pathogenicity": [0.0, 0.0],
            "CADD_phred": [10.0, 10.0],
            "CADD_reg_phred": [10.0, 10.0],
            "am_class": ["", ""],
            "vep_consq_lof": [True, False],
            "splicing_lof": [False, False],
            "5UTR_lof": [False, False],
            "CLNSIG": ["", ""],
            "CLNREVSTAT": ["", ""],
        }
    )

    pp3, bp4 = PP3_BP4_criteria(frame)

    assert pp3.tolist() == [1, 0]
    assert bp4.tolist() == [0, 1]


def test_variant_level_mechanism_contract(tmp_path: Path) -> None:
    sources = _write_fixture_sources(tmp_path)
    frame = pd.DataFrame(
        [
            {
                "SYMBOL": "TESTREC",
                "Gene": "ENSG1",
                "HPO_gene_inheritance": "Autosomal recessive inheritance",
                "Consequence": "missense_variant",
                "LoF": "",
                "NMD": "",
                "vep_consq_lof": False,
                "variant_gof_tag": "",
            },
            {
                "SYMBOL": "TESTLOF",
                "Gene": "ENSG2",
                "HPO_gene_inheritance": "Autosomal recessive inheritance",
                "Consequence": "splice_donor_5th_base_variant&intron_variant",
                "LoF": "OS",
                "NMD": "",
                "vep_consq_lof": True,
                "variant_gof_tag": "",
            },
            {
                "SYMBOL": "TESTGOF",
                "Gene": "ENSG3",
                "HPO_gene_inheritance": "Autosomal recessive inheritance",
                "Consequence": "missense_variant",
                "LoF": "",
                "NMD": "",
                "vep_consq_lof": False,
                "variant_gof_tag": "GOF",
                "gofcards_variant_id": "SNV|1|100|100|A|G",
                "gofcards_accession_id": "rs1",
            },
            {
                "SYMBOL": "TESTMONO",
                "Gene": "ENSG4",
                "HPO_gene_inheritance": "Autosomal dominant inheritance",
                "Consequence": "stop_gained",
                "LoF": "HC",
                "NMD": "NMD",
                "vep_consq_lof": True,
                "variant_gof_tag": "",
            },
        ]
    )

    annotated = annotate_gene_mechanism_categories(
        frame,
        clinvar_pathogenic_genes=set(),
        gene_to_am_score_map={},
        use_hgnc_package=False,
        **sources,
    ).set_index("SYMBOL")

    hpo_only = annotated.loc["TESTREC"]
    assert "gene_mech_inher_history" not in annotated.columns
    assert hpo_only["var_plausible_patho_mechs"] == ""
    assert hpo_only["variant_mechanism_applicable"] == ""

    biallelic_lof = annotated.loc["TESTLOF"]
    assert biallelic_lof["var_plausible_patho_mechs"] == "recessive_LOF"
    assert biallelic_lof["variant_effect"] == "predicted_LOF_high_confidence"
    assert "LOFTEE_OS" in biallelic_lof["variant_effect_evidence"]
    assert "recessive_LOF" in biallelic_lof["variant_mechanism_applicable"]

    biallelic_gof = annotated.loc["TESTGOF"]
    assert biallelic_gof["var_plausible_patho_mechs"] == "recessive_GOF"
    assert biallelic_gof["variant_effect"] == "exact_known_GOF"
    assert "recessive_GOF" in biallelic_gof["variant_mechanism_applicable"]
    assert biallelic_gof["clinvar_vcv_accessions"] == "VCV000000001"
    assert biallelic_gof["clinvar_rcv_conditions"] == "biallelic GOF disorder"
    assert biallelic_gof["clinvar_vcv_max_review_stars"] == "2"
    assert "NM_000001.1:c.1A>G" in biallelic_gof["clinvar_vcv_hgvs"]

    monoallelic_lof = annotated.loc["TESTMONO"]
    assert monoallelic_lof["var_plausible_patho_mechs"] == "dominant_LOF"
    assert "dominant_LOF" in monoallelic_lof["variant_mechanism_applicable"]

    masks = _variant_mechanism_masks(annotated.reset_index())
    assert masks["has_recessive_compatible"].tolist() == [False, True, True, False]
    assert masks["has_dominant_compatible"].tolist() == [False, False, False, True]
    assert masks["has_applicable_lof_assertion"].tolist() == [False, True, False, True]


def test_gene_wide_lof_signals_remain_audit_only_without_condition_history(
    tmp_path: Path,
) -> None:
    sources = _write_fixture_sources(tmp_path)
    frame = pd.DataFrame(
        [
            {
                "SYMBOL": "TESTREC",
                "Gene": "GENE_NONE",
                "HPO_gene_inheritance": "Autosomal recessive inheritance",
                "LOEUF": 1.0,
                "Consequence": "stop_gained",
                "LoF": "HC",
                "NMD": "NMD",
                "vep_consq_lof": True,
                "variant_gof_tag": "",
            },
            {
                "SYMBOL": "TESTREC",
                "Gene": "GENE_CLINVAR",
                "HPO_gene_inheritance": "Autosomal recessive inheritance",
                "LOEUF": 1.0,
                "Consequence": "stop_gained",
                "LoF": "HC",
                "NMD": "NMD",
                "vep_consq_lof": True,
                "variant_gof_tag": "",
            },
            {
                "SYMBOL": "TESTREC",
                "Gene": "GENE_LOEUF",
                "HPO_gene_inheritance": "Autosomal recessive inheritance",
                "LOEUF": 0.349,
                "Consequence": "stop_gained",
                "LoF": "HC",
                "NMD": "NMD",
                "vep_consq_lof": True,
                "variant_gof_tag": "",
            },
            {
                "SYMBOL": "TESTREC",
                "Gene": "GENE_AM",
                "HPO_gene_inheritance": "Autosomal recessive inheritance",
                "LOEUF": 1.0,
                "Consequence": "stop_gained",
                "LoF": "HC",
                "NMD": "NMD",
                "vep_consq_lof": True,
                "variant_gof_tag": "",
            },
            {
                "SYMBOL": "TESTREC",
                "Gene": "GENE_LOEUF_BOUNDARY",
                "HPO_gene_inheritance": "Autosomal recessive inheritance",
                "LOEUF": 0.35,
                "Consequence": "stop_gained",
                "LoF": "HC",
                "NMD": "NMD",
                "vep_consq_lof": True,
                "variant_gof_tag": "",
            },
            {
                "SYMBOL": "TESTREC",
                "Gene": "GENE_AM_BOUNDARY",
                "HPO_gene_inheritance": "Autosomal recessive inheritance",
                "LOEUF": 1.0,
                "Consequence": "stop_gained",
                "LoF": "HC",
                "NMD": "NMD",
                "vep_consq_lof": True,
                "variant_gof_tag": "",
            },
        ]
    )

    annotated = annotate_gene_mechanism_categories(
        frame,
        clinvar_pathogenic_genes={"GENE_CLINVAR"},
        gene_to_am_score_map={"GENE_AM": 0.565, "GENE_AM_BOUNDARY": 0.564},
        use_hgnc_package=False,
        **sources,
    ).set_index("Gene")

    for gene in annotated.index:
        assert annotated.loc[gene, "var_plausible_patho_mechs"] == ""
        assert annotated.loc[gene, "variant_mechanism_applicable"] == ""
    assert annotated.loc["GENE_CLINVAR", "gene_lof_evidence"] == "ClinVar_pathogenic_2plus"
    assert annotated.loc["GENE_LOEUF", "gene_lof_evidence"] == "LOEUF_lt_0.35"
    assert annotated.loc["GENE_AM", "gene_lof_evidence"] == "GeneAvgAM_gt_0.564"


def test_acmg_masks_reject_legacy_gene_level_fallback() -> None:
    legacy = pd.DataFrame(
        {
            "gene_mech_inher_history": ["LoF_recessive"],
            "Consequence": ["stop_gained"],
            "LoF": ["HC"],
            "variant_gof_tag": [""],
        }
    )

    with pytest.raises(KeyError, match="variant-level mechanism annotations are required"):
        _variant_mechanism_masks(legacy)


def test_clinvar_gene_lof_gate_requires_pathogenic_and_two_stars() -> None:
    clinvar = {
        "GENE_PATH": {
            1: {
                "p.Ala1Val": {
                    "CLNSIG": ["Likely_pathogenic"],
                    "CLNREVSTAT": ["criteria_provided,_multiple_submitters,_no_conflicts"],
                }
            }
        },
        "GENE_ONE_STAR": {
            1: {
                "p.Ala1Val": {
                    "CLNSIG": ["Pathogenic"],
                    "CLNREVSTAT": ["criteria_provided,_single_submitter"],
                }
            }
        },
        "GENE_CONFLICT": {
            1: {
                "p.Ala1Val": {
                    "CLNSIG": ["Conflicting_classifications_of_pathogenicity"],
                    "CLNREVSTAT": ["reviewed_by_expert_panel"],
                }
            }
        },
    }

    assert summarize_clinvar_gene_pathogenicity(clinvar) == {"GENE_PATH"}


def _modern_mechanism_rows() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "HPO_IDs": ["", "", ""],
            "var_plausible_patho_mechs": [
                "recessive_GOF",
                "recessive_GOF",
                "dominant_GOF",
            ],
            "variant_mechanism_applicable": [
                "recessive_GOF",
                "recessive_GOF",
                "dominant_GOF",
            ],
            "variant_mechanism_uncertain": ["", "", ""],
            "variant_effect": ["exact_known_GOF"] * 3,
            "variant_gof_tag": ["GOF"] * 3,
        }
    )


def test_bs1_and_bs2_use_variant_level_allelic_requirement() -> None:
    frame = _modern_mechanism_rows().assign(
        gnomAD_joint_AF_max=[0.1, 0.1, 0.1],
        gnomAD_joint_AN_max=[200, 200, 200],
        gnomAD_nhomalt_max=[0, 12, 0],
        gnomAD_nhomalt_XX=[0, 12, 0],
        gnomAD_nhomalt_XY=[0, 0, 0],
        gnomAD_joint_AF=[0.1, 0.1, 0.1],
        gnomAD_joint_AN=[200, 200, 200],
        gnomAD_joint_AF_XY=[0.0, 0.0, 0.0],
        clinvar_patho_gene_max_af=[0.0, 0.0, 0.0],
    )
    false = np.zeros(len(frame), dtype=bool)

    bs1 = BS1_criteria(
        frame,
        pm2_criteria=np.zeros(len(frame), dtype=int),
        non_monogenic=false,
        non_mendelian=false,
        incomplete_penetrance=false,
    )
    assert bs1.tolist() == [0, 3, 3]

    bs2 = BS2_criteria(
        frame,
        non_monogenic=false,
        non_mendelian=false,
        incomplete_penetrance=false,
        pm2_criteria=np.zeros(len(frame), dtype=int),
    )
    assert bs2.tolist() == [0, 3, 3]


def test_bs4_uses_biallelic_gof_as_recessive_requirement() -> None:
    frame = _modern_mechanism_rows().iloc[[0, 2]].copy()
    frame["PROBAND"] = "0/1:20"
    frame["CONTROL"] = "0/1:20"
    pedigree = pd.DataFrame(
        {
            "#FamilyID": ["F1", "F1"],
            "IndividualID": ["PROBAND", "CONTROL"],
            "Phenotype": ["2", "1"],
            "Sex": ["2", "2"],
        }
    )
    false = np.zeros(len(frame), dtype=bool)

    bs4 = BS4_criteria(
        frame,
        pedigree,
        "F1",
        non_monogenic=false,
        non_mendelian=false,
        incomplete_penetrance=false,
    )

    assert bs4.tolist() == [0, 1]
