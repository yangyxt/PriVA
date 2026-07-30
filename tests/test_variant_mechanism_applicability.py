import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from acmg_bs2_bp5_observation import BS2_criteria  # noqa: E402
from acmg_consequence import vep_consq_interpret_per_row  # noqa: E402
from acmg_pm2_bs1_ba1_frequency import BS1_criteria  # noqa: E402
from acmg_pp1_bs4_bp2_pm3_family import BS4_criteria  # noqa: E402
from acmg_pp3_bp4_bp7_insilico import PP3_BP4_criteria  # noqa: E402
from acmg_pvs1_null_variant import (  # noqa: E402
    summarize_clinvar_gene_pathogenicity,
)
from acmg_variant_mechanism import _variant_mechanism_masks  # noqa: E402
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
                                            "gofcards_variant_id": "loc_1:100:A->G_grch37",
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

    def cache_condition(
        condition_id: str,
        label: str,
        mechanism: str,
        allelic_requirement: str,
        inheritance_mode: str,
        source_record_id: str,
    ) -> dict:
        hpo_id = (
            "HP:0000007"
            if inheritance_mode == "autosomal_recessive"
            else "HP:0000006"
        )
        return {
            "label": label,
            "identifiers": {"OMIM": [condition_id]},
            "priva_scope": {
                "decision": "include",
                "category": "mendelian_non_neoplastic",
                "review_status": "auto_supported",
            },
            "inheritance": {
                "modes": [inheritance_mode],
                "assertions": [
                    {
                        "hpo_id": hpo_id,
                        "frequency": "-",
                        "evidence": "TAS",
                        "reference": condition_id,
                    }
                ],
            },
            "penetrance": {"statuses": [], "assertions": []},
            "onset": {"terms": [], "assertions": []},
            "hpo_assertion_count": 1,
            "pathogenic_mechanisms": {
                mechanism: {
                    "allelic_requirements": [allelic_requirement],
                    "evidence": [
                        {
                            "source": "G2P_DDG2P",
                            "source_record_id": source_record_id,
                            "condition_identifiers": [condition_id],
                            "condition_label": label,
                            "mechanism": mechanism,
                            "mechanism_raw": (
                                "gain of function"
                                if mechanism == "GOF"
                                else "loss of function"
                            ),
                            "allelic_requirement": allelic_requirement,
                            "mechanism_confidence": "high",
                            "disease_confidence": "definitive",
                            "pmids": [source_record_id],
                        }
                    ],
                    "variants": {},
                }
            },
        }

    lof_condition = cache_condition(
        "OMIM:1",
        "biallelic LOF disorder",
        "LOF",
        "biallelic_autosomal",
        "autosomal_recessive",
        "1",
    )
    mono_condition = cache_condition(
        "OMIM:2",
        "monoallelic LOF disorder",
        "LOF",
        "monoallelic_autosomal",
        "autosomal_dominant",
        "2",
    )
    gof_condition = cache_condition(
        "OMIM:3",
        "biallelic GOF disorder",
        "GOF",
        "biallelic_autosomal",
        "autosomal_recessive",
        "3",
    )
    gof_condition["pathogenic_mechanisms"]["GOF"]["variants"] = {
        "GOFCARDS:loc_1:100:A->G_grch37": {
            "mechanism": "GOF",
            "symbol": "GENE1",
            "gofcards_variant_id": "loc_1:100:A->G_grch37",
            # Transcript with its version, and the HGVS notations, per view.
            "transcripts": [
                {
                    "assembly": "hg19",
                    "transcript": "ENST00000000001.1",
                    "hgvsc": "NM_000001.1:c.1A>G",
                    "hgvsp": "NP_000001.1:p.Met1Val",
                }
            ],
            "clinvar_links": [
                {
                    "vcv_accession": "VCV000000001",
                    "rcv_accession": "RCV000000001",
                    "condition_identifiers": ["OMIM:3"],
                    "condition_names": ["biallelic GOF disorder"],
                    "hgvs": [
                        "NM_000001.1:c.1A>G",
                        "NP_000001.1:p.Met1Val",
                    ],
                    "clinical_significance": "Pathogenic",
                    "review_stars": 2,
                }
            ],
            "condition_link": {"status": "exact", "condition_key": "OMIM:3"},
        }
    }
    condition_cache = tmp_path / "condition-cache.json"
    condition_cache.write_text(
        json.dumps(
            {
                "_meta": {"schema_version": "1.0"},
                "genes": {
                    "TESTLOF": {
                        "conditions": {"OMIM:1": lof_condition},
                        "summary": {"pathogenic_mechanisms": ["LOF"]},
                        "unmapped_evidence": {"mechanisms": [], "variants": {}},
                    },
                    "TESTMONO": {
                        "conditions": {"OMIM:2": mono_condition},
                        "summary": {"pathogenic_mechanisms": ["LOF"]},
                        "unmapped_evidence": {"mechanisms": [], "variants": {}},
                    },
                    "TESTGOF": {
                        "conditions": {"OMIM:3": gof_condition},
                        "summary": {"pathogenic_mechanisms": ["GOF"]},
                        "unmapped_evidence": {"mechanisms": [], "variants": {}},
                    },
                },
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
        )
        + "\n"
        + "TESTLOF\tG2P_DDG2P\t1\tOMIM:1\tMONDO:1\t"
        "mendelian_non_neoplastic\tinclude\tauto_supported\t"
        "biallelic LOF disorder\t"
        "biallelic_autosomal\tloss of function\tLOF\thigh\tdefinitive\t1\n"
        + "TESTMONO\tG2P_DDG2P\t2\tOMIM:2\tMONDO:2\t"
        "mendelian_non_neoplastic\tinclude\tauto_supported\t"
        "monoallelic LOF disorder\t"
        "monoallelic_autosomal\tloss of function\tLOF\thigh\tdefinitive\t2\n"
        + "TESTGOF\tG2P_DDG2P\t3\tOMIM:3\tMONDO:3\t"
        "mendelian_non_neoplastic\tinclude\tauto_supported\t"
        "biallelic GOF disorder\t"
        "biallelic_autosomal\tgain of function\tGOF\thigh\tdefinitive\t3\n",
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
    # mechanism_json and ddg2p_evidence are no longer inputs: the condition
    # cache is the only evidence source the chain reads.
    return {
        "condition_cache": condition_cache,
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
                "gofcards_variant_id": "loc_1:100:A->G_grch37",
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
        use_hgnc_package=False,
        **sources,
    ).set_index("SYMBOL")

    hpo_only = annotated.loc["TESTREC"]
    assert "gene_mech_inher_history" not in annotated.columns
    assert hpo_only["var_plausible_patho_mechs"] == ""
    assert hpo_only["variant_mechanism_applicable"] == ""

    # A fifth-base splice-donor variant with LOFTEE "OS" is plausible loss of
    # function, not established: nothing says the transcript is destroyed, so
    # it scores 1 and its history is POSSIBLE rather than applicable.
    biallelic_lof = annotated.loc["TESTLOF"]
    assert biallelic_lof["var_plausible_patho_mechs"] == "recessive_LOF"
    assert biallelic_lof["variant_effect"] == "predicted_LOF_high_confidence"
    assert biallelic_lof["variant_lof_score"] == 1
    assert biallelic_lof["variant_mechanism_applicable"] == ""
    assert "recessive_LOF" in biallelic_lof["variant_mechanism_uncertain"]
    # The three facts the chain delivers, at this variant's own resolution.
    assert biallelic_lof["variant_condition_ids"] == "OMIM:1"
    assert biallelic_lof["variant_inheritance"] == "recessive"
    assert biallelic_lof["variant_x_linked"] == "false"
    assert biallelic_lof["variant_penetrance"] == "unknown"

    biallelic_gof = annotated.loc["TESTGOF"]
    assert biallelic_gof["var_plausible_patho_mechs"] == "recessive_GOF"
    assert biallelic_gof["variant_effect"] == "exact_known_GOF"
    assert "recessive_GOF" in biallelic_gof["variant_mechanism_applicable"]
    assert biallelic_gof["variant_condition_ids"] == "OMIM:3"
    assert biallelic_gof["variant_inheritance"] == "recessive"

    monoallelic_lof = annotated.loc["TESTMONO"]
    assert monoallelic_lof["var_plausible_patho_mechs"] == "dominant_LOF"
    assert "dominant_LOF" in monoallelic_lof["variant_mechanism_applicable"]

    masks = _variant_mechanism_masks(annotated.reset_index())
    assert masks["has_recessive_compatible"].tolist() == [False, True, True, False]
    assert masks["has_dominant_compatible"].tolist() == [False, False, False, True]
    # The masks answer two separate questions and never conjoin them. What the
    # gene's history says stays in the *_history masks; what this variant does
    # stays in the is_* masks. TESTMONO is a nonsense variant triggering decay,
    # so it reaches loss of function at full grade; TESTLOF escapes decay with
    # only LOFTEE "OS" behind it, so it reaches loss of function at the lower
    # grade. Both are is_predicted_lof; variant_lof_score separates them.
    assert masks["is_predicted_lof"].tolist() == [False, True, False, True]
    assert masks["has_rec_lof_history"].tolist() == [False, True, False, False]
    assert masks["has_dom_lof_history"].tolist() == [False, False, False, True]
    assert "has_established_lof_mechanism" not in masks
    assert "has_lof_mechanism_history" not in masks


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
        use_hgnc_package=False,
        **sources,
    ).set_index("Gene")

    # Gene-wide loss-of-function signals -- a ClinVar pathogenic history, a
    # constrained LOEUF, a high average AlphaMissense score -- say nothing about
    # which condition a variant acts on or by what mechanism. None of them
    # creates a mechanism, a plausible-mechanism tag, or a condition.
    for gene in annotated.index:
        assert annotated.loc[gene, "var_plausible_patho_mechs"] == ""
        assert annotated.loc[gene, "variant_mechanism_applicable"] == ""
        assert annotated.loc[gene, "variant_condition_ids"] == ""

    # Inheritance is the one exception, and only as a last resort. These genes
    # have no condition stating an inheritance, so rather than deliver nothing
    # the constraint data decides how many copies must be affected. The basis
    # column says so, and no mechanism accompanies it.
    for gene in annotated.index:
        assert annotated.loc[gene, "variant_inheritance"] in {"recessive", "dominant"}
        assert annotated.loc[gene, "variant_inheritance_basis"] == "gene_constraint"


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
    # Both ClinVar structures are keyed by transcript, so a transcript-to-gene
    # map is what makes the answer gene-level.
    tx2gene = {
        "ENST_PATH": "GENE_PATH",
        "ENST_ONE_STAR": "GENE_ONE_STAR",
        "ENST_CONFLICT": "GENE_CONFLICT",
        "ENST_SPLICE": "GENE_SPLICE_ONLY",
        "ENST_SPLICE_ONE_STAR": "GENE_SPLICE_ONE_STAR",
    }
    aa_dict = {
        "ENST_PATH": {
            1: {
                "p.Ala1Val": {
                    "CLNSIG": ["Likely_pathogenic"],
                    "CLNREVSTAT": ["criteria_provided,_multiple_submitters,_no_conflicts"],
                }
            }
        },
        "ENST_ONE_STAR": {
            1: {
                "p.Ala1Val": {
                    "CLNSIG": ["Pathogenic"],
                    "CLNREVSTAT": ["criteria_provided,_single_submitter"],
                }
            }
        },
        "ENST_CONFLICT": {
            1: {
                "p.Ala1Val": {
                    "CLNSIG": ["Conflicting_classifications_of_pathogenicity"],
                    "CLNREVSTAT": ["reviewed_by_expert_panel"],
                }
            }
        },
        # A transcript nothing maps to contributes nothing.
        "ENST_UNMAPPED": {
            1: {
                "p.Ala1Val": {
                    "CLNSIG": ["Pathogenic"],
                    "CLNREVSTAT": ["practice_guideline"],
                }
            }
        },
    }
    # A canonical splice-site variant carries no HGVSp, so it can only ever
    # reach the gate through this second structure.
    splice_dict = {
        "ENST_SPLICE": [
            {
                "hgvsc": "c.1234+1G>A",
                "consequence": "splice_donor_variant",
                "clinvar_sig": "Pathogenic",
                "clinvar_review": "reviewed_by_expert_panel",
            }
        ],
        "ENST_SPLICE_ONE_STAR": [
            {
                "hgvsc": "c.1234+1G>A",
                "consequence": "splice_donor_variant",
                "clinvar_sig": "Pathogenic",
                "clinvar_review": "criteria_provided,_single_submitter",
            }
        ],
    }

    # The amino-acid structure alone cannot see the splice-only gene.
    assert summarize_clinvar_gene_pathogenicity(
        tx2gene, clinvar_aa_dict=aa_dict
    ) == {"GENE_PATH"}

    # Adding the second structure recovers it, and the two-star bar and the
    # conflicting-classification exclusion still hold across both.
    assert summarize_clinvar_gene_pathogenicity(
        tx2gene, clinvar_aa_dict=aa_dict, clinvar_splice_dict=splice_dict
    ) == {"GENE_PATH", "GENE_SPLICE_ONLY"}


def test_scope_review_blocks_automatic_bs2() -> None:
    frame = pd.DataFrame(
        [
            {
                "chrom": "chr1",
                "HPO_IDs": "",
                "HPO_scope_review_required": 0,
                "var_plausible_patho_mechs": "dominant",
                "variant_effect": "uncertain",
                "variant_mechanism_applicable": "",
                "variant_mechanism_uncertain": "dominant",
                "Consequence": "missense_variant",
                "NMD": "",
                "LoF": "",
                "LoF_filter": "",
                "variant_gof_tag": "",
                "gnomAD_joint_AF": 0.01,
                "gnomAD_joint_AN": 2000,
                "gnomAD_nhomalt_XX": 0,
                "gnomAD_nhomalt_XY": 0,
            },
            {
                "chrom": "chr1",
                "HPO_IDs": "",
                "HPO_scope_review_required": 1,
                "var_plausible_patho_mechs": "dominant",
                "variant_effect": "uncertain",
                "variant_mechanism_applicable": "",
                "variant_mechanism_uncertain": "dominant",
                "Consequence": "missense_variant",
                "NMD": "",
                "LoF": "",
                "LoF_filter": "",
                "variant_gof_tag": "",
                "gnomAD_joint_AF": 0.01,
                "gnomAD_joint_AN": 2000,
                "gnomAD_nhomalt_XX": 0,
                "gnomAD_nhomalt_XY": 0,
            },
        ]
    )
    false = np.array([False, False])
    result = BS2_criteria(
        frame,
        non_monogenic=false,
        non_mendelian=false,
        incomplete_penetrance=false,
        pm2_criteria=np.array([0, 0]),
    )

    assert result.tolist() == [3, 0]


def test_gene_hub_does_not_reintroduce_excluded_hpo_inheritance(
    tmp_path: Path,
) -> None:
    sources = _write_fixture_sources(tmp_path)
    pd.DataFrame(
        [
            {
                "gene_symbol": "TESTREC",
                "disease_id": "OMIM:INCLUDED",
                "hpo_id": "HP:0000006",
                "frequency": "-",
                "evidence": "TAS",
                "reference": "OMIM:INCLUDED",
                "priva_scope": "include",
            },
            {
                "gene_symbol": "TESTREC",
                "disease_id": "OMIM:EXCLUDED",
                "hpo_id": "HP:0000007",
                "frequency": "-",
                "evidence": "TAS",
                "reference": "OMIM:EXCLUDED",
                "priva_scope": "exclude",
            },
            {
                "gene_symbol": "TESTREC",
                "disease_id": "OMIM:REVIEW",
                "hpo_id": "HP:0000007",
                "frequency": "-",
                "evidence": "TAS",
                "reference": "OMIM:REVIEW",
                "priva_scope": "review",
            },
        ]
    ).to_csv(sources["hpo_collapsed"], sep="\t", index=False)

    summary = GeneMechanismHub(
        use_hgnc_package=False,
        **sources,
    ).known_inheritance_mode("TESTREC")

    assert summary["dominant"] is True
    assert summary["recessive"] is False
    assert summary["hpo_scope_review_required"] is True
    assert summary["hpo_scope_review_disease_ids"] == "OMIM:REVIEW"
    assert summary["hpo_scope_excluded_disease_ids"] == "OMIM:EXCLUDED"


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
