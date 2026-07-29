from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gene_mechanism_hub import GeneMechanismHub  # noqa: E402
from acmg_criteria_assign import annotate_exact_nonlof_variants  # noqa: E402


def _exact_row(
    symbol: str,
    hgvsp: str,
    mechanism: str,
    allele_key: str,
    *,
    eligibility: str = "eligible",
    chrom: str = "10",
    hg19_pos: int = 100,
    hg38_pos: int = 200,
) -> dict[str, Any]:
    """One variant in the nested shape the canonical non-LOF JSON embeds.

    ``hgvsp`` arrives fully qualified (``ENSP...:p.Pro253Arg``); the cache
    stores the bare protein change, which is what the matcher normalizes and
    indexes on.
    """
    protein = hgvsp.split(":")[-1]
    coding = "c.758C>G"
    variant_id = (f"loc_{chrom}:{hg19_pos}:A->G_grch37"
                  if eligibility == "eligible"
                  else f"loc_{chrom}:{hg19_pos + 1}:A->G_grch37")

    def _assembly(pos: int, transcript: str, status: str) -> dict[str, Any]:
        return {
            "genomic": {"chrom": chrom, "pos": pos, "ref": "A", "alt": "G",
                        "status": status},
            "transcripts": {transcript: {
                "by_hgvsc": {coding: {"hgvsp": protein,
                                      "consequence": "missense_variant",
                                      "canonical": True, "mane_select": None}},
                "by_hgvsp": {protein: [coding]},
            }},
        }

    return {
        "source": "GoFCards" if mechanism == "GOF" else "CuratedDN",
        "mechanism": mechanism,
        "HGNC_Symbol": symbol,
        "variant_id": variant_id,
        "match_eligibility": eligibility,
        "record": {
            "source": {"gofcards_allele_key": allele_key, "variant_type_label": "SNV",
                       "assembly": "hg19", "chrom": chrom, "start": str(hg19_pos),
                       "ref": "A", "alt": "G"},
            "eligibility": {
                "status": eligibility,
                "gene_match_status": ("gene_concordant" if eligibility == "eligible"
                                      else "gene_discordant"),
                "vep_symbol": symbol,
                "reason": None if eligibility == "eligible" else "gene_discordant",
            },
            "liftover_status": "mapped",
            "annotations": {"clinvar_variation_id": f"ID-{allele_key}"},
            "evidence": [{"pmid": "1", "disease": "Test condition"}],
        },
        "assemblies": {
            "hg19": _assembly(hg19_pos, "ENST00000000001.1", "raw_ref_alt"),
            "hg38": _assembly(hg38_pos, "ENST00000000001.2", "lifted_ref_match"),
        },
    }


def _write_mechanism_cache(tmp_path: Path) -> Path:
    fgfr2 = _exact_row(
        "FGFR2",
        "ENSP00000338548.7:p.Pro253Arg",
        "GOF",
        "SNV|10|100|100|A|G",
    )
    quarantined = _exact_row(
        "FGFR2",
        "ENSP00000338548.7:p.Ala1Val",
        "GOF",
        "SNV|10|101|101|A|G",
        eligibility="quarantined_gene_discordance",
    )
    dual_gof = _exact_row(
        "DUAL1",
        "ENSP_DUAL.1:p.Arg10His",
        "GOF",
        "SNV|1|10|10|G|A",
    )
    dual_dn = _exact_row(
        "DUAL1",
        "ENSP_DUAL.1:p.Arg10His",
        "DOMINANT_NEGATIVE",
        "SNV|1|10|10|G|A",
    )
    payload = {
        "_meta": {"version": "2.0"},
        "HGNC:3689": {
            "symbol": "FGFR2",
            "mechanisms": ["GOF"],
            "variant_level": [
                {
                    "GoFCards": {
                        "mechanism": "GOF",
                        "exact_normalization_status": "matched_gene_concordant",
                        "exact_normalized_variants": [fgfr2],
                    }
                },
                {
                    "GoFCards": {
                        "mechanism": "GOF",
                        "exact_normalization_status": (
                            "quarantined_upstream_gene_discordance"
                        ),
                        "quarantined_exact_normalized_variants": [quarantined],
                    }
                },
            ],
        },
        "HGNC:99998": {
            "symbol": "DUAL1",
            "mechanisms": ["GOF", "DOMINANT_NEGATIVE"],
            "variant_level": [
                {
                    "CuratedGOF": {
                        "mechanism": "GOF",
                        "exact_normalized_variants": [dual_gof],
                    }
                },
                {
                    "CuratedDN": {
                        "mechanism": "DOMINANT_NEGATIVE",
                        "exact_normalized_variants": [dual_dn],
                    }
                },
            ],
        },
    }
    path = tmp_path / "canonical-nonlof.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_hgnc_gene_plus_hgvsp_is_an_exact_exclusive_match(tmp_path: Path) -> None:
    mechanism_json = _write_mechanism_cache(tmp_path)
    hub = GeneMechanismHub(
        mechanism_json=mechanism_json,
        use_hgnc_package=False,
    )

    match = hub.match_curated_nonlof_variant("FGFR2", hgvsp="p.Pro253Arg")
    wrong_gene = hub.match_curated_nonlof_variant("INPP5F", hgvsp="p.Pro253Arg")
    quarantined = hub.match_curated_nonlof_variant("FGFR2", hgvsp="p.Ala1Val")

    assert match["match_route"] == "hgvsp"
    assert match["mechanism_scores"] == {
        "LOF": 0,
        "GOF": 2,
        "DOMINANT_NEGATIVE": 0,
    }
    assert match["mechanism_exclusive"] is True
    assert match["exclusive_mechanisms"] == ["GOF"]
    assert wrong_gene["matches"] == []
    assert quarantined["matches"] == []


def test_compatibility_matcher_ignores_the_compact_tsv_argument(tmp_path: Path) -> None:
    mechanism_json = _write_mechanism_cache(tmp_path)
    hub = GeneMechanismHub(
        mechanism_json=mechanism_json,
        use_hgnc_package=False,
    )

    match = hub.match_gofcards_variant_gof(
        "FGFR2",
        "ENSP_OTHER:p.Pro253Arg",
        gofcards_exact_hgvsp_tsv=tmp_path / "does-not-exist.tsv.gz",
    )

    assert match["variant_gof_tag"] == "GOF"
    assert match["gof_score"] == 2
    assert match["matches"][0]["canonical_mechanism_json"] == str(mechanism_json)


def test_normalized_genomic_allele_is_an_exact_fallback(tmp_path: Path) -> None:
    mechanism_json = _write_mechanism_cache(tmp_path)
    hub = GeneMechanismHub(
        mechanism_json=mechanism_json,
        use_hgnc_package=False,
    )

    match = hub.match_curated_nonlof_variant(
        "FGFR2",
        chrom="chr10",
        pos=200,
        ref="A",
        alt="G",
        assembly="GRCh38",
    )

    assert match["match_route"] == "genomic"
    assert match["gof_score"] == 2
    assert match["matched_key_type"] in {"vcf", "genomic"}


def test_dual_exact_scores_require_two_explicit_allele_assertions(tmp_path: Path) -> None:
    mechanism_json = _write_mechanism_cache(tmp_path)
    hub = GeneMechanismHub(
        mechanism_json=mechanism_json,
        use_hgnc_package=False,
    )

    match = hub.match_curated_nonlof_variant("DUAL1", hgvsp="p.Arg10His")

    assert match["mechanism_scores"] == {
        "LOF": 0,
        "GOF": 2,
        "DOMINANT_NEGATIVE": 2,
    }
    assert match["exclusive_mechanisms"] == ["DOMINANT_NEGATIVE", "GOF"]
    assert {record["canonical_assertion_source"] for record in match["matches"]} == {
        "CuratedDN",
        "CuratedGOF",
    }


def test_dataframe_annotator_writes_exclusive_scores_and_audit(tmp_path: Path) -> None:
    mechanism_json = _write_mechanism_cache(tmp_path)
    frame = pd.DataFrame(
        [
            {
                "SYMBOL": "FGFR2",
                "HGVSp": "p.Pro253Arg",
                "chrom": "chr10",
                "pos": "200",
                "ref": "A",
                "alt": "G",
                "assembly": "GRCh38",
            },
            {
                "SYMBOL": "DUAL1",
                "HGVSp": "p.Arg10His",
                "chrom": "chr1",
                "pos": "10",
                "ref": "G",
                "alt": "A",
                "assembly": "GRCh38",
            },
        ]
    )

    annotated = annotate_exact_nonlof_variants(
        frame,
        mechanism_json=mechanism_json,
        context="unit_test",
    )

    assert annotated["variant_lof_score"].tolist() == [0, 0]
    assert annotated["variant_gof_score"].tolist() == [2, 2]
    assert annotated["variant_dn_score"].tolist() == [0, 2]
    assert annotated["variant_mechanism_exclusive"].tolist() == [True, True]
    assert annotated["variant_exact_mechanisms"].tolist() == [
        "GOF",
        "DOMINANT_NEGATIVE;GOF",
    ]
    assert annotated["variant_mechanism_match_route"].tolist() == [
        "hgvsp",
        "hgvsp",
    ]
    assert "canonical_mechanism_json" in annotated.loc[
        0, "variant_exact_mechanism_evidence"
    ]
