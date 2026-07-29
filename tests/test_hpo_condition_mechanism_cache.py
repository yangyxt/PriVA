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
from clinvar_vcv import partition_gofcards_variants  # noqa: E402


HPO_HEADER = (
    "gene_symbol\tdisease_id\thpo_id\tfrequency\tevidence\treference\t"
    "mondo_id\tmondo_name\tdisease_scope\tpriva_scope\tscope_evidence\t"
    "scope_reference\tscope_review_status\n"
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
    gofcards = _write_gofcards_cache(tmp_path / "gofcards.json", {})
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
