from __future__ import annotations

import json
import hashlib
import sys
from pathlib import Path

import pandas as pd
import pytest
from jsonschema import Draft202012Validator
from jsonschema.exceptions import ValidationError


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from clinvar_vcv import (  # noqa: E402
    load_exact_gofcards_lookup,
    review_stars,
    stream_parse_clinvar_vcv,
)
from build_gene_nonlof_mechanism_cache import (  # noqa: E402
    GOFCARDS_EXACT_COLUMNS,
    _latest_clinvar_xsd_url,
    _parse_remote_md5,
    build_nonlof_assertions_json,
    fetch_clinvar_vcv,
    inject_clinvar_matches,
    load_gofcards_exact_records,
    load_hgnc_mapping,
    make_unified_row,
    validate_canonical_json,
)
from audit_clinvar_gofcards_review_tiers import summarize_review_tiers  # noqa: E402


def _simple_allele(variation_id: int, pos: int, symbol: str = "MEFV") -> str:
    return f"""
      <SimpleAllele AlleleID="{variation_id + 100}" VariationID="{variation_id}">
        <GeneList><Gene Symbol="{symbol}" HGNC_ID="HGNC:6998" GeneID="4210"/></GeneList>
        <Name>NM_TEST.1({symbol}):c.1A&gt;G</Name>
        <CanonicalSPDI>NC_TEST:{pos - 1}:A:G</CanonicalSPDI>
        <VariantType>single nucleotide variant</VariantType>
        <Location>
          <SequenceLocation Assembly="GRCh37" Chr="16" positionVCF="{pos}" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
          <SequenceLocation Assembly="GRCh38" Chr="16" positionVCF="{pos + 100}" referenceAlleleVCF="A" alternateAlleleVCF="G"/>
        </Location>
        <HGVSlist><HGVS Type="coding"><NucleotideExpression sequenceAccessionVersion="NM_TEST.1"><Expression>NM_TEST.1:c.1A&gt;G</Expression></NucleotideExpression></HGVS></HGVSlist>
      </SimpleAllele>
    """


def _rcv(accession: str, cui: str, status: str, significance: str = "Pathogenic") -> str:
    return f"""
      <RCVAccession Accession="{accession}" Version="2" Title="test condition">
        <ClassifiedConditionList TraitSetID="1"><ClassifiedCondition DB="MedGen" ID="{cui}">Test disease {cui}</ClassifiedCondition></ClassifiedConditionList>
        <RCVClassifications><GermlineClassification><ReviewStatus>{status}</ReviewStatus><Description DateLastEvaluated="2026-01-01" SubmissionCount="2">{significance}</Description></GermlineClassification></RCVClassifications>
      </RCVAccession>
    """


def _clinical_assertions() -> str:
    return """
    <ClinicalAssertionList>
      <ClinicalAssertion ID="10" ContributesToAggregateClassification="true">
        <ClinVarSubmissionID localKey="x"/>
        <ClinVarAccession Accession="SCV000010" Version="1" SubmitterName="Lab A" OrgID="1" OrganizationCategory="laboratory"/>
        <RecordStatus>current</RecordStatus>
        <Classification DateLastEvaluated="2026-01-01"><ReviewStatus>criteria provided, single submitter</ReviewStatus><GermlineClassification>Pathogenic</GermlineClassification></Classification>
        <Assertion>variation to disease</Assertion>
        <AttributeSet><Attribute Type="ModeOfInheritance">Autosomal recessive inheritance</Attribute></AttributeSet>
        <AttributeSet><Attribute Type="Penetrance">complete</Attribute></AttributeSet>
        <ObservedInList><ObservedIn>
          <Sample><Origin>germline</Origin><AffectedStatus>yes</AffectedStatus><NumberTested>2</NumberTested></Sample>
          <ObservedData><Attribute Type="VariantAlleles">2</Attribute></ObservedData>
          <Co-occurrenceSet><Zygosity>CompoundHeterozygote</Zygosity><AlleleDescSet><Name>second allele</Name><RelativeOrientation>trans</RelativeOrientation><Zygosity>SingleHeterozygote</Zygosity></AlleleDescSet><Count>2</Count></Co-occurrenceSet>
        </ObservedIn></ObservedInList>
        <SimpleAllele><GeneList><Gene Symbol="MEFV"/></GeneList></SimpleAllele>
        <TraitSet Type="Disease"><Trait Type="Disease"><XRef DB="MedGen" ID="C-ELIGIBLE"/></Trait></TraitSet>
      </ClinicalAssertion>
      <ClinicalAssertion ID="11" ContributesToAggregateClassification="true">
        <ClinVarSubmissionID localKey="y"/>
        <ClinVarAccession Accession="SCV000011" Version="1" SubmitterName="Lab B"/>
        <RecordStatus>current</RecordStatus>
        <Classification><ReviewStatus>criteria provided, single submitter</ReviewStatus><GermlineClassification>Uncertain significance</GermlineClassification></Classification>
        <Assertion>variation to disease</Assertion>
        <SimpleAllele/>
        <TraitSet Type="Disease"><Trait Type="Disease"><XRef DB="MedGen" ID="C-LOW"/></Trait></TraitSet>
      </ClinicalAssertion>
      <ClinicalAssertion ID="12" ContributesToAggregateClassification="true">
        <ClinVarSubmissionID localKey="z"/>
        <ClinVarAccession Accession="SCV000012" Version="1" SubmitterName="Old Lab"/>
        <RecordStatus>replaced</RecordStatus>
        <Classification><ReviewStatus>criteria provided, single submitter</ReviewStatus><GermlineClassification>Pathogenic</GermlineClassification></Classification>
        <Assertion>variation to disease</Assertion><SimpleAllele/>
        <TraitSet Type="Disease"><Trait Type="Disease"><XRef DB="MedGen" ID="C-ELIGIBLE"/></Trait></TraitSet>
      </ClinicalAssertion>
      <ClinicalAssertion ID="13" ContributesToAggregateClassification="true">
        <ClinVarSubmissionID localKey="sparse"/>
        <ClinVarAccession Accession="SCV000013" Version="1" SubmitterName="Lab C"/>
        <RecordStatus>current</RecordStatus>
        <Classification><ReviewStatus>criteria provided, single submitter</ReviewStatus><GermlineClassification>Likely pathogenic</GermlineClassification></Classification>
        <Assertion>variation to disease</Assertion>
        <ObservedInList><ObservedIn><Sample><Origin>germline</Origin><AffectedStatus>unknown</AffectedStatus></Sample><ObservedData><Attribute Type="Description">not provided</Attribute></ObservedData></ObservedIn></ObservedInList>
        <SimpleAllele/>
        <TraitSet Type="Disease"><Trait Type="Disease"><XRef DB="MedGen" ID="C-ELIGIBLE"/></Trait></TraitSet>
      </ClinicalAssertion>
    </ClinicalAssertionList>
    <TraitMappingList><TraitMapping ClinicalAssertionID="10" TraitType="Disease" MappingType="XRef" MappingValue="C-ELIGIBLE"><MedGen Name="Test disease C-ELIGIBLE" CUI="C-ELIGIBLE"/></TraitMapping></TraitMappingList>
    """


def _archive(
    variation_id: int,
    body: str,
    rcvs: str,
    assertions: str = "",
    status: str = "current",
    record_type: str = "classified",
) -> str:
    record = (
        f"<ClassifiedRecord>{body}<RCVList>{rcvs}</RCVList>{assertions}</ClassifiedRecord>"
        if record_type == "classified"
        else f"<IncludedRecord>{body}<Classifications/><SubmittedClassificationList><SCV Accession='SCV1' Version='1'/></SubmittedClassificationList><ClassifiedVariationList><ClassifiedVariation VariationID='1' Version='1'/></ClassifiedVariationList></IncludedRecord>"
    )
    return f"""
    <VariationArchive VariationID="{variation_id}" VariationName="test {variation_id}" VariationType="single nucleotide variant" Accession="VCV{variation_id:09d}" Version="1" RecordType="{record_type}" NumberOfSubmissions="2" NumberOfSubmitters="2" DateCreated="2026-01-01">
      <RecordStatus>{status}</RecordStatus><Species>Homo sapiens</Species>{record}
    </VariationArchive>
    """


def _write_fixture(path: Path) -> None:
    eligible = _rcv(
        "RCV000001",
        "C-ELIGIBLE",
        "criteria provided, multiple submitters, no conflicts",
    )
    low = _rcv("RCV000002", "C-LOW", "criteria provided, single submitter")
    expert = _rcv("RCV000003", "C-HAP", "reviewed by expert panel")
    guideline = _rcv("RCV000004", "C-GENO", "practice guideline")
    somatic = "<RCVAccession Accession='RCV000005' Version='1'><ClassifiedConditionList TraitSetID='2'/><RCVClassifications><SomaticClinicalImpact><ReviewStatus>reviewed by expert panel</ReviewStatus><Description SubmissionCount='1'>Tier I</Description></SomaticClinicalImpact></RCVClassifications></RCVAccession>"
    xml = "<ClinVarVariationRelease ReleaseDate='2026-07-21'>" + "".join(
        [
            _archive(1, _simple_allele(1, 100), eligible + low, _clinical_assertions()),
            _archive(2, f"<Haplotype VariationID='200'>{_simple_allele(2, 200)}<Name>hap</Name><VariationType>Haplotype</VariationType></Haplotype>", expert),
            _archive(3, f"<Genotype VariationID='300'>{_simple_allele(3, 300)}<Name>genotype</Name><VariationType>Diplotype</VariationType></Genotype>", guideline),
            _archive(4, _simple_allele(4, 400), eligible, status="replaced"),
            _archive(5, _simple_allele(5, 500), eligible, record_type="included"),
            _archive(6, _simple_allele(6, 600), somatic),
            _archive(7, f"<Genotype VariationID='700'><Haplotype VariationID='701'>{_simple_allele(7, 700)}<Name>nested haplotype</Name><VariationType>Haplotype</VariationType></Haplotype><Name>complex genotype</Name><VariationType>Diplotype</VariationType></Genotype>", expert),
        ]
    ) + "</ClinVarVariationRelease>"
    path.write_text(xml, encoding="utf-8")


def _write_exact(path: Path) -> None:
    rows = []
    for pos in range(100, 701, 100):
        rows.append(
            {
                "HGNC_Symbol": "MEFV",
                "gofcards_accession_id": f"GF{pos}",
                "gofcards_variant_id": f"SNV|16|{pos}|{pos}|A|G",
                "allele_key": f"SNV|16|{pos}|{pos}|A|G",
                "hg19_vcf_key": f"chr16|{pos}|a|g",
                "hg38_vcf_key": f"16|{pos + 100}|A|G",
                "disease": "Familial Mediterranean fever",
                "pmids": "1",
                "pscore": "5",
            }
        )
    rows.append(
        {
            **rows[0],
            "HGNC_Symbol": "INPP5F",
            "gofcards_accession_id": "coordinate-collision",
        }
    )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_review_star_mapping_is_conservative() -> None:
    assert review_stars("criteria provided, single submitter") == 1
    assert review_stars("criteria provided, multiple submitters, no conflicts") == 2
    assert review_stars("reviewed by expert panel") == 3
    assert review_stars("practice guideline") == 4
    assert review_stars("future unknown status") == 0


def test_review_tier_summary_separates_lower_review_vcvs() -> None:
    def matched_vcv(vcv: str, gene: str, stars: list[int]) -> dict:
        return {
            "gene_symbols": [gene],
            "ClinVar_VCV": {
                "variation": {"vcv_accession": vcv},
                "condition_assertions": [
                    {
                        "rcv_accession": f"RCV{vcv[-1]}{index}",
                        "conditions": ["test"],
                        "germline_classification": {
                            "review_status": f"status {star}",
                            "review_stars": star,
                            "clinical_significance": "Pathogenic",
                        },
                    }
                    for index, star in enumerate(stars)
                ],
            },
        }

    summary, rows = summarize_review_tiers(
        [
            matched_vcv("VCV1", "GENE1", [0, 1]),
            matched_vcv("VCV2", "GENE2", [1, 2]),
            matched_vcv("VCV3", "GENE2", [3]),
        ],
        {"exact_variant_archives": 3},
    )

    assert summary["gene_concordant_vcvs_with_any_germline_rcv"] == 3
    assert summary["gene_concordant_vcvs_with_any_rcv_at_least_2_stars"] == 2
    assert summary["gene_concordant_vcvs_with_rcvs_but_none_at_least_2_stars"] == 1
    assert summary["genes_with_any_vcv_at_least_2_stars"] == 1
    assert summary["genes_with_only_lower_review_matched_vcvs"] == 1
    assert summary["rcv_review_star_counts"] == {"0": 1, "1": 2, "2": 1, "3": 1}
    assert len(rows) == 5


def test_release_metadata_parsing() -> None:
    assert (
        _parse_remote_md5(
            "6648ff42b4e03a881c3e5fb04931c31a *ClinVarVCVRelease.xml.gz"
        )
        == "6648ff42b4e03a881c3e5fb04931c31a"
    )
    version, url = _latest_clinvar_xsd_url(
        '<a href="ClinVar_VCV_2.5.xsd">old</a>'
        '<a href="ClinVar_VCV_2.6.xsd">current</a>'
    )
    assert version == "2.6"
    assert url.endswith("/ClinVar_VCV_2.6.xsd")


def test_vcv_fetch_uses_remote_md5_before_large_download(
    tmp_path: Path, monkeypatch
) -> None:
    payload = b"fake compressed VCV release"
    expected_md5 = hashlib.md5(payload).hexdigest()  # nosec: fixture checksum
    large_downloads = []

    def fake_text(url, *_args, **_kwargs):
        if url.endswith(".md5"):
            return f"{expected_md5} *ClinVarVCVRelease.xml.gz\n"
        if url.endswith("_README"):
            return "ClinVar XML README\n"
        if url.endswith("xsd_public/"):
            return '<a href="ClinVar_VCV_2.6.xsd">schema</a>'
        if url.endswith("ClinVar_VCV_2.6.xsd"):
            return "<xs:schema xmlns:xs='http://www.w3.org/2001/XMLSchema'/>"
        raise AssertionError(url)

    def fake_large_download(_url, out_path, **_kwargs):
        large_downloads.append(str(out_path))
        out_path.write_bytes(payload)

    monkeypatch.setattr(
        "build_gene_nonlof_mechanism_cache._download_text", fake_text
    )
    monkeypatch.setattr(
        "build_gene_nonlof_mechanism_cache.download_url_with_curl",
        fake_large_download,
    )

    first = fetch_clinvar_vcv(
        raw_dir=tmp_path,
        previous_meta={},
        force=True,
        timeout=1,
        retries=1,
        proxy_url="",
        download_tool="curl",
        max_download_seconds=60,
    )
    assert first["status"] == "downloaded"
    assert first["md5"] == expected_md5
    assert first["format_metadata"]["xsd_version"] == "2.6"
    assert len(large_downloads) == 1

    second = fetch_clinvar_vcv(
        raw_dir=tmp_path,
        previous_meta={"clinvar_vcv_weekly": first},
        force=True,
        timeout=1,
        retries=1,
        proxy_url="",
        download_tool="curl",
        max_download_seconds=60,
    )
    assert second["status"] == "checked_same_remote_md5"
    assert len(large_downloads) == 1


def test_streaming_parser_filters_and_preserves_edge_context(tmp_path: Path) -> None:
    xml_path = tmp_path / "fixture.xml"
    exact_path = tmp_path / "exact.tsv"
    _write_fixture(xml_path)
    _write_exact(exact_path)

    lookup = load_exact_gofcards_lookup(exact_path)
    matches, stats = stream_parse_clinvar_vcv(xml_path, lookup, min_review_stars=2)

    assert len(matches) == 4
    by_vcv = {
        match["ClinVar_VCV"]["variation"]["vcv_accession"]: match
        for match in matches
    }
    assert set(by_vcv) == {
        "VCV000000001",
        "VCV000000002",
        "VCV000000003",
        "VCV000000007",
    }
    assert by_vcv["VCV000000001"]["gene_symbols"] == ["MEFV"]

    simple = by_vcv["VCV000000001"]["ClinVar_VCV"]
    assert simple["match"]["gene_concordance"] is True
    assert simple["match"]["discarded_gene_discordant_symbols"] == ["INPP5F"]
    assert {
        item["HGNC_Symbol"] for item in simple["match"]["matched_gofcards_records"]
    } == {"MEFV"}
    assert simple["variation"]["classification_scope"] == "simple_allele"
    assert len(simple["condition_assertions"]) == 1
    condition = simple["condition_assertions"][0]
    assert condition["germline_classification"]["review_stars"] == 2
    assert [item["scv_accession"] for item in condition["matched_scvs"]] == [
        "SCV000010",
        "SCV000013",
    ]
    scv = condition["matched_scvs"][0]
    assert scv["submitted_mode_of_inheritance"] == [
        "Autosomal recessive inheritance"
    ]
    assert scv["penetrance"] == ["complete"]
    assert scv["observed_zygosity_counts"] == {"CompoundHeterozygote": 2}
    cooccurrence = scv["observations"][0]["cooccurrence_sets"][0]
    assert cooccurrence["cooccurring_alleles"][0]["relative_orientation"] == "trans"
    assert "observed_zygosity_counts" not in condition["matched_scvs"][1]
    assert condition["matched_scvs"][1]["observations"][0]["affected_status"] == "unknown"
    assert simple["allelic_requirement"]["value"] == "unresolved"

    assert (
        by_vcv["VCV000000002"]["ClinVar_VCV"]["variation"][
            "matched_component_context"
        ]
        == "haplotype_component"
    )
    assert (
        by_vcv["VCV000000003"]["ClinVar_VCV"]["variation"][
            "matched_component_context"
        ]
        == "genotype_component"
    )
    assert (
        by_vcv["VCV000000007"]["ClinVar_VCV"]["variation"][
            "matched_component_context"
        ]
        == "genotype_haplotype_component"
    )
    assert stats["skipped_noncurrent"] == 1
    assert stats["skipped_included_record"] == 1
    assert stats["exact_archives_without_eligible_rcv"] == 1


def test_clinvar_matches_are_injected_into_existing_canonical_json(
    tmp_path: Path,
) -> None:
    xml_path = tmp_path / "fixture.xml"
    exact_path = tmp_path / "exact.tsv"
    _write_fixture(xml_path)
    _write_exact(exact_path)
    matches, _ = stream_parse_clinvar_vcv(
        xml_path,
        load_exact_gofcards_lookup(exact_path),
        min_review_stars=2,
    )
    canonical = {
        "_meta": {
            "version": "2.0",
            "built_at": "2026-07-23T00:00:00+00:00",
            "total_genes": 1,
            "sources": {},
        },
        "HGNC:6998": {
            "symbol": "MEFV",
            "mechanisms": ["GOF"],
            "variant_level": [
                {
                    "GoFCards": {
                        "mechanism": "GOF",
                        "source_record_id": "1",
                        "allele_key": "SNV|16|100|100|A|G",
                        "exact_normalization_status": "unmatched_public_source_allele",
                    }
                }
            ],
        },
    }
    hgnc = {
        "MEFV": {
            "hgnc_id": "HGNC:6998",
            "symbol": "MEFV",
            "entrez_id": "4210",
            "ensembl_id": "ENSG00000103313",
        }
    }

    stats = inject_clinvar_matches(canonical, matches, hgnc)

    assert stats["injected_gene_matches"] == 4
    assert len(canonical["HGNC:6998"]["variant_level"]) == 5
    assert canonical["HGNC:6998"]["variant_level"][1]["ClinVar_VCV"][
        "allelic_requirement"
    ]["value"] == "unresolved"
    schema = json.loads(
        (
            ROOT
            / "data/gene_pathogenic_mechanism/schema/"
            "gene_nonlof_mechanism_curated_assertions.schema.json"
        ).read_text(encoding="utf-8")
    )
    Draft202012Validator(schema).validate(canonical)
    validate_canonical_json(
        canonical,
        ROOT
        / "data/gene_pathogenic_mechanism/schema/"
        "gene_nonlof_mechanism_curated_assertions.schema.json",
    )


def test_scv_linked_to_multiple_eligible_rcvs_is_marked_ambiguous(
    tmp_path: Path,
) -> None:
    xml_path = tmp_path / "ambiguous.xml"
    exact_path = tmp_path / "exact.tsv"
    eligible_a = _rcv(
        "RCV000010",
        "C-ELIGIBLE",
        "criteria provided, multiple submitters, no conflicts",
    )
    eligible_b = _rcv(
        "RCV000011",
        "C-ELIGIBLE",
        "reviewed by expert panel",
    )
    xml_path.write_text(
        "<ClinVarVariationRelease ReleaseDate='2026-07-21'>"
        + _archive(
            1,
            _simple_allele(1, 100),
            eligible_a + eligible_b,
            _clinical_assertions(),
        )
        + "</ClinVarVariationRelease>",
        encoding="utf-8",
    )
    _write_exact(exact_path)

    matches, _ = stream_parse_clinvar_vcv(
        xml_path,
        load_exact_gofcards_lookup(exact_path),
        min_review_stars=2,
    )

    conditions = matches[0]["ClinVar_VCV"]["condition_assertions"]
    for condition in conditions:
        scv = next(
            item
            for item in condition["matched_scvs"]
            if item["scv_accession"] == "SCV000010"
        )
        assert scv["trait_linkage_ambiguous_across_eligible_rcvs"] is True
        assert scv["linked_eligible_rcv_accessions"] == [
            "RCV000010",
            "RCV000011",
        ]

def test_exact_coordinate_with_contradictory_gene_is_rejected(tmp_path: Path) -> None:
    xml_path = tmp_path / "gene_collision.xml"
    exact_path = tmp_path / "exact.tsv"
    xml_path.write_text(
        "<ClinVarVariationRelease ReleaseDate='2026-07-21'>"
        + _archive(
            1,
            _simple_allele(1, 100, symbol="THPO"),
            _rcv(
                "RCV000010",
                "C-ELIGIBLE",
                "criteria provided, multiple submitters, no conflicts",
            ),
        )
        + "</ClinVarVariationRelease>",
        encoding="utf-8",
    )
    _write_exact(exact_path)

    matches, stats = stream_parse_clinvar_vcv(
        xml_path,
        load_exact_gofcards_lookup(exact_path),
        min_review_stars=2,
    )

    assert matches == []
    assert stats["skipped_gene_discordant_component_matches"] == 1


def test_non_clinvar_source_assertion_schema_contracts() -> None:
    schema = json.loads(
        (
            ROOT
            / "data/gene_pathogenic_mechanism/schema/"
            "gene_nonlof_mechanism_curated_assertions.schema.json"
        ).read_text(encoding="utf-8")
    )
    canonical = {
        "_meta": {
            "version": "2.0",
            "built_at": "2026-07-23T00:00:00+00:00",
            "total_genes": 1,
            "sources": {},
        },
        "HGNC:1": {
            "symbol": "A1BG",
            "mechanisms": [
                "GOF",
                "PANELAPP_GREEN_NON_LOF_PATHO_HISTORY",
                "TRIPLOSENSITIVITY",
            ],
            "gene_level": [
                {
                    "G2P_DDG2P": {
                        "mechanism": "GOF",
                        "mechanism_raw": "gain of function",
                        "disease": "test disease",
                        "inheritance": "monoallelic",
                        "confidence": "high",
                        "pmids": ["1"],
                    }
                },
                {
                    "PanelApp": {
                        "mechanism": "PANELAPP_GREEN_NON_LOF_PATHO_HISTORY",
                        "confidence": "high",
                        "panel": "test panel",
                    }
                },
                {
                    "ClinGen_Dosage": {
                        "mechanism": "TRIPLOSENSITIVITY",
                        "score": "2",
                    }
                },
            ],
            "variant_level": [
                {
                    "GoFCards": {
                        "mechanism": "GOF",
                        "source_record_id": "1",
                        "allele_key": "SNV|1|1|1|A|G",
                        "exact_normalization_status": "unmatched_public_source_allele",
                        "chr": "chr1",
                        "pos": "1",
                        "ref": "A",
                        "alt": "G",
                        "consequence": "missense_variant",
                        "animal_model": True,
                    }
                }
            ],
        },
    }
    validator = Draft202012Validator(schema)
    validator.validate(canonical)

    canonical["HGNC:1"]["gene_level"][1]["PanelApp"]["confidence"] = "moderate"
    with pytest.raises(ValidationError):
        validator.validate(canonical)

    canonical["HGNC:1"]["gene_level"][1]["PanelApp"]["confidence"] = "high"
    canonical["HGNC:1"]["variant_level"][0]["GoFCards"]["disease"] = ""
    with pytest.raises(ValidationError):
        validator.validate(canonical)


def test_gofcards_assertion_preserves_all_gene_concordant_exact_rows(
    tmp_path: Path,
) -> None:
    exact_path = tmp_path / "gofcards_exact.tsv"
    base = {column: "" for column in GOFCARDS_EXACT_COLUMNS}
    base.update(
        {
            "source": "GoFCards",
            "mechanism": "GOF",
            "build": "hg19_and_hg38",
            "HGNC_Symbol": "MEFV",
            "feature_type": "Transcript",
            "consequence": "missense_variant",
            "match_status": "both_cdna_protein_match",
            "raw_GoFCards_HGVS": "MEFV:NM_TEST.1:c.1A>G:p.Lys1Arg",
            "GoFCards_transcript": "NM_TEST.1",
            "hg19_chrom": "16",
            "hg19_pos": "100",
            "hg19_ref": "A",
            "hg19_alt": "G",
            "hg19_vcf_status": "raw_ref_alt",
            "hg38_chrom": "16",
            "hg38_pos": "200",
            "hg38_ref": "A",
            "hg38_alt": "G",
            "hg38_refalt_status": "snv_ref_match",
            "gofcards_variant_id": "SNV|16|100|100|A|G",
            "disease": "test disease",
            "pmids": "1;2",
            "pscore": "3",
            "derived_on": "2026-07-07",
            "allele_key": "SNV|16|100|100|A|G",
            "hg19_genomic_key": "16|100|A|G",
            "hg19_vcf_pos": "100",
            "hg19_vcf_ref": "A",
            "hg19_vcf_alt": "G",
            "hg38_vcf_pos": "200",
            "hg38_vcf_ref": "A",
            "hg38_vcf_alt": "G",
            "hg19_vcf_key": "16|100|A|G",
            "hg38_genomic_key": "16|200|A|G",
            "hg38_vcf_key": "16|200|A|G",
            "match_key_types": "hgvsp;hg19_vcf;hg38_vcf",
        }
    )
    rows = []
    for assembly, transcript, hgvsc, hgvsp in (
        ("hg19", "ENST_TEST.1", "ENST_TEST.1:c.1A>G", "ENSP_TEST.1:p.Lys1Arg"),
        ("hg38", "ENST_TEST.2", "ENST_TEST.2:c.1A>G", "ENSP_TEST.2:p.Lys1Arg"),
    ):
        rows.append(
            {
                **base,
                "VEP_assembly": assembly,
                "VEP_transcript": transcript,
                "HGVSc": hgvsc,
                "HGVSp": hgvsp,
                "hgvsp_key": "LYS1ARG",
            }
        )
    pd.DataFrame(rows, columns=GOFCARDS_EXACT_COLUMNS).to_csv(
        exact_path, sep="\t", index=False
    )
    exact_by_allele, stats = load_gofcards_exact_records(exact_path)
    hgnc = {
        "MEFV": {
            "hgnc_id": "HGNC:6998",
            "symbol": "MEFV",
            "entrez_id": "4210",
            "ensembl_id": "ENSG00000103313",
        }
    }
    unified = pd.DataFrame(
        [
            make_unified_row(
                gene_symbol="MEFV",
                mechanism=["GOF"],
                assertion_level="variant_level",
                source="GoFCards",
                source_record_id="1",
                assembly="hg19",
                chromosome="chr16",
                position="100",
                ref="A",
                alt="G",
                allele_key="SNV|16|100|100|A|G",
                consequence="missense_variant",
            )
        ]
    )

    canonical = build_nonlof_assertions_json(unified, hgnc, exact_by_allele)
    assertion = canonical["HGNC:6998"]["variant_level"][0]["GoFCards"]
    gofcards_meta = canonical["_meta"]["sources"]["GoFCards"]

    assert stats["cache_rows"] == 2
    assert "assembly" not in gofcards_meta
    assert gofcards_meta["raw_public_allele_fields"] == {
        "assembly": "hg19",
        "keys": ["chr", "pos", "ref", "alt"],
        "note": (
            "Applies only to the legacy top-level GoFCards source fields; "
            "exact_normalized_variants identify hg19 and hg38 fields explicitly"
        ),
    }
    assert assertion["exact_normalization_status"] == "matched_gene_concordant"
    assert len(assertion["exact_normalized_variants"]) == 2
    assert {
        row["VEP_assembly"] for row in assertion["exact_normalized_variants"]
    } == {"hg19", "hg38"}
    assert assertion["exact_normalized_variants"][0]["pmids"] == ["1", "2"]
    assert assertion["exact_normalized_variants"][0]["match_key_types"] == [
        "hg19_vcf",
        "hg38_vcf",
        "hgvsp",
    ]

    schema = json.loads(
        (
            ROOT
            / "data/gene_pathogenic_mechanism/schema/"
            "gene_nonlof_mechanism_curated_assertions.schema.json"
        ).read_text(encoding="utf-8")
    )
    Draft202012Validator(schema).validate(canonical)


def test_hgnc_mapping_splits_comma_delimited_aliases(tmp_path: Path) -> None:
    mapping_path = tmp_path / "hgnc_gene_id_map.txt"
    mapping_path.write_text(
        "HGNC ID\tApproved symbol\tAlias symbols\tPrevious symbols\t"
        "NCBI Gene ID\tEnsembl ID\n"
        "HGNC:4764\tH3-3A\tH3.3A\tH3F3, H3F3A\t3020\tENSG00000163041\n",
        encoding="utf-8",
    )

    mapping = load_hgnc_mapping(mapping_path)

    assert mapping["H3F3A"]["hgnc_id"] == "HGNC:4764"
    assert mapping["H3F3"]["symbol"] == "H3-3A"
