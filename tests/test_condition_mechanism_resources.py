from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_gene_pathogenic_mechanism_cache import (  # noqa: E402
    attach_mondo_condition_identity,
    parse_g2p,
    parse_orphadata,
)


def test_parse_g2p_preserves_assertion_and_condition_identifiers(
    tmp_path: Path,
) -> None:
    source = tmp_path / "AllG2P.csv"
    pd.DataFrame(
        [
            {
                "g2p id": "G2P00999",
                "gene symbol": "TEST1",
                "hgnc id": "1234",
                "disease name": "TEST1-related disorder",
                "disease mim": "612345",
                "disease MONDO": "MONDO:0012345",
                "allelic requirement": "monoallelic_autosomal",
                "confidence": "definitive",
                "molecular mechanism": "gain of function",
                "publications": "12345678",
            }
        ]
    ).to_csv(source, index=False)

    result = parse_g2p(source).iloc[0]

    assert result["source_record_id"] == "G2P00999"
    assert result["source_condition_id"] == "OMIM:612345"
    assert result["mondo_id"] == "MONDO:0012345"
    assert result["mechanism"] == "GOF"


def test_parse_orphadata_keeps_only_assessed_explicit_germline_mechanisms(
    tmp_path: Path,
) -> None:
    source = tmp_path / "en_product6.xml"
    source.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<JDBOR>
  <DisorderList>
    <Disorder>
      <OrphaCode>123</OrphaCode><Name>Condition one</Name>
      <DisorderGeneAssociationList>
        <DisorderGeneAssociation>
          <SourceOfValidation>111[PMID];222[PMID]</SourceOfValidation>
          <Gene><Symbol>GENE1</Symbol></Gene>
          <DisorderGeneAssociationType><Name>Disease-causing germline mutation(s) (gain of function) in</Name></DisorderGeneAssociationType>
          <DisorderGeneAssociationStatus><Name>Assessed</Name></DisorderGeneAssociationStatus>
        </DisorderGeneAssociation>
        <DisorderGeneAssociation>
          <Gene><Symbol>GENE2</Symbol></Gene>
          <DisorderGeneAssociationType><Name>Disease-causing germline mutation(s) in</Name></DisorderGeneAssociationType>
          <DisorderGeneAssociationStatus><Name>Assessed</Name></DisorderGeneAssociationStatus>
        </DisorderGeneAssociation>
        <DisorderGeneAssociation>
          <Gene><Symbol>GENE3</Symbol></Gene>
          <DisorderGeneAssociationType><Name>Disease-causing somatic mutation(s) in</Name></DisorderGeneAssociationType>
          <DisorderGeneAssociationStatus><Name>Assessed</Name></DisorderGeneAssociationStatus>
        </DisorderGeneAssociation>
      </DisorderGeneAssociationList>
    </Disorder>
    <Disorder>
      <OrphaCode>456</OrphaCode><Name>Condition two</Name>
      <DisorderGeneAssociationList>
        <DisorderGeneAssociation>
          <Gene><Symbol>GENE4</Symbol></Gene>
          <DisorderGeneAssociationType><Name>Disease-causing germline mutation(s) (loss of function) in</Name></DisorderGeneAssociationType>
          <DisorderGeneAssociationStatus><Name>Assessed</Name></DisorderGeneAssociationStatus>
        </DisorderGeneAssociation>
        <DisorderGeneAssociation>
          <Gene><Symbol>GENE5</Symbol></Gene>
          <DisorderGeneAssociationType><Name>Disease-causing germline mutation(s) (gain of function) in</Name></DisorderGeneAssociationType>
          <DisorderGeneAssociationStatus><Name>Not yet assessed</Name></DisorderGeneAssociationStatus>
        </DisorderGeneAssociation>
      </DisorderGeneAssociationList>
    </Disorder>
  </DisorderList>
</JDBOR>
""",
        encoding="utf-8",
    )

    result = parse_orphadata(source).set_index("gene_symbol")

    assert set(result.index) == {"GENE1", "GENE4"}
    assert result.loc["GENE1", "source_condition_id"] == "ORPHA:123"
    assert result.loc["GENE1", "mechanism"] == "GOF"
    assert result.loc["GENE1", "pmids"] == "111;222"
    assert result.loc["GENE4", "source_condition_id"] == "ORPHA:456"
    assert result.loc["GENE4", "mechanism"] == "LOF"


def test_attach_mondo_condition_identity_maps_source_ids_without_overwrite(
    tmp_path: Path,
) -> None:
    registry = tmp_path / "disease_scope.tsv"
    pd.DataFrame(
        [
            {
                "disease_id": "ORPHA:123",
                "mondo_id": "MONDO:0000123",
                "disease_scope": "mendelian_non_neoplastic",
                "priva_scope": "include",
                "scope_review_status": "auto_supported",
            },
            {
                "disease_id": "OMIM:456",
                "mondo_id": "MONDO:REGISTRY",
                "disease_scope": "neoplastic_uncertain",
                "priva_scope": "review",
                "scope_review_status": "manual_review_unresolved",
            },
        ]
    ).to_csv(registry, sep="\t", index=False)
    source = pd.DataFrame(
        [
            {
                "source_condition_id": "ORPHA:123",
                "mondo_id": "",
            },
            {
                "source_condition_id": "OMIM:456",
                "mondo_id": "MONDO:SOURCE",
            },
            {
                "source_condition_id": "ORPHA:999",
                "mondo_id": "",
            },
        ]
    )

    result = attach_mondo_condition_identity(source, registry)

    assert result.loc[0, "mondo_id"] == "MONDO:0000123"
    assert result.loc[0, "priva_scope"] == "include"
    assert result.loc[1, "mondo_id"] == "MONDO:SOURCE"
    assert result.loc[1, "priva_scope"] == "review"
    assert result.loc[2, "mondo_id"] == ""
    assert result.loc[2, "priva_scope"] == ""
