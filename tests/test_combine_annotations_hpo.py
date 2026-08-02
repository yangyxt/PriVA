import json
import sys
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from combine_annotations import (  # noqa: E402
    HPO_ANNOTATION_COLUMNS,
    HPO_BS1_BS2_CONTEXT_IDS,
    build_hpo_symbol_map,
    hpo_annotation_for_symbol,
    select_record_hpo_annotations,
)
from hpo_penetrance import HPO_PENETRANCE_STATUS_BY_TERM  # noqa: E402


def _hpo_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "SYMBOL": "HRG",
                "HPO_IDs": "HP:0000006",
                "HPO_terms": "Autosomal dominant inheritance",
                "HPO_sources": "OMIM:613116",
                "HPO_gene_inheritance": "Autosomal dominant inheritance",
            },
            {
                "SYMBOL": "P2RX2",
                "HPO_IDs": "HP:0000006",
                "HPO_terms": "Autosomal dominant inheritance",
                "HPO_sources": "OMIM:608224",
                "HPO_gene_inheritance": "Autosomal dominant inheritance",
            },
        ]
    )


def test_hpo_map_is_keyed_by_gene_symbol() -> None:
    hpo_by_symbol = build_hpo_symbol_map(_hpo_table())

    assert set(hpo_by_symbol) == {"HRG", "P2RX2"}
    assert (
        hpo_by_symbol["HRG"]["HPO_gene_inheritance"]
        == "Autosomal dominant inheritance"
    )


def test_multi_gene_record_keeps_separate_hpo_annotations() -> None:
    hpo_by_symbol = build_hpo_symbol_map(_hpo_table())
    csq_entries = (
        "A|HRG|ENST00000232003",
        "A|P2RX2|ENST00000343948",
        "A|UNKNOWN|ENST00000999999",
    )

    record_hpo = select_record_hpo_annotations(
        csq_entries,
        symbol_field_index=1,
        hpo_by_symbol=hpo_by_symbol,
    )

    assert set(record_hpo) == {"HRG", "P2RX2"}
    assert hpo_annotation_for_symbol(record_hpo, "HRG")["HPO_sources"] == "OMIM:613116"
    assert hpo_annotation_for_symbol(record_hpo, "P2RX2")["HPO_sources"] == "OMIM:608224"

    missing = hpo_annotation_for_symbol(record_hpo, "UNKNOWN")
    assert set(missing) == set(HPO_ANNOTATION_COLUMNS)
    assert all(pd.isna(value) for value in missing.values())


def test_duplicate_gene_symbols_are_rejected() -> None:
    duplicated = pd.concat([_hpo_table(), _hpo_table().iloc[[0]]], ignore_index=True)

    try:
        build_hpo_symbol_map(duplicated)
    except ValueError as error:
        assert "HRG" in str(error)
    else:
        raise AssertionError("duplicate HPO gene symbols should fail explicitly")


def test_assertion_rows_for_one_gene_are_grouped_without_losing_disease() -> None:
    assertions = pd.DataFrame(
        [
            {
                "gene_symbol": "KCNH5",
                "disease_id": "OMIM:620537",
                "hpo_id": "HP:0003584",
                "frequency": "1/17",
                "evidence": "PCS",
                "reference": "PMID:36307226",
            },
            {
                "gene_symbol": "KCNH5",
                "disease_id": "OMIM:OTHER",
                "hpo_id": "HP:0003584",
                "frequency": "-",
                "evidence": "TAS",
                "reference": "OMIM:OTHER",
            },
            {
                "gene_symbol": "KCNH5",
                "disease_id": "OMIM:620537",
                "hpo_id": "HP:0000006",
                "frequency": "-",
                "evidence": "TAS",
                "reference": "OMIM:620537",
            },
            {
                "gene_symbol": "KCNH5",
                "disease_id": "OMIM:620537",
                "hpo_id": "HP:0001250",
                "frequency": "17/17",
                "evidence": "PCS",
                "reference": "PMID:36307226",
            },
        ]
    )

    hpo_by_symbol = build_hpo_symbol_map(assertions)
    annotation = hpo_by_symbol["KCNH5"]
    structured = json.loads(annotation["HPO_assertions"])

    assert annotation["HPO_IDs"] == "HP:0003584;HP:0000006"
    assert set(annotation["HPO_sources"].split(";")) == {
        "OMIM:620537",
        "OMIM:OTHER",
    }
    assert annotation["HPO_gene_inheritance"] == "Autosomal dominant inheritance"
    assert len(structured) == 3
    assert {item["disease_id"] for item in structured} == {
        "OMIM:620537",
        "OMIM:OTHER",
    }


def test_scope_filters_flat_hpo_context_but_preserves_audit_assertions() -> None:
    assertions = pd.DataFrame(
        [
            {
                "gene_symbol": "SCOPE1",
                "disease_id": "OMIM:INCLUDED",
                "hpo_id": "HP:0000006",
                "frequency": "-",
                "evidence": "TAS",
                "reference": "OMIM:INCLUDED",
                "mondo_id": "MONDO:1",
                "mondo_name": "included disease",
                "disease_scope": "mendelian_non_neoplastic",
                "priva_scope": "include",
                "scope_evidence": "HPO_mendelian_inheritance",
                "scope_reference": "OMIM:INCLUDED",
                "scope_review_status": "auto_supported",
            },
            {
                "gene_symbol": "SCOPE1",
                "disease_id": "OMIM:EXCLUDED",
                "hpo_id": "HP:0003829",
                "frequency": "-",
                "evidence": "TAS",
                "reference": "OMIM:EXCLUDED",
                "mondo_id": "MONDO:2",
                "mondo_name": "complex disease",
                "disease_scope": "complex_or_non_monogenic",
                "priva_scope": "exclude",
                "scope_evidence": "HPO_non_monogenic_inheritance",
                "scope_reference": "OMIM:EXCLUDED",
                "scope_review_status": "auto_supported",
            },
            {
                "gene_symbol": "SCOPE1",
                "disease_id": "ORPHA:REVIEW",
                "hpo_id": "HP:0003584",
                "frequency": "-",
                "evidence": "PCS",
                "reference": "PMID:3",
                "mondo_id": "MONDO:3",
                "mondo_name": "uncertain neoplasm",
                "disease_scope": "neoplastic_uncertain",
                "priva_scope": "review",
                "scope_evidence": "MONDO_ancestor:MONDO:0005070",
                "scope_reference": "MONDO:test",
                "scope_review_status": "needs_review",
            },
        ]
    )

    annotation = build_hpo_symbol_map(assertions)["SCOPE1"]
    structured = json.loads(annotation["HPO_assertions"])

    assert annotation["HPO_IDs"] == "HP:0000006"
    assert annotation["HPO_sources"] == "OMIM:INCLUDED"
    assert annotation["HPO_gene_inheritance"] == "Autosomal dominant inheritance"
    assert annotation["HPO_scope_review_required"] == 1
    assert annotation["HPO_scope_review_disease_ids"] == "ORPHA:REVIEW"
    assert annotation["HPO_scope_excluded_disease_ids"] == "OMIM:EXCLUDED"
    assert len(structured) == 3
    assert {row["priva_scope"] for row in structured} == {
        "include",
        "exclude",
        "review",
    }


def test_flat_hpo_context_uses_the_shared_penetrance_vocabulary() -> None:
    assert set(HPO_PENETRANCE_STATUS_BY_TERM) <= HPO_BS1_BS2_CONTEXT_IDS
    assert "HP:0003587" in HPO_BS1_BS2_CONTEXT_IDS  # Insidious / gradual onset
    assert "HP:0034857" in HPO_BS1_BS2_CONTEXT_IDS  # Variable age of onset

    # These assertions remain deliberately outside the penetrance gate.
    assert "HP:0001470" not in HPO_BS1_BS2_CONTEXT_IDS  # Sex-limited expression
    assert "HP:0003677" not in HPO_BS1_BS2_CONTEXT_IDS  # Slowly progressive
    assert "HP:0031785" not in HPO_BS1_BS2_CONTEXT_IDS  # Stale insidious ID
