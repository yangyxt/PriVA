import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_mondo_disease_scope import (  # noqa: E402
    annotate_hpo_assertions,
    build_disease_scope_registry,
)


def _write_mondo(tmp_path: Path) -> Path:
    mondo = tmp_path / "mondo.obo"
    mondo.write_text(
        """format-version: 1.2
data-version: mondo/releases/test/mondo-simple.owl

[Term]
id: MONDO:0000001
name: disease

[Term]
id: MONDO:0003847
name: hereditary disease
is_a: MONDO:0000001 ! disease

[Term]
id: MONDO:0005070
name: neoplasm
is_a: MONDO:0000001 ! disease

[Term]
id: MONDO:1000001
name: inherited test disease
xref: OMIM:1 {source="MONDO:equivalentTo"}
is_a: MONDO:0003847 ! hereditary disease

[Term]
id: MONDO:1000002
name: hereditary cancer test disease
xref: OMIM:2 {source="MONDO:equivalentTo"}
is_a: MONDO:0003847 ! hereditary disease
is_a: MONDO:0005070 ! neoplasm

[Term]
id: MONDO:1000003
name: neoplasm with unresolved etiology
xref: Orphanet:3 {source="MONDO:equivalentTo"}
is_a: MONDO:0005070 ! neoplasm

[Term]
id: MONDO:1000004
name: susceptibility candidate
subset: omim_susceptibility
xref: OMIM:4 {source="MONDO:equivalentTo"}
is_a: MONDO:0000001 ! disease

[Term]
id: MONDO:1000006
name: inherited susceptibility without neoplasm ancestry
subset: omim_susceptibility
xref: OMIM:6 {source="MONDO:equivalentTo"}
is_a: MONDO:0003847 ! hereditary disease
""",
        encoding="utf-8",
    )
    return mondo


def _write_assertions(tmp_path: Path) -> Path:
    assertions = pd.DataFrame(
        [
            ["GENE1", "OMIM:1", "HP:0000006", "-", "TAS", "OMIM:1"],
            ["GENE2", "OMIM:2", "HP:0003829", "-", "PCS", "PMID:2"],
            ["GENE3", "ORPHA:3", "HP:0003584", "1/2", "PCS", "PMID:3"],
            ["GENE4", "OMIM:4", "HP:0010982", "-", "TAS", "OMIM:4"],
            ["GENE5", "OMIM:5", "HP:0001250", "1/1", "PCS", "PMID:5"],
            ["GENE6", "OMIM:6", "HP:0001250", "1/1", "PCS", "PMID:6"],
        ],
        columns=[
            "gene_symbol",
            "disease_id",
            "hpo_id",
            "frequency",
            "evidence",
            "reference",
        ],
    )
    path = tmp_path / "assertions.tsv"
    assertions.to_csv(path, sep="\t", index=False)
    return path


def test_builds_conservative_mondo_hpo_scope_registry(tmp_path: Path) -> None:
    registry = build_disease_scope_registry(
        _write_mondo(tmp_path),
        _write_assertions(tmp_path),
    ).set_index("disease_id")

    assert registry.loc["OMIM:1", "disease_scope"] == "mendelian_non_neoplastic"
    assert registry.loc["OMIM:1", "priva_scope"] == "include"
    assert registry.loc["OMIM:2", "disease_scope"] == "mendelian_cancer_predisposition"
    assert registry.loc["OMIM:2", "priva_scope"] == "include"
    assert registry.loc["ORPHA:3", "mondo_id"] == "MONDO:1000003"
    assert registry.loc["ORPHA:3", "disease_scope"] == "neoplastic_uncertain"
    assert registry.loc["ORPHA:3", "priva_scope"] == "review"
    assert registry.loc["OMIM:4", "disease_scope"] == "complex_or_non_monogenic"
    assert registry.loc["OMIM:4", "priva_scope"] == "exclude"
    assert registry.loc["OMIM:5", "disease_scope"] == "uncertain"
    assert (
        registry.loc["OMIM:6", "disease_scope"]
        == "mendelian_susceptibility_uncertain"
    )
    assert registry.loc["OMIM:6", "priva_scope"] == "review"


def test_manual_override_can_confirm_somatic_only_disease(tmp_path: Path) -> None:
    overrides = tmp_path / "overrides.tsv"
    overrides.write_text(
        "disease_id\tdisease_scope\tpriva_scope\tscope_review_status\t"
        "scope_evidence\tscope_reference\tscope_note\n"
        "ORPHA:3\tsomatic_oncogenic\texclude\tmanually_confirmed\t"
        "manual_review\tPMID:99\tSomatic-only condition\n",
        encoding="utf-8",
    )
    registry = build_disease_scope_registry(
        _write_mondo(tmp_path),
        _write_assertions(tmp_path),
        manual_overrides=overrides,
    ).set_index("disease_id")

    assert registry.loc["ORPHA:3", "disease_scope"] == "somatic_oncogenic"
    assert registry.loc["ORPHA:3", "priva_scope"] == "exclude"
    assert registry.loc["ORPHA:3", "scope_review_status"] == "manually_confirmed"


def test_scope_is_attached_to_every_assertion(tmp_path: Path) -> None:
    assertion_path = _write_assertions(tmp_path)
    assertions = pd.read_csv(assertion_path, sep="\t", dtype=str)
    registry = build_disease_scope_registry(_write_mondo(tmp_path), assertion_path)

    annotated = annotate_hpo_assertions(assertions, registry)

    assert len(annotated) == len(assertions)
    assert annotated["priva_scope"].isna().sum() == 0
    assert set(annotated["priva_scope"]) == {"include", "exclude", "review"}
