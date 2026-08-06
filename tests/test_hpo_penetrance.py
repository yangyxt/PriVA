import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from acmg_inheritance import hpo_onset_modes  # noqa: E402
from hpo_penetrance import (  # noqa: E402
    HPO_COMPLETE_PENETRANCE_TERMS,
    HPO_HIGH_PENETRANCE_TERMS,
    HPO_INCOMPLETE_PENETRANCE_EQUIVALENT_TERMS,
    HPO_PENETRANCE_STATUS_BY_TERM,
    normalize_penetrance_evidence,
    normalize_penetrance_hpo_ids,
)


INCOMPLETE_EQUIVALENT_TERMS = {
    "HP:0003829",  # Incomplete penetrance
    "HP:4000159",  # Moderate penetrance
    "HP:4000160",  # Low penetrance
    "HP:0003584",  # Late onset
}

NONBLOCKING_CONTEXT_TERMS = {
    "HP:0003831",  # Age-dependent penetrance / age-related onset
    "HP:0034857",  # Highly variable age of onset
    "HP:0003581",  # Adult onset
    "HP:0011462",  # Young adult onset
    "HP:0003596",  # Middle age onset
    "HP:0003587",  # Insidious / gradual onset
    "HP:0003828",  # Variable expressivity
}


@pytest.mark.parametrize("hpo_id", sorted(INCOMPLETE_EQUIVALENT_TERMS))
def test_incomplete_equivalent_terms_normalize_conservatively(hpo_id: str) -> None:
    assert normalize_penetrance_hpo_ids([hpo_id]) == "incomplete"
    assert hpo_onset_modes(hpo_id) is False


@pytest.mark.parametrize("hpo_id", sorted(NONBLOCKING_CONTEXT_TERMS))
def test_other_onset_and_expression_context_does_not_block(hpo_id: str) -> None:
    assert normalize_penetrance_hpo_ids([hpo_id]) == ""
    assert hpo_onset_modes(hpo_id) is True


def test_penetrance_vocabulary_has_an_explicit_boundary() -> None:
    assert HPO_INCOMPLETE_PENETRANCE_EQUIVALENT_TERMS == INCOMPLETE_EQUIVALENT_TERMS
    assert HPO_COMPLETE_PENETRANCE_TERMS == {"HP:0034950"}
    assert HPO_HIGH_PENETRANCE_TERMS == {"HP:4000158"}

    # High is complete in the binary contract. Sex-limited expression and slow
    # progression do not assert that some carriers fail to manifest.
    assert normalize_penetrance_hpo_ids(["HP:4000158"]) == "complete"
    assert normalize_penetrance_hpo_ids(["HP:0034950"]) == "complete"
    assert normalize_penetrance_hpo_ids(["HP:0001470"]) == ""
    assert normalize_penetrance_hpo_ids(["HP:0003677"]) == ""
    assert "HP:0003831" not in HPO_PENETRANCE_STATUS_BY_TERM
    assert "HP:0003828" not in HPO_PENETRANCE_STATUS_BY_TERM
    assert "HP:0001470" not in HPO_PENETRANCE_STATUS_BY_TERM
    assert "HP:0003677" not in HPO_PENETRANCE_STATUS_BY_TERM


def test_cautious_penetrance_precedence_handles_conflicting_assertions() -> None:
    assert (
        normalize_penetrance_hpo_ids("HP:0034950;HP:0003584")
        == "incomplete"
    )
    assert normalize_penetrance_hpo_ids("HP:0034950,HP:4000158") == "complete"
    assert hpo_onset_modes("HP:4000158;HP:0001470;HP:0003677") is True


def test_g2p_modifier_and_phenotype_evidence_share_binary_precedence() -> None:
    assert normalize_penetrance_evidence(
        raw_values="typically de novo; typified by incomplete penetrance"
    ) == "incomplete"
    assert normalize_penetrance_evidence(
        hpo_ids="HP:4000158;HP:0003829"
    ) == "incomplete"
    assert normalize_penetrance_evidence(
        hpo_ids="HP:0003828;HP:0003831"
    ) == ""
