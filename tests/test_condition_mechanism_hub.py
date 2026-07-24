from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gene_mechanism_hub import build_hpo_condition_index  # noqa: E402


def test_build_hpo_condition_index_preserves_assertion_provenance() -> None:
    assertions = pd.DataFrame(
        [
            {
                "gene_symbol": "GENE1",
                "disease_id": "OMIM:123",
                "mondo_id": "MONDO:0000123",
                "hpo_id": "HP:0000006",
                "frequency": "-",
                "evidence": "TAS",
                "reference": "OMIM:123",
                "disease_scope": "mendelian_non_neoplastic",
                "priva_scope": "include",
                "scope_review_status": "auto_supported",
            },
            {
                "gene_symbol": "GENE1",
                "disease_id": "OMIM:123",
                "mondo_id": "MONDO:0000123",
                "hpo_id": "HP:0003829",
                "frequency": "2/10",
                "evidence": "PCS",
                "reference": "PMID:111",
                "disease_scope": "mendelian_non_neoplastic",
                "priva_scope": "include",
                "scope_review_status": "auto_supported",
            },
            {
                "gene_symbol": "GENE1",
                "disease_id": "OMIM:123",
                "mondo_id": "MONDO:0000123",
                "hpo_id": "HP:0003584",
                "frequency": "1/10",
                "evidence": "PCS",
                "reference": "PMID:222",
                "disease_scope": "mendelian_non_neoplastic",
                "priva_scope": "include",
                "scope_review_status": "auto_supported",
            },
        ]
    )

    index = build_hpo_condition_index(assertions)
    by_omim = index[("GENE1", "OMIM:123")]
    by_mondo = index[("GENE1", "MONDO:0000123")]

    assert by_mondo is by_omim
    assert by_omim["inheritance_modes"] == ["Autosomal dominant inheritance"]
    assert by_omim["penetrance_hpo_ids"] == ["HP:0003829"]
    assert by_omim["onset_hpo_ids"] == ["HP:0003584"]
    assert by_omim["hpo_assertions"][1] == {
        "hpo_id": "HP:0003829",
        "frequency": "2/10",
        "evidence": "PCS",
        "reference": "PMID:111",
    }
