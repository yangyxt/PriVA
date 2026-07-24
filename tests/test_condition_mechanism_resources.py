from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_gene_pathogenic_mechanism_cache import parse_g2p  # noqa: E402


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
