from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_gofcards_source_coverage import load_gofcards_gene_counts  # noqa: E402


def test_load_gofcards_gene_counts_normalizes_and_counts(tmp_path: Path) -> None:
    cache = tmp_path / "gofcards.tsv"
    cache.write_text(
        "HGNC_Symbol\tgofcards_variant_id\n"
        "OLD1\tv1\n"
        "GENE1\tv2\n"
        "GENE2\tv3\n"
        "\tv4\n",
        encoding="utf-8",
    )
    aliases = {"OLD1": "GENE1"}

    counts = load_gofcards_gene_counts(
        cache,
        resolve_symbol=lambda symbol: aliases.get(symbol, symbol),
    )

    assert counts == {"GENE1": 2, "GENE2": 1}
