from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_gofcards_source_coverage import (  # noqa: E402
    load_condition_mechanism_gene_sets,
    load_gofcards_gene_counts,
)


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


def test_load_condition_mechanism_gene_sets_separates_scope(tmp_path: Path) -> None:
    evidence = tmp_path / "evidence.tsv"
    evidence.write_text(
        "gene_symbol\tsource\tnormalized_mechanisms\tpriva_scope\n"
        "OLD1\tG2P_DDG2P\tLOF\tinclude\n"
        "GENE2\tG2P_DDG2P\t\tinclude\n"
        "GENE3\tG2P_DDG2P\tGOF\treview\n"
        "GENE4\tOrphadata\tGOF\tinclude\n"
        "GENE5\tPanelApp\tGOF\tinclude\n",
        encoding="utf-8",
    )
    aliases = {"OLD1": "GENE1"}

    result = load_condition_mechanism_gene_sets(
        evidence,
        resolve_symbol=lambda symbol: aliases.get(symbol, symbol),
    )

    assert result["G2P_DDG2P"]["all_rows"] == {"GENE1", "GENE2", "GENE3"}
    assert result["G2P_DDG2P"]["canonical_mechanism"] == {"GENE1", "GENE3"}
    assert result["G2P_DDG2P"]["priva_included_mechanism"] == {"GENE1"}
    assert result["Orphadata"]["canonical_mechanism"] == {"GENE4"}
