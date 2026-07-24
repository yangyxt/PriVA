import sys
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from collapse_HPO_anno import (  # noqa: E402
    ASSERTION_COLUMNS,
    build_gene_disease_hpo_assertions,
    write_gene_disease_hpo_assertions,
)


def _write_sources(tmp_path: Path) -> tuple[Path, Path]:
    genes = pd.DataFrame(
        [
            {
                "ncbi_gene_id": "1",
                "gene_symbol": "GENE1",
                "hpo_id": "HP:0003829",
                "hpo_name": "Typified by incomplete penetrance",
                "frequency": "9/9",
                "disease_id": "OMIM:1",
            },
            {
                "ncbi_gene_id": "1",
                "gene_symbol": "GENE1",
                "hpo_id": "HP:0003829",
                "hpo_name": "Typified by incomplete penetrance",
                "frequency": "-",
                "disease_id": "OMIM:2",
            },
            {
                "ncbi_gene_id": "2",
                "gene_symbol": "GENE2",
                "hpo_id": "HP:0001250",
                "hpo_name": "Seizure",
                "frequency": "-",
                "disease_id": "OMIM:3",
            },
        ]
    )
    hpoa = pd.DataFrame(
        [
            {
                "database_id": "OMIM:1",
                "disease_name": "Disease one",
                "qualifier": "",
                "hpo_id": "HP:0003829",
                "reference": "PMID:1",
                "evidence": "PCS",
                "onset": "",
                "frequency": "1/3",
                "sex": "",
                "modifier": "",
                "aspect": "I",
                "biocuration": "HPO:test",
            },
            {
                "database_id": "OMIM:1",
                "disease_name": "Disease one",
                "qualifier": "",
                "hpo_id": "HP:0003829",
                "reference": "PMID:2",
                "evidence": "PCS",
                "onset": "",
                "frequency": "2/4",
                "sex": "",
                "modifier": "",
                "aspect": "I",
                "biocuration": "HPO:test",
            },
            {
                "database_id": "OMIM:1",
                "disease_name": "Disease one",
                "qualifier": "NOT",
                "hpo_id": "HP:0003829",
                "reference": "PMID:NEGATED",
                "evidence": "PCS",
                "onset": "",
                "frequency": "0/3",
                "sex": "",
                "modifier": "",
                "aspect": "I",
                "biocuration": "HPO:test",
            },
            {
                "database_id": "OMIM:2",
                "disease_name": "Disease two",
                "qualifier": "",
                "hpo_id": "HP:0003829",
                "reference": "OMIM:2",
                "evidence": "TAS",
                "onset": "",
                "frequency": "",
                "sex": "",
                "modifier": "",
                "aspect": "I",
                "biocuration": "HPO:test",
            },
            {
                "database_id": "OMIM:3",
                "disease_name": "Disease three",
                "qualifier": "NOT",
                "hpo_id": "HP:0001250",
                "reference": "OMIM:3",
                "evidence": "TAS",
                "onset": "",
                "frequency": "",
                "sex": "",
                "modifier": "",
                "aspect": "P",
                "biocuration": "HPO:test",
            },
        ]
    )
    genes_path = tmp_path / "genes_to_phenotype.txt"
    hpoa_path = tmp_path / "phenotype.hpoa"
    genes.to_csv(genes_path, sep="\t", index=False)
    with hpoa_path.open("w", encoding="ascii") as handle:
        handle.write("#version: test\n")
        hpoa.to_csv(handle, sep="\t", index=False)
    return genes_path, hpoa_path


def test_builds_one_row_per_full_assertion(tmp_path: Path) -> None:
    genes_path, hpoa_path = _write_sources(tmp_path)

    result = build_gene_disease_hpo_assertions(genes_path, hpoa_path)

    assert list(result.columns) == ASSERTION_COLUMNS
    assert len(result) == 3
    assert set(result["disease_id"]) == {"OMIM:1", "OMIM:2"}
    assert set(result.loc[result["disease_id"].eq("OMIM:1"), "frequency"]) == {
        "1/3",
        "2/4",
    }
    assert "PMID:NEGATED" not in set(result["reference"])
    omim_two = result.loc[result["disease_id"].eq("OMIM:2")].iloc[0]
    assert omim_two["frequency"] == "-"
    assert omim_two["evidence"] == "TAS"


def test_writes_gzip_assertion_table(tmp_path: Path) -> None:
    genes_path, hpoa_path = _write_sources(tmp_path)
    output = tmp_path / "assertions.tsv.gz"

    write_gene_disease_hpo_assertions(genes_path, hpoa_path, output)

    written = pd.read_csv(output, sep="\t", dtype=str, keep_default_na=False)
    assert list(written.columns) == ASSERTION_COLUMNS
    assert len(written) == 3
