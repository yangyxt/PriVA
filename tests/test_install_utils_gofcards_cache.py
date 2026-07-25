from __future__ import annotations

import gzip
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INSTALL_UTILS = ROOT / "scripts" / "install_utils.sh"

PROVENANCE_COLUMNS = [
    "source",
    "mechanism",
    "HGNC_Symbol",
    "GoFCards_HGNC_Symbol",
    "VEP_HGNC_Symbol",
    "gene_match_status",
    "match_eligibility",
    "HGVSp",
    "hgvsp_key",
    "hg19_vcf_key",
    "hg38_vcf_key",
    "hg19_genomic_key",
    "hg38_genomic_key",
    "gofcards_accession_id",
    "gofcards_variant_id",
]


def _write_cache(path: Path, row: dict[str, str], *, include_provenance: bool = True) -> None:
    columns = list(PROVENANCE_COLUMNS)
    if not include_provenance:
        columns = [
            column
            for column in columns
            if column
            not in {
                "GoFCards_HGNC_Symbol",
                "VEP_HGNC_Symbol",
                "gene_match_status",
                "match_eligibility",
            }
        ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(columns) + "\n")
        handle.write("\t".join(row.get(column, "") for column in columns) + "\n")


def _validate_cache(path: Path) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env.update(
        {
            "PRIVA_INSTALL_UTILS": str(INSTALL_UTILS),
            "GOFCARDS_CACHE": str(path),
        }
    )
    return subprocess.run(
        [
            "bash",
            "-c",
            'source "$PRIVA_INSTALL_UTILS" >/dev/null; '
            'validate_gofcards_exact_gof_cache "$GOFCARDS_CACHE"',
        ],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )


def _eligible_row() -> dict[str, str]:
    return {
        "source": "GoFCards",
        "mechanism": "GOF",
        "HGNC_Symbol": "RIT1",
        "GoFCards_HGNC_Symbol": "RIT1",
        "VEP_HGNC_Symbol": "RIT1",
        "gene_match_status": "gene_concordant",
        "match_eligibility": "eligible",
        "HGVSp": "ENSP00000499386.2:p.Met90Ile",
        "hgvsp_key": "M90I",
        "hg19_vcf_key": "1|155874747|G|A",
        "hg38_vcf_key": "1|155904618|G|A",
        "hg19_genomic_key": "1|155874747|G|A",
        "hg38_genomic_key": "1|155904618|G|A",
        "gofcards_accession_id": "rs730880487",
        "gofcards_variant_id": "SNV|1|155874747|155874747|G|A",
    }


def test_validator_accepts_source_preserving_compact_cache(tmp_path: Path) -> None:
    cache = tmp_path / "gofcards.tsv.gz"
    _write_cache(cache, _eligible_row())

    result = _validate_cache(cache)

    assert result.returncode == 0, result.stderr
    assert "eligible_rows=1" in result.stdout
    assert "quarantined_rows=0" in result.stdout


def test_validator_rejects_legacy_cache_without_provenance(tmp_path: Path) -> None:
    cache = tmp_path / "legacy.tsv.gz"
    _write_cache(cache, _eligible_row(), include_provenance=False)

    result = _validate_cache(cache)

    assert result.returncode != 0
    assert "missing required columns" in result.stderr


def test_validator_rejects_gene_discordance_marked_eligible(tmp_path: Path) -> None:
    cache = tmp_path / "unsafe.tsv.gz"
    row = _eligible_row()
    row.update(
        {
            "VEP_HGNC_Symbol": "INPP5F",
            "gene_match_status": "gene_discordant",
        }
    )
    _write_cache(cache, row)

    result = _validate_cache(cache)

    assert result.returncode != 0
    assert "invalid gene provenance" in result.stderr
