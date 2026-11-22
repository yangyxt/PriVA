from __future__ import annotations

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_production_config_names_both_step3_runtime_caches() -> None:
    config = _read(ROOT / "config.yaml")

    assert re.search(
        r"^gene_nonlof_mechanism_json: .*/"
        r"gene_nonlof_mechanism_curated_assertions\.json\.gz$",
        config,
        re.MULTILINE,
    )
    assert re.search(
        r"^hpo_condition_mechanism_json: .*/"
        r"hpo_condition_mechanism_cache\.json\.gz$",
        config,
        re.MULTILINE,
    )


