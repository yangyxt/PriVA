from __future__ import annotations

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_production_shell_forwards_the_configured_nonlof_cache() -> None:
    shell = _read(SCRIPTS / "prioritization_vcf_per_fam.sh")

    assert 'read_yaml ${config_file} "gene_nonlof_mechanism_json"' in shell
    assert (
        'check_path ${gene_nonlof_mechanism_json} "file" '
        '"gene_nonlof_mechanism_json"'
    ) in shell
    assert "--gene_nonlof_mechanism_json ${gene_nonlof_mechanism_json}" in shell


def test_nonlof_installer_builds_the_path_that_production_consumes() -> None:
    installer = _read(SCRIPTS / "install_utils.sh")
    builder = _read(SCRIPTS / "build_gene_nonlof_mechanism_cache.py")

    assert '--output-json "${cache_json}"' in installer
    assert '"--output-json"' in builder
    assert '"gene_nonlof_mechanism_json" "${cache_json}"' in installer


def test_obsolete_combined_mechanism_cache_is_not_a_production_interface() -> None:
    config = _read(ROOT / "config.yaml")
    acmg = _read(SCRIPTS / "acmg_criteria_assign.py")
    installer = _read(SCRIPTS / "install_utils.sh")
    evidence_builder = _read(
        SCRIPTS / "build_gene_pathogenic_mechanism_cache.py"
    )

    for source in (config, acmg, installer, evidence_builder):
        assert "gene_mechanism_curated_assertions.json" not in source
        assert re.search(r"(?<!nonlof_)gene_mechanism_json", source) is None

    assert "--ddg2p_mechanism_evidence" not in acmg
    assert "ddg2p_mechanism_evidence" in config
    assert "ddg2p_mechanism_evidence" in installer
    assert '"unified_json"' not in evidence_builder
