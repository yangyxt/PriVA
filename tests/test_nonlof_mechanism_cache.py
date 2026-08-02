import argparse
import json
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import build_gene_nonlof_mechanism_cache as nonlof  # noqa: E402
from clinvar_vcv import gofcards_variant_is_eligible  # noqa: E402


def test_portable_defaults_stay_inside_priva() -> None:
    expected_cache = ROOT / "data" / "patho_mechanism"

    assert nonlof.DEFAULT_CACHE_DIR == expected_cache
    assert nonlof.DEFAULT_SHARED_RAW_DIR == expected_cache / "raw"
    # Compressed on disk: the cache is tracked in git and gzip takes it from
    # 96.6 MiB to 5.7 MiB. Readers open it through the shared opener, which
    # decides by the .gz suffix, so nothing else changes.
    assert nonlof.DEFAULT_NONLOF_ASSERTIONS_JSON == (
        expected_cache / "gene_nonlof_mechanism_curated_assertions.json.gz"
    )
    assert nonlof.DEFAULT_GOFCARDS_EXACT_VARIANTS == (
        ROOT / "data" / "gofcards" / "gofcards_exact_gof.json.gz"
    )
    assert nonlof.DEFAULT_HGNC_TABLE == (
        ROOT / "data" / "hgnc" / "non_alt_loci_set.tsv"
    )
    assert "llm_gene_reranker" not in str(nonlof.DEFAULT_CACHE_DIR)


def test_nonlof_outputs_do_not_replace_broad_cache_manifests() -> None:
    assert nonlof.NONLOF_SOURCE_MANIFEST_FILENAME == "nonlof_source_manifest.json"
    assert (
        nonlof.NONLOF_SOURCE_MANIFEST_TSV_FILENAME
        == "nonlof_source_manifest.tsv"
    )
    assert nonlof.NONLOF_RUN_SUMMARY_FILENAME == "nonlof_run_summary.json"


def _variant_with_status(status: str) -> dict:
    """The smallest variant shape the quarantine reader looks at."""
    return {"record": {"eligibility": {"status": status}}}


def test_reviewed_non_gof_variant_is_never_runtime_eligible() -> None:
    # Eligibility is decided once by the normalizer and stored on the variant.
    # Every quarantine state, whatever its cause, fails the gate.
    for status in (
        "quarantined_reviewed_lof",
        "quarantined_reviewed_mixed",
        "quarantined_gene_discordance",
    ):
        assert not gofcards_variant_is_eligible(_variant_with_status(status))
    assert gofcards_variant_is_eligible(_variant_with_status("eligible"))


def test_upstream_quarantine_status_preserves_the_failed_gate() -> None:
    assert nonlof.gofcards_upstream_quarantine_status(
        _variant_with_status("quarantined_reviewed_lof")
    ) == "quarantined_upstream_mechanism_review"
    assert nonlof.gofcards_upstream_quarantine_status(
        _variant_with_status("quarantined_gene_discordance")
    ) == "quarantined_upstream_gene_discordance"


def test_hgnc_loader_accepts_priva_complete_set_columns(tmp_path: Path) -> None:
    table = tmp_path / "hgnc.tsv"
    table.write_text(
        "hgnc_id\tsymbol\talias_symbol\tprev_symbol\tentrez_id\t"
        "ensembl_gene_id\n"
        "HGNC:4764\tH3-3A\tH3.3A\tH3F3|H3F3A\t3020.0\t"
        "ENSG00000163041\n",
        encoding="utf-8",
    )

    mapping = nonlof.load_hgnc_mapping(table)

    assert mapping["H3-3A"]["hgnc_id"] == "HGNC:4764"
    assert mapping["H3.3A"]["symbol"] == "H3-3A"
    assert mapping["H3F3A"]["entrez_id"] == "3020"


def test_deployed_nonlof_json_validates_against_priva_schema() -> None:
    canonical_path = nonlof.DEFAULT_NONLOF_ASSERTIONS_JSON
    schema_path = nonlof.DEFAULT_OUTPUT_SCHEMA
    canonical = nonlof.read_json(canonical_path)

    nonlof.validate_canonical_json(canonical, schema_path)
    assert canonical["_meta"]["version"] == "2.0"
    assert canonical["_meta"]["total_genes"] == len(canonical) - 1


def test_failed_rebuild_does_not_advance_source_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cache_dir = tmp_path / "cache"
    raw_dir = cache_dir / "raw"
    manifest_path = cache_dir / nonlof.NONLOF_SOURCE_MANIFEST_FILENAME
    previous_manifest = {
        "schema_version": "2.0",
        "marker": "last_successful_build",
        "build_inputs": {"builder_sha256": "previous"},
    }
    manifest_path.parent.mkdir(parents=True)
    manifest_path.write_text(json.dumps(previous_manifest), encoding="utf-8")

    exact_gofcards = tmp_path / "gofcards.tsv"
    hgnc = tmp_path / "hgnc.tsv"
    schema = tmp_path / "schema.json"
    for path in (exact_gofcards, hgnc, schema):
        path.write_text("fixture\n", encoding="utf-8")

    args = argparse.Namespace(
        validate_only=None,
        cache_dir=cache_dir,
        output_json=None,
        shared_raw_dir=str(raw_dir),
        gofcards_exact_variants=exact_gofcards,
        hgnc_table=hgnc,
        output_schema=schema,
        clinvar_min_review_stars=2,
        clinvar_max_download_seconds=60,
        clinvar_only_refresh=False,
        fetch_sources_only=False,
        force=False,
        max_panelapp_panels=None,
        timeout=1,
        retries=1,
        proxy_url="",
        download_tool="auto",
        stale_lock_hours=1.0,
    )
    source_meta = {"sha256": "fixture", "status": "skipped_fresh"}
    monkeypatch.setattr(nonlof, "parse_args", lambda: args)
    monkeypatch.setattr(
        nonlof, "fetch_static_source", lambda **_kwargs: dict(source_meta)
    )
    monkeypatch.setattr(
        nonlof, "fetch_panelapp", lambda **_kwargs: dict(source_meta)
    )
    monkeypatch.setattr(
        nonlof, "fetch_clinvar_vcv", lambda **_kwargs: dict(source_meta)
    )

    def fail_during_rebuild(**_kwargs: object) -> None:
        raise RuntimeError("simulated parse failure")

    monkeypatch.setattr(nonlof, "parse_all_sources", fail_during_rebuild)

    with pytest.raises(RuntimeError, match="simulated parse failure"):
        nonlof.main()

    assert json.loads(manifest_path.read_text(encoding="utf-8")) == previous_manifest


def test_fetch_sources_only_stops_before_the_build(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The installer runs this phase before ClinVar injection exists, so the
    # GoFCards cache the ordinary build reads is still missing here. The phase
    # must place the raw sources without touching that file, and must leave the
    # previous build_inputs alone so the next full run can still tell whether a
    # rebuild is needed.
    cache_dir = tmp_path / "cache"
    raw_dir = cache_dir / "raw"
    manifest_path = cache_dir / nonlof.NONLOF_SOURCE_MANIFEST_FILENAME
    previous_build_inputs = {"builder_sha256": "previous"}
    manifest_path.parent.mkdir(parents=True)
    manifest_path.write_text(
        json.dumps({"schema_version": "2.0", "build_inputs": previous_build_inputs}),
        encoding="utf-8",
    )

    missing_gofcards = tmp_path / "gofcards_not_written_yet.json.gz"
    hgnc = tmp_path / "hgnc.tsv"
    schema = tmp_path / "schema.json"
    for path in (hgnc, schema):
        path.write_text("fixture\n", encoding="utf-8")

    args = argparse.Namespace(
        validate_only=None,
        cache_dir=cache_dir,
        output_json=None,
        shared_raw_dir=str(raw_dir),
        gofcards_exact_variants=missing_gofcards,
        hgnc_table=hgnc,
        output_schema=schema,
        clinvar_min_review_stars=2,
        clinvar_max_download_seconds=60,
        clinvar_only_refresh=False,
        fetch_sources_only=True,
        force=False,
        max_panelapp_panels=None,
        timeout=1,
        retries=1,
        proxy_url="",
        download_tool="auto",
        stale_lock_hours=1.0,
    )
    source_meta = {"sha256": "fixture", "status": "skipped_fresh"}
    monkeypatch.setattr(nonlof, "parse_args", lambda: args)
    monkeypatch.setattr(
        nonlof, "fetch_static_source", lambda **_kwargs: dict(source_meta)
    )
    monkeypatch.setattr(nonlof, "fetch_panelapp", lambda **_kwargs: dict(source_meta))
    monkeypatch.setattr(
        nonlof, "fetch_clinvar_vcv", lambda **_kwargs: dict(source_meta)
    )

    def fail_if_reached(**_kwargs: object) -> None:
        raise AssertionError("the build ran during a sources-only phase")

    monkeypatch.setattr(nonlof, "parse_all_sources", fail_if_reached)

    assert nonlof.main() == 0

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    # Downloads recorded, so the full run afterwards skips them.
    assert set(manifest["sources"]) == {
        source.name for source in nonlof.SOURCES
    } | {nonlof.PANELAPP_SPEC.name, nonlof.CLINVAR_VCV_SPEC.name}
    # Carried forward, not replaced by a partial record.
    assert manifest["build_inputs"] == previous_build_inputs
    # No cache published, and the missing GoFCards input was never required.
    assert not (cache_dir / nonlof.NONLOF_ASSERTIONS_FILENAME).exists()
    assert not missing_gofcards.exists()
