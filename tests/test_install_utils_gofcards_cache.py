"""Tests for the shell-level GoFCards cache validator.

The cache is keyed ``genes -> <SYMBOL> -> variants -> <loc_...._grch37> ->
{record, assemblies -> hg19|hg38 -> {genomic, transcripts}}``. Assembly sits
above everything build-dependent, and the build-independent verdict and
literature live once per variant. These tests pin the guarantees the validator
must enforce on that shape.
"""

from __future__ import annotations

import copy
import gzip
import json
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INSTALL_UTILS = ROOT / "scripts" / "install_utils.sh"

VARIANT_ID = "loc_1:155874747:G->A_grch37"


def _variant() -> dict:
    """A well-formed, runtime-eligible variant present on both builds."""
    return {
        "record": {
            "eligibility": {
                "status": "eligible", "gene_match_status": "gene_concordant",
                "reason": None, "vep_symbol": "RIT1",
            },
            "source": {
                "gofcards_allele_key": "SNV|1|155874747|155874747|G|A",
                "variant_type_label": "SNV", "assembly": "hg19",
                "chrom": "1", "start": "155874747", "ref": "G", "alt": "A",
                "protein_change": "RIT1:NM_006912:exon5:c.G270A:p.M90I",
            },
            "liftover_status": "mapped",
            "annotations": {"rsid": "rs730880487"},
            "evidence": [{"pmid": "23791108", "pscore": "3", "disease": "Noonan syndrome"}],
        },
        "assemblies": {
            "hg19": {
                "genomic": {"chrom": "1", "pos": 155874747, "ref": "G", "alt": "A",
                            "status": "raw_ref_alt"},
                "transcripts": {
                    "ENST00000368323.4": {
                        "by_hgvsc": {"c.270G>A": {
                            "hgvsp": "p.Met90Ile", "consequence": "missense_variant",
                            "canonical": True, "mane_select": "NM_006912.6"}},
                        "by_hgvsp": {"p.Met90Ile": ["c.270G>A"]},
                    }
                },
            },
            "hg38": {
                "genomic": {"chrom": "1", "pos": 155904618, "ref": "G", "alt": "A",
                            "status": "lifted_ref_match"},
                "transcripts": {
                    "ENST00000368323.5": {
                        "by_hgvsc": {"c.270G>A": {
                            "hgvsp": "p.Met90Ile", "consequence": "missense_variant",
                            "canonical": True, "mane_select": "NM_006912.6"}},
                        "by_hgvsp": {"p.Met90Ile": ["c.270G>A"]},
                    }
                },
            },
        },
    }


def _cache(variants: dict | None = None, *, metadata: dict | None = None,
           gene: str = "RIT1") -> dict:
    return {
        "metadata": {
            "source": "GoFCards", "mechanism": "GOF", "derived_on": "2026-07-28",
            **(metadata or {}),
        },
        "genes": {
            gene: {
                "hgnc_id": "HGNC:10023",
                "variants": variants if variants is not None else {VARIANT_ID: _variant()},
            }
        },
    }


def _validate(path: Path, *provenance: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env.update({"PRIVA_INSTALL_UTILS": str(INSTALL_UTILS), "GOFCARDS_CACHE": str(path)})
    return subprocess.run(
        ["bash", "-c",
         'source "$PRIVA_INSTALL_UTILS" >/dev/null; '
         'validate_gofcards_exact_gof_cache "$GOFCARDS_CACHE" "$@"',
         "bash", *provenance],
        check=False, capture_output=True, text=True, env=env,
    )


def _run(tmp_path: Path, cache: dict, name: str = "gofcards.json.gz"):
    path = tmp_path / name
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        json.dump(cache, handle)
    return _validate(path)


def test_validator_accepts_well_formed_cache(tmp_path: Path) -> None:
    result = _run(tmp_path, _cache())

    assert result.returncode == 0, result.stderr
    assert "variants=1" in result.stdout
    assert "eligible=1" in result.stdout
    assert "with_protein_key=1" in result.stdout
    assert "evidence_entries=1" in result.stdout


def test_validator_rejects_cache_without_genes_block(tmp_path: Path) -> None:
    cache = _cache()
    del cache["genes"]

    result = _run(tmp_path, cache)

    assert result.returncode != 0
    assert "no top-level 'genes' block" in result.stderr


def test_validator_rejects_cache_not_declaring_gofcards(tmp_path: Path) -> None:
    result = _run(tmp_path, _cache(metadata={"source": "SomewhereElse"}))

    assert result.returncode != 0
    assert "does not declare GoFCards" in result.stderr


def test_validator_rejects_non_gof_mechanism(tmp_path: Path) -> None:
    result = _run(tmp_path, _cache(metadata={"mechanism": "LOF"}))

    assert result.returncode != 0
    assert "does not declare GOF" in result.stderr


def test_validator_rejects_malformed_variant_identifier(tmp_path: Path) -> None:
    # The identifier is the variant's identity, so a malformed one must never
    # be accepted silently.
    result = _run(tmp_path, _cache(variants={"SNV|1|155874747|G|A": _variant()}))

    assert result.returncode != 0
    assert "not a well-formed variant identifier" in result.stderr


def test_validator_rejects_unknown_eligibility_state(tmp_path: Path) -> None:
    variant = _variant()
    variant["record"]["eligibility"]["status"] = "probably_fine"

    result = _run(tmp_path, _cache(variants={VARIANT_ID: variant}))

    assert result.returncode != 0
    assert "unknown eligibility" in result.stderr


def test_validator_rejects_variant_present_on_no_assembly(tmp_path: Path) -> None:
    variant = _variant()
    variant["assemblies"] = {}

    result = _run(tmp_path, _cache(variants={VARIANT_ID: variant}))

    assert result.returncode != 0
    assert "present on no assembly" in result.stderr


def test_validator_rejects_quarantine_without_a_reason(tmp_path: Path) -> None:
    # A quarantined variant must say why, or the decision cannot be audited.
    variant = _variant()
    variant["record"]["eligibility"].update(
        {"status": "quarantined_gene_discordance", "reason": None,
         "gene_match_status": "gene_discordant", "vep_symbol": "INPP5F"}
    )

    result = _run(tmp_path, _cache(variants={VARIANT_ID: variant}))

    assert result.returncode != 0
    assert "quarantined without a reason" in result.stderr


def test_validator_accepts_reviewed_mechanism_quarantine(tmp_path: Path) -> None:
    # The reviewed mechanism is carried in the state itself, so a reader can
    # tell loss of function from a mixed effect without opening the review table.
    reviewed = copy.deepcopy(_variant())
    reviewed["record"]["eligibility"].update(
        {"status": "quarantined_reviewed_lof", "reason": "article_mechanism_leakage:LOF"}
    )

    result = _run(tmp_path, _cache(variants={
        VARIANT_ID: _variant(),                       # a cache must retain eligible variants
        "loc_1:155874999:C->T_grch37": reviewed,
    }))

    assert result.returncode == 0, result.stderr
    assert "eligible=1" in result.stdout
    assert "quarantined=1" in result.stdout


def test_validator_rejects_assembly_block_without_coordinates(tmp_path: Path) -> None:
    # Absence of an assembly is how a failed liftover is recorded, so an empty
    # assembly block would be ambiguous and must not be accepted.
    variant = _variant()
    variant["assemblies"]["hg38"]["genomic"] = {}

    result = _run(tmp_path, _cache(variants={VARIANT_ID: variant}))

    assert result.returncode != 0
    assert "has no usable coordinates" in result.stderr


def test_validator_rejects_protein_key_pointing_at_unknown_coding_change(tmp_path: Path) -> None:
    # by_hgvsp maps a protein change to the coding changes that produce it, so a
    # dangling reference means the two views disagree.
    variant = _variant()
    variant["assemblies"]["hg19"]["transcripts"]["ENST00000368323.4"]["by_hgvsp"] = {
        "p.Met90Ile": ["c.999G>A"]
    }

    result = _run(tmp_path, _cache(variants={VARIANT_ID: variant}))

    assert result.returncode != 0
    assert "points at unknown" in result.stderr


def test_validator_accepts_variant_that_failed_liftover(tmp_path: Path) -> None:
    # CRLF2 sits in the pseudoautosomal region, whose boundaries moved between
    # builds. It keeps its protein key and simply offers no GRCh38 coordinate.
    variant = copy.deepcopy(_variant())
    variant["record"]["liftover_status"] = "unmapped"
    variant["assemblies"].pop("hg38", None)

    result = _run(tmp_path, _cache(variants={VARIANT_ID: variant}))

    assert result.returncode == 0, result.stderr
    assert "eligible=1" in result.stdout
    assert "on_hg38=0" in result.stdout


def test_validator_rejects_cache_with_no_eligible_variant(tmp_path: Path) -> None:
    variant = _variant()
    variant["record"]["eligibility"].update(
        {"status": "quarantined_gene_discordance",
         "reason": "curated_gene_absent_from_locus",
         "gene_match_status": "gene_discordant", "vep_symbol": "INPP5F"}
    )

    result = _run(tmp_path, _cache(variants={VARIANT_ID: variant}))

    assert result.returncode != 0
    assert "no runtime-eligible variants" in result.stderr


# The re-injection gate asks the validator whether the deployed cache was made
# by the injector, input and settings in use now. Those facts are recorded in a
# nested block, so the check has to reach into it, and it must be able to tell
# "recorded and different" from "never recorded at all".

def _injected(**clinvar: object) -> dict:
    return _cache(metadata={"clinvar": {
        "injector_sha256": "aaa", "source_cache_sha256": "bbb",
        "min_review_stars": 2, **clinvar,
    }})


def _write(tmp_path: Path, cache: dict) -> Path:
    path = tmp_path / "gofcards.json.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        json.dump(cache, handle)
    return path


def test_validator_accepts_matching_nested_provenance(tmp_path: Path) -> None:
    result = _validate(
        _write(tmp_path, _injected()),
        "clinvar.injector_sha256=aaa",
        "clinvar.source_cache_sha256=bbb",
        "clinvar.min_review_stars=2",
    )

    assert result.returncode == 0, result.stderr


def test_validator_rejects_a_changed_nested_value(tmp_path: Path) -> None:
    # A real code change to the injector: the cache must be rebuilt at once
    # rather than waiting out any refresh interval.
    result = _validate(
        _write(tmp_path, _injected()), "clinvar.injector_sha256=a_new_hash"
    )

    assert result.returncode != 0
    assert "clinvar.injector_sha256='aaa'" in result.stderr


def test_validator_rejects_a_cache_that_recorded_no_provenance(tmp_path: Path) -> None:
    # A cache built before provenance was recorded cannot prove it is current,
    # so it is treated as stale rather than assumed good.
    result = _validate(_write(tmp_path, _cache()), "clinvar.injector_sha256=aaa")

    assert result.returncode != 0
    assert "clinvar.injector_sha256=None" in result.stderr


def test_validator_ignores_a_provenance_pair_with_no_expected_value(tmp_path: Path) -> None:
    # An unset configuration value says nothing about fitness, so it is skipped
    # rather than treated as a mismatch.
    result = _validate(_write(tmp_path, _injected()), "clinvar.injector_sha256=")

    assert result.returncode == 0, result.stderr
