#!/usr/bin/env python3
"""Audit GoFCards and ClinVar VCV evidence in the canonical non-LOF JSON."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from build_gene_nonlof_mechanism_cache import gofcards_allele_identity


ALLOWED_REVIEW_STATUSES = {
    "criteria provided, multiple submitters, no conflicts": 2,
    "reviewed by expert panel": 3,
    "practice guideline": 4,
}

GOFCARDS_NORMALIZATION_STATUSES = {
    "matched_gene_concordant",
    "gene_discordant_coordinate_collision",
    "quarantined_upstream_gene_discordance",
    "quarantined_upstream_mechanism_review",
    "unmatched_public_source_allele",
}
GOFCARDS_QUARANTINE_ELIGIBILITIES = {
    "quarantined_gene_discordance",
    "quarantined_allele_gene_discordance",
    "quarantined_reviewed_lof",
    "quarantined_reviewed_mixed",
    "quarantined_reviewed_dominant_negative",
    "quarantined_reviewed_uncertain",
    "quarantined_reviewed_exclude",
    "quarantined_mechanism_review_required",
}
GOFCARDS_MECHANISM_QUARANTINE_ELIGIBILITIES = {
    value
    for value in GOFCARDS_QUARANTINE_ELIGIBILITIES
    if value.startswith(("quarantined_reviewed_", "quarantined_mechanism_"))
}

PRIVA_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_NONLOF_ASSERTIONS_JSON = (
    PRIVA_ROOT
    / "data"
    / "gene_pathogenic_mechanism"
    / "prepared"
    / "gene_nonlof_mechanism_curated_assertions.json"
)


def audit_compact_record_provenance(
    record: dict,
    *,
    parent_gene: str,
    expected_eligibility: str,
    label: str,
) -> list[str]:
    """Return violations in one canonicalized compact GoFCards record.

    Eligible rows must preserve the curated source gene and may have either a
    concordant VEP gene or no VEP gene. A directly discordant row must record a
    different VEP gene; an allele-quarantined sibling retains its concordant or
    source-only row provenance. Keeping these checks in one helper makes the
    GoFCards assertion and ClinVar-link audits identical.
    """
    violations: list[str] = []
    exported = str(record.get("HGNC_Symbol", "")).upper()
    source = str(record.get("GoFCards_HGNC_Symbol", "")).upper()
    vep = str(record.get("VEP_HGNC_Symbol", "")).upper()
    status = str(record.get("gene_match_status", "")).lower()
    eligibility = str(record.get("match_eligibility", "")).lower()

    if not exported or exported != parent_gene.upper():
        violations.append(f"{label}: compact gene differs from parent gene")
    if not source or source != exported:
        violations.append(f"{label}: curated source gene was not preserved")
    if eligibility != expected_eligibility:
        violations.append(
            f"{label}: expected eligibility {expected_eligibility!r}, "
            f"observed {eligibility!r}"
        )

    if expected_eligibility == "eligible":
        if status not in {"gene_concordant", "source_gene_only"}:
            violations.append(f"{label}: eligible row has invalid gene status")
        if status == "gene_concordant" and vep != source:
            violations.append(f"{label}: false gene-concordant label")
        if status == "source_gene_only" and vep:
            violations.append(f"{label}: source-gene-only row contains a VEP gene")
        if str(record.get("mechanism", "GOF")).upper() != "GOF":
            violations.append(f"{label}: eligible exact row is not reviewed GOF")
    elif expected_eligibility == "quarantined_gene_discordance":
        if status != "gene_discordant":
            violations.append(f"{label}: quarantined row lacks discordant status")
        if not vep or vep == source:
            violations.append(f"{label}: quarantined row lacks a different VEP gene")
    elif expected_eligibility == "quarantined_allele_gene_discordance":
        if status not in {"gene_concordant", "source_gene_only"}:
            violations.append(
                f"{label}: allele-quarantined sibling has invalid gene status"
            )
        if status == "gene_concordant" and vep != source:
            violations.append(
                f"{label}: allele-quarantined concordant sibling is discordant"
            )
        if status == "source_gene_only" and vep:
            violations.append(
                f"{label}: allele-quarantined source-only sibling has a VEP gene"
            )
    elif expected_eligibility in GOFCARDS_MECHANISM_QUARANTINE_ELIGIBILITIES:
        if status not in {"gene_concordant", "source_gene_only"}:
            violations.append(
                f"{label}: mechanism-quarantined row has invalid gene status"
            )
        if status == "gene_concordant" and vep != source:
            violations.append(
                f"{label}: mechanism-quarantined concordant row is gene-discordant"
            )
        if status == "source_gene_only" and vep:
            violations.append(
                f"{label}: mechanism-quarantined source-only row contains a VEP gene"
            )
        mechanism_eligibility = str(
            record.get("mechanism_match_eligibility", "")
        ).lower()
        if mechanism_eligibility != expected_eligibility:
            violations.append(
                f"{label}: mechanism eligibility does not match quarantine reason"
            )
        if expected_eligibility.startswith("quarantined_reviewed_") and (
            str(record.get("mechanism_review_status", "")).lower() != "reviewed"
        ):
            violations.append(
                f"{label}: reviewed mechanism quarantine lacks review status"
            )
    else:
        violations.append(f"{label}: unsupported eligibility {expected_eligibility!r}")
    return violations


def audit_compact_record_allele_identity(
    record: dict,
    *,
    parent_allele_key: str,
    label: str,
) -> list[str]:
    """Return a violation unless parent and compact keys identify one allele.

    The GoFCards backend calls some multi-base substitutions ``SNV`` while the
    public workbook parser derives ``Indel`` from allele length. The canonical
    builder intentionally joins those records by the type-independent
    chromosome/start/end/ref/alt identity. The audit must enforce that same
    contract instead of treating the source vocabulary prefix as identity.
    """
    nested_allele_key = str(record.get("allele_key", ""))
    if (
        not nested_allele_key
        or gofcards_allele_identity(nested_allele_key)
        != gofcards_allele_identity(parent_allele_key)
    ):
        return [f"{label}: nested allele identity disagrees"]
    return []


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--nonlof-assertions-json",
        "--json",
        dest="nonlof_assertions_json",
        type=Path,
        default=DEFAULT_NONLOF_ASSERTIONS_JSON,
        help="Canonical non-LOF assertions JSON; --json is a compatibility alias.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    canonical = json.loads(
        args.nonlof_assertions_json.read_text(encoding="utf-8")
    )
    gene_count = sum(key != "_meta" for key in canonical)
    entries = 0
    vcvs: set[str] = set()
    rcvs: set[tuple[str, str]] = set()
    scv_attachments = 0
    scvs: set[tuple[str, str]] = set()
    review_stars: Counter[int] = Counter()
    scopes: Counter[str] = Counter()
    contexts: Counter[str] = Counter()
    zygosity_observations: Counter[str] = Counter()
    with_moi = 0
    with_penetrance = 0
    ambiguous_scv_attachments = 0
    discarded_gene_symbols = 0
    gofcards_assertions = 0
    gofcards_statuses: Counter[str] = Counter()
    gofcards_exact_records = 0
    gofcards_exact_record_keys: set[str] = set()
    gofcards_quarantined_records = 0
    gofcards_quarantined_record_keys: set[str] = set()
    gofcards_with_hg19_vcf = 0
    gofcards_with_hg38_vcf = 0
    gofcards_exact_assertions: set[tuple[str, str]] = set()
    gofcards_quarantined_assertions: set[tuple[str, str]] = set()
    gofcards_clinvar_match_assertions: set[tuple[str, str]] = set()
    violations: list[str] = []

    for gene_key, gene in canonical.items():
        if gene_key == "_meta":
            continue
        gene_symbol = str(gene.get("symbol", ""))
        for assertion in gene.get("variant_level", []):
            gofcards = assertion.get("GoFCards")
            if isinstance(gofcards, dict):
                gofcards_assertions += 1
                source_record_id = str(gofcards.get("source_record_id", ""))
                allele_key = str(gofcards.get("allele_key", ""))
                status = str(gofcards.get("exact_normalization_status", ""))
                gofcards_statuses[status] += 1
                label = f"{gene_symbol}/GoFCards:{source_record_id or allele_key}"
                if not source_record_id or not allele_key:
                    violations.append(f"{label}: missing source record or allele key")
                if status not in GOFCARDS_NORMALIZATION_STATUSES:
                    violations.append(f"{label}: invalid normalization status {status!r}")

                exact_records = gofcards.get("exact_normalized_variants", [])
                if status == "matched_gene_concordant" and not exact_records:
                    violations.append(f"{label}: matched status without normalized records")
                if status != "matched_gene_concordant" and exact_records:
                    violations.append(f"{label}: non-matched status has normalized records")
                quarantined_records = gofcards.get(
                    "quarantined_exact_normalized_variants", []
                )
                if (
                    status in {
                        "quarantined_upstream_gene_discordance",
                        "quarantined_upstream_mechanism_review",
                    }
                    and not quarantined_records
                ):
                    violations.append(f"{label}: quarantine status without audit records")
                collision_symbols = gofcards.get("exact_cache_gene_symbols", [])
                if (
                    status in {
                        "gene_discordant_coordinate_collision",
                        "quarantined_upstream_gene_discordance",
                    }
                    and not collision_symbols
                ):
                    violations.append(f"{label}: gene collision lacks cache symbols")

                has_hg19 = False
                has_hg38 = False
                for record in exact_records:
                    if not isinstance(record, dict):
                        violations.append(f"{label}: non-object normalized record")
                        continue
                    gofcards_exact_records += 1
                    gofcards_exact_record_keys.update(record)
                    gofcards_exact_assertions.add(
                        (
                            gene_symbol.upper(),
                            gofcards_allele_identity(
                                record.get("gofcards_variant_id")
                                or record.get("allele_key")
                            ),
                        )
                    )
                    violations.extend(
                        audit_compact_record_provenance(
                            record,
                            parent_gene=gene_symbol,
                            expected_eligibility="eligible",
                            label=label,
                        )
                    )
                    violations.extend(
                        audit_compact_record_allele_identity(
                            record,
                            parent_allele_key=allele_key,
                            label=label,
                        )
                    )
                    has_hg19 = has_hg19 or bool(record.get("hg19_vcf_key"))
                    has_hg38 = has_hg38 or bool(record.get("hg38_vcf_key"))
                if has_hg19:
                    gofcards_with_hg19_vcf += 1
                if has_hg38:
                    gofcards_with_hg38_vcf += 1
                if status == "matched_gene_concordant" and not has_hg19:
                    violations.append(f"{label}: no normalized hg19 VCF key")
                if status == "matched_gene_concordant" and not has_hg38:
                    violations.append(f"{label}: no normalized hg38 VCF key")

                for record in quarantined_records:
                    if not isinstance(record, dict):
                        violations.append(f"{label}: non-object quarantine record")
                        continue
                    gofcards_quarantined_records += 1
                    gofcards_quarantined_record_keys.update(record)
                    gofcards_quarantined_assertions.add(
                        (
                            gene_symbol.upper(),
                            gofcards_allele_identity(
                                record.get("gofcards_variant_id")
                                or record.get("allele_key")
                            ),
                        )
                    )
                    record_eligibility = str(
                        record.get("match_eligibility", "")
                    ).lower()
                    if record_eligibility not in GOFCARDS_QUARANTINE_ELIGIBILITIES:
                        violations.append(
                            f"{label}: audit-only record is not quarantined"
                        )
                    else:
                        violations.extend(
                            audit_compact_record_provenance(
                                record,
                                parent_gene=gene_symbol,
                                expected_eligibility=record_eligibility,
                                label=label,
                            )
                        )
                    violations.extend(
                        audit_compact_record_allele_identity(
                            record,
                            parent_allele_key=allele_key,
                            label=label,
                        )
                    )

            payload = assertion.get("ClinVar_VCV")
            if not isinstance(payload, dict):
                continue
            entries += 1
            variation = payload.get("variation", {})
            vcv = str(variation.get("vcv_accession", ""))
            vcvs.add(vcv)
            scopes[str(variation.get("classification_scope", ""))] += 1
            contexts[str(variation.get("matched_component_context", ""))] += 1
            if variation.get("record_type") != "classified":
                violations.append(f"{gene_symbol}/{vcv}: non-classified record")
            if variation.get("record_status") != "current":
                violations.append(f"{gene_symbol}/{vcv}: non-current VCV")
            if payload.get("allelic_requirement", {}).get("value") != "unresolved":
                violations.append(f"{gene_symbol}/{vcv}: resolved allelic requirement")

            match = payload.get("match", {})
            clinvar_symbols = set(match.get("clinvar_gene_symbols", []))
            discarded_gene_symbols += len(
                match.get("discarded_gene_discordant_symbols", [])
            )
            for row in match.get("matched_gofcards_records", []):
                symbol = str(row.get("HGNC_Symbol", ""))
                gofcards_clinvar_match_assertions.add(
                    (
                        gene_symbol.upper(),
                        gofcards_allele_identity(
                            row.get("gofcards_variant_id")
                            or row.get("allele_key")
                        ),
                    )
                )
                violations.extend(
                    audit_compact_record_provenance(
                        row,
                        parent_gene=gene_symbol,
                        expected_eligibility="eligible",
                        label=f"{gene_symbol}/{vcv}",
                    )
                )
                if clinvar_symbols and symbol not in clinvar_symbols:
                    violations.append(
                        f"{gene_symbol}/{vcv}: retained discordant GoFCards gene {symbol}"
                    )

            for condition in payload.get("condition_assertions", []):
                rcv = str(condition.get("rcv_accession", ""))
                rcvs.add((vcv, rcv))
                classification = condition.get("germline_classification", {})
                status = str(classification.get("review_status", ""))
                stars = classification.get("review_stars")
                if ALLOWED_REVIEW_STATUSES.get(status) != stars:
                    violations.append(
                        f"{gene_symbol}/{vcv}/{rcv}: invalid review status/stars"
                    )
                if isinstance(stars, int):
                    review_stars[stars] += 1
                for scv in condition.get("matched_scvs", []):
                    scv_attachments += 1
                    scv_accession = str(scv.get("scv_accession", ""))
                    scvs.add((vcv, scv_accession))
                    if scv.get("record_status") != "current":
                        violations.append(
                            f"{gene_symbol}/{vcv}/{rcv}/{scv_accession}: non-current SCV"
                        )
                    if scv.get("contributes_to_aggregate_classification") is not True:
                        violations.append(
                            f"{gene_symbol}/{vcv}/{rcv}/{scv_accession}: non-contributing SCV"
                        )
                    if scv.get("submitted_mode_of_inheritance"):
                        with_moi += 1
                    if scv.get("penetrance"):
                        with_penetrance += 1
                    if scv.get("trait_linkage_ambiguous_across_eligible_rcvs"):
                        ambiguous_scv_attachments += 1
                    zygosity_observations.update(
                        {
                            str(key): int(value)
                            for key, value in scv.get(
                                "observed_zygosity_counts", {}
                            ).items()
                        }
                    )

    expected_gene_count = canonical.get("_meta", {}).get("total_genes")
    if expected_gene_count != gene_count:
        violations.append(
            f"meta total_genes={expected_gene_count}, observed={gene_count}"
        )

    exact_quarantine_overlap = (
        gofcards_exact_assertions & gofcards_quarantined_assertions
    ) - {("", "")}
    clinvar_quarantine_overlap = (
        gofcards_clinvar_match_assertions & gofcards_quarantined_assertions
    ) - {("", "")}
    if exact_quarantine_overlap:
        violations.append(
            "quarantined gene-allele assertions leaked into canonical exact "
            "records: "
            + ", ".join(
                f"{gene}:{allele}"
                for gene, allele in sorted(exact_quarantine_overlap)
            )
        )
    if clinvar_quarantine_overlap:
        violations.append(
            "quarantined gene-allele assertions leaked into ClinVar matches: "
            + ", ".join(
                f"{gene}:{allele}"
                for gene, allele in sorted(clinvar_quarantine_overlap)
            )
        )

    gofcards_meta = (
        canonical.get("_meta", {}).get("sources", {}).get("GoFCards", {})
    )
    expected_raw_fields = {
        "assembly": "hg19",
        "keys": ["chr", "pos", "ref", "alt"],
        "note": (
            "Applies only to the legacy top-level GoFCards source fields; "
            "exact_normalized_variants identify hg19 and hg38 fields explicitly"
        ),
    }
    if "assembly" in gofcards_meta:
        violations.append("GoFCards metadata has ambiguous top-level assembly")
    if gofcards_meta.get("raw_public_allele_fields") != expected_raw_fields:
        violations.append("GoFCards raw public allele metadata is missing or invalid")

    print(f"canonical_nonlof_json={args.nonlof_assertions_json.resolve()}")
    print(f"schema_version={canonical.get('_meta', {}).get('version', '')}")
    print(f"canonical_genes={gene_count}")
    print(f"gofcards_assertions={gofcards_assertions}")
    print(
        "gofcards_normalization_statuses="
        f"{json.dumps(dict(sorted(gofcards_statuses.items())))}"
    )
    print(f"gofcards_exact_normalized_records={gofcards_exact_records}")
    print(f"gofcards_quarantined_records={gofcards_quarantined_records}")
    print(
        "gofcards_exact_unique_gene_alleles="
        f"{len(gofcards_exact_assertions - {('', '')})}"
    )
    print(
        "gofcards_quarantined_unique_gene_alleles="
        f"{len(gofcards_quarantined_assertions - {('', '')})}"
    )
    print(
        "gofcards_exact_quarantine_overlap="
        f"{len(exact_quarantine_overlap)}"
    )
    print(
        "gofcards_clinvar_quarantine_overlap="
        f"{len(clinvar_quarantine_overlap)}"
    )
    print(f"gofcards_matched_assertions_with_hg19_vcf={gofcards_with_hg19_vcf}")
    print(f"gofcards_matched_assertions_with_hg38_vcf={gofcards_with_hg38_vcf}")
    print(
        "gofcards_exact_normalized_record_keys="
        f"{json.dumps(sorted(gofcards_exact_record_keys))}"
    )
    print(
        "gofcards_quarantined_record_keys="
        f"{json.dumps(sorted(gofcards_quarantined_record_keys))}"
    )
    print(f"clinvar_vcv_entries={entries}")
    print(f"unique_vcv_accessions={len(vcvs)}")
    print(f"unique_vcv_rcv_pairs={len(rcvs)}")
    print(f"rcv_review_stars={json.dumps(dict(sorted(review_stars.items()))) }")
    print(f"classification_scopes={json.dumps(dict(sorted(scopes.items()))) }")
    print(f"component_contexts={json.dumps(dict(sorted(contexts.items()))) }")
    print(f"scv_attachments={scv_attachments}")
    print(f"unique_vcv_scv_pairs={len(scvs)}")
    print(
        "zygosity_observation_counts="
        f"{json.dumps(dict(sorted(zygosity_observations.items())))}"
    )
    print(f"scv_attachments_with_submitted_moi={with_moi}")
    print(f"scv_attachments_with_penetrance={with_penetrance}")
    print(f"ambiguous_scv_attachments={ambiguous_scv_attachments}")
    print(f"discarded_gene_discordant_symbols={discarded_gene_symbols}")
    print(f"violations={len(violations)}")
    for violation in violations[:20]:
        print(f"violation={violation}")
    if violations:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
