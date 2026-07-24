from collections import Counter
import json
from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_gofcards_source_coverage import (  # noqa: E402
    DEFAULT_EVIDENCE,
    DEFAULT_GOFCARDS,
    build_gofcards_coverage_audit,
    load_condition_mechanism_gene_sets,
    load_exact_clinvar_linked_gene_counts,
    load_gofcards_gene_counts,
    main,
    parse_args,
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


def test_load_exact_clinvar_linked_gene_counts_requires_exact_allele(
    tmp_path: Path,
) -> None:
    mechanism_json = tmp_path / "mechanisms.json"
    mechanism_json.write_text(
        json.dumps(
            {
                "_meta": {},
                "HGNC:1": {
                    "symbol": "OLD1",
                    "variant_level": [
                        {
                            "ClinVar_VCV": {
                                "match": {
                                    "method": "exact_normalized_vcf_allele",
                                    "matched_gofcards_records": [{"id": "v1"}],
                                }
                            }
                        },
                        {
                            "ClinVar_VCV": {
                                "match": {
                                    "method": "exact_normalized_vcf_allele",
                                    "matched_gofcards_records": [{"id": "v2"}],
                                }
                            }
                        },
                    ],
                },
                "HGNC:2": {
                    "symbol": "GENE2",
                    "variant_level": [
                        {
                            "ClinVar_VCV": {
                                "match": {
                                    "method": "disease_name",
                                    "matched_gofcards_records": [{"id": "v3"}],
                                }
                            }
                        }
                    ],
                },
            }
        ),
        encoding="utf-8",
    )

    counts = load_exact_clinvar_linked_gene_counts(
        mechanism_json,
        resolve_symbol=lambda symbol: "GENE1" if symbol == "OLD1" else symbol,
    )

    assert counts == {"GENE1": 2}


def test_build_gofcards_coverage_audit_exposes_policy_boundaries() -> None:
    audit = build_gofcards_coverage_audit(
        Counter({"GENE1": 2, "GENE2": 1, "GENE3": 1, "GENE4": 1}),
        Counter({"GENE1": 1}),
        {
            "G2P_DDG2P": {
                "all_rows": {"GENE2", "GENE3"},
                "canonical_mechanism": {"GENE2"},
                "priva_included_mechanism": {"GENE2"},
            },
            "Orphadata": {
                "all_rows": {"GENE3"},
                "canonical_mechanism": {"GENE3"},
                "priva_included_mechanism": set(),
            },
        },
    ).set_index("gene_symbol")

    assert audit.loc["GENE1", "only_gofcards_vs_explicit_sources"] == 0
    assert audit.loc["GENE2", "only_gofcards_vs_explicit_sources"] == 0
    assert audit.loc["GENE3", "only_gofcards_vs_explicit_sources"] == 0
    assert audit.loc["GENE3", "only_gofcards_vs_priva_included"] == 1
    assert audit.loc["GENE3", "only_gofcards_vs_any_source_record"] == 0
    assert audit.loc["GENE4", "only_gofcards_vs_explicit_sources"] == 1
    assert audit.loc["GENE4", "only_gofcards_vs_any_source_record"] == 1


def test_parse_args_has_deployed_defaults_and_allows_archives(tmp_path: Path) -> None:
    defaults = parse_args([])
    archived = parse_args(
        [
            "--gofcards",
            str(tmp_path / "old_gofcards.tsv.gz"),
            "--condition-evidence",
            str(tmp_path / "old_evidence.tsv"),
        ]
    )

    assert defaults.gofcards == DEFAULT_GOFCARDS
    assert defaults.condition_evidence == DEFAULT_EVIDENCE
    assert archived.gofcards == tmp_path / "old_gofcards.tsv.gz"
    assert archived.condition_evidence == tmp_path / "old_evidence.tsv"


def test_main_writes_gene_audit_and_summary(tmp_path: Path) -> None:
    gofcards = tmp_path / "gofcards.tsv"
    gofcards.write_text(
        "HGNC_Symbol\tgofcards_variant_id\nGENE1\tv1\nGENE2\tv2\n",
        encoding="utf-8",
    )
    evidence = tmp_path / "evidence.tsv"
    evidence.write_text(
        "gene_symbol\tsource\tnormalized_mechanisms\tpriva_scope\n"
        "GENE1\tG2P_DDG2P\tLOF\tinclude\n",
        encoding="utf-8",
    )
    mechanism_json = tmp_path / "mechanisms.json"
    mechanism_json.write_text(json.dumps({"_meta": {}}), encoding="utf-8")
    hgnc = tmp_path / "hgnc.tsv"
    hgnc.write_text(
        "hgnc_id\tsymbol\tensembl_gene_id\tentrez_id\tprev_symbol\t"
        "alias_symbol\trefseq_accession\tuniprot_ids\tmane_select\n"
        "HGNC:1\tGENE1\tENSG1\t1\t\t\t\t\t\n"
        "HGNC:2\tGENE2\tENSG2\t2\t\t\t\t\t\n",
        encoding="utf-8",
    )
    output = tmp_path / "coverage.tsv"
    summary = tmp_path / "coverage.json"

    exit_code = main(
        [
            "--gofcards",
            str(gofcards),
            "--condition-evidence",
            str(evidence),
            "--mechanism-json",
            str(mechanism_json),
            "--hgnc-table",
            str(hgnc),
            "--output",
            str(output),
            "--summary-output",
            str(summary),
        ]
    )

    audit = pd.read_csv(output, sep="\t")
    summary_payload = json.loads(summary.read_text(encoding="utf-8"))
    assert exit_code == 0
    assert audit.loc[audit["gene_symbol"].eq("GENE2"), "only_gofcards_vs_explicit_sources"].item() == 1
    assert summary_payload["gofcards_exact_genes"] == 2
    assert summary_payload["genes_only_gofcards_vs_explicit_sources"] == 1
