import csv
import sys

evidence_tsv = sys.argv[1]

required = {
    "gene_symbol",
    "source",
    "source_record_id",
    "source_condition_id",
    "mondo_id",
    "disease_scope",
    "priva_scope",
    "scope_review_status",
    "disease_label",
    "inheritance",
    "penetrance_raw",
    "penetrance_hpo_ids",
    "normalized_penetrance",
    "patho_mode_raw",
    "normalized_mechanisms",
    "mechanism_confidence",
    "disease_confidence",
    "pmids",
}
row_count = 0
g2p_rows = 0
g2p_lof_rows = 0
g2p_penetrance_rows = 0
orphadata_rows = 0
orphadata_mondo_rows = 0
unexpected_sources = set()
with open(evidence_tsv, encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    missing = sorted(required - set(reader.fieldnames or []))
    if missing:
        raise SystemExit(f"{evidence_tsv} missing required columns: {', '.join(missing)}")
    for row in reader:
        row_count += 1
        source = row.get("source", "")
        if source not in {"G2P_DDG2P", "Orphadata"}:
            unexpected_sources.add(source or "<blank>")
        if source == "G2P_DDG2P":
            g2p_rows += 1
            normalized = row.get("normalized_mechanisms", "").upper()
            raw = row.get("patho_mode_raw", "").lower()
            if "LOF" in normalized or "loss of function" in raw:
                g2p_lof_rows += 1
            penetrance = row.get("normalized_penetrance", "").strip().lower()
            if penetrance:
                if penetrance not in {"complete", "incomplete"}:
                    raise SystemExit(
                        f"{evidence_tsv} has unsupported penetrance {penetrance!r}"
                    )
                if not (row.get("source_condition_id") or row.get("mondo_id")):
                    raise SystemExit(
                        f"{evidence_tsv} has penetrance without a condition identifier"
                    )
                g2p_penetrance_rows += 1
        if source == "Orphadata":
            orphadata_rows += 1
            if row.get("mondo_id"):
                orphadata_mondo_rows += 1
if row_count == 0:
    raise SystemExit(f"{evidence_tsv} contains no evidence rows")
if unexpected_sources:
    raise SystemExit(
        f"{evidence_tsv} contains non-condition sources: "
        + ", ".join(sorted(unexpected_sources))
    )
if g2p_rows == 0:
    raise SystemExit(f"{evidence_tsv} contains no G2P_DDG2P rows")
if g2p_lof_rows == 0:
    raise SystemExit(f"{evidence_tsv} contains no DDG2P/G2P LoF rows")
if g2p_penetrance_rows == 0:
    raise SystemExit(f"{evidence_tsv} contains no condition-linked G2P penetrance rows")
if orphadata_rows == 0:
    raise SystemExit(f"{evidence_tsv} contains no explicit Orphadata mechanism rows")
if orphadata_mondo_rows == 0:
    raise SystemExit(f"{evidence_tsv} contains no MONDO-mapped Orphadata rows")
print(
    f"evidence_rows={row_count}; g2p_rows={g2p_rows}; "
    f"g2p_lof_rows={g2p_lof_rows}; "
    f"g2p_penetrance_rows={g2p_penetrance_rows}; "
    f"orphadata_rows={orphadata_rows}; "
    f"orphadata_mondo_rows={orphadata_mondo_rows}"
)
