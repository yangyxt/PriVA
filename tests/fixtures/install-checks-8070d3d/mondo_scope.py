import csv
import gzip
import sys

registry_path, assertions_path = sys.argv[1:]

def inspect(path, required):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise SystemExit(f"{path} missing required columns: {', '.join(missing)}")
        count = sum(1 for _ in reader)
    if count == 0:
        raise SystemExit(f"{path} contains no data rows")
    return count

registry_count = inspect(registry_path, {
    "disease_id", "mondo_id", "disease_scope", "priva_scope",
    "scope_evidence", "scope_review_status",
})
assertion_count = inspect(assertions_path, {
    "gene_symbol", "disease_id", "hpo_id", "frequency", "evidence",
    "reference", "mondo_id", "disease_scope", "priva_scope",
    "scope_evidence", "scope_review_status",
})
print(f"disease_scope_rows={registry_count}; scoped_hpo_assertions={assertion_count}")
