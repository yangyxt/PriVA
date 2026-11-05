#!/usr/bin/env python3
"""Compare `priva hpo-condition-mechanism` with build_hpo_condition_mechanism_cache.py.

  python hpo_condition_mechanism_compare.py check --binary PRIVA --repo CHECKOUT --work-root NEW_DIR [--jobs N]
  python hpo_condition_mechanism_compare.py real --binary PRIVA --repo CHECKOUT --work-root NEW_DIR \
      --hpo-assertions TSV --mechanism-evidence TSV --clingen-dosage TSV --gofcards-variants JSON_GZ \
      --hpo-release TEXT --mondo-release TEXT

Each case writes the same inputs into two folders and runs its steps in both:
the script in one and the native command in the other. A step is a command
line or an edit of the files left by the previous steps. Exit status, the
failure category and message, the printed counts and every file left in the
folder (with its mode) are compared. The folder's own path and the build
time are normalized; gzip files are compared by their decompressed text.
`real` builds once with each from existing inputs read in place, compares
the two caches and checks that each one passes the other's provenance check.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import copy
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

STAMP = re.compile(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(?:\.\d+)?\+00:00")
SCRIPT = "scripts/build_hpo_condition_mechanism_cache.py"

HPO_COLUMNS = ["gene_symbol", "disease_id", "hpo_id", "frequency", "evidence", "reference",
               "mondo_id", "mondo_name", "disease_scope", "priva_scope", "scope_evidence",
               "scope_reference", "scope_review_status"]
MECHANISM_COLUMNS = ["gene_symbol", "source", "source_record_id", "source_condition_id", "mondo_id",
                     "disease_scope", "priva_scope", "scope_review_status", "disease_label",
                     "inheritance", "penetrance_raw", "penetrance_hpo_ids", "normalized_penetrance",
                     "patho_mode_raw", "normalized_mechanisms", "mechanism_confidence",
                     "disease_confidence", "pmids", "evidence_url"]
CLINGEN_COLUMNS = ["#Gene Symbol", "Genomic Location", "Haploinsufficiency Score",
                   "Haploinsufficiency Description", "Haploinsufficiency Disease ID",
                   "Haploinsufficiency PMID1", "Haploinsufficiency PMID2"]


def hpo(gene, disease, term, mondo, name, scope="include", category="germline_disease",
        evidence="TAS", reference="PMID:1", frequency="-", scope_evidence="", scope_reference="",
        review="automatic"):
    return [gene, disease, term, frequency, evidence, reference, mondo, name, category, scope,
            scope_evidence, scope_reference, review]


HPO_ROWS = [
    hpo("GENE1", "OMIM:100", "HP:0000007", "MONDO:0000001", "disease one",
        scope_evidence="rule a;rule b;-", scope_reference="PMID:5; PMID:5", review="manually_confirmed"),
    hpo("GENE1", "OMIM:100", "HP:0003829", "MONDO:0000001", "disease one", reference="PMID:2"),
    hpo("GENE1", "OMIM:100", "HP:0003829", "MONDO:0000001", "disease one", reference="PMID:2"),
    hpo("GENE1", "OMIM:100", "HP:0003581", "MONDO:0000001", "disease one", frequency="HP:0040281"),
    hpo("GENE1", "ORPHA:200", "HP:0003584", "MONDO:0000001", "", scope="review", review="pending"),
    hpo("GENE1", "OMIM:101", "HP:0000006", "MONDO:0000002", "dominant one"),
    hpo("gene1", "OMIM:101", "HP:0034950", "MONDO:0000002", "dominant one", review="pending"),
    hpo("GENE1", "OMIM:102", "HP:0000118", "MONDO:0000010", "\"quoted\tlabel\""),
    hpo("GENE2", "OMIM:300", "HP:0000007", "MONDO:0000003", "recessive two"),
    hpo("GENE2", "OMIM:301", "HP:0000006", "MONDO:0000004", "excluded two", scope="exclude",
        category="somatic"),
    hpo("GENE3", "OMIM:400", "HP:0001419", "MONDO:0000005", "x recessive"),
    hpo("GENE3", "", "HP:0001417", "MONDO:0000011", "mondo keyed", scope=""),
    hpo("GENE4", "OMIM:500", "HP:0000118", "MONDO:0000006", "x dosage"),
    hpo("GENE4", "OMIM:501", "HP:0001423", "MONDO:0000007", "x dominant"),
    hpo("GENE5", "OMIM:600", "HP:0001450", "MONDO:0000008", "y linked"),
    hpo("GENE6", "OMIM:700", "HP:0000006", "MONDO:0000009", "only excluded", scope="exclude"),
    hpo("", "OMIM:800", "HP:0000006", "MONDO:0000012", "no gene"),
    hpo("GENE7", "", "HP:0000006", "", "no disease"),
    hpo("GENE8", "OMIM:900", "HP:0000006", "MONDO:0000013", "unscoped", scope=""),
    hpo("GENE8", "ORPHA:901", "HP:0000007", "MONDO:0000013", "merged review", scope="review"),
    hpo("GENE8", "OMIM:902", "HP:0003577", "MONDO:0000014", "null text", reference="NULL",
        evidence=" nan "),
]


def mechanism(gene, source, record, condition, mondo, mechanisms, inheritance="", penetrance="",
              hpo_ids="", raw="", pmids="1;2;2", scope="include", label="label"):
    return [gene, source, record, condition, mondo, "germline_disease", scope, "automatic", label,
            inheritance, raw, hpo_ids, penetrance, "raw mode", mechanisms, "definitive", "strong",
            pmids, "https://example.org/" + record]


MECHANISM_ROWS = [
    mechanism("GENE1", "G2P_DDG2P", "g2p-1", "OMIM:100", "", "LOF;GOF", "biallelic_autosomal",
              "complete", "HP:0034950;hp:0003829;xx", "complete penetrance"),
    mechanism("GENE1", "Orphadata", "orpha-1", "ORPHA:200", "", "DOMINANT_NEGATIVE"),
    mechanism("GENE1", "G2P_DDG2P", "g2p-2", "OMIM:999", "", "lof|unknown", "monoallelic_autosomal",
              "incomplete"),
    mechanism("GENE1", "G2P_DDG2P", "g2p-2", "OMIM:999", "", "lof|unknown", "monoallelic_autosomal",
              "incomplete"),
    mechanism("GENE1", "PanelApp", "panel-1", "OMIM:100", "", "LOF"),
    mechanism("GENE2", "G2P_DDG2P", "g2p-3", "OMIM:300", "MONDO:0000004", "GOF"),
    mechanism("GENE9", "Orphadata", "orpha-9", "ORPHA:1", "", "GOF", penetrance="Complete"),
    mechanism("", "G2P_DDG2P", "g2p-blank", "OMIM:100", "", "LOF"),
    mechanism("GENE8", "G2P_DDG2P", "g2p-8", "MONDO:0000013", "", "", penetrance="incomplete",
              label="from g2p"),
]

CLINGEN_ROWS = [
    ["GENE1", "chr7:1-2", "3", "disease one HI", "MONDO:0000001", "11", ""],
    ["GENE2", "chr1:1-2", "30", "recessive", "MONDO:0000003", "", ""],
    ["GENE4", "chrX:1-2", "3", "x dosage HI", "MONDO:0000006", "12", "13"],
    ["GENE4", "X:1-2", "3", "x dominant HI", "MONDO:0000007", "", ""],
    ["GENE5", "chrY:1-2", "3", "y HI", "MONDO:0000008", "", ""],
    ["GENE3", "chrX:1-2", "2", "x some", "MONDO:0000005", "", ""],
    ["GENEX", "chr2:1-2", "3", "absent", "MONDO:0000099", "", ""],
    ["GENE2", "chr1:1-2", "1", "absent disease", "MONDO:0000999", "", ""],
    ["GENE8", "chr23:1-2", "3", "unusual chromosome", "MONDO:0000013", "", ""],
    ["GENE1", "chr7:1-2", "0", "no evidence", "MONDO:0000001", "", ""],
]


def variant(pos, eligibility="eligible", clinvar=None, additional=None):
    record = {
        "record": {"eligibility": {"status": eligibility}, "source": {"start": str(pos)}},
        "assemblies": {
            "hg38": {"genomic": {"chrom": "1", "pos": pos + 10},
                     "transcripts": {"ENST1.2": {"by_hgvsc": {f"c.{pos}A>G": {"hgvsp": f"p.K{pos}R"}}}}},
            "hg19": {"genomic": {"chrom": "1", "pos": pos},
                     "transcripts": {"ENST1.1": {"by_hgvsc": {f"c.{pos}A>G": {"hgvsp": f"p.K{pos}R"},
                                                              f"c.{pos}del": None}},
                                     "ENST0.1": {"by_hgvsc": {f"c.{pos + 1}A>G": {}}}}},
        },
    }
    if clinvar:
        record["clinvar"] = clinvar
    if additional:
        record["clinvar_additional"] = additional
    return record


def clinvar_block(vcv, assertions, hgvs=("NM_1.1:c.1A>G", "", None)):
    return {"vcv_accession": vcv, "hgvs": list(hgvs), "condition_assertions": assertions}


def assertion(rcv, conditions=(), mappings=(), stars=2):
    return {"rcv_accession": rcv,
            "germline_classification": {"clinical_significance": "Pathogenic", "review_stars": stars},
            "conditions": list(conditions),
            "matched_scvs": [{"trait_mappings": list(mappings)}, "not a record"]}


def gofcards_cache():
    omim101 = {"database": "OMIM", "id": "101", "name": "dominant one"}
    medgen = {"database": "MedGen", "id": "C1", "name": "medgen only"}
    return {
        "metadata": {"source": "GoFCards", "mechanism": "GOF"},
        "genes": {
            "GENE1": {"variants": {
                "loc_1:10_grch37": variant(10, clinvar=clinvar_block("VCV1", [
                    assertion("RCV1", [omim101, medgen, "text"]),
                    assertion("RCV1", [omim101, medgen]),
                ]), additional=[clinvar_block("VCV1b", [assertion(
                    "RCV1b", [], [{"mapping_ref": "Orphanet", "mapping_value": "ORPHA:200"},
                                  {"mapping_ref": "OMIM", "mapping_value": "MONDO:1"},
                                  {"mapping_ref": "HP", "mapping_value": "HP:1"}], stars=1.5)])]),
                "loc_1:20_grch37": variant(20, clinvar=clinvar_block("VCV2", [assertion(
                    "RCV2", [omim101, {"database": "omim", "id": "OMIM:100"}])])),
                "loc_1:30_grch37": variant(30),
                "loc_1:40_grch37": variant(40, clinvar=clinvar_block("VCV4", [assertion(
                    "RCV4", [{"database": "OMIM", "id": 999}])])),
                "loc_1:50_grch37": variant(50, eligibility="quarantined_gene_discordance"),
            }},
            "gene7": {"variants": {
                "loc_7:10_grch37": variant(70, clinvar=clinvar_block("VCV7", [assertion(
                    "RCV7", [{"database": "MONDO", "id": "MONDO:0000001"}])])),
            }},
            "GENE2": {"variants": {
                "loc_2:10_grch37": variant(80, clinvar=clinvar_block("VCV8", [assertion(
                    "RCV8", [{"database": "ORPHANET", "id": "ORPHANET:9"}])])),
            }},
        },
    }


def tsv(rows, columns):
    return "\t".join(columns) + "\n" + "".join("\t".join(row) + "\n" for row in rows)


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.name.endswith(".gz"):
        with path.open("wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as out:
            out.write(text.encode())
    else:
        path.write_bytes(text.encode())


def standard_inputs(folder, hpo_rows=None, mechanism_rows=None, clingen_rows=None, cache=None,
                    hpo_name="hpo.tsv", hpo_text=None, mechanism_text=None, clingen_text=None,
                    cache_text=None):
    write(folder / "in" / hpo_name, hpo_text if hpo_text is not None
          else tsv(HPO_ROWS if hpo_rows is None else hpo_rows, HPO_COLUMNS))
    write(folder / "in/mechanism.tsv", mechanism_text if mechanism_text is not None
          else tsv(MECHANISM_ROWS if mechanism_rows is None else mechanism_rows, MECHANISM_COLUMNS))
    write(folder / "in/clingen.tsv", clingen_text if clingen_text is not None
          else tsv(CLINGEN_ROWS if clingen_rows is None else clingen_rows, CLINGEN_COLUMNS))
    write(folder / "in/gofcards.json.gz", cache_text if cache_text is not None
          else json.dumps(gofcards_cache() if cache is None else cache))


def build_args(output="out/cache.json.gz", hpo_name="hpo.tsv", extra=()):
    return ["--hpo-assertions", f"in/{hpo_name}", "--mechanism-evidence", "in/mechanism.tsv",
            "--clingen-dosage", "in/clingen.tsv", "--gofcards-variants", "in/gofcards.json.gz",
            "--hpo-release", "2026-10-09", "--mondo-release", "v2026-10-01", "--output", output, *extra]


def provenance_args(cache="out/cache.json.gz", hpo_release="2026-10-09", mondo_release="v2026-10-01"):
    return ["--validate-only", cache, "--hpo-assertions", "in/hpo.tsv", "--mechanism-evidence",
            "in/mechanism.tsv", "--clingen-dosage", "in/clingen.tsv", "--gofcards-variants",
            "in/gofcards.json.gz", "--hpo-release", hpo_release, "--mondo-release", mondo_release]


def read_cache(path):
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.name.endswith(".gz") else data)


def edit_cache(change, name="out/cache.json.gz"):
    def apply(folder):
        path = folder / name
        payload = read_cache(path)
        change(payload)
        write(path, json.dumps(payload, sort_keys=True))
    return apply


def append(name, text):
    def apply(folder):
        with (folder / name).open("a") as handle:
            handle.write(text)
    return apply


def condition(payload, gene="GENE1", key="OMIM:100"):
    return payload["genes"][gene]["conditions"][key]


def set_counts(payload, name, value):
    payload["_meta"]["counts"][name] = value


CASES = {}


def case(name, setup, *steps):
    CASES[name] = (setup, list(steps))


def run_cases():
    base = standard_inputs
    case("build_gzip", base, build_args())
    case("build_plain", base, build_args("out/cache.json"))
    case("build_pretty_nested", base, build_args("out/a/b/cache.json", extra=["--pretty"]))
    case("build_gzip_inputs", lambda f: base(f, hpo_name="hpo.tsv.gz"), build_args(hpo_name="hpo.tsv.gz"))
    case("build_default_clingen", base,
         [a for a in build_args() if a not in ("--clingen-dosage", "in/clingen.tsv")])
    case("build_empty_clingen_argument", base,
         [("" if a == "in/clingen.tsv" else a) for a in build_args()])
    case("build_no_releases", base, build_args()[:8] + ["--output", "out/cache.json"])
    case("build_replaces_existing", lambda f: (base(f), write(f / "out/cache.json", "old\n")),
         build_args("out/cache.json"))
    case("validate_plain", base, build_args(), ["--validate-only", "out/cache.json.gz"])
    case("validate_provenance", base, build_args(), provenance_args())
    case("validate_provenance_pretty_plain", base, build_args("out/c.json", extra=["--pretty"]),
         provenance_args("out/c.json"))
    for name, target in (("hpo", "in/hpo.tsv"), ("mechanism", "in/mechanism.tsv"),
                         ("clingen", "in/clingen.tsv")):
        case(f"provenance_changed_{name}", base, build_args(), append(target, "\n"), provenance_args())
    case("provenance_changed_gofcards", base, build_args(),
         lambda f: write(f / "in/gofcards.json.gz", json.dumps(gofcards_cache(), indent=1)),
         provenance_args())
    case("provenance_releases", base, build_args(), provenance_args(hpo_release="other", mondo_release=" "))
    case("provenance_missing_source", base, build_args(),
         edit_cache(lambda p: p["_meta"]["sources"].pop("clingen_dosage")),
         edit_cache(lambda p: p["_meta"]["sources"]["hpo_assertions"].update(size_bytes=1.0)),
         provenance_args())
    case("provenance_null_releases", base, build_args(),
         edit_cache(lambda p: p["_meta"].update(releases={"HPO": None})), provenance_args())
    case("provenance_missing_input", base, build_args(),
         lambda f: (f / "in/clingen.tsv").unlink(), provenance_args())
    case("provenance_partial", base, build_args(),
         ["--validate-only", "out/cache.json.gz", "--hpo-release", ""])
    case("provenance_partial_input", base, build_args(),
         ["--validate-only", "out/cache.json.gz", "--hpo-assertions", "in/hpo.tsv"])

    edits = {
        "schema_version": lambda p: p["_meta"].update(schema_version="1.1"),
        "schema_version_number": lambda p: p["_meta"].update(schema_version=1.2),
        "no_genes": lambda p: p.update(genes={}),
        "genes_list": lambda p: p.update(genes=[1]),
        "gene_not_object": lambda p: p["genes"].update(GENE1=[]),
        "empty_symbol": lambda p: p["genes"].update({"": p["genes"]["GENE2"]}),
        "missing_sections": lambda p: p["genes"]["GENE1"].pop("unmapped_evidence"),
        "missing_scope_gate": lambda p: p["genes"]["GENE1"].pop("scope_gate"),
        "bool_excluded_count": lambda p: p["genes"]["GENE1"]["scope_gate"].update(excluded_condition_count=True),
        "negative_excluded_count": lambda p: p["genes"]["GENE1"]["scope_gate"].update(excluded_condition_count=-1),
        "survived_exclusion": lambda p: condition(p)["priva_scope"].update(decision=" Exclude "),
        "condition_list": lambda p: p["genes"]["GENE1"]["conditions"].update({"OMIM:100": []}),
        "key_absent": lambda p: condition(p).update(identifiers={"OMIM": ["OMIM:1"]}),
        "identifier_text": lambda p: condition(p).update(identifiers={"X": "OMIM:100"}),
        "identifier_number": lambda p: condition(p).update(identifiers={"X": 5}),
        "identifier_unhashable": lambda p: condition(p).update(identifiers={"X": [["OMIM:100"]]}),
        "identifiers_list": lambda p: condition(p).update(identifiers=["OMIM:100"]),
        "no_penetrance": lambda p: condition(p).pop("penetrance"),
        "status_unknown": lambda p: condition(p)["penetrance"].update(statuses=["partial"]),
        "status_unhashable": lambda p: condition(p)["penetrance"].update(statuses=[["complete"]]),
        "assertions_text": lambda p: condition(p)["penetrance"].update(assertions="x"),
        "assertion_extra_key": lambda p: condition(p)["penetrance"]["assertions"][0].update(extra=1),
        "assertion_value": lambda p: condition(p)["penetrance"]["assertions"][0].update(normalized_value="high"),
        "assertion_absent": lambda p: condition(p)["penetrance"].update(statuses=[]),
        "unsupported_mechanism": lambda p: condition(p)["pathogenic_mechanisms"].update(
            HAPLO={"variants": {}}),
        "mechanisms_null": lambda p: condition(p).update(pathogenic_mechanisms=None),
        "variant_link": lambda p: condition(p, key="OMIM:101")["pathogenic_mechanisms"]["GOF"]["variants"][
            "GOFCARDS:loc_1:10_grch37"]["condition_link"].update(condition_key="OMIM:100"),
        "unmapped_exact": lambda p: p["genes"]["GENE1"]["unmapped_evidence"]["variants"][
            "GOFCARDS:loc_1:30_grch37"]["condition_link"].update(status="exact"),
        "unmapped_penetrance": lambda p: p["genes"]["GENE1"]["unmapped_evidence"]["penetrance"].append({"a": 1}),
        "unmapped_penetrance_text": lambda p: p["genes"]["GENE1"]["unmapped_evidence"].update(penetrance="ab"),
        "counts": lambda p: set_counts(p, "conditions", 1),
        "counts_float": lambda p: set_counts(p, "genes", float(p["_meta"]["counts"]["genes"])),
        "no_inheritance_key": lambda p: condition(p).pop("inheritance"),
        "inheritance_list": lambda p: condition(p).update(inheritance=[]),
        "inheritance_text": lambda p: condition(p).update(inheritance="x"),
        "inheritance_null": lambda p: condition(p).update(inheritance=None),
        "no_mechanisms_key": lambda p: condition(p).pop("pathogenic_mechanisms"),
        "block_no_variants": lambda p: condition(p)["pathogenic_mechanisms"]["LOF"].pop("variants"),
        "unmapped_mechanisms_number": lambda p: p["genes"]["GENE1"]["unmapped_evidence"].update(mechanisms=3),
        "unmapped_no_variants": lambda p: p["genes"]["GENE1"]["unmapped_evidence"].pop("variants"),
        "scope_stats_missing": lambda p: p["_meta"]["build_statistics"].pop("scope_gate"),
        "build_statistics_list": lambda p: p["_meta"].update(build_statistics=[1]),
        "scope_stats_count": lambda p: p["_meta"]["build_statistics"]["scope_gate"].update(conditions_excluded="3"),
        "payload_list": lambda p: (p.clear(), p.update(x=1)),
    }
    for name, change in edits.items():
        case(f"invalid_{name}", base, build_args(), edit_cache(change), ["--validate-only", "out/cache.json.gz"])
    case("invalid_payload_array", base, lambda f: write(f / "out/cache.json", "[1]\n"),
         ["--validate-only", "out/cache.json"])
    case("invalid_json_text", base, lambda f: write(f / "out/cache.json", "{\"_meta\": \n"),
         ["--validate-only", "out/cache.json"])
    case("invalid_not_gzip", base, lambda f: (f / "out").mkdir() or (f / "out/cache.json.gz").write_text("{}"),
         ["--validate-only", "out/cache.json.gz"])
    case("invalid_missing_cache", base, ["--validate-only", "out/none.json"])

    case("missing_all_arguments", base, [])
    case("missing_output", base, build_args()[:-2])
    case("empty_output", base, build_args(""))
    case("missing_hpo_file", base, build_args(hpo_name="absent.tsv"))
    case("missing_hpo_columns", lambda f: base(f, hpo_text="gene_symbol\tdisease_id\nA\tB\n"), build_args())
    case("empty_hpo_file", lambda f: base(f, hpo_text=""), build_args())
    case("blank_header", lambda f: base(f, hpo_text="\n" + tsv(HPO_ROWS, HPO_COLUMNS)), build_args())
    case("hpo_header_only", lambda f: base(f, hpo_text=tsv([], HPO_COLUMNS), mechanism_rows=[],
                                         cache={"metadata": {"source": "GoFCards", "mechanism": "GOF"},
                                                "genes": {}}), build_args())
    case("mechanism_only_gene", lambda f: base(f, hpo_text=tsv([], HPO_COLUMNS)), build_args())
    case("hpo_crlf_short_rows_blank_lines", lambda f: base(f, hpo_text=(
        "\t".join(HPO_COLUMNS) + "\r\n" + "\r\n".join("\t".join(r) for r in HPO_ROWS[:6])
        + "\r\n\r\nGENE1\tOMIM:100\r\n\n" + "\t".join(HPO_ROWS[8]) + "\textra\tmore\r"
        + "\t".join(HPO_ROWS[9]))), build_args())
    case("hpo_multiline_quote", lambda f: base(f, hpo_text=tsv(HPO_ROWS[:3], HPO_COLUMNS) + "\t".join(
        hpo("GENE1", "OMIM:100", "HP:0000007", "MONDO:0000001", '"line\none""x"')) + "\n"), build_args())
    case("hpo_unterminated_quote", lambda f: base(f, hpo_text=tsv(HPO_ROWS[:3], HPO_COLUMNS) + "GENE1\t\"OMIM:100\n"),
         build_args())
    case("hpo_quote_inside_field", lambda f: base(f, hpo_text=tsv(HPO_ROWS[:3], HPO_COLUMNS) + "\t".join(
        hpo("GENE1", "OMIM:100", "HP:0000007", "MONDO:0000001", 'a"b')) + "\n" + "\t".join(
        hpo("GENE1", "OMIM:100", "HP:0000007", "MONDO:0000001", '"a"b')) + "\n"), build_args())
    case("hpo_field_limit", lambda f: base(f, hpo_text=tsv(HPO_ROWS[:3], HPO_COLUMNS) + "\t".join(
        hpo("GENE1", "OMIM:100", "HP:0000007", "MONDO:0000001", "x" * 131073)) + "\n"), build_args())
    case("hpo_field_at_limit", lambda f: base(f, hpo_text=tsv(HPO_ROWS[:3], HPO_COLUMNS) + "\t".join(
        hpo("GENE1", "OMIM:100", "HP:0000007", "MONDO:0000001", "x" * 131072)) + "\n"), build_args())
    case("hpo_invalid_utf8", lambda f: (base(f), (f / "in/hpo.tsv").write_bytes(
        tsv(HPO_ROWS[:3], HPO_COLUMNS).encode() + b"GENE1\t\xff\n")), build_args())
    case("hpo_duplicate_header", lambda f: base(f, hpo_text=tsv(
        [r + [r[0]] for r in HPO_ROWS[:4]] + [HPO_ROWS[5][:2]], HPO_COLUMNS + ["gene_symbol"])), build_args())
    case("conflicting_scope", lambda f: base(f, hpo_rows=HPO_ROWS + [
        hpo("GENE1", "OMIM:100", "HP:0000006", "MONDO:0000001", "x", scope="exclude")]), build_args())
    case("conflicting_category", lambda f: base(f, hpo_rows=HPO_ROWS + [
        hpo("GENE1", "OMIM:100", "HP:0000006", "MONDO:0000001", "x", category="other")]), build_args())
    case("unknown_scope_in_merge", lambda f: base(f, hpo_rows=HPO_ROWS + [
        hpo("GENE1", "ORPHA:201", "HP:0000006", "MONDO:0000002", "x", scope="maybe")]), build_args())
    case("unknown_scope_alone", lambda f: base(f, hpo_rows=HPO_ROWS + [
        hpo("GENE1", "ORPHA:202", "HP:0000006", "MONDO:0000077", "x", scope="maybe")]), build_args())
    case("disease_without_colon", lambda f: base(f, hpo_rows=HPO_ROWS + [
        hpo("GENE1", "12345", "HP:0000006", "", "x")]), build_args())
    case("unsupported_penetrance", lambda f: base(f, mechanism_rows=MECHANISM_ROWS + [
        mechanism("GENE1", "G2P_DDG2P", "", "OMIM:100", "", "LOF", penetrance="Variable")]), build_args())
    case("missing_mechanism_columns", lambda f: base(f, mechanism_text="gene_symbol\tsource\n"), build_args())
    case("missing_clingen_columns", lambda f: base(f, clingen_text="#Gene Symbol\n"), build_args())
    case("gofcards_wrong_source", lambda f: base(f, cache={"metadata": {"source": "x"}, "genes": {}}),
         build_args())
    case("gofcards_no_genes", lambda f: base(f, cache={"metadata": {}}), build_args())
    case("gofcards_list", lambda f: base(f, cache_text="[\"genes\"]"), build_args())
    case("gofcards_not_json", lambda f: base(f, cache_text="{"), build_args())
    case("gofcards_missing", lambda f: (base(f), (f / "in/gofcards.json.gz").unlink()), build_args())
    case("gofcards_not_gzip", lambda f: (base(f), (f / "in/gofcards.json.gz").write_text("{}")), build_args())
    bad = copy.deepcopy(gofcards_cache())
    bad["genes"]["GENE1"]["variants"]["loc_1:10_grch37"]["clinvar"]["condition_assertions"] = ["text"]
    case("gofcards_assertion_text", lambda f: base(f, cache=bad), build_args())
    bad2 = copy.deepcopy(gofcards_cache())
    bad2["genes"]["GENE1"]["variants"]["loc_1:20_grch37"]["clinvar_additional"] = 5
    case("gofcards_additional_number", lambda f: base(f, cache=bad2), build_args())
    bad3 = copy.deepcopy(gofcards_cache())
    bad3["genes"]["GENE1"]["variants"]["loc_1:20_grch37"]["clinvar_additional"] = {"vcv": 1}
    case("gofcards_additional_object", lambda f: base(f, cache=bad3), build_args())
    bad4 = copy.deepcopy(gofcards_cache())
    bad4["genes"]["GENE1"]["variants"][" "] = variant(90)
    case("gofcards_blank_variant_id", lambda f: base(f, cache=bad4), build_args())
    bad5 = copy.deepcopy(gofcards_cache())
    bad5["genes"]["GENE1"]["variants"]["loc_1:30_grch37"]["assemblies"]["hg19"] = []
    case("gofcards_assembly_list", lambda f: base(f, cache=bad5), build_args())
    bad6 = copy.deepcopy(gofcards_cache())
    bad6["genes"]["GENE1"]["variants"]["loc_1:30_grch37"]["record"] = "eligible"
    case("gofcards_record_text", lambda f: base(f, cache=bad6), build_args())
    case("unknown_option", base, build_args(extra=["--extra", "1"]))
    case("flag_with_value", base, build_args(extra=["--pretty=1"]))
    case("option_without_value", base, ["--output"])
    case("prefix_options", base, ["--hpo-a", "in/hpo.tsv", "--mech", "in/mechanism.tsv", "--clin",
                                  "in/clingen.tsv", "--gof", "in/gofcards.json.gz", "--out", "out/x.json"])
    case("help", base, ["--help"])


def outcome(command, folder):
    result = subprocess.run([str(part) for part in command], cwd=folder, capture_output=True, text=True)
    lines = result.stderr.splitlines()
    value = {"exit": result.returncode}
    if result.returncode == 1 and lines:
        last = lines[-1]
        match = re.match(r"priva: ([A-Za-z]+): (.*)", last) or re.match(
            r"(?:[\w.]+\.)?([A-Za-z]+(?:Error|Exception|Exit)?): (.*)", last)
        if match:
            value["category"], value["message"] = match[1], match[2]
        else:
            value["category"], value["message"] = "SystemExit", last
        value["message"] = STAMP.sub("<TIME>", value["message"].replace(str(folder), "<ARM>"))
    if result.returncode == 0 and result.stdout.strip() and not result.stdout.startswith("usage"):
        value["stdout"] = result.stdout
    return value


def normalize(value, folder):
    if isinstance(value, dict):
        return {key: normalize(item, folder) for key, item in value.items()}
    if isinstance(value, list):
        return [normalize(item, folder) for item in value]
    if isinstance(value, str):
        return STAMP.sub("<TIME>", value.replace(str(folder), "<ARM>"))
    return value


def snapshot(folder):
    state = {}
    for path in sorted(folder.rglob("*")):
        if path.is_dir():
            continue
        relative = str(path.relative_to(folder))
        data = path.read_bytes()
        if path.name.endswith(".gz"):
            try:
                data = gzip.decompress(data)
            except OSError:
                pass
        mode = oct(path.stat().st_mode & 0o777)
        try:
            body = normalize(json.loads(data), folder)
            # The compact and indented forms are compared byte for byte too.
            text = data.decode().replace(str(folder), "<ARM>")
            state[relative] = [mode, body, STAMP.sub("<TIME>", text)]
        except (json.JSONDecodeError, UnicodeDecodeError):
            state[relative] = [mode, hashlib.sha256(data).hexdigest()]
    return state


def run_case(name, binary, repo, work_root):
    setup, steps = CASES[name]
    arms = {}
    for arm in ("python", "native"):
        folder = work_root / name / arm
        folder.mkdir(parents=True)
        setup(folder)
        results = []
        for step in steps:
            if callable(step):
                step(folder)
                continue
            if arm == "python":
                command = [sys.executable, repo / SCRIPT, *step]
            else:
                command = [binary, "hpo-condition-mechanism", "--source-root", repo, *step]
            results.append(outcome(command, folder))
        arms[arm] = {"steps": results, "files": snapshot(folder)}
    same = arms["python"] == arms["native"]
    return name, same, len([s for s in steps if not callable(s)]), arms


def check(args):
    run_cases()
    work_root = Path(args.work_root)
    work_root.mkdir(parents=True)
    repo = Path(args.repo).resolve()
    failed, steps = [], 0
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(run_case, name, Path(args.binary), repo, work_root) for name in CASES]
        for future in futures:
            name, same, count, arms = future.result()
            steps += count
            if not same:
                failed.append(name)
                print(f"DIFFERENT {name}", file=sys.stderr)
                for arm in ("python", "native"):
                    print(f"  {arm}: {json.dumps(arms[arm]['steps'])[:2000]}", file=sys.stderr)
                py_files, na_files = arms["python"]["files"], arms["native"]["files"]
                for key in sorted(set(py_files) | set(na_files)):
                    if py_files.get(key) != na_files.get(key):
                        print(f"  file {key} differs", file=sys.stderr)
    print(json.dumps({"cases": len(CASES), "steps": steps, "failed": failed}))
    return 1 if failed else 0


def real(args):
    work_root = Path(args.work_root)
    work_root.mkdir(parents=True)
    repo = Path(args.repo).resolve()
    inputs = ["--hpo-assertions", args.hpo_assertions, "--mechanism-evidence", args.mechanism_evidence,
              "--clingen-dosage", args.clingen_dosage, "--gofcards-variants", args.gofcards_variants,
              "--hpo-release", args.hpo_release, "--mondo-release", args.mondo_release]
    commands = {"python": [sys.executable, repo / SCRIPT],
                "native": [args.binary, "hpo-condition-mechanism", "--source-root", repo]}
    report, caches = {}, {}
    for arm, command in commands.items():
        output = work_root / arm / "hpo_condition_mechanism_cache.json.gz"
        start = time.monotonic()
        result = subprocess.run([str(p) for p in [*command, *inputs, "--output", output]],
                                capture_output=True, text=True)
        report[f"{arm}_seconds"] = round(time.monotonic() - start, 1)
        report[f"{arm}_exit"] = result.returncode
        report[f"{arm}_counts"] = result.stdout.strip()
        if result.returncode:
            report[f"{arm}_stderr"] = result.stderr[-2000:]
            continue
        caches[arm] = output
    if len(caches) == 2:
        bodies = {arm: gzip.decompress(path.read_bytes()).decode() for arm, path in caches.items()}
        report["same"] = STAMP.sub("<TIME>", bodies["python"]) == STAMP.sub("<TIME>", bodies["native"])
        report["genes"] = len(json.loads(bodies["native"])["genes"])
        for checker, cache in (("python", "native"), ("native", "python")):
            result = subprocess.run([str(p) for p in [*commands[checker], "--validate-only", caches[cache],
                                                       *inputs]], capture_output=True, text=True)
            report[f"{checker}_accepts_{cache}"] = result.returncode == 0
    else:
        report["same"] = False
    print(json.dumps(report))
    return 0 if report["same"] and all(v for k, v in report.items() if "_accepts_" in k) else 1


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)
    for mode in ("check", "real"):
        part = sub.add_parser(mode)
        part.add_argument("--binary", required=True)
        part.add_argument("--repo", required=True)
        part.add_argument("--work-root", required=True)
        if mode == "check":
            part.add_argument("--jobs", type=int, default=4)
        else:
            for name in ("hpo-assertions", "mechanism-evidence", "clingen-dosage", "gofcards-variants",
                         "hpo-release", "mondo-release"):
                part.add_argument("--" + name, required=True)
    args = parser.parse_args()
    raise SystemExit(check(args) if args.mode == "check" else real(args))


if __name__ == "__main__":
    main()
