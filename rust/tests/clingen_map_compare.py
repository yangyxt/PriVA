#!/usr/bin/env python3
"""Differential fixtures for the opt-in native ClinGen producer.

  python clingen_map_compare.py edges --binary PRIVA \
      --python-producer scripts/prepare_clingen_map.py --work-root NEW_DIRECTORY

The reference producer is executed unchanged. Each successful Python pickle is
converted to schema 10 from those same fixture inputs; native output is compared
against that pickle with CPython's independent full-value comparator. Text/NaN
maps retain schema 10 v1; surviving numeric/Boolean codes require typed v2.
Logs and fixtures are written under the supplied directory, never beside source.
"""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import cache_values_compare as oracle


HEADER = ("##fileformat=VCFv4.2\n##contig=<ID=1,length=10000>\n"
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
TABLE = "ClinVar Variation Id\tApplied Evidence Codes (Met)\n"
RECORDS = ("1\t10\t1\tA\tT,C\t.\t.\t.\n"
           "1\t11\t2\tG\tC\t.\t.\t.\n"
           "1\t12\t0\tC\tA\t.\t.\t.\n"
           "1\t13\t-9223372036854775808\tT\tG\t.\t.\t.\n"
           "1\t14\t1000\tA\tG\t.\t.\t.\n")


def cases() -> dict[str, dict]:
    def case(rows="1\tPS3\n2\tPM2\n", records=RECORDS, **options):
        return {"table": TABLE + rows, "vcf": HEADER + records, **options}

    result = {
        "ordinary_first_alt": case(),
        "duplicate_vcf_ids": case(records="1\t10\t1; 2 ;.; ;1\tA\tT,C\n1\t11\t1\tG\tC\n"),
        "duplicate_variant_rows_last_nan": case(rows="1\tPS3\n2\tPM2\n1\tBA1\n3\tNA\n4\tPVS1\n",
              records="1\t10\t1;2;3\tA\tT,C\n1\t11\t1\tG\tC\n"),
        "missing_codes": case(rows="1\tNA\n2\t\n"),
        "float_id_truncation": case(rows="1.9\tPS3\n\tPM2\n2.0\tNA\n"),
        "boolean_ids": case(rows="True\tPS3\nFalse\tPM2\n"),
        "boolean_ids_with_na": case(rows="True\tPS3\n\tPM2\nFalse\tPVS1\n"),
        "na_markers": case(rows="NA\tPS3\nNULL\tPM2\nNone\tBA1\n1\tNaN\n"),
        "finite_float_overflow": case(rows="1e30\tPS3\n1\tPM2\n"),
        "uint64_wrap": case(rows="9223372036854775808\tPS3\n1\tPM2\n"),
        "uint64_with_na_failure": case(rows="9223372036854775808\tPS3\n\tPM2\n1\tBA1\n", fails=True),
        "int64_min_missing_sentinel": case(rows="-9223372036854775808\tPS3\n\tPM2\n1\tBA1\n"),
        "signed_whitespace": case(rows=" +001 \tPS3\n-2\tPM2\n"),
        "object_integer_underscores": case(rows="1_000\tPS3\n1\tPM2\n"),
        "unicode_integer": case(rows="１\tPS3\n２\tPM2\n"),
        "unicode_integer_whitespace": case(rows="\u00a0+１\u2003\tPS3\n２\tPM2\n"),
        "object_float_text_failure": case(rows="1.0\tPS3\nx\tPM2\n", fails=True),
        "infinite_id_failure": case(rows="inf\tPS3\n1\tPM2\n", fails=True),
        "huge_object_id_failure": case(rows="18446744073709551616\tPS3\n1\tPM2\n", fails=True),
        "numeric_codes_refused": case(rows="1\t1\n2\t2\n", version=2),
        "boolean_codes_refused": case(rows="1\tTrue\n2\tFalse\n", version=2),
        "float_codes_refused": case(rows="1\t1.5\n2\t2.5\n", version=2),
        "uint64_codes": case(rows="1\t18446744073709551615\n2\t9223372036854775808\n", version=2),
        "signed_integer_codes": case(rows="1\t-9223372036854775808\n2\t-2\n", version=2),
        "signed_zero_codes": case(rows="1\t-0.0\n2\t0.0\n", version=2),
        "infinite_codes": case(rows="1\tinf\n2\t-inf\n", version=2),
        "boolean_codes_missing": case(rows="1\tFalse\n2\tNA\n", version=2),
        "float_codes_missing": case(rows="1\t2.5\n2\tNA\n", version=2),
        "numeric_codes_last_missing": case(rows="1\t2\n1\tNA\n", version=1),
        "unmapped_numeric_codes": case(rows="9\t1\n8\t2\n"),
        "mapped_nan_unmapped_numeric": case(rows="1\tNA\n9\t2\n"),
        "last_nan_overwrites_numeric": case(rows="1\t2\n1\tNA\n"),
        "object_codes_preserve_text": case(rows="1\t1.5\n2\tPS3\n"),
        "uint_codes_preserve_na": case(rows="9\t9223372036854775808\n8\t1\n1\tNA\n"),
        "uint_negative_text_codes_preserve_na": case(rows="9\t9223372036854775808\n8\t-x\n1\tNA\n"),
        "uint_negative_overflow_codes_preserve_na": case(rows="9\t9223372036854775808\n8\t-9223372036854775809\n1\tNA\n"),
        "negative_overflow_codes_drop_na": case(rows="9\t-9223372036854775809\n8\t1\n1\tNA\n"),
        "uint_text_codes_drop_na": case(rows="9\t9223372036854775808\n8\tx\n1\tNA\n"),
        "all_missing_ids": case(rows="NA\tPS3\n\tPM2\n"),
        "empty_table": case(rows=""),
        "empty_vcf": case(records=""),
        "quoted_codes": case(rows='1\t"PS3\tPM2\nPVS1"\n2\t"A""B"\n'),
        "implicit_index": {"table": TABLE + "row\t1\tPS3\nrow2\t2\tPM2\n", "vcf": HEADER + RECORDS},
        "mangled_headers": {"table": "\ufeff" + TABLE.rstrip("\n") + "\tClinVar Variation Id\r\n1\tPS3\tbad\r\n", "vcf": HEADER + RECORDS},
        "extra_late_table_field": case(rows="1\tPS3\n2\tPM2\textra\n", fails=True),
        "missing_table_column": {"table": "ClinVar Variation Id\n1\n", "vcf": HEADER + RECORDS, "fails": True},
        "alt_dot": case(records="1\t10\t1\tA\t.\n", fails=True),
        "unused_alt_dot": case(records="1\t10\t.\tA\t.\n", fails=True),
        "empty_table_alt_dot": case(rows="", records="1\t10\t.\tA\t.\n", fails=True),
        "alt_dot_first": case(records="1\t10\t1\tA\t.,T\n"),
        "empty_alleles_pos_zero": case(records="1\t0\t1\t\t\n"),
        "symbolic_first_alt": case(records="1\t10\t1\tA\t<DEL>,C\n"),
        "five_fields": case(records="1\t10\t1\tA\tT\n"),
        "seven_fields": case(records="1\t10\t1\tA\tT\tbad\t.\n"),
        "four_fields_failure": case(records="1\t10\t1\tA\n", fails=True),
        "negative_pos_failure": case(records="1\t-1\t1\tA\tT\n", fails=True),
        "noninteger_pos_failure": case(records="1\tx\t1\tA\tT\n", fails=True),
        "float_pos_failure": case(records="1\t1.5\t1\tA\tT\n", fails=True),
        "blank_record_failure": case(records="\n" + RECORDS, fails=True),
        "record_comment_failure": case(records="#comment\n" + RECORDS, fails=True),
        "bad_header_failure": {"table": TABLE + "1\tPS3\n", "vcf": HEADER.replace("#CHROM", "#CHR") + RECORDS, "fails": True},
        "missing_fileformat_failure": {"table": TABLE + "1\tPS3\n", "vcf": HEADER.split("\n", 1)[1] + RECORDS, "fails": True},
    }
    result["gzip_table"] = {**case(), "gzip_table": True}
    result["bgzf_vcf"] = {**case(), "bgzf": True}
    result["ordinary_gzip_vcf_failure"] = {**case(), "gzip_vcf": True, "fails": True}
    result["crlf"] = {"table": (TABLE + "1\tPS3\n").replace("\n", "\r\n"),
                      "vcf": (HEADER + RECORDS).replace("\n", "\r\n")}
    result["empty_position"] = case(records="1\t\t1\tA\tT\n")
    result["positive_signed_position"] = case(records="1\t+10\t1\tA\tT\n")
    for kind, values in {"GT": [("", False), (".", True), ("0//1", False), ("0|", False), ("/1", False), ("0a", False), ("+1", True), ("-1", False), ("1/99", True)],
                         "DP": [("", False), (".", True), ("-2", True), ("+2", True), ("1x", False), ("1.0", False), (" 1", False), ("1 ", False), ("2147483648", True), ("999999999999999999999999999", True)],
                         "AF": [("", False), (".", True), ("nan", True), ("NaN", True), ("inf", True), ("Infinity", True), ("1e999", True), (" 1", True), ("1 ", False), ("x", False)]}.items():
        value_type = {"GT": "String", "DP": "Integer", "AF": "Float"}[kind]
        typed = HEADER.replace("#CHROM", f'##FORMAT=<ID={kind},Number=1,Type={value_type},Description="Value">\n#CHROM').rstrip("\n") + "\tFORMAT\tS1\n"
        for i, (value, succeeds) in enumerate(values):
            result[f"sample_{kind}_{i}"] = {"table": TABLE + "1\tPS3\n", "vcf": typed + f"1\t10\t1\tA\tT\t.\t.\t.\t{kind}\t{value}\n", "fails": not succeeds}
    result["format_without_sample"] = {"table": TABLE + "1\tPS3\n", "vcf": typed + "1\t10\t1\tA\tT\t.\t.\t.\tAF\n", "fails": True}
    result["header_samples_short_record"] = {"table": TABLE + "1\tPS3\n", "vcf": typed + "1\t10\t1\tA\tT\n"}
    return result


def identity(path: Path) -> dict:
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": oracle.sha256(path)}


def invoke(command: list[str], root: Path, name: str) -> int:
    with (root / f"{name}.stdout.log").open("x") as out, (root / f"{name}.stderr.log").open("x") as err:
        return subprocess.run(command, stdout=out, stderr=err, check=False).returncode


def edges(args: argparse.Namespace) -> None:
    import gzip
    args.work_root.mkdir(parents=True, exist_ok=False)
    sources = {"binary": identity(args.binary), "producer": identity(args.python_producer),
               "comparator": identity(Path(oracle.__file__).resolve()), "fixtures": identity(Path(__file__).resolve())}
    results = {}
    for name, fixture in cases().items():
        root = args.work_root / name
        root.mkdir()
        table, vcf = root / "curated.tsv", root / "clinvar.vcf"
        if fixture.get("gzip_table"):
            table = Path(str(table) + ".gz")
            with gzip.open(table, "wt", newline="") as stream:
                stream.write(fixture["table"])
        else:
            table.write_text(fixture["table"])
        if fixture.get("gzip_vcf"):
            vcf = Path(str(vcf) + ".gz")
            with gzip.open(vcf, "wt", newline="") as stream:
                stream.write(fixture["vcf"])
        elif fixture.get("bgzf"):
            import pysam
            plain = vcf
            plain.write_text(fixture["vcf"])
            vcf = Path(str(plain) + ".gz")
            pysam.tabix_compress(str(plain), str(vcf))
        else:
            vcf.write_text(fixture["vcf"])
        reference, converted, native = root / "reference.pkl.gz", root / "reference.rkyv", root / "native.rkyv"
        # Input/producer errors must leave a previous published archive intact.
        sentinel = b"previous published archive"
        native.write_bytes(sentinel)
        python_status = invoke([sys.executable, str(args.python_producer), "-v", str(vcf), "-t", str(table), "-o", str(reference)], root, "python")
        native_status = invoke([str(args.binary), "clingen-map", "build", "--clinvar-vcf", str(vcf), "--clingen-tab", str(table), "--output", str(native)], root, "native")
        row = {"python_status": python_status, "native_status": native_status, "passed": False}
        if fixture.get("fails"):
            row["passed"] = python_status == native_status == 1 and not reference.exists() and native.read_bytes() == sentinel
        elif python_status == 0:
            conversion = invoke([str(args.binary), "cache-convert", "clingen_map", "--pickle", str(reference), "--output", str(converted)], root, "conversion")
            row["conversion_status"] = conversion
            if fixture.get("unsupported"):
                row["passed"] = conversion != 0 and native_status == 1 and native.read_bytes() == sentinel and not converted.exists()
            elif conversion == native_status == 0:
                with redirect_stdout(io.StringIO()):
                    comparison = oracle.compare(SimpleNamespace(binary=args.binary, cache="clingen_map", pickle=reference,
                        archive=native, report=root / "compare.json", mismatches=root / "mismatches.jsonl"))
                    converted_comparison = oracle.compare(SimpleNamespace(binary=args.binary, cache="clingen_map", pickle=reference,
                        archive=converted, report=root / "converted_compare.json", mismatches=root / "converted_mismatches.jsonl"))
                row.update(passed=comparison["passed"] and converted_comparison["passed"],
                           compared_keys=comparison["exported_entries"], mismatches=comparison["mismatches"])
                versions = [int.from_bytes(path.read_bytes()[16:20], "little") for path in (native, converted)]
                row["schema_versions"] = versions
                row["passed"] = row["passed"] and versions == [fixture.get("version", 1)] * 2
        results[name] = row
    stable = {name: oracle.sha256(Path(record["path"])) == record["sha256"] for name, record in sources.items()}
    passed = all(row["passed"] for row in results.values()) and all(stable.values())
    report = {"passed": passed, "cases": results, "sources": sources, "source_stable": stable}
    (args.work_root / "edges.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"passed": passed, "cases": len(results), "failed_cases": [name for name, row in results.items() if not row["passed"]],
                      "report": str(args.work_root / "edges.json")}))
    if not passed:
        raise SystemExit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("edges")
    for name in ("binary", "python-producer", "work-root"):
        command.add_argument("--" + name, type=Path, required=True)
    command.set_defaults(function=edges)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
