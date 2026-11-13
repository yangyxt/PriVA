#!/usr/bin/env python3
"""Compare AA/splice producers against the unchanged Python implementation.

probe measures source-supported shapes and existing schema refusals; check runs
the native builder and compares all values with an independent Python decoder.
All fixtures, outputs, logs and mismatch ledgers stay in --work-root.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import pickle
import re
import shutil
import subprocess
import sys
from types import SimpleNamespace

PRODUCER_SHA256 = "c222bf399eb3b7d1d948c733d6b26d5d17a110d5460478ae082adce3215e0b89"
SENTINEL = b"prior output survives a processing failure\n"
BASE_FIELDS = ["Protein_position", "HGVSp", "Consequence", "Feature_type", "Feature", "EXON", "INTRON", "HGVSc"]
AI_FIELDS = ["SpliceAI_pred_" + name for name in ["SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"]]
FIELDS = BASE_FIELDS + AI_FIELDS + ["TAIL"]
DEFAULTS = dict(zip(FIELDS, ["2", "p.A2T", "missense_variant", "Transcript", "ENST_ONE", "1/2", "", "ENST_ONE:c.2+1A>T",
                            "GENE", "0.2", "0.3", "0.4", "0.5", "1", "2", "3", "4", "end"]))


def annotation(fields=FIELDS, **values):
    data = dict(DEFAULTS, **values)
    return "|".join(data.get(name, "value") for name in fields)


def row(pos=10, chrom="chr1", csq=None, alt="T", ref="A", sig="Pathogenic", review="reviewed_by_expert_panel", extra=""):
    if csq is None:
        csq = annotation()
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tCLNSIG={sig};CLNREVSTAT={review};CSQ={csq}{extra}\n"


def cases():
    result = {"normal": {}, "all_splice_empty": {"values": {name: "" for name in AI_FIELDS}},
              "partial_splice": {"values": {"SpliceAI_pred_DS_AL": ""}},
              "partial_splice_diff_rows": {"rows": [row(csq=annotation(SpliceAI_pred_DS_AL="")), row(pos=11, csq=annotation(SpliceAI_pred_DS_AG=""))]},
              "extra_splice_field": {"fields": FIELDS[:-1] + ["SpliceAI_pred_NEW", "TAIL"]},
              "trailing_splice_header": {"fields": FIELDS[:-1]},
              "alt_missing": {"rows": [row(alt=".")]},
              "empty_ref_alt": {"mutate": ("\tA\tT\t", "\t\t\t")},
              "no_feature_header": {"fields": [v for v in FIELDS if v != "Feature"]},
              "no_aa_headers": {"fields": [v for v in FIELDS if v not in ["Protein_position", "HGVSp"]]},
              "no_splice_headers": {"fields": [v for v in FIELDS if v not in ["EXON", "INTRON"]]},
              "no_format": {"description": "no format annotation"},
              "other_last_format": {"extra_header": '##INFO=<ID=OTHER,Number=.,Type=String,Description="Format: Feature|EXON|TAIL">\n'},
              "duplicate_header_names": {"fields": FIELDS[:-1] + ["Feature", "SpliceAI_pred_DS_AG", "TAIL"]},
              "csi": {"csi": True}, "no_index": {"index": False}, "no_eof": {"no_eof": True},
              "ordinary_gzip": {"compression": "gzip", "index": False}, "plain_unindexed": {"compression": "plain", "index": False},
              "empty": {"rows": []}, "no_header_contigs": {"contigs": []},
              "empty_no_contigs_no_index": {"rows": [], "contigs": [], "index": False},
              "header_order": {"contigs": ["chr2", "chr1"], "rows": [row(), row(chrom="chr2",csq=annotation(HGVSp="p.A2G"))]},
              "duplicates": {"rows": [row(csq=annotation()+","+annotation()), row(pos=11)]},
              "multiallelic": {"rows": [row(alt="T,G")]},
              "short_annotation": {"rows": [row(csq="short")]},
              "unused_invalid_gt": {"sample": "GT", "sample_value": "bad", "rows": [row(sig="Benign")]},
              "unused_invalid_format": {"sample": "DP", "sample_value": "bad", "format_type": "Integer", "rows": [row(sig="Benign")]},
              "sig_scalar": {"sig_number": "1"}, "review_scalar": {"review_number": "1"},
              "sig_integer": {"sig_type": "Integer"}, "sig_flag": {"sig_type": "Flag", "sig_number": "0"},
              "csq_scalar": {"csq_number": "1"}, "csq_integer": {"csq_type": "Integer"},
              "csq_flag": {"csq_type": "Flag", "csq_number": "0"},
              "csq_number_g": {"csq_number": "G"},
              "csq_missing": {"rows": [row().split(";CSQ=")[0]+"\n"]},
              "csq_empty": {"rows": [row(csq="")]},
              "csq_dot": {"rows": [row(csq=".")]},
              "no_clinical_info": {"rows": ["chr1\t10\t.\tA\tT\t.\tPASS\t.\n"]},
              "zero_threads": {"threads": "0"}, "negative_threads": {"threads": "-1"},
              "bad_threads": {"threads": "bad"}, "two_threads": {"threads": "2"},
              "crlf": {"crlf": True}, "unicode": {"values": {"Feature": "ENST_α", "HGVSp": "p.Δ"}},
              "literal_paths": {"filename": "input ;literal $(unused).vcf.gz"}}
    for field in ["EXON", "INTRON", "HGVSc", "Consequence"]:
        result["missing_" + field.lower()] = {"fields": [v for v in FIELDS if v != field]}
    result.update({
        "missing_exon_present_intron": {"fields": [v for v in FIELDS if v != "EXON"], "values": {"INTRON": "1/2"}},
        "empty_header_contig": {"contigs": ["chr1", "chr2"]},
        "negative_pos": {"mutate": ("chr1\t10\t", "chr1\t-1\t")},
        "invalid_pos": {"mutate": ("chr1\t10\t", "chr1\tbad\t")},
        "empty_pos": {"mutate": ("chr1\t10\t", "chr1\t\t")},
        "first_alt_dot": {"rows": [row(alt=".,T")]},
        "first_alt_empty": {"rows": [row(alt=",T")]},
        "empty_reference": {"mutate": ("\tA\tT\t", "\t\tT\t")},
        "csq_none_after_valid": {"rows": [row(csq=annotation()+",,"+annotation())]},
        "csq_no_raw_value": {"rows": [row().split(";CSQ=")[0]+";CSQ\n"]},
        "csq_empty_scalar": {"csq_number": "1", "rows": [row(csq="")]},
        "review_number_g": {"review_number": "G"},
        "sig_number_g": {"sig_number": "G"},
        "duplicate_info_first": {"rows": [row(extra=";CLNSIG=Benign;CSQ=short")]},
        "duplicate_csq_header": {"extra_header": '##INFO=<ID=CSQ,Number=1,Type=Integer,Description="Format: OTHER|TAIL">\n'},
        "sig_undefined": {"omit_headers": ["CLNSIG"]},
        "review_undefined": {"omit_headers": ["CLNREVSTAT"]},
        "csq_undefined": {"omit_headers": ["CSQ"]},
        "two_threads_spaced": {"threads": " +2 "},
        "two_threads_underscore": {"threads": "0_2"},
        "two_threads_unicode": {"threads": "٢"},
        "thread_invalid_whitespace": {"threads": "\x1c2"},
        "aa_only": {"splice_output": "empty"},
        "splice_only": {"aa_output": "empty"},
        "neither_output": {"aa_output": "empty", "splice_output": "empty"},
        "aa_bad_parent": {"aa_output": "bad_parent"},
        "splice_bad_parent": {"splice_output": "bad_parent"},
    })
    for n in [5, 6, 7]:
        result["short_record_"+str(n)] = {"mutate": (row(), "\t".join(row().rstrip().split("\t")[:n])+"\n")}
    for key in ["CLNSIG", "CLNREVSTAT"]:
        for raw in [".", "", "Benign", "Pathogenic,Likely_pathogenic"]:
            result[key.lower()+"_"+str(len(result))] = {"rows": [row(**{"sig" if key == "CLNSIG" else "review": raw})]}
    return result


def make_case(folder, spec):
    import pysam
    fields = spec.get("fields", FIELDS)
    lines = ["##fileformat=VCFv4.2\n"] + [f"##contig=<ID={c},length=1000000>\n" for c in spec.get("contigs", ["chr1"])]
    for key, prefix in [("CLNSIG", "sig"), ("CLNREVSTAT", "review"), ("CSQ", "csq")]:
        if key in spec.get("omit_headers", []):
            continue
        description = spec.get("description", "Format: " + "|".join(fields)) if key == "CSQ" else key
        lines.append(f'##INFO=<ID={key},Number={spec.get(prefix+"_number", ".")},Type={spec.get(prefix+"_type", "String")},Description="{description}">\n')
    lines.append(spec.get("extra_header", ""))
    if spec.get("sample"):
        lines.append(f'##FORMAT=<ID={spec["sample"]},Number=1,Type={spec.get("format_type", "String")},Description="sample">\n')
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" + ("\tFORMAT\tSAMPLE" if spec.get("sample") else "") + "\n")
    for item in spec.get("rows", [row(csq=annotation(fields, **spec.get("values", {})))]):
        if spec.get("sample"):
            item = item.rstrip("\n") + "\t" + spec["sample"] + "\t" + spec["sample_value"] + "\n"
        lines.append(item)
    data = "".join(lines)
    if spec.get("crlf"):
        data = data.replace("\n", "\r\n")
    plain = folder / "input.vcf"
    plain.write_text(data, newline="")
    if spec.get("compression") == "plain":
        return plain
    output = folder / spec.get("filename", "input.vcf.gz")
    if spec.get("compression") == "gzip":
        with gzip.open(output, "wb") as stream:
            stream.write(data.encode())
        return output
    pysam.tabix_compress(str(plain), str(output), force=True)
    if spec.get("index", True):
        pysam.tabix_index(str(output), preset="vcf", force=True, csi=spec.get("csi", False))
    if spec.get("mutate"):
        plain.write_text(data.replace(*spec["mutate"]), newline="")
        pysam.tabix_compress(str(plain), str(output), force=True)
    if spec.get("no_eof"):
        output.write_bytes(output.read_bytes()[:-28])
    return output


def run(command, folder, name):
    with (folder / (name + ".stdout.log")).open("x") as out, (folder / (name + ".stderr.log")).open("x") as err:
        return subprocess.run(list(map(str, command)), stdout=out, stderr=err).returncode


def failure_category(folder, name):
    """Compare actual exception types, without requiring Python tracebacks."""
    text = (folder / (name + ".stderr.log")).read_text()
    if name == "native":
        found = re.findall(r"^priva: ([A-Za-z_][A-Za-z_0-9]*):", text, re.M)
    else:
        found = re.findall(r"^([A-Za-z_][A-Za-z_0-9]*(?:Error|Exception)|error):", text, re.M)
    return found[-1] if found else None


def consumer_probe(records, scripts):
    sys.path.insert(0, str(scripts))
    from acmg_ps1_pm5_known_variant import check_splice_pathogenic
    target = {"INTRON": "", "EXON": "1/2", "SpliceAI_pred_DS_AG": .9, "SpliceAI_pred_DS_AL": .2,
              "SpliceAI_pred_DS_DG": .1, "SpliceAI_pred_DS_DL": .1, "HGVSc": "ENST_ONE:c.2+1A>T",
              "STRAND": 1, "pos": 10, "ref": "A", "alt": "T", "chrom": "chr1"}
    result = {}
    for name, update in [("matching", {}), ("nonmatching", {"EXON": "other", "INTRON": "other"}),
                         ("below_score", {"SpliceAI_pred_DS_AG": .01, "SpliceAI_pred_DS_AL": .01,
                                          "SpliceAI_pred_DS_DG": .01, "SpliceAI_pred_DS_DL": .01})]:
        try:
            result[name] = {"result": check_splice_pathogenic(dict(target, **update), records, 1, .5)}
        except Exception as error:
            result[name] = {"exception": type(error).__name__, "message": str(error)}
    return result


def check(args):
    if hashlib.sha256(args.producer.read_bytes()).hexdigest() != PRODUCER_SHA256:
        raise RuntimeError("unchanged Python producer identity differs")
    args.work_root.mkdir(parents=True, exist_ok=False)
    reports = {}
    if args.command == "check":
        import cache_values_compare as comparator
    for name, case in cases().items():
        folder = args.work_root / name
        folder.mkdir()
        vcf = make_case(folder, case)
        aa, splice = folder / "python.aa.pkl", folder / "python.splice.pkl"
        for path in [aa, splice]:
            path.write_bytes(SENTINEL)
        def target(path, kind):
            policy = case.get(kind+"_output")
            return "" if policy == "empty" else folder/"absent"/path.name if policy == "bad_parent" else path
        status = run([sys.executable, args.producer, vcf, target(aa,"aa"), target(splice,"splice"), case.get("threads", "1")], folder, "python")
        report = {"python_exit": status, "written": [p.read_bytes() != SENTINEL for p in [aa, splice]]}
        if not status and all(report["written"]):
            aa_value, splice_value = [pickle.loads(p.read_bytes()) for p in [aa, splice]]
            records = [v for values in splice_value.values() for v in values]
            report.update(aa_transcripts=len(aa_value), splice_transcripts=len(splice_value), splice_records=len(records),
                          nullable_fields=sorted({k for r in records for k, v in r.items() if v is None}),
                          splice_key_counts=sorted({len(r["splice_ai"]) for r in records}))
            if args.command == "probe":
                report["consumer"] = consumer_probe(records, args.scripts)
                report["conversion"] = {}
                for cache, source in [("clinvar_aa_stat", aa), ("clinvar_splice_stat", splice)]:
                    report["conversion"][cache] = run([args.binary, "cache-convert", cache, "--pickle", source,
                        "--output", folder / (cache + ".rkyv")], folder, cache)
        if args.command == "check":
            native = [folder / "native.aa.rkyv", folder / "native.splice.rkyv"]
            for path in native:
                path.write_bytes(SENTINEL)
            native_status = run([args.binary, "clinvar-aa-splice", "build", "--vcf", vcf,
                "--aa-output", target(native[0],"aa"), "--splice-output", target(native[1],"splice"), "--threads", case.get("threads", "1")], folder, "native")
            report.update(native_exit=native_status, native_written=[p.read_bytes() != SENTINEL for p in native])
            report["python_category"] = failure_category(folder, "python") if status else None
            report["native_category"] = failure_category(folder, "native") if native_status else None
            report["failure_category_match"] = (status == 0 or
                (report["python_category"] is not None and report["python_category"] == report["native_category"]))
            report["passed"] = (status == native_status and report["failure_category_match"] and
                                report["written"] == report["native_written"])
            if any(report["written"]) and any(report["native_written"]):
                report["values"] = {}
                for cache, source, archive in zip(["clinvar_aa_stat", "clinvar_splice_stat"], [aa, splice], native):
                    if source.read_bytes() == SENTINEL or archive.read_bytes() == SENTINEL:
                        continue
                    result = comparator.compare(SimpleNamespace(binary=args.binary, cache=cache, pickle=source, archive=archive,
                        report=folder/(cache+".values.json"), mismatches=folder/(cache+".mismatches.jsonl")))
                    report["values"][cache] = result
                    report["passed"] &= result["passed"]
                import struct
                if native[1].read_bytes() != SENTINEL:
                    report["splice_schema_version"] = struct.unpack("<I", native[1].read_bytes()[16:20])[0]
            if name in ["normal", "csi", "header_order", "duplicates", "partial_splice_diff_rows", "alt_missing", "missing_hgvsc", "no_header_contigs"]:
                import pysam
                with pysam.VariantFile(str(vcf)) as handle:
                    chromosomes = list(handle.header.contigs)
                merge_command = [args.binary, "clinvar-aa-splice", "merge"]
                shard_exits = []
                for i, chrom in enumerate(chromosomes):
                    shard_aa, shard_splice = folder/f"chr{i}.aa.rkyv", folder/f"chr{i}.splice.rkyv"
                    shard_exits.append(run([args.binary, "clinvar-aa-splice", "chromosome", "--vcf", vcf, "--chrom", chrom,
                        "--aa-output", shard_aa, "--splice-output", shard_splice], folder, f"chr{i}"))
                    merge_command += ["--aa-input", shard_aa, "--splice-input", shard_splice]
                merged_aa, merged_splice = folder/"merged.aa.rkyv", folder/"merged.splice.rkyv"
                code = run(merge_command + ["--aa-output", merged_aa, "--splice-output", merged_splice], folder, "merge")
                pipeline = {"shard_exits": shard_exits, "merge_exit": code, "passed": not any(shard_exits) and code == 0}
                if pipeline["passed"]:
                    for cache, source, archive in zip(["clinvar_aa_stat", "clinvar_splice_stat"], [aa, splice], [merged_aa, merged_splice]):
                        result = comparator.compare(SimpleNamespace(binary=args.binary, cache=cache, pickle=source, archive=archive,
                            report=folder/(cache+".merged.values.json"), mismatches=folder/(cache+".merged.mismatches.jsonl")))
                        pipeline[cache] = result
                        pipeline["passed"] &= result["passed"]
                report["pipeline"] = pipeline
                report["passed"] &= pipeline["passed"]
        reports[name] = report
    result = {"producer_sha256": PRODUCER_SHA256, "count": len(reports), "cases": reports}
    result["failed"] = [name for name, case in reports.items() if not case.get("passed", True)]
    (args.work_root / "comparison.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"count": len(reports), "python_success": sum(r["python_exit"] == 0 for r in reports.values()), "failed": result["failed"]}))
    if result["failed"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for command in ["probe", "check"]:
        p = commands.add_parser(command)
        for name in ["producer", "scripts", "binary", "work-root"]:
            p.add_argument("--" + name, type=Path, required=True)
        p.set_defaults(function=check)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
