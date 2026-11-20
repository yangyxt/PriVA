#!/usr/bin/env python3
"""Exact exon-AF producer fixtures against the unchanged Python procedure.

python clinvar_exon_af_compare.py check --producer ORIGINAL.py --binary PRIVA --work-root NEW_DIR
The reference source controls success, exit code and exception category. Shared
VCF fixture construction and independent cache decoding are reused from the gene
and archive tests; exon policy expectations are produced only by Python.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
from types import SimpleNamespace
sys.dont_write_bytecode = True

PRODUCER_SHA256 = "9387b36b80b57606447744c4d4cb7aac73d85e8b83585a90a554a58b5f27c861"
SENTINEL = b"previous output must survive pre-publication failure\n"


def module(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(name + ".py"))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def row(pos=10, chrom="chr1", af="0.1", joint="0.01", csq="ENST_ONE|1/3", sig="Pathogenic", review="reviewed_by_expert_panel"):
    return f"{chrom}\t{pos}\t.\tA\tT\t.\tPASS\tCLNSIG={sig};CLNREVSTAT={review};AF_grpmax_joint={af};AF_joint={joint};CSQ={csq}\n"


def cases():
    result = {
        "normal": {}, "csi": {"csi": True}, "two_threads": {"threads": "2"},
        "zero_threads": {"threads": "0"}, "negative_threads": {"threads": "-1"}, "bad_threads": {"threads": "bad"},
        "no_index": {"index": False}, "ordinary_gzip": {"compression": "gzip", "index": False},
        "plain_unindexed": {"compression": "plain", "index": False},
        "empty": {"rows": []}, "no_contigs": {"contigs": []},
        "empty_no_contigs_no_index": {"rows": [], "contigs": [], "index": False},
        "empty_header_contig": {"contigs": ["chr1", "chr2"]},
        "header_order": {"contigs": ["chr2", "chr1"], "rows": [row(csq="T|1/2"), row(chrom="chr2", csq="T|2/2", af="0.2")]},
        "global_stats": {"contigs": ["chr2", "chr1"], "rows": [row(csq="T|1/2"), row(pos=11,csq="T|1/2",af="0.3"),row(chrom="chr2",csq="T|1/2",af="0.7")]},
        "same_record_last_exon": {"rows": [row(csq="T|1/3,U|1/2,T|2/3,T|,V|1/1,U|2/2")]},
        "exon_insertion_order": {"rows": [row(csq="T|2/9"),row(pos=11,csq="T|1/7"),row(pos=12,csq="T|2/9",af="0.2")]},
        "duplicates": {"rows": [row(),row(),row(pos=11)]},
        "unicode": {"rows": [row(csq="ENST_α|β/2")]},
        "empty_exon": {"rows": [row(csq="T|")]}, "empty_transcript": {"rows": [row(csq="|1/3")]},
        "short_csq": {"rows": [row(csq="T")]}, "csq_dot": {"rows": [row(csq=".")]},
        "csq_empty": {"rows": [row(csq="")]}, "csq_missing_tuple_member": {"rows": [row(csq="T|1/2,,U|1/1")]},
        "csq_scalar": {"csq_number": "1"}, "csq_integer": {"csq_type": "Integer"},
        "csq_number_g": {"csq_number": "G"}, "csq_flag": {"csq_type": "Flag", "csq_number": "0"},
        "feature_alias": {"description": "Format: Feature|EXON"},
        "lowercase_exon": {"description": "Format: Transcript|exon"},
        "feature_overrides_transcript": {"description": "Format: Transcript|Feature|exon", "rows": [row(csq="OLD|NEW|2/4")]},
        "literal_suffix": {"description": "Format: Transcript|EXON tail"},
        "format_no_space": {"description": "Format:Transcript|EXON"},
        "no_format": {"description": "No annotation fields"},
        "missing_transcript": {"description": "Format: Gene|EXON"},
        "missing_exon": {"description": "Format: Transcript|Other"},
        "duplicate_transcript_header": {"description": "Format: Transcript|EXON|Transcript", "rows": [row(csq="OLD|1/2|NEW")]},
        "duplicate_exon_header": {"description": "Format: Transcript|EXON|EXON", "rows": [row(csq="T|1/3|2/3")]},
        "multiple_format_markers": {"description": "Format: Transcript|EXON|Format: Other"},
        "af_integer": {"af_type": "Integer", "rows": [row(af="3",joint="1")]},
        "af_integer_scalar": {"af_type": "Integer", "af_number": "1", "rows": [row(af="3",joint="1")]},
        "af_float_scalar": {"af_number": "1"},
        "af_string_tuple": {"af_type": "String"},
        "af_string_scalar": {"af_type": "String", "af_number": "1"},
        "af_flag": {"af_type": "Flag", "af_number": "0", "rows": [row().replace("AF_grpmax_joint=0.1", "AF_grpmax_joint")]},
        "af_number_g_skip": {"af_number": "G"},
        "af_empty_tuple": {"rows": [row().replace("AF_grpmax_joint=0.1", "AF_grpmax_joint")]},
        "af_empty_scalar": {"af_number": "1", "rows": [row().replace("AF_grpmax_joint=0.1", "AF_grpmax_joint")]},
        "af_string_empty_tuple": {"af_type": "String", "rows": [row(af="")]},
        "af_scalar_missing_fallback": {"af_number": "1", "rows": [row(af=".")]},
        "undefined_af": {"omit_af_header": True},
        "undefined_af_absent": {"omit_af_header": True, "rows": [row().replace(";AF_grpmax_joint=0.1;AF_joint=0.01", "")]},
        "missing_both_af": {"rows": [row().replace(";AF_grpmax_joint=0.1;AF_joint=0.01", "")]},
        "missing_group_af": {"rows": [row().replace(";AF_grpmax_joint=0.1", "")]},
        "missing_joint_af": {"rows": [row().replace(";AF_joint=0.01", "")]},
        "missing_csq": {"rows": [row().split(";CSQ=")[0] + "\n"]},
        "missing_sig": {"rows": [row().replace("CLNSIG=Pathogenic;", "")]},
        "missing_review": {"rows": [row().replace("CLNREVSTAT=reviewed_by_expert_panel;", "")]},
        "duplicate_info_first": {"rows": [row().rstrip()+";AF_grpmax_joint=0.9;CSQ=OTHER|4/4\n"]},
        "invalid_unused_gt": {"sample": "GT", "sample_value": "x", "rows": [row(sig="Benign")]},
        "valid_gt": {"sample": "GT", "sample_value": "0/1"},
        "invalid_unused_integer_format": {"sample": "DP", "sample_value": "oops", "format_type": "Integer", "rows": [row(sig="Benign")]},
        "invalid_pos": {"mutate": ("chr1\t10\t", "chr1\tbad\t")},
        "negative_pos": {"mutate": ("chr1\t10\t", "chr1\t-1\t")},
        "empty_pos": {"mutate": ("chr1\t10\t", "chr1\t\t")},
        "truncated_record": {"rows": ["chr1\t10\t.\tA\tT\n"]},
        "crlf": {"crlf": True}, "literal_paths": {"filename": "input ;literal $(unused).vcf.gz"},
        "missing_file": {"missing_file": True}, "output_directory": {"output_directory": True},
        "output_missing_parent": {"output_missing_parent": True},
        "crc_failure": {"damage": "crc"}, "truncated_bgzf": {"damage": "truncate"},
        "missing_bgzf_eof": {"damage": "no_eof"}, "trailing_bgzf_garbage": {"damage": "trailing"},
        "skip_bad_then_valid": {"af_type": "String", "rows": [row(af="bad"), row(pos=11,af="0.3")]},
        "fatal_after_prefix": {"rows": [row(),row(pos=11,af=".")]},
    }
    for type_ in ["Float", "String", "Integer"]:
        for token in ["0", "-0", "-1", ".", "NaN", "-NaN", "Inf", "-Inf", "1e-50", "1e40", "x", "0.123456789", "0.1junk", " 0.2", ".2", "0x1p-2", "1_2", "2147483648"]:
            result[f"af_{type_}_{len(result)}"] = {"af_type": type_, "rows": [row(af=token)]}
    for token in ["0", "-1", ".", "NaN", "-NaN", "Inf", "-Inf"]:
        result[f"joint_{len(result)}"] = {"rows": [row(af="0", joint=token)]}
    for sig in ["Benign", "pathogenic", "PATHOGENIC", "Likely_pathogenic", "not_Pathogenic", "."]:
        result[f"significance_{len(result)}"] = {"rows": [row(sig=sig)]}
    for review in ["no_assertion_criteria_provided", "REVIEWED_BY_EXPERT_PANEL", "criteria_provided,_multiple_submitters,_no_conflicts", "xpractice_guidelinex"]:
        result[f"review_{len(result)}"] = {"rows": [row(review=review)]}
    for name in ["CLNSIG", "CLNREVSTAT"]:
        for kind, number in [("String", "1"), ("Integer", "."), ("String", "G")]:
            result[f"{name}_{kind}_{number}"] = {"header_replace": (f"ID={name},Number=.,Type=String", f"ID={name},Number={number},Type={kind}")}
    for n in [7, 8, 9, 127, 128, 129, 1803, 8193]:
        result[f"many_samples_{n}"] = {"rows": [row(pos=i+1, af=str((i % 31)/31)) for i in range(n)]}
    result["mixed_nan"] = {"rows": [row(pos=i+1, af=value) for i, value in enumerate(["NaN", "-NaN", "0.1", "NaN", "-NaN", "Inf", "NaN", "0.2"])]}
    return result


def make_case(folder, case, shared):
    import pysam
    spec = {"description": "Format: Transcript|EXON", "rows": [row()], **case}
    vcf = shared.make_case(folder, spec)
    if "header_replace" in case:
        plain = folder / "input.vcf"
        plain.write_text(plain.read_text().replace(*case["header_replace"]))
        pysam.tabix_compress(str(plain), str(vcf), force=True)
        pysam.tabix_index(str(vcf), preset="vcf", force=True)
    if case.get("damage"):
        content = bytearray(vcf.read_bytes())
        block_size = int.from_bytes(content[16:18], "little") + 1
        if case["damage"] == "crc": content[block_size - 8] ^= 1
        elif case["damage"] == "no_eof": content = content[:-28]
        elif case["damage"] == "trailing": content.extend(b"trailing garbage")
        else: content = content[:block_size - 10]
        vcf.write_bytes(content)
    if case.get("missing_file"):
        return folder / "absent.vcf.gz"
    return vcf


def run(command, folder, name):
    with (folder / (name + ".stdout.log")).open("x") as out, (folder / (name + ".stderr.log")).open("x") as err:
        status = subprocess.run(command, stdout=out, stderr=err, timeout=90).returncode
    text = (folder / (name + ".stderr.log")).read_text(errors="replace")
    matches = re.findall(r"^(?:priva: )?([A-Za-z_][A-Za-z0-9_.]*Error|error):", text, flags=re.M)
    return {"exit": status, "exception": matches[-1].split(".")[-1] if matches else None}


def check(args):
    if hashlib.sha256(args.producer.read_bytes()).hexdigest() != PRODUCER_SHA256:
        raise RuntimeError("reference producer identity differs")
    args.work_root.mkdir(parents=True, exist_ok=False)
    shared, comparator = module("clinvar_gene_af_compare"), module("cache_values_compare")
    reports = {}
    for name, case in cases().items():
        folder = args.work_root / name
        folder.mkdir()
        vcf = make_case(folder, case, shared)
        reference, archive = folder / "python.pkl", folder / "native.rkyv"
        if case.get("output_directory"):
            reference.mkdir(); archive.mkdir()
        elif case.get("output_missing_parent"):
            reference, archive = folder / "absent/python.pkl", folder / "absent/native.rkyv"
        else:
            reference.write_bytes(SENTINEL)
            archive.write_bytes(SENTINEL)
        threads = case.get("threads", "1")
        py = run([sys.executable, str(args.producer), "--vcf_path", str(vcf), "--output", str(reference), "--threads", threads], folder, "python")
        native = run([str(args.binary), "clinvar-exon-af", "build", "--vcf", str(vcf), "--output", str(archive), "--threads", threads], folder, "native")
        py["written"], native["written"] = [p.is_file() and p.read_bytes() != SENTINEL for p in [reference, archive]]
        report = {"python": py, "native": native, "passed": py == native}
        if py["exit"] == 0 and native["exit"] == 0:
            comparison = comparator.compare(SimpleNamespace(binary=args.binary, cache="clinvar_patho_exon_af_stat", pickle=reference,
                archive=archive, report=folder/"values.json", mismatches=folder/"mismatches.jsonl"))
            report["comparison"] = comparison
            report["passed"] &= comparison["passed"]
        if name in ["normal", "csi", "header_order", "global_stats", "same_record_last_exon", "exon_insertion_order", "empty_header_contig", "many_samples_129", "mixed_nan"]:
            import pysam
            with pysam.VariantFile(str(vcf)) as source:
                contigs = list(source.header.contigs)
            command = [str(args.binary), "clinvar-exon-af", "merge"]
            statuses = []
            for at, contig in enumerate(contigs):
                shard = folder / f"chr_{at}.rkyv"
                statuses.append(run([str(args.binary), "clinvar-exon-af", "chromosome", "--vcf", str(vcf), "--chrom", contig, "--output", str(shard)], folder, f"chr_{at}"))
                command += ["--input", str(shard)]
            merged = folder / "merged.rkyv"
            merged.write_bytes(SENTINEL)
            status = run(command + ["--output", str(merged)], folder, "merge")
            pipeline = {"shards": statuses, "merge": status, "passed": not any(s["exit"] for s in statuses) and status["exit"] == 0}
            if pipeline["passed"]:
                comparison = comparator.compare(SimpleNamespace(binary=args.binary, cache="clinvar_patho_exon_af_stat", pickle=reference,
                    archive=merged, report=folder/"merged.values.json", mismatches=folder/"merged.mismatches.jsonl"))
                pipeline["comparison"] = comparison
                pipeline["passed"] &= comparison["passed"]
            report["pipeline"] = pipeline
            report["passed"] &= pipeline["passed"]
        reports[name] = report
    admissions = {}
    for name, extra in [("missing_required", []), ("missing_output", ["--vcf", "missing.vcf"]),
                        ("unknown_option", ["--unknown"]), ("help", ["--help"])]:
        folder = args.work_root / ("admission_" + name)
        folder.mkdir()
        py_args = ["--vcf_path" if x == "--vcf" else x for x in extra]
        py = run([sys.executable, str(args.producer)] + py_args, folder, "python")
        native = run([str(args.binary), "clinvar-exon-af", "build"] + extra, folder, "native")
        admissions[name] = {"python": py, "native": native, "passed": py == native}
    result = {"count": len(reports), "cases": reports, "admission": admissions, "producer_sha256": PRODUCER_SHA256,
              "failed": [name for name, r in reports.items() if not r["passed"]]}
    result["failed"] += ["admission_" + name for name, r in admissions.items() if not r["passed"]]
    result["passed"] = not result["failed"]
    (args.work_root / "comparison.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k:v for k,v in result.items() if k != "cases"}), flush=True)
    return 0 if result["passed"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("check")
    for name in ["producer", "binary", "work-root"]:
        command.add_argument("--" + name, type=Path, required=True)
    command.set_defaults(function=check)
    args = parser.parse_args()
    return args.function(args)


if __name__ == "__main__":
    raise SystemExit(main())
