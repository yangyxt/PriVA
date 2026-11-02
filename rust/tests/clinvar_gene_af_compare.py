#!/usr/bin/env python3
"""Differential gene-AF producer fixtures against a pinned unchanged Python source.

  python clinvar_gene_af_compare.py probe --producer ORIGINAL.py --work-root NEW_DIR
  python clinvar_gene_af_compare.py check --producer ORIGINAL.py --binary PRIVA --work-root NEW_DIR
The source-matched producer controls outcomes. Every surviving value is compared
with CPython's independent decoder, including float bits and record field order.
All raw data and complete mismatch ledgers remain in --work-root.
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

PRODUCER_SHA256 = "0e9dce6d1cd9e193e1549b74ec3ccf83786342c51e80536c050aab816b9710bb"
SENTINEL = b"previous output must survive a pre-publication failure\n"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def row(pos=10, gene="ENSG_ONE", af="0.1", joint="0.01", chrom="chr1", csq=None, sig="Pathogenic", review="reviewed_by_expert_panel", alt="T", ref="A", extras=""):
    if csq is None:
        csq = f"{gene}|p.A1T|MODERATE|missense_variant|GENE"
    info = f"CLNSIG={sig};CLNREVSTAT={review};AF_grpmax_joint={af};AF_joint={joint};CSQ={csq}"
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}{extras}\n"


def cases():
    result = {"normal": {}, "csi": {"csi": True}, "csi_large_position": {"csi": True, "rows": [row(pos=600000000)]}, "no_index": {"index": False},
              "ordinary_gzip": {"compression": "gzip", "index": False},
              "plain_unindexed": {"compression": "plain", "index": False},
              "no_csq_format": {"description": "no field description", "index": False},
              "no_csq_format_threads_zero": {"description": "no field description", "threads": "0", "index": False},
              "no_csq_format_bad_threads": {"description": "no field description", "threads": "bad", "index": False},
              "zero_threads": {"threads": "0"}, "negative_threads": {"threads": "-1"},
              "two_threads": {"threads": "2"}, "empty_records": {"rows": []},
              "no_header_contigs": {"contigs": []}, "empty_header_contig": {"contigs": ["chr1","chr2"]},
              "empty_no_contigs_no_index": {"rows": [], "contigs": [], "index": False},
              "unicode_fields": {"rows": [row(gene="ENSG_α", csq="ENSG_α|p.Δ|MODERATE|missense_variant|μ")]},
              "multiallelic": {"rows": [row(alt="T,G", af="0.1,0.8", joint="0.01,0.9")]},
              "duplicates_partial_sum": {"rows": [row(pos=10), row(pos=11, af="0.2"), row(pos=12, af="0.1",joint="0.5")]},
              "duplicate_csq": {"rows": [row(csq="ENSG_ONE|p.A1T|MODERATE|missense_variant|FIRST,ENSG_ONE|p.A2T|HIGH|stop_gained|SECOND"), row(pos=11,af="0.3")]},
              "same_chr_tie": {"rows": [row(), row(pos=11, joint="0.9")]},
              "header_order_tie": {"contigs": ["chr2","chr1"], "rows": [row(), row(chrom="chr2",joint="0.9")]},
              "cross_chr_larger": {"rows": [row(),row(chrom="chr2",af="0.2")], "contigs": ["chr1","chr2"]},
              "merge_none_first": {"rows": [row(af="0",joint="0"),row(chrom="chr2")],"contigs": ["chr1","chr2"]},
              "merge_none_second": {"rows": [row(),row(chrom="chr2",af="0",joint="0")],"contigs": ["chr1","chr2"]},
              "merge_none_both": {"rows": [row(af="0",joint="0"),row(chrom="chr2",af="0",joint="0")],"contigs": ["chr1","chr2"]},
              "none_separate_genes": {"rows": [row(af="0",joint="0"),row(chrom="chr2",gene="ENSG_TWO",af="0",joint="0")],"contigs": ["chr1","chr2"]},
              "missing_alt_winner": {"rows": [row(alt=".")]},
              "missing_alt_not_winner": {"rows": [row(alt=".",af="0",joint="0")]},
              "missing_alt_tied": {"rows": [row(),row(pos=11,alt=".")]},
              "empty_ref_alt": {"mutate": ("\tA\tT\t", "\t\t\t")},
              "short_csq": {"rows": [row(csq="ENSG_ONE|short")]},
              "csq_missing_entry": {"rows": [row(csq=".,ENSG_ONE|p.A1T|MODERATE|missense_variant|GENE")]},
              "csq_dot": {"rows": [row(csq=".")]}, "csq_empty": {"rows": [row(csq="")]},
              "csq_flag": {"rows": [row().split(";CSQ=")[0]+";CSQ\n"]},
              "csq_scalar": {"csq_number": "1"}, "csq_integer": {"csq_type": "Integer"},
              "csq_duplicate_field": {"description": "Format: Gene|HGVSp|IMPACT|Consequence|SYMBOL|Gene", "rows": [row(csq="WRONG|p.A1T|HIGH|stop_gained|S|ENSG_RIGHT")]},
              "missing_csq_key": {"description": "Format: Gene|HGVSp|IMPACT|Consequence|OTHER"},
              "no_info": {"rows": ["chr1\t10\t.\tA\tT\t.\tPASS\t.\n"]},
              "invalid_unused_gt": {"sample": "GT", "sample_value": "x", "rows": [row(sig="Benign")]},
              "invalid_unused_integer_format": {"sample": "DP", "sample_value": "oops", "format_type": "Integer", "rows": [row(sig="Benign")]},
              "valid_samples": {"sample": "GT", "sample_value": "0/1"},
              "af_integer": {"af_type": "Integer", "rows": [row(af="2",joint="1")]},
              "af_integer_overflow": {"af_type": "Integer", "rows": [row(af="2147483648",joint="1")]},
              "af_integer_prefix": {"af_type": "Integer", "rows": [row(af="2.3x",joint="1")]},
              "af_string": {"af_type": "String"},
              "af_flag": {"af_type": "Flag", "af_number": "0", "rows": [row().replace("AF_grpmax_joint=0.1", "AF_grpmax_joint")]},
              "af_number_g": {"af_number": "G"},
              "af_empty_tuple": {"rows": [row().replace("AF_grpmax_joint=0.1", "AF_grpmax_joint")]},
              "af_equals_empty": {"rows": [row(af="")]},
              "af_empty_scalar": {"af_number": "1", "rows": [row().replace("AF_grpmax_joint=0.1", "AF_grpmax_joint")]},
              "duplicate_info_first": {"rows": [row(extras=";AF_grpmax_joint=0.9;CSQ=ENSG_TWO|p.X|HIGH|stop_gained|OTHER")]},
              "undefined_af": {"omit_af_header": True},
              "undefined_af_absent": {"omit_af_header": True, "rows": [row().replace(";AF_grpmax_joint=0.1;AF_joint=0.01", "")]},
              "missing_joint": {"rows": [row().replace(";AF_joint=0.01", "")]},
              "missing_both_af": {"rows": [row().replace(";AF_grpmax_joint=0.1;AF_joint=0.01", "")]},
              "gene_prefix_only": {"rows": [row(gene="ENSG")]},
              "non_gene": {"rows": [row(gene="OTHER")]},
              "invalid_pos": {"mutate": ("chr1\t10\t", "chr1\txyz\t")},
              "negative_pos": {"mutate": ("chr1\t10\t", "chr1\t-1\t")},
              "empty_pos": {"mutate": ("chr1\t10\t", "chr1\t\t")},
              "truncated_no_info": {"rows": ["chr1\t10\t.\tA\tT\n"]},
              "crlf": {"crlf": True}, "literal_paths": {"filename": "input ;literal $(unused).vcf.gz"},
              }
    # Boolean fields are source-supported; v2 is needed only when one survives.
    result.update({
        "flag_group_float_joint": {"af_type":"Flag", "af_number":"0", "joint_type":"Float", "joint_number":"A"},
        "flag_group_absent_falls_back": {"af_type":"Flag", "af_number":"0", "joint_type":"Float", "joint_number":"A", "rows":[row().replace(";AF_grpmax_joint=0.1", "")]},
        "flag_joint_first_snapshot": {"joint_type":"Flag", "joint_number":"0"},
        "flag_joint_added_twice": {"joint_type":"Flag", "joint_number":"0", "rows":[row(),row(pos=11,af="0.2")]},
        "flag_joint_absent_then_true": {"joint_type":"Flag", "joint_number":"0", "rows":[row().replace(";AF_joint=0.01", ""),row(pos=11,af="0.2")]},
        "flag_joint_true_then_absent": {"joint_type":"Flag", "joint_number":"0", "rows":[row(),row(pos=11,af="0.2").replace(";AF_joint=0.01", "")]},
        "flag_group_later_float_winner": {"af_type":"Flag", "af_number":"0", "joint_type":"Float", "joint_number":"A", "rows":[row(),row(pos=11,joint="2").replace(";AF_grpmax_joint=0.1", "")]},
        "flag_unused_record": {"af_type":"Flag", "af_number":"0", "rows":[row(sig="Benign")]},
        "flag_empty_contig": {"af_type":"Flag", "af_number":"0", "contigs":["chr1","chr2"]},
        "flag_header_order_tie": {"af_type":"Flag", "af_number":"0", "contigs":["chr2","chr1"], "rows":[row(),row(chrom="chr2",ref="C")]},
        "flag_none_merge_failure": {"af_type":"Flag", "af_number":"0", "contigs":["chr1","chr2"], "rows":[row().replace(";AF_grpmax_joint=0.1;AF_joint=0.01", ""),row(chrom="chr2")]},
        "flag_all_absent": {"af_type":"Flag", "af_number":"0", "rows":[row().replace(";AF_grpmax_joint=0.1;AF_joint=0.01", "")]},
    })
    for value in ["0", "-0", "-1", ".", "NaN", "-NaN", "Inf", "-Inf", "1e-50", "1e40", "x", "0.123456789", "0.1junk", " 0.2", ".2", "0x1p-2"]:
        name = "af_token_" + str(len(result))
        result[name] = {"rows": [row(af=value)], "token": value}
    for value in ["0", "-1", ".", "NaN", "Inf", "-Inf"]:
        result["joint_token_"+str(len(result))] = {"rows": [row(joint=value),row(pos=11,af="0.2",joint="0.1")],"token": value}
    for sig in ["Benign", "PATHOGENIC", "pathogenic", "Likely_pathogenic", "not_athogenic", "."]:
        result["sig_"+str(len(result))] = {"rows": [row(sig=sig)]}
    for review in ["no_assertion_criteria_provided", "REVIEWED_BY_EXPERT_PANEL", "criteria_provided,_multiple_submitters,_no_conflicts", "xpractice_guidelinex"]:
        result["review_"+str(len(result))] = {"rows": [row(review=review)]}
    # Append these without renaming the earlier numbered fixture catalogue.
    result.update({"whitespace_threads": {"threads": " 1 "}, "underscore_threads": {"threads": "0_1"},
        "unicode_threads": {"threads": "١"}, "unicode_zero_threads": {"threads": "٠"},
        "no_csq_unicode_zero_threads": {"description": "no field description", "threads": "٠", "index": False}})
    return result


def make_case(folder, spec):
    import pysam
    contigs = spec.get("contigs", ["chr1"])
    lines = ["##fileformat=VCFv4.2\n"] + [f"##contig=<ID={c},length=1000000>\n" for c in contigs]
    lines += ['##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">\n',
              '##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="Review status">\n']
    if not spec.get("omit_af_header"):
        for key in ["AF_grpmax_joint", "AF_joint"]:
            number = spec.get("joint_number", spec.get("af_number", "A")) if key == "AF_joint" else spec.get("af_number", "A")
            kind = spec.get("joint_type", spec.get("af_type", "Float")) if key == "AF_joint" else spec.get("af_type", "Float")
            lines.append(f'##INFO=<ID={key},Number={number},Type={kind},Description="AF">\n')
    lines.append(f'##INFO=<ID=CSQ,Number={spec.get("csq_number", ".")},Type={spec.get("csq_type", "String")},Description="{spec.get("description", "Format: Gene|HGVSp|IMPACT|Consequence|SYMBOL")}">\n')
    if spec.get("sample"):
        lines.append(f'##FORMAT=<ID={spec["sample"]},Number=1,Type={spec.get("format_type", "String")},Description="sample">\n')
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" + ("\tFORMAT\tSAMPLE" if spec.get("sample") else "") + "\n")
    for item in spec.get("rows", [row()]):
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
        data = data.replace(*spec["mutate"])
        plain.write_text(data, newline="")
        pysam.tabix_compress(str(plain), str(output), force=True)
    return output


def run(command, folder, name):
    with (folder / (name+".stdout.log")).open("x") as out, (folder / (name+".stderr.log")).open("x") as err:
        result = subprocess.run(command, stdout=out, stderr=err)
    return result.returncode


def failure_category(path):
    # Multiprocessing tracebacks may repeat the worker exception. The final
    # exception is the actual command failure, not the RemoteTraceback wrapper.
    matches = re.findall(r"(?m)^(?:priva: )?([A-Za-z_][A-Za-z_.0-9]*(?:Error|Exception)): (.*)$", path.read_text())
    return matches[-1][0].rsplit(".", 1)[-1] if matches else None


def matching_outcome(python_exit, native_exit, python_category, native_category):
    return (python_exit == native_exit == 0 or
        python_exit == native_exit == 1 and python_category is not None and python_category == native_category)


def verify_pin(record):
    path = Path(record["path"])
    if path.stat().st_size != record["bytes"] or digest(path) != record["sha256"]:
        raise RuntimeError("saved reference changed: " + str(path))
    return path


def admission(args):
    """Native named-option admission stays exit 2, separate from source calls."""
    args.work_root.mkdir()
    absent = str(args.work_root / "absent.vcf.gz")
    output = str(args.work_root / "must_not_be_written.rkyv")
    common = ["--vcf", absent, "--output", output]
    cases = {
        "unknown_mode": ["unknown"],
        "build_missing_all": ["build"],
        "build_missing_output": ["build", "--vcf", absent],
        "build_missing_vcf": ["build", "--output", output],
        "chromosome_missing_contig": ["chromosome", *common],
        "merge_missing_inputs": ["merge", "--output", output],
        "unknown_option": ["build", *common, "--unknown", "x"],
        "missing_option_value": ["build", *common, "--threads"],
    }
    reports = {}
    for name, options in cases.items():
        status = run([str(args.binary), "clinvar-gene-af", *options], args.work_root, name)
        reports[name] = {"exit": status, "passed": status == 2}
    result = {"cases": reports, "passed": not Path(output).exists() and all(row["passed"] for row in reports.values()),
        "scope": "native named-option admission; the legacy Python positional usage interface is separate"}
    (args.work_root / "admission.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if not result["passed"]:
        raise SystemExit(1)


def check(args):
    if digest(args.producer) != PRODUCER_SHA256:
        raise RuntimeError("the reference producer does not match its pinned source")
    args.work_root.mkdir(parents=True, exist_ok=False)
    comparator_path = Path(__file__).with_name("cache_values_compare.py")
    if args.command == "check":
        spec = importlib.util.spec_from_file_location("compare_values", comparator_path)
        comparator = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(comparator)
    reports = {}
    saved = json.loads(args.reference_report.read_text()) if args.reference_report else None
    if saved:
        verify_pin(saved["producer"])
        if saved["producer"]["sha256"] != PRODUCER_SHA256:
            raise RuntimeError("saved oracle producer differs")
    reused = 0
    for name, case in cases().items():
        folder = args.work_root / name
        folder.mkdir()
        if saved and name in saved["cases"]:
            old = saved["cases"][name]
            if old["spec_sha256"] != hashlib.sha256(json.dumps(case, sort_keys=True).encode()).hexdigest():
                raise RuntimeError("saved fixture specification differs: " + name)
            vcf = verify_pin(old["vcf"])
            for record in old["sidecars"]:
                verify_pin(record)
            reference = verify_pin(old["python"]["output"])
            python_stderr = verify_pin(old["python"]["stderr"])
            verify_pin(old["python"]["stdout"])
            python_status = old["python"]["returncode"]
            reused += 1
        else:
            vcf = make_case(folder, case)
            reference = folder / "python.pkl"
            reference.write_bytes(SENTINEL)
            python_status = run([sys.executable, str(args.producer), str(vcf), str(reference), case.get("threads", "1")], folder, "python")
            python_stderr = folder / "python.stderr.log"
        python_written = reference.read_bytes() != SENTINEL
        python_category = failure_category(python_stderr) if python_status else None
        report = {"python_exit": python_status, "python_output_written": python_written,
                  "python_exception_category": python_category, "python_stderr": str(python_stderr),
                  "saved_reference_reused": bool(saved and name in saved["cases"])}
        if python_status == 0 and python_written:
            value = pickle.loads(reference.read_bytes())
            report.update(genes=len(value), none=sum(v is None for v in value.values()))
        if args.command == "check":
            archive = folder / "native.rkyv"
            archive.write_bytes(SENTINEL)
            status = run([str(args.binary), "clinvar-gene-af", "build", "--vcf", str(vcf), "--output", str(archive), "--threads", case.get("threads", "1")], folder, "native")
            written = archive.read_bytes() != SENTINEL
            native_category = failure_category(folder / "native.stderr.log") if status else None
            report.update(native_exit=status, native_output_written=written, native_exception_category=native_category)
            report["passed"] = matching_outcome(python_status, status, python_category, native_category) and written == python_written
            if status == 0 and python_status == 0 and written and python_written:
                from types import SimpleNamespace
                comparison = comparator.compare(SimpleNamespace(binary=args.binary, cache="clinvar_patho_af_stat", pickle=reference,
                    archive=archive, report=folder/"values.json", mismatches=folder/"mismatches.jsonl"))
                report["comparison"] = comparison
                report["passed"] &= comparison["passed"]
        if args.command == "check" and name in ["flag_empty_contig", "flag_header_order_tie", "flag_none_merge_failure", "normal", "csi", "header_order_tie", "cross_chr_larger", "merge_none_first", "merge_none_second", "merge_none_both", "none_separate_genes", "empty_header_contig", "duplicate_csq",
            "no_index", "ordinary_gzip", "plain_unindexed", "invalid_unused_gt", "invalid_unused_integer_format",
            "invalid_pos", "negative_pos", "empty_pos", "af_string", "af_number_g", "af_empty_tuple",
            "undefined_af", "undefined_af_absent", "csq_integer", "missing_alt_winner"]:
            contigs = case.get("contigs", ["chr1"])
            merge_command = [str(args.binary), "clinvar-gene-af", "merge"]
            shard_status = []
            shard_categories = []
            shard_preserved = []
            for index, chrom in enumerate(contigs):
                shard = folder / f"chr{index}.rkyv"
                shard.write_bytes(SENTINEL)
                shard_status.append(run([str(args.binary), "clinvar-gene-af", "chromosome", "--vcf", str(vcf), "--chrom", chrom, "--output", str(shard)],folder,f"chr{index}"))
                shard_categories.append(failure_category(folder / f"chr{index}.stderr.log") if shard_status[-1] else None)
                shard_preserved.append(not shard_status[-1] or shard.read_bytes() == SENTINEL)
                merge_command += ["--input",str(shard)]
            merged = folder/"merged.rkyv"
            merged.write_bytes(SENTINEL)
            if any(shard_status):
                first_failure = next(i for i, status in enumerate(shard_status) if status)
                merged_status = None
                pipeline_status = shard_status[first_failure]
                pipeline_category = shard_categories[first_failure]
            else:
                merged_status = run(merge_command+["--output",str(merged)],folder,"merge")
                pipeline_status = merged_status
                pipeline_category = failure_category(folder / "merge.stderr.log") if merged_status else None
            merged_written = merged.read_bytes() != SENTINEL
            pipeline = {"shard_exits":shard_status,"merge_exit":merged_status,"merge_output_written":merged_written,
                        "shard_categories":shard_categories,"failed_shard_outputs_preserved":all(shard_preserved),
                        "pipeline_exit":pipeline_status,"pipeline_exception_category":pipeline_category,
                        "passed":matching_outcome(python_status,pipeline_status,python_category,pipeline_category)
                            and merged_written == python_written and all(shard_preserved)}
            if pipeline_status == 0 and python_status == 0 and merged_written:
                comparison = comparator.compare(SimpleNamespace(binary=args.binary,cache="clinvar_patho_af_stat",pickle=reference,
                    archive=merged,report=folder/"merged.values.json",mismatches=folder/"merged.mismatches.jsonl"))
                pipeline["comparison"] = comparison
                pipeline["passed"] &= comparison["passed"]
            report["chromosome_pipeline"] = pipeline
            report["passed"] &= pipeline["passed"]
        reports[name] = report
    result = {"producer_sha256": digest(args.producer), "count": len(reports), "cases": reports,
              "reused_python_cases":reused,"exact_failure_checks":sum(row["python_exit"] != 0 for row in reports.values())}
    result["failed"] = [name for name,r in reports.items() if not r.get("passed", True)]
    (args.work_root / "comparison.json").write_text(json.dumps(result, sort_keys=True, indent=2)+"\n")
    print(json.dumps({"count":len(reports),"failed":result["failed"],"python_success":sum(r["python_exit"] == 0 for r in reports.values())}))
    if result["failed"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ["probe", "check"]:
        p = commands.add_parser(name)
        p.add_argument("--producer", type=Path, required=True)
        p.add_argument("--work-root", type=Path, required=True)
        p.add_argument("--reference-report", type=Path, help="Verified saved Python fixtures; preserve their producing source and diagnostics")
        if name == "check":
            p.add_argument("--binary", type=Path, required=True)
        p.set_defaults(function=check)
    p = commands.add_parser("admission")
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--work-root", type=Path, required=True)
    p.set_defaults(function=admission)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
