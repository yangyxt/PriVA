#!/usr/bin/env python3
"""Compare versioned splice conversion and the unchanged Python splice consumer.

  python clinvar_splice_typed_compare.py check --scripts SCRIPTS --binary PRIVA \
      --cases-root COMPLETED_PRODUCER_CASES --work-root NEW_DIRECTORY

Reuses completed source-produced synthetic caches. Every cache value is checked;
the consumer comparisons cover complete result frames and conditional failures.
"""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import pickle
import struct
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import acmg_step_oracle_common as oracle
import cache_values_compare as compare
from clinvar_aa_splice_compare import run, SENTINEL

PINNED = {
    "acmg_ps1_pm5_known_variant.py": "f3b76ac982c8f0e1956528e844b521baae7da2433f4b780c61ff20f7998a6a48",
    "splicing_var_analysis.py": "d4ae0fec7e477ad020ca5157f4af4be2660ac737b8ce7f6946687b855b561196",
    "acmg_pvs1_null_variant.py": "3f77d412cfb388f53a93f05d768b6b08ba478a6dc503dee85bab98dc4a069ae9",
}


def convert(binary, source, output, folder, name):
    return run([binary, "cache-convert", "clinvar_splice_stat", "--pickle", source, "--output", output], folder, name)


def consumer(args, folder, value, aa_archive, splice_archive, module):
    records = [record for records in value.values() for record in records]
    # Use the first actual source transcript, and add absent/matched scope tests.
    transcript = next(iter(value), "ABSENT")
    original = {"Feature": transcript, "INTRON": "", "EXON": "1/2", "SpliceAI_pred_DS_AG": .9,
                "SpliceAI_pred_DS_AL": .2, "SpliceAI_pred_DS_DG": .1, "SpliceAI_pred_DS_DL": .1,
                "HGVSc": "ENST_ONE:c.2+1A>T", "STRAND": 1, "pos": 10, "ref": "A", "alt": "T", "chrom": "chr1"}
    if records and records[0]["exon"] is None and records[0]["intron"]:
        original["INTRON"] = records[0]["intron"]
    variants = {"matching": {}, "nonmatching": {"EXON": "other", "INTRON": "other"},
                "below_score": {name: .01 for name in original if name.startswith("SpliceAI_pred_DS_")},
                "missing_query_exon": {"EXON": None}, "missing_query_intron": {"INTRON": None},
                "nan_query_score": {"SpliceAI_pred_DS_AG": np.nan},
                "other_largest_score": {"SpliceAI_pred_DS_AL": .99},
                "absent_transcript": {"Feature": "ABSENT"}, "noncanonical_query": {"HGVSc": "ENST_ONE:c.2A>T"}}
    results = {}
    for name, update in variants.items():
        case = folder/name
        case.mkdir()
        rows = [dict(original, **update)] * 3
        frame = pd.DataFrame(rows)
        for field in ["Feature", "INTRON", "EXON", "HGVSc", "ref", "alt", "chrom"]:
            frame[field] = pd.Series([row[field] for row in rows], dtype=object)
        strengths = [0, 1, 4]
        oracle.write_frame(frame, case/"input.frame")
        oracle.capture(case, lambda: [module.check_splice_pathogenic(row, value.get(row["Feature"], []), strength, .5)
                                      for row, strength in zip(frame.to_dict("records"), strengths)],
                       lambda result: pd.DataFrame({"result": pd.Series([str(v) for v in result], dtype=object)}))
        (case/"observed.frame").write_bytes(SENTINEL)
        status, error = oracle.run_rust(args.binary, "acmg_ps1_pm5_known_variant", "check_splice_pathogenic", case,
            {"clinvar_aa_dict_pkl": str(aa_archive), "clinvar_splice_dict_pkl": str(splice_archive), "pvs1_criteria": strengths})
        (case/"native.stderr.log").write_text(error)
        expected_ok = (case/"expected.frame").is_file()
        record = {"python_success": expected_ok, "native_exit": status, "passed": expected_ok == (status == 0)}
        if expected_ok and status == 0:
            record["differences"] = oracle.compare_frames(case/"expected.frame", case/"observed.frame", show=1000000)
            record["passed"] &= not record["differences"]
        elif not expected_ok:
            record["python_error"] = json.loads((case/"expected_error.json").read_text())
            record["preserved"] = (case/"observed.frame").read_bytes() == SENTINEL
            record["passed"] &= record["preserved"] and record["python_error"]["type"] in error
        results[name] = record
    # Exercise the callers of the shared reader. The empty AA/PM5 inputs
    # isolate the changed splice path; PM5's independently configured residue
    # rules are outside this comparison and receive no candidate records.
    for name in ["ps1_splice_route", "pvs1_gene_summary"]:
        case = folder/name
        case.mkdir()
        options = {"clinvar_aa_dict_pkl": str(args.work_root/"empty.aa.rkyv"),
                   "clinvar_splice_dict_pkl": str(splice_archive)}
        if name == "ps1_splice_route":
            frame = pd.DataFrame([dict(original, Protein_position="", HGVSp="", Consequence="intron_variant",
                CLNSIG="Benign" if i < 3 else "likely_pathogenic", CLNREVSTAT="reviewed_by_expert_panel") for i in range(6)])
            strengths = [0,1,4,0,1,4]
            def expected():
                return module.PS1_PM5_criteria(frame.copy(deep=True), {}, value, np.zeros(6,dtype=bool),
                    np.array(strengths,dtype=np.int64), threads=1)
            def to_frame(result):
                return pd.DataFrame({"PS1":result[0],"PM5":result[1],"pm5_rule":pd.Series([""]*6,dtype=object)})
            rust_module, function = "acmg_ps1_pm5_known_variant", "PS1_PM5_criteria"
            options.update(clinvar_pm5_index=str(args.work_root/"empty.pm5.tsv"),pvs1_criteria=strengths,threads=1,gof_hgvsp=[])
        else:
            import acmg_pvs1_null_variant as pvs1
            keys = list(value)+["ABSENT"]
            frame = pd.DataFrame({"Feature":pd.Series(keys,dtype=object),"Gene":pd.Series(["G"+str(i) for i in range(len(keys))],dtype=object)})
            mapping = dict(zip(frame["Feature"],frame["Gene"]))
            expected = lambda: pvs1.summarize_clinvar_gene_pathogenicity(mapping,{},value)
            to_frame = lambda result: pd.DataFrame({"gene":pd.Series(sorted(result),dtype=object)})
            rust_module, function = "acmg_pvs1_null_variant", "summarize_clinvar_gene_pathogenicity"
        oracle.write_frame(frame,case/"input.frame")
        oracle.capture(case,expected,to_frame)
        (case/"observed.frame").write_bytes(SENTINEL)
        status,error = oracle.run_rust(args.binary,rust_module,function,case,options)
        (case/"native.stderr.log").write_text(error)
        py_ok = (case/"expected.frame").is_file()
        record = {"python_success":py_ok,"native_exit":status,"passed":py_ok == (status==0)}
        if py_ok and status==0:
            record["differences"] = oracle.compare_frames(case/"expected.frame",case/"observed.frame",show=1000000)
            record["passed"] &= not record["differences"]
        elif not py_ok:
            record["python_error"] = json.loads((case/"expected_error.json").read_text())
            record["preserved"] = (case/"observed.frame").read_bytes() == SENTINEL
            record["passed"] &= record["preserved"] and record["python_error"]["type"] in error
        results[name] = record
    return results


def check(args):
    args.work_root.mkdir(parents=True, exist_ok=False)
    for name, wanted in PINNED.items():
        if oracle.sha256_file(args.scripts/name) != wanted:
            raise RuntimeError("Python consumer identity differs: " + name)
    sys.path.insert(0, str(args.scripts))
    import acmg_ps1_pm5_known_variant as module
    (args.work_root/"empty.aa.pkl").write_bytes(pickle.dumps({},protocol=4))
    if run([args.binary,"cache-convert","clinvar_aa_stat","--pickle",args.work_root/"empty.aa.pkl",
            "--output",args.work_root/"empty.aa.rkyv"],args.work_root,"empty_aa"):
        raise RuntimeError("empty AA control conversion failed")
    (args.work_root/"empty.pm5.tsv").write_text("\t".join(["transcript","protein_position","hgvsp","consequence","nmd","symbol",
        "hgvsc","chrom","pos","ref","alt","clnsig","clnrevstat"])+"\n")
    result = {"provenance": oracle.manifest([args.scripts/name for name in PINNED]), "cases": {}}
    for name in ["normal", "all_splice_empty", "partial_splice", "partial_splice_diff_rows", "extra_splice_field",
                 "trailing_splice_header", "alt_missing", "missing_consequence", "missing_hgvsc", "missing_intron", "other_last_format",
                 "missing_exon_present_intron"]:
        folder = args.work_root/name
        folder.mkdir()
        source = args.cases_root/name/"python.splice.pkl"
        aa_archive = args.cases_root/name/"native.aa.rkyv"
        archive = folder/"converted.rkyv"
        status = convert(args.binary, source, archive, folder, "convert")
        record = {"conversion_exit": status, "passed": status == 0}
        if not status:
            record["values"] = compare.compare(SimpleNamespace(binary=args.binary, cache="clinvar_splice_stat", pickle=source,
                archive=archive, report=folder/"values.json", mismatches=folder/"mismatches.jsonl"))
            record["schema_version"] = struct.unpack("<I", archive.read_bytes()[16:20])[0]
            value = pickle.loads(source.read_bytes())
            record["consumers"] = consumer(args, folder, value, aa_archive, archive, module)
            record["passed"] &= record["values"]["passed"] and all(v["passed"] for v in record["consumers"].values())
        result["cases"][name] = record
    # An invalid score is source-supported text, and must fail only when read.
    base = pickle.loads((args.cases_root/"normal/python.splice.pkl").read_bytes())
    for name, changes in [("invalid_matched_score", {"SpliceAI_pred_DS_AG": "bad"}),
                          ("missing_selected_score", {"SpliceAI_pred_DS_AG": None}),
                          ("nullable_mixed_records", {})]:
        value = copy.deepcopy(base)
        record = next(iter(value.values()))[0]
        for key, text in changes.items():
            if text is None:
                record["splice_ai"].pop(key)
            else:
                record["splice_ai"][key] = text
        if name == "nullable_mixed_records":
            record["alt"] = None
            record["exon"], record["intron"] = "other", "other"
            next(iter(value.values())).append(copy.deepcopy(next(iter(base.values()))[0]))
        folder = args.work_root/name
        folder.mkdir()
        source, archive = folder/"input.pkl", folder/"converted.rkyv"
        source.write_bytes(pickle.dumps(value, protocol=4))
        status = convert(args.binary, source, archive, folder, "convert")
        r = {"conversion_exit": status, "passed": status == 0}
        if not status:
            r["values"] = compare.compare(SimpleNamespace(binary=args.binary, cache="clinvar_splice_stat", pickle=source,
                archive=archive, report=folder/"values.json", mismatches=folder/"mismatches.jsonl"))
            r["consumers"] = consumer(args, folder, value, args.cases_root/"normal/native.aa.rkyv", archive, module)
            r["passed"] &= r["values"]["passed"] and all(v["passed"] for v in r["consumers"].values())
        result["cases"][name] = r
    result["passed"] = all(v["passed"] for v in result["cases"].values())
    (args.work_root/"comparison.json").write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"passed": result["passed"], "cases": len(result["cases"]),
                      "failed": [k for k,v in result["cases"].items() if not v["passed"]]}))
    if not result["passed"]:
        raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    p = parser.add_subparsers(dest="command", required=True).add_parser("check")
    for name in ["scripts", "binary", "cases-root", "work-root"]:
        p.add_argument("--"+name, type=Path, required=True)
    p.set_defaults(function=check)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
