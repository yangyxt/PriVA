#!/usr/bin/env python3
"""Differential producer tests with explicit worker-completion traces.

  python clinvar_gene_stat_edges.py check --producer ORIGINAL.py --binary PRIVA \
      --work-root NEW_DIRECTORY

Reference runs call the pinned Python producer. For controlled-order cases only,
the Pool transport yields unchanged worker results in a declared permutation;
the original worker, whole-result collector and statistics are unmodified.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

from clinvar_gene_stat_compare import FIELDS, digest, producer, invoke


def annotation(index=1, **changes):
    values = dict(zip(FIELDS, ["SYMBOL", "ENSG_1", "protein_coding", f"ENST_{index}",
        "Pfam:PF1&InterPro:IPR1", "missense_variant", "", "", "NA", "0.0", "0.0"]))
    values.update(changes)
    return "|".join(values[name] for name in FIELDS)


def record(index=1, *, csq=None, **changes):
    value = dict(chrom="chr1", pos=str(index), ref="A", alt="T", allele=str(index),
                 review="reviewed_by_expert_panel", significance="Pathogenic",
                 csq=csq if csq is not None else annotation(index), extra="")
    value.update(changes)
    return value


def cases():
    """Small VCFs cover source values and exceptions, including unused data."""
    yield "one_record", {}
    yield "empty_records", {"records": []}
    yield "no_ensembl", {"records": [record(csq=annotation(Feature="NM_1"))]}
    for name, changes in [
        ("review_low", {"review": "criteria_provided,_single_submitter"}),
        ("review_two_star", {"review": "criteria_provided,_multiple_submitters,_no_conflicts"}),
        ("review_four_star", {"review": "practice_guideline"}),
        ("review_missing", {"review": "."}),
        ("pathogenic_case", {"significance": "not_pathogenic"}),
        ("benign", {"significance": "Benign"}),
        ("first_alt", {"alt": "T,G"}),
        ("empty_ref", {"ref": ""}), ("empty_alt", {"alt": ""}),
        ("dot_alt", {"alt": "."}), ("empty_pos", {"pos": ""}),
        ("bad_pos", {"pos": "bad"}), ("zero_pos", {"pos": "0"}),
        ("negative_pos", {"pos": "-1"}),
        ("missing_alleleid", {"allele": None}), ("null_alleleid", {"allele": "."}),
        ("alleleid_tuple", {"allele": "1,2"}),
        ("numeric_info_bad", {"extra": ";AF_joint=bad"}),
        ("numeric_info_tuple", {"extra": ";AF_joint=0.3,0.1"}),
        ("duplicate_info_first", {"extra": ";ALLELEID=2"}),
        ("unknown_info", {"extra": ";UNKNOWN=some_text"}),
    ]:
        yield name, {"records": [record(**changes)]}
    for name, changes in [
        ("noncoding", {"BIOTYPE": "lncRNA"}),
        ("gene_missing", {"Gene": ""}),
        ("domain_missing", {"DOMAINS": ""}),
        ("domain_numeric", {"DOMAINS": "1"}),
        ("domain_nan_substring", {"DOMAINS": "domain_nan_fragment&OTHER"}),
        ("domain_empty_component", {"DOMAINS": "&PF&&PF&"}),
        ("consequence_missing", {"Consequence": ""}),
        ("consequence_numeric", {"Consequence": "1"}),
        ("synonymous", {"Consequence": "synonymous_variant"}),
        ("missense_synonymous", {"Consequence": "missense_variant&synonymous_variant"}),
        ("lof_hc", {"LoF": "HC"}), ("lof_os", {"LoF": "OS"}),
        ("lof_normalization", {"LoF": "lc& hc &os"}),
        ("nmd_escape", {"LoF": "HC", "NMD": "NMD_escaping_variant"}),
        ("nmd_case", {"LoF": "HC", "NMD": "ESCAPE"}),
        ("donor_below", {"SpliceAI_pred_DS_DL": "0.249999999"}),
        ("donor_equal", {"SpliceAI_pred_DS_DL": "0.25"}),
        ("donor_above", {"SpliceAI_pred_DS_DL": "0.250001"}),
        ("donor_float32_boundary", {"SpliceAI_pred_DS_DL": "0.250000001"}),
        ("donor_missing", {"SpliceAI_pred_DS_DL": ""}),
        ("donor_bad", {"SpliceAI_pred_DS_DL": "bad"}),
        ("splice_missing_delta", {"SpliceVault_top_events": "Top1:CD:+4:40%:Inframe", "SpliceVault_SpliceAI_delta": ""}),
        ("splice_below_ignored_bad", {"SpliceVault_top_events": "bad", "SpliceVault_SpliceAI_delta": "0.1"}),
        ("splice_bad", {"SpliceVault_top_events": "bad", "SpliceVault_SpliceAI_delta": "0.5"}),
        ("splice_bad_frequency", {"SpliceVault_top_events": "Top1:ES:1:0.4:Inframe", "SpliceVault_SpliceAI_delta": "0.5"}),
        ("splice_bad_offset", {"SpliceVault_top_events": "Top1:CD:bad:40%:Inframe", "SpliceVault_SpliceAI_delta": "0.5"}),
        ("splice_delta_bad", {"SpliceVault_SpliceAI_delta": "bad"}),
    ]:
        yield name, {"records": [record(csq=annotation(**changes))]}
    for length in [0, 3, 9, 10, 14, 15, 16]:
        yield f"inframe_{length}", {"records": [record(alt="T"+"A"*length,
            csq=annotation(Consequence="inframe_insertion"))]}
    for offset in ["2", "3", "9", "10", "15", "16", "+-+4", "٤", "1_6", "9"*30]:
        yield "offset_"+hashlib.sha256(offset.encode()).hexdigest()[:8], {"records": [record(csq=annotation(
            SpliceVault_top_events=f"Top1:CD:{offset}:40%:Inframe", SpliceVault_SpliceAI_delta="0.5"))]}
    for event in ["Top1:ES:2-3:20%:Frameshift", "Top1:ES:2-3:#40:Inframe",
                  "Top1:CD:+4:20%%:Inframe", "Top1:CA:-20:40%:Inframe",
                  "Top1:OTHER:bad:40%:Inframe", "Top1:CD:bad:40%:Frameshift",
                  "Top1:ES:2-3:10%:Inframe&Top2:CD:4:40%:Inframe"]:
        yield "event_"+hashlib.sha256(event.encode()).hexdigest()[:8], {"records": [record(csq=annotation(
            SpliceVault_top_events=event, SpliceVault_SpliceAI_delta="0.5"))]}
    yield "duplicate_csq_last", {"records": [record(csq=annotation(SYMBOL="FIRST")+","+annotation(SYMBOL="LAST", LoF="HC"))]}
    yield "short_csq", {"records": [record(csq="SYMBOL|ENSG_1|protein_coding|ENST_1")]}
    yield "long_csq", {"records": [record(csq=annotation()+"|extra")]}
    yield "domain_mixed", {"records": [record(), record(2, csq=annotation(2, DOMAINS="1"))]}
    yield "gene_numeric", {"records": [record(csq=annotation(Gene="1"))]}
    yield "gene_numeric_noncoding", {"records": [record(csq=annotation(Gene="1", BIOTYPE="lncRNA"))]}
    yield "gene_numeric_lowreview", {"records": [record(csq=annotation(Gene="1"), review="none")]}
    yield "empty_ref_mixed", {"records": [record(ref=""), record(2, csq=annotation(2))]}
    yield "empty_alt_mixed", {"records": [record(alt=""), record(2, csq=annotation(2))]}
    yield "duplicate_variant_last_allele", {"records": [record(), record(2, pos="1", csq=annotation(2))]}
    yield "global_cross_chrom_gene", {"records": [record(), record(2, chrom="chr2", csq=annotation(2))]}
    for order in [[0, 1], [1, 0]]:
        yield "prefilter_symbol_"+"".join(map(str, order)), {"records": [
            record(csq=annotation(SYMBOL="FIRST", BIOTYPE="lncRNA")),
            record(2, csq=annotation(2, SYMBOL="SECOND"))], "order": order}
        # One worker returns 20,000 rows. The threshold applies to the whole
        # result, not an arbitrary split within it. Only ENST_1 is the target.
        yield "chunk_coercion_"+"".join(map(str, order)), {"records": [
            record(csq=",".join(annotation(index, SpliceAI_pred_DS_DL="0.250000001") for index in range(1, 20001))),
            record(2, csq=annotation(2, BIOTYPE="lncRNA", SpliceAI_pred_DS_DL="bad"))], "order": order}
    ordered = [record(index, csq=annotation(index, DOMAINS=f"D{index}&COMMON&D{index%3}")) for index in range(1, 25)]
    for seed in [0, 1, 42, 4294967295]:
        yield f"nested_order_seed_{seed}", {"records": ordered, "seed": seed}
    for name in ["normal_gzip", "missing_eof", "plain", "crlf", "blank_record", "bad_unused_format", "valid_sample",
                 "missing_required_info", "duplicate_required_info", "clncsq_bad", "clncsq_good", "duplicate_csq_field",
                 "summary_duplicate_last", "mane_duplicate_last", "summary_missing_name", "summary_numeric_name",
                 "summary_missing_id", "summary_float_id", "summary_two_nan_ids", "summary_quoted", "summary_grch37",
                 "mane_numeric_target", "mane_mixed_target", "mane_missing_refseq",
                 "summary_two_nan_object_ids", "mane_two_nan_refseq", "mane_two_nan_object_refseq"]:
        yield name, {"mutation": name}


def write_case(folder, spec):
    import pysam
    from clinvar_annotation_contract import TABLE_INFO_FIELDS
    folder.mkdir()
    header = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1000000>"]
    for name in sorted(TABLE_INFO_FIELDS):
        header.append(f'##INFO=<ID={name},Number=.,Type=Float,Description="Fixture frequency">')
    for name in ["CLNSIG", "CLNREVSTAT", "CLNDN"]:
        header.append(f'##INFO=<ID={name},Number=.,Type=String,Description="Fixture annotation">')
    header += ['##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="Fixture identifier">',
               '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: '+"|".join(FIELDS)+'">',
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    records = spec.get("records", [record()])
    lines = []
    for row in records:
        info = ("" if row["allele"] is None else "ALLELEID="+row["allele"]+";")
        info += "CLNSIG="+row["significance"]+";CLNREVSTAT="+row["review"]+";CSQ="+row["csq"]+row["extra"]
        lines.append("\t".join([row["chrom"], row["pos"], ".", row["ref"], row["alt"], ".", "PASS", info]))
    count = max(2, len(records))
    summary = ["#AlleleID\tAssembly\tName"]+[f"{i}\tGRCh38\tNM_{i:06d}.1:c.1A>T" for i in range(1,count+1)]
    mane = ["RefSeq_nuc\tEnsembl_nuc"]+[f"NM_{i:06d}.1\tENST_{i}.1" for i in range(1,count+1)]
    mutation = spec.get("mutation")
    if mutation in ["bad_unused_format", "valid_sample"]:
        header.insert(-1, '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">')
        header[-1] += "\tFORMAT\tSAMPLE"
        lines = [line+"\tDP\t"+("bad" if mutation.startswith("bad") else "12") for line in lines]
    if mutation == "blank_record": lines.append("")
    if mutation == "missing_required_info": header = [x for x in header if not x.startswith("##INFO=<ID=AF_joint,")]
    if mutation == "duplicate_required_info": header.insert(-1, '##INFO=<ID=AF_joint,Number=1,Type=String,Description="Second">')
    if mutation in ["clncsq_bad", "clncsq_good"]:
        description = "not a format" if mutation.endswith("bad") else "Format: SYMBOL|Consequence"
        header.insert(-1, f'##INFO=<ID=CLNCSQ,Number=.,Type=String,Description="{description}">')
    if mutation == "duplicate_csq_field": header = [x.replace("Format: SYMBOL|", "Format: SYMBOL|SYMBOL|") for x in header]
    if mutation == "summary_duplicate_last": summary.append("1\tGRCh38\tNM_000002.1:c.1A>T")
    if mutation == "mane_duplicate_last": mane.append("NM_000001.1\tENST_2.1")
    if mutation == "summary_missing_name": summary = [summary[0], "1\tGRCh38\t", summary[2]]
    if mutation == "summary_numeric_name": summary = [summary[0], "1\tGRCh38\t123", "2\tGRCh38\t456"]
    if mutation == "summary_missing_id": summary[1] = summary[1][1:]
    if mutation == "summary_float_id": summary[1] = "1.0"+summary[1][1:]
    if mutation == "summary_two_nan_ids": summary = [summary[0], summary[1][1:], summary[2][1:]]
    if mutation == "summary_two_nan_object_ids": summary = [summary[0], summary[1][1:], summary[2][1:], "text\tGRCh38\tNM_000001.1:c.1A>T"]
    if mutation == "summary_quoted": summary[1] = '1\tGRCh38\t"NM_000001.1:c.1A>T\tquoted"'
    if mutation == "summary_grch37": summary = [x.replace("GRCh38", "GRCh37") for x in summary]
    if mutation == "mane_numeric_target": mane = [mane[0], "NM_000001.1\t1", "NM_000002.1\t2"]
    if mutation == "mane_mixed_target": mane[2] = "NM_000002.1\t2"
    if mutation == "mane_missing_refseq": mane[1] = "\tENST_1.1"
    if mutation == "mane_two_nan_refseq": mane = [mane[0], "\tENST_1.1", "\tENST_2.1"]
    if mutation == "mane_two_nan_object_refseq": mane += ["\tENST_1.1", "\tENST_2.1"]
    separator = "\r\n" if mutation == "crlf" else "\n"
    payload = (separator.join(header+lines)+separator).encode()
    plain = folder/"input.vcf"
    plain.write_bytes(payload)
    vcf = folder/"input.vcf.gz"
    if mutation == "plain": vcf = plain
    elif mutation == "normal_gzip": vcf.write_bytes(gzip.compress(payload))
    else:
        pysam.tabix_compress(str(plain), str(vcf), force=False)
        if mutation == "missing_eof": vcf.write_bytes(vcf.read_bytes()[:-28])
    (folder/"summary.tsv").write_text("\n".join(summary)+"\n")
    (folder/"mane.tsv").write_text("\n".join(mane)+"\n")
    if "order" in spec:
        (folder/"order.json").write_text(json.dumps({"vcf_sha256": digest(vcf), "order": spec["order"]})+"\n")
    return vcf


class ReplayPool:
    """Controlled transport only; source worker and collector stay unchanged."""
    order = None
    def __init__(self, threads): pass
    def imap_unordered(self, function, args):
        results = [function(arg) for arg in args]
        if sorted(self.order) != list(range(len(results))):
            raise ValueError("fixture completion order is not a permutation")
        yield from (results[index] for index in self.order)
    def close(self): pass
    def join(self): pass
    def terminate(self): pass


_WORKER = None
_POOL = None


def indexed_worker(item):
    index, args = item
    return index, _WORKER(args)


class ObservedPool:
    order = []
    def __init__(self, *args, **kwargs): self.pool = _POOL(*args, **kwargs)
    def imap_unordered(self, function, args):
        for index, result in self.pool.imap_unordered(indexed_worker, enumerate(args)):
            self.order.append(index)
            yield result
    def __getattr__(self, name): return getattr(self.pool, name)


def reference(args):
    original = producer(args.producer)
    # producer() imports the source module from its pinned script directory.
    import combine_annotations
    mapping = None
    if args.target_map:
        import pickle
        if not args.target_map_sha256 or digest(args.target_map) != args.target_map_sha256:
            raise ValueError("saved transcript-map identity differs")
        def prepared_map(*unused_args, **unused_keywords):
            # This is the explicitly supplied trusted result of the source's
            # summary/MANE procedure, not a guessed or rebuilt mapping.
            with args.target_map.open("rb") as stream:
                return pickle.load(stream)
        original.offer_clinvar_alleleid_to_transcript_map = prepared_map
        mapping = {"path": str(args.target_map), "sha256": args.target_map_sha256}
    if args.order:
        ReplayPool.order = json.loads(args.order.read_text())["order"]
        combine_annotations.mp.Pool = ReplayPool
    elif args.observe:
        global _WORKER, _POOL
        _WORKER, _POOL = combine_annotations.convert_record_to_tab, combine_annotations.mp.Pool
        ObservedPool.order = []
        combine_annotations.mp.Pool = ObservedPool
    report = {"success": False, "prepared_target_map": mapping}
    try:
        value = original.main_stat_variants(str(args.vcf), str(args.output), args.threads if args.order or args.observe else 1,
            clinvar_tab=str(args.summary), mane_tranx_tab=str(args.mane), donor_loss_cutoff=.25)
        report.update(success=True, genes=len(value),
            set_counts={key: sum(len(gene[key]) for gene in value.values()) for key in [
                "clinvar_variants", "putative_nmd_variants", "small_aachange", "large_aachange", "no_aachange"]})
    except Exception as error:
        report.update(exception=type(error).__name__, message=str(error))
    if args.observe:
        args.observe.write_text(json.dumps({"vcf_sha256": digest(args.vcf), "order": ObservedPool.order})+"\n")
    args.report.write_text(json.dumps(report, indent=2)+"\n")


def check(args):
    import cache_values_compare as compare
    producer(args.producer)  # add pinned imports for fixture construction only
    args.work_root.mkdir(parents=True, exist_ok=False)
    reports = {}
    for name, spec in cases():
        if args.case and name not in args.case: continue
        folder = args.work_root/name
        inputs = args.reference_fixtures/name if args.reference_fixtures else folder
        if args.reference_fixtures:
            folder.mkdir()
            vcf = inputs/"input.vcf.gz"
            if not vcf.is_file(): vcf = inputs/"input.vcf"
        else:
            vcf = write_case(folder, spec)
        reference_path = inputs/"reference.pkl.gz"
        command = [sys.executable, Path(__file__), "reference", "--producer", args.producer,
            "--vcf", vcf, "--summary", inputs/"summary.tsv", "--mane", inputs/"mane.tsv",
            "--output", reference_path, "--report", inputs/"reference.json"]
        if "order" in spec: command += ["--order", inputs/"order.json"]
        seed = str(spec.get("seed", 0))
        if not args.reference_fixtures:
            with (folder/"python.log").open("x") as log:
                code = subprocess.run(list(map(str, command)), stdout=log, stderr=subprocess.STDOUT,
                    env=dict(os.environ, PYTHONHASHSEED=seed, PYTHONDONTWRITEBYTECODE="1")).returncode
            if code:
                raise RuntimeError(f"reference harness failed for {name}, exit {code}")
        reference_result = json.loads((inputs/"reference.json").read_text())
        archive = folder/"native.rkyv"
        sentinel = b"previous generation\n"
        archive.write_bytes(sentinel)
        command = [args.binary, "clinvar-gene-stats", "build", "--vcf", vcf,
            "--clinvar-summary", inputs/"summary.tsv", "--mane-summary", inputs/"mane.tsv",
            "--output", archive, "--threads", "2" if "order" in spec else "1",
            "--python-hash-seed", seed, "--donor-loss-cutoff", ".25"]
        if "order" in spec: command += ["--completion-order", inputs/"order.json"]
        status = invoke(command, folder, "native")
        row = {"python": reference_result, "native_exit": status,
               "passed": reference_result["success"] == (status == 0), "seed": int(seed)}
        if status:
            row["output_preserved"] = archive.read_bytes() == sentinel
            row["passed"] &= row["output_preserved"]
            row["native_error"] = (folder/"native.stderr.log").read_text().splitlines()[-2:]
        if reference_result["success"] and not status:
            try:
                row["comparison"] = compare.compare(SimpleNamespace(binary=args.binary,
                    cache="clinvar_gene_stat", pickle=reference_path, archive=archive,
                    report=folder/"values.json", mismatches=folder/"mismatches.jsonl"))
                row["passed"] &= row["comparison"]["passed"]
            except Exception as error:
                row["comparison_error"] = str(error)
                row["passed"] = False
        reports[name] = row
        print(json.dumps({"case": name, "passed": row["passed"], "python_success": reference_result["success"], "native_exit": status}), flush=True)
    result = {"passed": all(row["passed"] for row in reports.values()), "cases": reports,
              "producer_sha256": digest(args.producer), "binary_sha256": digest(args.binary),
              "scope": "complete producer fixtures; explicit completion traces only where declared; no universal unobserved-trace equivalence"}
    (args.work_root/"comparison.json").write_text(json.dumps(result, indent=2)+"\n")
    if not result["passed"]: raise SystemExit(1)


def transport(args):
    """Actual pool arrivals in both languages, plus publication/trace controls."""
    import cache_values_compare as compare
    producer(args.producer)
    args.work_root.mkdir(parents=True, exist_ok=False)
    reports = {}
    records = [record(index, csq=",".join(annotation(index+j*100, DOMAINS=f"D{index}&COMMON")
        for j in range(1+index%5))) for index in range(1, 31)]
    for side in ["python_arrival", "native_arrival"]:
        folder = args.work_root/side
        vcf = write_case(folder, {"records": records})
        trace = folder/"order.json"
        archive = folder/"native.rkyv"
        native = [args.binary, "clinvar-gene-stats", "build", "--vcf", vcf,
            "--clinvar-summary", folder/"summary.tsv", "--mane-summary", folder/"mane.tsv",
            "--output", archive, "--threads", "3", "--python-hash-seed", "42", "--donor-loss-cutoff", ".25"]
        reference_command = [sys.executable, Path(__file__), "reference", "--producer", args.producer,
            "--vcf", vcf, "--summary", folder/"summary.tsv", "--mane", folder/"mane.tsv",
            "--output", folder/"reference.pkl.gz", "--report", folder/"reference.json", "--threads", "3"]
        if side == "native_arrival":
            status = invoke(native+["--completion-trace", trace], folder, "native")
            if status: raise RuntimeError("native default transport failed")
            reference_command += ["--order", trace]
        else:
            reference_command += ["--observe", trace]
        with (folder/"python.log").open("x") as log:
            code = subprocess.run(list(map(str, reference_command)), stdout=log, stderr=subprocess.STDOUT,
                env=dict(os.environ, PYTHONHASHSEED="42", PYTHONDONTWRITEBYTECODE="1")).returncode
        reference_result = json.loads((folder/"reference.json").read_text())
        if code or not reference_result["success"]: raise RuntimeError("reference transport failed")
        if side == "python_arrival":
            status = invoke(native+["--completion-order", trace], folder, "native")
        order = json.loads(trace.read_text())["order"]
        row = {"native_exit": status, "records": len(records), "observed_source_order": order == list(range(len(records))),
               "trace_is_permutation": sorted(order) == list(range(len(records))), "trace_sha256": digest(trace)}
        if not status:
            row["comparison"] = compare.compare(SimpleNamespace(binary=args.binary, cache="clinvar_gene_stat",
                pickle=folder/"reference.pkl.gz", archive=archive, report=folder/"values.json", mismatches=folder/"mismatches.jsonl"))
        row["passed"] = status == 0 and row["trace_is_permutation"] and row.get("comparison", {}).get("passed", False)
        reports[side] = row
    folder = args.work_root/"controls"
    vcf = write_case(folder, {})
    sentinel = b"preserved generation\n"
    archive = folder/"native.rkyv"
    base = [args.binary, "clinvar-gene-stats", "build", "--vcf", vcf,
        "--clinvar-summary", folder/"summary.tsv", "--mane-summary", folder/"mane.tsv",
        "--output", archive, "--threads", "1", "--python-hash-seed", "0"]
    for name, trace in [("duplicate", {"vcf_sha256": digest(vcf), "order": [0, 0]}),
                        ("short", {"vcf_sha256": digest(vcf), "order": []}),
                        ("long", {"vcf_sha256": digest(vcf), "order": [0, 1]}),
                        ("negative", {"vcf_sha256": digest(vcf), "order": [-1]}),
                        ("float", {"vcf_sha256": digest(vcf), "order": [0.0]}),
                        ("wrong_input", {"vcf_sha256": "0"*64, "order": [0]})]:
        path = folder/(name+".json")
        path.write_text(json.dumps(trace)+"\n")
        archive.write_bytes(sentinel)
        status = invoke(base+["--completion-order", path], folder, name)
        reports[name] = {"native_exit": status, "passed": status != 0 and archive.read_bytes() == sentinel}
    for name, destination in [("input_alias", vcf), ("trace_alias", archive)]:
        archive.write_bytes(sentinel)
        original_hash = digest(vcf)
        status = invoke(base+["--completion-trace", destination], folder, name)
        reports[name] = {"native_exit": status, "passed": status != 0 and digest(vcf) == original_hash and archive.read_bytes() == sentinel}
    result = {"passed": all(row["passed"] for row in reports.values()), "cases": reports,
              "scope": "observed default worker arrivals replayed without changing source worker/collector; explicit trace and publication rejection controls"}
    (args.work_root/"comparison.json").write_text(json.dumps(result, indent=2)+"\n")
    if not result["passed"]: raise SystemExit(1)


def merge_check(args):
    import pickle
    import cache_values_compare as compare
    args.work_root.mkdir(parents=True, exist_ok=False)
    reports = {}
    for name in ["text", "mixed_numeric_text", "all_missing"]:
        folder = args.work_root/name
        folder.mkdir()
        reference_path = args.source_probes/name/"reference.pkl.gz"
        with gzip.open(reference_path, "rb") as stream: reference_value = pickle.load(stream)
        archives = []
        for index, (key, value) in enumerate(reference_value.items()):
            source = folder/f"part_{index}.pkl"
            source.write_bytes(pickle.dumps({key: value}, protocol=5))
            archive = folder/f"part_{index}.rkyv"
            status = invoke([args.binary, "cache-convert", "clinvar_gene_stat", "--pickle", source, "--output", archive], folder, f"convert_{index}")
            if status: raise RuntimeError("shard conversion failed")
            archives.append(archive)
        merged = folder/"merged.rkyv"
        command = [args.binary, "clinvar-gene-stats", "merge", "--output", merged]
        for path in archives: command += ["--input", path]
        status = invoke(command, folder, "merge")
        row = {"native_exit": status}
        if not status:
            row["comparison"] = compare.compare(SimpleNamespace(binary=args.binary, cache="clinvar_gene_stat",
                pickle=reference_path, archive=merged, report=folder/"values.json", mismatches=folder/"mismatches.jsonl"))
        row["passed"] = status == 0 and row.get("comparison", {}).get("passed", False)
        reports[name] = row
        before = merged.read_bytes()
        status = invoke(command+["--input", archives[0]], folder, "duplicate")
        reports[name+"_overlap_refused"] = {"native_exit": status, "passed": status != 0 and merged.read_bytes() == before}
        before_hash = digest(archives[0])
        status = invoke([args.binary, "clinvar-gene-stats", "merge", "--input", archives[0], "--output", archives[0]], folder, "alias")
        reports[name+"_alias_refused"] = {"native_exit": status, "passed": status != 0 and digest(archives[0]) == before_hash}
        bad = folder/"malformed.rkyv"
        bad.write_bytes(b"not an archive")
        status = invoke(command+["--input", bad], folder, "malformed")
        reports[name+"_malformed_refused"] = {"native_exit": status, "passed": status != 0 and merged.read_bytes() == before}
    result = {"passed": all(row["passed"] for row in reports.values()), "cases": reports,
              "scope": "disjoint version-1/version-2 partition union only; overlap refusal does not restrict the whole-input builder"}
    (args.work_root/"comparison.json").write_text(json.dumps(result, indent=2)+"\n")
    if not result["passed"]: raise SystemExit(1)


def operational(args):
    """Retain actual source-only estimator failures as explicit divergences."""
    import cache_values_compare as compare
    producer(args.producer)
    args.work_root.mkdir(parents=True, exist_ok=False)
    reports = {}
    for name in ["input_path_space", "estimator_path_empty"]:
        folder = args.work_root/name
        folder.mkdir()
        inputs = folder/("input with space" if name == "input_path_space" else "inputs")
        vcf = write_case(inputs, {})
        reference = folder/"reference.pkl.gz"
        command = [sys.executable, Path(__file__), "reference", "--producer", args.producer,
            "--vcf", vcf, "--summary", inputs/"summary.tsv", "--mane", inputs/"mane.tsv",
            "--output", reference, "--report", folder/"reference.json"]
        environment = dict(os.environ, PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1")
        if name == "estimator_path_empty": environment["PATH"] = ""
        with (folder/"python.log").open("x") as log:
            code = subprocess.run(list(map(str, command)), stdout=log, stderr=subprocess.STDOUT, env=environment).returncode
        if code: raise RuntimeError("operational reference harness failed")
        result = json.loads((folder/"reference.json").read_text())
        archive = folder/"native.rkyv"
        status = invoke([args.binary, "clinvar-gene-stats", "build", "--vcf", vcf,
            "--clinvar-summary", inputs/"summary.tsv", "--mane-summary", inputs/"mane.tsv",
            "--output", archive, "--threads", "1", "--python-hash-seed", "0", "--donor-loss-cutoff", ".25"], folder, "native")
        row = {"python": result, "reference_output_exists": reference.exists(), "native_exit": status,
               "native_output_exists": archive.exists(), "source_failure_equivalent": result["success"] == (status == 0),
               "role": "explicit operational divergence, excluded from producer failure-equivalence counts"}
        if not status:
            row["comparison_to_same_valid_data"] = compare.compare(SimpleNamespace(binary=args.binary, cache="clinvar_gene_stat",
                pickle=args.reference_fixtures/"one_record/reference.pkl.gz", archive=archive,
                report=folder/"values.json", mismatches=folder/"mismatches.jsonl"))
        row["passed"] = not result["success"] and not status and row.get("comparison_to_same_valid_data", {}).get("passed", False)
        reports[name] = row
    result = {"passed": all(row["passed"] for row in reports.values()), "cases": reports,
              "meaning_of_passed": "the declared mismatch was reproduced and native values equal the same valid-data reference; this is not failure equivalence",
              "scientific_or_caller_waiver": False}
    (args.work_root/"comparison.json").write_text(json.dumps(result, indent=2)+"\n")
    if not result["passed"]: raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    sub = commands.add_parser("reference")
    for name in ["producer", "vcf", "summary", "mane", "output", "report"]:
        sub.add_argument("--"+name, type=Path, required=True)
    sub.add_argument("--order", type=Path)
    sub.add_argument("--observe", type=Path)
    sub.add_argument("--threads", type=int, default=2)
    sub.add_argument("--target-map", type=Path)
    sub.add_argument("--target-map-sha256")
    sub.set_defaults(function=reference)
    sub = commands.add_parser("check")
    for name in ["producer", "binary", "work-root"]:
        sub.add_argument("--"+name, type=Path, required=True)
    sub.add_argument("--case", action="append")
    sub.add_argument("--reference-fixtures", type=Path)
    sub.set_defaults(function=check)
    sub = commands.add_parser("transport")
    for name in ["producer", "binary", "work-root"]:
        sub.add_argument("--"+name, type=Path, required=True)
    sub.set_defaults(function=transport)
    sub = commands.add_parser("merge")
    for name in ["producer", "binary", "source-probes", "work-root"]:
        sub.add_argument("--"+name, type=Path, required=True)
    sub.set_defaults(function=merge_check)
    sub = commands.add_parser("operational")
    for name in ["producer", "binary", "reference-fixtures", "work-root"]:
        sub.add_argument("--"+name, type=Path, required=True)
    sub.set_defaults(function=operational)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
