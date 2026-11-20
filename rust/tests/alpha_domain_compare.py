#!/usr/bin/env python3
"""Independent AlphaMissense domain producer/Parquet comparisons.

check runs bounded unchanged Python procedures and native commands; oracle is
the subprocess boundary that preserves their runtime exceptions. No DAS policy
or VEP is executed. All generated data/logs require an explicit work directory.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import importlib
import json
import math
from pathlib import Path
import pickle
import re
import os
import shutil
import struct
import subprocess
import sys

sys.dont_write_bytecode = True
FIELDS = ["Gene", "Feature_type", "Feature", "EXON", "Protein_position", "DOMAINS", "SYMBOL", "TAIL"]
DEFAULT = dict(zip(FIELDS, ["ENSG_ONE", "Transcript", "ENST_ONE", "2/4", "2/20", "Pfam:PF1", "GENE", "end"]))
PIN = "cc7b2d60665e8052e6362107e359cc133e31c398f243589ee63ad7e8f90341c4"
SENTINEL = b"previous complete publication\n"


def annotation(fields=FIELDS, **values):
    data = dict(DEFAULT, **values)
    return "|".join(data.get(name, "unknown") for name in fields)


def row(chrom="chr1", pos=10, score="0.25", af="0.1", variant="A2T", csq=None, extra=""):
    return f"{chrom}\t{pos}\t.\tA\tT\t.\tPASS\tAM_PATHOGENICITY={score};AF_grpmax_joint={af};PVAR={variant};CSQ={csq if csq is not None else annotation()}{extra}\n"


def cases():
    items = {
        "normal": {}, "empty": {"rows": []}, "csi": {"csi": True},
        "plain": {"compression": "plain", "index": False}, "plain_indexed_route": {"compression": "plain", "index": False, "threads": "2"},
        "gzip": {"compression": "gzip", "index": False}, "no_index_sequential": {"index": False}, "no_index_parallel": {"index": False, "threads": "2"},
        "two_threads": {"threads": "2"}, "zero_threads": {"threads": "0"}, "negative_threads": {"threads": "-1"},
        "space_threads": {"threads": " 2 "}, "unicode_threads": {"threads": "٢"},
        "underscore_threads": {"threads": "0_2"}, "large_threads": {"threads": "9"*100},
        "empty_chromosome": {"contigs": ["chr1", "chr2"]},
        "no_header_contigs_sequential": {"contigs": [], "index": False},
        "empty_no_header_contigs_parallel": {"contigs": [], "rows": [], "index": False, "threads": "2"},
        "header_order": {"contigs": ["chr2", "chr1"], "threads": "2", "rows": [row(), row(chrom="chr2", score="0.75", csq=annotation(EXON="1/4", Protein_position="9/20", SYMBOL="FIRST"))]},
        "nonadjacent_gene": {"contigs": ["chr3", "chr1", "chr2"], "threads": "3", "rows": [row(), row(chrom="chr2", score="0.75"), row(chrom="chr3", score="0.125", af="0.3")]},
        "parent_multiplicity": {"values": {"DOMAINS": "PANTHER:A&B:X:P&X:Y&PANTHER:A:B:C&PANTHER:A"}},
        "duplicates": {"rows": [row(csq=annotation()+","+annotation()), row(), row(pos=11)]},
        "first_score_max_af": {"rows": [row(), row(pos=11, score="0.75", af="0.2"), row(pos=12, score="0.125", af="0.05")]},
        "residue_min_max": {"rows": [row(), row(pos=11, variant="A2G", score="0.75", af="0.2"), row(pos=12, variant="A2C", score="0.125", af="0.05")]},
        "transcript_order": {"rows": [row(csq=annotation(Feature="ENST_B", EXON="4/9")+","+annotation(Feature="ENST_A", EXON="1/4")), row(pos=11, csq=annotation(Feature="ENST_B", EXON="2/9"))]},
        "symbol_first_length_max": {"rows": [row(csq=annotation(SYMBOL="first", Protein_position="3/20")), row(pos=11, csq=annotation(SYMBOL="second", Protein_position="2/50"))]},
        "skipped_sources": {"values": {"DOMAINS": "PDB-ENSP_mappings:A&AFDB-ENSP_mappings:B&MobiDB_lite:C&Pfam:PF1"}},
        "only_skipped_sources": {"values": {"DOMAINS": "PDB-ENSP_mappings:A"}},
        "empty_domain": {"values": {"DOMAINS": ""}}, "empty_gene": {"values": {"Gene": ""}},
        "empty_exon": {"values": {"EXON": ""}}, "empty_protein_position": {"values": {"Protein_position": ""}},
        "nontranscript": {"values": {"Feature_type": "RegulatoryFeature"}}, "non_enst": {"values": {"Feature": "NM_1"}},
        "hierarchy_empty_ids": {"values": {"DOMAINS": "Pfam::X&:A&DB:"}},
        "unicode": {"values": {"Gene": "ENSG_α", "Feature": "ENST_β", "DOMAINS": "DB:δ", "EXON": "ε", "SYMBOL": "ζ"}, "variant": "α2δ"},
        "no_format": {"description": "none"}, "format_no_space": {"description": "Format:"+"|".join(FIELDS)},
        "last_required_field": {"fields": FIELDS[:-1]},
        "last_other_format": {"extra_header": '##INFO=<ID=OTHER,Number=.,Type=String,Description="Format: OTHER|TAIL">\n'},
        "duplicate_header": {"extra_header": '##INFO=<ID=CSQ,Number=1,Type=Integer,Description="Format: OTHER|TAIL">\n'},
        "duplicate_field": {"fields": FIELDS[:-1]+["Gene", "DOMAINS", "TAIL"]},
        "short_csq": {"rows": [row(csq="short")]}, "csq_empty": {"rows": [row(csq="")]},
        "csq_none_after_prefix": {"rows": [row(csq=annotation()+",,"+annotation())]},
        "csq_scalar": {"csq_number": "1"}, "csq_integer": {"csq_type": "Integer"},
        "csq_g": {"csq_number": "G"}, "csq_flag": {"csq_type": "Flag", "csq_number": "0"},
        "csq_zero_integer": {"csq_number": "1", "csq_type": "Integer", "rows": [row(csq="0")]},
        "csq_zero_float": {"csq_number": "1", "csq_type": "Float", "rows": [row(csq="-0")]},
        "score_tuple": {"score_number": "."}, "score_flag": {"score_type": "Flag", "score_number": "0"},
        "af_scalar_float": {"af_number": "1"}, "af_scalar_string": {"af_number": "1", "af_type": "String"},
        "af_tuple_string": {"af_type": "String"}, "af_g": {"af_number": "G"},
        "variant_tuple_layout": {"variant_number": "."},
        "variant_integer": {"variant_type": "Integer"}, "variant_empty": {"variant": ""},
        "invalid_unused_gt": {"sample": "GT", "sample_value": "oops"},
        "invalid_unused_dp": {"sample": "DP", "sample_value": "oops", "format_type": "Integer"},
        "valid_gt": {"sample": "GT", "sample_value": "0/1"},
        "negative_pos": {"mutate": ("chr1\t10\t", "chr1\t-1\t")},
        "bad_pos": {"mutate": ("chr1\t10\t", "chr1\tbad\t")},
        "empty_pos": {"mutate": ("chr1\t10\t", "chr1\t\t")},
        "missing_file": {"missing_file": True}, "no_eof": {"damage": "no_eof"},
        "crc": {"damage": "crc"}, "truncated": {"damage": "truncated"}, "trailing": {"damage": "trailing"},
        "crlf": {"crlf": True}, "output_directory": {"output_directory": True},
        "output_missing_parent": {"output_missing_parent": True},
        "literal_paths": {"filename": "input ;literal $(unused).vcf.gz"},
    }
    for name in ["AM_PATHOGENICITY", "AF_grpmax_joint", "PVAR", "CSQ"]:
        items["missing_"+name] = {"rows": [re.sub(r"(?:^|(?<=\t)|;)"+name+r"=[^;\n]*", "", row())]}
        items["undefined_"+name] = {"omit_headers": [name]}
    for field in FIELDS[:-1]:
        items["header_missing_"+field] = {"fields": [v for v in FIELDS if v != field]}
    for value in ["2-3/20", "-2/20", "bad/20", "2/bad", "2/0", "2/-9", "2", " 3 / 40", "+4/50", "1_2/30", "٣/٤٠"]:
        items["position_"+str(len(items))] = {"values": {"Protein_position": value}}
    for prefix in ["score", "af"]:
        for kind in ["Float", "Integer", "String"]:
            for value in ["0", "-0", ".", "NaN", "-NaN", "Inf", "-Inf", "1e-50", "1e40", "x", "0.123456789", "0.1junk", " 0.2", "0x1p-2", "1_2"]:
                items[f"{prefix}_{kind}_{len(items)}"] = {prefix+"_type": kind, prefix: value}
    for cut in [*range(33), 40, 64, 80, 100, 128, 160, 192, 224, 256]:
        items["truncated_at_"+str(cut)] = {"damage": cut}
    items.update(tuple_cases())
    items["tuple_nan_parent_multiplicity"] = dict(variant_type="Float", variant_number=".",
        rows=[row(variant="NaN,1", csq=annotation(DOMAINS="DB:A:B:C&DB:A")), row(pos=11,variant="NaN,1",csq=annotation(DOMAINS="DB:A:B:C&DB:A"))])
    items["tuple_nan_chromosomes"] = dict(variant_type="Float", variant_number=".", threads="2", contigs=["chr2","chr1"],
        rows=[row(variant="NaN,1"), row(chrom="chr2",variant="NaN,1",score="0.75")])
    return items


def make_case(folder, spec):
    import pysam
    fields = spec.get("fields", FIELDS)
    lines = ["##fileformat=VCFv4.2\n"]+[f"##contig=<ID={c},length=1000000>\n" for c in spec.get("contigs", ["chr1"])]
    for name, prefix, number, kind in [("AM_PATHOGENICITY", "score", "1", "Float"), ("AF_grpmax_joint", "af", ".", "Float"),
                                       ("PVAR", "variant", "1", "String"), ("CSQ", "csq", ".", "String")]:
        if name in spec.get("omit_headers", []): continue
        description = spec.get("description", "Format: "+"|".join(fields)) if name == "CSQ" else name
        lines.append(f'##INFO=<ID={name},Number={spec.get(prefix+"_number", number)},Type={spec.get(prefix+"_type", kind)},Description="{description}">\n')
    lines.append(spec.get("extra_header", ""))
    if spec.get("sample"):
        lines.append(f'##FORMAT=<ID={spec["sample"]},Number=1,Type={spec.get("format_type", "String")},Description="sample">\n')
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"+("\tFORMAT\tSAMPLE" if spec.get("sample") else "")+"\n")
    rows = spec.get("rows", [row(score=spec.get("score", "0.25"), af=spec.get("af", "0.1"), variant=spec.get("variant", "A2T"),
                                csq=annotation(fields, **spec.get("values", {})))])
    for value in rows:
        if spec.get("sample"): value = value.rstrip("\n")+"\t"+spec["sample"]+"\t"+spec["sample_value"]+"\n"
        lines.append(value)
    data = "".join(lines)
    if spec.get("crlf"): data = data.replace("\n", "\r\n")
    plain = folder / "input.vcf"
    plain.write_text(data, newline="")
    if spec.get("compression") == "plain": return plain
    path = folder / spec.get("filename", "input.vcf.gz")
    if spec.get("compression") == "gzip":
        with gzip.open(path, "wb") as out: out.write(data.encode())
        return path
    pysam.tabix_compress(str(plain), str(path), force=True)
    if spec.get("index", True): pysam.tabix_index(str(path), preset="vcf", force=True, csi=spec.get("csi", False))
    if spec.get("mutate"):
        plain.write_text(data.replace(*spec["mutate"]), newline=""); pysam.tabix_compress(str(plain), str(path), force=True)
    damage = spec.get("damage")
    if damage is not None:
        payload = bytearray(path.read_bytes())
        if isinstance(damage, int): payload = payload[:damage]
        elif damage == "no_eof": payload = payload[:-28]
        elif damage == "trailing": payload.extend(b"junk")
        elif damage == "truncated": payload = payload[:14]
        elif damage == "crc": payload[-36] ^= 1
        path.write_bytes(payload)
    if spec.get("missing_file"): return folder / "missing.vcf.gz"
    return path


def oracle(args):
    sys.path.insert(0, str(args.source.parent))
    source = importlib.import_module(args.source.stem)
    if args.procedure == "statistics":
        if args.raw:
            collector = source.DomainAMScoreCollector(str(args.vcf)); collector.collect_scores()
            with args.output.open("wb") as out: pickle.dump(source.to_plain_dict(collector.scores_dict), out)
        else:
            source.main(str(args.vcf), str(args.output), int(args.threads))
    else:
        targets = None if args.target_domain is None else set(args.target_domain)
        if args.procedure == "associated":
            wrapper = importlib.import_module("domain_metadata_cache")
            calls = 0
            def extract():
                nonlocal calls
                calls += 1
                return source.extract_domain_metadata_parallel(str(args.vcf), target_domains=targets, threads=int(args.threads))
            wrapper.load_or_extract_domain_metadata(cache_path=args.output, vcf_path=args.vcf, assembly=args.assembly,
                extractor=extract, source_paths=[args.source])
            print(json.dumps({"extractor_calls":calls}))
        else:
            value = source.extract_domain_metadata_parallel(str(args.vcf), target_domains=targets, threads=int(args.threads))
            with args.output.open("wb") as out: pickle.dump(value, out)


def tuple_cases():
    variants = {}
    for kind, values in [("String", ["A2T", "A2T,A2G", ".", ".,A2T", "A2T,.", "A2T,,A2G", ""]),
                         ("Integer", ["1", "1,2", ".", "1,.", "-1,0"]),
                         ("Float", ["1", "1,2", ".", "1,.", "-0,0", "NaN,1", "1,NaN", "NaN", "Inf,-Inf"]),
                         ("Flag", ["1", "."])]:
        for value in values:
            variants[f"tuple_{kind}_{len(variants)}"] = dict(variant_number=".", variant_type=kind,
                rows=[row(variant=value), row(pos=11, variant=value, score="0.75", af="0.2")])
    return variants


def keys_probe(args):
    if hashlib.sha256(args.producer.read_bytes()).hexdigest() != PIN: raise RuntimeError("original producer pin differs")
    args.work_root.mkdir(parents=True, exist_ok=False)
    records = []
    def key_shape(key):
        if isinstance(key, tuple): return {"type": "tuple", "items": [key_shape(item) for item in key]}
        if isinstance(key, float): return {"type": "float", "bits": struct.pack(">d", key).hex()}
        return {"type": type(key).__name__, "value": key}
    for name, spec in tuple_cases().items():
        folder = args.work_root/name; folder.mkdir()
        vcf = make_case(folder, spec)
        output = folder/"raw.pkl"
        code = run([sys.executable, __file__, "oracle", "--procedure", "statistics", "--source", args.producer,
                    "--vcf", vcf, "--output", output, "--raw"], folder, "python")
        record = dict(case=name, python_exit=code, exception=category(folder/"python.stderr.log"))
        if code == 0:
            with output.open("rb") as stream: value = pickle.load(stream)
            domain = value["ENSG_ONE"]["Pfam"]["PF1"]
            record["keys"] = {name:[key_shape(k) for k in domain[name]] for name in
                              ["distribution", "gnomAD_distribution", "gnomAD_residue_distribution", "min_distribution", "max_distribution"]}
            record["sample_lengths"] = list(map(len, domain["gnomAD_residue_distribution"].values()))
        records.append(record); print(json.dumps(record), flush=True)
    (args.work_root/"summary.json").write_text(json.dumps({"cases":records,"count":len(records)},indent=2)+"\n")


class SampleRuns:
    def __init__(self, values, counts): self.values,self.counts=values,counts


def scientific_key(value):
    if value["kind"] == "text": return value["text"]
    if value["kind"] != "tuple": raise ValueError("unknown scientific key tag")
    return tuple(None if item["kind"] == "none" else item[{"text":"text", "int":"integer", "float":"float", "bool":"boolean"}[item["kind"]]] for item in value["elements"])


def iter_genes(path, expand_samples=True):
    import numpy as np
    import pyarrow.parquet as pq
    file = pq.ParquetFile(path)
    current = None
    for batch in file.iter_batches(batch_size=4096):
        for row_ in batch.to_pylist():
            g, kind, key = row_["gene"], row_["kind"], row_["key"]
            if kind == "metadata":
                yield g,row_["gene_order"],dict(gene_symbol=key, aa_start=row_["minimum_int"], aa_end=row_["maximum"], prot_length=row_["auxiliary_int"])
                continue
            if kind == "gene":
                if current is not None: yield current
                current = (g,row_["gene_order"],{});continue
            if current is None or current[:2] != (g,row_["gene_order"]): raise ValueError("noncontiguous model gene")
            node = current[2]
            path_ = row_["path"][:-1] if kind in ["source", "domain"] else row_["path"]
            for part in path_: node = node[part]
            if kind in ["source", "domain"]: node[key] = {}
            elif kind == "transcript":
                node[key] = dict(exons=set(row_["keys"]), aa_min=row_["minimum_int"] if row_["minimum_int"] is not None else row_["minimum_float"], aa_max=row_["maximum"])
            elif kind == "scores": node[key] = dict(zip(map(scientific_key, row_["typed_keys"]), row_["values"]))
            elif kind == "array": node[key] = np.array(row_["values"], dtype=np.float64)
            elif kind == "residues": node[key] = {}
            elif kind == "residue":
                items = SampleRuns(row_["values"],row_["counts"])
                if expand_samples:
                    items = []
                    for value, count in zip(row_["values"], row_["counts"]): items.extend([value]*count)
                node["gnomAD_residue_distribution"][scientific_key(row_["typed_keys"][0])] = items
            else: raise ValueError("Unknown native row kind: "+kind)
    if current is not None: yield current


def decode(path):
    values = sorted(iter_genes(path),key=lambda item:item[1])
    return {gene:value for gene,_,value in values}


def exact(expected, actual, path="root"):
    import numpy as np
    if type(expected) is not type(actual): raise AssertionError(f"{path}: type {type(expected).__name__} != {type(actual).__name__}")
    if isinstance(expected, dict):
        if len(expected) != len(actual): raise AssertionError(f"{path}: dictionary length differs")
        for index, ((left_key,left), (right_key,right)) in enumerate(zip(expected.items(), actual.items())):
            exact(left_key, right_key, path+f"/key-{index}")
            exact(left, right, path+"/"+str(left_key))
    elif isinstance(expected, np.ndarray):
        if expected.dtype != actual.dtype or expected.shape != actual.shape or expected.tobytes() != actual.tobytes(): raise AssertionError(f"{path}: exact array bits/dtype/shape differ")
    elif isinstance(expected, (tuple, list)):
        if len(expected) != len(actual): raise AssertionError(f"{path}: multiplicity differs")
        for i, (left, right) in enumerate(zip(expected, actual)): exact(left, right, path+f"/{i}")
    elif isinstance(expected, float):
        if struct.pack(">d", expected) != struct.pack(">d", actual): raise AssertionError(f"{path}: float bits differ")
    elif expected != actual: raise AssertionError(f"{path}: scalar/set differs")


def run(command, folder, name):
    with (folder/(name+".stdout.log")).open("x") as out, (folder/(name+".stderr.log")).open("x") as err:
        return subprocess.run(list(map(str, command)), stdout=out, stderr=err, timeout=90).returncode


def category(path, native=False):
    text = path.read_text(errors="replace")
    if native:
        match = re.search(r"priva: ([A-Za-z_][A-Za-z_0-9]*):", text)
        return match[1] if match else None
    matches = re.findall(r"^([A-Za-z_][A-Za-z_0-9]*):", text, re.M)
    return matches[-1] if matches else None


def cache_checks(args):
    import pyarrow as pa
    import pyarrow.parquet as pq
    args.work_root.mkdir(parents=True, exist_ok=False)
    results = []
    def stamp(path):
        s = path.stat(); return s.st_ino,s.st_mtime_ns,s.st_size,s.st_mode & 0o777
    def pair(folder, vcf, pyout, rsout, name, assembly="hg38", threads="1", targets=()):
        target_args = [part for target in targets for part in ["--target-domain",target]]
        py = run([sys.executable,__file__,"oracle","--procedure","associated","--source",args.producer.with_name("am_pick_intolerant_domains.py"),
                  "--vcf",vcf,"--output",pyout,"--assembly",assembly,"--threads",threads,*target_args],folder,name+".python")
        rs = run([args.binary,"domain-metadata","build","--vcf",vcf,"--output",rsout,"--assembly",assembly,"--threads",threads,*target_args],folder,name+".native")
        record = dict(case=name,python_exit=py,native_exit=rs,
            python_exception=category(folder/(name+".python.stderr.log")),native_exception=category(folder/(name+".native.stderr.log"),True))
        try:
            if py != rs: raise AssertionError("exit codes differ")
            if py:
                if record["python_exception"] != record["native_exception"]: raise AssertionError("exception categories differ")
            else:
                with pyout.open("rb") as stream: expected = pickle.load(stream)
                exact(expected,decode(rsout))
                messages = [json.loads(v) for v in (folder/(name+".python.stdout.log")).read_text().splitlines() if v.startswith('{')]
                record["python_extractor_calls"] = messages[-1]["extractor_calls"]
                messages = [json.loads(v) for v in (folder/(name+".native.stdout.log")).read_text().splitlines() if v.startswith('{')]
                record["native_reused"] = messages[-1].get("reused",False)
                if record["native_reused"] != (record["python_extractor_calls"] == 0): raise AssertionError("reuse decisions differ")
            record["matched"] = True
        except Exception as error:
            record.update(matched=False,mismatch=str(error))
        return record
    for scenario in ["cold_missing_parent","warm","alias_assembly","changed_assembly","index_removed","changed_vcf_metadata",
                     "changed_targets_warm","empty_regeneration","zero_length_reuse","changed_output","missing_receipt","invalid_source",
                     "existing_mode","receipt_directory","merged_cache_reuse"]:
        folder = args.work_root/scenario; folder.mkdir()
        spec = {"rows":[]} if scenario == "empty_regeneration" else ({"values":{"Protein_position":"2"}} if scenario == "zero_length_reuse" else {})
        vcf = make_case(folder,spec)
        pyout,rsout = folder/"cache/reference.pkl",folder/"cache/native.parquet"
        if scenario != "cold_missing_parent": pyout.parent.mkdir()
        if scenario == "existing_mode":
            pyout.write_bytes(SENTINEL);rsout.write_bytes(SENTINEL)
            pyout.chmod(0o660);rsout.chmod(0o660)
        first = pair(folder,vcf,pyout,rsout,"setup")
        if not first["matched"]:
            first["case"] = scenario+"/setup";results.append(first);continue
        before = [stamp(p) for p in [pyout,rsout]]
        for path in [pyout,rsout,Path(str(pyout)+".install.json"),Path(str(rsout)+".source.json")]:
            if path.is_file(): shutil.copy2(path,folder/("before."+path.name))
        assembly = "GRCh38" if scenario == "alias_assembly" else "hg19" if scenario == "changed_assembly" else "hg38"
        threads = "1"; targets = ()
        if scenario == "index_removed": Path(str(vcf)+".tbi").unlink(); threads="2"
        if scenario == "changed_vcf_metadata":
            old=vcf.stat();os.utime(vcf,ns=(old.st_atime_ns,old.st_mtime_ns+1000000))
        if scenario == "changed_targets_warm": targets=("not_a_matching_domain",)
        if scenario == "changed_output": pyout.write_bytes(b"bad pickle");rsout.write_bytes(b"bad parquet")
        if scenario == "missing_receipt": Path(str(pyout)+".install.json").unlink();Path(str(rsout)+".source.json").unlink()
        if scenario == "receipt_directory":
            for path in [Path(str(pyout)+".install.json"),Path(str(rsout)+".source.json")]:path.unlink();path.mkdir()
        if scenario == "invalid_source": vcf.write_bytes(b"\x1f")
        if scenario == "merged_cache_reuse":
            shard=folder/"metadata.raw.parquet"
            if run([args.binary,"domain-metadata","chromosome","--vcf",vcf,"--chrom","chr1","--output",shard,"--assembly","GRCh38"],folder,"shard"):
                raise RuntimeError("metadata shard setup failed")
            if run([args.binary,"domain-metadata","merge","--input",shard,"--output",rsout,"--assembly","hg38"],folder,"merge"):
                raise RuntimeError("metadata association merge failed")
            before=[stamp(p) for p in [pyout,rsout]]
            Path(str(vcf)+".tbi").unlink();threads="2"
        record = pair(folder,vcf,pyout,rsout,"operation",assembly,threads,targets);record["case"] = scenario
        if record["matched"] and scenario in ["warm","alias_assembly","index_removed","changed_targets_warm","zero_length_reuse","invalid_source","merged_cache_reuse"]:
            if before != [stamp(p) for p in [pyout,rsout]]:record.update(matched=False,mismatch="cache changed on reuse/failure")
        if record["matched"] and scenario == "existing_mode" and any(stamp(p)[3] != 0o660 for p in [pyout,rsout]):
            record.update(matched=False,mismatch="existing output permission mode changed")
        results.append(record);print(json.dumps(record),flush=True)
    # Native-format integrity and ordered-shard admission are separate from
    # source Python exceptions. Each invalid artifact retains its own evidence.
    folder=args.work_root/"native_guards";folder.mkdir()
    vcf=make_case(folder,{"contigs":["chr1","chr2"]}); shards=[]
    for i,chrom in enumerate(["chr1","chr2"]):
        shard=folder/f"chr{i}.parquet";shards.append(shard)
        if run([args.binary,"alphamissense-domain","chromosome","--vcf",vcf,"--chrom",chrom,"--output",shard],folder,f"create{i}"):
            raise RuntimeError("native guard setup failed")
    for name,inputs in [("missing_shard",shards[:1]),("reordered_shards",list(reversed(shards))),("duplicate_shards",[shards[0],shards[0]])]:
        output=folder/(name+".parquet");output.write_bytes(SENTINEL)
        code=run([args.binary,"alphamissense-domain","merge",*[v for p in inputs for v in ["--input",p]],"--output",output],folder,name)
        record=dict(case=name,native_exit=code,matched=code==1 and output.read_bytes()==SENTINEL)
        results.append(record);print(json.dumps(record),flush=True)
    source=pq.read_table(shards[0]); rows=source.to_pylist()
    for name in ["zero_sample_count","unknown_key_tag","duplicate_gene","wrong_empty_schema"]:
        altered=folder/(name+".parquet")
        values=json.loads(json.dumps(rows))
        if name=="zero_sample_count":next(row for row in values if row["kind"]=="residue")["counts"][0]=0
        elif name=="unknown_key_tag":next(row for row in values if row["kind"]=="residue")["typed_keys"][0]["kind"]="invalid"
        elif name=="duplicate_gene":values.extend(values)
        if name=="wrong_empty_schema": table=pa.table({"unexpected":pa.array([],type=pa.int64())}).replace_schema_metadata(source.schema.metadata)
        else: table=pa.Table.from_pylist(values,schema=source.schema)
        pq.write_table(table,altered,compression="zstd")
        code=run([args.binary,"alphamissense-domain","check","--input",altered],folder,name)
        record=dict(case=name,native_exit=code,matched=code==1)
        results.append(record);print(json.dumps(record),flush=True)
    summary=dict(cases=len(results),results=results,mismatches=[v for v in results if not v["matched"]],passed=all(v["matched"] for v in results))
    (args.work_root/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    if not summary["passed"]:raise SystemExit(1)


def check(args):
    if hashlib.sha256(args.producer.read_bytes()).hexdigest() != PIN: raise RuntimeError("original producer pin differs")
    args.work_root.mkdir(parents=True, exist_ok=False)
    results, mismatches, boundaries, pipelines = [], [], [], []
    for procedure in ["statistics", "metadata"]:
        source = args.producer if procedure == "statistics" else args.producer.with_name("am_pick_intolerant_domains.py")
        for name, spec in cases().items():
            folder = args.work_root/(procedure+"_"+name); folder.mkdir()
            vcf = make_case(folder, spec)
            reference, native = folder/"reference.pkl", folder/"native.parquet"
            reference.write_bytes(SENTINEL); native.write_bytes(SENTINEL)
            if spec.get("output_directory"):
                reference.unlink(); native.unlink(); reference.mkdir(); native.mkdir()
            if spec.get("output_missing_parent"):
                reference, native = folder/"absent/reference.pkl", folder/"absent/native.parquet"
            command = "alphamissense-domain" if procedure == "statistics" else "domain-metadata"
            threads = spec.get("threads", "1")
            py = run([sys.executable, __file__, "oracle", "--procedure", procedure, "--source", source, "--vcf", vcf, "--output", reference, "--threads", threads], folder, "python")
            rs = run([args.binary, command, "build", "--vcf", vcf, "--output", native, "--threads", threads], folder, "native")
            record = dict(procedure=procedure, case=name, python_exit=py, native_exit=rs,
                          python_exception=category(folder/"python.stderr.log"), native_exception=category(folder/"native.stderr.log", True))
            problem = None
            if spec.get("unsupported") and procedure == "statistics" and py == 0 and record["native_exception"] == "UnsupportedLayout":
                record["qualification"] = "unqualified_successful_python_layout"; boundaries.append(record)
            else:
                if py != rs: problem = "exit codes differ"
                elif py:
                    if record["python_exception"] != record["native_exception"]: problem = "exception categories differ"
                    if native.is_file() and native.read_bytes() != SENTINEL: problem = "previous output changed after processing failure"
                else:
                    try:
                        with reference.open("rb") as stream: expected = pickle.load(stream)
                        exact(expected, decode(native))
                    except Exception as error: problem = str(error)
                record["matched"] = problem is None
                if problem: record["mismatch"] = problem; mismatches.append(record)
            results.append(record)
            print(json.dumps(record), flush=True)
            if record.get("matched") and py == 0 and (name.startswith("tuple_") or name in ["normal", "header_order", "nonadjacent_gene", "empty_chromosome", "parent_multiplicity", "duplicates", "first_score_max_af", "residue_min_max", "transcript_order", "symbol_first_length_max", "only_skipped_sources"]):
                shards = []
                for ordinal, chrom in enumerate(spec.get("contigs", ["chr1"])):
                    shard = folder/f"chromosome-{ordinal}.parquet"
                    status = run([args.binary, command, "chromosome", "--vcf", vcf, "--chrom", chrom, "--output", shard], folder, f"chromosome-{ordinal}")
                    if status: raise RuntimeError("chromosome failed: "+procedure+"/"+name)
                    shards += ["--input", shard]
                merged = folder/"merged.parquet"
                status = run([args.binary, command, "merge", *shards, "--output", merged], folder, "merge")
                pipeline = dict(procedure=procedure, case=name, merge_exit=status)
                try:
                    if status: raise AssertionError("ordered merge failed")
                    with reference.open("rb") as stream: expected = pickle.load(stream)
                    exact(expected, decode(merged))
                    pipeline["matched"] = True
                except Exception as error:
                    pipeline["matched"] = False; pipeline["mismatch"] = str(error); mismatches.append(pipeline)
                pipelines.append(pipeline)
                print(json.dumps(dict(pipeline=pipeline)), flush=True)
    summary = dict(cases=len(results), qualified_cases=len(results)-len(boundaries), matched=sum(r.get("matched", False) for r in results),
                   mismatches=mismatches, unqualified_layouts=boundaries, results=results, pipelines=pipelines, passed=not mismatches)
    (args.work_root/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    if mismatches: raise SystemExit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    test = commands.add_parser("check")
    test.add_argument("--producer", type=Path, required=True); test.add_argument("--binary", type=Path, required=True)
    test.add_argument("--work-root", type=Path, required=True)
    original = commands.add_parser("oracle")
    original.add_argument("--procedure", choices=["statistics", "metadata", "associated"], required=True)
    original.add_argument("--source", type=Path, required=True); original.add_argument("--vcf", type=Path, required=True)
    original.add_argument("--output", type=Path, required=True); original.add_argument("--threads", default="1")
    original.add_argument("--raw", action="store_true")
    original.add_argument("--assembly", default="hg38")
    original.add_argument("--target-domain", action="append")
    probe = commands.add_parser("keys_probe")
    probe.add_argument("--producer", type=Path, required=True); probe.add_argument("--work-root", type=Path, required=True)
    cached = commands.add_parser("cache_checks")
    cached.add_argument("--producer", type=Path, required=True); cached.add_argument("--binary", type=Path, required=True)
    cached.add_argument("--work-root", type=Path, required=True)
    args = parser.parse_args(); globals()[args.command](args)


if __name__ == "__main__": main()
