#!/usr/bin/env python3
"""Compare native domain Parquet against already completed Python checkpoints.

compare accepts a JSON list of {path, sha256, genes} reference entries. It
deserializes one original chromosome and one native gene at a time. Optional
--finalize-source applies only the pinned original dictionary-to-array finalizer;
it never collects variants, runs VEP, or executes the DAS procedure.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import pickle
import resource
import struct
import sys
import time

sys.dont_write_bytecode = True
from alpha_domain_compare import SampleRuns, iter_genes, PIN


def digest(path):
    with Path(path).open("rb") as stream: return hashlib.file_digest(stream,"sha256").hexdigest()


def exact(expected, actual, counts, path="root"):
    import numpy as np
    if isinstance(actual,SampleRuns):
        if type(expected) is not list or any(type(v) is not float for v in expected):
            raise AssertionError(path+": residue sample Python types differ")
        if len(actual.values) != len(actual.counts) or any(type(n) is not int or n <= 0 for n in actual.counts) or sum(actual.counts) != len(expected):
            raise AssertionError(path+": sample multiplicity differs")
        bits=np.asarray(expected,dtype=np.float64).view(np.uint64);at=0
        for value,count in zip(actual.values,actual.counts):
            word=struct.unpack("=Q",struct.pack("=d",value))[0]
            if np.any(bits[at:at+count] != word): raise AssertionError(f"{path}: sample bits differ in run starting {at}")
            at+=count
        counts["sample_lists"]+=1;counts["samples"]+=len(expected);return
    if type(expected) is not type(actual):raise AssertionError(f"{path}: type {type(expected).__name__} != {type(actual).__name__}")
    counts["nodes"]+=1
    if isinstance(expected,dict):
        counts["dicts"]+=1
        if len(expected) != len(actual):raise AssertionError(path+": dictionary length differs")
        for i,((ek,ev),(ak,av)) in enumerate(zip(expected.items(),actual.items())):
            exact(ek,ak,counts,f"{path}/key-{i}");exact(ev,av,counts,path+"/"+str(ek))
    elif isinstance(expected,np.ndarray):
        counts["arrays"]+=1;counts["array_values"]+=expected.size
        if expected.dtype != actual.dtype or expected.shape != actual.shape or not np.array_equal(expected.view(np.uint64),actual.view(np.uint64)):
            raise AssertionError(path+": array dtype/shape/bits differ")
    elif isinstance(expected,(tuple,list)):
        if len(expected) != len(actual):raise AssertionError(path+": sequence length differs")
        for i,(ev,av) in enumerate(zip(expected,actual)):exact(ev,av,counts,f"{path}/{i}")
    elif isinstance(expected,set):
        counts["sets"]+=1;counts["set_members"]+=len(expected)
        if expected != actual:raise AssertionError(path+": set differs")
    elif isinstance(expected,float):
        counts["float_scalars"]+=1
        if struct.pack("=d",expected) != struct.pack("=d",actual):raise AssertionError(path+": float bits differ")
    elif expected != actual:raise AssertionError(path+": scalar differs")


def compare(args):
    import pyarrow.parquet as pq
    entries=json.loads(args.reference_list.read_text())
    producer=None
    if args.finalize_source:
        if digest(args.finalize_source) != PIN:raise RuntimeError("original finalizer source differs")
        spec=importlib.util.spec_from_file_location("domain_reference_finalizer",args.finalize_source)
        producer=importlib.util.module_from_spec(spec);spec.loader.exec_module(producer)
    header=pq.ParquetFile(args.native).metadata.metadata
    if header[b"priva.kind"] != b"statistics" or header[b"priva.stage"] != (b"final" if producer else b"raw"):
        raise RuntimeError("native artifact kind/stage differs")
    started=time.monotonic();counts=Counter();references=[];mismatches=[];ordinal=0
    native=iter_genes(args.native,expand_samples=False)
    with args.mismatches.open("x") as ledger:
        for entry in entries:
            path=Path(entry["path"]);before=path.stat()
            sha=digest(path)
            if sha != entry["sha256"]:raise RuntimeError("reference checkpoint checksum changed: "+str(path))
            with path.open("rb") as stream:reference=pickle.load(stream)
            if list(reference) != entry["genes"]:raise RuntimeError("reference gene order disagrees with saved inventory")
            after=path.stat()
            if (before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns) != (after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns):
                raise RuntimeError("reference changed while loading")
            references.append(dict(path=str(path),sha256=sha,genes=len(reference)))
            for gene in entry["genes"]:
                expected=reference.pop(gene)
                if producer:
                    value={gene:expected};producer._finalize_scores_dict(value);expected=value[gene]
                try:
                    actual_gene,actual_ordinal,actual=next(native)
                    if (actual_gene,actual_ordinal) != (gene,ordinal):raise AssertionError("gene identity/order differs")
                    exact(expected,actual,counts)
                except Exception as error:
                    record=dict(gene_ordinal=ordinal,mismatch=str(error),exception=type(error).__name__)
                    mismatches.append(record);ledger.write(json.dumps(record)+"\n");ledger.flush()
                ordinal+=1;counts["genes"]+=1
                del expected
            del reference
        for _,extra_ordinal,_ in native:
            record=dict(gene_ordinal=extra_ordinal,mismatch="extra native gene")
            mismatches.append(record);ledger.write(json.dumps(record)+"\n")
    result=dict(passed=not mismatches,counts=dict(counts),mismatches=mismatches,references=references,
        native_sha256=digest(args.native),native_bytes=args.native.stat().st_size,elapsed_seconds=time.monotonic()-started,
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,original_collector_rerun=False,
        finalizer_sha256=PIN if producer else None)
    with args.report.open("x") as out:json.dump(result,out,indent=2);out.write("\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ["references","mismatches"]}),flush=True)
    if mismatches:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest="command",required=True)
    check=commands.add_parser("compare")
    for name in ["reference-list","native","report","mismatches"]:check.add_argument("--"+name,type=Path,required=True)
    check.add_argument("--finalize-source",type=Path)
    args=parser.parse_args();globals()[args.command](args)


if __name__=="__main__":main()
