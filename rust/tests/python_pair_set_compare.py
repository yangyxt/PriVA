#!/usr/bin/env python3
"""Compare exact two-text tuple hashes/set order with installed CPython 3.12.

compare --probe PROGRAM --work-root NEW_DIRECTORY
Also reruns the unchanged string-set oracle against the shared table. Seeds are
local subprocess profiles, never a production seed or a historical-seed claim.
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import random
import subprocess
import sys

SEEDS=[0,1,42,4294967295]


def encoded(value):return value.encode().hex() or "-"
def line(pair):return " ".join(map(encoded,pair))+"\n"


def cases():
    for size in [0,1,4,5,18,19,76,77,306,307,1228,1229,50000,50001,78643]:
        values=[("PFAM",f"PF{index:08d}") for index in range(size)]
        yield f"forward_{size}",values
        yield f"reverse_{size}",values[::-1]
    chars=["","a","é","α","𐀀","😀","\0","\n","\r","\t","Aé","Aα","A😀"]
    yield "unicode",list(itertools.product(chars,repeat=2))*2
    yield "duplicate_resize",[(f"DB{index%3}",f"ID{index%97}") for index in range(4000)]
    yield "component_order",[(str(index),str(index+1)) for index in range(97)]+[(str(index+1),str(index)) for index in range(97)]
    # These collisions are selected by the actual independent Python hash.
    selected=[]
    for index in range(20000):
        pair=("DB",f"collision{index}")
        if hash(pair)&255==0:selected.append(pair)
        if len(selected)==40:break
    if len(selected)!=40:raise RuntimeError("could not construct bounded collision case")
    yield "collisions",selected
    yield "collisions_reversed",selected[::-1]
    rng=random.Random(9091)
    for ordinal in range(256):
        values=[(rng.choice(chars)+str(rng.randrange(150)),rng.choice(chars)+str(rng.randrange(150)))
            for _ in range(rng.randrange(350))]
        yield f"random_{ordinal}",values


def reference(args):
    if not (sys.version_info[:2]==(3,12) and sys.byteorder=="little" and sys.hash_info.algorithm=="siphash13"
            and sys.hash_info.width==64 and sys.hash_info.cutoff==0):
        raise RuntimeError("reference requires CPython 3.12, 64-bit little-endian SipHash13")
    args.work_root.mkdir(parents=True,exist_ok=False);manifest=[]
    for name,values in cases():
        (args.work_root/(name+".input")).write_text("".join(map(line,values)))
        (args.work_root/(name+".hashes")).write_text("".join(str(hash(value))+"\n" for value in values))
        ordered=set()
        for value in values:ordered.add(value)
        (args.work_root/(name+".order")).write_text("".join(map(line,ordered)))
        manifest.append(dict(case=name,values=len(values),unique=len(ordered)))
    (args.work_root/"manifest.json").write_text(json.dumps(manifest))


def compare(args):
    args.work_root.mkdir(parents=True,exist_ok=False);rows=[]
    for seed in SEEDS:
        root=args.work_root/str(seed)
        with (args.work_root/f"reference_{seed}.log").open("x") as log:
            subprocess.run([sys.executable,__file__,"reference","--work-root",str(root)],
                env=dict(os.environ,PYTHONHASHSEED=str(seed)),stdout=log,stderr=subprocess.STDOUT,check=True)
        for case in json.loads((root/"manifest.json").read_text()):
            name=case["case"];output=root/(name+".native.order");hashes=root/(name+".native.hashes")
            with (root/(name+".stderr.log")).open("xb") as error:
                status=subprocess.run([str(args.probe),"pairs","--input",str(root/(name+".input")),
                    "--output",str(output),"--hashes",str(hashes),"--seed",str(seed)],
                    stdout=subprocess.DEVNULL,stderr=error,timeout=60).returncode
            rows.append(dict(case,seed=seed,returncode=status,
                equal_hashes=status==0 and hashes.read_bytes()==(root/(name+".hashes")).read_bytes(),
                equal_order=status==0 and output.read_bytes()==(root/(name+".order")).read_bytes()))
    with (args.work_root/"string_oracle.log").open("x") as log:
        status=subprocess.run([sys.executable,str(Path(__file__).with_name("python_set_order_compare.py")),
            "check","--binary",str(args.probe),"--work-root",str(args.work_root/"strings")],
            stdout=log,stderr=subprocess.STDOUT,timeout=600).returncode
    strings=json.loads((args.work_root/"strings/comparison.json").read_text())
    mismatch=[row for row in rows if not(row["equal_hashes"] and row["equal_order"])]
    (args.work_root/"mismatches.json").write_text(json.dumps(mismatch,indent=2)+"\n")
    result=dict(passed=not mismatch and status==0 and strings["mismatches"]==0,
        pair_cases=len(rows),pair_hash_values=sum(row["values"] for row in rows),pair_mismatches=len(mismatch),
        string_cases=strings["cases"],string_hash_values=strings["hash_values"],string_mismatches=strings["mismatches"],
        seeds=SEEDS,python_version=sys.version,probe_sha256=hashlib.sha256(args.probe.read_bytes()).hexdigest(),results=rows)
    (args.work_root/"report.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({key:value for key,value in result.items() if key!="results"}),flush=True)
    if not result["passed"]:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest="command",required=True)
    for name,function in [("reference",reference),("compare",compare)]:
        cmd=commands.add_parser(name);cmd.add_argument("--work-root",type=Path,required=True)
        if name=="compare":cmd.add_argument("--probe",type=Path,required=True)
        cmd.set_defaults(function=function)
    args=parser.parse_args();args.function(args)


if __name__=="__main__":main()
