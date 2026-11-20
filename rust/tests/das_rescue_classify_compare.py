#!/usr/bin/env python3
"""Compare DAS candidate decisions with unchanged sorted-reference ASTs.

compare --probe PROGRAM --reference SOURCE --work-root NEW_DIR
The caller-normalized first InterPro entry is supplied identically. Raw InterPro
lookup, whole-stage direct-rescue early returns and full DAS remain separate.
"""
import argparse
import ast
import copy
import hashlib
import itertools
import json
import logging
from pathlib import Path
import random
import struct
import subprocess
import time
import warnings
import numpy as np
warnings.filterwarnings("ignore",category=RuntimeWarning)

NAMES={"go":"rescue_bimodal_go_mf_domains","non":"rescue_non_gomf_domains"}
PIN="1a2b4efc312c41c6fb20afa3c1d90566da7cd0772d2f261c2f7f09d094fae834"
def bits(value):return struct.pack(">d",value).hex()
def hextext(value):return value.encode().hex() or "-"


def references(path):
    text=path.read_text()
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:raise RuntimeError("explicit sorted oracle source changed")
    tree=ast.parse(text);functions={};pins={}
    for family,name in NAMES.items():
        original=next(node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name==name)
        function=copy.deepcopy(original);prefix=[]
        for node in function.body:
            if isinstance(node,ast.Expr) and isinstance(node.value,ast.Call) and ast.unparse(node.value.func)=="logger.info":break
            prefix.append(node)
        returned="[candidates,candidate_composites,"+("directly_rescued,n_go_mf" if family=="go" else "set(),n_eligible_type")+"]"
        function.body=prefix+[ast.Return(value=ast.parse(returned,mode="eval").body)]
        scope={"np":np,"logger":logging.getLogger("oracle")}
        exec(compile(ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[])),str(path),"exec"),scope)
        functions[family]=scope[name]
        pins[name]=hashlib.sha256(ast.get_source_segment(text,original).encode()).hexdigest()
    return functions,pins


def cases():
    paths=["G:"+source+":A" for source in ["Gene3D","Superfamily","PANTHER","Pfam","PIRSF","SMART","PROSITE_profiles","Prints",
        "Coiled-coils_(Ncoils)","Low_complexity_(Seg)","Transmembrane_helices","Cleavage_site_(Signalp)","PROSITE_patterns","Other","pfam",""]]
    paths+=["Pfam:A","Pfam","",":Pfam:A:B","G:Other:Pfam:A","G:Pfam:","α:Pfam:δ"]
    kinds=["Domain","Conserved_site","Binding_site","Active_site","PTM","Homologous_superfamily","Repeat","Family","Unknown","domain","","δ"]
    arrays=[[],[.564],[np.nextafter(.564,0.)],[np.nextafter(.564,1.)],[.75],[np.nextafter(.75,0.)],
        [np.nextafter(.75,1.)],[0.],[-0.],[1.],[np.nan],[np.inf],[-np.inf],[0.,1.],[0.,.564,1.],
        [np.nan,.8],[.8,np.nan],[np.inf,-np.inf],[1.,1.,0.,0.],[-1.,2.,3.]]
    for family,path,kind,mf in itertools.product(NAMES,paths,kinds,[False,True]):
        for scores in arrays:
            yield (family,path,True,True,kind,mf,.564 if family=="go" else .75,scores)
    for family,path,present,entry,threshold,scores in itertools.product(NAMES,paths,[False,True],[False,True],
        [-np.inf,-0.,.564,.75,np.inf,np.nan],arrays[:10]):
        yield (family,path,present,entry,"Family",True,threshold,scores)
    rng=random.Random(0xD45)
    for _ in range(4000):
        family=rng.choice(list(NAMES));scores=[rng.random() for _ in range(rng.randrange(65))]
        yield (family,rng.choice(paths),True,rng.choice([False,True]),rng.choice(kinds),rng.choice([False,True]),
            rng.choice([rng.random(),.564,.75,np.nan]),scores)


def compare(args):
    args.work_root.mkdir(parents=True,exist_ok=False)
    functions,pins=references(args.reference);start=time.monotonic();count=values=0;by_family={}
    class Normalizer:
        def __init__(self,entry):self.entry=entry
        def query_interpro_entry_vep_anno(self,*args):return {"interpro_entries":self.entry}
    with (args.work_root/"input.bits").open("x") as incoming,(args.work_root/"python.bits").open("x") as truth:
        for family,path,present,entry,kind,mf,threshold,scores in cases():
            values+=len(scores);count+=1;by_family[family]=by_family.get(family,0)+1
            incoming.write(" ".join([family,hextext(path),str(int(present)),str(int(entry)),hextext(kind),str(int(mf)),bits(threshold),str(len(scores)),*map(bits,scores)])+"\n")
            # Tuple GO terms are intentional: DAS accesses g[1] directly.
            entries=[("IPR",kind,"short","name",[("GO","molecular_function" if mf else "biological_process","description")],"")] if entry else []
            arrays={path:np.asarray(scores,dtype=np.float64)} if present else {}
            comps={kind:kind,"Domain":"Domain"}
            candidate,selected,direct,eligible=functions[family]({path},arrays,"GLOBAL",Normalizer(entries),{},type_composites=comps,min_median_am=threshold)
            action="test" if candidate else "direct" if direct else "skip"
            composite=selected[0] if selected else ""
            truth.write(f"{int(bool(eligible))} {action} {hextext(composite)}\n")
    python_seconds=time.monotonic()-start;start=time.monotonic()
    with (args.work_root/"input.bits").open("rb") as incoming,(args.work_root/"native.bits").open("xb") as out,(args.work_root/"native.stderr.log").open("xb") as err:
        status=subprocess.run([str(args.probe),"classify"],stdin=incoming,stdout=out,stderr=err,timeout=300).returncode
    mismatches=0
    with (args.work_root/"python.bits").open() as expected,(args.work_root/"native.bits").open() as actual,(args.work_root/"mismatches.jsonl").open("x") as ledger:
        for index,(left,right) in enumerate(itertools.zip_longest(expected,actual)):
            if left!=right:mismatches+=1;ledger.write(json.dumps(dict(case=index,python=left,native=right))+"\n")
    result=dict(passed=status==0 and not mismatches,cases=count,values=values,by_family=by_family,mismatches=mismatches,
        native_returncode=status,python_seconds=python_seconds,native_seconds=time.monotonic()-start,
        numpy_version=np.__version__,reference_sha256=PIN,procedure_source_sha256=pins,
        probe_sha256=hashlib.sha256(args.probe.read_bytes()).hexdigest(),
        scope="Candidate action, composite type selection and pre-median eligibility counter for typed float64 arrays and an already-resolved ordered first InterPro entry; whole DAS remains unqualified")
    (args.work_root/"report.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result),flush=True)
    if not result["passed"]:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);commands=parser.add_subparsers(dest="command",required=True)
    cmd=commands.add_parser("compare")
    for name in ["probe","work-root"]:cmd.add_argument("--"+name,type=Path,required=True)
    cmd.add_argument("--reference",type=Path,default=Path(__file__).parents[2]/"scripts/am_pick_intolerant_domains.py")
    cmd.set_defaults(function=compare);args=parser.parse_args();args.function(args)


if __name__=="__main__":main()

