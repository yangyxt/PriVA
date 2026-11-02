#!/usr/bin/env python3
"""Compare the optional read-only resource inventory against its Python source.

check --producer SCRIPT --binary PRIVA --registry JSON --work-root NEW_DIRECTORY
     [--extra-registry JSON ...]
oracle observes original exception events without changing the producer's code.
All fixtures are metadata only. Unsupported profiles remain measured mismatches.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import time
sys.dont_write_bytecode=True
SOURCE_SHA="dd0e1bc61466068c6f8ae772c39560771aa5d1722c4dce16fc7e3254c8e79b33"
SEEDS=(0,1,42,4294967295)

def pin(path):
    path=Path(path);data=path.read_bytes()
    return dict(path=str(path),bytes=len(data),sha256=hashlib.sha256(data).hexdigest())

def save(path,value):
    with Path(path).open("x") as out:json.dump(value,out,indent=2,sort_keys=True);out.write("\n")

def oracle(args):
    spec=importlib.util.spec_from_file_location("resource_oracle",args.producer)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    last={}
    def observe(frame,event,arg):
        if event=="exception" and frame.f_code.co_filename==str(args.producer):
            category,exception,_=arg
            if category is not SystemExit:last.update(category=category.__name__,message=str(exception))
        return observe
    sys.settrace(observe)
    try:returncode=module.main(args.arguments)
    except Exception as exc:
        # An interpreter-depth failure can disable the trace callback before
        # it records an exception; observe it after stack unwinding as well.
        sys.settrace(None);last.update(category=type(exc).__name__,message=str(exc));raise
    finally:
        sys.settrace(None)
        if last:save(args.exception,last)
    raise SystemExit(returncode)

def entry(*parents,lifecycle="derived",functions=()):
    return dict(lifecycle=lifecycle,config_keys=[],bundled_paths=[],depends_on=list(parents),
        refresh=dict(automatic=False,keep_function=bool(functions),entrypoints=list(functions)))

def registry(resources=None):
    return dict(schema_version=1,resources=resources or dict(a=entry(lifecycle="snapshot"),b=entry("a",functions=("build_b",)),c=entry("b")))

def catalogue():
    cases={}
    def add(name,document=None,**kw):cases[name]=dict(document=registry() if document is None else document,**kw)
    for command,options in (("validate",[]),("list",[]),("list",["--refreshable"]),("show",["--resource","b"]),("refresh-plan",["--resource","a"])):
        add("basic_"+command+("_filtered" if options==["--refreshable"] else ""),command=command,options=options)
    rich=registry();rich["resources"]["b"]["extras"]={"z":"α\t\n💙","a":[-0.0,1.0,1e-7,1e17,float('nan'),float('inf')],"n":None}
    add("show_unknown_fields",rich,command="show",options=["--resource","b"])
    for name,value in (("true",True),("float",1.0),("zero",0),("null",None),("string","1")):
        document=registry();document["schema_version"]=value;add("schema_"+name,document)
    document=registry();del document["schema_version"];add("schema_missing",document)
    for name,raw in (("null","null"),("list","[]"),("number","1"),("string",'"registry"'),("empty",""),
        ("space"," \r\n\t"),("unfinished",'{"schema_version":1'),("extra",'{} {}'),("colon",'{"schema_version" 1}'),
        ("duplicate",'{"schema_version":1,"schema_version":2}'),("nested_duplicate",'{"x":{"a":1,"a":2}}'),
        ("bad_escape",'{"a":"\\q"}'),("control",'{"a":"x\ty"}'),("bom",'\ufeff{}')):
        add("json_"+name,raw=raw)
    for name,value in (("empty",{}),("list",[]),("null",None)):
        document=registry();document["resources"]=value;add("resources_"+name,document)
    for name in ("A","a-","a\n","á","9a",""):
        add("id_"+str(len(cases)),registry({name:entry()}))
    for name,value in (("null",None),("list",[]),("number",1),("text","entry")):
        add("entry_"+name,registry(dict(a=value)))
    for name,value in (("bad","fixed"),("none",None),("list",[]),("object",{}),("bool",True)):
        document=registry();document["resources"]["a"]["lifecycle"]=value;add("lifecycle_"+name,document)
    for key in ("config_keys","bundled_paths","depends_on"):
        for name,value in (("none",None),("text","a"),("nonstrings",["a",1])):
            document=registry();document["resources"]["a"][key]=value;add(key+"_"+name,document)
        document=registry();del document["resources"]["a"][key];add(key+"_missing",document)
    for name,value in (("none",None),("list",[]),("bad_value",{"hg38":"a"}),("bad_member",{"hg38":[1]})):
        document=registry();document["resources"]["b"]["assembly_dependencies"]=value;add("assemblies_"+name,document)
    for name,value in (("none",None),("list",[]),("number",1),("text","refresh")):
        document=registry();document["resources"]["a"]["refresh"]=value;add("refresh_"+name,document)
    for key in ("automatic","keep_function"):
        for name,value in (("one",1),("none",None),("text","false")):
            document=registry();document["resources"]["a"]["refresh"][key]=value;add(key+"_"+name,document)
        document=registry();del document["resources"]["a"]["refresh"][key];add(key+"_missing",document)
    for name,value in (("none",None),("text","build_b"),("bad_member",["build_b",1])):
        document=registry();document["resources"]["b"]["refresh"]["entrypoints"]=value;add("entrypoints_"+name,document)
    document=registry();document["resources"]["b"]["refresh"]["keep_function"]=False;add("inconsistent_functions",document)
    document=registry();document["resources"]["a"]["refresh"]["keep_function"]=True;add("inconsistent_empty",document)
    add("missing_dependency",registry(dict(a=entry("unknown"))))
    add("seeded_unknown_dependencies",registry(dict(a=entry(*[f"missing_{i}" for i in range(60)]))))
    document=registry(dict(a=entry("missing_2","missing_1","missing_2")));document["resources"]["a"]["assembly_dependencies"]={"hg19":["missing_0","missing_9"],"hg38":["missing_1","missing_8"]};add("seeded_assembly_updates",document)
    add("self_cycle",registry(dict(a=entry("a"))))
    add("two_cycle",registry(dict(a=entry("b"),b=entry("a"))))
    add("seeded_cycle",registry(dict(a=entry("b","c"),b=entry("d"),c=entry("a"),d=entry("b"))))
    add("diamond_plan",registry(dict(a=entry(),d=entry("b","c"),c=entry("a"),b=entry("a"),z=entry("d"))),command="refresh-plan",options=["--resource","a"])
    document=registry();document["resources"]["c"]["assembly_dependencies"]={"hg38":["a"]};add("assembly_plan",document,command="refresh-plan",options=["--resource","a"])
    add("unknown_resource",command="show",options=["--resource","absent"])
    document=registry();document["resources"]["a"]["config_keys"]=None;add("validation_precedes_unknown",document,command="show",options=["--resource","absent"])
    for label,shell in (("simple","function build_b() {}\n"),("newline","function\n build_b \n (\n"),("unicode_space","function\u2003build_b\u001c(\n"),
        ("indented"," function build_b() {}\n"),("comment","# function build_b() {}\n"),("prefix","functionbuild_b() {}\n"),("invalid_name","function build_b-extra() {}\n")):
        add("routes_"+label,shell=shell,options=["--check-routes"])
    add("routes_missing_file",shell=None,options=["--check-routes"])
    add("routes_directory",shell_kind="directory",options=["--check-routes"])
    add("routes_utf8",shell_bytes=b"function build_b() {}\n\xff".hex(),options=["--check-routes"])
    document=registry();document["resources"]["b"]["install_entrypoints"]=["missing_"+str(i) for i in range(50)]+["build_b"]
    add("seeded_routes",document,shell="function build_b() {}\n",options=["--check-routes"])
    for label,value in (("none",None),("text","build_b"),("int",1),("object",{})):
        document=registry();document["resources"]["b"]["install_entrypoints"]=value
        add("install_routes_"+label,document,shell="function build_b() {}\n",options=["--check-routes"])
    for name,kind in (("missing","missing"),("directory","directory"),("not_directory","not_directory"),("denied","denied")):
        add("registry_"+name,registry_kind=kind)
    add("registry_utf8",raw_bytes=b'{"a":"\xff"}'.hex())
    add("registry_utf8_late",raw_bytes=(b" "*9000+b"\xff").hex())
    add("registry_cr",raw=json.dumps(registry(),indent=2).replace("\n","\r"))
    add("registry_symlink_routes",symlink=True,shell="function build_b() {}\n",options=["--check-routes"])
    add("registry_relative",relative=True)
    # Object hooks run only after an object has completely parsed.
    for name,raw in (("duplicate_then_syntax",'{"a":1,"a":2,"b":}'),("duplicate_then_nested_duplicate",'{"a":1,"a":2,"b":{"x":1,"x":2}}')):
        add("hook_"+name,raw=raw)
    for name,raw in {
        "child_before_syntax":'{"child":{"a":1,"a":2},"broken":}',
        "syntax_before_child":'{"first":1 "child":{"a":1,"a":2}}',
        "parent_not_closed":'{"a":1,"a":2',
        "child_closed_parent_not":'{"child":{"a":1,"a":2}',
        "empty_child_before_syntax":'{"a":1,"a":2,"child":{},"broken":}',
        "array_child_before_syntax":'{"a":1,"a":2,"child":[{"x":1,"x":2}],"broken":}',
        "unicode_offset":'{"α":1,"α":2,"broken":}',
        "newlines":'{\r\n"a":1,\r"a":2,\n"broken":}',
        "trailing_data":'{"a":1,"a":2} trailing',
        "quoted_braces":'{"a":1,"a":2,"child":{"x\\\"}":1,"x\\\"}":2}}',
        "multiple_duplicates":'{"a":1,"a":2,"b":1,"b":2}',
        "array_order":'[{"a":1,"a":2},{"x":1,"x":2}]',
        "nonfinite":'{"a":NaN,"a":Infinity,"child":{"x":-Infinity,"x":0}}',
        "wide_flat_prefix":json.dumps({"items":[{"v":i} for i in range(1000)]})[:-1]+',"items":0}',
        "deep_valid_prefix":'{"x":'*120+'{"a":1,"a":2}'+'}'*120,
    }.items():add("hook_"+name,raw=raw)
    for label,value in (("null",None),("int",1),("list",[]),("object",{})):
        document=registry();document["resources"]["b"]["install_entrypoints"]=["build_b",value]
        add("mixed_route_"+label,document,shell="function build_b() {}\n",options=["--check-routes"])
    for label,value in (("true",True),("false",False),("decimal",1.25),("minus_zero",-0.0),("nan",float('nan')),
        ("inf",float('inf')),("minus_inf",float('-inf')),("wide_integer",2**100)):
        document=registry();document["resources"]["b"]["install_entrypoints"]=["build_b",value]
        add("single_nonstring_"+label,document,shell="function build_b() {}\n",options=["--check-routes"])
    for label,values in (("missing_before_list",["missing",[]]),("null_before_list",[None,[]]),
        ("list_before_dict",[[],{}]),("dict_before_list",[{},[]]),("mixed_before_dict",[1,"missing",{}])):
        document=registry();document["resources"]["b"]["install_entrypoints"]=values
        add("unhashable_"+label,document,shell="function build_b() {}\n",options=["--check-routes"])
    document=registry();document["resources"]["a"]["install_entrypoints"]=["missing"];document["resources"]["b"]["install_entrypoints"]=[[]]
    add("earlier_resource_precedes_unhashable",document,shell="function build_b() {}\n",options=["--check-routes"])
    # These remain measured, never counted as matched failures or successes.
    for label,values in (("missing_string",["missing",None]),("two_nonstring",[None,1]),("equal_numbers",[1,True]),
        ("negative_zero",[-0.0,0.0]),("nan_duplicates",[float('nan'),float('nan')])):
        document=registry();document["resources"]["b"]["install_entrypoints"]=values
        add("unqualified_ambiguous_route_"+label,document,shell="function build_b() {}\n",options=["--check-routes"],
            unqualified="More than one possible missing candidate remains outside the string-set/single-candidate profile.")
    add("unqualified_duplicate_with_wide_integer",raw='{"a":1,"a":2,"large":'+str(2**160)+'}',
        unqualified="A competing duplicate hook with an unsupported JSON integer remains outside the represented JSON profile.")
    for label,raw in (("wide_int",'{"schema_version":1,"resources":'+json.dumps(registry()["resources"])+',"extra":'+str(2**160)+'}'),
        ("surrogate",json.dumps(registry())[:-1]+',"extra":"\\ud800"}'),
        ("nesting",json.dumps(registry())[:-1]+',"extra":'+"["*405+"0"+"]"*405+'}')):
        add("unqualified_json_"+label,raw=raw,unqualified="Existing JSON representation limit; Python success is retained as a mismatch.")
    chain={f"r{i}":entry(f"r{i+1}") if i<1099 else entry() for i in range(1100)}
    add("unqualified_python_recursion",registry(chain),unqualified="Native iterative DFS preserves graph order but does not emulate interpreter-frame recursion limits.")
    add("unqualified_random_seed",seed="random",unqualified="An independently randomized Python hash secret cannot be inferred.")
    add("unqualified_unset_seed",seed=None,unqualified="An unspecified Python hash secret cannot be inferred.")
    return cases

def prepare(folder,spec):
    folder.mkdir();root=folder/"input";root.mkdir();(root/"scripts").mkdir();path=root/"cache_resources.json"
    raw=bytes.fromhex(spec["raw_bytes"]) if "raw_bytes" in spec else spec.get("raw",json.dumps(spec["document"],ensure_ascii=True)).encode()
    path.write_bytes(raw)
    if "shell_bytes" in spec:(root/"scripts/install_utils.sh").write_bytes(bytes.fromhex(spec["shell_bytes"]))
    elif spec.get("shell_kind")=="directory":(root/"scripts/install_utils.sh").mkdir()
    elif spec.get("shell") is not None:(root/"scripts/install_utils.sh").write_text(spec["shell"])
    kind=spec.get("registry_kind")
    if kind in ("missing","directory"):path.unlink()
    if kind=="directory":path.mkdir()
    if kind=="not_directory":path=path/"child.json"
    if kind=="denied":path.chmod(0)
    if spec.get("symlink"):
        link=folder/"selected.json";link.symlink_to(path);path=link
    if spec.get("relative"):path=path.relative_to(folder)
    return path

def run(command,folder,label,env):
    started=time.monotonic()
    with (folder/(label+".stdout.log")).open("xb") as stdout,(folder/(label+".stderr.log")).open("xb") as stderr:
        child=subprocess.run(list(map(str,command)),cwd=folder,env=env,stdout=stdout,stderr=stderr,timeout=20)
    return dict(command=list(map(str,command)),returncode=child.returncode,wall_seconds=time.monotonic()-started,
        stdout=pin(folder/(label+".stdout.log")),stderr=pin(folder/(label+".stderr.log")))

def compare(args,folder,path,spec,seed):
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE="1")
    if seed is None:env.pop("PYTHONHASHSEED",None)
    else:env["PYTHONHASHSEED"]=str(seed)
    tail=[spec.get("command","validate"),"--registry",str(path),*spec.get("options",[])]
    expected=folder/"exception.json"
    def states():
        result={}
        for item in sorted((folder/"input").rglob("*")):
            value=item.lstat();record=dict(mode=stat.S_IMODE(value.st_mode),bytes=value.st_size,mtime_ns=value.st_mtime_ns)
            if item.is_symlink():record.update(kind="symlink",target=os.readlink(item))
            elif item.is_file():record.update(kind="file",sha256=pin(item)["sha256"] if record["mode"]&0o444 else None)
            else:record.update(kind="directory")
            result[str(item.relative_to(folder))]=record
        return result
    before=states()
    measurements={"python":run([sys.executable,Path(__file__).resolve(),"oracle","--producer",args.producer,"--exception",expected,"--",*tail],folder,"python",env),
        "native":run([args.binary,"cache-resources",*tail],folder,"native",env)}
    outcomes={}
    for arm,measurement in measurements.items():
        outcomes[arm]=dict(exit=measurement["returncode"],stdout_sha256=measurement["stdout"]["sha256"])
        if measurement["returncode"]:
            if arm=="python":outcomes[arm].update(json.loads(expected.read_text()) if expected.exists() else dict(category=None,message=Path(measurement["stderr"]["path"]).read_text()))
            else:
                text=Path(measurement["stderr"]["path"]).read_text();match=re.fullmatch(r"priva: (\w+): (.*)\n",text,re.S)
                outcomes[arm].update(category=match[1] if match else None,message=match[2] if match else text)
    preserved=before==states()
    result=dict(seed=seed,measurements=measurements,outcomes=outcomes,exact=outcomes["python"]==outcomes["native"] and preserved,inputs_preserved=preserved,
        stderr_bytes_equal=measurements["python"]["stderr"]["sha256"]==measurements["native"]["stderr"]["sha256"],
        unqualified=spec.get("unqualified"))
    save(folder/"comparison.json",result);return result

def check(args):
    if pin(args.producer)["sha256"]!=SOURCE_SHA:raise RuntimeError("Original source differs")
    args.work_root.mkdir(parents=True);cases={};real={}
    for name,spec in catalogue().items():
        for seed in (SEEDS if "seed" not in spec else (spec["seed"],)):
            folder=args.work_root/(name+"_"+str(seed));path=prepare(folder,spec)
            cases[folder.name]=compare(args,folder,path,spec,seed)
            if spec.get("registry_kind")=="denied":(folder/"input/cache_resources.json").chmod(0o644)
    for i,registry_path in enumerate([args.registry,*args.extra_registry]):
        source_root=registry_path.resolve().parent;root=args.work_root/("real_"+str(i));root.mkdir();(root/"scripts").mkdir()
        shutil.copyfile(registry_path,root/"cache_resources.json");shutil.copyfile(source_root/"scripts/install_utils.sh",root/"scripts/install_utils.sh")
        original={"registry":pin(registry_path),"installer":pin(source_root/"scripts/install_utils.sh")}
        resources=json.loads(registry_path.read_text())["resources"]
        specs=[dict(command="validate"),dict(command="validate",options=["--check-routes"]),dict(command="list"),dict(command="list",options=["--refreshable"])]
        specs += [dict(command=command,options=["--resource",name]) for name in resources for command in ("show","refresh-plan")]
        results={}
        for j,spec in enumerate(specs):
            folder=root/("command_"+str(j));folder.mkdir();results[str(j)]=compare(args,folder,root/"cache_resources.json",spec,42)
        if original!={"registry":pin(registry_path),"installer":pin(source_root/"scripts/install_utils.sh")}:raise RuntimeError("Reference inputs changed")
        if pin(root/"cache_resources.json")["sha256"]!=original["registry"]["sha256"] or pin(root/"scripts/install_utils.sh")["sha256"]!=original["installer"]["sha256"]:raise RuntimeError("Read-only inputs changed")
        real[str(i)]=dict(original=original,resources=len(resources),cases=results)
    admission={};selected=args.work_root/'real_0/cache_resources.json'
    for name,tail,expected in [('missing_command',[],2),('unknown_command',['unknown'],2),('help',['--help'],0),
        ('missing_resource',['show','--registry',selected],2),('unknown_flag',['validate','--registry',selected,'--refreshable'],2),
        ('flag_value',['list','--registry',selected,'--refreshable=yes'],2),('missing_registry_value',['list','--registry'],2)]:
        folder=args.work_root/('admission_'+name);folder.mkdir();env=dict(os.environ,PYTHONHASHSEED='42')
        reference=run([sys.executable,args.producer,*tail],folder,'python',env)
        native=run([args.binary,'cache-resources',*tail],folder,'native',env)
        admission[name]=dict(python=reference,native=native,expected=expected,passed=reference['returncode']==native['returncode']==expected)
    folder=args.work_root/'admission_explicit_registry_required';folder.mkdir()
    measured=run([args.binary,'cache-resources','validate'],folder,'native',dict(os.environ,PYTHONHASHSEED='42'))
    admission['explicit_registry_required']=dict(native=measured,expected=2,passed=measured['returncode']==2,
        interface_difference='Native requires the registry path; original Python has a source-file-relative default.')
    qualified={name:row for name,row in cases.items() if not row["unqualified"]}
    unqualified={name:row for name,row in cases.items() if row["unqualified"]}
    failed=[name for name,row in qualified.items() if not row["exact"]]
    failed += ["real_"+i+"_"+j for i,group in real.items() for j,row in group["cases"].items() if not row["exact"]]
    failed += ['admission_'+name for name,row in admission.items() if not row['passed']]
    report=dict(producer=pin(args.producer),binary=pin(args.binary),comparator=pin(Path(__file__)),cases=cases,real=real,admission=admission,
        qualified_count=len(qualified),qualified_successes=sum(row["outcomes"]["python"]["exit"]==0 for row in qualified.values()),
        real_count=sum(len(group["cases"]) for group in real.values()),unqualified_count=len(unqualified),
        unqualified_mismatches=[name for name,row in unqualified.items() if not row["exact"]],failed=failed,
        stderr_policy="Exit, exception category/message and stdout bytes are exact gates; stderr prefixes/tracebacks are separately recorded.",passed=not failed)
    save(args.work_root/"comparison.json",report)
    print(json.dumps({key:report[key] for key in ("qualified_count","qualified_successes","real_count","unqualified_count","unqualified_mismatches","failed","passed")}))
    if failed:raise SystemExit(1)

def main():
    parser=argparse.ArgumentParser(description=__doc__);subs=parser.add_subparsers(dest="procedure",required=True)
    p=subs.add_parser("check");p.add_argument("--producer",type=Path,required=True);p.add_argument("--binary",type=Path,required=True);p.add_argument("--registry",type=Path,required=True);p.add_argument("--extra-registry",type=Path,action="append",default=[]);p.add_argument("--work-root",type=Path,required=True);p.set_defaults(function=check)
    p=subs.add_parser("oracle");p.add_argument("--producer",type=Path,required=True);p.add_argument("--exception",type=Path,required=True);p.add_argument("arguments",nargs=argparse.REMAINDER);p.set_defaults(function=oracle)
    args=parser.parse_args()
    if args.procedure=="oracle" and args.arguments[:1]==["--"]:args.arguments=args.arguments[1:]
    args.function(args)

if __name__=="__main__":main()
