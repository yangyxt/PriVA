#!/usr/bin/env python3
"""Compare the optional Conservation guard with its unchanged Python source.

check --producer SCRIPT --binary PRIVA --work-root NEW_DIRECTORY
oracle is the internal reference launcher. It records exception type without
altering the original main(), validation, publication or subprocess behavior.
Scientific inputs are tiny real BigWigs made by the existing Perl fixture.
Tool-failure cases use explicitly identified fault-injection executables.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
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
PRODUCER_SHA="b742dafd0f961aba170bc1d637b0911914cf05702e2cbc8dd9e2c7a9ba7bfd73"
SENTINEL=b"existing target must be retained\n"

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def pin(path):
    path=Path(path)
    return dict(path=str(path),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())

def write_json(path,value):
    Path(path).write_text(json.dumps(value,sort_keys=True,indent=2)+"\n")

def oracle(args):
    module=load(args.producer,"original_conservation")
    if args.tool_prefix:
        # Explicit discovery fixture: emulate the source's selected interpreter
        # path without replacing its interpreter or installed scientific tools.
        module.sys.executable=str(args.tool_prefix/"python")
    def observed(function):
        def call(*positional,**keywords):
            try:return function(*positional,**keywords)
            except Exception as exc:
                write_json(args.exception,dict(category=type(exc).__name__,message=str(exc)))
                raise
        return call
    module.validate=observed(module.validate);module.publish=observed(module.publish)
    if args.procedure == "oracle-sum":
        try: print(json.dumps(module.bsd_sum(args.input)))
        except Exception as exc:
            write_json(args.exception,dict(category=type(exc).__name__,message=str(exc)))
            raise
        return
    raise SystemExit(module.main(args.arguments))

def stat_identity(path):
    value=Path(path).stat()
    return dict(device=value.st_dev,inode=value.st_ino,bytes=value.st_size,
        mtime_ns=value.st_mtime_ns,ctime_ns=value.st_ctime_ns)

def file_state(path):
    path=Path(path)
    if path.is_symlink():return dict(kind="symlink",target=os.readlink(path))
    if not path.exists():return dict(kind="absent")
    if path.is_dir():return dict(kind="directory",mode=stat.S_IMODE(path.stat().st_mode))
    return dict(kind="file",mode=stat.S_IMODE(path.stat().st_mode),bytes=path.stat().st_size,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest())

def run(command,folder,label,env):
    started=time.monotonic()
    with (folder/(label+".stdout.log")).open("wb") as out,(folder/(label+".stderr.log")).open("wb") as err:
        process=subprocess.run(list(map(str,command)),cwd=folder,stdout=out,stderr=err,env=env,timeout=150)
    return dict(command=list(map(str,command)),returncode=process.returncode,wall_seconds=time.monotonic()-started,
        stdout=pin(folder/(label+".stdout.log")),stderr=pin(folder/(label+".stderr.log")))

def executable(path,text,mode=0o755):
    path.write_text(text);path.chmod(mode);return path

def catalogue():
    cases={
        "valid":{},"unicode_path":{"filename":"selected_α.bw"},"quoted_path":{"filename":"selected \"quoted\" 'track'.bw"},"hg19":{"assembly":"hg19","track":"hg19"},"grch37":{"assembly":"GRCh37","track":"hg19"},
        "grch38":{"assembly":"GRCh38"},"extra_contig":{"track":"extra"},"alias_equal":{"track":"alias_equal"},
        "alias_conflict":{"track":"alias_conflict"},"mitochondrial_alias":{"track":"mito"},
        "wrong_track":{"track":"hg19"},"wrong_reference":{"reference_assembly":"hg19"},
        "shared_length":{"track":"wrong_chr2"},"source_empty":{"url":""},"custom_url":{"url":"https://example.invalid/custom.bw"},
        "unknown_assembly":{"assembly":"hg20"},"missing_input":{"input":"missing"},"empty_input":{"input":"empty"},
        "directory_input":{"input":"directory"},"corrupt_input":{"input":"corrupt"},"truncated_input":{"input":"truncated"},
        "missing_reference":{"reference":"missing"},"empty_reference":{"reference":"empty"},"directory_reference":{"reference":"directory"},
        "missing_index":{"fai":"missing"},"empty_index":{"fai":"empty"},"directory_index":{"fai":"directory"},
        "no_chr1":{"fai":"2\t100\t6\t1\t2\n"},"fai_duplicate":{"fai":"chr1\t248956422\t6\t1\t2\nchr1\t248956422\t6\t1\t2\n"},
        "fai_alias_conflict":{"fai":"1\t248956422\t6\t1\t2\nchr1\t248956421\t6\t1\t2\n"},
        "fai_alias_equal":{"fai":"1\t248956422\t6\t1\t2\nchr1\t248956422\t6\t1\t2\n"},
        "fai_cr":{"fai":"chr1\t248956422\t6\t1\t2\rchr2\t242193529\t6\t1\t2\r"},
        "fai_crlf":{"fai":"chr1\t248956422\t6\t1\t2\r\nchr2\t242193529\t6\t1\t2\r\n"},
        "fai_no_newline":{"fai":"chr1\t248956422\t6\t1\t2"},
        "fai_unicode_decimal":{"fai":"chr1\t٢٤٨٩٥٦٤٢٢\t６\t１\t２\n"},
        "fai_offset_nondecimal_digit":{"fai":"chr1\t248956422\t²\t1\t2\n"},
        "fai_length_nondecimal_digit":{"fai":"chr1\t²\t6\t1\t2\n"},
        "fai_unbounded_offset":{"fai":"chr1\t248956422\t"+"9"*5000+"\t1\t2\n"},
        "fai_big_integer":{"fai":"chr1\t248956422\t6\t1\t2\nextra\t"+"9"*200+"\t6\t1\t2\n"},
        "fai_integer_limit":{"fai":"chr1\t"+"9"*4301+"\t6\t1\t2\n"},
        "fai_bad_utf8":{"fai_bytes":b"chr1\t248956422\t6\t1\t2\n\xff".hex()},
        "fai_truncated_utf8":{"fai_bytes":b"chr1\t248956422\t6\t1\t2\n\xe2\x82".hex()},
        "fai_utf8_second_block":{"fai_bytes":(b"x"*8191+b"\n\xff").hex()},
        "existing_symlink":{"input":"symlink"},"relative_paths":{"relative":True},
        "default_size":{"url":"default"},"default_wrong_assembly":{"assembly":"hg19","track":"hg19","url":"default"},
        "publish":{"publish":True},"publish_existing":{"publish":True,"destination":"existing"},"publish_existing_corrupt":{"publish":True,"destination":"existing","input":"corrupt"},
        "publish_dangling":{"publish":True,"destination":"dangling"},"publish_directory":{"publish":True,"destination":"directory"},
        "publish_same":{"publish":True,"destination":"same"},"publish_symlink":{"publish":True,"input":"symlink"},
        "publish_missing_parent":{"publish":True,"destination":"missing_parent"},
        "publish_no_url":{"publish":True,"url":None},"publish_empty_url":{"publish":True,"url":""},
        "publish_link_denied":{"publish":True,"deny":"link"},"publish_unlink_denied":{"publish":True,"deny":"unlink"},
        "prefix_preferred":{"tools":"preferred"},"path_fallback":{"tools":"fallback"},"missing_perl":{"tools":"missing"},
        "nonexecutable_perl":{"tools":"nonexecutable"},"perl_exit":{"tools":"exit"},
        "perl_bad_json":{"tools":"json"},"perl_json_list":{"tools":"list"},"perl_json_null":{"tools":"null"},"perl_empty_dictionary":{"tools":"empty"},"perl_no_reader":{"tools":"no_reader"},
        "perl_unicode_stderr":{"tools":"unicode"},"perl_large_pipes":{"tools":"large"},
        "perl_timeout":{"tools":"timeout"},"track_changes_during_reader":{"tools":"change"},
    }
    for label,row in {"short":"chr1\t1", "empty_name":"\t248956422\t6\t1\t2", "zero":"chr1\t0\t6\t1\t2",
        "negative":"chr1\t-248956422\t6\t1\t2", "signed":"chr1\t+248956422\t6\t1\t2",
        "underscore":"chr1\t248_956_422\t6\t1\t2", "line_zero":"chr1\t248956422\t6\t0\t2",
        "line_short":"chr1\t248956422\t6\t2\t1", "blank":"", "extra_columns":"chr1\t248956422\t6\t1\t2\tignored"}.items():
        cases["fai_"+label]={"fai":row+"\n"}
    for name,details in {"actual":{},"locale":{"sum_tool":"locale"},"missing_input":{"input":"missing"},
        "directory":{"input":"directory"},"missing_tool":{"sum_tool":"missing"},"tool_exit":{"sum_tool":"exit"},
        "tool_signal":{"sum_tool":"signal"},"malformed":{"sum_tool":"malformed"},"unicode_digits":{"sum_tool":"unicode"},
        "nondecimal_digit":{"sum_tool":"nondecimal"},"invalid_utf8":{"sum_tool":"utf8"}}.items():
        cases["bsd_"+name]=dict(details,procedure="bsd-sum")
    return cases

def setup(folder,case,fixtures,guard):
    folder.mkdir();stage=folder/"staging";stage.mkdir();target_folder=folder/"target";target_folder.mkdir()
    path=stage/case.get("filename","selected.bw");shutil.copyfile(fixtures[case.get("track","valid")],path)
    reference=folder/"reference.fa";reference.write_bytes(b">chr1\nA\n")
    assembly=case.get("assembly","hg38");reference_assembly=case.get("reference_assembly",assembly)
    length=249250621 if reference_assembly in ("hg19","GRCh37") else 248956422
    fai=Path(str(reference)+".fai");fai.write_text(f"chr1\t{length}\t6\t1\t2\nchr2\t242193529\t6\t1\t2\n")
    if "fai_bytes" in case:fai.write_bytes(bytes.fromhex(case["fai_bytes"]))
    if "fai" in case:
        value=case["fai"]
        if value in ("missing","directory"):fai.unlink()
        if value=="directory":fai.mkdir()
        if value=="empty":fai.write_bytes(b"")
        elif value not in ("missing","directory"):fai.write_text(value)
    value=case.get("input")
    if value in ("missing","directory","symlink"):path.unlink()
    if value=="directory":path.mkdir()
    if value=="symlink":path.symlink_to(fixtures["valid"])
    if value=="empty":path.write_bytes(b"")
    if value=="corrupt":path.write_bytes(b"not a BigWig\n")
    if value=="truncated":path.write_bytes(path.read_bytes()[:path.stat().st_size//2])
    value=case.get("reference")
    if value in ("missing","directory"):reference.unlink()
    if value=="directory":reference.mkdir()
    if value=="empty":reference.write_bytes(b"")
    destination=target_folder/"published.bw"
    value=case.get("destination")
    if value=="existing":destination.write_bytes(SENTINEL)
    if value=="dangling":destination.symlink_to("missing.bw")
    if value=="directory":destination.mkdir()
    if value=="same":destination=path
    if value=="missing_parent":destination=folder/"missing"/"published.bw"
    if case.get("deny")=="link":target_folder.chmod(0o555)
    if case.get("deny")=="unlink":stage.chmod(0o555)
    prefix=Path(sys.executable).parent;env=dict(os.environ,PYTHONDONTWRITEBYTECODE="1")
    if case.get("tools"):
        prefix=folder/"prefix";prefix.mkdir();path_bin=folder/"path";path_bin.mkdir();env["PATH"]=str(path_bin)
        tool=case["tools"];perl=prefix/"perl";real=Path(sys.executable).parent/"perl"
        if tool=="preferred":
            perl.symlink_to(real);executable(path_bin/"perl","#!/bin/sh\nprintf WRONG >&2\nexit 77\n")
        elif tool=="fallback":(path_bin/"perl").symlink_to(real)
        elif tool=="nonexecutable":executable(perl,"not executable\n",0o644)
        elif tool=="exit":executable(perl,"#!/bin/sh\nprintf 'reader fixture failed\\n' >&2\nexit 73\n")
        elif tool=="json":executable(perl,"#!/bin/sh\nprintf 'not JSON\\n'\n")
        elif tool=="list":executable(perl,"#!/bin/sh\nprintf '[]\\n'\n")
        elif tool=="null":executable(perl,"#!/bin/sh\nprintf 'null\\n'\n")
        elif tool=="empty":executable(perl,"#!/bin/sh\nprintf '{\"chroms\":{}}\\n'\n")
        elif tool=="no_reader":executable(perl,"#!/bin/sh\nprintf '{\"chroms\":{\"1\":248956422}}\\n'\n")
        elif tool=="unicode":executable(perl,"#!/bin/sh\nprintf '\\377' >&2\nexit 73\n")
        elif tool=="change":
            executable(perl,f'#!{real}\nsystem {{q({real})}} q({real}), @ARGV; my $status=$? >> 8; utime(time()+1,time()+1,$ARGV[-1]); exit $status;\n')
        elif tool in ("timeout","large"):
            payload='sleep 130;' if tool=="timeout" else 'print STDERR "w" x 200000; print " " x 200000; print q({"chroms":{"1":248956422},"reader":"fixture","reader_module":"fixture"});'
            executable(perl,f'#!{real}\n{payload}\n')
    if case.get("sum_tool"):
        path_bin=folder/"sum_path";path_bin.mkdir();env["PATH"]=str(path_bin);env["LC_ALL"]="C.UTF-8"
        mode=case["sum_tool"];tool=path_bin/"sum"
        if mode=="locale":executable(tool,f'#!/bin/sh\n[ "$LC_ALL" = C ] || exit 74\nexec {shutil.which("sum")} "$@"\n')
        elif mode=="exit":executable(tool,"#!/bin/sh\nexit 73\n")
        elif mode=="signal":executable(tool,"#!/bin/sh\nkill -TERM $$\n")
        elif mode=="malformed":executable(tool,"#!/bin/sh\nprintf 'wrong\\n'\n")
        elif mode=="unicode":executable(tool,"#!/bin/sh\nprintf '١٢٣ ٤\\n'\n")
        elif mode=="nondecimal":executable(tool,"#!/bin/sh\nprintf '² 1\\n'\n")
        elif mode=="utf8":executable(tool,"#!/bin/sh\nprintf '\\377'\n")
    url=case.get("url","https://example.invalid/custom.bw" if case.get("publish") else None)
    if url=="default":url=guard.GRCH38_DEFAULT_URL
    supplied_input=str(path.relative_to(folder)) if case.get("relative") else str(path)
    supplied_reference=str(reference.relative_to(folder)) if case.get("relative") else str(reference)
    arguments=["--input",supplied_input,"--assembly",assembly,"--reference",supplied_reference]
    if url is not None:arguments += ["--source-url",url]
    if case.get("publish"):arguments += ["--publish-to",str(destination)]
    return dict(path=path,reference=reference,index=fai,destination=destination,prefix=prefix,env=env,arguments=arguments)

def normalize(value,folder):
    if isinstance(value,str):return value.replace(str(folder),"<ARM>")
    if isinstance(value,list):return [normalize(v,folder) for v in value]
    if isinstance(value,dict):return {k:normalize(v,folder) for k,v in value.items()}
    return value

def outcome(folder,measured,settings,reference_arm):
    stderr=Path(measured["stderr"]["path"]).read_text()
    result=dict(exit=measured["returncode"],stage=file_state(settings["path"]),destination=file_state(settings["destination"]))
    if measured["returncode"]==0:
        value=json.loads(Path(measured["stdout"]["path"]).read_text())
        if not isinstance(value,dict):
            result["receipt"]=value
            return result
        published=Path(value["path"])
        result["identities_exact"]=(value["identity"]==stat_identity(published) and value["reference"]["identity"]==stat_identity(settings["reference"]))
        result["identities_exact"] &= value["reference"]["fai_sha256"]==pin(settings["index"])["sha256"]
        value["identity"]="<VERIFIED_ACTUAL_STAT>";value["reference"]["identity"]="<VERIFIED_ACTUAL_STAT>"
        result["receipt"]=normalize(value,folder)
    elif reference_arm:
        result.update(json.loads((folder/"exception.json").read_text()))
    else:
        matched=re.fullmatch(r"priva: ([A-Za-z_][A-Za-z_0-9]*): (.*)\n",stderr,re.S)
        result.update(category=matched[1] if matched else None,message=matched[2] if matched else stderr)
    if "message" in result:result["message"]=normalize(result["message"],folder)
    result["stage"]=normalize(result["stage"],folder);result["destination"]=normalize(result["destination"],folder)
    return result

def check(args):
    if pin(args.producer)["sha256"]!=PRODUCER_SHA:raise RuntimeError("Oracle source differs")
    args.work_root.mkdir(parents=True,exist_ok=False);fixture_root=args.work_root/"fixtures";fixture_root.mkdir()
    helper=load(args.producer.parent.parent/"tests/test_conservation_resource_check.py","conservation_fixtures")
    guard=helper.guard;fixtures={}
    for name,chroms in {"valid":{"1":248956422,"2":242193529},"hg19":{"chr1":249250621},
        "wrong_chr2":{"chr1":248956422,"chr2":242193528},"extra":{"chr1":248956422,"extra":1000},
        "alias_equal":{"chr1":248956422,"1":248956422},"alias_conflict":{"chr1":248956422,"1":248956421},
        "mito":{"1":248956422,"chrM":16569}}.items():
        fixtures[name]=helper.make_bigwig(fixture_root,name,chroms)
    cases={}
    for name,case in catalogue().items():
        folder=args.work_root/name;folder.mkdir();settings={arm:setup(folder/arm,case,fixtures,guard) for arm in ("python","native")}
        commands={}
        for arm,s in settings.items():
            if case.get("procedure")=="bsd-sum":
                if arm=="python":commands[arm]=[sys.executable,str(Path(__file__).resolve()),"oracle-sum","--producer",str(args.producer),"--exception",str(folder/arm/"exception.json"),"--input",str(s["path"])]
                else:commands[arm]=[str(args.binary),"conservation-resource","bsd-sum","--input",str(s["path"])]
            elif arm=="python":commands[arm]=[sys.executable,str(Path(__file__).resolve()),"oracle","--producer",str(args.producer),"--exception",str(folder/arm/"exception.json"),"--tool-prefix",str(s["prefix"]),"--",*s["arguments"]]
            else:commands[arm]=[str(args.binary),"conservation-resource","publish" if case.get("publish") else "validate",*s["arguments"],"--tool-prefix",str(s["prefix"])]
        # Timeout probes run together; both retain the original 120-second limit.
        with ThreadPoolExecutor(max_workers=2) as pool:
            pending={arm:pool.submit(run,commands[arm],folder/arm,arm,settings[arm]["env"]) for arm in settings}
            measurements={arm:future.result() for arm,future in pending.items()}
        outcomes={arm:outcome(folder/arm,measurements[arm],settings[arm],arm=="python") for arm in settings}
        exact=outcomes["python"]==outcomes["native"]
        if outcomes["python"].get("identities_exact") is False:exact=False
        cases[name]=dict(spec=case,measurements=measurements,outcomes=outcomes,exact=exact,
            diagnostic_message_exact=outcomes["python"].get("message")==outcomes["native"].get("message"),
            stderr_bytes_equal=Path(measurements["python"]["stderr"]["path"]).read_bytes()==Path(measurements["native"]["stderr"]["path"]).read_bytes())
        write_json(folder/"comparison.json",cases[name])
        for arm in settings:
            for directory in (folder/arm/"staging",folder/arm/"target"):directory.chmod(0o755)
    # Native command admission is distinct from the Python helper's legacy
    # top-level interface. Missing named options remain exit 2.
    admission={}
    for name,tail,expected in [("no_procedure",[],2),("no_input",["validate"],2),("missing_prefix",["validate","--input","x","--assembly","hg38","--reference","y"],2),
        ("unknown",["unknown"],2),("help",["--help"],0)]:
        folder=args.work_root/("admission_"+name);folder.mkdir()
        measured=run([args.binary,"conservation-resource",*tail],folder,"native",dict(os.environ))
        admission[name]=dict(measurement=measured,expected=expected,passed=measured["returncode"]==expected)
    report=dict(producer=pin(args.producer),comparator=pin(Path(__file__)),binary=pin(args.binary),cases=cases,
        admission=admission,count=len(cases),successes=sum(c["outcomes"]["python"]["exit"]==0 for c in cases.values()),
        failed=[name for name,c in cases.items() if not c["exact"]]+["admission_"+name for name,c in admission.items() if not c["passed"]],
        stderr_format="Native prefix/category and reference traceback formatting are reported separately from exact exception message/category.")
    write_json(args.work_root/"comparison.json",report)
    print(json.dumps({key:report[key] for key in ("count","successes","failed")}))
    if report["failed"]:raise SystemExit(1)

def main():
    parser=argparse.ArgumentParser(description=__doc__);subs=parser.add_subparsers(dest="procedure",required=True)
    p=subs.add_parser("check");p.add_argument("--producer",type=Path,required=True);p.add_argument("--binary",type=Path,required=True);p.add_argument("--work-root",type=Path,required=True);p.set_defaults(function=check)
    p=subs.add_parser("oracle-sum");p.add_argument("--producer",type=Path,required=True);p.add_argument("--exception",type=Path,required=True);p.add_argument("--input",type=Path,required=True);p.set_defaults(function=oracle,tool_prefix=None)
    p=subs.add_parser("oracle");p.add_argument("--producer",type=Path,required=True);p.add_argument("--exception",type=Path,required=True);p.add_argument("--tool-prefix",type=Path);p.add_argument("arguments",nargs=argparse.REMAINDER);p.set_defaults(function=oracle)
    args=parser.parse_args()
    if args.procedure=="oracle" and args.arguments[:1]==["--"]:args.arguments=args.arguments[1:]
    args.function(args)

if __name__=="__main__":main()
