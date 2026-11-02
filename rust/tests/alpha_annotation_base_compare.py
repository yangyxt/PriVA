#!/usr/bin/env python3
"""Differential annotation-base procedures using tiny real BGZF/TBI fixtures.

  python alpha_annotation_base_compare.py check --repo SOURCE --binary PRIVA --work-root EMPTY_DIRECTORY

VEP and Perl runtime discovery use explicit fixtures. No annotation, model
generation, package transaction or download is performed. Retain all failures.
"""
from __future__ import annotations
import argparse
import fcntl
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from unittest import mock
from types import SimpleNamespace


def write(path,value):
    path.write_text(json.dumps(value,indent=2,sort_keys=True)+"\n")


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def observe(function):
    try:
        return dict(success=True,value=function())
    except Exception as error:
        return dict(success=False,category=type(error).__name__,message=str(error))


class Fixture:
    def __init__(self,root,repo):
        self.root=root;root.mkdir()
        self.cache=root/"cache";self.cache.mkdir()
        self.source=root/"science";self.scripts=self.source/"scripts";self.scripts.mkdir(parents=True)
        for name in ["alphamissense_annotation_base.py","alphamissense_cache.py","mapped_artifact_publish.py","biomart_nm_to_enst.py","biomart_ftp_nm_to_enst.py"]:
            shutil.copy2(repo/"scripts"/name,self.scripts/name)
        (self.scripts/"alphmis_tsv2vcf.py").write_text("# scientific converter fixture\n")
        (self.scripts/"am_pick_missense_constrained_segments.py").write_text("# unchanged motif producer fixture\n")
        (self.scripts/"install_utils.sh").write_text(
            "function basic_vep_annotation() {\n    vep --fork 22 --buffer_size 10000 --use_transcript_ref --total_length --domains\n}\n"
            'function AlphaMissense_anno() {\n    python "$convert_py" "$alphamissense_prescore" "$converted" && \\\n'
            '    basic_vep_annotation \\\n    -i "$converted" -a "$assembly" -v "$conservation_file" && \\\n    :\n}\n')
        self.bin=root/"vep_api";self.bin.mkdir()
        loaded={}
        for name,relative in [("Bio/EnsEMBL/VEP/Runner.pm","modules/Bio/EnsEMBL/VEP/Runner.pm"),("Bio/EnsEMBL/Registry.pm","Bio/EnsEMBL/Registry.pm")]:
            path=self.bin/relative;path.parent.mkdir(parents=True,exist_ok=True);path.write_text("# pinned API\n");loaded[name]=str(path)
        self.vep=self.bin/"vep";self.vep.write_text("#!/bin/sh\nprintf 'ensembl-vep : 116.3\\n'\n");self.vep.chmod(0o755)
        perl=self.bin/"perl";perl.write_text("#!"+sys.executable+"\nprint("+repr(json.dumps(dict(version="5.032001",loaded=loaded)))+")\n");perl.chmod(0o755)
        plugins=self.cache/"plugins";plugins.mkdir();(plugins/"Conservation.pm").write_text("# plugin\n")
        tree=self.cache/"VEP/homo_sapiens_merged/116_GRCh38";tree.mkdir(parents=True);(tree/"info.txt").write_text("cache metadata\n")
        self.cfg=dict(cache_root=str(self.cache),assembly="hg38",vep_installed_version="116",threads=22,
                      alphamissense_vcf=str(self.cache/"converted.vcf.gz"),
                      alphamissense_vep_base_vcf=str(self.cache/"complete.base.vcf.gz"),
                      alphamissense_vep_vcf=str(self.cache/"runtime.vcf.gz"),
                      alphamissense_prescore=str(self.cache/"publisher.tsv.gz"),
                      ref_genome=str(self.cache/"reference.fa"),conservation_file=str(self.cache/"gerp.bw"),
                      vep_plugins_dir=str(plugins),vep_cache_dir=str(self.cache/"VEP"))
        for name in ["publisher.tsv.gz","publisher.tsv.gz.tbi","reference.fa","reference.fa.fai","gerp.bw"]:(self.cache/name).write_bytes(b"unchanged fixture\n")
        self.config=root/"config.json";self.save()
        self.env=dict(os.environ,PATH=str(self.bin)+os.pathsep+os.environ["PATH"],PYTHONDONTWRITEBYTECODE="1")
        self.module=load(self.scripts/"alphamissense_annotation_base.py","reference_"+root.name)

    def save(self):
        write(self.config,self.cfg)

    def ref(self,function,*args,**kwargs):
        with mock.patch.dict(os.environ,self.env,clear=True):
            return getattr(self.module,function)(*args,**kwargs)

    def common(self):
        return self.config,self.source

    def vcf(self,path,annotated=True,mitochondrial=True,header=None):
        text=header or ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##contig=<ID=chrM,length=100>\n'
             +('##VEP="v116.3" assembly="GRCh38.p14"\n' if annotated else '')
             +'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        text+='chr1\t10\t.\tA\tC\t.\tPASS\t.\n'
        if mitochondrial:text+='chrM\t20\t.\tG\tT\t.\tPASS\t.\n'
        with path.open("wb") as stream:subprocess.run(["bgzip","-c"],input=text.encode(),stdout=stream,check=True)
        subprocess.run(["tabix","-p","vcf",str(path)],check=True,capture_output=True)
        return path

    def stage(self,mitochondrial=True,header=None):
        root=self.cache/("stage_"+str(len(list(self.cache.glob("stage_*")))));root.mkdir()
        return self.vcf(root/"converted.vcf.gz",False),self.vcf(root/"annotation.vcf.gz",mitochondrial=mitochondrial,header=header)

    def published(self):
        converted,annotation=self.stage()
        expected=self.ref("capture",*self.common())
        self.ref("publish",*self.common(),converted,annotation,expected)
        return expected


class Harness:
    def __init__(self,args):
        self.args=args;self.root=args.work_root;self.root.mkdir();self.rows=[];self.controls=[];self.serial=0
        # Actual importlib discovery on deliberately small six-package fixtures.
        # The full retained environment is qualified separately.
        runtime=self.root/"python_runtime";runtime.mkdir()
        for module,dist in [("numpy","numpy"),("pandas","pandas"),("scipy","scipy"),("sklearn","scikit-learn"),("pysam","pysam"),("joblib","joblib")]:
            (runtime/(module+".py")).write_text("__version__ = '1.2.3'\n")
            folder=runtime/(dist.replace("-","_")+"-1.2.3.dist-info");folder.mkdir()
            (folder/"METADATA").write_text("Metadata-Version: 2.1\nName: "+dist+"\nVersion: 1.2.3\n")
            (folder/"RECORD").write_text(module+".py,,\n"+folder.name+"/METADATA,,\n"+folder.name+"/RECORD,,\n")
        sys.path.insert(0,str(runtime));os.environ["PYTHONPATH"]=str(runtime)+os.pathsep+os.environ.get("PYTHONPATH","")
        os.environ["PYTHONDONTWRITEBYTECODE"]="1"

    def fixture(self,name):
        self.serial+=1
        return Fixture(self.root/(f"{self.serial:03}_"+name),self.args.repo)

    def native(self,f,command,extra=()):
        index=len(self.rows)
        argv=[str(self.args.binary),"alpha-annotation-base",command,"--config",str(f.config),"--source-root",str(f.source),
              "--python",sys.executable,"--helper-root",str(f.source),*map(str,extra)]
        result=subprocess.run(argv,capture_output=True,text=True,env=f.env)
        logs=self.root/"commands";logs.mkdir(exist_ok=True)
        label=f"{index:04}_{len(list(logs.glob('*.stderr.log'))):04}"
        (logs/(label+".stdout.log")).write_text(result.stdout);(logs/(label+".stderr.log")).write_text(result.stderr)
        value=dict(success=result.returncode==0,returncode=result.returncode,argv=argv,stdout=result.stdout)
        if result.returncode:
            message=result.stderr.strip().rsplit("priva: ",1)[-1]
            category,sep,message=message.partition(": ")
            value.update(category=category if sep else "unknown",message=message)
        return value

    def record(self,name,reference,native,*,value=False,message=True):
        passed=reference["success"]==native["success"] and native["returncode"]==(0 if reference["success"] else 1)
        if reference["success"] and value:
            expected=reference["value"]
            emitted=native["stdout"]
            passed=passed and emitted==json.dumps(expected,sort_keys=True)+"\n"
        if not reference["success"]:
            passed=passed and reference["category"]==native.get("category")
            if message:passed=passed and reference["message"]==native.get("message")
        self.rows.append(dict(name=name,passed=passed,reference=reference,native=native,
                              exact_message_checked=message and not reference["success"],value_checked=value))
        write(self.root/"progress.json",dict(cases=len(self.rows),passed=sum(r["passed"] for r in self.rows),rows=self.rows))
        return passed

    def compare(self,name,f,command,function,extra=(),value=False,message=True):
        reference=observe(function);native=self.native(f,command,extra)
        return self.record(name,reference,native,value=value,message=message)

    def success(self,name,f,command,extra=(),verify=None):
        native=self.native(f,command,extra)
        reference=observe(verify or (lambda:None))
        self.record(name,reference,native)
        self.rows[-1]["passed"] &= native["success"] and reference["success"]

    def path_and_capture(self):
        f=self.fixture("ordinary")
        self.compare("capture_exact",f,"capture",lambda:f.ref("capture",*f.common()),value=True)
        ref=observe(lambda:str(f.ref("settings",*f.common())[2]["base"]))
        native=self.native(f,"path");native["stdout"]=json.dumps(native["stdout"].rstrip("\n"))+"\n"
        self.record("path_exact",ref,native,value=True)
        output=f.root/"capture.json"
        self.success("capture_file",f,"capture",["--output",output],lambda:f.ref("capture",*f.common()))
        self.rows[-1]["passed"] &= json.loads(output.read_text())==f.ref("capture",*f.common()) and output.stat().st_mode&0o777==0o600
        self.compare("recover_absent",f,"recover",lambda:f.ref("recover",*f.common()))
        mutations={
            "root_relative":("cache_root","relative"),"root_empty":("cache_root",""),
            "converted_relative":("alphamissense_vcf","relative.vcf.gz"),"converted_suffix":("alphamissense_vcf",str(f.cache/"bad.vcf")),
            "destinations_equal":("alphamissense_vep_base_vcf",f.cfg["alphamissense_vcf"]),
            "outside_root":("alphamissense_vep_base_vcf",str(f.root/"outside.vcf.gz")),
            "parent_component":("alphamissense_vep_base_vcf",str(f.cache/"a/../bad.vcf.gz")),
            "null_path":("alphamissense_vcf",None),"numeric_path":("alphamissense_vcf",7),
            "assembly_unknown":("assembly","GRCh39"),"major_bad":("vep_installed_version","0116"),
            "threads_text_bad":("threads","many"),"threads_float":("threads",2.75),
            "conservation_none":("conservation_file",None),"conservation_null":("conservation_file","null"),
        }
        for name,(key,value)in mutations.items():
            old=f.cfg[key];f.cfg[key]=value;f.save()
            self.compare(name,f,"capture",lambda:f.ref("capture",*f.common()),value=True)
            f.cfg[key]=old;f.save()
        for content in ["[]","null","7","{}"]:
            f.config.write_text(content)
            self.compare("config_"+content,f,"path",lambda:f.ref("settings",*f.common()))
        f.save()
        target=Path(f.cfg["alphamissense_vcf"]);target.symlink_to(f.cache/"missing")
        self.compare("dangling_target",f,"path",lambda:f.ref("settings",*f.common()));target.unlink()
        target.write_text("same");other=Path(f.cfg["alphamissense_vep_base_vcf"]);os.link(target,other)
        self.compare("hardlinked_targets",f,"path",lambda:f.ref("settings",*f.common()));other.unlink();target.unlink()
        f.vep.write_text("#!/bin/sh\nexit 7\n")
        self.compare("vep_failure",f,"capture",lambda:f.ref("capture",*f.common()))

    def publication(self):
        f=self.fixture("publication");converted,annotation=f.stage()
        expected_bytes={key:path.read_bytes() for key,path in [
            ("converted",converted),("converted_index",Path(str(converted)+".tbi")),
            ("base",annotation),("base_index",Path(str(annotation)+".tbi"))]}
        expected=f.ref("capture",*f.common());capture=f.root/"expected.json";write(capture,expected)
        extra=["--staged-converted",converted,"--staged-annotation",annotation,"--expected-inputs",capture]
        self.success("native_publish",f,"publish",extra,lambda:f.ref("check",*f.common()))
        self.rows[-1]["passed"] &= all(Path(value["path"]).read_bytes()==expected_bytes[key]
            for key,value in f.module.output_state(f.module.settings(*f.common())[2]).items())
        self.compare("native_published_check",f,"check",lambda:f.ref("check",*f.common()))
        self.success("native_restore",f,"restore",verify=lambda:f.ref("check",*f.common()))
        base=Path(f.cfg["alphamissense_vep_base_vcf"]);runtime=Path(f.cfg["alphamissense_vep_vcf"])
        self.rows[-1]["passed"] &= base.read_bytes()==runtime.read_bytes() and base.stat().st_ino!=runtime.stat().st_ino
        self.compare("runtime_without_lineage",f,"runtime-associated",lambda:f.ref("runtime_associated",*f.common()))
        for field,value in [("threads",11),("assembly","hg19")]:
            old=f.cfg[field];f.cfg[field]=value;f.save()
            self.compare("changed_"+field,f,"check",lambda:f.ref("check",*f.common()))
            f.cfg[field]=old;f.save()
        receipt=f.module.receipt_path(base);original=receipt.read_bytes()
        for key,value in [("base_sha256","BAD"),("outputs",{}),("schema","unknown"),("science",{})]:
            data=json.loads(original);data[key]=value;write(receipt,data)
            self.compare("receipt_"+key,f,"check",lambda:f.ref("check",*f.common()))
        receipt.write_bytes(original)
        marker=f.module.pending_path(base);marker.write_text("{}")
        self.compare("pending_check",f,"check",lambda:f.ref("check",*f.common()));marker.unlink()
        for count in range(7):
            f=self.fixture("recover_prefix_"+str(count));converted,annotation=f.stage();expected=f.ref("capture",*f.common())
            with mock.patch.object(f.module,"commit",side_effect=RuntimeError("retained before commit")):
                result=observe(lambda:f.ref("publish",*f.common(),converted,annotation,expected))
            marker=f.module.pending_path(Path(f.cfg["alphamissense_vep_base_vcf"]));data=json.loads(marker.read_text())
            for op in data["operations"][:count]:os.replace(op["source"],op["target"])
            self.success("recover_prefix_"+str(count),f,"recover",verify=lambda:f.ref("check",*f.common()))
        f=self.fixture("recover_bad_order");converted,annotation=f.stage();expected=f.ref("capture",*f.common())
        with mock.patch.object(f.module,"commit",side_effect=RuntimeError("retained before commit")):
            observe(lambda:f.ref("publish",*f.common(),converted,annotation,expected))
        marker=f.module.pending_path(Path(f.cfg["alphamissense_vep_base_vcf"]));data=json.loads(marker.read_text())
        op=data["operations"][-1];os.replace(op["source"],op["target"])
        self.compare("recover_bad_order",f,"recover",lambda:f.ref("recover",*f.common()))
        f=self.fixture("lost_chromosome");converted,annotation=f.stage(mitochondrial=False);expected=f.ref("capture",*f.common());capture=f.root/"expected.json";write(capture,expected)
        self.compare("lost_chromosome",f,"publish",lambda:f.ref("publish",*f.common(),converted,annotation,expected),
                     ["--staged-converted",converted,"--staged-annotation",annotation,"--expected-inputs",capture])

    def headers(self):
        variants={
            "missing_vep":"##fileformat=VCFv4.2\n",
            "duplicate_vep":'##VEP="v116.3" assembly="GRCh38"\n'*2,
            "patch_mismatch":'##VEP="v116.4" assembly="GRCh38"\n',
            "major_mismatch":'##VEP="v115.3" assembly="GRCh38"\n',
            "wrong_assembly":'##VEP="v116.3" assembly="GRCh37"\n',
            "bad_patch_assembly":'##VEP="v116.3" assembly="GRCh38.p0"\n',
            "quoted_unclosed":'##VEP="v116.3 assembly=GRCh38\n',
            "duplicate_assembly":'##VEP="v116.3" assembly=GRCh38 assembly=GRCh38\n',
        }
        for name,header in variants.items():
            f=self.fixture(name);f.published()
            path=Path(f.cfg["alphamissense_vep_base_vcf"])
            path.write_text("##fileformat=VCFv4.2\n"+header+"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            receipt=f.module.receipt_path(path);data=json.loads(receipt.read_text())
            data["outputs"]=f.module.output_state(f.module.settings(*f.common())[2]);write(receipt,data)
            self.compare(name,f,"check",lambda:f.ref("check",*f.common()))
        normal=b'##fileformat=VCFv4.2\n##VEP="v116.3" assembly="GRCh38"\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        for name,payload in {
            "ordinary_gzip_header":gzip.compress(normal),
            "invalid_body_ignored":normal+b'not a record\n',
            "no_header_newline":normal[:-1],
            "crlf_header":normal.replace(b'\n',b'\r\n'),
            "invalid_utf8_vep":normal.replace(b'v116.3',b'v116.\xff'),
            "escaped_assembly":normal.replace(b'GRCh38',b'GRCh\\38'),
            "header_byte_limit":b'##'+b'a'*(4*1024*1024)+b'\n'+normal,
            "header_line_limit":b'##x\n'*50000+normal,
        }.items():
            f=self.fixture(name);f.published();path=Path(f.cfg["alphamissense_vep_base_vcf"]);path.write_bytes(payload)
            receipt=f.module.receipt_path(path);data=json.loads(receipt.read_text());data["outputs"]=f.module.output_state(f.module.settings(*f.common())[2]);write(receipt,data)
            self.compare(name,f,"check",lambda:f.ref("check",*f.common()))

    def adoption(self):
        for native_runtime in [False,True]:
            f=self.fixture("adopt_"+str(native_runtime));converted,annotation=f.stage();expected=f.ref("capture",*f.common())
            if native_runtime:
                for source,key in [(converted,"alphamissense_vcf"),(annotation,"alphamissense_vep_vcf")]:
                    shutil.copy2(source,f.cfg[key]);shutil.copy2(str(source)+".tbi",f.cfg[key]+".tbi")
                expected=f.ref("capture",*f.common());converted=Path(f.cfg["alphamissense_vcf"])
            record=dict(schema=f.module.SCHEMA,science=expected["science"],counts=f.ref("index_counts",annotation),provenance={"fixture":"completed"},
                        converted=dict(path=str(converted),sha256=f.module.digest(converted),index_sha256=f.module.digest(str(converted)+".tbi")),
                        annotation=dict(path=str(annotation),sha256=f.module.digest(annotation),index_sha256=f.module.digest(str(annotation)+".tbi")))
            if native_runtime:
                proof=f.root/"proof.json";write(proof,{"closed":"fixture"})
                record["native_runtime"]=dict(vcf=f.module.identity(f.cfg["alphamissense_vep_vcf"]),index=f.module.identity(f.cfg["alphamissense_vep_vcf"]+".tbi"),
                    annotation_base_sha256=f.module.digest(annotation),science_sha256=f.module.json_digest(expected["science"]),
                    closed_native_proof=dict(path=str(proof),sha256=f.module.digest(proof)))
            lineage=f.root/"lineage.json";write(lineage,record)
            self.success("adopt_"+str(native_runtime),f,"adopt",["--lineage",lineage],lambda:f.ref("check",*f.common()))
            if native_runtime:self.compare("associated_runtime",f,"runtime-associated",lambda:f.ref("runtime_associated",*f.common()))
            if native_runtime:
                proof.write_text("changed closed proof\n")
                self.compare("associated_proof_changed",f,"runtime-associated",lambda:f.ref("runtime_associated",*f.common()))

    def motifs(self):
        f=self.fixture("motifs");f.published()
        self.compare("motif_runtime",f,"motif-runtime",lambda:f.ref("motif_runtime"),value=True)
        self.compare("motif_binding",f,"motif-binding",lambda:f.ref("motif_binding",*f.common()),value=True)
        artifact=f.cache/"historical.pkl";artifact.write_bytes(b"opaque complete historical motif bytes\n")
        mapping=f.root/"custom.tsv";mapping.write_text("custom mapping fixture\n")
        biomart=load(f.scripts/"biomart_nm_to_enst.py","fixture_biomart")
        record=dict(mapping_identity=biomart.mapping_identity(mapping,[]),artifact_bytes=artifact.stat().st_size,artifact_sha256=f.module.digest(artifact),
                    scientific_input=f.ref("motif_binding",*f.common()))
        receipt=Path(str(artifact)+".nm_to_enst.json");write(receipt,record)
        self.compare("motif_check",f,"motif-check",lambda:f.ref("motif_check",*f.common(),artifact),["--artifact",artifact])
        lineage=dict(schema=f.module.MOTIF_SCHEMA,provenance={"fixture":"retained"},scientific_input=record["scientific_input"],
                     artifact_path=str(artifact),receipt_sha256=f.module.digest(receipt),**{k:record[k] for k in ["mapping_identity","artifact_bytes","artifact_sha256"]})
        lineage_file=f.root/"motif_lineage.json";write(lineage_file,lineage)
        target=f.cache/"adopted.pkl"
        extra=["--artifact",target,"--mapping",mapping,"--lineage",lineage_file]
        self.success("adopt_motif",f,"adopt-motif",extra,lambda:biomart.check_receipt(target,mapping,[]))
        self.rows[-1]["passed"] &= target.read_bytes()==artifact.read_bytes() and target.stat().st_ino!=artifact.stat().st_ino
        def unchanged_state(path):
            s=path.stat()
            return s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns,s.st_mode
        before=(unchanged_state(target),Path(str(target)+".nm_to_enst.json").read_bytes())
        self.success("adopt_motif_reuse",f,"adopt-motif",extra,lambda:biomart.check_receipt(target,mapping,[]))
        self.rows[-1]["passed"] &= before==(unchanged_state(target),Path(str(target)+".nm_to_enst.json").read_bytes())
        old=receipt.read_bytes();receipt.write_text("{}")
        target2=f.cache/"second.pkl"
        self.compare("adopt_motif_stale_historical_receipt",f,"adopt-motif",
                     lambda:f.ref("adopt_motif",*f.common(),target2,lineage_file,mapping,[]),
                     ["--artifact",target2,"--mapping",mapping,"--lineage",lineage_file])
        receipt.write_bytes(old)

    def boundaries(self):
        f=self.fixture("profile_edges")
        original=f.vep.read_text()
        for name,body in {
            "vep_no_colon":"printf 'ensembl-vep 116.3\\n'\n",
            "vep_colon_no_space":"printf 'ensembl-vep :116.3\\n'\n",
            "vep_missing_space":"printf 'ensembl-vep116.3\\n'\n",
            "vep_stderr_first":"printf 'ensembl-vep : 115.0\\n' >&2\nprintf 'ensembl-vep : 116.3\\n'\n",
            "vep_stdout_first":"printf 'ensembl-vep : 116.3\\n'\nprintf 'ensembl-vep : 115.0\\n' >&2\n",
            "vep_patch_prefix":"printf 'ensembl-vep : 116.3suffix\\n'\n",
        }.items():
            f.vep.write_text("#!/bin/sh\n"+body)
            self.compare(name,f,"capture",lambda:f.ref("capture",*f.common()),value=True)
        f.vep.write_text(original)
        script=f.scripts/"install_utils.sh";original=script.read_text()
        for name,body in {
            "missing_basic":original.replace("function basic_vep_annotation","function other"),
            "missing_recipe":original.replace("function AlphaMissense_anno","function other"),
            "missing_converter_call":original.replace('python "$convert_py"','python "$other"'),
            "missing_vep_call":original.replace("    vep --fork","    other --fork"),
            "multiline_blank_before_vep":original.replace("    vep --fork","\n \t\n    vep --fork"),
            "plugin_execution_start":original.replace("    vep --fork","    if [[ -f ${spliceai_snv_prescore} ]]; then\n        :\n    fi\n    log ignored\n    vep --fork"),
            "crlf_shell_source":original.replace("\n","\r\n"),
        }.items():
            script.write_bytes(body.encode())
            self.compare(name,f,"capture",lambda:f.ref("capture",*f.common()),value=True)
        script.write_text(original)
        for name,path in [
            ("prescore_missing",Path(f.cfg["alphamissense_prescore"])),
            ("reference_index_missing",Path(f.cfg["ref_genome"]+".fai")),
            ("plugin_missing",Path(f.cfg["vep_plugins_dir"])/"Conservation.pm"),
        ]:
            saved=path.read_bytes();path.unlink()
            self.compare(name,f,"capture",lambda:f.ref("capture",*f.common()))
            path.write_bytes(saved)
        tree=Path(f.cfg["vep_cache_dir"])/"homo_sapiens_merged/116_GRCh38"
        (tree/"info.txt").unlink()
        self.compare("cache_empty",f,"capture",lambda:f.ref("capture",*f.common()))
        (tree/"info.txt").write_text("restored")
        hidden=tree/".git";hidden.mkdir();(hidden/"ignore").write_text("unrelated")
        linked=tree/"linked";linked.symlink_to(f.root,target_is_directory=True)
        self.compare("tree_ignored_and_linked_dirs",f,"capture",lambda:f.ref("capture",*f.common()),value=True)
        linked.unlink()
        target=Path(f.cfg["alphamissense_vcf"]);target.write_bytes(b"")
        self.compare("empty_prior_output_allowed",f,"capture",lambda:f.ref("capture",*f.common()),value=True)
        target.unlink()
        for marker_kind in ["symlink","empty","list","oversize","wrong_schema"]:
            f=self.fixture("pending_"+marker_kind);base=Path(f.cfg["alphamissense_vep_base_vcf"])
            marker=f.module.pending_path(base)
            if marker_kind=="symlink":
                ordinary=f.root/"marker.json";ordinary.write_text("{}");marker.symlink_to(ordinary)
            else:
                marker.write_text({"empty":"","list":"[]","oversize":" "*4194304,"wrong_schema":'{"schema":"wrong"}'}[marker_kind])
            self.compare("pending_"+marker_kind,f,"recover",lambda:f.ref("recover",*f.common()),message=marker_kind!="empty")
        f=self.fixture("index_missing");converted,annotation=f.stage();expected=f.ref("capture",*f.common());capture=f.root/"expected.json";write(capture,expected)
        Path(str(converted)+".tbi").unlink()
        self.compare("index_missing",f,"publish",lambda:f.ref("publish",*f.common(),converted,annotation,expected),
                     ["--staged-converted",converted,"--staged-annotation",annotation,"--expected-inputs",capture])
        f=self.fixture("lock_contention");f.published()
        held=os.open(f.cache,os.O_RDONLY|os.O_DIRECTORY)
        try:
            fcntl.flock(held,fcntl.LOCK_EX|fcntl.LOCK_NB)
            self.compare("restore_lock_contention",f,"restore",lambda:f.ref("restore",*f.common()))
        finally:os.close(held)
        for change in ["stage_changed","prior_changed","source_escape","target_symlink"]:
            f=self.fixture(change);converted,annotation=f.stage();expected=f.ref("capture",*f.common())
            with mock.patch.object(f.module,"commit",side_effect=RuntimeError("retained before commit")):
                observe(lambda:f.ref("publish",*f.common(),converted,annotation,expected))
            marker=f.module.pending_path(Path(f.cfg["alphamissense_vep_base_vcf"]));data=json.loads(marker.read_text())
            op=data["operations"][0]
            if change=="stage_changed":Path(op["source"]).write_bytes(b"changed")
            elif change=="prior_changed":Path(op["target"]).write_bytes(b"changed")
            elif change=="source_escape":op["source"]=str(f.root/"outside");write(marker,data)
            else:
                outside=f.root/"linked";outside.write_text("kept");Path(op["target"]).symlink_to(outside)
            self.compare(change,f,"recover",lambda:f.ref("recover",*f.common()))

    def motif_recovery(self):
        for prefix in [0,1,2]:
            f=self.fixture("motif_prefix_"+str(prefix));f.published()
            mapping=f.source/"data/hgnc/ensembl_GRCh37_mart116_nm_to_enst.tsv";mapping.parent.mkdir(parents=True)
            # Bundled reuse has only existence/nonzero admission, even for this
            # deliberately non-table fixture. Mapping identity still binds its
            # bytes; it does not validate against a publisher checksum.
            mapping.write_text("deliberately tiny bundled fixture\n")
            biomart=load(f.scripts/"biomart_nm_to_enst.py","motif_biomart_"+str(prefix))
            artifact=f.cache/"historical.pkl";artifact.write_bytes(b"opaque historical motif")
            record=dict(mapping_identity=biomart.mapping_identity(mapping,[]),artifact_bytes=artifact.stat().st_size,
                        artifact_sha256=f.module.digest(artifact),scientific_input=f.ref("motif_binding",*f.common()))
            receipt=Path(str(artifact)+".nm_to_enst.json");write(receipt,record)
            lineage=dict(schema=f.module.MOTIF_SCHEMA,scientific_input=record["scientific_input"],provenance={"fixture":"closed"},
                         artifact_path=str(artifact),receipt_sha256=f.module.digest(receipt),**{k:record[k] for k in ["mapping_identity","artifact_bytes","artifact_sha256"]})
            lineage_file=f.root/"lineage.json";write(lineage_file,lineage)
            target=f.cache/"adopted.pkl";publisher=load(f.scripts/"mapped_artifact_publish.py","publisher_"+str(prefix))
            args=SimpleNamespace(artifact=target,mapping=mapping,input=[],cache_root=f.cache,output_dir=target.parent,
                assembly="hg38",config=f.config,producer=f.scripts/"am_pick_missense_constrained_segments.py",kind="motifs",
                annotation_base=Path(f.cfg["alphamissense_vep_base_vcf"]),native_lineage=lineage_file)
            with mock.patch.dict(os.environ,f.env,clear=True),mock.patch.object(publisher,"finish_pending",side_effect=RuntimeError("retain completed stage")):
                result=observe(lambda:publisher.adopt_native(args,artifact,lineage))
            marker=publisher.pending_path(target)
            if not marker.exists():raise RuntimeError("reference motif stage did not close: "+str(result))
            pending=json.loads(marker.read_text());staged=Path(pending["stage"])/target.name
            for src,dst in [(staged,target),(publisher.receipt_path(staged),publisher.receipt_path(target))][:prefix]:os.replace(src,dst)
            self.success("motif_recover_"+str(prefix),f,"adopt-motif",["--artifact",target,"--lineage",lineage_file,"--mapping",mapping],
                         lambda:biomart.check_receipt(target,mapping,[]))
            self.rows[-1]["passed"] &= not marker.exists() and target.read_bytes()==artifact.read_bytes()

    def finish(self):
        result=dict(cases=len(self.rows),successes=sum(r["reference"]["success"] for r in self.rows),
                    failures=sum(not r["reference"]["success"] for r in self.rows),
                    mismatch_count=sum(not r["passed"] for r in self.rows),rows=self.rows,
                    admission_controls=self.controls)
        write(self.root/"complete.json",result)
        print(json.dumps({k:v for k,v in result.items() if k not in ("rows","admission_controls")}))
        return result

    def admission(self):
        import yaml
        f=self.fixture("admission");f.published()
        original=f.config.read_text()
        for name,content in [
            ("yaml_alias", "unused: &alias 1\nunused_alias: *alias\n"+yaml.safe_dump(f.cfg)),
            ("yaml_flow",yaml.safe_dump(f.cfg,default_flow_style=True)),
        ]:
            f.config.write_text(content)
            reference=observe(lambda:f.ref("check",*f.common()));native=self.native(f,"check")
            self.controls.append(dict(name=name,reference=reference,native=native,
                                      source_supported=reference["success"],native_admitted=native["success"]))
        f.config.write_text(original)
        receipt=f.module.receipt_path(Path(f.cfg["alphamissense_vep_base_vcf"]));original=receipt.read_text()
        for name,value in [("json_extra_bigint",2**200),("json_extra_surrogate","\ud800")]:
            data=json.loads(original);data["unused_extra"]=value;write(receipt,data)
            reference=observe(lambda:f.ref("check",*f.common()));native=self.native(f,"check")
            self.controls.append(dict(name=name,reference=reference,native=native,
                                      source_supported=reference["success"],native_admitted=native["success"]))
        receipt.write_text(original)


def check(args):
    harness=Harness(args)
    try:
        harness.path_and_capture();harness.publication();harness.headers();harness.adoption();harness.motifs()
        harness.boundaries();harness.motif_recovery()
        harness.admission()
    finally:
        result=harness.finish()
    if result["mismatch_count"]:raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest="command",required=True)
    p=sub.add_parser("check")
    p.add_argument("--repo",type=Path,required=True);p.add_argument("--binary",type=Path,required=True);p.add_argument("--work-root",type=Path,required=True)
    args=parser.parse_args();check(args)


if __name__=="__main__":main()
