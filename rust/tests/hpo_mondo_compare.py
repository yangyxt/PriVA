#!/usr/bin/env python3
"""Compare opt-in builders with unchanged Python/pandas source.

Run `compare --native PATH --source REPO --inputs-report JSON --work-root NEW`.
The input report names already acquired, frozen release files and SHA256 pins.
No source download or installer call occurs. `oracle --request JSON` runs one
unchanged Python procedure; source-frame observation captures coverage metadata
without replacing its joins, classifiers, parsers or publishers. Every failure,
input, command, complete output and mismatch is retained in the work directory.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import resource
import signal
import struct
import subprocess
import sys
import time
import traceback


def identity(path):
    with path.open("rb") as stream:
        return {"path": str(path), "bytes": path.stat().st_size,
                "sha256": hashlib.file_digest(stream, "sha256").hexdigest()}


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True) + "\n")


def oracle(args):
    request = json.loads(args.request.read_text())
    sys.path.insert(0, str(Path(request["source"]) / "scripts"))
    import pandas as pd
    import collapse_HPO_anno as hpo
    import build_mondo_disease_scope as mondo
    captured = {}

    def profile(frame, event, value):
        if event != "return" or value is None:
            return
        name, state = frame.f_code.co_name, frame.f_locals
        if frame.f_code.co_filename == hpo.__file__:
            if name == "_read_tsv" and str(state["path"]) == request["paths"].get("genes-to-phenotype"):
                captured["input_gene_rows"] = len(value)
            if name == "build_gene_disease_hpo_assertions":
                captured.update({"input_hpoa_rows": len(state["hpoa_df"]),
                    "unmapped_gene_rows": int(state["unmapped_gene"].sum()),
                    "positive_hpoa_assertions": len(state["hpoa_assertions"]), "joined_rows": len(state["joined"]),
                    "matched_rows": int(state["matched"].sum()), "negative_only_rows": int(state["negative_only"].sum()),
                    "unknown_rows": int(state["unknown"].sum()), "assertion_rows": len(value),
                    "genes": int(value["gene_symbol"].nunique()), "diseases": int(value["disease_id"].nunique()),
                    "match_rate": state["match_rate"]})
        if frame.f_code.co_filename == mondo.__file__ and name == "build_disease_scope_registry":
            captured.update({"assertion_rows": len(state["assertions"]), "diseases": len(value),
                "mondo_terms": len(state["terms"]), "mapped_diseases": int(value["mondo_id"].ne("").sum()),
                "unmapped_diseases": int(value["mondo_id"].eq("").sum()),
                "override_diseases": sum(key in state["overrides"] for key in value["disease_id"]),
                "loaded_overrides": len(state["overrides"]), "mondo_version": state["mondo_version"],
                "scope_counts": {str(key): int(count) for key, count in value["disease_scope"].value_counts().items()},
                "priva_counts": {str(key): int(count) for key, count in value["priva_scope"].value_counts().items()},
                "review_counts": {str(key): int(count) for key, count in value["scope_review_status"].value_counts().items()}})

    paths = request["paths"]
    sys.setprofile(profile)
    try:
        if request["procedure"] == "hpo":
            result = hpo.write_gene_disease_hpo_assertions(paths["genes-to-phenotype"], paths["phenotype-hpoa"],
                paths["output"], minimum_match_rate=float(request.get("threshold", ".999")),
                disease_scope_registry=paths.get("disease-scope-registry"))
            if "disease-scope-registry" in paths:
                captured["annotated_rows"] = len(result)
        elif request["procedure"] == "mondo":
            registry, annotated = mondo.build_and_write_resources(paths["mondo-obo"], paths["hpo-assertions"],
                paths["registry-output"], manual_overrides=paths.get("manual-overrides"),
                annotated_assertions_output=paths.get("annotated-assertions-output"))
            if annotated is not None:
                captured["annotated_rows"] = len(annotated)
        else:
            assertions, registry = mondo._read_tsv(paths["hpo-assertions"]), mondo._read_tsv(paths["registry"])
            annotated = mondo.annotate_hpo_assertions(assertions, registry)
            mondo._write_tsv_atomic(annotated, paths["output"])
            captured.update(assertion_rows=len(annotated), registry_rows=len(registry), annotated_rows=len(annotated))
        sys.setprofile(None)
        print(json.dumps({"success": True, "metadata": captured, "pandas": pd.__version__}), flush=True)
    except Exception as error:
        sys.setprofile(None)
        traceback.print_exc()
        print(json.dumps({"success": False, "error_type": type(error).__name__, "error": str(error)}), flush=True)
        raise SystemExit(1)


def write_tsv(path, names, rows, *, comments="", newline="\n"):
    stream = io.StringIO(newline="")
    stream.write(comments)
    writer = csv.writer(stream, delimiter="\t", lineterminator=newline)
    writer.writerow(names)
    writer.writerows(rows)
    path.write_bytes(stream.getvalue().encode())
    return path


def payload(path):
    data = path.read_bytes()
    return gzip.decompress(data) if data.startswith(b"\x1f\x8b") else data


def gzip_header(path):
    data = path.read_bytes()
    if not data.startswith(b"\x1f\x8b"):
        return None
    if data[2] != 8 or data[3] & ~8:
        raise AssertionError("unexpected gzip header flags")
    name = data[10:data.index(b"\0", 10)].decode("latin1") if data[3] & 8 else None
    return {"mtime": struct.unpack_from("<I", data, 4)[0], "xfl": data[8], "os": data[9], "filename": name}


def child_limits(file_limit, memory_limit):
    if file_limit:
        signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
        resource.setrlimit(resource.RLIMIT_FSIZE, (64, 64))
    if memory_limit:
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
        resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
        resource.setrlimit(resource.RLIMIT_CPU, (10, 15))


class Comparison:
    def __init__(self, args):
        self.args = args
        self.root = args.work_root
        self.root.mkdir(parents=True)
        self.cases = []
        self.mismatches = []
        self.byte_comparisons = 0
        self.header_checks = 0
        self.compressed_differences = 0
        self.release_checks = {}

    def case(self, label, procedure, paths, outputs, *, threshold=".999", failure=False,
             exact_error=False, old=False, directory_output=None, file_limit=False, release=False,
             resource_probe=False):
        root = self.root / label
        root.mkdir()
        request = {"procedure": procedure, "source": str(self.args.source), "threshold": str(threshold)}
        record = {"label": label, "procedure": procedure, "expected_failure": failure, "issues": [],
                  "resource_probe": resource_probe}
        paths = {name: str(path) for name, path in paths.items()}
        results, actual_paths, commands = {}, {}, {}
        start = time.time()
        for arm in ("python", "native"):
            folder = root / arm
            folder.mkdir()
            arm_paths = paths | {key: str(folder / name) for key, name in outputs.items()}
            actual_paths[arm] = arm_paths
            for key, name in outputs.items():
                output = folder / name
                output.parent.mkdir(parents=True, exist_ok=True)
                if key == directory_output:
                    output.mkdir()
                elif old:
                    output.write_bytes(b"prior-complete-output\n")
                    output.chmod(0o664)
                # A previous attempt is outside this invocation's ownership.
                stale = output.parent / ("." + output.name + ".write-preserved")
                stale.mkdir()
                (stale / "marker").write_bytes(b"preserved")
            if arm == "python":
                request_path = root / "request.json"
                save(request_path, request | {"paths": arm_paths})
                command = [sys.executable, str(Path(__file__).resolve()), "oracle", "--request", str(request_path)]
            else:
                head = {"hpo": ["hpo-assertions", "build"], "mondo": ["mondo-scope", "build"], "annotate": ["mondo-scope", "annotate"]}[procedure]
                command = [str(self.args.native), *head]
                for key, value in arm_paths.items():
                    command.extend(["--" + key, value])
                if procedure == "hpo":
                    command.append("--minimum-match-rate=" + str(threshold))
            commands[arm] = command
            previous_umask = os.umask(0o027)
            try:
                result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=600,
                    env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1", POLARS_MAX_THREADS="2"),
                    preexec_fn=(lambda: child_limits(file_limit, resource_probe)) if file_limit or resource_probe else None)
                (folder / "stdout.log").write_bytes(result.stdout)
                (folder / "stderr.log").write_bytes(result.stderr)
            finally:
                os.umask(previous_umask)
            raw_stdout = (folder / "stdout.log").read_text()
            try:
                metadata = json.loads(raw_stdout) if raw_stdout.strip() else None
            except json.JSONDecodeError:
                metadata = None
            results[arm] = {"returncode": result.returncode, "metadata": metadata,
                            "stderr": (folder / "stderr.log").read_text()}
            if result.returncode not in ((0, 1) if arm == "python" else (0, 2)):
                record["issues"].append(arm + " abnormal process exit " + str(result.returncode))
            if arm == "python" and result.returncode == 1 and (metadata is None or metadata.get("success") is not False):
                record["issues"].append("oracle exception metadata missing")
            if arm == "python" and metadata and metadata.get("error_type") == "MemoryError" and not resource_probe:
                record["issues"].append("oracle resource failure is not a parser-equivalence result")
            if arm == "python" and failure is None:
                failure = result.returncode != 0
                record["expected_failure"] = failure
            if (result.returncode != 0) != failure:
                record["issues"].append(arm + " unexpected exit " + str(result.returncode))
            for key, name in outputs.items():
                stale = folder / ("." + name + ".write-preserved")
                if not (stale / "marker").is_file() or (stale / "marker").read_bytes() != b"preserved":
                    record["issues"].append(arm + " modified prior staging marker")
                expected_files = {"stdout.log", "stderr.log", *outputs.values()}
                expected_files.update("." + value + ".write-preserved" for value in outputs.values())
                if {entry.name for entry in folder.iterdir()} - expected_files:
                    record["issues"].append(arm + " left unexpected staging products")
        record["commands"] = commands
        record["results"] = results
        if resource_probe:
            # A resource-constrained reference failure is recorded, not promoted
            # to source-language parser parity or a caller activation gate.
            record["parity_qualified"] = False
        elif not failure:
            py = results["python"]["metadata"]
            native = results["native"]["metadata"]
            if py is None or native is None or not py.get("success") or py.get("metadata") != native:
                record["issues"].append("structured coverage metadata differs")
            if release and py is not None:
                self.release_checks[label] = {"python": py.get("metadata"), "native": native}
        elif exact_error:
            py = results["python"]["metadata"]
            error = results["native"]["stderr"].strip().removeprefix("priva: ")
            if py is None or py.get("error", "").rstrip("\n") != error:
                record["issues"].append("failure message differs")
        output_records = {}
        for key in outputs:
            py, native = (Path(actual_paths[arm][key]) for arm in ("python", "native"))
            state = {arm: ("file" if path.is_file() else "directory" if path.is_dir() else "absent")
                     for arm, path in (("python", py), ("native", native))}
            if state["python"] != state["native"]:
                record["issues"].append(key + " publication state differs")
            if not py.is_file() or not native.is_file():
                output_records[key] = state
                continue
            expected, observed = payload(py), payload(native)
            self.byte_comparisons += 1
            state.update(python=identity(py), native=identity(native), decompressed_bytes=len(expected),
                         decompressed_sha256=hashlib.sha256(expected).hexdigest())
            if expected != observed:
                at = next((index for index, (a, b) in enumerate(zip(expected, observed)) if a != b), min(len(expected), len(observed)))
                mismatch = {"first_different_byte": at, "python_bytes": len(expected), "native_bytes": len(observed),
                            "python_sha256": hashlib.sha256(expected).hexdigest(), "native_sha256": hashlib.sha256(observed).hexdigest()}
                save(root / (key + ".mismatch.json"), mismatch)
                record["issues"].append(key + " complete decompressed output differs")
            if py.stat().st_mode & 0o7777 != native.stat().st_mode & 0o7777:
                record["issues"].append(key + " output mode differs")
            ph, nh = gzip_header(py), gzip_header(native)
            state["headers"] = {"python": ph, "native": nh}
            if (ph is None) != (nh is None):
                record["issues"].append(key + " compression inference differs")
            if ph is not None and nh is not None:
                self.header_checks += 1
                if not (int(start) <= nh["mtime"] <= int(time.time()) and nh["xfl"] == ph["xfl"] == 2 and nh["os"] == ph["os"] == 255):
                    record["issues"].append(key + " gzip time/level/OS header differs")
                if procedure == "hpo":
                    if nh["filename"] != ph["filename"]:
                        record["issues"].append(key + " gzip filename policy differs")
                else:
                    prefix = "." + py.name + "."
                    if (ph["filename"] is None) != (nh["filename"] is None):
                        record["issues"].append(key + " gzip filename encoding policy differs")
                    for header in (ph, nh):
                        name = header["filename"]
                        if name is not None and (not name.startswith(prefix) or not name.endswith(".tmp") or not name[len(prefix):-4].isdigit()):
                            record["issues"].append(key + " gzip staged filename policy differs")
                self.compressed_differences += py.read_bytes() != native.read_bytes()
            output_records[key] = state
        record["outputs"] = output_records
        record["seconds"] = time.time() - start
        save(root / "case.json", record)
        self.cases.append({"label": label, "issues": record["issues"], "expected_failure": failure, "case": str(root / "case.json")})
        if record["issues"]:
            self.mismatches.append(self.cases[-1])
        return {name: Path(value) for name, value in actual_paths["python"].items()}, {name: Path(value) for name, value in actual_paths["native"].items()}


def fixtures(root):
    root.mkdir()
    gene_columns = ["gene_symbol", "disease_id", "hpo_id", "frequency", "extra"]
    hpoa_columns = ["database_id", "qualifier", "hpo_id", "frequency", "evidence", "reference"]
    genes = [
        [" G2 ", " D2 ", " HP:2 ", "summary", ""], ["G1", "D1", "HP:1", "4/5", ""],
        ["G1", "D1", "HP:1", "4/5", ""], ["G1", "D1", "HP:1", "6/7", ""],
        ["GNEG", "DNEG", "HP:3", "rare", ""], ["GBOTH", "D1", "HP:1", "", ""],
        ["", "D1", "HP:1", "", ""], ["-", "D1", "HP:1", "", ""],
        ["GUNKNOWN", "DMISSING", "HP:x", " \t ", ""], ["NA", "NA", "NaN", "N/A", ""],
        ["\u001c ÅΩ😀\u00a0", "DQUOTE", "HP:q", "", ""], ["G\nline", "DL", "HP:l", "", ""],
    ]
    hpoa = [
        ["D1", "", "HP:1", "1/2", "TAS", "PMID:1"], ["D1", "", "HP:1", "1/3", "PCS", "PMID:2"],
        ["D1", "", "HP:1", "1/2", "TAS", "PMID:1"], ["D1", "NOT", "HP:1", "", "", ""],
        ["D2", "", "HP:2", "", "", ""], ["DNEG", "NOT", "HP:3", "", "", ""],
        ["NA", "", "NaN", "N/A", "NaN", "NULL"], ["DQUOTE", "", "HP:q", "50%", "A\tB", '"PMID:#quoted"'],
        ["DL", "", "HP:l", "", "", "line1\nline2"],
    ]
    paths = {"genes-to-phenotype": write_tsv(root / "genes.tsv", gene_columns, genes),
             "phenotype-hpoa": write_tsv(root / "hpoa.tsv", hpoa_columns, hpoa, comments="#release fixture\n#second comment\n")}
    obo = ["format-version: 1.2", "data-version: fixture-v1", "[Term]", "id: MONDO:0003847", "name: hereditary", "[Term]", "id: MONDO:0005070", "name: neoplasm"]
    assertions = []
    assertion_columns = ["gene_symbol", "disease_id", "hpo_id", "frequency", "evidence", "reference"]
    for bits in range(64):
        disease = "OMIM:" + str(100000 + bits)
        obo.extend(["[Term]", "id: MONDO:" + str(9000000 + bits), "name: scope " + str(bits) + " Ω", "xref: " + disease])
        if bits & 4: obo.append("is_a: MONDO:0003847 ! hereditary")
        if bits & 8: obo.append("is_a: MONDO:0005070 ! neoplasm")
        if bits & 16: obo.append("subset: predisposition")
        if bits & 32: obo.append("subset: omim_susceptibility")
        hpos = ["HP:other"] + (["HP:0000006"] if bits & 1 else []) + (["HP:0001426"] if bits & 2 else [])
        for index, hpo in enumerate(hpos):
            assertions.append(["G" + str(bits), disease, hpo, "", "TAS", "PMID:" + str(index)])
    obo.extend(["[Term]", "id: MONDO:9900000", "name: cycle", "is_a: MONDO:9900001", "xref: Orphanet:cycle",
                "[Term]", "id: MONDO:9900001", "name: cycle", "is_a: MONDO:9900000", "is_a: MONDO:0003847",
                "xref: ORPHA:cycle", "subset: omim_susceptibility", "[Term]", "id: MONDO:9900002",
                "name: broad", "xref: OMIM:broad {source=MONDO:relatedTo}", "xref: OMIM:eq {source=MONDO:equivalentTo}",
                "xref: UMLS:ignored", "[Term]", "id: MONDO:9900003", "name: obsolete", "is_a: MONDO:0005070",
                "xref: OMIM:obsolete", "is_obsolete: true", "[Term]", "id: MONDO:9900004", "name: child", "is_a: MONDO:9900003",
                "xref: OMIM:child", "[Term]", "id: MONDO:9900005", "name: replaced", "xref: OMIM:old",
                "[Term]", "id: MONDO:9900005", "name: replacement", "xref: OMIM:new", "is_a: MONDO:unknown",
                "[Typedef]", "id: ignored", "data-version: fixture-v2", "[Term]", "id: NOTMONDO:1", "xref: OMIM:notmondo"])
    for disease in ["ORPHA:cycle", "OMIM:broad", "OMIM:eq", "OMIM:obsolete", "OMIM:child", "OMIM:old", "OMIM:new", "OMIM:notmondo", "", "NA"]:
        assertions.append(["EXTRA", disease, "HP:other", "", "", ""])
    for hpo in ["HP:0000007", "HP:0001417", "HP:0001419", "HP:0001423", "HP:0001427", "HP:0001450", "HP:0012275", "HP:0034341", "HP:0010982", "HP:0010983", "HP:0010984"]:
        assertions.extend([["MORE", "OMIM:100000", hpo, "", "", "R;R"], ["MORE", "OMIM:100000", hpo, "", "", "-"], ["MORE", "OMIM:100000", hpo, "", "", ""]])
    mondo_path = root / "mondo.obo"
    mondo_path.write_text("\n".join(obo) + "\n")
    assertion_path = write_tsv(root / "assertions.tsv", assertion_columns, list(reversed(assertions)))
    scopes = ["mendelian_non_neoplastic", "mendelian_cancer_predisposition", "mendelian_susceptibility_uncertain", "complex_or_non_monogenic", "complex_susceptibility_candidate", "somatic_oncogenic", "neoplastic_uncertain", "conflicting_inheritance", "uncertain"]
    override_columns = ["disease_id", "disease_scope", "priva_scope", "scope_review_status", "scope_evidence", "scope_reference", "scope_note"]
    overrides = [["OMIM:" + str(100000 + i), scope, ["include", "exclude", "review"][i % 3], "" if i % 2 else "reviewed", "" if i % 2 else "manual", "", "note Ω"] for i, scope in enumerate(scopes)]
    overrides.extend([[" OMIM:100009 ", "uncertain", "review", "", "", "", "before"], ["OMIM:100009", "somatic_oncogenic", "exclude", "final", "override", 'PMID:"x"', "after"], ["", "bad-but-skipped", "bad-but-skipped", "", "", "", ""], ["OMIM:unused", "uncertain", "review", "", "", "", ""]])
    override_path = write_tsv(root / "overrides.tsv", override_columns, overrides, comments="#manual overrides\n")
    return paths, {"mondo-obo": mondo_path, "hpo-assertions": assertion_path}, override_path, (gene_columns, genes, hpoa_columns, hpoa, assertion_columns, assertions, override_columns, overrides)


def compare(args):
    inventory = json.loads(args.inputs_report.read_text())
    inputs = {name: Path(pin["frozen"]["path"]) for name, pin in inventory["inputs"].items()}
    for name, pin in inventory["inputs"].items():
        if identity(inputs[name]) != pin["frozen"]:
            raise RuntimeError("frozen input identity changed: " + name)
    for name, pin in inventory["python_sources"].items():
        if identity(args.source / "scripts" / name)["sha256"] != pin["frozen"]["sha256"]:
            raise RuntimeError("unchanged oracle source differs: " + name)
    comparison = Comparison(args)
    base_py, base_native = comparison.case("release_hpo", "hpo", {"genes-to-phenotype": inputs["genes"], "phenotype-hpoa": inputs["hpoa"]}, {"output": "assertions.tsv.gz"}, release=True)
    scope_py, scope_native = comparison.case("release_mondo", "mondo", {"mondo-obo": inputs["mondo"], "hpo-assertions": base_py["output"], "manual-overrides": inputs["overrides"]}, {"registry-output": "scope.tsv.gz", "annotated-assertions-output": "scoped.tsv.gz"}, release=True)
    comparison.case("release_hpo_scoped", "hpo", {"genes-to-phenotype": inputs["genes"], "phenotype-hpoa": inputs["hpoa"], "disease-scope-registry": scope_py["registry-output"]}, {"output": "scoped.tsv.gz"}, release=True)
    comparison.case("release_native_chain", "mondo", {"mondo-obo": inputs["mondo"], "hpo-assertions": base_native["output"], "manual-overrides": inputs["overrides"]}, {"registry-output": "scope.tsv", "annotated-assertions-output": "scoped.tsv"}, release=True)
    fixture_root = args.work_root / "fixtures"
    hp, mp, overrides, raw = fixtures(fixture_root)
    gc, genes, hc, hpoa, ac, assertions, oc, override_rows = raw
    for suffix in (".tsv", ".tsv.gz", ".tsv.gz.tmp", ".tsv.GZ"):
        comparison.case("hpo_suffix_" + suffix.replace(".", "_"), "hpo", hp, {"output": "assertions" + suffix}, threshold=0, old=True)
    comparison.case("hpo_unicode_filename", "hpo", hp, {"output": "Ω😀.tsv.gz"}, threshold="nan")
    comparison.case("hpo_latin1_filename", "hpo", hp, {"output": "Å.tsv.gz"}, threshold="-inf")
    comparison.case("hpo_default_mismatch", "hpo", hp, {"output": "out.tsv.gz"}, failure=True, exact_error=True, old=True)
    comparison.case("hpo_infinite_threshold", "hpo", hp, {"output": "out.tsv"}, threshold="inf", failure=True, exact_error=True)
    empty_genes = write_tsv(fixture_root / "empty_genes.tsv", gc, [])
    empty_hpoa = write_tsv(fixture_root / "empty_hpoa.tsv", hc, [])
    comparison.case("hpo_empty", "hpo", {"genes-to-phenotype": empty_genes, "phenotype-hpoa": empty_hpoa}, {"output": "empty.tsv.gz"})
    negative_genes = write_tsv(fixture_root / "negative_genes.tsv", gc, [["G", "DNEG", "HP:3", "", ""]])
    comparison.case("hpo_negative_only", "hpo", hp | {"genes-to-phenotype": negative_genes}, {"output": "negative.tsv"})
    threshold_genes = write_tsv(fixture_root / "threshold_genes.tsv", gc, [["G", "D1", "HP:1", "", ""], ["U", "UNKNOWN", "HP:x", "summary", ""]])
    rate = 1.0 - 1.0 / 3.0
    for label, threshold, failed in (("equal", rate, False), ("below", math.nextafter(rate, -math.inf), False), ("above", math.nextafter(rate, math.inf), True)):
        comparison.case("hpo_threshold_" + label, "hpo", hp | {"genes-to-phenotype": threshold_genes}, {"output": "out.tsv"}, threshold=threshold, failure=failed, exact_error=failed)
    for label, rows in (("blank_disease", [["", " ", "HP:1", "", ""]]), ("dash_hpo", [["-", "D1", "-", "", ""]])):
        path = write_tsv(fixture_root / (label + ".tsv"), gc, rows)
        comparison.case("hpo_" + label, "hpo", hp | {"genes-to-phenotype": path}, {"output": "out.tsv.gz"}, failure=True, exact_error=True, old=True)
    for key, columns in (("genes-to-phenotype", ["gene_symbol"]), ("phenotype-hpoa", ["database_id"])):
        path = write_tsv(fixture_root / (key + "_missing.tsv"), columns, [])
        comparison.case("hpo_missing_" + key, "hpo", hp | {key: path}, {"output": "out.tsv"}, failure=True, exact_error=True)
    cr_genes = write_tsv(fixture_root / "cr_genes.tsv", gc, genes, newline="\r")
    comparison.case("hpo_cr_input", "hpo", hp | {"genes-to-phenotype": cr_genes}, {"output": "out.tsv"}, threshold=0)
    parser_rows = [["G0", "D1", "HP:1", "", ""], ["G1", "D1", "HP:1", "", ""], ["G2", "D1", "HP:1", "", ""]]
    for label, count, spaces, ending in (("first", 1, [0], "\r"), ("second", 2, [1], "\r"),
            ("both", 2, [0, 1], "\r"), ("third", 3, [2], "\r"), ("lf", 3, [0, 1, 2], "\n"),
            ("crlf", 3, [0, 1, 2], "\r\n")):
        rows = [list(row) for row in parser_rows[:count]]
        for index in spaces: rows[index][0] = " " + rows[index][0]
        path = write_tsv(fixture_root / ("parser_" + label + ".tsv"), gc, rows, newline=ending)
        comparison.case("hpo_parser_" + label, "hpo", hp | {"genes-to-phenotype": path}, {"output": "out.tsv"}, threshold=0, failure=None, exact_error=True)
    # Chunk boundaries and the two-record header call jointly bound the pinned
    # C tokenizer's backtracking. These inputs are ASCII so the byte offsets are
    # explicit; complete source outputs/errors are retained for every offset.
    boundary_paths = {}
    for offset in (-2, -1, 0, 1, 2):
        prefix = "\t".join(gc) + "\r" + "\t".join(parser_rows[0]) + "\r" + "G1\tD1\tHP:1\t\t"
        padding = 256 * 1024 + offset - len(prefix.encode()) - 1
        path = fixture_root / ("chunk_" + str(offset) + ".tsv")
        path.write_bytes((prefix + "x" * padding + "\r G2\tD1\tHP:1\t\t\r").encode())
        if offset <= 0:
            comparison.case("hpo_parser_chunk_" + str(offset), "hpo", hp | {"genes-to-phenotype": path}, {"output": "out.tsv"}, threshold=0, failure=None)
        else:
            boundary_paths[offset] = path
    commented = fixture_root / "comments.tsv"
    commented.write_bytes(hp["phenotype-hpoa"].read_bytes() + b"Dinline\t\tHP:i\tvalue#discard\tE\tR\n")
    comparison.case("hpo_inline_comments", "hpo", hp | {"phenotype-hpoa": commented}, {"output": "out.tsv"}, threshold=0)
    zipped = fixture_root / "members.tsv.gz"
    content = hp["phenotype-hpoa"].read_bytes()
    midpoint = len(content) // 2
    zipped.write_bytes(gzip.compress(content[:midpoint], mtime=0) + gzip.compress(content[midpoint:], mtime=0))
    comparison.case("hpo_gzip_members", "hpo", hp | {"phenotype-hpoa": zipped}, {"output": "out.tsv.gz"}, threshold=0)
    for label, data in (("truncated", zipped.read_bytes()[:-5]), ("crc", zipped.read_bytes()[:-8] + b"\x00\x00\x00\x00" + zipped.read_bytes()[-4:]), ("invalid_utf8", b"gene_symbol\tfrequency\tdisease_id\thpo_id\n\xff\t\tD\tH\n"), ("unclosed_quote", b'gene_symbol\tfrequency\tdisease_id\thpo_id\n"G\t\tD\tH\n')):
        path = fixture_root / (label + (".gz" if label in ("truncated", "crc") else ".tsv"))
        path.write_bytes(data)
        key = "phenotype-hpoa" if label in ("truncated", "crc") else "genes-to-phenotype"
        comparison.case("hpo_" + label, "hpo", hp | {key: path}, {"output": "out.tsv.gz"}, failure=True, old=True)
    comparison.case("hpo_missing_input", "hpo", hp | {"phenotype-hpoa": fixture_root / "absent.tsv"}, {"output": "out.tsv"}, failure=True, old=True)
    comparison.case("hpo_partial_write", "hpo", hp, {"output": "out.tsv"}, threshold=0, failure=True, old=True, file_limit=True)
    comparison.case("hpo_gzip_partial_write", "hpo", hp, {"output": "out.tsv.gz"}, threshold=0, failure=True, old=True, file_limit=True)
    comparison.case("hpo_rename_failure", "hpo", hp, {"output": "out.tsv.gz"}, threshold=0, failure=True, directory_output="output")
    registry_py, _ = comparison.case("mondo_all_64_classifications", "mondo", mp, {"registry-output": "registry.tsv", "annotated-assertions-output": "scoped.tsv"})
    for suffix in (".tsv", ".tsv.gz", ".tsv.gz.tmp", ".tsv.GZ", ".tsv.bz2"):
        comparison.case("mondo_suffix_" + suffix.replace(".", "_"), "mondo", mp | {"manual-overrides": overrides}, {"registry-output": "registry" + suffix, "annotated-assertions-output": "scoped" + suffix}, old=True)
    comparison.case("mondo_missing_overrides", "mondo", mp | {"manual-overrides": fixture_root / "absent-overrides.tsv"}, {"registry-output": "registry.tsv"})
    comparison.case("mondo_unicode_filename", "mondo", mp, {"registry-output": "Ω😀.tsv.gz"})
    zero_overrides = fixture_root / "zero-overrides.tsv"
    zero_overrides.write_bytes(b"")
    comparison.case("mondo_empty_overrides", "mondo", mp | {"manual-overrides": zero_overrides}, {"registry-output": "registry.tsv"})
    loop = fixture_root / "override-loop.tsv"
    loop.symlink_to(loop.name)
    comparison.case("mondo_loop_overrides", "mondo", mp | {"manual-overrides": loop}, {"registry-output": "registry.tsv"})
    comparison.case("mondo_not_directory_overrides", "mondo", mp | {"manual-overrides": overrides / "child"}, {"registry-output": "registry.tsv"})
    blocked = fixture_root / "blocked_override_parent"
    blocked.mkdir()
    (blocked / "overrides.tsv").write_bytes(overrides.read_bytes())
    blocked.chmod(0)
    try:
        comparison.case("mondo_override_stat_permission", "mondo", mp | {"manual-overrides": blocked / "overrides.tsv"}, {"registry-output": "registry.tsv.gz"}, failure=True, old=True)
    finally:
        blocked.chmod(0o700)
    empty_assertions = write_tsv(fixture_root / "empty_assertions.tsv", ac, [])
    comparison.case("mondo_empty_assertions", "mondo", mp | {"hpo-assertions": empty_assertions}, {"registry-output": "registry.tsv.gz", "annotated-assertions-output": "scoped.tsv.gz"})
    for label, rows in (("duplicate", override_rows + [override_rows[0]]), ("unused_duplicate", override_rows + [override_rows[-1]]), ("bad_scope", [["OMIM:unused", "unknown", "include", "", "", "", ""]]), ("bad_priva", [["OMIM:unused", "uncertain", "invalid", "", "", "", ""]])):
        path = write_tsv(fixture_root / (label + "_overrides.tsv"), oc, rows)
        comparison.case("mondo_overrides_" + label, "mondo", mp | {"manual-overrides": path}, {"registry-output": "registry.tsv.gz"}, failure=True, exact_error=True, old=True)
    path = write_tsv(fixture_root / "missing_override_columns.tsv", ["disease_id"], [])
    comparison.case("mondo_override_columns", "mondo", mp | {"manual-overrides": path}, {"registry-output": "out.tsv"}, failure=True, exact_error=True)
    path = write_tsv(fixture_root / "missing_assertion_columns.tsv", ["disease_id"], [])
    comparison.case("mondo_assertion_columns", "mondo", mp | {"hpo-assertions": path}, {"registry-output": "out.tsv"}, failure=True, exact_error=True)
    obo = mp["mondo-obo"].read_text()
    for label, data in (("no_version", obo.replace("data-version:", "other-version:")), ("missing_root", obo.replace("id: MONDO:0003847", "id: OTHER:root")), ("blank_parent", obo + "\n[Term]\nid: MONDO:BAD\nis_a: \n"), ("blank_subset", obo + "\n[Term]\nid: MONDO:BAD\nsubset: \n"), ("blank_xref", obo + "\n[Term]\nid: MONDO:BAD\nxref: \n")):
        path = fixture_root / (label + ".obo")
        path.write_text(data)
        comparison.case("mondo_" + label, "mondo", mp | {"mondo-obo": path}, {"registry-output": "out.tsv.gz"}, failure=True, exact_error=True, old=True)
    cr_obo = fixture_root / "cr.obo"
    cr_obo.write_bytes(obo.replace("\n", "\r").encode())
    comparison.case("mondo_cr_obo", "mondo", mp | {"mondo-obo": cr_obo}, {"registry-output": "registry.tsv"})
    comparison.case("mondo_partial_write", "mondo", mp, {"registry-output": "registry.tsv"}, failure=True, old=True, file_limit=True)
    comparison.case("mondo_gzip_partial_write", "mondo", mp, {"registry-output": "registry.tsv.gz"}, failure=True, old=True, file_limit=True)
    comparison.case("mondo_second_publication_failure", "mondo", mp, {"registry-output": "registry.tsv.gz", "annotated-assertions-output": "scoped.tsv.gz"}, failure=True, old=True, directory_output="annotated-assertions-output")
    scope_columns = ["mondo_id", "mondo_name", "disease_scope", "priva_scope", "scope_evidence", "scope_reference", "scope_review_status"]
    extra_assertions = write_tsv(fixture_root / "stale_scope.tsv", ac + scope_columns + ["extra"], [row + ["stale"] * 7 + ["extra\tΩ"] for row in assertions + assertions[:3]])
    ap = {"hpo-assertions": extra_assertions, "registry": registry_py["registry-output"]}
    comparison.case("annotate_stale_columns_multiplicity", "annotate", ap, {"output": "scoped.tsv.gz"})
    with registry_py["registry-output"].open(newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        columns, rows = next(reader), list(reader)
    for label, altered in (("unused_duplicate", rows + [["UNUSED"] + rows[0][1:]] * 2), ("missing", rows[1:]), ("blank_scope", [row[:11] + [""] + row[12:] if index == 0 else row for index, row in enumerate(rows)])):
        path = write_tsv(fixture_root / (label + "_registry.tsv"), columns, altered)
        comparison.case("annotate_" + label, "annotate", ap | {"registry": path}, {"output": "out.tsv.gz"}, failure=True, exact_error=True, old=True)
    literal_nan = write_tsv(fixture_root / "nan_scope.tsv", columns, [row[:11] + ["nan"] + row[12:] for row in rows])
    comparison.case("annotate_literal_nan_scope", "annotate", ap | {"registry": literal_nan}, {"output": "out.tsv"})
    bounded = Comparison(argparse.Namespace(**(vars(args) | {"work_root": args.work_root / "bounded_parser"})))
    for offset, path in boundary_paths.items():
        bounded.case("hpo_parser_chunk_" + str(offset), "hpo", hp | {"genes-to-phenotype": path}, {"output": "out.tsv"}, threshold=0, failure=None, resource_probe=True)
    save(args.work_root / "bounded_parser/cases.json", bounded.cases)
    for name, pin in inventory["inputs"].items():
        if identity(inputs[name]) != pin["frozen"]:
            raise RuntimeError("frozen input changed during comparison: " + name)
    summary = {"passed": not comparison.mismatches, "cases": len(comparison.cases),
               "qualified_cases_passed": not comparison.mismatches,
               "expected_failures": sum(case["expected_failure"] for case in comparison.cases),
               "complete_output_comparisons": comparison.byte_comparisons, "gzip_header_checks": comparison.header_checks,
               "permitted_compressed_byte_differences": comparison.compressed_differences,
               "mismatches": len(comparison.mismatches), "release_checks": comparison.release_checks,
               "classification_combinations": 64, "input_identities_rechecked": len(inputs),
               "caller_general_parity": False,
               "unqualified_resource_boundaries": {"cases": len(bounded.cases), "offsets_after_256_kib": [1, 2],
                    "address_space_limit_bytes": 2 * 1024**3, "cpu_soft_seconds": 10, "cpu_hard_seconds": 15,
                    "ledger": str(args.work_root / "bounded_parser/cases.json"), "abnormal_cases": bounded.mismatches,
                    "reason": "Unbounded pandas children were OOM-killed at these offsets. Bounded resource behavior is retained separately and is not parser exception equivalence."},
               "python": sys.version, "native": identity(args.native),
               "scope": "Fixed-release and defined edge cases: complete decompressed bytes, headers, row order/multiplicity, source-executed metadata, errors and atomic publication. Two resource-dependent parser boundaries remain unqualified; installer callers remain unchanged."}
    save(args.work_root / "cases.json", comparison.cases)
    save(args.work_root / "mismatches.json", comparison.mismatches)
    save(args.work_root / "summary.json", summary)
    print(json.dumps({key: value for key, value in summary.items() if key not in ("release_checks", "native")}), flush=True)
    if comparison.mismatches:
        raise RuntimeError(str(len(comparison.mismatches)) + " cases differ; complete evidence retained")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("compare")
    for name in ("native", "source", "inputs-report", "work-root"):
        command.add_argument("--" + name, type=Path, required=True)
    command.set_defaults(function=compare)
    command = commands.add_parser("oracle")
    command.add_argument("--request", type=Path, required=True)
    command.set_defaults(function=oracle)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
