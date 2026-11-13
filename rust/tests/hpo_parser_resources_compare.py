#!/usr/bin/env python3
"""Bounded differential checks for pandas' bare-CR tokenizer rewind.

Run `compare --native PATH --source REPO --work-root NEW`. Both real command
lines run under the same 512 MiB address-space and 10/15 second CPU limits.
ParserError is required for allocation-cycle cases; MemoryError, a signal, a
timeout or an arbitrary nonzero exit is never accepted as parser equivalence.
Complete inputs, outputs, diagnostics and wait4 resource measurements remain
in the new work directory. No release input or installed environment changes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import sys
import time

MEMORY = 512 * 1024**2
ALLOCATION_ERROR = "Error tokenizing data. C error: out of memory"
BUFFER_ERROR = "Error tokenizing data. C error: Buffer overflow caught - possible malformed input file."


def pin(path):
    with path.open("rb") as stream:
        return {"path": str(path), "bytes": path.stat().st_size,
                "sha256": hashlib.file_digest(stream, "sha256").hexdigest()}


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def limits():
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    resource.setrlimit(resource.RLIMIT_AS, (MEMORY, MEMORY))
    resource.setrlimit(resource.RLIMIT_CPU, (10, 15))
    os.umask(0o027)


def measured(command, folder):
    started = time.monotonic()
    with (folder / "stdout.log").open("xb") as out, (folder / "stderr.log").open("xb") as err:
        process = subprocess.Popen(command, stdout=out, stderr=err, cwd=folder,
            preexec_fn=limits, env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
                OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", POLARS_MAX_THREADS="1"))
        _, status, usage = os.wait4(process.pid, 0)
        process.returncode = os.waitstatus_to_exitcode(status)
    diagnostic = (folder / "stderr.log").read_text()
    matches = re.findall(r"(?m)^(?:priva: )?([A-Za-z_][A-Za-z_.0-9]*Error): ([^\n]*)$", diagnostic)
    category, message = matches[-1] if matches else (None, None)
    return {"command": list(map(str, command)), "returncode": process.returncode,
        "error_type": category.rsplit(".", 1)[-1] if category else None, "error": message,
        "wall_seconds": time.monotonic() - started, "peak_rss_kib": usage.ru_maxrss,
        "user_seconds": usage.ru_utime, "system_seconds": usage.ru_stime,
        "stdout": pin(folder / "stdout.log"), "stderr": pin(folder / "stderr.log")}


def publication(path):
    return ({"state": "file", "bytes": path.stat().st_size,
        "sha256": pin(path)["sha256"], "mode": path.stat().st_mode & 0o7777}
        if path.is_file() else {"state": "absent"})


def payload(names, values, offset, *, newline="\r", bom=False):
    """Place the space after row two at a specified 256 KiB buffer offset."""
    width = len(names)
    first = ["G0", *values] + [""] * (width - len(values) - 1)
    second = ["G1", *values] + [""] * (width - len(values) - 1)
    third = [" G2", *values] + [""] * (width - len(values) - 1)
    prefix = ("\ufeff" if bom else "") + "\t".join(names) + "\r" + "\t".join(first) + "\r" + "\t".join(second)
    padding = 256 * 1024 + offset - len(prefix.encode()) - 1
    result = prefix + "x" * padding + "\r" + "\t".join(third) + "\r"
    return result.replace("\r", newline).encode()


class Comparison:
    def __init__(self, args):
        self.args = args
        self.cases = []
        self.previous = json.loads(args.reuse_summary.read_text()) if args.reuse_summary else None
        if self.previous:
            for name, source in self.previous["original_sources"].items():
                if pin(args.source / "scripts" / name)["sha256"] != source["sha256"]:
                    raise RuntimeError("saved Python source differs")

    def case(self, label, procedure, paths, *, expected_error=None, old=False):
        root = self.args.work_root / label
        root.mkdir()
        outputs = {"output": "assertions.tsv"} if procedure == "hpo" else {
            "registry-output": "registry.tsv", "annotated-assertions-output": "annotated.tsv"}
        arms = {}
        issues = []
        for arm in ("python", "native"):
            reused = None
            if arm == "python" and self.previous:
                selected_inputs = {key: pin(path)["sha256"] for key, path in paths.items()}
                candidates = [case for case in self.previous["cases"] if case["procedure"] == procedure
                    and {key: value["sha256"] for key, value in case["input_files"].items()} == selected_inputs]
                for candidate in candidates:
                    expected_before = {key: ({"state": "file", "bytes": len(b"previous complete output\n"),
                        "sha256": hashlib.sha256(b"previous complete output\n").hexdigest(), "mode": 0o664}
                        if old else {"state": "absent"}) for key in outputs}
                    if candidate["arms"]["python"]["before"] == expected_before:
                        reused = candidate
                        break
            if reused:
                prior = reused["arms"]["python"]
                for record in reused["input_files"].values():
                    if pin(Path(record["path"])) != record:
                        raise RuntimeError("saved Python input changed")
                for key in ("stdout", "stderr"):
                    record = prior["measurement"][key]
                    if pin(Path(record["path"])) != record:
                        raise RuntimeError("saved Python diagnostic changed")
                original_folder = Path(prior["measurement"]["stdout"]["path"]).parent
                retained_modes = {}
                for key, name in outputs.items():
                    stored = publication(original_folder / name)
                    producing = dict(prior["after"][key])
                    # Retention systems may make immutable artifacts read-only.
                    # The pinned receipt carries their producing publication
                    # mode; the complete retained content must still match.
                    retained_modes[key] = stored.pop("mode", None)
                    producing.pop("mode", None)
                    if stored != producing:
                        raise RuntimeError("saved Python publication changed")
                result = prior["measurement"]
                if result["returncode"] != (1 if expected_error else 0):
                    issues.append("reused python exact exit differs")
                if expected_error and (result["error_type"], result["error"]) != ("ParserError", expected_error):
                    issues.append("reused python exact exception category/message differs")
                if expected_error and prior["before"] != prior["after"]:
                    issues.append("reused python modified output on parser failure")
                arms[arm] = {**prior, "reused_summary": pin(self.args.reuse_summary),
                    "original_label": reused["label"], "retained_modes": retained_modes}
                continue
            folder = root / arm
            folder.mkdir()
            selected = {**paths, **{key: folder / name for key, name in outputs.items()}}
            before = {}
            for key, name in outputs.items():
                output = folder / name
                if old:
                    output.write_bytes(b"previous complete output\n")
                    output.chmod(0o664)
                stale = folder / ("." + name + ".write-preserved")
                stale.mkdir()
                (stale / "marker").write_bytes(b"preserve prior attempt\n")
                before[key] = publication(output)
            if arm == "python":
                if procedure == "hpo":
                    command = [sys.executable, self.args.source / "scripts/collapse_HPO_anno.py",
                        selected["genes-to-phenotype"], selected["phenotype-hpoa"], selected["output"],
                        "--minimum-match-rate=0"]
                else:
                    command = [sys.executable, self.args.source / "scripts/build_mondo_disease_scope.py"]
                    for key, value in selected.items():
                        command.extend(["--" + key, value])
            else:
                command = [self.args.native, "hpo-assertions" if procedure == "hpo" else "mondo-scope", "build"]
                for key, value in selected.items():
                    command.extend(["--" + key, value])
                if procedure == "hpo":
                    command.append("--minimum-match-rate=0")
            result = measured(list(map(str, command)), folder)
            after = {key: publication(folder / name) for key, name in outputs.items()}
            expected_exit = 1 if expected_error else 0
            if result["returncode"] != expected_exit:
                issues.append(arm + " exact exit differs")
            if expected_error and (result["error_type"], result["error"]) != ("ParserError", expected_error):
                issues.append(arm + " exact exception category/message differs")
            if expected_error and before != after:
                issues.append(arm + " modified output on parser failure")
            expected_entries = {"stdout.log", "stderr.log"}
            for key, name in outputs.items():
                stale = folder / ("." + name + ".write-preserved")
                expected_entries.add(stale.name)
                if (stale / "marker").read_bytes() != b"preserve prior attempt\n":
                    issues.append(arm + " changed another attempt's staging file")
                if after[key]["state"] == "file":
                    expected_entries.add(name)
            if {entry.name for entry in folder.iterdir()} != expected_entries:
                issues.append(arm + " staging cleanup differs")
            arms[arm] = {"measurement": result, "before": before, "after": after}
        if arms["python"]["after"] != arms["native"]["after"]:
            issues.append("complete output bytes, mode or publication state differs")
        value = {"label": label, "procedure": procedure, "expected_error": expected_error,
            "input_files": {key: pin(path) for key, path in paths.items()},
            "arms": arms, "issues": issues, "passed": not issues}
        save(root / "case.json", value)
        self.cases.append(value)


def compare(args):
    args.work_root.mkdir()
    fixtures = args.work_root / "inputs"
    fixtures.mkdir()
    hpoa = fixtures / "hpoa.tsv"
    hpoa.write_text("database_id\tqualifier\thpo_id\tfrequency\tevidence\treference\nD1\t\tHP:1\t1/2\tTAS\tPMID:1\n")
    obo = fixtures / "mondo.obo"
    obo.write_text("format-version: 1.2\ndata-version: fixture\n\n[Term]\nid: MONDO:0003847\nname: hereditary\nxref: D1\n"
                   "\n[Term]\nid: MONDO:0005070\nname: neoplasm\n")
    comparison = Comparison(args)
    for width in (4, 5, 8):
        names = ["gene_symbol", "disease_id", "hpo_id", "frequency"] + ["extra" + str(n) for n in range(width - 4)]
        for offset in sorted({1, 2, width}):
            path = fixtures / (f"width_{width}_offset_{offset}.tsv")
            path.write_bytes(payload(names, ["D1", "HP:1"], offset))
            comparison.case(path.stem, "hpo", {"genes-to-phenotype": path, "phenotype-hpoa": hpoa},
                expected_error=ALLOCATION_ERROR, old=(offset == width))
        for ending in ("\n", "\r\n"):
            name = "lf" if ending == "\n" else "crlf"
            path = fixtures / (f"width_{width}_{name}.tsv")
            path.write_bytes(payload(names, ["D1", "HP:1"], 2, newline=ending))
            comparison.case(path.stem, "hpo", {"genes-to-phenotype": path, "phenotype-hpoa": hpoa}, old=True)
    names = ["gene_symbol", "disease_id", "hpo_id", "frequency", "extra"]
    for label, data, error in [
        ("bom_boundary_success", payload(names, ["D1", "HP:1"], 2, bom=True), None),
        ("bom_short_cycle", payload(names, ["D1", "HP:1"], 3, bom=True), ALLOCATION_ERROR),
        ("full_width_cycle", ("\t".join(names) + "\r G0\tD1\tHP:1\t\t\r G1\tD1\tHP:1\t\t\r").encode(), BUFFER_ERROR),
        ("chunk_boundary_no_rewind", payload(names, ["D1", "HP:1"], 0), None),
    ]:
        path = fixtures / (label + ".tsv")
        path.write_bytes(data)
        comparison.case(label, "hpo", {"genes-to-phenotype": path, "phenotype-hpoa": hpoa}, expected_error=error, old=True)
    for label, suffix, extra, bom in (("accent", "é", 1, False), ("emoji", "😀", 3, False),
                                      ("bom_emoji", "😀", 5, True)):
        names = ["gene_symbol", "disease_id", "hpo_id", "frequency", "extra" + suffix]
        for after in (0, 1):
            path = fixtures / (label + ("_success.tsv" if after == 0 else "_cycle.tsv"))
            path.write_bytes(payload(names, ["D1", "HP:1"], extra + after, bom=bom))
            comparison.case(path.stem, "hpo", {"genes-to-phenotype": path, "phenotype-hpoa": hpoa},
                expected_error=ALLOCATION_ERROR if after else None, old=True)
    names = ["gene_symbol", "disease_id", "hpo_id", "frequency", "evidence", "reference", "extra"]
    for label, ending, error in (("mondo_short_cycle", "\r", ALLOCATION_ERROR),
                                  ("mondo_lf_control", "\n", None)):
        path = fixtures / (label + ".tsv")
        path.write_bytes(payload(names, ["D1", "HP:1", "1/2", "TAS", "PMID:1"], 2, newline=ending))
        comparison.case(label, "mondo", {"hpo-assertions": path, "mondo-obo": obo}, expected_error=error, old=True)
    report = {"passed": all(case["passed"] for case in comparison.cases), "cases": comparison.cases,
        "case_count": len(comparison.cases), "failure_cases": sum(bool(case["expected_error"]) for case in comparison.cases),
        "reused_python_cases": sum("reused_summary" in case["arms"]["python"] for case in comparison.cases),
        "address_space_limit_bytes": MEMORY, "cpu_limit_seconds": [10, 15], "native": pin(args.native),
        "original_sources": {name: pin(args.source / "scripts" / name) for name in
            ("collapse_HPO_anno.py", "build_mondo_disease_scope.py")}, "comparison_source": pin(Path(__file__))}
    save(args.work_root / "summary.json", report)
    print(json.dumps({key: report[key] for key in ("passed", "case_count", "failure_cases")}))
    if not report["passed"]:
        raise RuntimeError("parser resource differential mismatch; retained complete case receipts")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="procedure", required=True)
    sub = commands.add_parser("compare")
    sub.add_argument("--native", type=Path, required=True)
    sub.add_argument("--source", type=Path, required=True)
    sub.add_argument("--work-root", type=Path, required=True)
    sub.add_argument("--reuse-summary", type=Path, help="Verify and reuse saved Python measurements with identical inputs, source and prior outputs")
    sub.set_defaults(run=compare)
    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
