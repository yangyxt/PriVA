#!/usr/bin/env python3
"""Compare `priva gene-nonlof-mechanism` with build_gene_nonlof_mechanism_cache.py.

  python gene_nonlof_mechanism_compare.py check --binary PRIVA --repo CHECKOUT --work-root NEW_DIR [--jobs N]
  python gene_nonlof_mechanism_compare.py real --binary PRIVA --repo CHECKOUT --work-root NEW_DIR \
      --raw-inputs DIR --gofcards-cache JSON_GZ

Each case prepares two identical copies of its inputs and runs the script
and the native command on them. Exit status, the failure category and
message, printed JSON, the progress lines on stderr and every file left in
the copy are compared. Timestamps and the copy's own path are normalized;
gzip files are compared by their decompressed text. Sources are generated
(G2P CSV, ClinGen TSV, PanelApp JSON, the GoFCards workbook written with
openpyxl, small ClinVar VCV XML files) or taken from the checkout (HGNC
table, schema, GoFCards exact cache). Selected-URL cases are served by a
local HTTP server; no case reaches a publisher. `real` runs one build on
existing raw sources read in place, made fresh by a prepared manifest.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import copy
import gzip
import hashlib
import http.server
import json
import os
import re
import shutil
import socketserver
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

STAMP = re.compile(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(?:\.\d+)?\+00:00")
PROGRESS = ("checking_source=", "carrying_source=", "filter:", "hgnc_mapping:", "clinvar_vcv:",
            "json_schema_validation=", "json_output:", "panelapp_details_fetched=")
G2P_URL = "https://www.ebi.ac.uk/gene2phenotype/api/panel/all/download"
CLINGEN_URL = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
GOFCARDS_URL = "https://download.genemed.tech/upload/GainFunCards/gofcards_data_download.xlsx"
PANELAPP_URL = "https://panelapp.genomicsengland.co.uk/api/v1/panels/"
CLINVAR_RAW = "clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def normalize(value, arm):
    if isinstance(value, dict):
        return {key: normalize(item, arm) for key, item in value.items()}
    if isinstance(value, list):
        return [normalize(item, arm) for item in value]
    if isinstance(value, str):
        return STAMP.sub("<TIME>", value.replace(str(arm), "<ARM>"))
    return value


def outcome(command, arm, env):
    result = subprocess.run([str(part) for part in command], cwd=arm, env=env,
                            capture_output=True, text=True)
    lines = result.stderr.splitlines()
    value = {"exit": result.returncode}
    if result.returncode == 1 and lines:
        last = lines[-1]
        match = re.match(r"priva: ([A-Za-z]+): (.*)", last) or re.match(r"(?:[\w.]+\.)?([A-Za-z]+(?:Error|Exception|Exit)?): (.*)", last)
        if match:
            value["category"], value["message"] = match[1], match[2]
        else:
            value["category"], value["message"] = "SystemExit", last
        value["message"] = normalize(value["message"], arm)
    progress = []
    for line in lines:
        if line.startswith(PROGRESS):
            progress.append(normalize(line, arm))
        elif line.startswith("parse_error source="):
            progress.append(line.split(":", 1)[0])
    value["progress"] = progress
    if result.returncode == 0 and result.stdout.strip():
        try:
            value["stdout"] = normalize(json.loads(result.stdout), arm)
        except json.JSONDecodeError:
            value["stdout"] = normalize(result.stdout, arm)
    return value


def snapshot(arm, folder=None):
    state = {}
    for path in sorted((folder or arm).rglob("*")):
        if path.is_dir() or path.name.endswith(".log"):
            continue
        relative = str(path.relative_to(folder or arm))
        data = path.read_bytes()
        if path.name.endswith(".gz"):
            try:
                data = gzip.decompress(data)
            except OSError:
                pass
        if path.name.endswith((".json", ".json.gz")):
            try:
                body = json.loads(data)
                if isinstance(body, dict) and isinstance(body.get("sources"), dict):
                    for name, meta in body["sources"].items():
                        if isinstance(meta, dict) and str(meta.get("status", "")).endswith("_after_download_error"):
                            meta["error"] = "<ERROR>"
                        # A fetched PanelApp file records its own download time.
                        if name == "panelapp_all_panels" and isinstance(meta, dict) and \
                                meta.get("status") in ("downloaded", "updated", "checked_same_checksum"):
                            meta.update(sha256="<PANELAPP>", md5="<PANELAPP>")
                            hashes = (body.get("build_inputs") or {}).get("source_sha256") or {}
                            if name in hashes:
                                hashes[name] = "<PANELAPP>"
                state[relative] = normalize(body, arm)
                continue
            except (json.JSONDecodeError, UnicodeDecodeError):
                pass
        if path.name.endswith(".tsv"):
            text = data.decode(errors="replace").replace(str(arm), "<ARM>")
            rows = [line.split("\t") for line in text.split("\n")]
            if rows and rows[0] and rows[0][0] == "source" and "error" in rows[0]:
                where = rows[0].index("error")
                for row in rows[1:]:
                    if len(row) > where and row[1].endswith("_after_download_error"):
                        row[where] = "<ERROR>"
                    if len(row) > where and row[0] == "panelapp_all_panels" and \
                            row[1] in ("downloaded", "updated", "checked_same_checksum"):
                        for column in ("sha256", "md5"):
                            row[rows[0].index(column)] = "<PANELAPP>"
            state[relative] = STAMP.sub("<TIME>", "\n".join("\t".join(row) for row in rows))
            continue
        state[relative] = hashlib.sha256(data).hexdigest()
    return state


class Server:
    """A local HTTP server over one folder, for selected-URL cases."""

    def __init__(self, folder):
        handler = lambda *a, **k: http.server.SimpleHTTPRequestHandler(*a, directory=str(folder), **k)
        http.server.SimpleHTTPRequestHandler.log_message = lambda *a, **k: None
        self.server = socketserver.ThreadingTCPServer(("127.0.0.1", 0), handler)
        self.port = self.server.server_address[1]
        threading.Thread(target=self.server.serve_forever, daemon=True).start()

    def url(self, path):
        return f"http://127.0.0.1:{self.port}/{path}"

    def close(self):
        self.server.shutdown()
        self.server.server_close()


def gofcards_genes(repo, count=8):
    cache = json.loads(gzip.decompress((repo / "data/gofcards/gofcards_exact_gof.json.gz").read_bytes()))
    picked = []
    for symbol, gene in cache["genes"].items():
        for variant_id, variant in gene["variants"].items():
            genomic = ((variant.get("assemblies") or {}).get("hg19") or {}).get("genomic") or {}
            status = ((variant.get("record") or {}).get("eligibility") or {}).get("status")
            if genomic.get("pos") and (status == "eligible" and len(picked) < count or status != "eligible"
                                       and not any(p[4] != "eligible" for p in picked)):
                picked.append((symbol, genomic["chrom"], genomic["pos"], genomic["ref"], genomic["alt"], status))
                break
        if len(picked) > count:
            break
    return cache, picked


def write_workbook(path, picked, variant="plain"):
    import openpyxl
    from openpyxl.styles import NamedStyle
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = "total3161"
    header = ["Order numbe", "genesymbol", "chr", "hg19start", "hg19end", "ref", "alt", "transcript",
              "Function", "Disorder involved", "PMID", "Pscore", "Pathways proteins involved",
              "Animal model", "Cell model", "Function"]
    if variant == "unnamed":
        header[7] = None
    sheet.append(header)
    for index, (symbol, chrom, pos, ref, alt, _) in enumerate(picked, 1):
        sheet.append([index, symbol, chrom if index % 2 else f"chr{chrom}", pos, float(pos) if index % 3 == 0 else pos,
                      ref, alt, "NM_000001.1" if index % 2 else None, "gain of function" if index % 4 else "_",
                      "Disorder " + str(index) if index != 2 else "NA", 12345 + index if index % 2 else f"PMID:{index}; 99",
                      [6.5, 3, 1, "nan", "x", 5.0, 4.25, " 7 "][index % 8], "pathway A;B" if index % 2 else "-",
                      "Y" if index % 3 else "N", "y" if index % 2 else "", "dup"])
    sheet.append([99, "UNMATCHED1", "X", 1000, 1000, "A", "T", None, None, None, None, None, None, None, None])
    sheet.append([None, None, None])
    sheet.append([100, "lowercase_gene", "7", "2000", "2001", "AC", "-", "NM_2", "activating", "D", "", "N/A"])
    sheet.append([101, picked[0][0], picked[0][1], picked[0][2], picked[0][2], picked[0][3], picked[0][4], "T", "f", "d", "1", "#N/A"])
    late = sheet.max_row + 2
    sheet.cell(row=late, column=1, value=102)
    sheet.cell(row=late, column=2, value="LATEGENE")
    sheet.cell(row=late, column=3, value=True)
    if variant == "dates":
        cell = sheet.cell(row=2, column=11, value=45000.5)
        cell.number_format = "yyyy-mm-dd hh:mm"
        cell = sheet.cell(row=3, column=11, value=0.25)
        cell.number_format = "h:mm:ss"
    if variant == "wrong_sheet":
        sheet.title = "other"
    workbook.save(path)


def write_sources(raw, picked, repo, strict=False):
    raw.mkdir(parents=True, exist_ok=True)
    g2p = raw / "AllG2P.csv"
    genes = [p[0] for p in picked]
    lines = ["g2p id,gene symbol,hgnc id,disease name,allelic requirement,confidence,molecular mechanism,publications,gene symbol"]
    lines += [
        f"G2P1,{genes[0]},HGNC:1,\"Disease, one\",monoallelic,definitive,gain of function,123;456,x",
        f"G2P2,{genes[1]},HGNC:2,Disease two,biallelic,limited,loss of function,,y",
        f"G2P3,{genes[2]},HGNC:3,Disease three,monoallelic,strong,dominant negative,PMID:9,z",
        ",,,,,,,,",
        "G2P4,ALIASGENE,,NA,monoallelic,moderate,\"activating, increased dosage\",1 2 3,q",
        *([] if strict else ["G2P5,EXTRA,a,b,c,d,e,f,g,h,i"]),
        f"G2P6,{genes[3].lower()},HGNC:4,D4,monoallelic,refuted,dominant-negative effect,,",
    ]
    g2p.write_text("\ufeff" + "\n".join(lines) + "\n")
    clingen = raw / "ClinGen_gene_curation_list_GRCh38.tsv"
    header = ["#Gene Symbol", "Gene ID", "Haploinsufficiency Score", "Haploinsufficiency Description",
              "Haploinsufficiency PMID1", "Triplosensitivity Score", "Triplosensitivity Description",
              "Triplosensitivity PMID1", "Triplosensitivity PMID2"]
    rows = [[genes[4], "100", "3", "Sufficient evidence", "111", "3", "Sufficient evidence for dosage pathogenicity", "222", "111"],
            [genes[5], "101", "", "", "", "2", "Emerging", "", ""],
            ["NOTS", "102", "1", "Little evidence", "", "0", "No evidence", "", ""],
            ["", "103", "3", "x", "", "3", "y", "", ""]]
    clingen.write_text("#ClinGen dosage\n#Generated\n" + "\t".join(header) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")
    panels = {"downloaded_at_utc": now(), "list_url": PANELAPP_URL, "panel_count": 2, "panels": [
        {"id": 1, "name": "Panel one", "version": "1.2", "genes": [
            {"gene_data": {"gene_symbol": genes[6]}, "mode_of_pathogenicity":
             "Loss-of-function variants (as defined in pop up message) DO NOT cause this phenotype - please provide details in the comments",
             "mode_of_inheritance": "MONOALLELIC", "publications": ["1", "", "PMID 2"], "confidence_level": "3"},
            {"gene_data": None, "entity_name": genes[7] if len(genes) > 7 else "EXTRA2", "mode_of_pathogenicity": "lof do not cause",
             "confidence_level": 3, "publications": "12"},
            {"gene_data": {"hgnc_symbol": "AMBER"}, "mode_of_pathogenicity": "do not cause loss of function", "confidence_level": "2"},
            {"gene_data": {}, "mode_of_pathogenicity": "x", "confidence_level": "3"}]},
        {"id": 2, "name": 7, "version": None, "genes": "not a list"}]}
    (raw / "panelapp_all_panels.json").write_text(json.dumps(panels, sort_keys=True) + "\n")
    (raw / "gofcards").mkdir(exist_ok=True)
    write_workbook(raw / "gofcards/gofcards_data_download.xlsx", picked)
    (raw / "clinvar_vcv").mkdir(exist_ok=True)
    (raw / CLINVAR_RAW).write_bytes(gzip.compress(b"<ClinVarVariationRelease/>"))


def fresh_manifest(cache_dir, raw, urls=None, clinvar=True, checked=None):
    urls = urls or {}
    checked = checked or now()
    sources = {}
    for name, url, filename in [("g2p_ddg2p", G2P_URL, "AllG2P.csv"),
                                ("clingen_dosage_gene_grch38", CLINGEN_URL, "ClinGen_gene_curation_list_GRCh38.tsv"),
                                ("gofcards_gof_variants", GOFCARDS_URL, "gofcards/gofcards_data_download.xlsx"),
                                ("panelapp_all_panels", PANELAPP_URL, "panelapp_all_panels.json")]:
        path = raw / filename
        sources[name] = {"name": name, "url": urls.get(name, url), "raw_path": str(path), "last_checked_at_utc": checked,
                         "sha256": sha256(path) if path.exists() else "", "status": "downloaded"}
    if clinvar:
        path = raw / CLINVAR_RAW
        sources["clinvar_vcv_weekly"] = {"name": "clinvar_vcv_weekly", "raw_path": str(path), "sha256": sha256(path),
                                         "md5": "0" * 32, "last_remote_md5_checked_at_utc": checked, "status": "downloaded"}
    cache_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / "nonlof_source_manifest.json").write_text(json.dumps({"sources": sources}, indent=2) + "\n")


def vcv_xml(kind="valid"):
    records = '<VariationArchive VariationID="12" Accession="VCV1"><X/></VariationArchive>' \
              '<VariationArchive VariationID="34"/>'
    body = {
        "valid": f'<?xml version="1.0"?><ClinVarVariationRelease ReleaseDate="2026-10-01">{records}</ClinVarVariationRelease>',
        "namespaced": f'<ClinVarVariationRelease xmlns="urn:x" ReleaseDate="d">{records}</ClinVarVariationRelease>',
        "wrong_root": f"<Other>{records}</Other>",
        "no_id": '<ClinVarVariationRelease><VariationArchive VariationID="x1"/></ClinVarVariationRelease>',
        "no_records": "<ClinVarVariationRelease><Other/></ClinVarVariationRelease>",
        "broken": f"<ClinVarVariationRelease>{records}<VariationArchive",
        "entity": "<ClinVarVariationRelease>&undefined;</ClinVarVariationRelease>",
    }[kind]
    return body.encode()


class Case:
    def __init__(self, name, prepare, steps):
        self.name, self.prepare, self.steps = name, prepare, steps


def run_case(args, case, root, env):
    folder = root / case.name
    arms = {"python": folder / "python", "native": folder / "native"}
    arms["python"].mkdir(parents=True)
    case.prepare(arms["python"])
    # Identical input bytes and times in both copies; only the copy's own
    # path inside the prepared JSON files differs.
    shutil.copytree(arms["python"], arms["native"], symlinks=True)
    for path in arms["native"].rglob("*.json*"):
        if path.is_file() and path.name.endswith((".json", ".json.gz")):
            data = path.read_bytes()
            compressed = path.name.endswith(".gz")
            text = (gzip.decompress(data) if compressed else data).decode()
            if str(arms["python"]) in text:
                text = text.replace(str(arms["python"]), str(arms["native"])).encode()
                stamp = path.stat()
                path.write_bytes(gzip.compress(text, mtime=0) if compressed else text)
                os.utime(path, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))
    results = []
    for step in case.steps:
        pair = {}
        for label, arm in arms.items():
            arguments = step(arm)
            if label == "python":
                command = [sys.executable, args.repo / "scripts/build_gene_nonlof_mechanism_cache.py", *arguments]
            else:
                command = [args.binary, "gene-nonlof-mechanism", "--source-root", args.repo, *arguments]
            pair[label] = outcome(command, arm, env)
            pair[label]["files"] = snapshot(arm)
        python, native = pair["python"], pair["native"]
        keys = ["exit", "progress", "stdout", "files"] + (["category", "message"] if python["exit"] == 1 else [])
        differences = [key for key in keys if python.get(key) != native.get(key)]
        differences += [f"file:{name}" for name in sorted(set(python["files"]) | set(native["files"]))
                        if python["files"].get(name) != native["files"].get(name)]
        results.append({"equal": not differences, "differences": differences,
                        "python": {k: v for k, v in python.items() if k != "files"},
                        "native": {k: v for k, v in native.items() if k != "files"}})
    (folder / "comparison.json").write_text(json.dumps(results, indent=2, default=str) + "\n")
    return {"name": case.name, "same": all(r["equal"] for r in results),
            "steps": [{"exit": r["python"]["exit"], "differences": r["differences"]} for r in results]}


class Fixtures:
    """The inputs shared by every case and the two argument helpers."""

    def __init__(self, repo):
        self.repo = repo
        self.cache, self.picked = gofcards_genes(repo)
        self.schema = repo / "data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.schema.json"
        self.hgnc = repo / "data/hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv"

    def base(self, arm, manifest=True, vcv="valid", inject=False, workbook="plain"):
        raw = arm / "raw"
        write_sources(raw, self.picked, self.repo)
        if workbook != "plain":
            write_workbook(raw / "gofcards/gofcards_data_download.xlsx", self.picked, workbook)
        (arm / "vcv.xml.gz").write_bytes(gzip.compress(vcv_xml(vcv)))
        local = copy.deepcopy(self.cache)
        if inject:
            # Relative to the copy, where both commands run: identical cache bytes.
            local["metadata"]["clinvar"]["source"] = "vcv.xml.gz"
            local["metadata"]["clinvar"]["source_sha256"] = sha256(arm / "vcv.xml.gz")
        (arm / "gofcards.json.gz").write_bytes(gzip.compress(json.dumps(local).encode(), mtime=0))
        shutil.copy2(self.schema, arm / "schema.json")
        if manifest:
            fresh_manifest(arm / "cache", raw)

    def common(self, arm, *extra):
        return ["--cache-dir", arm / "cache", "--shared-raw-dir", arm / "raw", "--gofcards-exact-variants",
                arm / "gofcards.json.gz", "--hgnc-table", self.hgnc, "--output-schema", arm / "schema.json",
                "--output-json", arm / "out/nonlof.json.gz", *extra]


def local_cases(f):
    base, common = f.base, f.common
    vcv = lambda a, *extra: common(a, "--clinvar-vcv", a / "vcv.xml.gz", *extra)

    def cache_local(arm):
        base(arm, manifest=False)
        shutil.copytree(arm / "raw", arm / "cache/raw")
        shutil.rmtree(arm / "raw")
        fresh_manifest(arm / "cache", arm / "cache/raw")

    def replaced(arm, relative, body):
        base(arm)
        path = arm / relative
        path.write_bytes(body) if isinstance(body, bytes) else path.write_text(body)

    def schema_change(arm, change):
        base(arm)
        (arm / "schema.json").write_text(change(f.schema.read_text()))

    def stale_lock(arm):
        base(arm)
        (arm / "raw/.build.lock").write_text("pid=1\n")
        os.utime(arm / "raw/.build.lock", (1, 1))

    missing_property = lambda text: json.dumps({**json.loads(text), "required": ["_meta", "MISSING"],
                                                "properties": {"_meta": {"$ref": "#/$defs/meta"},
                                                               "MISSING": {"type": "string"}}})
    deeper = lambda text: text.replace('"minimum": 0', '"minimum": 100000').replace('"^SYMBOL:.+$"', '"^SYMBOL:ZZZ.+$"')
    result = [
        Case("build_then_reuse", base, [common, common]),
        Case("force_rebuild_explicit_vcv", lambda a: base(a, inject=True),
             [lambda a: vcv(a), lambda a: vcv(a, "--force", "--download-tool", "curl", "--timeout", "2",
                                             "--retries", "1", "--proxy-url", "http://127.0.0.1:9")]),
        Case("cache_local_raw", cache_local,
             [lambda a: ["--cache-dir", a / "cache", "--shared-raw-dir", " ", "--gofcards-exact-variants",
                         a / "gofcards.json.gz", "--hgnc-table", f.hgnc, "--output-schema", a / "schema.json"]]),
        Case("explicit_vcv_fetch_only", base, [lambda a: vcv(a, "--fetch-sources-only")]),
        Case("explicit_vcv_build", lambda a: base(a, inject=True), [vcv, vcv]),
        Case("explicit_vcv_not_injected", base, [vcv]),
        Case("fetch_only_keeps_build_inputs", base, [common, lambda a: common(a, "--fetch-sources-only"), common]),
        Case("validate_only", base, [common, lambda a: common(a, "--validate-only", a / "out/nonlof.json.gz")]),
        Case("validate_only_missing", base, [lambda a: common(a, "--validate-only", a / "absent.json")]),
        Case("dated_workbook", lambda a: base(a, workbook="dates"), [common]),
        Case("unnamed_header", lambda a: base(a, workbook="unnamed"), [common]),
        Case("wrong_sheet", lambda a: base(a, workbook="wrong_sheet"), [common]),
        Case("not_a_workbook", lambda a: (replaced(a, "raw/gofcards/gofcards_data_download.xlsx", "plain"),
                                          fresh_manifest(a / "cache", a / "raw")), [common]),
        Case("missing_gofcards_cache", lambda a: (base(a), (a / "gofcards.json.gz").unlink()), [common]),
        Case("bad_gofcards_cache", lambda a: replaced(a, "gofcards.json.gz", gzip.compress(
            b'{"genes": {}, "metadata": {"source": "x"}}')), [common]),
        Case("lock_present", lambda a: replaced(a, "raw/.build.lock", "pid=1\n"), [common]),
        Case("stale_lock", stale_lock, [common]),
        Case("clinvar_only_refresh", base, [common, lambda a: common(a, "--clinvar-only-refresh")]),
        Case("clinvar_only_refresh_changed_url", base,
             [lambda a: common(a, "--clinvar-only-refresh", "--g2p-ddg2p-url", "https://example.invalid/g2p.csv")]),
        Case("schema_failure", lambda a: schema_change(a, missing_property), [common]),
        Case("schema_failure_deep", lambda a: schema_change(a, deeper), [common]),
        Case("schema_missing", lambda a: replaced(a, "schema.json", "{}"), [common]),
        Case("stale_sources_unreachable_curl", lambda a: base(a, manifest=False),
             [lambda a: vcv(a, "--fetch-sources-only", "--download-tool", "curl", "--timeout", "2", "--retries", "1",
                            "--proxy-url", "http://127.0.0.1:9")]),
        Case("vcv_empty", lambda a: replaced(a, "vcv.xml.gz", b""), [lambda a: vcv(a, "--fetch-sources-only")]),
    ]
    for kind in ("namespaced", "wrong_root", "no_id", "no_records", "broken", "entity"):
        result.append(Case(f"vcv_{kind}", lambda a, kind=kind: base(a, vcv=kind),
                           [lambda a: vcv(a, "--fetch-sources-only")]))
    for name, extra in [("bad_int", ["--timeout", "x"]), ("bad_choice", ["--download-tool", "wget"]),
                        ("bad_stars", ["--clinvar-min-review-stars", "5"]), ("bad_url", ["--g2p-ddg2p-url", "ftp://x"]),
                        ("unknown_option", ["--nope"]), ("bad_float", ["--stale-lock-hours", "soon"])]:
        result.append(Case(f"argument_{name}", base, [lambda a, extra=extra: common(a, *extra)]))
    return result


def served_cases(f, server, root):
    """Selected URLs served locally, through urllib's route in the script."""
    write_sources(root / "src", f.picked, f.repo, strict=True)
    shutil.copy2(root / "src/AllG2P.csv", root / "g2p.csv")
    shutil.copy2(root / "src/ClinGen_gene_curation_list_GRCh38.tsv", root / "clingen.tsv")
    (root / "g2p_bad.csv").write_text("a,b\n1,2\n")
    gene = {"gene_data": {"gene_symbol": f.picked[0][0]}, "mode_of_pathogenicity": "LOF do not cause",
            "confidence_level": "3", "publications": ["7"]}
    pages = {"panelapp": {"results": [{"id": 2}], "next": server.url("panelapp/page2/")},
             "panelapp/page2": {"results": [{"id": 3}], "next": None},
             "panelapp/2": {"id": 2, "name": "P2", "version": "1", "genes": [gene]},
             "panelapp/3": {"id": 3, "name": "P3", "version": "1", "genes": []},
             "badpanel": {"results": "x"}}
    for path, body in pages.items():
        (root / path).mkdir(parents=True, exist_ok=True)
        (root / path / "index.html").write_text(json.dumps(body))
    selected = ["--g2p-ddg2p-url", server.url("g2p.csv"), "--clingen-dosage-grch38-url", server.url("clingen.tsv"),
                "--panelapp-panels-url", server.url("panelapp/")]
    vcv = lambda a, *extra: f.common(a, "--clinvar-vcv", a / "vcv.xml.gz", *extra)
    fresh = lambda a: f.base(a, manifest=False)
    return [
        Case("served_sources", fresh, [lambda a: vcv(a, "--fetch-sources-only", *selected)]),
        Case("served_then_build_then_reuse", lambda a: f.base(a, manifest=False, inject=True),
             [lambda a: vcv(a, *selected), lambda a: vcv(a, *selected)]),
        Case("served_404", fresh, [lambda a: vcv(a, "--fetch-sources-only", "--retries", "1",
                                                 "--g2p-ddg2p-url", server.url("absent.csv"))]),
        Case("served_incompatible", fresh, [lambda a: vcv(a, "--fetch-sources-only",
                                                          "--g2p-ddg2p-url", server.url("g2p_bad.csv"))]),
        Case("served_panelapp_schema", fresh, [lambda a: vcv(a, "--fetch-sources-only", "--retries", "1",
                                                             "--panelapp-panels-url", server.url("badpanel/"))]),
        Case("served_panelapp_cap", fresh, [lambda a: vcv(a, "--fetch-sources-only", "--max-panelapp-panels", "1",
                                                          "--panelapp-panels-url", server.url("panelapp/"))]),
    ]


def check(args):
    args.work_root.mkdir(parents=True)
    env = dict(os.environ, no_proxy="127.0.0.1,localhost", NO_PROXY="127.0.0.1,localhost", PYTHONHASHSEED="0")
    for key in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY", "all_proxy", "ALL_PROXY"):
        env.pop(key, None)
    fixtures = Fixtures(args.repo)
    served = args.work_root / "served"
    served.mkdir()
    server = Server(served)
    selected = local_cases(fixtures) + served_cases(fixtures, server, served)

    def one(case):
        row = run_case(args, case, args.work_root / "cases", env)
        print(json.dumps(row), file=sys.stderr, flush=True)
        return row
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            rows = list(pool.map(one, selected))
    finally:
        server.close()
    failed = [row["name"] for row in rows if not row["same"]]
    print(json.dumps({"cases": len(rows), "steps": sum(len(r["steps"]) for r in rows), "failed": failed}))
    return 1 if failed else 0


def real(args):
    args.work_root.mkdir(parents=True)
    env = dict(os.environ, PYTHONHASHSEED="0")
    repo = args.repo
    rows = []
    for label in ("python", "native"):
        arm = args.work_root / label
        raw = arm / "raw"
        (raw / "gofcards").mkdir(parents=True)
        (raw / "clinvar_vcv").mkdir(parents=True)
        for name in ("AllG2P.csv", "ClinGen_gene_curation_list_GRCh38.tsv", "panelapp_all_panels.json",
                     "gofcards/gofcards_data_download.xlsx"):
            (raw / name).symlink_to((args.raw_inputs / name).resolve())
        (raw / CLINVAR_RAW).write_bytes(gzip.compress(b"<ClinVarVariationRelease/>"))
        fresh_manifest(arm / "cache", raw)
    for label in ("python", "native"):
        arm = args.work_root / label
        arguments = ["--cache-dir", arm / "cache", "--shared-raw-dir", arm / "raw", "--gofcards-exact-variants",
                     args.gofcards_cache, "--hgnc-table", repo / "data/hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv",
                     "--output-schema", repo / "data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.schema.json"]
        if label == "python":
            command = [sys.executable, repo / "scripts/build_gene_nonlof_mechanism_cache.py", *arguments]
        else:
            command = [args.binary, "gene-nonlof-mechanism", "--source-root", repo, *arguments]
        start = time.time()
        value = outcome(command, arm, env)
        value["seconds"] = round(time.time() - start, 1)
        value["files"] = snapshot(arm, arm / "cache")
        rows.append(value)
    python, native = rows
    same = all(python.get(k) == native.get(k) for k in ("exit", "progress", "stdout", "files"))
    output = python["files"].get("gene_nonlof_mechanism_curated_assertions.json.gz") or {}
    report = {"same": same, "python_seconds": python["seconds"], "native_seconds": native["seconds"],
              "genes": len(output) - 1 if isinstance(output, dict) else None,
              "differences": sorted(k for k in set(python) | set(native) if k != "seconds" and python.get(k) != native.get(k))}
    (args.work_root / "comparison.json").write_text(json.dumps({"python": python, "native": native}, indent=1, default=str))
    print(json.dumps(report))
    return 0 if same else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("check", "real"):
        command = commands.add_parser(name)
        command.add_argument("--binary", type=Path, required=True)
        command.add_argument("--repo", type=Path, required=True)
        command.add_argument("--work-root", type=Path, required=True)
        command.add_argument("--jobs", type=int, default=1)
        if name == "real":
            command.add_argument("--raw-inputs", type=Path, required=True)
            command.add_argument("--gofcards-cache", type=Path, required=True)
    args = parser.parse_args()
    args.repo = args.repo.resolve()
    sys.exit(check(args) if args.command == "check" else real(args))


if __name__ == "__main__":
    main()
