#!/usr/bin/env python3
"""Compare `priva gofcards-exact-gof` with build_gofcards_exact_gof_cache.py.

  python gofcards_exact_gof_compare.py check --binary PRIVA --repo CHECKOUT --work-root NEW_DIR [--jobs N]
  python gofcards_exact_gof_compare.py real --binary PRIVA --repo CHECKOUT --work-root NEW_DIR \
      --workbook XLSX --summary-cache JSONL --hg19-fasta FA --hg38-fasta FA --chain CHAIN \
      --vep-cache-dir DIR --vep-cache-version N --hgnc-table TSV --mechanism-review-tsv TSV

Each case writes the same inputs into two folders and runs the builder in
both: the script through a small driver that points its summary endpoint at a
local server, and the native command with --summary-endpoint. Stand-in
bcftools, vep and CrossMap programs read a control file, so no case reaches a
publisher or a real annotation. Exit status, failure category and message, the
`[gofcards]` log lines, the requests the server saw and every file left in the
folder are compared. gzip files are compared by decompressed text and by their
bytes apart from the time in the header.

`real` serves an existing workbook from this machine, seeds each folder with an
existing summary cache so that no allele is requested, and runs the real
bcftools, VEP and CrossMap with both. VCF files are compared by hash only.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import copy
import datetime
import gzip
import hashlib
import http.server
import json
import os
import random
import re
import shutil
import stat
import subprocess
import sys
import threading
import time
from pathlib import Path
from urllib.parse import parse_qsl, urlsplit

SCRIPT_NAME = "build_gofcards_exact_gof_cache"
COLUMNS = ["GeneSymbol", "Transcript", "Chr", "hg19Start", "hg19End", "Ref", "Alt", "Function",
           "Pathways Proteins Involved", "Disorder Involved", "PMID", "Animal Model", "Cell Model",
           "Pscore", "Order Numbe"]
REVIEW_COLUMNS = ["source_order", "source_gene", "allele_key", "source_mechanism", "reviewed_mechanism",
                  "review_status", "gof_eligibility", "reason_code", "evidence_summary", "pmids",
                  "reviewer", "reviewed_at", "review_version"]
LENGTH = 3000

DRIVER = '''
import importlib, sys
sys.path.insert(0, sys.argv[1])
module = importlib.import_module("%s")
module.SUMMARY_ENDPOINT = sys.argv[2]
sys.argv = ["builder"] + sys.argv[3:]
sys.exit(module.main())
''' % SCRIPT_NAME

TOOL = r'''#!%s
import json, os, sys
from pathlib import Path
name = Path(sys.argv[0]).name
args = sys.argv[1:]
control = json.loads(Path("control.json").read_text())
with open("tool_calls.log", "a") as calls:
    calls.write(json.dumps([name, args]) + "\n")
print(name, "says hello")
print(name, "complains", file=sys.stderr)
if name in control.get("fail", []):
    sys.exit(7)
if name in control.get("silent", []):
    sys.exit(0)

def value(option):
    return args[args.index(option) + 1]

def records(path):
    for line in Path(path).read_text().splitlines():
        if not line.startswith("#"):
            yield line.split("\t")

if name == "bcftools":
    output, source = value("-o"), args[-1]
    text = Path(source).read_text()
    Path(output).write_text(text)
    for extra in control.get("extra_lines", []):
        with open(output, "a") as handle:
            handle.write(extra + "\n")
elif name == "CrossMap":
    chain, source, fasta, output = args[1:5]
    lines, unmap = [], []
    for line in Path(source).read_text().splitlines():
        if line.startswith("#"):
            lines.append(line)
            continue
        fields = line.split("\t")
        if fields[2] in control.get("unmapped", []):
            unmap.append(line)
            continue
        fields[0] = control.get("contig_prefix", "chr") + fields[0].removeprefix("chr")
        fields[1] = str(int(fields[1]) + control.get("shift", 2))
        lines.append("\t".join(fields))
    Path(output).write_text("\n".join(lines) + "\n")
    if unmap:
        Path(output + ".unmap").write_text("#header\n" + "\n".join(unmap) + "\n")
elif name == "vep":
    assembly = value("--assembly")
    rows = control["vep"].get(assembly, {})
    header = ("#Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
              "HGVSc", "HGVSp", "MANE_SELECT", "MANE_PLUS_CLINICAL", "CANONICAL", "SYMBOL")
    out = ["## ENSEMBL VARIANT EFFECT PREDICTOR", "## a comment", "\t".join(header)]
    for fields in records(value("--input_file")):
        for row in rows.get(fields[2], []):
            out.append("\t".join([fields[2], f"{fields[0]}:{fields[1]}", fields[4], "ENSG1", row["feature"],
                                  "Transcript", row.get("consequence", "missense_variant"), row.get("hgvsc", "-"),
                                  row.get("hgvsp", "-"), row.get("mane", "-"), "-", row.get("canonical", "-"),
                                  row.get("symbol", "-")]))
    if not control.get("vep_no_output"):
        Path(value("--output_file")).write_text("\n".join(out) + "\n")
''' % sys.executable


def sequence(chrom):
    generator = random.Random("gofcards" + chrom)
    return "".join(generator.choice("ACGT") for _ in range(LENGTH))


def bases(chrom, start, length=1):
    return sequence(chrom)[start - 1:start - 1 + length]


def write_fasta(path, records, index=True):
    offset, fai, lines = 0, [], []
    for name, text in records:
        header = f">{name}\n"
        lines.append(header)
        offset += len(header)
        fai.append(f"{name}\t{len(text)}\t{offset}\t60\t61\n")
        wrapped = "".join(text[at:at + 60] + "\n" for at in range(0, len(text), 60))
        lines.append(wrapped)
        offset += len(wrapped)
    path.write_text("".join(lines))
    if index:
        Path(f"{path}.fai").write_text("".join(fai))


def row(gene, chrom, start, ref=None, alt=None, length=1, vep="same", **overrides):
    """One workbook row at a real position; `ref` defaults to the genome there."""
    if ref is None:
        ref = bases(chrom, start, length)
    if alt is None:
        alt = {"A": "C", "C": "G", "G": "T", "T": "A"}.get(ref[0], "T")
    values = {"GeneSymbol": gene, "Transcript": "NM_000001.1", "Chr": chrom, "hg19Start": start,
              "hg19End": start + max(len(ref), 1) - 1 if start is not None else None, "Ref": ref, "Alt": alt, "Function": "gain of function",
              "Pathways Proteins Involved": "pathway A", "Disorder Involved": "disease A", "PMID": 1000 + (start or 0),
              "Animal Model": "mouse", "Cell Model": "_", "Pscore": 0.5, "Order Numbe": start or 0}
    values.update(overrides)
    values["vep"] = gene if vep == "same" else vep
    return values


def base_spec():
    return {
        "rows": [
            row("GENE1", "1", 100),
            row("GENE1", "1", 100, PMID=4242, **{"Disorder Involved": "disease B"}),
            row("GENE1", "1", 200, alt="-", length=3),
            row("GENE2", "2", 300, ref="-", alt="TTA"),
            row("gene3", "X", 400, length=2, alt="GG"),
            row("OLDNAME", "1", 500, vep="GENE4"),
            row("GENE5", "1", 600, ref="N"),
            row("GENE6", "Z", 700, ref="A", alt="C"),
            row("GENE7", "1", None, ref="A", alt="C", hg19End=None),
            row("GENE8", "2", 800, vep="OTHER"),
            row("GENE9", "2", 900, vep=None),
            row("GENE10", "2", 1000, vep=["OTHER", ""]),
        ],
        "columns": list(COLUMNS),
        "vep_extra": {},
        "control": {"unmapped": []},
        "summary": {},
        "reviews": [],
        "hgnc": [("HGNC:1", "GENE1", "", ""), ("HGNC:2", "GENE2", "", ""), ("HGNC:3", "GENE3", "", ""),
                 ("HGNC:4", "GENE4", "OLDNAME|OTHERNAME", ""), ("HGNC:5", "GENE5", "", ""),
                 ("HGNC:6", "GENE6", "", ""), ("HGNC:7", "GENE7", "", ""), ("HGNC:8", "GENE8", "", ""),
                 ("HGNC:9", "GENE9", "", ""), ("HGNC:10", "GENE10", "", ""), ("", "", "", ""),
                 ("HGNC:11", "OTHER", "", "PREV1,PREV2")],
        "args": [],
        "drop_args": [],
        "workbook": None,
        "server": {},
        "hg38_prefix": "chr",
        "fai": True,
        "hide": [],
        "setup": None,
    }


def vcf_ids(rows):
    ids, seen = {}, {}
    for at, values in enumerate(rows):
        key = tuple(str(values.get(name)) for name in ("Chr", "hg19Start", "hg19End", "Ref", "Alt"))
        ids[at] = seen.setdefault(key, f"gof{len(seen):06d}")
    return ids


def transcript_rows(values, assembly):
    gene = values["vep"]
    if gene is None:
        return []
    genes = gene if isinstance(gene, list) else [gene]
    position = values["hg19Start"]
    out = []
    for number, symbol in enumerate(genes):
        version = 1 if assembly == "GRCh37" else 2
        out.append({"feature": f"ENST{number + 1:05d}.{version}", "symbol": symbol,
                    "hgvsc": f"ENST{number + 1:05d}.{version}:c.{position}A>G",
                    "hgvsp": f"ENSP{number + 1:05d}.{version}:p.Lys{position}Glu",
                    "canonical": "YES" if number == 0 else "-", "mane": "NM_000001.1" if number == 0 else "-"})
    return out


def write_workbook(path, spec):
    import openpyxl
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.append(spec["columns"])
    for values in spec["rows"]:
        sheet.append([values.get(name) for name in spec["columns"]])
    for name, rows in spec.get("extra_sheets", {}).items():
        extra = book.create_sheet(name)
        for values in rows:
            extra.append(values)
    book.save(path)


def materialize(folder, spec, base_url, tools):
    folder.mkdir(parents=True)
    (folder / "server").mkdir()
    workbook = folder / "server" / "public.xlsx"
    if spec["workbook"] is None:
        write_workbook(workbook, spec)
    else:
        workbook.write_bytes(spec["workbook"])
    control = dict(spec["control"])
    ids = vcf_ids(spec["rows"])
    control["vep"] = {}
    for assembly in ("GRCh37", "GRCh38"):
        table = {}
        for at, values in enumerate(spec["rows"]):
            if isinstance(values.get("hg19Start"), int):
                table.setdefault(ids[at], transcript_rows(values, assembly))
        table.update(spec["vep_extra"].get(assembly, {}))
        control["vep"][assembly] = table
    (folder / "control.json").write_text(json.dumps(control))
    (folder / "server" / "control.json").write_text(json.dumps({**spec["server"], "summary": spec["summary"]}))
    records19 = [(f"chr{name}", sequence(name)) for name in ("1", "2", "X")]
    records38 = [(f"{spec['hg38_prefix']}{name}", "GG" + sequence(name)) for name in ("1", "2", "X")]
    write_fasta(folder / "hg19.fa", records19, spec["fai"])
    write_fasta(folder / "hg38.fa", records38, spec["fai"])
    (folder / "chain.txt").write_text("chain\n")
    (folder / "vep_cache").mkdir()
    lines = ["hgnc_id\tsymbol\talias_symbol\tprev_symbol"]
    lines += ["\t".join(item) for item in [(a, b, f'"{c}"' if c else "", d) for a, b, c, d in spec["hgnc"]]]
    (folder / "hgnc.tsv").write_text("\n".join(lines) + "\n")
    if spec["reviews"] is not None:
        review_lines = ["\t".join(spec.get("review_columns", REVIEW_COLUMNS))]
        for values in spec["reviews"]:
            review_lines.append("\t".join(values.get(name, "") for name in spec.get("review_columns", REVIEW_COLUMNS)))
        (folder / "reviews.tsv").write_text("\n".join(review_lines) + "\n")
    if spec["setup"]:
        spec["setup"](folder)
    arguments = {"--workdir": "work", "--out-json": "out/cache.json.gz", "--hg19-fasta": "hg19.fa",
                 "--hg38-fasta": "hg38.fa", "--chain": "chain.txt", "--vep-cache-dir": "vep_cache",
                 "--vep-cache-version": "116", "--hgnc-table": "hgnc.tsv",
                 "--mechanism-review-tsv": "reviews.tsv", "--public-excel-url": f"{base_url}/public.xlsx",
                 "--api-timeout": "20", "--stats-json": "stats.json"}
    command = []
    for option, value in arguments.items():
        if option not in spec["drop_args"]:
            command += [option, value]
    return command + spec["args"]


class Handler(http.server.BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    folders: dict = {}
    hits: dict = {}
    lock = threading.Lock()

    def log_message(self, *_):
        pass

    def reply(self, status, body, headers=()):
        self.send_response(status)
        self.send_header("Content-Length", str(len(body)))
        for name, value in headers:
            self.send_header(name, value)
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        parts = urlsplit(self.path)
        _, token, resource = parts.path.split("/", 2)
        folder = self.folders[token]
        control = json.loads((folder / "server" / "control.json").read_text())
        names = ("User-Agent", "Accept", "Referer", "Origin", "Admin-Work-Version", "Content-Type")
        with open(folder / "requests.log", "a") as log:
            log.write(json.dumps([resource, sorted(parse_qsl(parts.query, keep_blank_values=True)), parts.query,
                                  [self.headers.get(name) for name in names]]) + "\n")
        if resource == "public.xlsx":
            status = control.get("workbook_status", 200)
            if status != 200:
                self.reply(status, b"no")
            elif control.get("redirect"):
                self.reply(302, b"", [("Location", f"/{token}/moved.xlsx")])
            else:
                self.send_workbook(folder)
            return
        if resource == "moved.xlsx":
            self.send_workbook(folder)
            return
        query = dict(parse_qsl(parts.query, keep_blank_values=True))
        key = f"{query.get('chr')}:{query.get('start')}"
        answer = control["summary"].get(key, control["summary"].get("*", {}))
        with self.lock:
            count = self.hits[(token, key)] = self.hits.get((token, key), 0) + 1
        if count <= answer.get("fail_first", 0):
            self.reply(answer.get("fail_status", 503), b"busy")
            return
        if "raw" in answer:
            self.reply(answer.get("status", 200), answer["raw"].encode())
            return
        data = answer.get("data", {"Accession": f"{query.get('start')}5", "CLNSIG": "x"})
        self.reply(answer.get("status", 200), json.dumps({"code": 0, "data": data}).encode())

    def send_workbook(self, folder):
        self.reply(200, (folder / "server" / "public.xlsx").read_bytes(),
                   [("Last-Modified", "Fri, 27 Sep 2024 03:01:42 GMT")])


class Server(http.server.ThreadingHTTPServer):
    daemon_threads = True


def start_server():
    server = Server(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server


CASES = {}


def case(function):
    CASES[function.__name__] = function
    return function


def set_column(spec, name, values):
    for item, value in zip(spec["rows"], values):
        item[name] = value


def date_rows():
    return [datetime.datetime(2020, 1, 2 + at, 3, 4, 5) for at in range(12)]


@case
def base(spec): pass
@case
def out_tsv_option(spec):
    spec["drop_args"] = ["--out-json"]
    spec["args"] = ["--out-tsv", "out/other.json.gz"]
@case
def output_without_gz_suffix(spec):
    spec["drop_args"] = ["--out-json"]
    spec["args"] = ["--out-json", "plain_cache.json"]
@case
def no_stats(spec): spec["drop_args"] = ["--stats-json"]
@case
def no_vep_version(spec): spec["drop_args"] = ["--vep-cache-version"]
@case
def unnormalized_paths(spec):
    spec["drop_args"] = ["--workdir", "--vep-cache-dir"]
    spec["args"] = ["--workdir", "./work//deep/", "--vep-cache-dir", "./vep_cache/"]
@case
def stats_parent_missing(spec):
    spec["drop_args"] = ["--stats-json"]
    spec["args"] = ["--stats-json", "missing/stats.json"]
@case
def redirected_workbook(spec): spec["server"] = {"redirect": True}
@case
def workbook_404(spec): spec["server"] = {"workbook_status": 404}
@case
def workbook_500(spec): spec["server"] = {"workbook_status": 500}
@case
def workbook_not_zip(spec): spec["workbook"] = b"this is not a workbook"
@case
def workbook_empty_file(spec): spec["workbook"] = b""
@case
def workbook_empty_sheet(spec):
    spec["columns"], spec["rows"] = [], []
@case
def workbook_header_only(spec): spec["rows"] = []
@case
def workbook_second_sheet_ignored(spec):
    spec["extra_sheets"] = {"Other": [["Chr"], ["9"]]}
@case
def header_spelling(spec):
    renames = {"Chr": " CHR ", "hg19Start": "HG19START", "Ref": "ref", "PMID": "pmid "}
    for values in spec["rows"]:
        for old, new in renames.items():
            values[new] = values.pop(old)
    spec["columns"] = [renames.get(name, name) for name in spec["columns"]]
@case
def extra_and_unnamed_columns(spec):
    spec["columns"] = spec["columns"] + ["Extra", None, "Extra"]
    set_column(spec, "Extra", ["x"] * 12)
@case
def missing_columns(spec):
    spec["columns"] = [name for name in spec["columns"] if name not in ("Function", "Order Numbe", "Transcript")]
@case
def missing_gene_column(spec):
    spec["columns"] = [name for name in spec["columns"] if name != "GeneSymbol"]
@case
def missing_coordinate_columns(spec):
    spec["columns"] = [name for name in spec["columns"] if name not in ("Chr", "hg19Start")]
@case
def pmid_numbers_with_blank(spec):
    set_column(spec, "PMID", [111, None] + list(range(10)))
@case
def pmid_numbers_and_text(spec):
    set_column(spec, "PMID", [111, "222;333", None, 4.5] + list(range(8)))
@case
def pmid_numeric_text(spec):
    set_column(spec, "PMID", ["111", 222, " 333 ", "4e2", "+5", ".5", "6.", None, "7", "8", "9", "10"])
@case
def pmid_na_spellings(spec):
    set_column(spec, "PMID", ["NA", "N/A", "null", "None", "nan", "NaN", "n/a", "<NA>", "#N/A", "x", "-", "_"])
@case
def order_with_blank(spec):
    set_column(spec, "Order Numbe", [1, None] + list(range(2, 12)))
@case
def pscore_mixed_floats(spec):
    set_column(spec, "Pscore", [0.1, 1, 2.0, 1e20, 1e-7, None, 3.25, 7, 8, 9, 10, 11])
@case
def pscore_all_whole_floats(spec):
    set_column(spec, "Pscore", [float(n) for n in range(12)])
@case
def boolean_column(spec):
    set_column(spec, "Animal Model", [True, False] * 6)
@case
def boolean_column_with_blank(spec):
    set_column(spec, "Animal Model", [True, None] * 6)
@case
def boolean_and_numbers(spec):
    set_column(spec, "Animal Model", [True, 5] * 6)
@case
def boolean_and_text(spec):
    set_column(spec, "Animal Model", [True, "mouse"] * 6)
@case
def boolean_text(spec):
    set_column(spec, "Animal Model", ["true", "FALSE", "True", None] * 3)
@case
def date_column(spec):
    set_column(spec, "Function", date_rows())
@case
def date_column_with_blank(spec):
    set_column(spec, "Function", [None] + date_rows()[1:])
@case
def date_and_text_column(spec):
    set_column(spec, "Function", ["text"] + date_rows()[1:])
@case
def date_in_pmid_at_midnight(spec):
    set_column(spec, "PMID", [datetime.datetime(2021, 5, 6)] * 12)
@case
def numeric_chromosomes_only(spec):
    spec["rows"] = [row("GENE1", "1", 100), row("GENE2", "2", 300)]
    set_column(spec, "Chr", [1, 2])
@case
def numeric_chromosomes_with_blank(spec):
    spec["rows"] = [row("GENE1", "1", 100), row("GENE2", "2", 300), row("GENE5", "2", 310)]
    set_column(spec, "Chr", [1, 2, None])
@case
def mixed_chromosomes(spec):
    set_column(spec, "Chr", [1, 1, 1, 2, "X", "chr1", "1", "Z", "1", 2, "CHR2", 2])
@case
def float_coordinates(spec):
    for item in spec["rows"]:
        if isinstance(item["hg19Start"], int):
            item["hg19Start"] = float(item["hg19Start"])
    spec["rows"][0]["hg19End"] = 100.5
@case
def text_coordinates(spec):
    for item in spec["rows"]:
        if isinstance(item["hg19Start"], int):
            item["hg19Start"] = f"{item['hg19Start']}.0"
@case
def bad_coordinate_text(spec):
    spec["rows"][0]["hg19Start"] = "abc"
    spec["rows"][2]["hg19Start"] = "-5"
    spec["rows"][3]["hg19Start"] = "²"
@case
def lower_case_alleles(spec):
    spec["rows"][0]["Ref"] = spec["rows"][0]["Ref"].lower()
    spec["rows"][1]["Ref"] = spec["rows"][1]["Ref"].lower()
    spec["rows"][4]["Alt"] = "gg"
@case
def placeholder_alleles(spec):
    spec["rows"] = [row("GENE1", "1", 100, alt="."), row("GENE2", "1", 200, ref="_", alt="AC"),
                    row("GENE5", "1", 220, ref="-", alt="-"), row("GENE6", "1", 240, ref=None, alt=None, Ref=None, Alt=None),
                    row("GENE8", "1", 260)]
@case
def text_with_quotes_and_tabs(spec):
    spec["rows"][7]["Chr"] = 'Z"\tQ'
    spec["rows"][0]["Disorder Involved"] = 'a "quoted" disease\twith tab'
    spec["rows"][1]["GeneSymbol"] = "GENE1"
    spec["rows"][2]["Function"] = "café é"
@case
def deletion_at_first_base(spec):
    spec["rows"].append(row("GENE1", "1", 1, alt="-"))
@case
def start_zero(spec):
    spec["rows"].append(row("GENE1", "1", 0, ref="A", alt="C", hg19End=0))
@case
def insertion_past_contig_end(spec):
    spec["rows"].append(row("GENE1", "1", LENGTH + 500, ref="-", alt="AC", vep=None))
@case
def reference_past_contig_end(spec):
    spec["rows"].append(row("GENE1", "1", LENGTH - 1, ref="AAAAAA", alt="C", vep=None))
@case
def no_rejects(spec):
    spec["rows"] = [item for item in spec["rows"] if item["GeneSymbol"] not in ("GENE5", "GENE6", "GENE7")]
@case
def every_allele_rejected(spec):
    spec["rows"] = [row("GENE5", "1", 600, ref="N"), row("GENE6", "Z", 700, ref="A", alt="C")]
@case
def no_curated_gene(spec):
    for item in spec["rows"]:
        item["GeneSymbol"] = None
        item["vep"] = "GENE1"
@case
def multiple_source_genes(spec):
    spec["rows"][1]["GeneSymbol"] = "OTHER"
    spec["rows"][1]["Transcript"] = "NM_000002.1"
@case
def alleles_sharing_identifier(spec):
    spec["rows"].append(row("GENE1", "1", 100, hg19End=101))
    spec["rows"].append(row("GENE9", "2", 900, hg19End=950, vep="GENE9"))
@case
def summary_cache_resume(spec):
    def setup(folder):
        work = folder / "work"
        work.mkdir()
        (work / "gofcards_summary_cache.jsonl").write_text(
            '{"allele_key": "SNV|2|800|800|%s|%s", "data": {"Accession": 77}}\n\n'
            '  {"allele_key": "unused", "data": null}  \n' % (bases("2", 800), row("X", "2", 800)["Alt"]))
    spec["setup"] = setup
@case
def summary_cache_bad_line(spec):
    def setup(folder):
        (folder / "work").mkdir()
        (folder / "work" / "gofcards_summary_cache.jsonl").write_text("not json\n")
    spec["setup"] = setup
@case
def summary_cache_no_key(spec):
    def setup(folder):
        (folder / "work").mkdir()
        (folder / "work" / "gofcards_summary_cache.jsonl").write_text('{"data": {}}\n')
    spec["setup"] = setup
@case
def summary_cache_not_dict(spec):
    def setup(folder):
        (folder / "work").mkdir()
        (folder / "work" / "gofcards_summary_cache.jsonl").write_text('{"allele_key": "x", "data": "text"}\n')
    spec["setup"] = setup
@case
def summary_retry_then_ok(spec):
    spec["rows"] = spec["rows"][:2]
    spec["summary"] = {"1:100": {"fail_first": 1}}
@case
def summary_fails_four_times(spec):
    spec["rows"] = spec["rows"][:1]
    spec["summary"] = {"*": {"fail_first": 9}}
@case
def summary_404_four_times(spec):
    spec["rows"] = spec["rows"][:1]
    spec["summary"] = {"*": {"fail_first": 9, "fail_status": 404}}
@case
def summary_not_json(spec):
    spec["rows"] = spec["rows"][:1]
    spec["summary"] = {"*": {"raw": "<html>"}}
@case
def summary_list_body(spec):
    spec["rows"] = spec["rows"][:1]
    spec["summary"] = {"*": {"raw": "[1]"}}
@case
def summary_null_body(spec):
    spec["summary"] = {"*": {"raw": "null"}}
@case
def summary_data_text(spec):
    spec["summary"] = {"*": {"raw": '{"data": "text"}'}}
@case
def summary_data_null(spec):
    spec["summary"] = {"*": {"raw": '{"data": null}'}}
@case
def accession_kinds(spec):
    spec["summary"] = {"1:100": {"data": {"Accession": "12345"}}, "1:200": {"data": {"Accession": 7.5}},
                       "2:300": {"data": {"Accession": 0}}, "X:400": {"data": {"Accession": "nan"}},
                       "1:500": {"data": {"Accession": True}}, "2:800": {"data": {"Accession": [1, 2]}},
                       "2:900": {"data": {"Accession": " 9 "}}}
def recover(spec, position, **extra):
    spec["rows"] = [row("GENE1", "1", 100), row("GENE2", "1", 200, vep=None)]
    spec["control"]["unmapped"] = ["gof000001"]
    spec["summary"] = {"1:200": {"data": {"hg38_start": position, **extra}}}
@case
def endpoint_position_accepted(spec): recover(spec, 202)
@case
def endpoint_position_float(spec): recover(spec, 202.0)
@case
def endpoint_position_text(spec): recover(spec, " 202.0 ")
@case
def endpoint_position_rejected(spec): recover(spec, 1234)
@case
def endpoint_position_past_end(spec): recover(spec, 99999)
@case
def endpoint_position_bad(spec): recover(spec, "unknown")
@case
def endpoint_position_absent(spec):
    recover(spec, 202)
    spec["summary"] = {}
@case
def endpoint_position_zero_kept_out(spec): recover(spec, "0")
@case
def endpoint_accepted_with_transcripts(spec):
    recover(spec, 202)
    spec["rows"][1]["vep"] = "GENE2"
@case
def endpoint_contig_missing_on_hg38(spec):
    recover(spec, 202)
    def setup(folder):
        write_fasta(folder / "hg38.fa", [("chr2", "GG" + sequence("2"))])
    spec["setup"] = setup
@case
def unmapped_without_endpoint(spec):
    spec["control"]["unmapped"] = ["gof000000", "gof000003"]
@case
def lifted_contig_without_prefix(spec):
    spec["hg38_prefix"] = ""
    spec["control"]["contig_prefix"] = ""
@case
def everything_unmapped(spec):
    spec["control"]["unmapped"] = [f"gof{at:06d}" for at in range(12)]
@case
def duplicate_id_in_normalized_vcf(spec):
    spec["control"]["extra_lines"] = ["chr1\t105\tgof000000\tA\tC\t.\t.\t."]
@case
def blank_id_and_short_vcf_lines(spec):
    spec["control"]["extra_lines"] = ["chr1\t105\t\tA\tC", "chr1\t106", "chr1\t107\tother\ta\tc\t.\t.\t."]
@case
def non_numeric_position(spec):
    spec["control"]["extra_lines"] = ["chr1\tabc\tgof000000\tA\tC\t.\t.\t."]
@case
def several_transcripts(spec):
    spec["vep_extra"] = {"GRCh38": {"gof000000": [
        {"feature": "ENST1.2", "symbol": "GENE1", "hgvsc": "ENST1.2:c.1A>G", "hgvsp": "ENSP1:p.K1E"},
        {"feature": "ENST1.2", "symbol": "GENE1", "hgvsc": "ENST1.2:c.2A>G", "hgvsp": "ENSP1:p.K1E", "canonical": "YES"},
        {"feature": "ENST1.2", "symbol": "GENE1", "hgvsc": "ENST1.2:c.2A>G", "hgvsp": "ENSP1:p.K1E"},
        {"feature": "ENST2.2", "symbol": "GENE1", "hgvsc": "-", "hgvsp": "ENSP2:p.K1E"},
        {"feature": "ENST3.2", "symbol": "GENE1", "hgvsc": "ENST3.2:n.5A>G", "hgvsp": "-", "mane": "NM_1.1"},
        {"feature": "-", "symbol": "GENE1"},
        {"feature": "ENST4", "symbol": "OTHER", "hgvsc": "ENST4:c.9A>G"},
        {"feature": "ENSR5", "symbol": "-", "consequence": "regulatory_region_variant"},
    ]}}
@case
def discordant_with_symbol_free_rows(spec):
    spec["vep_extra"] = {"GRCh37": {"gof000000": [
        {"feature": "ENST1.1", "symbol": "OTHER", "hgvsc": "ENST1.1:c.1A>G"},
        {"feature": "ENST2.1", "symbol": "-", "hgvsc": "ENST2.1:c.1A>G"}]}, "GRCh38": {"gof000000": []}}
@case
def vep_empty_for_grch37(spec):
    spec["vep_extra"] = {"GRCh37": {f"gof{at:06d}": [] for at in range(12)}}
@case
def vep_empty_for_both(spec):
    empty = {f"gof{at:06d}": [] for at in range(12)}
    spec["vep_extra"] = {"GRCh37": empty, "GRCh38": empty}
@case
def vep_writes_nothing(spec): spec["control"]["vep_no_output"] = True
@case
def bcftools_fails(spec): spec["control"]["fail"] = ["bcftools"]
@case
def vep_fails(spec): spec["control"]["fail"] = ["vep"]
@case
def crossmap_fails(spec): spec["control"]["fail"] = ["CrossMap"]
@case
def bcftools_writes_nothing(spec): spec["control"]["silent"] = ["bcftools"]
@case
def crossmap_writes_nothing(spec): spec["control"]["silent"] = ["CrossMap"]
@case
def bcftools_missing(spec): spec["hide"] = ["bcftools"]
@case
def vep_missing(spec): spec["hide"] = ["vep"]
@case
def crossmap_missing(spec): spec["hide"] = ["CrossMap"]
@case
def hgnc_missing(spec):
    spec["drop_args"] = ["--hgnc-table"]
    spec["args"] = ["--hgnc-table", "absent.tsv"]
@case
def hgnc_ambiguous_alias(spec):
    spec["hgnc"] = spec["hgnc"] + [("HGNC:12", "SECOND", "OLDNAME", "")]
@case
def hgnc_previous_symbol(spec):
    spec["rows"][0]["GeneSymbol"] = "prev2"
    spec["rows"][1]["GeneSymbol"] = "PREV1"
    spec["rows"][0]["vep"] = "other"
@case
def hgnc_unknown_symbols(spec):
    spec["hgnc"] = [("", "GENE1", "", "")]
@case
def reviews_absent(spec): spec["reviews"] = None
@case
def reviews_not_supplied(spec): spec["drop_args"] = ["--mechanism-review-tsv"]
@case
def reviews_missing_columns(spec):
    spec["review_columns"] = ["source_gene", "allele_key", "gof_eligibility"]
@case
def reviews_header_only(spec): spec["reviews"] = []
def review(gene, position, chrom="1", **values):
    item = row(gene, chrom, position)
    key = f"SNV|{chrom}|{position}|{position}|{item['Ref']}|{item['Alt']}"
    return {"source_gene": gene, "allele_key": key, "gof_eligibility": "quarantined", **values}, item
@case
def reviews_every_state(spec):
    pairs = [review("GENE1", 1500 + at, reviewed_mechanism=name, reason_code=reason)
             for at, (name, reason) in enumerate([
                 ("LOF", "lof_reason"), ("mixed", "mix"), ("DOMINANT_NEGATIVE", ""), ("DN", "dn"),
                 ("UNCERTAIN", "u"), ("EXCLUDE", "e"), ("SOMETHING", "s"), ("", ""), ("NA", "NA")])]
    pairs.append(review("GENE1", 1520, reviewed_mechanism="LOF", gof_eligibility="Quarantined"))
    pairs.append(review("GENE1", 1521, reviewed_mechanism="LOF", gof_eligibility="eligible"))
    pairs.append(review("OTHER", 1522, reviewed_mechanism="LOF"))
    spec["rows"] = spec["rows"] + [item for _, item in pairs[:-1]]
    spec["rows"].append({**pairs[-1][1], "GeneSymbol": "GENE1", "vep": "GENE1"})
    spec["reviews"] = [values for values, _ in pairs]
@case
def review_of_discordant_gene(spec):
    values, item = review("GENE8", 1600, chrom="2", reviewed_mechanism="LOF")
    spec["rows"].append({**item, "vep": "OTHER"})
    spec["reviews"] = [values]
@case
def hg19_fasta_missing(spec):
    spec["drop_args"] = ["--hg19-fasta"]
    spec["args"] = ["--hg19-fasta", "absent.fa"]
@case
def hg38_fasta_missing(spec):
    spec["drop_args"] = ["--hg38-fasta"]
    spec["args"] = ["--hg38-fasta", "absent.fa"]
@case
def fasta_without_index(spec): spec["fai"] = False
@case
def fasta_with_duplicate_contig(spec):
    def setup(folder):
        path = folder / "hg19.fa.fai"
        path.write_text(path.read_text() + path.read_text().splitlines()[0] + "\n")
    spec["setup"] = setup
@case
def argument_missing(spec): spec["drop_args"] = ["--workdir", "--chain"]
@case
def argument_out_missing(spec): spec["drop_args"] = ["--out-json"]
@case
def argument_bad_timeout(spec): spec["args"] = ["--api-timeout", "soon"]
@case
def argument_unknown(spec): spec["args"] = ["--unknown", "x"]


def outcome(command, folder, env):
    result = subprocess.run([str(part) for part in command], cwd=folder, capture_output=True, text=True, env=env)
    # pyfaidx's __del__ complains after a failed open; the failure is what came before.
    lines = result.stderr.split("Exception ignored in")[0].splitlines()
    value = {"exit": result.returncode, "stdout": result.stdout,
             "log": [line for line in lines if line.startswith("[gofcards]")]}
    if result.returncode == 2 and lines:
        value["message"] = lines[-1].split("error: ", 1)[-1]
        value["message"] = value["message"].removeprefix("priva: ")
    if result.returncode == 1 and lines:
        last = lines[-1]
        match = re.match(r"priva: ([A-Za-z]+): (.*)", last) or re.match(
            r"(?:[\w.]+\.)?([A-Za-z]+(?:Error|Exception|Exit)?): (.*)", last)
        value["category"], value["message"] = (match[1], match[2]) if match else ("?", last)
    return value


def normalize(value, folder, token):
    if isinstance(value, dict):
        return {key: normalize(item, folder, token) for key, item in value.items()}
    if isinstance(value, list):
        return [normalize(item, folder, token) for item in value]
    if isinstance(value, str):
        value = value.replace(str(folder), "<FOLDER>").replace(token, "<TOKEN>")
        return re.sub(r"(127\.0\.0\.1):\d+", r"\1:PORT", value)
    return value


def snapshot(folder, token):
    files = {}
    for path in sorted(folder.rglob("*")):
        key = str(path.relative_to(folder))
        if key.startswith("server") or path.is_dir():
            continue
        raw = path.read_bytes()
        if raw[:2] == b"\x1f\x8b":
            text = gzip.decompress(raw).decode()
            files[key] = {"text": normalize(text, folder, token),
                          "header": raw[:4].hex() + raw[8:raw.index(b"\0", 10) + 1 if raw[3] else 10].hex(),
                          "stream_sha256": hashlib.sha256(raw[raw.index(b"\0", 10) + 1 if raw[3] else 10:]).hexdigest()
                          if token not in text else "url in text"}
        else:
            text = raw.decode("utf-8", "replace")
            if key.endswith(".xlsx"):
                text = "served bytes" if raw == (folder / "server" / "public.xlsx").read_bytes() else "other bytes"
            files[key] = normalize(text, folder, token)
        files[key + ":mode"] = stat.S_IMODE(path.stat().st_mode)
    return files


def tool_folder(root, hide):
    folder = root / ("tools_" + "_".join(hide))
    folder.mkdir(parents=True, exist_ok=True)
    program = root / "stand_in_tool"
    if not program.exists():
        program.write_text(TOOL)
        program.chmod(0o755)
    for name in ("bcftools", "vep", "CrossMap"):
        if name not in hide and not (folder / name).exists():
            (folder / name).symlink_to(program)
    for name in ("curl", "samtools"):
        found = shutil.which(name)
        if found and not (folder / name).exists():
            (folder / name).symlink_to(found)
    return folder


def run_case(name, binary, repo, root, port):
    spec = base_spec()
    CASES[name](spec)
    tools = tool_folder(root / "tools", spec["hide"])
    driver = root / "driver.py"
    arms = {}
    for arm in ("python", "native"):
        folder = root / name / arm
        token = f"{name}-{arm}"
        base_url = f"http://127.0.0.1:{port}/{token}"
        Handler.folders[token] = folder
        command = materialize(folder, copy.deepcopy(spec), base_url, tools)
        endpoint = f"{base_url}/summary"
        env = dict(os.environ, PATH=str(tools), PYTHONDONTWRITEBYTECODE="1", TZ="UTC",
                   no_proxy="*", NO_PROXY="*")
        env.pop("http_proxy", None)
        if arm == "python":
            full = [sys.executable, driver, repo / "scripts", endpoint, *command]
        else:
            full = [binary, "gofcards-exact-gof", "--summary-endpoint", endpoint, *command]
        result = outcome(full, folder, env)
        arms[arm] = {"result": normalize(result, folder, token), "files": snapshot(folder, token)}
    python, native = arms["python"], arms["native"]
    return name, python == native, python["result"]["exit"], arms


def check(args):
    root = Path(args.work_root)
    root.mkdir(parents=True)
    (root / "driver.py").write_text(DRIVER)
    tool_folder(root / "tools", [])
    server = start_server()
    port = server.server_address[1]
    binary = Path(args.binary).resolve()
    repo = Path(args.repo).resolve()
    names = [name for name in CASES if not args.cases or name in args.cases.split(",")]
    failed, exits = [], {}
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(run_case, name, binary, repo, root, port) for name in names]
        for future in futures:
            name, same, status, arms = future.result()
            exits[status] = exits.get(status, 0) + 1
            if same:
                continue
            failed.append(name)
            print(f"DIFFERENT {name}", file=sys.stderr)
            python, native = arms["python"], arms["native"]
            for field in sorted(set(python["result"]) | set(native["result"])):
                left, right = python["result"].get(field), native["result"].get(field)
                if left != right:
                    if isinstance(left, list) and isinstance(right, list):
                        pairs = [(a, b) for a, b in zip(left, right) if a != b][:3]
                        print(f"  {field} lengths {len(left)} {len(right)}: {pairs}"[:1500], file=sys.stderr)
                    else:
                        print(f"  {field}: python {str(left)[:600]!r} native {str(right)[:600]!r}", file=sys.stderr)
            for key in sorted(set(python["files"]) | set(native["files"])):
                left, right = python["files"].get(key), native["files"].get(key)
                if left != right:
                    at = next((i for i, (a, b) in enumerate(zip(str(left), str(right))) if a != b), 0)
                    print(f"  file {key} differs at {at}: python {str(left)[max(at - 60, 0):at + 120]!r} "
                          f"native {str(right)[max(at - 60, 0):at + 120]!r}", file=sys.stderr)
    print(json.dumps({"cases": len(names), "exit_statuses": exits, "failed": failed}))
    return 1 if failed else 0


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def real(args):
    root = Path(args.work_root)
    root.mkdir(parents=True)
    (root / "driver.py").write_text(DRIVER)
    server = start_server()
    port = server.server_address[1]
    binary = Path(args.binary).resolve()
    repo = Path(args.repo).resolve()
    folders = {arm: root / arm for arm in ("python", "native")}
    runs = {}

    def run(arm):
        folder = folders[arm]
        (folder / "server").mkdir(parents=True)
        (folder / "work").mkdir()
        shutil.copy(args.workbook, folder / "server" / "public.xlsx")
        (folder / "server" / "control.json").write_text('{"summary": {}}')
        shutil.copy(args.summary_cache, folder / "work" / "gofcards_summary_cache.jsonl")
        Handler.folders[arm] = folder
        base_url = f"http://127.0.0.1:{port}/{arm}"
        command = ["--workdir", "work", "--out-json", "work/cache.json.gz", "--hg19-fasta", args.hg19_fasta,
                   "--hg38-fasta", args.hg38_fasta, "--chain", args.chain, "--vep-cache-dir", args.vep_cache_dir,
                   "--vep-cache-version", args.vep_cache_version, "--hgnc-table", args.hgnc_table,
                   "--mechanism-review-tsv", args.mechanism_review_tsv, "--public-excel-url",
                   "http://127.0.0.1:%d/%s/public.xlsx" % (port, arm), "--stats-json", "work/stats.json"]
        endpoint = f"{base_url}/summary"
        if arm == "python":
            full = [sys.executable, root / "driver.py", repo / "scripts", endpoint, *command]
        else:
            full = [binary, "gofcards-exact-gof", "--summary-endpoint", endpoint, *command]
        started = time.time()
        env = dict(os.environ, no_proxy="*", NO_PROXY="*")
        result = outcome(full, folder, env)
        runs[arm] = {"seconds": round(time.time() - started, 1), "result": normalize(result, folder, arm)}

    with concurrent.futures.ThreadPoolExecutor(2) as pool:
        for future in [pool.submit(run, arm) for arm in folders]:
            future.result()

    def view(folder, arm):
        work = folder / "work"
        stats = json.loads((work / "stats.json").read_text())
        cache = json.loads(gzip.decompress((work / "cache.json.gz").read_bytes()))
        return {
            "stats": normalize(stats, folder, arm),
            "cache_sha256": hashlib.sha256(json.dumps(normalize(cache, folder, arm), sort_keys=True).encode()).hexdigest(),
            "rejects": (work / "gofcards_reference_rejects.tsv").read_text(),
            "vcf_sha256": {name: sha256(work / name) for name in
                           ("gofcards.hg19.vcf", "gofcards.hg19.norm.vcf", "gofcards.hg38.vcf",
                            "gofcards.hg38.norm.vcf", "gofcards.hg38.vcf.unmap")},
            "summary_sha256": sha256(work / "gofcards_summary_cache.jsonl"),
            "requests": (folder / "requests.log").read_text().count("\n"),
        }
    views = {arm: view(folder, arm) for arm, folder in folders.items()}
    raw = {arm: (folders[arm] / "work" / "cache.json.gz").read_bytes() for arm in folders}
    same_log = runs["python"]["result"]["log"] == runs["native"]["result"]["log"]
    report = {"same": views["python"] == views["native"] and same_log, "same_log": same_log,
              "python_seconds": runs["python"]["seconds"], "native_seconds": runs["native"]["seconds"],
              "exits": [runs["python"]["result"]["exit"], runs["native"]["result"]["exit"]],
              "stats": views["native"]["stats"], "requests": views["native"]["requests"],
              "same_vcf": views["python"]["vcf_sha256"] == views["native"]["vcf_sha256"],
              "same_cache": views["python"]["cache_sha256"] == views["native"]["cache_sha256"],
              "same_rejects": views["python"]["rejects"] == views["native"]["rejects"],
              "gzip_size": [len(raw["python"]), len(raw["native"])]}
    print(json.dumps(report))
    return 0 if report["same"] else 1


def main():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("check", "real"):
        part = commands.add_parser(name)
        part.add_argument("--binary", required=True)
        part.add_argument("--repo", required=True)
        part.add_argument("--work-root", required=True)
        if name == "check":
            part.add_argument("--jobs", type=int, default=4)
            part.add_argument("--cases", default="", help="comma-separated case names (default: all)")
        else:
            for option in ("workbook", "summary-cache", "hg19-fasta", "hg38-fasta", "chain", "vep-cache-dir",
                           "vep-cache-version", "hgnc-table", "mechanism-review-tsv"):
                part.add_argument("--" + option, required=True)
    args = parser.parse_args()
    return check(args) if args.command == "check" else real(args)


if __name__ == "__main__":
    sys.exit(main())
